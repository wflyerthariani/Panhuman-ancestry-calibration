"""Stage 4: reference and annotation resources (scope-aware).

Scope-independent (always the same, about 15 MB):
  - GENCODE v50 basic, MANE Select transcripts and CDS only, genome-wide
    -> reference/gencode.v50.basic.mane_select.gtf.gz and tables of transcripts and CDS

Scope-dependent (depend on resources/regions/<scope>.bed):
  - regions BED: MANE Select CDS of the scope's genes
  - GRCh38 FASTA for the scope's contigs only (uppercased, bgzipped, faidx-indexed)
  - AlphaMissense hg38 rows inside the regions (protein -> genome lookup only; D-008)
  - ClinVar records inside the regions, from the dated release, chromosomes renamed to chrN
"""

from __future__ import annotations

import bisect
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from .. import budget, config, provenance, status
from ..settings import Settings

PRIMARY_CONTIGS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
GRCH38_LENGTHS = {  # GRCh38 primary assembly (checked against the FASTA and the VCF header)
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555, "chr5": 181538259,
    "chr6": 170805979, "chr7": 159345973, "chr8": 145138636, "chr9": 138394717, "chr10": 133797422,
    "chr11": 135086622, "chr12": 133275309, "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
    "chr16": 90338345, "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895,
}
MANE_MIN_TRANSCRIPTS = 19000  # MANE v1.x has about 19,300 Select transcripts
_ATTR = re.compile(r'(\w+) "([^"]*)"')


class ScopeBlocked(RuntimeError):
    pass


# ------------------------------------------------------------ GENCODE / MANE


def gtf_keep(line: str) -> bool:
    """Keep MANE Select transcript and CDS rows on the primary contigs."""
    if line.startswith("#"):
        return False
    parts = line.split("\t", 3)
    return (
        len(parts) == 4 and parts[0] in GRCH38_LENGTHS and parts[2] in ("transcript", "CDS")
        and 'tag "MANE_Select"' in line
    )


def parse_mane_gtf(path: Path) -> tuple[pl.DataFrame, pl.DataFrame]:
    rows = []
    import gzip

    with gzip.open(path, "rt") as fh:
        for line in fh:
            chrom, _, feature, start, end, _, strand, frame, attrs = line.rstrip("\n").split("\t")
            a = dict(_ATTR.findall(attrs))
            rows.append({
                "chrom": chrom, "feature": feature, "start": int(start), "end": int(end), "strand": strand,
                "frame": frame, "gene_name": a.get("gene_name"), "gene_id": a.get("gene_id"),
                "transcript_id": a.get("transcript_id"), "protein_id": a.get("protein_id"),
                "gene_type": a.get("gene_type"), "exon_number": a.get("exon_number"),
            })
    df = pl.DataFrame(rows)
    cds = df.filter(pl.col("feature") == "CDS").drop("feature")
    tx = (
        df.filter(pl.col("feature") == "transcript")
        .select("gene_name", "gene_id", "transcript_id", "chrom", "strand", "gene_type",
                pl.col("start").alias("tx_start"), pl.col("end").alias("tx_end"))
        .join(
            cds.group_by("transcript_id").agg(
                pl.col("protein_id").first(), pl.col("start").min().alias("cds_start"),
                pl.col("end").max().alias("cds_end"), pl.len().alias("n_cds_segments"),
                (pl.col("end") - pl.col("start") + 1).sum().alias("cds_length"),
            ),
            on="transcript_id", how="left",
        )
        .with_columns((pl.col("cds_length") % 3 == 0).alias("cds_length_multiple_of_3"))
        .sort("chrom", "tx_start")
    )
    return tx, cds


# --------------------------------------------------------------------- scope


def scope_genes(ds: config.DataSources, tx: pl.DataFrame, root: Path) -> list[str]:
    genes = ds.scopes[ds.scope].genes
    if isinstance(genes, list):
        return genes
    if genes == "all_mane_select":
        return sorted(tx["gene_name"].unique().to_list())
    if genes == "from_dms_inventory":
        inv = root / "tables" / "dms_candidate_inventory.tsv"
        if not inv.exists():
            raise ScopeBlocked(
                "scope dms_genes takes its genes from tables/dms_candidate_inventory.tsv (Stage 5), which does "
                "not exist yet. Run Stage 5 once (the inventory is scope-independent), then re-run Stage 4."
            )
        cand = pl.read_csv(inv, separator="\t", infer_schema_length=0).filter(pl.col("prescreen_status") == "candidate")
        return sorted(cand["gene_symbol"].drop_nulls().unique().to_list())
    raise ValueError(f"unknown gene spec {genes!r} for scope {ds.scope}")


def regions_bed(cds: pl.DataFrame, genes: list[str]) -> tuple[pl.DataFrame, list[str]]:
    """Merged 0-based BED of MANE Select CDS for `genes`; also returns genes not found."""
    sel = cds.filter(pl.col("gene_name").is_in(genes))
    missing = sorted(set(genes) - set(sel["gene_name"].unique().to_list()))
    iv = sel.select("chrom", (pl.col("start") - 1).alias("start"), "end", "gene_name").sort("chrom", "start")
    merged = []
    for chrom, start, end, gene in iv.iter_rows():
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            m = merged[-1]
            merged[-1] = (chrom, m[1], max(m[2], end), m[3] if gene in m[3].split(",") else f"{m[3]},{gene}")
        else:
            merged.append((chrom, start, end, gene))
    bed = pl.DataFrame(merged, schema=["chrom", "start", "end", "name"], orient="row")
    order = {c: i for i, c in enumerate(PRIMARY_CONTIGS)}
    bed = bed.with_columns(pl.col("chrom").replace_strict(order, return_dtype=pl.Int32).alias("_o")).sort("_o", "start").drop("_o")
    return bed, missing


class IntervalIndex:
    """Point-in-interval lookup on a merged, sorted BED (0-based half-open)."""

    def __init__(self, bed: pl.DataFrame):
        self._starts: dict[str, list[int]] = {}
        self._ends: dict[str, list[int]] = {}
        for chrom, start, end, _ in bed.iter_rows():
            self._starts.setdefault(chrom, []).append(start)
            self._ends.setdefault(chrom, []).append(end)

    def contains(self, chrom: str, pos1: int) -> bool:
        starts = self._starts.get(chrom)
        if not starts:
            return False
        i = bisect.bisect_right(starts, pos1 - 1) - 1
        return i >= 0 and pos1 - 1 < self._ends[chrom][i]


# ----------------------------------------------------------------- resources


def fetch_reference(settings: Settings, uri_template: str, contig: str) -> Path:
    dest = settings.data_dir("reference") / f"GRCh38.{contig}.fa.gz"
    uri = uri_template.format(contig=contig)
    if not provenance.recorded_row(settings, dest, uri):
        provenance.stream_filter(
            settings, uri, dest, transform=lambda l: l if l.startswith(">") else l.upper(),
            stage="stage4", source_name=f"UCSC hg38 {contig}", reference_build="GRCh38",
            database_version="hg38 (UCSC per-chromosome)", filter_description="uppercase soft-masked bases; bgzip",
            bgzip_output=True,
        )
    fai = Path(f"{dest}.fai")
    if not fai.exists():
        subprocess.run(["samtools", "faidx", str(dest)], check=True)
        for idx in (fai, Path(f"{dest}.gzi")):
            provenance.record(settings, idx, stage="stage4", source_name=f"faidx index {contig}",
                              command=f"samtools faidx {dest.name}", software_name="samtools",
                              software_version=_version("samtools"), input_uri=dest.name, reference_build="GRCh38")
    return dest


def read_contig(fasta: Path, contig: str) -> str:
    out = subprocess.run(["samtools", "faidx", str(fasta), contig], capture_output=True, text=True, check=True).stdout
    return "".join(out.split("\n")[1:])


def subset_alphamissense(settings: Settings, uri: str, index: IntervalIndex, scope: str) -> Path:
    dest = settings.data_dir("reference") / f"alphamissense_hg38.{scope}.tsv.gz"
    if provenance.recorded_row(settings, dest, uri):
        return dest

    def keep(line: str) -> bool:
        if line.startswith("#"):
            return True
        chrom, pos, _ = line.split("\t", 2)
        return index.contains(chrom, int(pos))

    provenance.stream_filter(
        settings, uri, dest, keep, stage="stage4", source_name="AlphaMissense hg38", reference_build="GRCh38",
        database_version="AlphaMissense 2023 (CC BY-NC-SA 4.0)",
        filter_description=f"rows inside resources/regions/{scope}.bed", gzip_output=True,
    )
    return dest


@dataclass
class ClinVarRelease:
    date: str        # YYYYMMDD
    uri: str
    info: provenance.RemoteInfo


def resolve_clinvar(src: config.Source) -> ClinVarRelease:
    date = src.release
    if not date or date == "latest":
        header = subprocess.run(["bcftools", "view", "-h", src.uri], capture_output=True, text=True, check=True).stdout
        m = re.search(r"^##fileDate=(\d{4})-?(\d{2})-?(\d{2})", header, re.M)
        if not m:
            raise RuntimeError("could not read ##fileDate from the latest ClinVar VCF")
        date = "".join(m.groups())
    for template in src.dated_uris or []:
        uri = template.format(date=date, year=date[:4])
        try:
            info = provenance.head(uri)
            provenance.head(uri + (src.index_suffix or ".tbi"))
            return ClinVarRelease(date, uri, info)
        except Exception:
            continue
    raise RuntimeError(f"no dated ClinVar file found for {date}")


def subset_clinvar(settings: Settings, src: config.Source, bed: pl.DataFrame, scope: str) -> tuple[Path, ClinVarRelease]:
    rel = resolve_clinvar(src)
    d = settings.data_dir("clinvar")
    tbi = d / f"clinvar_{rel.date}.vcf.gz.tbi"
    provenance.fetch(settings, rel.uri + ".tbi", tbi, stage="stage4", source_name="ClinVar tabix index",
                     database_version=rel.date, reference_build="GRCh38")

    # ClinVar uses 1..22,X; the rest of the pipeline uses chr1..chr22,chrX.
    nochr_bed = d / f"regions.{scope}.nochr.bed"
    bed.with_columns(pl.col("chrom").str.strip_prefix("chr")).write_csv(nochr_bed, separator="\t", include_header=False)
    rename = d / "chr_rename.txt"
    rename.write_text("".join(f"{c.removeprefix('chr')}\t{c}\n" for c in PRIMARY_CONTIGS))

    dest = d / f"clinvar_{rel.date}.{scope}.vcf.gz"
    src_arg = f"{rel.uri}##idx##{tbi}"
    cmd = (f"bcftools view --no-version -R {nochr_bed} '{src_arg}' | "
           f"bcftools annotate --no-version --rename-chrs {rename} -Oz -o {dest} && tabix -f -p vcf {dest}")
    if not provenance.recorded_row(settings, dest, rel.uri):
        subprocess.run(cmd, shell=True, check=True, cwd=d)
        provenance.record(settings, dest, stage="stage4", source_name="ClinVar GRCh38 (region subset)",
                          command=cmd, software_name="bcftools", software_version=_version("bcftools"),
                          input_uri=rel.uri, access_method="remote_tabix", database_version=rel.date,
                          reference_build="GRCh38", input_info=rel.info)
        provenance.record(settings, Path(f"{dest}.tbi"), stage="stage4", source_name="tabix index",
                          command=f"tabix -p vcf {dest.name}", software_name="tabix",
                          software_version=_version("tabix"), input_uri=dest.name)
    after = provenance.head(rel.uri)
    if after.etag != rel.info.etag or after.last_modified != rel.info.last_modified:
        raise RuntimeError("ClinVar dated file changed while it was being read")
    return dest, rel


def _version(tool: str) -> str:
    return subprocess.run([tool, "--version"], capture_output=True, text=True).stdout.splitlines()[0]


# ---------------------------------------------------------------- validation


def check_alphamissense_ref(path: Path, seqs: dict[str, str]) -> tuple[int, int, int]:
    """(rows, rows whose REF matches the reference base, rows on contigs not loaded)."""
    df = pl.read_csv(path, separator="\t", comment_prefix="#", has_header=False,
                     new_columns=["chrom", "pos", "ref", "alt", "genome", "uniprot_id", "transcript_id",
                                  "protein_variant", "am_pathogenicity", "am_class"],
                     schema_overrides={"pos": pl.Int64})
    ok = skipped = 0
    for chrom, pos, ref in df.select("chrom", "pos", "ref").iter_rows():
        seq = seqs.get(chrom)
        if seq is None:
            skipped += 1
        elif seq[pos - 1] == ref:
            ok += 1
    return df.height, ok, skipped


def check_clinvar(path: Path, seqs: dict[str, str]) -> dict:
    out = subprocess.run(["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", str(path)],
                         capture_output=True, text=True, check=True).stdout
    n = snv = ok = 0
    chroms = set()
    for line in out.splitlines():
        chrom, pos, ref, alt = line.split("\t")
        n += 1
        chroms.add(chrom)
        seq = seqs.get(chrom)
        if seq is not None and seq[int(pos) - 1 : int(pos) - 1 + len(ref)] == ref.upper():
            ok += 1
        if len(ref) == 1 and len(alt) == 1:
            snv += 1
    return {"records": n, "snv": snv, "ref_match": ok, "chroms": sorted(chroms)}


# ----------------------------------------------------------------------- run


def run(settings: Settings) -> status.StageStatus:
    budget.check(settings)
    ds = config.load_data_sources(settings.root)
    checks: list[status.Check] = []

    # --- GENCODE MANE Select (scope-independent)
    g = ds.sources["gene_annotation"]
    gtf = settings.data_dir("reference") / f"gencode.{g.version}.basic.mane_select.gtf.gz"
    if not provenance.recorded_row(settings, gtf, g.uri):
        provenance.stream_filter(settings, g.uri, gtf, gtf_keep, stage="stage4", source_name="GENCODE basic",
                                 database_version=g.version, reference_build="GRCh38",
                                 filter_description="MANE_Select transcript+CDS rows on chr1-22,X", gzip_output=True)
    tx, cds = parse_mane_gtf(gtf)
    tx_path = settings.data_dir("reference") / "mane_select_transcripts.parquet"
    cds_path = settings.data_dir("reference") / "mane_select_cds.parquet"
    tx.write_parquet(tx_path)
    cds.write_parquet(cds_path)
    n_tx, n_genes = tx.height, tx["gene_name"].n_unique()
    frac3 = tx["cds_length_multiple_of_3"].mean()
    checks += [
        status.Check("mane_select_transcripts", n_tx >= MANE_MIN_TRANSCRIPTS, f"{n_tx} transcripts, {n_genes} genes"),
        status.Check("mane_one_transcript_per_gene", tx["gene_id"].n_unique() == n_tx, "one MANE Select per gene"),
        status.Check("mane_cds_multiple_of_3", frac3 >= 0.999, f"{frac3:.4%} of CDS lengths are multiples of 3"),
    ]

    # --- scope regions
    try:
        genes = scope_genes(ds, tx, settings.root)
    except ScopeBlocked as e:
        return status.write(settings, "stage4", checks, notes="Reference and annotation.", blocked_reason=str(e))
    bed, missing = regions_bed(cds, genes)
    bed_path = settings.root / "resources" / "regions" / f"{ds.scope}.bed"
    bed_path.parent.mkdir(parents=True, exist_ok=True)
    bed.write_csv(bed_path, separator="\t", include_header=False)
    contigs = bed["chrom"].unique(maintain_order=True).to_list()
    bp = int((bed["end"] - bed["start"]).sum())
    checks.append(status.Check("scope_genes_found", not missing,
                               f"{len(genes) - len(missing)}/{len(genes)} genes; {bed.height} intervals, {bp:,} bp on "
                               f"{', '.join(contigs)}" + (f"; missing {missing}" if missing else "")))

    # --- reference FASTA for the scope's contigs
    seqs = {}
    ref_src = ds.sources["reference_fasta"]
    for contig in contigs:
        fasta = fetch_reference(settings, ref_src.uri, contig)
        seqs[contig] = read_contig(fasta, contig)
        ok_len = len(seqs[contig]) == GRCH38_LENGTHS[contig]
        alphabet = set(seqs[contig]) <= set("ACGTN")
        checks.append(status.Check(f"reference_{contig}", ok_len and alphabet,
                                   f"length {len(seqs[contig]):,} (GRCh38 {GRCH38_LENGTHS[contig]:,}); "
                                   f"alphabet {''.join(sorted(set(seqs[contig])))}"))
    budget.check(settings)

    # --- AlphaMissense subset
    index = IntervalIndex(bed)
    am = subset_alphamissense(settings, ds.sources["alphamissense"].uri, index, ds.scope)
    n_am, am_ok, am_skip = check_alphamissense_ref(am, seqs)
    checks.append(status.Check("alphamissense_ref_matches_grch38", n_am > 0 and am_ok == n_am - am_skip,
                               f"{n_am:,} rows; REF matches reference for {am_ok:,}"))

    # --- ClinVar subset
    cv, rel = subset_clinvar(settings, ds.sources["clinvar"], bed, ds.scope)
    cvs = check_clinvar(cv, seqs)
    checks.append(status.Check("clinvar_subset", cvs["records"] > 0 and cvs["ref_match"] == cvs["records"]
                               and set(cvs["chroms"]) <= set(contigs),
                               f"release {rel.date}; {cvs['records']:,} records ({cvs['snv']:,} SNVs); REF matches "
                               f"for {cvs['ref_match']:,}; contigs {cvs['chroms']}"))

    u = budget.usage(settings)
    checks.append(status.Check("disk_budget", u.ok, f"{u.used_mb:.1f} MB of {settings.local_budget_gb} GB"))

    for p in (tx_path, cds_path, bed_path):
        provenance.record(settings, p, stage="stage4", source_name="derived", command="ancal stage4",
                          software_name="ancal.stages.reference", software_version=provenance.code_version(settings.root),
                          input_uri=gtf.name)
    write_report(settings, ds.scope, genes, bed, tx, n_am, am_ok, cvs, rel, u)
    return status.write(settings, "stage4", checks,
                        notes=f"Reference and annotation for scope {ds.scope}: ClinVar release {rel.date}.")


def write_report(settings, scope, genes, bed, tx, n_am, am_ok, cvs, rel, usage):
    by_gene = (tx.filter(pl.col("gene_name").is_in(genes))
               .select("gene_name", "gene_id", "transcript_id", "protein_id", "chrom", "strand",
                       "cds_start", "cds_end", "n_cds_segments", "cds_length",
                       (pl.col("cds_length") // 3).alias("protein_length_aa")))
    lines = "\n".join("| " + " | ".join(str(v) for v in r) + " |" for r in by_gene.iter_rows())
    text = f"""# Reference and annotation report (Stage 4)

Generated by `ancal stage4`. Do not edit. Scope: `{scope}`. Run times are in the manifest.

## Resources

| Resource | Version | Stored | Check |
|---|---|---|---|
| GENCODE basic, MANE Select transcript and CDS rows (genome-wide) | v50 | `reference/gencode.v50.basic.mane_select.gtf.gz` | {tx.height:,} transcripts; {tx['cds_length_multiple_of_3'].mean():.3%} of CDS lengths are multiples of 3 |
| GRCh38 FASTA, scope contigs only | UCSC hg38 | `reference/GRCh38.<contig>.fa.gz` (+ .fai, .gzi) | lengths equal GRCh38; uppercased |
| AlphaMissense hg38, rows in scope regions | 2023 release (**CC BY-NC-SA 4.0**) | `reference/alphamissense_hg38.{scope}.tsv.gz` | {n_am:,} rows; REF equals GRCh38 for {am_ok:,} |
| ClinVar GRCh38, records in scope regions | **{rel.date}** (dated file) | `clinvar/clinvar_{rel.date}.{scope}.vcf.gz` | {cvs['records']:,} records ({cvs['snv']:,} SNVs); REF equals GRCh38 for {cvs['ref_match']:,}; chromosomes renamed to `chrN` |

ClinVar was read from `{rel.uri}`. That dated file does not change after release. The
permanent copy will be under `archive_2.0/{rel.date[:4]}/`.

## Scope regions (`resources/regions/{scope}.bed`)

{bed.height} merged MANE Select CDS intervals covering {int((bed['end'] - bed['start']).sum()):,} bp.

| Gene | Ensembl gene | MANE Select transcript | Protein | Chrom | Strand | CDS start | CDS end | CDS segments | CDS length (nt) | Protein length (aa) |
|---|---|---|---|---|---|---|---|---|---|---|
{lines}

Local data used after this stage: {usage.used_mb:.1f} MB of {settings.local_budget_gb} GB.
"""
    (settings.root / "reports" / "qc" / "reference_report.md").write_text(text)
