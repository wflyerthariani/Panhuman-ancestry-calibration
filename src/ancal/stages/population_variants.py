"""Stage 7: normalized population variants (item 7, protocol Steps 6-8).

Per contig in scope, in chunks of intervals, over the gnomAD HGDP+1KG callset:

  region query                                     raw records to a temporary file (tabix for local
                                                   files, HTTP range requests for remote ones)
  view -S <post-QC samples>                        genotypes of eligible samples only
  annotate -x INFO,^FORMAT/GT                      drop everything except genotypes and FILTER
  norm -m -any -f GRCh38 -c w                      split multiallelics, left-align, check REF
  +fill-tags -S <sample groups>                    AN, AC, het/hom/hemi counts and NS for every
                                                   population, group and sample set, from genotypes
  view -G                                          drop genotypes; only site-level counts reach disk

Sample sets (Step 4): `unrelated` is primary; `all_eligible` (unrelated + related post-QC
samples) is kept for sensitivity analysis. Genotypes are never stored.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from .. import budget, config, consequence, provenance, status
from ..settings import Settings
from ..sources import remote_vcf
from . import common
from .reference import PRIMARY_CONTIGS, read_contig

TAGS = ("AN", "AC", "AC_Hom", "AC_Het", "AC_Hemi", "NS")
SETS = {"u": "unrelated", "a": "all_eligible"}
CHUNK_INTERVALS = 20   # scope intervals per bcftools call: a failed call repeats only its own chunk
WORKERS = 4
RETRIES = 4


@dataclass(frozen=True)
class Group:
    token: str          # suffix of the INFO tags, e.g. u_gAFR
    sample_set: str     # unrelated | all_eligible
    level: str          # all | group | population
    label: str
    population_group: str
    n_samples: int


def sample_groups(samples: pl.DataFrame) -> tuple[list[Group], dict[str, list[str]]]:
    """Group definitions and, per sample, the tokens it belongs to."""
    post = samples.filter(pl.col("post_qc") & pl.col("in_vcf")).sort("sample_id")
    pops = sorted(post["original_population"].unique().to_list())
    pop_tok = {p: f"p{i:03d}" for i, p in enumerate(pops, start=1)}
    pop_group = dict(zip(post["original_population"], post["population_group"]))
    membership: dict[str, list[str]] = {}
    counts: dict[str, int] = {}
    for sid, pop, grp, incl in post.select("sample_id", "original_population", "population_group", "inclusion_status").iter_rows():
        toks = []
        for s in ("u", "a") if incl == "primary" else ("a",):
            toks += [f"{s}_all", f"{s}_g{grp}", f"{s}_{pop_tok[pop]}"]
        membership[sid] = toks
        for t in toks:
            counts[t] = counts.get(t, 0) + 1
    groups = []
    for s, set_name in SETS.items():
        groups.append(Group(f"{s}_all", set_name, "all", "ALL", "", counts.get(f"{s}_all", 0)))
        for g in sorted(set(pop_group.values())):
            groups.append(Group(f"{s}_g{g}", set_name, "group", g, g, counts.get(f"{s}_g{g}", 0)))
        for pop in pops:
            groups.append(Group(f"{s}_{pop_tok[pop]}", set_name, "population", pop, pop_group[pop],
                                counts.get(f"{s}_{pop_tok[pop]}", 0)))
    return [g for g in groups if g.n_samples > 0], membership


def variant_type(ref: str, alt: str) -> str:
    if alt == "*" or alt.startswith("<"):
        return "star_or_symbolic"
    if len(ref) == 1 and len(alt) == 1:
        return "snv"
    return "mnv" if len(ref) == len(alt) else "indel"


# ------------------------------------------------------------------- bcftools


def extract_contig(settings: Settings, uri: str, index_suffix: str, contig: str, bed: pl.DataFrame,
                   samples_file: Path, groups_file: Path, scope: str) -> tuple[Path, str]:
    """Run the pipeline for one contig, in resumable chunks. Returns the sites-only VCF and the norm summary."""
    d = settings.data_dir("variants")
    d.mkdir(parents=True, exist_ok=True)
    dest = d / f"hgdp_1kg.{scope}.{contig}.sites.vcf.gz"
    log = d / f"hgdp_1kg.{scope}.{contig}.norm.log"
    if provenance.recorded_row(settings, dest, uri) and log.exists():
        return dest, log.read_text()
    tbi = d / f"{Path(uri).name}{index_suffix}"
    provenance.fetch(settings, uri + index_suffix, tbi, stage="stage7", source_name=f"gnomAD HGDP+1KG tabix index {contig}",
                     database_version="v3.1.2", reference_build="GRCh38")
    fasta = settings.data_dir("reference") / f"GRCh38.{contig}.fa.gz"
    chunks = d / "chunks" / f"{scope}.{contig}"
    chunks.mkdir(parents=True, exist_ok=True)
    iv = bed.filter(pl.col("chrom") == contig).select("chrom", "start", "end").sort("start")
    # Two steps per chunk: copy the raw records of the chunk's intervals to a temporary file, then
    # reduce them to site counts and delete the raw copy. Remote files are read with HTTP range
    # requests (ancal.sources.remote_vcf): htslib's own remote reads of this bucket are reset mid-region.
    local = "://" not in uri or uri.startswith("file://")
    path = uri.removeprefix("file://")
    fetch = f"tabix -h -R {{bed}} '{path}##idx##{tbi}' | gzip -1 > {{raw}}" if local else \
        f"ancal.sources.remote_vcf.write_regions {uri} {{bed}} > {{raw}}"
    head = "" if local else remote_vcf.header(uri)
    reduce = (
        f"bcftools view --no-version -S {samples_file} {{raw}} -Ou"
        " | bcftools annotate --no-version -x INFO,^FORMAT/GT -Ou"
        f" | bcftools norm --no-version -m -any --old-rec-tag OLD_REC -f {fasta} -c w -Ou 2> {{log}}"
        f" | bcftools +fill-tags --no-version -Ou -- -S {groups_file} -t {','.join(TAGS)}"
        " | bcftools view --no-version -G -Ob -o {out}"
    )
    cmd = f"{fetch} && {reduce}  # per chunk of {CHUNK_INTERVALS} intervals, then bcftools concat"

    def run_chunk(i: int) -> Path:
        """One chunk of intervals -> a sites-only BCF. Skipped if already complete; retried on network errors."""
        part, part_log, part_bed = chunks / f"{i:05d}.bcf", chunks / f"{i:05d}.log", chunks / f"{i:05d}.bed"
        if part.exists() and part_log.exists():
            return part
        iv.slice(i * CHUNK_INTERVALS, CHUNK_INTERVALS).write_csv(part_bed, separator="\t", include_header=False)
        raw, tmp, tmp_log = chunks / f"{i:05d}.raw.vcf.gz", part.with_suffix(".bcf.part"), part_log.with_suffix(".log.part")
        part_iv = iv.slice(i * CHUNK_INTERVALS, CHUNK_INTERVALS)
        err = ""
        try:
            for attempt in range(RETRIES):
                try:
                    if local:
                        subprocess.run(["bash", "-c", "set -o pipefail; " + fetch.format(bed=part_bed, raw=raw)],
                                       cwd=chunks, capture_output=True, text=True, check=True)
                    else:
                        remote_vcf.write_regions(uri, tbi, contig, list(zip(part_iv["start"], part_iv["end"])), raw, head)
                    subprocess.run(["bash", "-c", "set -o pipefail; " + reduce.format(raw=raw, log=tmp_log, out=tmp)],
                                   cwd=chunks, capture_output=True, text=True, check=True)
                    tmp_log.replace(part_log)
                    tmp.replace(part)
                    return part
                except subprocess.CalledProcessError as e:
                    err = e.stderr.strip()[-600:]
                except (OSError, ValueError) as e:   # network errors are OSError subclasses in requests
                    err = f"{type(e).__name__}: {e}"
                time.sleep(5 * 2 ** attempt)
        finally:
            raw.unlink(missing_ok=True)
        raise RuntimeError(f"chunk {i} of {contig} failed {RETRIES} times: {err}")

    n_chunks = -(-iv.height // CHUNK_INTERVALS)
    with ThreadPoolExecutor(WORKERS) as pool:
        parts = list(pool.map(run_chunk, range(n_chunks)))
    budget.check(settings)
    listing = chunks / "parts.txt"
    listing.write_text("".join(f"{x}\n" for x in parts))
    subprocess.run(["bcftools", "concat", "--no-version", "-f", str(listing), "-Oz", "-o", str(dest)], check=True,
                   capture_output=True)
    totals: dict[str, int] = {}
    for x in parts:
        for k, v in norm_counts(x.with_suffix(".log").read_text()).items():
            totals[k] = totals.get(k, 0) + v
    log.write_text(f"Lines   {'/'.join(totals)}:\t{'/'.join(map(str, totals.values()))}\n")
    shutil.rmtree(chunks)
    provenance.record(settings, dest, stage="stage7", source_name=f"gnomAD HGDP+1KG sites with per-population counts ({contig})",
                      command=cmd, software_name="bcftools", software_version=common.tool_version("bcftools"),
                      input_uri=uri, access_method="remote_tabix", database_version="v3.1.2", reference_build="GRCh38",
                      input_info=provenance.head(uri))
    budget.check(settings)
    return dest, log.read_text()


def read_sites_vcf(path: Path, groups: list[Group]) -> pl.DataFrame:
    """Wide table: one row per record, one column per (tag, group)."""
    cols = [f"{t}_{g.token}" for g in groups for t in TAGS]
    fmt = "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/OLD_REC" + "".join(f"\t%INFO/{c}" for c in cols) + "\n"
    tsv = path.with_suffix(".query.tsv")
    try:
        with open(tsv, "w") as fh:
            subprocess.run(["bcftools", "query", "-f", fmt, str(path)], stdout=fh, check=True)
        schema = {"chrom": pl.String, "pos": pl.Int64, "ref": pl.String, "alt": pl.String, "source_filter": pl.String,
                  "old_rec": pl.String, **{c: pl.Int64 for c in cols}}
        return pl.read_csv(tsv, separator="\t", has_header=False, new_columns=list(schema), schema_overrides=schema,
                           null_values=["."])
    finally:
        tsv.unlink(missing_ok=True)


def norm_counts(log: str) -> dict[str, int]:
    """Parse bcftools norm's `Lines   total/split/...: 1/2/...` summary."""
    m = re.search(r"^Lines\s+([\w/]+):\s*([\d/]+)", log, re.M)
    if not m:
        return {}
    return dict(zip(m.group(1).split("/"), map(int, m.group(2).split("/"))))


# --------------------------------------------------------------------- tables


def build_tables(wide: pl.DataFrame, groups: list[Group], seq: str, region_pos: pl.DataFrame, cons: pl.DataFrame,
                 cfg: config.AnalysisConfig) -> tuple[pl.DataFrame, pl.DataFrame]:
    """(sites, long) for one contig from the wide per-record table."""
    key = ["chrom", "pos", "ref", "alt"]
    count_cols = [c for c in wide.columns if c not in (*key, "source_filter", "old_rec")]
    wide = wide.with_columns(pl.col(count_cols).fill_null(0))
    # Exact duplicates (same alleles, same counts) are removed; same alleles with different counts conflict.
    dup = wide.group_by(key).agg(pl.len().alias("n_records"), pl.struct(count_cols).n_unique().alias("n_distinct"))
    wide = wide.unique(subset=key, keep="first", maintain_order=True).join(dup, on=key)

    refs = [seq[p - 1 : p - 1 + len(r)] == r.upper() for p, r in zip(wide["pos"], wide["ref"])]
    types = [variant_type(r, a) for r, a in zip(wide["ref"], wide["alt"])]
    sites = wide.select(
        pl.format("{}:{}:{}:{}", *key).alias("variant_id"), *key, pl.Series("variant_type", types),
        pl.col("source_filter").fill_null("."),
        pl.when(pl.col("old_rec").is_null()).then(pl.lit("unchanged"))
        .when(pl.col("old_rec").str.split("|").list.get(3, null_on_oob=True).str.contains(","))
        .then(pl.lit("split_multiallelic")).otherwise(pl.lit("realigned")).alias("normalization"),
        pl.Series("ref_matches_grch38", refs), (pl.col("n_records") - 1).alias("n_duplicate_records_removed"),
        (pl.col("n_distinct") > 1).alias("allele_conflict"),
        pl.col("AC_u_all").alias("ac_global"), pl.col("AN_u_all").alias("an_global"),
        pl.col("AC_a_all").alias("ac_all_eligible"), pl.col("AN_a_all").alias("an_all_eligible"),
    ).join(region_pos.with_columns(pl.lit(True).alias("in_scope_region")), on=["chrom", "pos"], how="left").with_columns(
        pl.col("in_scope_region").fill_null(False),
        pl.when(pl.col("an_global") > 0).then(pl.col("ac_global") / pl.col("an_global")).alias("af_global"),
        pl.when(pl.col("an_all_eligible") > 0).then(pl.col("ac_all_eligible") / pl.col("an_all_eligible")).alias("af_all_eligible"),
    ).with_columns(
        pl.when(~pl.col("ref_matches_grch38")).then(pl.lit("reference_allele_mismatch"))
        .when(pl.col("allele_conflict")).then(pl.lit("unresolved_allele_conflict"))
        .when(pl.col("variant_type") == "star_or_symbolic").then(pl.lit("alleles_cannot_be_harmonized"))
        .when(pl.col("source_filter") != "PASS").then(pl.lit("fails_source_filter"))
        .otherwise(pl.lit("")).alias("qc_exclusion_reason"),
        common.af_bin_expr(cfg),
    ).with_columns(
        pl.when(pl.col("qc_exclusion_reason") == "").then(pl.lit("pass")).otherwise(pl.lit("excluded")).alias("qc_status"),
    )
    ann = consequence.per_variant(cons)
    sites = sites.join(ann, on=key, how="left").with_columns(
        pl.col("consequence").is_not_null().alias("coding"),
        pl.when(pl.col("consequence").is_not_null()).then(pl.format("{}{}{}", "ref_aa", "protein_pos", "alt_aa")).alias("protein_change"),
    ).with_columns(
        ((pl.col("variant_type") == "snv") & (pl.col("qc_status") == "pass") & pl.col("in_scope_region")).alias("primary_snv"))

    ids = sites.select("variant_id", pl.col("qc_exclusion_reason").alias("_site_reason"))
    base = wide.select(pl.format("{}:{}:{}:{}", *key).alias("variant_id"), *count_cols)
    parts = []
    for g in groups:
        t = g.token
        parts.append(base.select(
            "variant_id", pl.lit(g.sample_set).alias("sample_set"), pl.lit(g.level).alias("level"),
            pl.lit(g.label).alias("label"), pl.lit(g.population_group).alias("population_group"),
            pl.lit(g.n_samples).alias("n_samples"), pl.col(f"AC_{t}").alias("ac"), pl.col(f"AN_{t}").alias("an"),
            pl.col(f"NS_{t}").alias("n_contributing_individuals"), pl.col(f"AC_Het_{t}").alias("n_het"),
            (pl.col(f"AC_Hom_{t}") // 2).alias("n_hom_alt"), pl.col(f"AC_Hemi_{t}").alias("n_hemi")))
    long = pl.concat(parts).join(ids, on="variant_id").with_columns(
        pl.when(pl.col("an") > 0).then(pl.col("ac") / pl.col("an")).alias("af"),
        (pl.col("n_contributing_individuals") / pl.col("n_samples")).alias("call_rate"),
        (pl.col("an") > 0).alias("callable"),
        pl.when(pl.col("_site_reason") != "").then(pl.col("_site_reason"))
        .when(pl.col("an") == 0).then(pl.lit("allele_number_zero")).otherwise(pl.lit("")).alias("qc_exclusion_reason"),
    ).with_columns(
        pl.when(pl.col("qc_exclusion_reason") == "").then(pl.lit("pass")).otherwise(pl.lit("excluded")).alias("qc_status"),
        ((pl.col("ac") >= 1) & (pl.col("qc_exclusion_reason") == "")).alias("observed"),
    ).drop("_site_reason")

    obs = long.filter(pl.col("observed") & (pl.col("sample_set") == "unrelated"))
    n_obs = obs.group_by("variant_id").agg(
        (pl.col("level") == "population").sum().alias("n_populations_observed"),
        (pl.col("level") == "group").sum().alias("n_groups_observed"))
    sites = sites.join(n_obs, on="variant_id", how="left").with_columns(
        pl.col("n_populations_observed", "n_groups_observed").fill_null(0))
    return sites, long


def summarize(sites: pl.DataFrame, long: pl.DataFrame, cfg: config.AnalysisConfig) -> pl.DataFrame:
    """Observed-variant counts per sample set and population (primary SNVs only beyond the first column)."""
    j = long.filter(pl.col("observed")).join(
        sites.select("variant_id", "primary_snv", "coding", "consequence", "af_bin"), on="variant_id")
    snv = pl.col("primary_snv")
    aggs = [pl.len().alias("n_observed_variants"), snv.sum().alias("n_observed_snvs"),
            (snv & pl.col("coding")).sum().alias("n_coding_snvs")]
    aggs += [(snv & (pl.col("consequence") == c)).sum().alias(f"n_{c}") for c in ("missense", "synonymous", "nonsense")]
    aggs += [(snv & (pl.col("af_bin") == b.name)).sum().alias(f"n_snvs_{b.name}") for b in cfg.allele_frequency_bins.bins]
    levels = long.select("sample_set", "level", "label", "population_group", "n_samples").unique()
    order = {"all": 0, "group": 1, "population": 2}
    return (levels.join(j.group_by("sample_set", "level", "label").agg(aggs), on=["sample_set", "level", "label"], how="left")
            .fill_null(0).with_columns(pl.col("level").replace_strict(order, return_dtype=pl.Int32).alias("_o"))
            .sort(pl.col("sample_set") != "unrelated", "_o", "population_group", "label").drop("_o"))


def consistency(long: pl.DataFrame) -> tuple[int, int]:
    """(rows with AN above 2 x samples, variants whose population ACs do not sum to the group/overall AC)."""
    over = long.filter(pl.col("an") > 2 * pl.col("n_samples")).height
    pop = long.filter(pl.col("level") == "population").group_by("variant_id", "sample_set", "population_group").agg(pl.col("ac").sum())
    grp = long.filter(pl.col("level") == "group").select("variant_id", "sample_set", "population_group", pl.col("ac").alias("ac_group"))
    bad = pop.join(grp, on=["variant_id", "sample_set", "population_group"]).filter(pl.col("ac") != pl.col("ac_group")).height
    tot = long.filter(pl.col("level") == "group").group_by("variant_id", "sample_set").agg(pl.col("ac").sum())
    allr = long.filter(pl.col("level") == "all").select("variant_id", "sample_set", pl.col("ac").alias("ac_all"))
    bad += tot.join(allr, on=["variant_id", "sample_set"]).filter(pl.col("ac") != pl.col("ac_all")).height
    return over, bad


# ------------------------------------------------------------------------- run


def run(settings: Settings) -> status.StageStatus:
    budget.check(settings)
    ds = config.load_data_sources(settings.root)
    cfg = config.load_analysis(settings.root)
    p = common.Paths(settings)
    notes = "Normalized population-variant table (item 7)."
    for need, stage in ((p.bed, 4), (p.samples, 3), (p.consequences, 6)):
        if not need.exists():
            return status.write(settings, "stage7", [], notes=notes,
                                blocked_reason=f"{need.name} is missing; run `make stage{stage}` first.")
    bed = common.read_bed(p.bed)
    region_pos = common.bed_positions(bed)
    cons = pl.read_parquet(p.consequences)
    samples = pl.read_parquet(p.samples)
    groups, membership = sample_groups(samples)
    d = settings.data_dir("variants")
    d.mkdir(parents=True, exist_ok=True)
    samples_file, groups_file = d / "samples.post_qc.txt", d / "sample_groups.txt"
    samples_file.write_text("".join(f"{s}\n" for s in membership))
    groups_file.write_text("".join(f"{s}\t{','.join(t)}\n" for s, t in membership.items()))

    src = ds.sources["population_vcf"]
    contigs = [c for c in PRIMARY_CONTIGS if c in set(bed["chrom"])]
    site_parts, long_parts, norm, n_records = [], [], {}, 0
    for contig in contigs:
        vcf, log = extract_contig(settings, src.uri.format(contig=contig), src.index_suffix or ".tbi", contig, bed,
                                  samples_file, groups_file, ds.scope)
        for k, v in norm_counts(log).items():
            norm[k] = norm.get(k, 0) + v
        wide = read_sites_vcf(vcf, groups)
        n_records += wide.height
        seq = read_contig(settings.data_dir("reference") / f"GRCh38.{contig}.fa.gz", contig)
        sites_c, long_c = build_tables(wide, groups, seq, region_pos.filter(pl.col("chrom") == contig),
                                       cons.filter(pl.col("chrom") == contig), cfg)
        sp, lp = d / f"sites.{ds.scope}.{contig}.part.parquet", d / f"long.{ds.scope}.{contig}.part.parquet"
        sites_c.write_parquet(sp)
        long_c.write_parquet(lp)
        site_parts.append(sp)
        long_parts.append(lp)
        del wide, seq, sites_c, long_c
        budget.check(settings)
    pl.scan_parquet(site_parts).sink_parquet(p.sites)
    pl.scan_parquet(long_parts).sink_parquet(p.long)
    for part in (*site_parts, *long_parts):
        part.unlink()

    sites, long = pl.read_parquet(p.sites), pl.read_parquet(p.long)
    summary = summarize(sites, long, cfg)
    excl = (sites.filter(pl.col("qc_status") == "excluded").group_by("qc_exclusion_reason", "variant_type").len()
            .rename({"len": "n_variants"}).sort(["n_variants", "qc_exclusion_reason"], descending=[True, False]))
    pop_excl = (long.filter((pl.col("qc_exclusion_reason") == "allele_number_zero")).group_by("sample_set", "level", "label").len()
                .rename({"len": "n_variants_allele_number_zero"}).sort("sample_set", "level", "label"))
    t_sum = common.write_tsv(summary, p.tables / "population_variant_summary.tsv")
    t_excl = common.write_tsv(excl, p.exclusions / "population_variant_exclusions.tsv")
    t_an0 = common.write_tsv(pop_excl, p.exclusions / "population_variant_allele_number_zero.tsv")
    common.record(settings, "stage7", "population_variants", "variants/hgdp_1kg.*.sites.vcf.gz", p.sites, p.long,
                  t_sum, t_excl, t_an0)

    over, bad = consistency(long)
    n_u = next(g.n_samples for g in groups if g.token == "u_all")
    n_a = next(g.n_samples for g in groups if g.token == "a_all")
    pass_snv = sites.filter(pl.col("primary_snv"))
    mismatch = sites.filter(~pl.col("ref_matches_grch38")).height
    checks = [
        status.Check("samples", n_a == len(membership) and n_u > 0, f"{n_u:,} unrelated (primary), {n_a:,} all eligible; "
                     f"{sum(g.level == 'population' and g.sample_set == 'all_eligible' for g in groups)} populations"),
        status.Check("sites_extracted", sites.height > 0, f"{n_records:,} records after normalization; {sites.height:,} unique "
                     f"variants on {', '.join(contigs)}"),
        status.Check("unique_variant_ids", sites["variant_id"].n_unique() == sites.height, "one row per chrom:pos:ref:alt"),
        status.Check("long_table_complete", long.height == sites.height * len(groups),
                     f"{long.height:,} rows = {sites.height:,} variants x {len(groups)} population/group/set levels"),
        status.Check("allele_number_within_bounds", over == 0, f"{over} rows with AN > 2 x samples"),
        status.Check("counts_consistent", bad == 0, f"{bad} variants whose population counts do not sum to the group or overall count"),
        status.Check("reference_alleles", mismatch == 0 or mismatch / sites.height < 0.001,
                     f"{mismatch} variants with a REF that differs from GRCh38 (excluded and listed)"),
        status.Check("primary_snvs", pass_snv.height > 0, f"{pass_snv.height:,} PASS SNVs in scope regions; "
                     f"{pass_snv.filter(pl.col('ac_global') >= 1).height:,} observed in the unrelated set"),
    ]
    u = budget.usage(settings)
    checks.append(status.Check("disk_budget", u.ok, f"{u.used_mb:.1f} MB of {settings.local_budget_gb} GB"))
    write_report(settings, p, ds.scope, sites, long, summary, excl, norm, n_records, groups, contigs, src.uri)
    return status.write(settings, "stage7", checks, notes=notes)


def write_report(settings, p, scope, sites, long, summary, excl, norm, n_records, groups, contigs, uri):
    by_type = sites.group_by("variant_type", "qc_status").len().sort("variant_type", "qc_status")
    by_norm = sites.group_by("normalization").len().sort("normalization")
    by_filter = sites.group_by("source_filter").len().sort("len", descending=True).head(15)
    conflicts = sites.filter(pl.col("allele_conflict") | ~pl.col("ref_matches_grch38")).select(
        "variant_id", "variant_type", "source_filter", "ref_matches_grch38", "allele_conflict")
    prim = sites.filter(pl.col("primary_snv") & (pl.col("ac_global") >= 1))
    by_bin = prim.group_by("af_bin").len().sort("af_bin")
    by_cons = prim.group_by(pl.col("consequence").fill_null("noncoding")).len().sort("len", descending=True)
    top = summary.filter((pl.col("sample_set") == "unrelated") & (pl.col("level") != "population")).drop("sample_set")
    text = f"""# Variant normalization and allele-conflict report (Stage 7, protocol Steps 6-8)

Generated by `ancal stage7`. Do not edit. Scope: `{scope}` ({', '.join(contigs)}).

{common.provisional_banner(scope)}Source: `{uri}`, read by region. Genotypes are used in the stream and never stored.

## Normalization (Step 6)

| | Count |
|---|---|
| Source records in the scope regions (bcftools norm `total`) | {common.fmt(norm.get('total'))} |
| Multiallelic records split | {common.fmt(norm.get('split'))} |
| Records left-aligned | {common.fmt(norm.get('realigned'))} |
| Records after normalization | {n_records:,} |
| Exact duplicate records removed | {int(sites['n_duplicate_records_removed'].sum()):,} |
| Unique variants (`chrom:pos:ref:alt`) | {sites.height:,} |
| REF differs from GRCh38 | {sites.filter(~pl.col('ref_matches_grch38')).height:,} |
| Same alleles with conflicting counts | {sites.filter(pl.col('allele_conflict')).height:,} |

Chromosome names already follow `chr1..chr22, chrX`. Variant IDs are `chromosome:position:reference:alternate`.

{common.md_table(by_norm, ['Normalization', 'Variants'])}

## Variant types and QC status

All types stay in the master table. Only SNVs that pass QC and lie inside the scope regions
(`primary_snv`) enter the primary DMS analysis. An SNV split from a multiallelic site counts as an SNV.

{common.md_table(by_type, ['Type', 'QC status', 'Variants'])}

### Site-level exclusions (Step 7)

{common.md_table(excl, ['Reason', 'Type', 'Variants']) if excl.height else 'None.'}

Source FILTER values (top 15):

{common.md_table(by_filter, ['FILTER', 'Variants'])}

### Allele-conflict report

{common.md_table(conflicts.head(50)) if conflicts.height else 'No reference-allele mismatches and no conflicting duplicate records.'}

## Population-level quantities (Step 7)

For each of {len(groups)} levels (overall, 7 groups and 80 populations, for the unrelated and the
all-eligible sample sets), `{p.long.name}` holds: alternate allele count, allele number, allele
frequency, number of contributing individuals, call rate, heterozygotes, alternate homozygotes,
hemizygotes, observed, callable, QC status and exclusion reason.

- **Observed** = alternate allele count of at least 1 after QC.
- **Callable** = allele number above 0. No per-population callability mask is distributed (L-03).
- **Call rate** = individuals with a non-missing genotype / individuals in the group.
- **Global allele frequency** (`af_global`, Step 16) = unrelated samples of the combined resource.

## Observed primary SNVs (unrelated set)

{common.md_table(top)}

By global allele-frequency bin:

{common.md_table(by_bin, ['Bin', 'SNVs'])}

By consequence (MANE Select):

{common.md_table(by_cons, ['Consequence', 'SNVs'])}

All 80 populations and the all-eligible sample set: `tables/population_variant_summary.tsv`.

## Outputs

- `data_store/variants/{p.sites.name}`: one row per variant (type, QC, global frequency and bin, consequence).
- `data_store/variants/{p.long.name}`: one row per variant, population level and sample set. The table for
  one population or group (Step 8) is the rows with that `level` and `label`.
"""
    (p.qc / "variant_normalization_report.md").write_text(text)
