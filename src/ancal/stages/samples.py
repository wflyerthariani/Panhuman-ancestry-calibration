"""Stage 3: sample metadata, population mapping and counts (items 4, 5, 6; protocol Steps 3-5).

Sample set: the Koenig et al. (2024) post-QC HGDP+1KG release, i.e. the 4,094 samples
in the global PCA "without outliers". gnomAD's own `high_quality` flag is NOT used: it
removes every Biaka, Mbuti, San, Bougainville and Papuan sample, which would delete the
OCE group entirely.

Relatedness: a post-QC sample is related if it is in Koenig's related_sample_ids table.
This reproduces the sample counts of their released unrelated (3,400) and related (694)
PCA sets.
"""

from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path

import polars as pl
import yaml

from .. import budget, config, provenance, status
from ..settings import Settings
from ..sources.hail_table import read_string_key_table

N_PCS = 20
PCS_FOR_DISCORDANCE = 10
# 1000 Genomes populations described as recently admixed. Flagged for review, never excluded (Step 5).
RECENTLY_ADMIXED = ("ACB", "ASW", "CLM", "MXL", "PEL", "PUR")
GRCH38_CONTIG_LENGTHS = {"chr1": 248956422, "chr17": 83257441, "chrX": 156040895}
EXPECTED_SAMPLES = 4094
EXPECTED_POPULATIONS = 80
META = {
    "sample_id": "s",
    "source_dataset": "hgdp_tgp_meta.Project",
    "original_population": "population",
    "gnomad_population_label": "hgdp_tgp_meta.Population",
    "population_group": "hgdp_tgp_meta.Genetic.region",
    "geographic_region": "hgdp_tgp_meta.Study.region",
    "latitude": "hgdp_tgp_meta.Latitude",
    "longitude": "hgdp_tgp_meta.Longitude",
    "sex_karyotype": "sex_imputation.sex_karyotype",
    "reported_sex": "sex",
    "hard_filtered": "sample_filters.hard_filtered",
    "hard_filters": "sample_filters.hard_filters",
}


# ------------------------------------------------------------------ pure logic


def read_metadata(path: Path) -> pl.DataFrame:
    raw = pl.read_csv(path, separator="\t", infer_schema_length=0, null_values=["NA", ""])
    df = raw.select([pl.col(src).alias(dst) for dst, src in META.items()])
    return df.with_columns(pl.col("latitude").cast(pl.Float64), pl.col("longitude").cast(pl.Float64))


def read_pcs(path: Path) -> pl.DataFrame:
    import gzip
    import io

    data = gzip.decompress(path.read_bytes()) if path.name.endswith((".gz", ".bgz")) else path.read_bytes()
    df = pl.read_csv(io.BytesIO(data), separator="\t", infer_schema_length=0)
    pcs = [f"PC{i}" for i in range(1, N_PCS + 1)]
    return df.select(
        pl.col("s").alias("sample_id"), *[pl.col(p).cast(pl.Float64).alias(p.lower()) for p in pcs]
    )


def build_sample_table(
    meta: pl.DataFrame, vcf_samples: list[str], pcs: pl.DataFrame, related: set[str], outliers: set[str]
) -> pl.DataFrame:
    """One row per individual in the VCF (Step 4), including excluded ones with a reason."""
    post_qc = set(pcs["sample_id"])
    vcf = pl.DataFrame({"sample_id": vcf_samples, "in_vcf": [True] * len(vcf_samples)})
    df = vcf.join(meta, on="sample_id", how="full", coalesce=True).with_columns(
        pl.col("in_vcf").fill_null(False),
        pl.col("sample_id").is_in(meta["sample_id"].to_list()).alias("in_metadata"),
        pl.col("sample_id").is_in(list(post_qc)).alias("post_qc"),
        pl.col("sample_id").is_in(list(related)).alias("in_related_list"),
    )
    df = df.with_columns(
        pl.when(~pl.col("post_qc")).then(pl.lit("not_assessed"))
        .when(pl.col("in_related_list")).then(pl.lit("related"))
        .otherwise(pl.lit("unrelated")).alias("relatedness_status"),
        pl.when(pl.col("post_qc") & ~pl.col("in_related_list")).then(pl.lit("primary"))
        .when(pl.col("post_qc")).then(pl.lit("sensitivity_only"))
        .otherwise(pl.lit("excluded")).alias("inclusion_status"),
        pl.when(pl.col("post_qc")).then(pl.lit(None, dtype=pl.String))
        .when(~pl.col("in_metadata")).then(pl.lit("not_in_sample_metadata"))
        .when(~pl.col("in_vcf")).then(pl.lit("not_in_vcf"))
        .when(pl.col("hard_filtered") == "true")
        .then(pl.concat_str(pl.lit("gnomad_hard_filter:"), pl.col("hard_filters").str.strip_chars("[]")))
        .when(pl.col("sample_id").is_in(list(outliers))).then(pl.lit("pca_outlier"))
        .otherwise(pl.lit("not_in_koenig_post_qc_release")).alias("exclusion_reason"),
    )
    df = df.join(pcs, on="sample_id", how="left")
    df = add_pc_discordance(df)
    df = df.with_columns(pl.col("original_population").is_in(RECENTLY_ADMIXED).alias("recently_admixed_population"))
    cols = [
        "sample_id", "source_dataset", "original_population", "gnomad_population_label", "population_group",
        "geographic_region", "latitude", "longitude", "sex_karyotype", "reported_sex", "in_vcf", "in_metadata",
        "post_qc", "relatedness_status", "inclusion_status", "exclusion_reason", "recently_admixed_population",
        "pc_nearest_group", "pc_group_discordant", *[f"pc{i}" for i in range(1, N_PCS + 1)],
    ]
    return df.select(cols).sort("population_group", "original_population", "sample_id", nulls_last=True)


def add_pc_discordance(df: pl.DataFrame) -> pl.DataFrame:
    """Flag post-QC samples whose nearest group centroid (PCs 1-10) is not their assigned group.

    Review aid only (Step 5); nobody is excluded on this basis.
    """
    pcs = [f"pc{i}" for i in range(1, PCS_FOR_DISCORDANCE + 1)]
    scored = df.filter(pl.col("post_qc"))
    centroids = scored.group_by("population_group").agg([pl.col(p).median() for p in pcs]).sort("population_group")
    groups = centroids["population_group"].to_list()
    dist = [
        sum(((pl.col(p) - centroids[p][i]) ** 2) for p in pcs).sqrt().alias(f"_d_{g}")
        for i, g in enumerate(groups)
    ]
    scored = scored.with_columns(dist)
    dcols = [f"_d_{g}" for g in groups]
    scored = scored.with_columns(
        pl.concat_list(dcols).list.arg_min().alias("_argmin")
    ).with_columns(
        pl.col("_argmin").map_elements(lambda i: groups[i], return_dtype=pl.String).alias("pc_nearest_group")
    ).with_columns(
        (pl.col("pc_nearest_group") != pl.col("population_group")).alias("pc_group_discordant")
    )
    return df.join(scored.select("sample_id", "pc_nearest_group", "pc_group_discordant"), on="sample_id", how="left")


def population_mapping(samples: pl.DataFrame) -> pl.DataFrame:
    """One row per original population for team review (Step 5)."""
    labelled = samples.filter(pl.col("in_metadata"))
    return (
        labelled.group_by("original_population")
        .agg(
            pl.col("source_dataset").unique().sort().str.join("|").alias("source_dataset"),
            pl.col("population_group").unique().sort().str.join("|").alias("population_group"),
            pl.col("population_group").n_unique().alias("n_groups"),
            pl.col("gnomad_population_label").unique().sort().str.join("|").alias("gnomad_population_label"),
            pl.col("geographic_region").unique().sort().str.join("|").alias("geographic_region"),
            pl.len().alias("n_samples_in_metadata"),
            pl.col("post_qc").sum().alias("n_post_qc"),
            (pl.col("inclusion_status") == "primary").sum().alias("n_unrelated_post_qc"),
            pl.col("pc_group_discordant").fill_null(False).sum().alias("n_pc_group_discordant"),
            pl.col("recently_admixed_population").first().alias("recently_admixed"),
        )
        .with_columns(
            (pl.col("n_groups") > 1).alias("ambiguous_group_label"),
            pl.col("original_population").is_null().alias("missing_label"),
            (pl.col("gnomad_population_label") != pl.col("original_population")).alias("label_differs_from_gnomad"),
        )
        .drop("n_groups")
        .sort("population_group", "original_population")
    )


def population_counts(samples: pl.DataFrame, koenig: pl.DataFrame) -> pl.DataFrame:
    """Counts per population, per group and in total (Step 3 / item 6), beside the published counts."""
    s = samples.filter(pl.col("in_metadata"))

    def agg(df: pl.DataFrame, keys: list[str]) -> pl.DataFrame:
        exprs = [
            pl.len().alias("n_in_metadata"),
            pl.col("in_vcf").sum().alias("n_in_vcf"),
            pl.col("post_qc").sum().alias("n_post_qc"),
            (pl.col("inclusion_status") == "primary").sum().alias("n_unrelated"),
            (pl.col("inclusion_status") == "sensitivity_only").sum().alias("n_related"),
            (pl.col("inclusion_status") == "excluded").sum().alias("n_excluded"),
            (pl.col("post_qc") & (pl.col("source_dataset") == "HGDP")).sum().alias("n_post_qc_hgdp"),
            (pl.col("post_qc") & (pl.col("source_dataset") == "1000 Genomes")).sum().alias("n_post_qc_1kg"),
            (pl.col("post_qc") & (pl.col("sex_karyotype") == "XX")).sum().alias("n_post_qc_xx"),
            (pl.col("post_qc") & (pl.col("sex_karyotype") == "XY")).sum().alias("n_post_qc_xy"),
        ]
        return df.group_by(keys).agg(exprs) if keys else df.select(exprs)

    pop = agg(s, ["population_group", "original_population"]).join(
        koenig, on=["population_group", "original_population"], how="full", coalesce=True
    ).with_columns(pl.lit("population").alias("level"))
    grp = agg(s, ["population_group"]).join(
        koenig.group_by("population_group").agg(pl.col("koenig_n_post_qc").sum(), pl.col("koenig_n_unrelated").sum()),
        on="population_group", how="left",
    ).with_columns(pl.lit("group").alias("level"), pl.lit(None, dtype=pl.String).alias("original_population"))
    tot = agg(s, []).with_columns(
        pl.lit("total").alias("level"), pl.lit(None, dtype=pl.String).alias("population_group"),
        pl.lit(None, dtype=pl.String).alias("original_population"),
        pl.lit(koenig["koenig_n_post_qc"].sum()).alias("koenig_n_post_qc"),
        pl.lit(koenig["koenig_n_unrelated"].sum()).alias("koenig_n_unrelated"),
    )
    order = ["level", "population_group", "original_population"]
    out = pl.concat([t.select(order + [c for c in pop.columns if c not in order]) for t in (tot, grp, pop)],
                    how="vertical_relaxed")
    out = out.with_columns(
        (pl.col("n_post_qc") - pl.col("koenig_n_post_qc")).alias("diff_post_qc_vs_koenig"),
        (pl.col("n_unrelated") - pl.col("koenig_n_unrelated")).alias("diff_unrelated_vs_koenig_summary"),
    )
    level_rank = pl.col("level").replace_strict({"total": 0, "group": 1, "population": 2}, return_dtype=pl.Int8)
    return out.sort(level_rank, "population_group", "original_population", nulls_last=True)


def read_koenig_summary(path: Path) -> pl.DataFrame:
    df = pl.read_csv(path, separator="\t", infer_schema_length=0)
    return df.select(
        pl.col("hgdp_tgp_meta.Genetic.region").alias("population_group"),
        pl.col("population").alias("original_population"),
        pl.col("n_snp_stats.n").cast(pl.Int64).alias("koenig_n_post_qc"),
        pl.col("n_unrelated").cast(pl.Int64).alias("koenig_n_unrelated"),
    )


# --------------------------------------------------------------------- inputs


@dataclass
class Inputs:
    meta: Path
    summary: Path
    pcs: Path
    outliers: Path
    related: Path
    vcf_header: Path
    unrelated_mt_meta: Path
    related_mt_meta: Path
    contig: str
    vcf_uri: str


def fetch_inputs(settings: Settings) -> Inputs:
    ds = config.load_data_sources(settings.root)
    src = ds.sources["sample_metadata"]
    d = settings.data_dir("population_metadata")
    kw = dict(stage="stage3", database_version="hgdp_1kg_v2")

    meta = d / "gnomad_meta_updated.tsv"
    provenance.fetch(settings, src.uri, meta, source_name="gnomAD HGDP+1KG sample metadata", **kw)
    summary = d / "post_qc_summary.tsv"
    provenance.fetch(settings, src.qc_summary_uri, summary, source_name="gnomAD HGDP+1KG post-QC summary", **kw)
    pcs = d / "GLOBAL_scores_without_outliers.txt.bgz"
    provenance.fetch(settings, src.pc_scores_uri, pcs, source_name="Koenig et al. global PC scores (post-QC)", **kw)
    outliers = d / "pca_outliers.txt"
    provenance.fetch(settings, src.pca_outliers_uri, outliers, source_name="Koenig et al. PCA outliers", **kw)
    unrel_meta = d / "unrelateds_without_outliers.cols.metadata.json.gz"
    provenance.fetch(settings, f"{src.unrelated_mt_uri}/cols/metadata.json.gz", unrel_meta,
                     source_name="Koenig et al. unrelated PCA set (column metadata)", **kw)
    rel_meta = d / "relateds_without_outliers.cols.metadata.json.gz"
    provenance.fetch(settings, f"{src.related_mt_uri}/cols/metadata.json.gz", rel_meta,
                     source_name="Koenig et al. related PCA set (column metadata)", **kw)

    related = d / "related_sample_ids.txt"
    ids, info = read_string_key_table(src.related_ids_uri)
    related.write_text("\n".join(sorted(ids)) + "\n")
    provenance.record(
        settings, related, stage="stage3", source_name="Koenig et al. related sample IDs (Hail Table)",
        command=f"ancal.sources.hail_table.read_string_key_table {src.related_ids_uri}",
        software_name="ancal.sources.hail_table", software_version=f"reads Hail {info['hail_version']}",
        input_uri=src.related_ids_uri, access_method="hail_table_decode", database_version="hgdp_1kg_v2",
        input_sha256=info["parts_sha256"],
    )

    contig = _first_contig(ds)
    vcf_uri = ds.sources["population_vcf"].uri.format(contig=contig)
    header = d / f"vcf_header_{contig}.txt"
    cmd = ["bcftools", "view", "--no-version", "-h", vcf_uri]
    header.write_text(subprocess.run(cmd, capture_output=True, text=True, check=True).stdout)
    provenance.record(
        settings, header, stage="stage3", source_name="gnomAD HGDP+1KG VCF header",
        command=" ".join(cmd), software_name="bcftools", software_version=_bcftools_version(),
        input_uri=vcf_uri, access_method="remote_header", database_version="v3.1.2",
        input_info=provenance.head(vcf_uri),
    )
    return Inputs(meta, summary, pcs, outliers, related, header, unrel_meta, rel_meta, contig, vcf_uri)


def _first_contig(ds: config.DataSources) -> str:
    contigs = ds.scopes[ds.scope].contigs
    return contigs[0] if isinstance(contigs, list) else "chr1"


def _bcftools_version() -> str:
    return subprocess.run(["bcftools", "--version"], capture_output=True, text=True).stdout.splitlines()[0]


def _mt_col_count(path: Path) -> int:
    import gzip
    import json

    return sum(json.loads(gzip.decompress(path.read_bytes()))["components"]["partition_counts"]["counts"])


def parse_vcf_header(path: Path) -> tuple[list[str], dict[str, int], str]:
    samples, contigs, reference = [], {}, ""
    for line in path.read_text().splitlines():
        if line.startswith("##contig=<"):
            fields = dict(kv.split("=", 1) for kv in line[len("##contig=<"):-1].split(",") if "=" in kv)
            if "ID" in fields and "length" in fields:
                contigs[fields["ID"]] = int(fields["length"])
        elif line.startswith("##reference="):
            reference = line.split("=", 1)[1]
        elif line.startswith("#CHROM"):
            samples = line.split("\t")[9:]
    return samples, contigs, reference


# ----------------------------------------------------------------------- run


def run(settings: Settings) -> status.StageStatus:
    budget.check(settings)
    inp = fetch_inputs(settings)

    meta = read_metadata(inp.meta)
    pcs = read_pcs(inp.pcs)
    related = set(inp.related.read_text().split())
    outliers = set(inp.outliers.read_text().split())
    vcf_samples, contig_lengths, reference = parse_vcf_header(inp.vcf_header)
    koenig = read_koenig_summary(inp.summary)
    n_unrel_mt, n_rel_mt = _mt_col_count(inp.unrelated_mt_meta), _mt_col_count(inp.related_mt_meta)

    samples = build_sample_table(meta, vcf_samples, pcs, related, outliers)
    mapping = population_mapping(samples)
    counts = population_counts(samples, koenig)

    # ---- write outputs
    tables = settings.root / "tables"
    samples.write_parquet(settings.data_dir("population_metadata") / "sample_metadata.parquet")
    samples.write_csv(tables / "sample_metadata.tsv", separator="\t", float_precision=6)
    mapping.write_csv(tables / "population_mapping.tsv", separator="\t")
    counts.write_csv(tables / "population_counts.tsv", separator="\t")
    groups_yaml = settings.root / "config" / "population_groups.yaml"
    groups = {r["original_population"]: r["population_group"] for r in mapping.iter_rows(named=True)}
    groups_yaml.write_text(
        "# GENERATED by Stage 3 (ancal stage3) from the HGDP+1KG metadata. Do not edit by hand.\n"
        "# Original population (80 labels) -> seven-group analysis label.\n"
        "# Review tables/population_mapping.tsv, then mark population_groups.mapping as\n"
        "# confirmed in config/analysis.yaml.\n"
        + yaml.safe_dump(dict(sorted(groups.items())), sort_keys=False)
    )

    post = samples.filter(pl.col("post_qc"))
    n_post, n_pops = post.height, post["original_population"].n_unique()
    n_unrel = post.filter(pl.col("inclusion_status") == "primary").height
    per_pop = counts.filter(pl.col("level") == "population")
    pop_mismatch = per_pop.filter(pl.col("diff_post_qc_vs_koenig").fill_null(1) != 0)
    unrel_summary_mismatch = per_pop.filter(pl.col("diff_unrelated_vs_koenig_summary") != 0)
    missing_labels = post.filter(pl.col("original_population").is_null() | pl.col("population_group").is_null())
    bad_groups = post.filter(~pl.col("population_group").is_in(config.GROUP_LABELS))
    multi_group = mapping.filter(pl.col("ambiguous_group_label"))
    not_in_vcf = post.filter(~pl.col("in_vcf"))
    expected_len = GRCH38_CONTIG_LENGTHS.get(inp.contig)

    checks = [
        status.Check("one_label_each", missing_labels.height == 0 and bad_groups.height == 0,
                     f"{missing_labels.height} post-QC samples missing a label; {bad_groups.height} with an unknown group"),
        status.Check("no_ambiguous_population", multi_group.height == 0,
                     f"{multi_group.height} populations map to more than one group"),
        status.Check("expected_sample_count", n_post == EXPECTED_SAMPLES, f"{n_post} post-QC samples (expected {EXPECTED_SAMPLES})"),
        status.Check("expected_population_count", n_pops == EXPECTED_POPULATIONS,
                     f"{n_pops} populations (expected {EXPECTED_POPULATIONS})"),
        status.Check("per_population_counts_match_koenig", pop_mismatch.height == 0,
                     f"{pop_mismatch.height} populations differ from the published post-QC counts"),
        status.Check("unrelated_matches_released_pca_sets", n_unrel == n_unrel_mt and n_post - n_unrel == n_rel_mt,
                     f"unrelated {n_unrel} vs {n_unrel_mt}; related {n_post - n_unrel} vs {n_rel_mt}"),
        status.Check("all_post_qc_samples_in_vcf", not_in_vcf.height == 0, f"{not_in_vcf.height} missing from the VCF header"),
        status.Check("grch38_contig_length", expected_len is None or contig_lengths.get(inp.contig) == expected_len,
                     f"{inp.contig} length {contig_lengths.get(inp.contig)} (GRCh38 {expected_len})"),
    ]

    review = []
    if unrel_summary_mismatch.height:
        review.append(
            f"Unrelated count {n_unrel} matches Koenig's released unrelated PCA set ({n_unrel_mt}) but not the "
            f"n_unrelated column of post_qc_summary.tsv ({koenig['koenig_n_unrelated'].sum()}); "
            f"{unrel_summary_mismatch.height} populations differ."
        )
    write_report(settings, samples, counts, mapping, inp, n_unrel_mt, n_rel_mt, reference, contig_lengths, review)

    for path in (tables / "sample_metadata.tsv", tables / "population_mapping.tsv", tables / "population_counts.tsv",
                 groups_yaml, settings.data_dir("population_metadata") / "sample_metadata.parquet",
                 settings.root / "reports" / "qc" / "sample_count_report.md"):
        provenance.record(settings, path, stage="stage3", source_name="derived", software_name="ancal.stages.samples",
                          software_version=provenance.code_version(settings.root), command="ancal stage3",
                          input_uri="population_metadata/{gnomad_meta_updated.tsv,post_qc_summary.tsv,"
                                    "GLOBAL_scores_without_outliers.txt.bgz,pca_outliers.txt,related_sample_ids.txt,"
                                    f"vcf_header_{inp.contig}.txt}}")

    notes = "Sample metadata, population mapping and counts (items 4-6). REVIEW STOP: see reports/qc/sample_count_report.md."
    if review:
        notes += " Review items: " + " | ".join(review)
    return status.write(settings, "stage3", checks, notes=notes)


# -------------------------------------------------------------------- report


def _md_table(df: pl.DataFrame) -> str:
    cols = df.columns
    rows = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for r in df.iter_rows():
        rows.append("| " + " | ".join("" if v is None else str(v) for v in r) + " |")
    return "\n".join(rows)


def write_report(settings, samples, counts, mapping, inp, n_unrel_mt, n_rel_mt, reference, contig_lengths, review):
    meta_rows = samples.filter(pl.col("in_metadata"))
    post = samples.filter(pl.col("post_qc"))
    excluded = samples.filter(pl.col("inclusion_status") == "excluded")
    total = counts.filter(pl.col("level") == "total").row(0, named=True)
    groups = counts.filter(pl.col("level") == "group").select(
        "population_group", "n_in_metadata", "n_post_qc", "n_unrelated", "n_related", "n_post_qc_hgdp",
        "n_post_qc_1kg", "koenig_n_post_qc", "koenig_n_unrelated")
    small = counts.filter((pl.col("level") == "population") & (pl.col("n_unrelated") < 10)).select(
        "population_group", "original_population", "n_post_qc", "n_unrelated")
    excl = excluded.group_by("exclusion_reason").agg(pl.len().alias("n"), pl.col("sample_id").sort().str.join(", ").alias("samples")).sort("n", descending=True)
    gnomad_hq_dropped = ["Biaka", "Mbuti", "San", "Bougainville", "PapuanHighlands", "PapuanSepik"]
    merged = mapping.filter(pl.col("label_differs_from_gnomad")).select("original_population", "gnomad_population_label", "population_group")
    disc = post.group_by("population_group").agg(pl.col("pc_group_discordant").sum().alias("n_discordant"), pl.len().alias("n")).sort("population_group")
    admixed = mapping.filter(pl.col("recently_admixed")).select("original_population", "population_group", "n_post_qc", "n_pc_group_discordant")

    text = f"""# Sample count report (protocol Step 3)

Generated by `ancal stage3`. Do not edit. Scope: `{settings.scope}`.

## Summary

| Check (Step 3) | Result |
|---|---|
| Samples in the VCF header ({inp.contig}) | {len([s for s in samples['in_vcf'].to_list() if s])} |
| Samples in the sample metadata | {meta_rows.height} |
| **Post-QC samples (Koenig et al. 2024 release)** | **{post.height}** (expected about 4,094) |
| HGDP post-QC samples | {total['n_post_qc_hgdp']} |
| 1000 Genomes post-QC samples | {total['n_post_qc_1kg']} |
| **Original populations** | **{post['original_population'].n_unique()}** (expected 80) |
| Seven-group categories | {post['population_group'].n_unique()} |
| Unrelated post-QC samples (primary set) | {total['n_unrelated']} |
| Related post-QC samples (sensitivity set only) | {total['n_related']} |
| Samples lacking population metadata | {samples.filter(~pl.col('in_metadata')).height} (in the VCF but not in the metadata) |
| Genome build | GRCh38 ({inp.contig} length {contig_lengths.get(inp.contig)} matches GRCh38){' - reference: ' + reference if reference else ''} |
| Variant representation | Per-chromosome bgzipped VCF with genotypes; normalization is checked in Stage 7 |
| Callability information | Not distributed with the callset. Stage 7 uses AN > 0 as a proxy (listed in item 15) |

## Counts by seven-group category

{_md_table(groups)}

Per-population counts are in `tables/population_counts.tsv`; the review table is `tables/population_mapping.tsv`.

## Differences from what was expected, and how they were resolved

1. **gnomAD's `high_quality` flag is not used.** In the metadata it marks only 3,942 samples, and it
   drops every sample from {', '.join(gnomad_hq_dropped)}. That would remove the whole OCE group.
   Koenig et al. (2024) re-ran QC for this callset to keep these populations. This audit uses their
   post-QC release: the {post.height} samples in their global PCA without outliers.
2. **80 vs 78 population labels.** The metadata field `hgdp_tgp_meta.Population` has 78 labels. The
   80-population label is `population`, which matches Koenig et al. and `post_qc_summary.tsv`. Both
   are kept in the sample table. Where they differ:

{_md_table(merged)}

3. **Unrelated set.** Unrelated = post-QC samples not in Koenig's `related_sample_ids` table. This gives
   {total['n_unrelated']} unrelated and {total['n_related']} related, which matches the sample counts of their
   released unrelated ({n_unrel_mt}) and related ({n_rel_mt}) PCA sets. The `n_unrelated` column of
   `post_qc_summary.tsv` sums to {total['koenig_n_unrelated']} instead. **For review:** it is unclear which
   relatedness pass that column comes from. This audit uses the released PCA sets.
4. **Excluded samples** ({excluded.height}):

{_md_table(excl)}

   `CHMI_CHMI3_WGS2` is the CHM13 synthetic-diploid control in the callset. The
   `not_in_koenig_post_qc_release` samples pass gnomAD's hard filters and are not PCA outliers,
   but are absent from Koenig's released post-QC set; their reason is not documented in the release.

## Populations with fewer than 10 unrelated samples

These stay in the descriptive 80-population analyses; they will contribute to inferential
results only through their seven-group category (Step 18).

{_md_table(small)}

## Population assignment review (Step 5)

- Every post-QC sample has exactly one original population label and one seven-group label.
- **Recently admixed populations (flagged, never excluded):**

{_md_table(admixed)}

- **PC discordance:** for each post-QC sample, the nearest seven-group centroid on Koenig's
  global PCs 1-{PCS_FOR_DISCORDANCE} is compared with its assigned group. This is a review aid only.

{_md_table(disc)}

- All 20 PCs are kept in `tables/sample_metadata.tsv` for the continuous-structure analysis (Step 28).
"""
    if review:
        text += "\n## Open review items\n\n" + "\n".join(f"- {r}" for r in review) + "\n"
    (settings.root / "reports" / "qc" / "sample_count_report.md").write_text(text)
