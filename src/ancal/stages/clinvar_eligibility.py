"""Stage 10: ClinVar coding-missense eligibility (item 12, protocol Step 29).

Classifies every ClinVar record in the scope regions for the clinical benchmark:

  eligible   an SNV whose REF equals GRCh38, missense in MANE Select, classified
             Pathogenic / Likely pathogenic (positive) or Benign / Likely benign (negative),
             with a valid AVI score
  excluded   with the first reason that applies (not an SNV, allele mismatch, not coding
             missense, VUS, conflicting, other classification, no AVI score)

VUS are never treated as benign. Review status and the submission fields the ClinVar VCF
carries are kept. This stage counts variants; it computes no classification metric.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import polars as pl

from .. import budget, config, consequence, status
from ..settings import Settings
from . import common

TAGS = ("ALLELEID", "CLNSIG", "CLNSIGCONF", "CLNREVSTAT", "CLNVC", "MC", "GENEINFO", "ORIGIN", "CLNDN", "CLNHGVS")
KEY = ["chrom", "pos", "ref", "alt"]


def read_clinvar(path: Path) -> pl.DataFrame:
    header = subprocess.run(["bcftools", "view", "-h", str(path)], capture_output=True, text=True, check=True).stdout
    present = [t for t in TAGS if re.search(rf"^##INFO=<ID={t},", header, re.M)]
    fmt = "%CHROM\t%POS\t%ID\t%REF\t%ALT" + "".join(f"\t%INFO/{t}" for t in present) + "\n"
    out = subprocess.run(["bcftools", "query", "-f", fmt, str(path)], capture_output=True, text=True, check=True).stdout
    cols = ["chrom", "pos", "clinvar_variation_id", "ref", "alt", *[t.lower() for t in present]]
    df = pl.read_csv(out.encode(), separator="\t", has_header=False, new_columns=cols, infer_schema_length=0,
                     null_values=["."], quote_char=None).with_columns(pl.col("pos").cast(pl.Int64))
    for t in TAGS:
        if t.lower() not in df.columns:
            df = df.with_columns(pl.lit(None, dtype=pl.String).alias(t.lower()))
    return df


def significance_class(clnsig: str | None, positive: list[str], negative: list[str]) -> str:
    """positive / negative / vus / conflicting / other / not_provided, from the primary CLNSIG term."""
    if not clnsig:
        return "not_provided"
    term = clnsig.split("|")[0].split(",")[0]
    if term in positive:
        return "positive"
    if term in negative:
        return "negative"
    if term == "Uncertain_significance":
        return "vus"
    if term.startswith("Conflicting"):
        return "conflicting"
    return "not_provided" if term == "not_provided" else "other"


def classify(cv: pl.DataFrame, cons: pl.DataFrame, cfg: config.AnalysisConfig, avi_scores: pl.DataFrame | None) -> pl.DataFrame:
    """One row per ClinVar record with eligibility and the first exclusion reason.

    `avi_scores` is Stage 8's table of AVI records in scope (None while Stage 8 is blocked).
    """
    pos_c, neg_c = cfg.clinvar.positive_class, cfg.clinvar.negative_class
    # REF is checked against the reference base at coding positions; Stage 4 checked every record.
    ref_at = cons.select("chrom", "pos", pl.col("ref").alias("_ref")).unique()
    cv = cv.join(ref_at, on=["chrom", "pos"], how="left").with_columns(
        pl.format("{}:{}:{}:{}", *KEY).alias("variant_id"),
        (pl.col("_ref").is_null() | (pl.col("ref").str.slice(0, 1) == pl.col("_ref"))).alias("ref_matches_grch38"),
        ((pl.col("ref").str.len_chars() == 1) & pl.col("alt").is_in(["A", "C", "G", "T"])).alias("is_snv"),
        pl.col("clnsig").map_elements(lambda s: significance_class(s, pos_c, neg_c), return_dtype=pl.String,
                                      skip_nulls=False).alias("significance_class"),
        pl.col("mc").fill_null("").str.contains("missense_variant").alias("clinvar_says_missense"),
    ).drop("_ref").join(consequence.per_variant(cons).select(*KEY, "gene_symbol", "ensembl_gene_id", "consequence", "protein_pos",
                                                              "ref_aa", "alt_aa"), on=KEY, how="left")
    if avi_scores is not None:
        valid = avi_scores.group_by(KEY).agg(pl.col("avi_score").drop_nulls().n_unique().alias("_n")).filter(pl.col("_n") == 1)
        cv = cv.join(valid, on=KEY, how="left")
        avi_state = pl.when(pl.col("_n").is_not_null()).then(pl.lit("valid_avi_score")).otherwise(pl.lit("no_valid_avi_score"))
    else:
        avi_state = pl.lit("avi_not_available")
    cv = cv.with_columns(avi_state.alias("avi_status"))
    reason = (
        pl.when(~pl.col("is_snv")).then(pl.lit("not_an_snv"))
        .when(~pl.col("ref_matches_grch38")).then(pl.lit("no_valid_allele_match"))
        .when(pl.col("consequence").is_null() | (pl.col("consequence") != "missense")).then(pl.lit("not_coding_missense"))
        .when(pl.col("significance_class") == "vus").then(pl.lit("uncertain_significance"))
        .when(pl.col("significance_class") == "conflicting").then(pl.lit("conflicting_classifications"))
        .when(~pl.col("significance_class").is_in(["positive", "negative"])).then(pl.lit("other_or_no_classification"))
        .otherwise(pl.lit(""))
    )
    return cv.with_columns(reason.alias("exclusion_reason")).with_columns(
        (pl.col("exclusion_reason") == "").alias("eligible_before_avi"),
        pl.when(pl.col("exclusion_reason") != "").then(pl.lit("excluded"))
        .when(pl.col("avi_status") == "valid_avi_score").then(pl.lit("eligible"))
        .when(pl.col("avi_status") == "no_valid_avi_score").then(pl.lit("excluded"))
        .otherwise(pl.lit("eligible_pending_avi")).alias("eligibility"),
    ).with_columns(
        pl.when((pl.col("exclusion_reason") == "") & (pl.col("avi_status") == "no_valid_avi_score"))
        .then(pl.lit("no_avi_score")).otherwise(pl.col("exclusion_reason")).alias("exclusion_reason"))


def by_population(cv: pl.DataFrame, long: pl.DataFrame) -> pl.DataFrame:
    """Eligible ClinVar variants observed in each population (unrelated and all-eligible sets)."""
    e = cv.filter(pl.col("eligibility") != "excluded").select("variant_id", "significance_class", "eligibility").unique(subset="variant_id")
    obs = long.filter(pl.col("observed")).select("variant_id", "sample_set", "level", "label", "population_group", "n_samples")
    levels = long.select("sample_set", "level", "label", "population_group", "n_samples").unique()
    agg = obs.join(e, on="variant_id").group_by("sample_set", "level", "label").agg(
        (pl.col("significance_class") == "positive").sum().alias("n_pathogenic_or_likely_pathogenic"),
        (pl.col("significance_class") == "negative").sum().alias("n_benign_or_likely_benign"),
        (pl.col("eligibility") == "eligible").sum().alias("n_with_valid_avi_score"),
        pl.len().alias("n_eligible_observed"))
    order = {"all": 0, "group": 1, "population": 2}
    return (levels.join(agg, on=["sample_set", "level", "label"], how="left").fill_null(0)
            .with_columns(pl.col("level").replace_strict(order, return_dtype=pl.Int32).alias("_o"))
            .sort(pl.col("sample_set") != "unrelated", "_o", "population_group", "label").drop("_o"))


def run(settings: Settings) -> status.StageStatus:
    budget.check(settings)
    ds = config.load_data_sources(settings.root)
    cfg = config.load_analysis(settings.root)
    p = common.Paths(settings)
    notes = "ClinVar coding-missense eligibility report (item 12)."
    vcfs = sorted(settings.data_dir("clinvar").glob(f"clinvar_*.{ds.scope}.vcf.gz"))
    if not vcfs or not p.consequences.exists():
        return status.write(settings, "stage10", [], notes=notes, blocked_reason=
                            "The ClinVar subset (Stage 4) or the consequence table (Stage 6) is missing for this scope.")
    # Use the release Stage 4 fetched for this run, not simply the newest file on disk.
    m = re.search(r"ClinVar release (\d{8})", (status.read(settings, "stage4") or {}).get("notes", ""))
    if m and (settings.data_dir("clinvar") / f"clinvar_{m.group(1)}.{ds.scope}.vcf.gz") in vcfs:
        vcfs = [settings.data_dir("clinvar") / f"clinvar_{m.group(1)}.{ds.scope}.vcf.gz"]
    release = vcfs[-1].name.split(".")[0].removeprefix("clinvar_")
    cons = pl.read_parquet(p.consequences)
    cv = read_clinvar(vcfs[-1])
    avi_scores = pl.read_parquet(p.avi_scores) if p.avi_scores.exists() else None
    long = pl.read_parquet(p.long) if p.long.exists() else None
    out = classify(cv, cons, cfg, avi_scores)
    out.write_parquet(p.clinvar)

    counts = out.group_by("eligibility", "exclusion_reason", "significance_class").len().rename({"len": "n_records"}).sort(
        "eligibility", "exclusion_reason", "significance_class")
    tables = [common.write_tsv(counts, p.tables / "clinvar_eligibility_counts.tsv")]
    elig = out.filter(pl.col("eligibility") != "excluded")
    keep = ["variant_id", "clinvar_variation_id", "alleleid", "gene_symbol", "ensembl_gene_id", "protein_pos", "ref_aa", "alt_aa",
            "clnsig", "significance_class", "clnrevstat", "clnsigconf", "origin", "clndn", "eligibility", "avi_status"]
    tables.append(common.write_tsv(elig.select(keep).sort("variant_id"), p.tables / "clinvar_eligible_variants.tsv"))
    pop = by_population(out, long) if long is not None else None
    if pop is not None:
        tables.append(common.write_tsv(pop, p.tables / "clinvar_eligibility_by_population.tsv"))
    tables.append(common.write_tsv(out.filter(pl.col("eligibility") == "excluded").group_by("exclusion_reason").len()
                                   .rename({"len": "n_records"}).sort("n_records", descending=True),
                                   p.exclusions / "clinvar_exclusions.tsv"))
    common.record(settings, "stage10", "clinvar_eligibility", f"clinvar/{vcfs[-1].name}", p.clinvar, *tables)

    snv_cds = out.filter(pl.col("is_snv") & pl.col("consequence").is_not_null())
    agree = snv_cds.filter(pl.col("clinvar_says_missense") == (pl.col("consequence") == "missense")).height
    dup = elig.height - elig["variant_id"].n_unique()
    checks = [
        status.Check("clinvar_records", out.height > 0, f"release {release}: {out.height:,} records in scope"),
        status.Check("every_record_classified", out.filter((pl.col("eligibility") == "excluded") & (pl.col("exclusion_reason") == "")).height == 0,
                     f"{elig.height:,} eligible or eligible pending AVI; {out.height - elig.height:,} excluded"),
        status.Check("vus_not_benign", out.filter((pl.col("significance_class") == "vus") & (pl.col("eligibility") != "excluded")).height == 0,
                     "no VUS in the eligible set"),
        status.Check("missense_agrees_with_clinvar", snv_cds.height == 0 or agree / snv_cds.height >= 0.95,
                     f"MANE Select consequence agrees with ClinVar's MC for {agree:,} of {snv_cds.height:,} coding SNV records"),
        status.Check("unique_eligible_variants", dup == 0, f"{dup} eligible variants with more than one ClinVar record"),
    ]
    u = budget.usage(settings)
    checks.append(status.Check("disk_budget", u.ok, f"{u.used_mb:.1f} MB of {settings.local_budget_gb} GB"))
    write_report(p, ds.scope, release, out, counts, pop, avi_scores is not None, agree, snv_cds.height)
    return status.write(settings, "stage10", checks, notes=notes + f" ClinVar release {release}.")


def write_report(p, scope, release, out, counts, pop, have_avi, agree, n_cds):
    elig = out.filter(pl.col("eligibility") != "excluded")
    funnel = [
        ("ClinVar records in the scope regions", out.height),
        ("SNVs", out.filter(pl.col("is_snv")).height),
        ("... with a REF equal to GRCh38", out.filter(pl.col("is_snv") & pl.col("ref_matches_grch38")).height),
        ("... missense in MANE Select", out.filter(pl.col("is_snv") & pl.col("ref_matches_grch38") & (pl.col("consequence") == "missense")).height),
        ("... classified P/LP or B/LB (eligible before AVI)", out.filter(pl.col("eligible_before_avi")).height),
        ("... with a valid AVI score (eligible)", elig.filter(pl.col("eligibility") == "eligible").height if have_avi else None),
    ]
    cls = elig.group_by("significance_class").len().sort("significance_class")
    rev = elig.group_by("clnrevstat").len().sort("len", descending=True)
    gene = elig.group_by("gene_symbol").agg(
        (pl.col("significance_class") == "positive").sum().alias("P/LP"), (pl.col("significance_class") == "negative").sum().alias("B/LB")
    ).sort("gene_symbol").head(100)
    excl = out.filter(pl.col("eligibility") == "excluded").group_by("exclusion_reason").len().sort("len", descending=True)
    caveats = [] if have_avi else ["AVI scores are not available (Stage 8 blocked), so `eligible_pending_avi` is an upper bound"]
    poptext = "Population variants (Stage 7) are not available, so observed counts were not computed."
    if pop is not None:
        poptext = common.md_table(pop.filter((pl.col("sample_set") == "unrelated") & (pl.col("level") != "population")).drop("sample_set", "population_group")) + \
            "\n\nAll 80 populations and the all-eligible sample set: `tables/clinvar_eligibility_by_population.tsv`."
    text = f"""# ClinVar coding-missense eligibility report (Stage 10, protocol Step 29)

Generated by `ancal stage10`. Do not edit. Scope: `{scope}`. ClinVar release: **{release}**.

{common.provisional_banner(scope, caveats)}This report counts variants. It contains no AUROC or other classification result (Step 30 is
part of the second deliverable).

## Eligibility

{common.md_table(pl.DataFrame(funnel, orient='row', schema=['Step', 'Records']))}

Positive class = Pathogenic, Likely pathogenic. Negative class = Benign, Likely benign.
Variants of uncertain significance and conflicting classifications are excluded; a VUS is never
counted as benign. "Missense" is the pipeline's MANE Select consequence; it agrees with ClinVar's
own molecular consequence for {agree:,} of {n_cds:,} coding SNV records.

### Eligible variants by class

{common.md_table(cls, ['Class', 'Variants'])}

### Eligible variants by ClinVar review status

{common.md_table(rev, ['Review status', 'Variants'])}

### Eligible variants by gene{' (first 100)' if gene.height == 100 else ''}

{common.md_table(gene)}

### Exclusions (first reason that applies)

{common.md_table(excl, ['Reason', 'Records'])}

## Eligible variants observed in each population

A ClinVar variant counts for a population when its alternate allele count there is at least 1
after QC. Most ClinVar variants are not observed in any of the 4,094 genomes.

{poptext}

## Fields kept

Review status (`CLNREVSTAT`), conflicting-classification detail (`CLNSIGCONF`), allele origin
(`ORIGIN`), condition (`CLNDN`), variation ID and allele ID. The ClinVar VCF does not carry
per-submitter records; if submitter-level detail is needed, it has to come from ClinVar's
`submission_summary.txt.gz` (listed under unresolved issues).

## Outputs

- `tables/clinvar_eligible_variants.tsv`: the eligible variants.
- `tables/clinvar_eligibility_counts.tsv`: every record by eligibility, reason and class.
- `data_store/clinvar/{p.clinvar.name}`: every ClinVar record in scope with its classification.
"""
    (p.qc / "clinvar_eligibility_report.md").write_text(text)
