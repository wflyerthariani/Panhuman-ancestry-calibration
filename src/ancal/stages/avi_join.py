"""Stage 8: AVI join and coverage (item 8, protocol Steps 9-11).

Blocked until the `avi:` source in config/data_sources.yaml is filled in (L-02).

1. Load AVI records inside the scope regions (ancal.sources.avi) and normalize them to the
   pipeline's variant identifier (chrN, upper-case alleles).
2. Join to the population variants on chrom, pos, ref, alt only (never rsID), and classify
   every attempted join into one of the six Step 10 categories.
3. Report coverage and no-evidence rate by population, group, allele-frequency bin, coding
   status, consequence and gene.

This stage reports score availability only. It computes no performance metric. The direction
check compares AVI with AlphaMissense (one of its own components) and with nonsense versus
synonymous variants; it confirms the sign of the score and says nothing about accuracy.
"""

from __future__ import annotations

import polars as pl

from .. import budget, config, status
from ..settings import Settings
from ..sources import avi
from . import common

VALID = ("exact_allele_match", "corrected_representation_match")
CLASSES = (*VALID, "allele_swapped_match", "multiple_possible_matches", "no_avi_score", "unresolved_allele_conflict")
KEY = ["chrom", "pos", "ref", "alt"]
MIN_REF_MATCH = 0.99


def normalize(raw: pl.DataFrame) -> pl.DataFrame:
    """AVI records with pipeline-style chromosome names and alleles; flags rows that had to be changed."""
    return raw.with_columns(
        pl.when(pl.col("chrom_raw").str.starts_with("chr")).then(pl.col("chrom_raw"))
        .otherwise(pl.concat_str(pl.lit("chr"), pl.col("chrom_raw"))).alias("chrom"),
        pl.col("ref_raw").str.to_uppercase().alias("ref"), pl.col("alt_raw").str.to_uppercase().alias("alt"),
    ).with_columns(
        (~pl.col("chrom_raw").str.starts_with("chr") | (pl.col("ref_raw") != pl.col("ref"))
         | (pl.col("alt_raw") != pl.col("alt"))).alias("representation_corrected"))


def classify_joins(sites: pl.DataFrame, scores: pl.DataFrame) -> pl.DataFrame:
    """One row per population variant: join class (Step 10) and, for valid joins, the AVI score.

    `sites` needs chrom, pos, ref, alt, variant_id, ref_matches_grch38. `scores` is normalize() output.
    """
    per_key = scores.group_by(KEY).agg(
        pl.len().alias("n_avi_records"), pl.col("avi_score").drop_nulls().n_unique().alias("n_distinct_scores"),
        pl.col("avi_score").drop_nulls().first(), pl.col("representation_corrected").any())
    swapped = scores.select("chrom", "pos", pl.col("alt").alias("ref"), pl.col("ref").alias("alt")).unique().with_columns(
        pl.lit(True).alias("_swapped"))
    # Reference alleles AVI uses at each position, to tell "not scored" from "scored against another reference".
    at_pos = scores.group_by("chrom", "pos").agg(pl.col("ref").unique().alias("_avi_refs"))
    j = (sites.select("variant_id", *KEY, "ref_matches_grch38").join(per_key, on=KEY, how="left")
         .join(swapped, on=KEY, how="left").join(at_pos, on=["chrom", "pos"], how="left"))
    matched = pl.col("n_avi_records").is_not_null()
    return j.with_columns(
        pl.when(~pl.col("ref_matches_grch38")).then(pl.lit("unresolved_allele_conflict"))
        .when(matched & (pl.col("n_distinct_scores") > 1)).then(pl.lit("multiple_possible_matches"))
        .when(matched & (pl.col("n_distinct_scores") == 1) & pl.col("representation_corrected"))
        .then(pl.lit("corrected_representation_match"))
        .when(matched & (pl.col("n_distinct_scores") == 1)).then(pl.lit("exact_allele_match"))
        .when(matched).then(pl.lit("no_avi_score"))
        .when(pl.col("_swapped").fill_null(False)).then(pl.lit("allele_swapped_match"))
        .when((pl.col("ref").str.len_chars() == 1) & pl.col("_avi_refs").is_not_null()
              & ~pl.col("_avi_refs").list.contains(pl.col("ref")))
        .then(pl.lit("unresolved_allele_conflict"))
        .otherwise(pl.lit("no_avi_score")).alias("avi_join_class"),
    ).with_columns(
        pl.col("avi_join_class").is_in(VALID).alias("avi_valid"),
        pl.when(pl.col("avi_join_class").is_in(VALID)).then(pl.col("avi_score")).alias("avi_score"),
    ).select("variant_id", "avi_join_class", "avi_valid", "avi_score", pl.col("n_avi_records").fill_null(0))


def coverage(long: pl.DataFrame, sites: pl.DataFrame, join: pl.DataFrame) -> pl.DataFrame:
    """AVI coverage and no-evidence rate (Step 11) per population level and stratum.

    Eligible observed variants = primary SNVs (PASS, in scope, REF confirmed) with an alternate
    allele count of at least 1 in that population.
    """
    s = sites.filter(pl.col("primary_snv")).join(join.select("variant_id", "avi_valid"), on="variant_id").select(
        "variant_id", "avi_valid", "af_bin",
        pl.when(pl.col("coding")).then(pl.lit("coding")).otherwise(pl.lit("noncoding")).alias("coding_status"),
        pl.col("consequence").fill_null("noncoding"), pl.col("gene_symbol").fill_null("(no MANE Select CDS)").alias("gene"))
    obs = long.filter(pl.col("observed")).select("variant_id", "sample_set", "level", "label", "population_group").join(s, on="variant_id")
    out = []
    for name, col in (("overall", None), ("af_bin", "af_bin"), ("coding_status", "coding_status"),
                      ("consequence", "consequence"), ("gene", "gene")):
        by = ["sample_set", "level", "label", "population_group"] + ([col] if col else [])
        g = obs.group_by(by).agg(pl.len().alias("n_eligible_observed"), pl.col("avi_valid").sum().alias("n_with_valid_avi_score"))
        g = g.with_columns(pl.lit(name).alias("stratifier"), (pl.col(col) if col else pl.lit("all")).alias("stratum"))
        out.append(g.select("sample_set", "level", "label", "population_group", "stratifier", "stratum",
                            "n_eligible_observed", "n_with_valid_avi_score"))
    order = {"all": 0, "group": 1, "population": 2}
    return (pl.concat(out).with_columns(
        (pl.col("n_with_valid_avi_score") / pl.col("n_eligible_observed")).alias("avi_coverage"),
        (1 - pl.col("n_with_valid_avi_score") / pl.col("n_eligible_observed")).alias("no_evidence_rate"),
        pl.col("level").replace_strict(order, return_dtype=pl.Int32).alias("_o"))
        .sort(pl.col("sample_set") != "unrelated", "stratifier", "_o", "population_group", "label", "stratum").drop("_o"))


def direction_evidence(scores: pl.DataFrame, cons: pl.DataFrame, am_path) -> dict:
    """Sign checks for the raw AVI score. Not a performance result."""
    one = scores.group_by(KEY).agg(pl.col("avi_score").drop_nulls().mean()).drop_nulls("avi_score")
    c = one.join(cons.select(*KEY, "consequence").unique(subset=KEY), on=KEY)
    med = dict(c.group_by("consequence").agg(pl.col("avi_score").median()).iter_rows())
    n = dict(c.group_by("consequence").len().iter_rows())
    am = pl.read_csv(am_path, separator="\t", comment_prefix="#", has_header=False,
                     new_columns=[*KEY, "genome", "uniprot_id", "transcript_id", "protein_variant", "am_pathogenicity", "am_class"],
                     schema_overrides={"pos": pl.Int64, "am_pathogenicity": pl.Float64})
    j = one.join(am.select(*KEY, "am_pathogenicity").unique(subset=KEY), on=KEY)
    rho = common.spearman(j, "avi_score", "am_pathogenicity")
    votes = []
    if rho is not None:
        votes.append("higher_is_more_damaging" if rho > 0 else "lower_is_more_damaging")
    if n.get("nonsense", 0) >= 10 and n.get("synonymous", 0) >= 10 and med["nonsense"] != med["synonymous"]:
        votes.append("higher_is_more_damaging" if med["nonsense"] > med["synonymous"] else "lower_is_more_damaging")
    inferred = votes[0] if votes and len(set(votes)) == 1 else None
    return {"rho_alphamissense": rho, "n_alphamissense": j.height, "median": med, "n": n, "inferred": inferred, "votes": votes}


def blocked(settings: Settings, p: common.Paths, reason: str, notes: str) -> status.StageStatus:
    for stale in (p.avi_scores, p.avi_join, p.tables / "avi_coverage.tsv", p.tables / "avi_join_classes.tsv",
                  p.exclusions / "avi_join_exclusions.tsv"):
        stale.unlink(missing_ok=True)
    (p.qc / "avi_join_report.md").write_text(
        "# AVI join and coverage report (Stage 8, protocol Steps 9-11)\n\nGenerated by `ancal stage8`. Do not edit. "
        f"Scope: `{settings.scope}`.\n\n**BLOCKED.** {reason}\n\nNo AVI score has been joined. Stages 9-11 report every "
        "AVI-dependent count as not available. See `SCALING_UP.md` section 3.\n")
    return status.write(settings, "stage8", [], notes=notes, blocked_reason=reason)


def run(settings: Settings) -> status.StageStatus:
    budget.check(settings)
    ds = config.load_data_sources(settings.root)
    p = common.Paths(settings)
    src = ds.sources["avi"]
    notes = "AVI join and coverage report (item 8)."
    if not src.uri:
        return blocked(settings, p, "The AVI source is not configured: `sources.avi.uri` in config/data_sources.yaml is "
                       "null. Access to the AlphaGenome Atlas AVI scores has not been arranged (L-02).", notes)
    if not p.sites.exists():
        return blocked(settings, p, "The population-variant table is missing; run `make stage7` first.", notes)
    bed = common.read_bed(p.bed)
    try:
        raw, info = avi.load(settings, src, bed, common.bed_positions(bed), ds.scope)
    except avi.AviConfigError as e:
        return blocked(settings, p, f"The AVI source is configured but cannot be read: {e}", notes)

    scores = normalize(raw)
    sites, long, cons = pl.read_parquet(p.sites), pl.read_parquet(p.long), pl.read_parquet(p.consequences)
    ev = direction_evidence(scores, cons, p.alphamissense)
    direction = src.direction or ev["inferred"]
    sign = -1.0 if direction == "lower_is_more_damaging" else 1.0
    join = classify_joins(sites, scores).with_columns((pl.col("avi_score") * sign).alias("avi_score_harmonized"))
    scores = scores.with_columns((pl.col("avi_score") * sign).alias("avi_score_harmonized"))
    cov = coverage(long, sites, join)

    for path in (p.avi_scores, p.avi_join):
        path.parent.mkdir(parents=True, exist_ok=True)
    scores.write_parquet(p.avi_scores)
    join.write_parquet(p.avi_join)
    t_cov = common.write_tsv(cov, p.tables / "avi_coverage.tsv")
    classes = (pl.DataFrame({"avi_join_class": list(CLASSES)}).join(
        join.join(sites.select("variant_id", "primary_snv", "variant_type"), on="variant_id")
        .group_by("avi_join_class").agg(pl.len().alias("n_variants"), pl.col("primary_snv").sum().alias("n_primary_snvs")),
        on="avi_join_class", how="left").fill_null(0))
    t_cls = common.write_tsv(classes, p.tables / "avi_join_classes.tsv")
    not_valid = join.filter(~pl.col("avi_valid")).join(sites.select("variant_id", "primary_snv"), on="variant_id").filter(
        pl.col("primary_snv")).group_by("avi_join_class").len().rename({"len": "n_primary_snvs_excluded"}).sort("avi_join_class")
    t_ex = common.write_tsv(not_valid, p.exclusions / "avi_join_exclusions.tsv")
    common.record(settings, "stage8", "avi_join", src.uri, p.avi_scores, p.avi_join, t_cov, t_cls, t_ex)

    snv = scores.filter((pl.col("ref").str.len_chars() == 1) & (pl.col("alt").str.len_chars() == 1))
    ref_ok = snv.join(cons.select("chrom", "pos", pl.col("ref").alias("_ref")).unique(), on=["chrom", "pos"])
    ref_rate = (ref_ok["ref"] == ref_ok["_ref"]).mean() if ref_ok.height else None
    multi = scores.group_by(KEY).agg(pl.len().alias("n"), pl.col("avi_score").n_unique().alias("d"))
    facts = {
        "n_records": scores.height, "n_variants": multi.height, "min": scores["avi_score"].min(), "max": scores["avi_score"].max(),
        "n_null": scores["avi_score"].null_count(), "multi_records": multi.filter(pl.col("n") > 1).height,
        "multi_scores": multi.filter(pl.col("d") > 1).height,
        "transcript_specific": scores["transcript"].drop_nulls().len() > 0 and multi.filter(pl.col("d") > 1).height > 0,
        "ref_rate": ref_rate, "n_ref_checked": ref_ok.height,
    }
    n_valid = join.filter(pl.col("avi_valid")).height
    agree = src.direction is None or not ev["votes"] or all(v == src.direction for v in ev["votes"])
    checks = [
        status.Check("avi_records_in_scope", scores.height > 0, f"{scores.height:,} records for {multi.height:,} variants "
                     f"({info['access']})"),
        status.Check("avi_reference_build", ref_rate is not None and ref_rate >= MIN_REF_MATCH,
                     "REF equals GRCh38 for " + (f"{ref_rate:.2%} of {ref_ok.height:,} coding SNV records" if ref_rate is not None
                                                 else "no checkable records")),
        status.Check("score_direction", direction is not None and agree,
                     f"configured: {src.direction or 'not set'}; evidence: {ev['votes'] or 'none'}"
                     + ("" if src.direction else f"; INFERRED {direction}, confirm and set avi.direction")),
        status.Check("every_join_classified", join.height == sites.height and join["avi_join_class"].is_in(CLASSES).all(),
                     f"{join.height:,} population variants; {n_valid:,} with a valid AVI score"),
        status.Check("coverage_in_range", cov.height > 0 and cov["avi_coverage"].is_between(0, 1).all(), f"{cov.height:,} coverage rows"),
    ]
    u = budget.usage(settings)
    checks.append(status.Check("disk_budget", u.ok, f"{u.used_mb:.1f} MB of {settings.local_budget_gb} GB"))
    write_report(settings, p, ds.scope, src, info, facts, ev, direction, classes, cov)
    return status.write(settings, "stage8", checks, notes=notes)


def write_report(settings, p, scope, src, info, facts, ev, direction, classes, cov):
    u = cov.filter(pl.col("sample_set") == "unrelated").drop("sample_set")

    def strat(name, levels=("all", "group")):
        return common.md_table(u.filter((pl.col("stratifier") == name) & pl.col("level").is_in(levels)).drop("stratifier", "population_group"))

    med = "; ".join(f"{k} {v:.4g} (n={ev['n'].get(k, 0):,})" for k, v in sorted(ev["median"].items()))
    text = f"""# AVI join and coverage report (Stage 8, protocol Steps 9-11)

Generated by `ancal stage8`. Do not edit. Scope: `{scope}`.

{common.provisional_banner(scope)}**This report is about score availability. It contains no accuracy or performance result.**
A high coverage or a high damaging-call rate is not evidence of accuracy (Step 11).

## The AVI source (Step 9)

| | |
|---|---|
| Source | `{src.uri}` |
| Version | {src.version or 'not stated in config'} |
| File last modified | {info.get('last_modified') or 'not available'} |
| Format and access | {src.format}; {info['access']} |
| Reference build | GRCh38 assumed; REF equals GRCh38 for {'—' if facts['ref_rate'] is None else format(facts['ref_rate'], '.2%')} of {facts['n_ref_checked']:,} coding SNV records checked |
| Score field | `{src.score_field}` |
| Score range in scope | {common.fmt(facts['min'])} to {common.fmt(facts['max'])} |
| Records / distinct variants in scope | {facts['n_records']:,} / {facts['n_variants']:,} |
| Records with a null score | {facts['n_null']:,} |
| Variants with more than one record | {facts['multi_records']:,} |
| Variants with more than one distinct score | {facts['multi_scores']:,} |
| Transcript-specific | {'yes' if facts['transcript_specific'] else 'no evidence of it in scope'} |
| Damaging threshold | {src.damaging_threshold if src.damaging_threshold is not None else 'none defined'} |
| Columns used | {info['columns']} |

## Score direction

- **Configured:** {src.direction or 'not set'}. **Used:** {direction or 'undetermined'}.
- **Check 1:** Spearman correlation of the raw AVI score with AlphaMissense pathogenicity over
  {ev['n_alphamissense']:,} missense SNVs in scope: {common.fmt(ev['rho_alphamissense'])}. AlphaMissense is a component of AVI,
  so this only confirms the sign.
- **Check 2:** median raw AVI score by consequence: {med or 'not available'}.
- `avi_score_harmonized` = the raw score{' multiplied by -1' if direction == 'lower_is_more_damaging' else ''}, so that higher
  means more predicted damage. `avi_score` is kept unchanged.

## Join classification (Step 10)

Joins use chromosome, position, reference allele and alternate allele. rsIDs are never used.
Only the first two classes count as a valid score.

{common.md_table(classes)}

- `corrected_representation_match`: matched after adding the `chr` prefix or upper-casing alleles.
- `allele_swapped_match`: AVI has the variant with REF and ALT exchanged. Not used.
- `multiple_possible_matches`: several AVI records with different scores for one variant. Not used.
- `unresolved_allele_conflict`: the population variant's REF, or AVI's REF at that position, differs from GRCh38.
- `no_avi_score`: no AVI record for the exact alleles, or a null score.

## Coverage and no-evidence rate (Step 11)

AVI coverage = eligible observed variants with a valid AVI score / all eligible observed variants.
No-evidence rate = 1 - coverage. Unrelated sample set; all strata, all 80 populations and the
all-eligible set are in `tables/avi_coverage.tsv`.

### Overall and by group

{strat('overall')}

### By global allele-frequency bin

{strat('af_bin')}

### By coding status

{strat('coding_status')}

### By consequence

{strat('consequence')}

### By gene (all samples combined)

{strat('gene', levels=('all',))}

### By population

{strat('overall', levels=('population',))}
"""
    (p.qc / "avi_join_report.md").write_text(text)
