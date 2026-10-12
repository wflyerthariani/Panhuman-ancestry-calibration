"""Stage 9: proposed DMS gene universe (item 10, protocol Steps 12 and 35).

Applies the Step 12 inclusion criteria to every candidate gene in scope:

  quantitative variant-level measurements   a candidate assay exists (Stage 5 pre-screen)
  effect can be consistently oriented       an included assay has a score direction (review, Stage 6)
  acceptable assay quality                  the reviewer included the assay
  identifiable on GRCh38                    Stage 6 placed its variants on GRCh38
  at least 20 eligible variants             counted as set in dms_assays.min_variants_basis
  not excluded for known leakage            recorded per gene; nothing is known yet (Step 35)

The result is a PROPOSAL. The universe is frozen only after the team approves it, by writing
it into config/analysis.yaml (`dms_gene_universe`) and running `ancal freeze-config`.
"""

from __future__ import annotations

import polars as pl

from .. import budget, config, status
from ..settings import Settings
from . import common

LEAKAGE_UNKNOWN = "overlap_unknown"


def build(inventory: pl.DataFrame, assays: pl.DataFrame, dms: pl.DataFrame, sites: pl.DataFrame | None,
          avi_join: pl.DataFrame | None, scope_genes: list[str], min_variants: int, basis: str) -> pl.DataFrame:
    cand = inventory.filter(pl.col("prescreen_status") == "candidate")
    genes = cand.group_by("gene_symbol", "ensembl_gene_id").agg(pl.len().alias("n_candidate_assays")).sort("gene_symbol")
    a = assays.group_by("gene_symbol").agg(
        (pl.col("decision") == "include").sum().alias("n_assays_passing_review"),
        (pl.col("harmonization_status") == "harmonized").sum().alias("n_assays_harmonized"),
        pl.col("score_set_urn").filter(pl.col("harmonization_status") == "harmonized").sort().str.join(";").alias("harmonized_assays"),
        ((pl.col("harmonization_status") == "harmonized") & (pl.col("direction_source") != "reviewer")).sum().alias("n_assays_unreviewed"),
        (pl.col("decision") == "pending").sum().alias("n_assays_awaiting_review"),
    )
    miss = dms.filter(pl.col("consequence") == "missense").select("gene_symbol", "variant_id").unique()
    counts = miss.group_by("gene_symbol").agg(pl.len().alias("n_harmonized_missense_snvs"))
    if sites is not None:
        obs = sites.filter(pl.col("primary_snv") & (pl.col("ac_global") >= 1)).select("variant_id")
        m_obs = miss.join(obs, on="variant_id")
        counts = counts.join(m_obs.group_by("gene_symbol").agg(pl.len().alias("n_population_observed_missense_snvs")),
                             on="gene_symbol", how="left").with_columns(pl.col("n_population_observed_missense_snvs").fill_null(0))
        if avi_join is not None:
            m_avi = m_obs.join(avi_join.filter(pl.col("avi_valid")).select("variant_id"), on="variant_id")
            counts = counts.join(m_avi.group_by("gene_symbol").agg(pl.len().alias("n_observed_with_valid_avi_score")),
                                 on="gene_symbol", how="left").with_columns(pl.col("n_observed_with_valid_avi_score").fill_null(0))
    for c in ("n_population_observed_missense_snvs", "n_observed_with_valid_avi_score"):
        if c not in counts.columns:
            counts = counts.with_columns(pl.lit(None, dtype=pl.UInt32).alias(c))

    g = genes.join(a, on="gene_symbol", how="left").join(counts, on="gene_symbol", how="left").with_columns(
        pl.col("gene_symbol").is_in(scope_genes).alias("in_scope"),
        pl.col("n_assays_passing_review", "n_assays_harmonized", "n_assays_unreviewed", "n_assays_awaiting_review",
               "n_harmonized_missense_snvs").fill_null(0),
        pl.col("harmonized_assays").fill_null(""))
    n_eligible = pl.col("n_harmonized_missense_snvs" if basis == "harmonized_missense_snvs" else "n_population_observed_missense_snvs")
    g = g.with_columns(
        n_eligible.alias("n_eligible_variants"),
        (pl.col("n_assays_passing_review") > 0).alias("orientable_and_quality_ok"),
        (pl.col("n_assays_harmonized") > 0).alias("identifiable_on_grch38"),
        (n_eligible.fill_null(0) >= min_variants).alias("meets_min_variants"),
        pl.lit(LEAKAGE_UNKNOWN).alias("leakage_status"),
    ).with_columns(
        pl.when(~pl.col("in_scope")).then(pl.lit("not_assessed_at_this_scope"))
        .when(pl.col("n_assays_awaiting_review") == pl.col("n_candidate_assays")).then(pl.lit("awaiting_review"))
        .when(~pl.col("orientable_and_quality_ok")).then(pl.lit("no_assay_with_a_direction_and_acceptable_quality"))
        .when(~pl.col("identifiable_on_grch38")).then(pl.lit("no_assay_identifiable_on_grch38"))
        .when(n_eligible.is_null()).then(pl.lit("population_variants_not_available"))
        .when(~pl.col("meets_min_variants")).then(pl.lit(f"fewer_than_{min_variants}_eligible_variants"))
        .otherwise(pl.lit("")).alias("exclusion_reason"),
    ).with_columns(
        pl.when(~pl.col("in_scope")).then(pl.lit("not_assessed"))
        .when(pl.col("exclusion_reason") == "").then(pl.lit("proposed_include"))
        .otherwise(pl.lit("proposed_exclude")).alias("inclusion_status"),
        pl.when(pl.col("n_assays_unreviewed") > 0).then(pl.lit("provisional_unreviewed_assays")).otherwise(pl.lit("")).alias("caveat"),
    )
    return g.select(
        "gene_symbol", "ensembl_gene_id", "in_scope", "inclusion_status", "exclusion_reason", "caveat", "n_candidate_assays",
        "n_assays_passing_review", "n_assays_harmonized", "n_assays_unreviewed", "harmonized_assays",
        "n_harmonized_missense_snvs", "n_population_observed_missense_snvs", "n_observed_with_valid_avi_score",
        "n_eligible_variants", "orientable_and_quality_ok", "identifiable_on_grch38", "meets_min_variants", "leakage_status")


def run(settings: Settings) -> status.StageStatus:
    budget.check(settings)
    ds = config.load_data_sources(settings.root)
    cfg = config.load_analysis(settings.root)
    p = common.Paths(settings)
    notes = "Proposed DMS gene universe (item 10). Not frozen: needs team approval."
    if not p.dms.exists():
        return status.write(settings, "stage9", [], notes=notes,
                            blocked_reason="The harmonized DMS table is missing; run `make stage6` first.")
    inventory = pl.read_csv(p.tables / "dms_candidate_inventory.tsv", separator="\t", infer_schema_length=0)
    assays = pl.read_csv(p.tables / "dms_assay_harmonization.tsv", separator="\t", infer_schema_length=0)
    dms = pl.read_parquet(p.dms)
    sites = pl.read_parquet(p.sites) if p.sites.exists() else None
    avi_join = pl.read_parquet(p.avi_join) if p.avi_join.exists() else None
    scope_genes = common.scope_genes(common.read_bed(p.bed))
    min_n, basis = cfg.dms_assays.min_variants_per_gene, cfg.dms_assays.min_variants_basis.value
    uni = build(inventory, assays, dms, sites, avi_join, scope_genes, min_n, basis)
    t = common.write_tsv(uni, p.tables / "proposed_dms_gene_universe.tsv")
    common.record(settings, "stage9", "gene_universe", "tables/dms_candidate_inventory.tsv, tables/dms_assay_harmonization.tsv, "
                  f"dms/{p.dms.name}", t)

    in_scope = uni.filter(pl.col("in_scope"))
    inc = in_scope.filter(pl.col("inclusion_status") == "proposed_include")
    missing = sorted(set(scope_genes) & set(inventory.filter(pl.col("prescreen_status") == "candidate")["gene_symbol"]) - set(in_scope["gene_symbol"]))
    checks = [
        status.Check("criteria_basis", basis in ("harmonized_missense_snvs", "population_observed_missense_snvs"),
                     f"at least {min_n} {basis}"),
        status.Check("every_scope_gene_assessed", not missing and in_scope.height > 0,
                     f"{in_scope.height} candidate genes in scope" + (f"; missing {missing}" if missing else "")),
        status.Check("every_exclusion_has_a_reason",
                     in_scope.filter((pl.col("inclusion_status") == "proposed_exclude") & (pl.col("exclusion_reason") == "")).height == 0,
                     f"{inc.height} proposed for inclusion, {in_scope.height - inc.height} proposed for exclusion"),
    ]
    u = budget.usage(settings)
    checks.append(status.Check("disk_budget", u.ok, f"{u.used_mb:.1f} MB of {settings.local_budget_gb} GB"))
    write_report(p, ds.scope, uni, min_n, basis, sites is not None, avi_join is not None, ds.review_gates.dms_review)
    return status.write(settings, "stage9", checks, notes=notes + (" PROVISIONAL at smoke scope." if ds.scope == "smoke" else ""))


def write_report(p, scope, uni, min_n, basis, have_sites, have_avi, mode):
    in_scope = uni.filter(pl.col("in_scope"))
    caveats = []
    if in_scope.filter(pl.col("caveat") != "").height:
        caveats.append(f"{in_scope.filter(pl.col('caveat') != '').height} gene(s) rest on assays without a reviewer decision")
    if not have_sites:
        caveats.append("population variants (Stage 7) are not available, so observed counts are missing")
    if not have_avi:
        caveats.append("AVI scores (Stage 8) are not available")
    by_reason = in_scope.group_by("inclusion_status", "exclusion_reason").len().sort("inclusion_status", "exclusion_reason")
    shown = in_scope.select("gene_symbol", "ensembl_gene_id", "inclusion_status", "exclusion_reason", "caveat",
                            "n_candidate_assays", "n_assays_harmonized", "n_harmonized_missense_snvs",
                            "n_population_observed_missense_snvs", "n_observed_with_valid_avi_score")
    text = f"""# Gene-universe report (Stage 9, protocol Step 12)

Generated by `ancal stage9`. Do not edit. Scope: `{scope}`. Review gate: `{mode}`.

{common.provisional_banner(scope, caveats)}**This is a proposal for team review. It is not frozen.** After approval, the approved Ensembl gene
IDs are written to `dms_gene_universe` in `config/analysis.yaml` and the config is frozen with
`ancal freeze-config`, before any population comparison is generated.

## Criteria applied (Step 12)

| Criterion | How it is tested |
|---|---|
| Variant-level quantitative measurements | the gene has a candidate assay after the Stage 5 pre-screen |
| Effect can be consistently oriented | at least one assay is included with a score direction |
| Assay quality is acceptable | the reviewer's `include` decision (or, provisionally, the pipeline's proposal) |
| Variants identifiable on GRCh38 | Stage 6 placed the assay's changes on MANE Select |
| At least {min_n} eligible variants after harmonization | `{basis}` (both counts are shown below) |
| Not excluded for known data leakage | no leakage information is available yet; every gene is `{LEAKAGE_UNKNOWN}` (Step 35) |

## Result

| | Genes |
|---|---|
| Candidate genes in MaveDB (all scopes) | {uni.height} |
| Assessed at this scope | {in_scope.height} |
| Proposed for inclusion | {in_scope.filter(pl.col('inclusion_status') == 'proposed_include').height} |
| Proposed for exclusion | {in_scope.filter(pl.col('inclusion_status') == 'proposed_exclude').height} |
| Not assessed at this scope | {uni.height - in_scope.height} |

{common.md_table(by_reason, ['Status', 'Reason', 'Genes'])}

## Genes assessed{' (first 200)' if shown.height > 200 else ''}

{common.md_table(shown.head(200))}

Every candidate gene, with the assays behind it: `tables/proposed_dms_gene_universe.tsv`. Per-assay
detail (phenotype, publication, replicates, direction, quality flags): `tables/dms_candidate_inventory.tsv`
and `tables/dms_assay_harmonization.tsv`.

## Leakage (Step 35)

No information is available on which genes, variants or DMS assays were used to develop AVI or
its components. Every gene is therefore recorded as `{LEAKAGE_UNKNOWN}`, and no gene is excluded
for leakage. This must be revisited before the benchmark is described as independent.
"""
    (p.qc / "gene_universe_report.md").write_text(text)
