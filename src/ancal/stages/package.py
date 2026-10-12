"""Stage 11: attrition, inferential forecast, unresolved issues and the package index
(items 13, 14, 15; protocol Steps 14, 15, 18, 39).

- Variant attrition table: per population, how many variants survive each step towards the
  primary evaluable set, plus the reason every observed variant was lost.
- Expected inferential groups: which groups reach the Step 18 count criteria. The Spearman
  and confidence-interval criteria cannot be evaluated before the freeze and are marked so.
- Unresolved issues: every stage that is not `pass`, every open review item and every known
  limitation in docs/DECISIONS.md.
- deliverables/d1/README.md: one row per item 1-15 with its files, status and scope.

No performance metric is computed here.
"""

from __future__ import annotations

import re

import polars as pl

from .. import budget, config, status
from ..settings import Settings
from . import common

SCOPE_DEPENDENT = {7, 8, 10, 11, 12, 13, 14}
NOT_EVALUABLE = "not evaluable before the config freeze"


def variant_flags(sites: pl.DataFrame, dms: pl.DataFrame, universe: pl.DataFrame, avi_join: pl.DataFrame | None) -> pl.DataFrame:
    """One row per variant: the flags that decide the primary evaluable set (Step 15) and the first reason it is lost."""
    d = dms.filter(pl.col("consequence") == "missense").group_by("variant_id").agg(
        pl.col("score_set_urn").n_unique().alias("n_dms_assays"), pl.col("score_set_urn").unique().sort().alias("dms_assays"))
    uni = universe.filter(pl.col("inclusion_status") == "proposed_include").select("gene_symbol", pl.lit(True).alias("in_universe"))
    f = sites.join(d, on="variant_id", how="left").join(uni, on="gene_symbol", how="left").with_columns(
        pl.col("n_dms_assays").fill_null(0), pl.col("in_universe").fill_null(False),
        (pl.col("consequence").fill_null("") == "missense").alias("missense"))
    if avi_join is not None:
        f = f.join(avi_join.select("variant_id", "avi_join_class", "avi_valid"), on="variant_id", how="left")
    else:
        f = f.with_columns(pl.lit(None, dtype=pl.String).alias("avi_join_class"), pl.lit(None, dtype=pl.Boolean).alias("avi_valid"))
    has_dms = pl.col("n_dms_assays") > 0
    upper = pl.col("primary_snv") & pl.col("missense") & has_dms & pl.col("in_universe")
    return f.with_columns(
        has_dms.alias("has_dms"), upper.alias("evaluable_ignoring_avi"),
        (upper & pl.col("avi_valid")).alias("evaluable"),
        pl.when(pl.col("qc_status") != "pass").then(pl.col("qc_exclusion_reason"))
        .when(pl.col("variant_type") != "snv").then(pl.lit("not_an_snv"))
        .when(~pl.col("in_scope_region")).then(pl.lit("outside_scope_regions"))
        .when(~pl.col("coding")).then(pl.lit("noncoding"))
        .when(~pl.col("missense")).then(pl.concat_str(pl.lit("coding_"), pl.col("consequence")))
        .when(~has_dms).then(pl.lit("no_eligible_dms_measurement"))
        .when(~pl.col("in_universe")).then(pl.lit("gene_not_in_proposed_universe"))
        .when(pl.col("avi_valid").is_null()).then(pl.lit(""))
        .when(~pl.col("avi_valid")).then(pl.concat_str(pl.lit("avi_"), pl.col("avi_join_class")))
        .otherwise(pl.lit("")).alias("lost_at"),
    )


def attrition(long: pl.DataFrame, flags: pl.DataFrame, have_avi: bool) -> pl.DataFrame:
    """Counts per sample set and population level at each step (Steps 14-15)."""
    f = flags.select("variant_id", "primary_snv", "coding", "missense", "has_dms", "avi_valid", "evaluable_ignoring_avi",
                     "evaluable", "gene_symbol", "dms_assays")
    j = long.filter(pl.col("ac") >= 1).join(f, on="variant_id")
    obs, snv = pl.col("observed"), pl.col("observed") & pl.col("primary_snv")
    avi = pl.col("avi_valid").fill_null(False)
    final = pl.col("evaluable") if have_avi else pl.col("evaluable_ignoring_avi")
    na = pl.lit(None, dtype=pl.UInt32)
    agg = j.group_by("sample_set", "level", "label").agg(
        pl.len().alias("n_alt_allele_count_ge_1"), obs.sum().alias("n_observed_after_qc"), snv.sum().alias("n_primary_snvs"),
        (snv & pl.col("coding")).sum().alias("n_coding"), (snv & pl.col("missense")).sum().alias("n_missense"),
        ((snv & avi).sum() if have_avi else na).alias("n_snvs_with_avi_score"),
        (snv & pl.col("has_dms")).sum().alias("n_snvs_with_dms_measurement"),
        ((snv & avi & pl.col("has_dms")).sum() if have_avi else na).alias("n_snvs_with_avi_and_dms"),
        (obs & pl.col("evaluable_ignoring_avi")).sum().alias("n_evaluable_ignoring_avi"),
        ((obs & pl.col("evaluable")).sum() if have_avi else na).alias("n_evaluable"),
        pl.col("gene_symbol").filter(obs & final).n_unique().alias("n_genes"),
        pl.col("dms_assays").filter(obs & final).explode(empty_as_null=True).drop_nulls().n_unique().alias("n_dms_assays"),
    )
    levels = long.select("sample_set", "level", "label", "population_group", "n_samples").unique()
    order = {"all": 0, "group": 1, "population": 2}
    out = levels.join(agg, on=["sample_set", "level", "label"], how="left")
    fill = [c for c in agg.columns if c.startswith("n_") and (have_avi or c not in ("n_snvs_with_avi_score", "n_snvs_with_avi_and_dms", "n_evaluable"))]
    return (out.with_columns(pl.col(fill).fill_null(0))
            .with_columns(pl.col("level").replace_strict(order, return_dtype=pl.Int32).alias("_o"))
            .sort(pl.col("sample_set") != "unrelated", "_o", "population_group", "label").drop("_o"))


def forecast(att: pl.DataFrame, counts: pl.DataFrame, cfg: config.AnalysisConfig, have_avi: bool) -> pl.DataFrame:
    """Step 18 count criteria per group and population, from the unrelated (primary) sample set."""
    pc = cfg.precision_criteria
    n = pl.col("n_evaluable" if have_avi else "n_evaluable_ignoring_avi")
    people = pl.concat([
        counts.filter(pl.col("level") == "group").select(pl.lit("group").alias("level"), pl.col("population_group").alias("label"),
                                                         "n_post_qc", "n_unrelated"),
        counts.filter(pl.col("level") == "population").select(pl.lit("population").alias("level"),
                                                              pl.col("original_population").alias("label"), "n_post_qc", "n_unrelated"),
    ]).with_columns(pl.col("n_post_qc", "n_unrelated").cast(pl.Int64))
    a = att.filter((pl.col("sample_set") == "unrelated") & (pl.col("level") != "all"))
    ok_v, ok_g = n >= pc.min_evaluable_variants, pl.col("n_genes") >= pc.min_represented_genes
    return a.join(people, on=["level", "label"], how="left").select(
        "level", "label", "population_group", pl.col("n_post_qc").alias("n_individuals"), pl.col("n_unrelated").alias("n_unrelated_individuals"),
        n.alias("n_evaluable_variants"), "n_genes",
        pl.lit("evaluable" if have_avi else "upper bound: AVI scores not yet joined").alias("count_basis"),
        ok_v.alias(f"meets_min_{pc.min_evaluable_variants}_variants"), ok_g.alias(f"meets_min_{pc.min_represented_genes}_genes"),
        pl.lit(NOT_EVALUABLE).alias("valid_spearman"), pl.lit(NOT_EVALUABLE).alias("bootstrap_ci_width_le_%.2f" % pc.max_bootstrap_ci_width),
        pl.when(ok_v & ok_g).then(pl.lit("expected_to_meet_count_criteria")).otherwise(pl.lit("expected_to_fail")).alias("expectation"),
        pl.when(ok_v & ok_g).then(pl.lit(""))
        .when(~ok_v & ~ok_g).then(pl.lit(f"fewer than {pc.min_evaluable_variants} evaluable variants; fewer than {pc.min_represented_genes} genes"))
        .when(~ok_v).then(pl.lit(f"fewer than {pc.min_evaluable_variants} evaluable variants"))
        .otherwise(pl.lit(f"fewer than {pc.min_represented_genes} genes")).alias("criteria_failed"),
        pl.when(pl.col("level") == "group").then(pl.lit("inferential (primary spread estimate)"))
        .otherwise(pl.lit("descriptive only")).alias("role"),
    )


# ------------------------------------------------------------ unresolved issues


def decisions_tables(text: str) -> dict[str, list[list[str]]]:
    """Rows of the first markdown table under each `## ` heading of docs/DECISIONS.md."""
    out: dict[str, list[list[str]]] = {}
    for section in re.split(r"^## ", text, flags=re.M)[1:]:
        title, _, body = section.partition("\n")
        rows = [[c.strip() for c in line.strip().strip("|").split("|")] for line in body.splitlines()
                if line.startswith("|") and not re.match(r"^\|[\s:|-]+\|$", line.strip())]
        if len(rows) > 1:
            out[title.strip()] = rows[1:]
    return out


def unresolved_issues(settings: Settings, statuses: list[dict], mode: str, scope: str) -> tuple[str, int]:
    tables = decisions_tables((settings.root / "docs" / "DECISIONS.md").read_text())
    review = next((v for k, v in tables.items() if "review" in k.lower()), [])
    limits = next((v for k, v in tables.items() if "limitation" in k.lower()), [])
    not_pass = [s for s in statuses if s["status"] != "pass"]
    lines = ["# Unresolved data-access, licensing and review issues (item 15)", "",
             "Generated by `ancal stage11` from the stage status files and `docs/DECISIONS.md`. Do not edit;",
             "change `docs/DECISIONS.md` and re-run. Scope: `" + scope + "`.", "", "## Stages that did not pass", ""]
    if not_pass:
        lines += ["| Stage | Status | Reason |", "|---|---|---|"]
        lines += [f"| {s['stage']} | {s['status']} | {s.get('blocked_reason') or '; '.join(c['name'] + ': ' + c['detail'] for c in s['checks'] if not c['ok'])} |"
                  for s in not_pass]
    else:
        lines.append("Every stage passed.")
    lines += ["", "## Results that are provisional", ""]
    prov = []
    if scope == "smoke":
        prov.append("- **Scope.** Items 7, 8, 10, 11, 12, 13 and 14 were produced at `smoke` scope (a few genes on one "
                    "chromosome). They must be re-run at `dms_genes` scope (`SCALING_UP.md`).")
    if mode != "strict":
        prov.append("- **DMS review.** `review_gates.dms_review` is `provisional`: assays without a reviewer decision use "
                    "the pipeline's proposed direction. Items 10, 11, 13 and 14 depend on it. The config cannot be frozen in this mode.")
    lines += prov or ["None."]
    lines += ["", "## Data-access and licensing limitations", "", "| ID | Limitation | Effect | Mitigation |", "|---|---|---|---|"]
    lines += ["| " + " | ".join(r) + " |" for r in limits]
    lines += ["", "## Open review items", "", "| ID | Item | Proposed resolution | Where | Status |", "|---|---|---|---|---|"]
    lines += ["| " + " | ".join(r) + " |" for r in review]
    return "\n".join(lines) + "\n", len(limits) + len(review) + len(not_pass)


# ---------------------------------------------------------------- package index

ITEMS = [
    (1, "Repository structure and README", "stage0", ["README.md", "SCALING_UP.md", "docs/DECISIONS.md"]),
    (2, "Versioned analysis configuration", "stage1", ["config/analysis.yaml", "reports/qc/config_status.md"]),
    (3, "Dataset and software manifest", "stage2", ["reports/manifests/data_manifest.tsv", "reports/manifests/software_manifest.tsv"]),
    (4, "Sample metadata table", "stage3", ["tables/sample_metadata.tsv"]),
    (5, "Population mapping table", "stage3", ["tables/population_mapping.tsv", "config/population_groups.yaml"]),
    (6, "Counts for all 80 populations and seven groups", "stage3", ["tables/population_counts.tsv", "reports/qc/sample_count_report.md"]),
    (7, "Normalized population-variant table", "stage7", ["data_store/variants/population_variants_long.{scope}.parquet",
        "data_store/variants/population_sites.{scope}.parquet", "tables/population_variant_summary.tsv",
        "reports/qc/variant_normalization_report.md", "reports/exclusions/population_variant_exclusions.tsv"]),
    (8, "AVI join and coverage report", "stage8", ["reports/qc/avi_join_report.md", "tables/avi_coverage.tsv", "tables/avi_join_classes.tsv"]),
    (9, "Candidate DMS gene and assay inventory", "stage5", ["tables/dms_candidate_inventory.tsv", "reports/qc/dms_inventory_report.md",
        "config/dms_curation.tsv", "reports/review/dms_direction_review.md"]),
    (10, "Proposed frozen DMS gene universe", "stage9", ["tables/proposed_dms_gene_universe.tsv", "reports/qc/gene_universe_report.md"]),
    (11, "DMS harmonization report", "stage6", ["reports/qc/dms_harmonization_report.md", "tables/dms_assay_harmonization.tsv",
        "data_store/dms/dms_harmonized.{scope}.parquet", "reports/exclusions/dms_harmonization_exclusions.tsv"]),
    (12, "ClinVar coding-missense eligibility report", "stage10", ["reports/qc/clinvar_eligibility_report.md",
        "tables/clinvar_eligible_variants.tsv", "tables/clinvar_eligibility_counts.tsv", "tables/clinvar_eligibility_by_population.tsv"]),
    (13, "Variant attrition table", "stage11", ["tables/variant_attrition.tsv", "reports/exclusions/variant_attrition_exclusions.tsv",
        "reports/qc/attrition_and_eligibility_report.md"]),
    (14, "Groups expected to meet inferential criteria", "stage11", ["tables/expected_inferential_groups.tsv"]),
    (15, "Unresolved data-access or licensing issues", "stage11", ["deliverables/d1/unresolved_issues.md"]),
]
QC_REPORTS = [  # protocol Step 39 -> where each report lives
    ("Sample count report", "reports/qc/sample_count_report.md"),
    ("Population mapping report", "tables/population_mapping.tsv"),
    ("Variant normalization report", "reports/qc/variant_normalization_report.md"),
    ("Allele-conflict report", "reports/qc/variant_normalization_report.md (section 'Allele-conflict report')"),
    ("AVI join report", "reports/qc/avi_join_report.md"),
    ("DMS join report", "reports/qc/attrition_and_eligibility_report.md"),
    ("ClinVar join report", "reports/qc/clinvar_eligibility_report.md"),
    ("Score-direction verification", "reports/review/dms_direction_review.md (DMS); reports/qc/avi_join_report.md (AVI)"),
    ("Gene-universe report", "reports/qc/gene_universe_report.md"),
    ("Frequency-bin report", "reports/qc/variant_normalization_report.md; tables/population_variant_summary.tsv"),
    ("Inferential eligibility report", "tables/expected_inferential_groups.tsv"),
    ("Full exclusion log", "reports/exclusions/; data_store/analysis_datasets/variant_exclusion_log.{scope}.parquet"),
    ("Data-version manifest", "reports/manifests/data_manifest.tsv"),
]


def package_index(settings: Settings, statuses: dict[str, dict], scope: str, mode: str, own_status: str) -> tuple[str, list[int]]:
    rows, missing_items = [], []
    for n, title, stage, files in ITEMS:
        st = own_status if stage == "stage11" else (statuses.get(stage) or {}).get("status", "not run")
        paths = [f.format(scope=scope) for f in files]
        present = [f for f in paths if (settings.root / f).exists()]
        notes = []
        if n in SCOPE_DEPENDENT:
            item_scope = (statuses.get(stage) or {}).get("scope", scope)
            notes.append(f"scope `{item_scope}`" + ("; rerun at `dms_genes`" if item_scope == "smoke" else ""))
        if n in (10, 11, 13, 14) and mode != "strict":
            notes.append("provisional: DMS review not complete")
        if n == 2:
            frozen = sorted((settings.root / "config" / "frozen").glob("analysis_*.yaml"))
            notes.append(f"frozen: `config/frozen/{frozen[-1].name}`" if frozen else "draft, not yet frozen")
        if st == "blocked":
            notes.append((statuses[stage].get("blocked_reason") or "")[:160])
        if st == "pass" and len(present) < len(paths):
            missing = sorted(set(paths) - set(present))
            if not (n == 12 and missing == ["tables/clinvar_eligibility_by_population.tsv"]):
                missing_items.append(n)
            notes.append("missing: " + ", ".join(f"`{m}`" for m in missing))
        rows.append(f"| {n} | {title} | {st} | " + "<br>".join(f"`{f}`" for f in present) + " | " + "; ".join(x for x in notes if x) + " |")
    qc = "\n".join(f"| {name} | {path.format(scope=scope)} |" for name, path in QC_REPORTS)
    text = f"""# Deliverable 1: data-readiness package

Generated by `ancal stage11`. Do not edit. Current scope: `{scope}`. DMS review gate: `{mode}`.

This package contains **no comparative performance result** (no AVI-DMS correlation, no CADD
comparison, no ClinVar AUROC). Those follow the Step 40 review and the config freeze.
Start with `docs/DECISIONS.md` for the reasoning behind every item, and
`deliverables/d1/unresolved_issues.md` for what is still open.

| # | Item | Status | Files | Notes |
|---|---|---|---|---|
{chr(10).join(rows)}

`pass` means the stage's gate checks passed. `blocked` means the stage cannot run yet for the
reason shown; it is not a failure. Files under `data_store/` are not in git; they are rebuilt by
the pipeline and recorded with checksums in `reports/manifests/data_manifest.tsv`.

## Quality-control reports required by protocol Step 39

| Report | Where |
|---|---|
{qc}
"""
    return text, missing_items


# ------------------------------------------------------------------------- run


def run(settings: Settings) -> status.StageStatus:
    budget.check(settings)
    ds = config.load_data_sources(settings.root)
    cfg = config.load_analysis(settings.root)
    p = common.Paths(settings)
    mode = ds.review_gates.dms_review
    notes = "Attrition table, inferential forecast, unresolved issues and package index (items 13-15)."
    statuses = {s["stage"]: s for s in status.read_all(settings) if s["stage"] != "stage11"}
    checks: list[status.Check] = []
    have_core = p.sites.exists() and p.long.exists() and p.dms.exists() and (p.tables / "proposed_dms_gene_universe.tsv").exists()
    have_avi = p.avi_join.exists()

    if have_core:
        sites, long, dms = pl.read_parquet(p.sites), pl.read_parquet(p.long), pl.read_parquet(p.dms)
        universe = pl.read_csv(p.tables / "proposed_dms_gene_universe.tsv", separator="\t", infer_schema_length=0)
        avi_join = pl.read_parquet(p.avi_join) if have_avi else None
        flags = variant_flags(sites, dms, universe, avi_join)
        att = attrition(long, flags, have_avi)
        counts = pl.read_csv(p.tables / "population_counts.tsv", separator="\t", infer_schema_length=0)
        fc = forecast(att, counts, cfg, have_avi)
        log = flags.filter(pl.col("ac_all_eligible") >= 1).select(
            "variant_id", "variant_type", "source_filter", "gene_symbol", "consequence", "af_bin", "n_dms_assays", "avi_join_class",
            "evaluable_ignoring_avi", "evaluable", "lost_at")
        p.exclusion_log.parent.mkdir(parents=True, exist_ok=True)
        log.write_parquet(p.exclusion_log)
        lost = (log.group_by("lost_at").len().rename({"len": "n_variants"})
                .with_columns(pl.col("lost_at").replace("", "retained (evaluable" + ("" if have_avi else ", AVI pending") + ")"))
                .sort("n_variants", descending=True))
        outs = [common.write_tsv(att, p.tables / "variant_attrition.tsv"),
                common.write_tsv(fc, p.tables / "expected_inferential_groups.tsv"),
                common.write_tsv(lost, p.exclusions / "variant_attrition_exclusions.tsv")]
        common.record(settings, "stage11", "package", f"variants/{p.sites.name}, dms/{p.dms.name}, tables/proposed_dms_gene_universe.tsv",
                      p.exclusion_log, *outs)
        mono = att.filter((pl.col("n_observed_after_qc") > pl.col("n_alt_allele_count_ge_1")) | (pl.col("n_primary_snvs") > pl.col("n_observed_after_qc"))
                          | (pl.col("n_coding") > pl.col("n_primary_snvs")) | (pl.col("n_missense") > pl.col("n_coding"))
                          | (pl.col("n_evaluable_ignoring_avi") > pl.col("n_snvs_with_dms_measurement"))).height
        groups = fc.filter(pl.col("level") == "group")
        checks += [
            status.Check("attrition_table", att.height > 0 and mono == 0, f"{att.height} rows; counts never increase along the funnel"),
            status.Check("exclusion_log_complete",
                         log.filter((pl.col("lost_at") == "") != pl.col("evaluable" if have_avi else "evaluable_ignoring_avi")).height == 0,
                         f"{log.height:,} variants with an alternate allele in an eligible sample; each is retained or has a reason"),
            status.Check("seven_groups_forecast", groups.height == 7, f"{groups.filter(pl.col('expectation') == 'expected_to_meet_count_criteria').height} "
                         f"of {groups.height} groups expected to meet the count criteria"
                         + ("" if have_avi else " (upper bound; AVI pending)")),
        ]
        write_report(p, ds.scope, att, fc, lost, have_avi, cfg, mode)
    else:
        checks.append(status.Check("inputs_for_attrition", False, "Stages 6, 7 and 9 must run before the attrition table can be built"))

    text, n_issues = unresolved_issues(settings, list(statuses.values()), mode, ds.scope)
    p.d1.mkdir(parents=True, exist_ok=True)
    (p.d1 / "unresolved_issues.md").write_text(text)
    checks.append(status.Check("unresolved_issues", n_issues > 0, f"{n_issues} entries (stages not passed, limitations, review items)"))
    failed = [s for s, v in statuses.items() if v["status"] == "fail"]
    not_run = [f"stage{i}" for i in range(0, 11) if f"stage{i}" not in statuses]
    checks.append(status.Check("earlier_stages", not failed and not not_run,
                               "every earlier stage is pass or blocked" if not failed and not not_run else f"failed: {failed}; not run: {not_run}"))
    u = budget.usage(settings)
    checks.append(status.Check("disk_budget", u.ok, f"{u.used_mb:.1f} MB of {settings.local_budget_gb} GB"))
    own = "pass" if all(c.ok for c in checks) else "fail"
    index, missing_items = package_index(settings, statuses, ds.scope, mode, own)
    checks.append(status.Check("package_files_present", not missing_items, "every passed item has its files"
                               if not missing_items else f"items with missing files: {missing_items}"))
    if missing_items:
        index, _ = package_index(settings, statuses, ds.scope, mode, "fail")
    (p.d1 / "README.md").write_text(index)
    common.record(settings, "stage11", "package", "reports/status/*.json, docs/DECISIONS.md", p.d1 / "README.md", p.d1 / "unresolved_issues.md")
    return status.write(settings, "stage11", checks, notes=notes)


def write_report(p, scope, att, fc, lost, have_avi, cfg, mode):
    caveats = []
    if not have_avi:
        caveats.append("AVI scores are not joined (Stage 8 blocked), so evaluable counts are upper bounds that ignore the AVI requirement")
    if mode != "strict":
        caveats.append("the DMS review is not complete")
    u = att.filter(pl.col("sample_set") == "unrelated").drop("sample_set", "population_group")
    pc = cfg.precision_criteria
    g = fc.filter(pl.col("level") == "group").drop("level", "population_group", "role")
    text = f"""# Variant attrition and inferential-eligibility report (Stage 11, protocol Steps 14, 15 and 18)

Generated by `ancal stage11`. Do not edit. Scope: `{scope}`.

{common.provisional_banner(scope, caveats)}This is the DMS join report: counts only, no correlation.

## Attrition by group (unrelated sample set)

Each column is a subset of the one before it, except the AVI and DMS columns, which are counted
independently among primary SNVs.

| Column | Meaning |
|---|---|
| `n_alt_allele_count_ge_1` | variants with at least one alternate allele in the population, before QC |
| `n_observed_after_qc` | ... that pass site QC and have allele number above 0 (Step 7) |
| `n_primary_snvs` | ... that are SNVs inside the scope regions |
| `n_coding` / `n_missense` | ... in a MANE Select CDS / missense |
| `n_snvs_with_avi_score` | primary SNVs with a valid AVI score{'' if have_avi else ' (not available)'} |
| `n_snvs_with_dms_measurement` | primary SNVs that are missense with a harmonized DMS measurement |
| `n_snvs_with_avi_and_dms` | both{'' if have_avi else ' (not available)'} |
| `n_evaluable_ignoring_avi` | missense, DMS measurement, gene in the proposed universe |
| `n_evaluable` | ... and a valid AVI score: the primary evaluable set (Step 15){'' if have_avi else ' (not available)'} |
| `n_genes`, `n_dms_assays` | genes and assays represented in the evaluable set |

{common.md_table(u.filter(pl.col('level') != 'population').drop('level'))}

All 80 populations and the all-eligible sample set: `tables/variant_attrition.tsv`.

## Where observed variants are lost

One row per reason, over variants with an alternate allele in at least one eligible sample. Each
variant is counted once, at the first step it fails. Per variant:
`data_store/analysis_datasets/{p.exclusion_log.name}`.

{common.md_table(lost, ['Lost at', 'Variants'])}

## Groups expected to meet the inferential criteria (Step 18)

Criteria: at least {pc.min_evaluable_variants} evaluable variants, at least {pc.min_represented_genes} represented genes, a valid Spearman
correlation, and a bootstrap 95% confidence-interval width of at most {pc.max_bootstrap_ci_width}. Only the two count
criteria can be checked now. The other two need the correlation itself and are {NOT_EVALUABLE}.

{common.md_table(g)}

The 80 original populations are in `tables/expected_inferential_groups.tsv` (descriptive only).
"""
    (p.qc / "attrition_and_eligibility_report.md").write_text(text)
