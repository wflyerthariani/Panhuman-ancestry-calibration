"""Stage 6: DMS harmonization (item 11, protocol Step 13).

For every in-scope assay that the review (config/dms_curation.tsv) includes:
  1. parse each measurement to one protein change, or one coding SNV where the assay
     reports single-nucleotide variants on the MANE Select CDS;
  2. confirm the reference amino acid / base against GRCh38 (MANE Select);
  3. resolve duplicate measurements of the same change (mean, count kept);
  4. orient the score so that higher = more damaging; keep the original score;
  5. compute the within-assay percentile over missense changes;
  6. place each change on GRCh38: every SNV that encodes it.

Review gate (D-035): with `review_gates.dms_review: strict`, the stage is blocked until every
in-scope assay has a reviewer decision. With `provisional`, unreviewed assays use the
pipeline's proposed direction and every output is marked PROVISIONAL.
"""

from __future__ import annotations

import gzip
import io
import re

import polars as pl

from .. import budget, config, consequence, provenance, status
from ..settings import Settings
from . import common
from .dms_inventory import DIRECTIONS, _SINGLE, _protein_column, classify_hgvs_pro
from .reference import read_contig

MIN_REF_MATCH = 0.95        # below this an assay's numbering is not MANE Select numbering
MIN_OFFSET_RESIDUES = 20    # distinct residues needed before a numbering offset is inferred
NT_MODE_MIN_FRACTION = 0.5  # share of scored rows written as a single coding SNV to use nucleotide mapping
AM_MIN_AGREEMENT = 0.98
AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H",
       "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
       "Tyr": "Y", "Val": "V", "Ter": "*", "*": "*"}
_NT = re.compile(r"^(?:[\w.]+:)?c\.(\d+)([ACGT])>([ACGT])$")

OUT_SCHEMA = {
    "variant_id": pl.String, "chrom": pl.String, "pos": pl.Int64, "ref": pl.String, "alt": pl.String,
    "gene_symbol": pl.String, "ensembl_gene_id": pl.String, "transcript_id": pl.String, "score_set_urn": pl.String,
    "mapping_method": pl.String, "protein_pos": pl.Int64, "ref_aa": pl.String, "alt_aa": pl.String,
    "consequence": pl.String, "dms_score_original": pl.Float64, "dms_score_harmonized": pl.Float64,
    "dms_percentile": pl.Float64, "n_measurements": pl.Int64, "n_snvs_for_change": pl.Int64,
    "raw_direction": pl.String, "direction_source": pl.String,
}


# ------------------------------------------------------------------ review gate


def assay_decisions(curation: pl.DataFrame, mode: str) -> tuple[pl.DataFrame, list[str]]:
    """Per in-scope candidate: decision, direction and where they came from. Also the urns still unreviewed."""
    rows, unreviewed = [], []
    cand = curation.filter((pl.col("in_scope") == "true") & (pl.col("prescreen_status") == "candidate"))
    for r in cand.iter_rows(named=True):
        decision, direction = (r["reviewer_decision"] or "").strip().lower(), (r["reviewer_direction"] or "").strip()
        if decision in ("include", "exclude") and (decision == "exclude" or direction in DIRECTIONS):
            source, reason = "reviewer", "" if decision == "include" else "excluded_by_reviewer"
        else:
            unreviewed.append(r["score_set_urn"])
            source, direction = "provisional_proposal", r["proposed_direction"]
            if mode != "provisional":
                decision, reason = "pending", "awaiting_review"
            elif direction in DIRECTIONS:
                decision, reason = "include", ""
            else:
                decision, reason = "exclude", f"no_direction_proposal ({direction}); needs review"
        rows.append({"score_set_urn": r["score_set_urn"], "gene_symbol": r["gene_symbol"], "decision": decision,
                     "raw_direction": direction if decision == "include" else "", "direction_source": source,
                     "direction_confidence": r["direction_confidence"], "assay_exclusion_reason": reason})
    return pl.DataFrame(rows, schema=["score_set_urn", "gene_symbol", "decision", "raw_direction", "direction_source",
                                      "direction_confidence", "assay_exclusion_reason"]), unreviewed


# ------------------------------------------------------------------- one assay


def parse_protein(h: str | None) -> tuple[int, str, str] | str:
    """(position, ref aa, alt aa) for a single substitution, synonymous or nonsense change; else a reason."""
    if not h or h in ("NA", "_wt"):
        return "no_protein_change"
    cls = classify_hgvs_pro(h)
    if cls == "frameshift":
        return "frameshift"
    h = h.split(":", 1)[1] if ":" in h else h
    m = _SINGLE.match(h)
    if cls == "other" or not m:
        return "multiple_changes" if (";" in h or "[" in h) else "other_notation"
    ref, alt = AA3.get(m.group(1)), AA3.get(m.group(1) if m.group(3) == "=" else m.group(3))
    if ref is None or alt is None:
        return "other_notation"
    return int(m.group(2)), ref, alt


def numbering_offset(residues: set[tuple[int, str]], ref_at: dict[int, str]) -> tuple[int, float]:
    """(offset, match rate): the shift that puts an assay's residue numbering on MANE Select numbering.

    Many assays number residues within the construct they tested (a domain, or a tagged
    protein). Offset 0 is used when it already matches. Otherwise the one offset at which at
    least MIN_REF_MATCH of the assay's distinct residues equal the MANE Select residue is
    used, provided there are at least MIN_OFFSET_RESIDUES of them. If none qualifies, the
    result is (0, rate at offset 0) and the caller rejects the assay.
    """
    if not residues or not ref_at:
        return 0, 0.0

    def rate(o: int) -> float:
        return sum(ref_at.get(pos + o) == aa for pos, aa in residues) / len(residues)

    r0 = rate(0)
    if r0 >= MIN_REF_MATCH or len(residues) < MIN_OFFSET_RESIDUES:
        return 0, r0
    lo, hi = min(p for p, _ in residues), max(p for p, _ in residues)
    best = sorted(((rate(o), o) for o in range(1 - lo, max(ref_at) - lo + 1) if o), reverse=True)[:2]
    if best and best[0][0] >= MIN_REF_MATCH and (len(best) == 1 or best[1][0] < MIN_REF_MATCH):
        return best[0][1], best[0][0]
    return 0, r0


def read_scores(path) -> pl.DataFrame:
    raw = gzip.decompress(path.read_bytes())
    return pl.read_csv(io.BytesIO(raw), infer_schema_length=0, null_values=["NA", ""])


def harmonize_assay(df: pl.DataFrame, urn: str, direction: str, cons: pl.DataFrame) -> tuple[pl.DataFrame, dict, dict]:
    """Harmonize one assay against its gene's consequence table.

    Returns (mapped rows, {exclusion reason: measurement count}, summary).
    """
    excl: dict[str, int] = {}
    summ = {"n_rows": df.height, "mapping_method": "", "ref_match_rate": None, "n_changes": 0, "n_missense_changes": 0}

    def drop(reason: str, n: int = 1):
        if n:
            excl[reason] = excl.get(reason, 0) + n

    empty = pl.DataFrame(schema=OUT_SCHEMA)
    if "score" not in df.columns:
        drop("no_score_column", df.height)
        return empty, excl, summ
    scores = df["score"].cast(pl.Float64, strict=False).to_list()
    col = _protein_column(df)
    pro = df[col].to_list() if col else [None] * df.height
    nt = df["hgvs_nt"].to_list() if "hgvs_nt" in df.columns else [None] * df.height
    scored = [(s, p, n) for s, p, n in zip(scores, pro, nt) if s is not None and s == s and abs(s) != float("inf")]
    drop("no_score", df.height - len(scored))
    summ["n_scored"] = len(scored)
    if not scored:
        return empty, excl, summ

    nt_hits = [(_NT.match(n) if n else None) for _, _, n in scored]
    use_nt = sum(m is not None for m in nt_hits) >= NT_MODE_MIN_FRACTION * len(scored)
    if use_nt:
        ref_at = dict(zip(cons["cds_pos"], cons["cds_ref"]))
        keys = [(int(m.group(1)), m.group(2), m.group(3)) for m in nt_hits if m]
        in_range = [k for k in keys if k[0] in ref_at]
        rate = sum(ref_at[k[0]] == k[1] for k in in_range) / len(in_range) if in_range else 0.0
        use_nt = rate >= MIN_REF_MATCH
    units: dict[tuple, list[float]] = {}
    if use_nt:
        summ.update(mapping_method="hgvs_nt_cds", ref_match_rate=rate)
        for (s, _, _), m in zip(scored, nt_hits):
            if m is None:
                drop("not_a_single_coding_snv")
            else:
                units.setdefault((int(m.group(1)), m.group(2), m.group(3)), []).append(s)
        key_cols, cons_keys, ref_at = ["cds_pos", "cds_ref", "cds_alt"], ["cds_pos", "cds_ref", "cds_alt"], ref_at
    else:
        for s, p, _ in scored:
            parsed = parse_protein(p)
            if isinstance(parsed, str):
                drop(parsed)
            else:
                units.setdefault(parsed, []).append(s)
        ref_at = dict(zip(cons["protein_pos"], cons["ref_aa"]))
        offset, rate = numbering_offset({(k[0], k[1]) for k in units}, ref_at)
        if offset:
            units = {(k[0] + offset, k[1], k[2]): v for k, v in units.items()}
        summ.update(mapping_method="hgvs_pro", ref_match_rate=rate if units else None, numbering_offset=offset)
        key_cols, cons_keys = ["protein_pos", "ref_aa", "alt_aa"], ["protein_pos", "ref_aa", "alt_aa"]
        if units and rate < MIN_REF_MATCH:
            drop("assay_numbering_mismatch", sum(len(v) for v in units.values()))
            return empty, excl, summ

    sign = -1.0 if direction == "lower_is_more_damaging" else 1.0
    ok = []
    for key, vals in units.items():
        if key[0] not in ref_at:
            drop("position_out_of_range", len(vals))
        elif ref_at[key[0]] != key[1]:
            drop("reference_mismatch", len(vals))
        else:
            mean = sum(vals) / len(vals)
            ok.append((*key, mean, sign * mean, len(vals)))
    if not ok:
        return empty, excl, summ
    u = pl.DataFrame(ok, orient="row", schema=[*key_cols, "dms_score_original", "dms_score_harmonized", "n_measurements"])
    mapped = u.join(cons, left_on=key_cols, right_on=cons_keys, how="left")
    unreachable = mapped.filter(pl.col("pos").is_null())
    drop("not_reachable_by_one_snv", int(unreachable["n_measurements"].sum()))
    mapped = mapped.filter(pl.col("pos").is_not_null())
    if mapped.height == 0:
        return empty, excl, summ

    # Percentile over missense changes: every scored missense change for protein-level assays
    # (reachable by one SNV or not), every mapped missense SNV for nucleotide-level assays.
    if use_nt:
        pool = mapped.filter(pl.col("consequence") == "missense").select(*key_cols, "dms_score_harmonized")
    else:
        pool = u.filter((pl.col("ref_aa") != pl.col("alt_aa")) & (pl.col("alt_aa") != "*") & (pl.col("protein_pos") > 1)
                        ).select(*key_cols, "dms_score_harmonized")
    pct = pool.with_columns((pl.col("dms_score_harmonized").rank("average") / pl.len()).alias("dms_percentile")
                            ).drop("dms_score_harmonized")
    summ.update(n_changes=u.height, n_missense_changes=pool.height)
    mapped = (
        mapped.join(pct, on=key_cols, how="left")
        .with_columns(
            pl.format("{}:{}:{}:{}", "chrom", "pos", "ref", "alt").alias("variant_id"),
            pl.lit(urn).alias("score_set_urn"), pl.lit(summ["mapping_method"]).alias("mapping_method"),
            pl.len().over(key_cols).alias("n_snvs_for_change"), pl.lit(direction).alias("raw_direction"),
            pl.lit("").alias("direction_source"))
        .select(list(OUT_SCHEMA)).cast(OUT_SCHEMA).sort("chrom", "pos", "alt")
    )
    return mapped, excl, summ


# ------------------------------------------------------------- consequence table


def build_consequences(settings: Settings, p: common.Paths) -> tuple[pl.DataFrame, list[str]]:
    bed = common.read_bed(p.bed)
    genes = common.scope_genes(bed)
    cache: dict[str, str] = {}

    def get_seq(chrom: str) -> str:
        if chrom not in cache:
            cache.clear()  # one contig in memory at a time
            cache[chrom] = read_contig(settings.data_dir("reference") / f"GRCh38.{chrom}.fa.gz", chrom)
        return cache[chrom]

    cons, skipped = consequence.build(pl.read_parquet(p.tx), pl.read_parquet(p.cds), genes, get_seq)
    p.consequences.parent.mkdir(parents=True, exist_ok=True)
    cons.write_parquet(p.consequences)
    return cons, skipped


def alphamissense_agreement(cons: pl.DataFrame, am_path) -> tuple[int, int]:
    """(AlphaMissense rows joined to a missense SNV of the same gene model, rows whose protein change agrees)."""
    am = pl.read_csv(am_path, separator="\t", comment_prefix="#", has_header=False,
                     new_columns=["chrom", "pos", "ref", "alt", "genome", "uniprot_id", "transcript_id",
                                  "protein_variant", "am_pathogenicity", "am_class"], schema_overrides={"pos": pl.Int64})
    j = am.select("chrom", "pos", "ref", "alt", "protein_variant").unique().join(
        cons.filter(pl.col("consequence").is_in(["missense", "start_lost"])), on=["chrom", "pos", "ref", "alt"])
    agree = j.filter(pl.col("protein_variant") == pl.format("{}{}{}", "ref_aa", "protein_pos", "alt_aa"))
    return j.select("chrom", "pos", "ref", "alt").n_unique(), agree.select("chrom", "pos", "ref", "alt").n_unique()


# ------------------------------------------------------------------------- run


def run(settings: Settings) -> status.StageStatus:
    budget.check(settings)
    ds = config.load_data_sources(settings.root)
    p = common.Paths(settings)
    mode = ds.review_gates.dms_review
    notes = "DMS harmonization (item 11)."
    if not p.bed.exists() or not p.tx.exists():
        return status.write(settings, "stage6", [], notes=notes,
                            blocked_reason="Stage 4 outputs for this scope are missing; run `make stage4` first.")
    curation = pl.read_csv(settings.root / "config" / "dms_curation.tsv", separator="\t", infer_schema_length=0)
    decisions, unreviewed = assay_decisions(curation, mode)
    if mode == "strict" and unreviewed:
        return status.write(
            settings, "stage6", [], notes=notes,
            blocked_reason=f"{len(unreviewed)} of {decisions.height} in-scope assays have no reviewer decision in "
                           "config/dms_curation.tsv (guide: reports/review/dms_direction_review.md). To run on the "
                           "pipeline's proposed directions instead, set review_gates.dms_review: provisional in "
                           "config/data_sources.yaml.")

    cons, skipped = build_consequences(settings, p)
    n_am, n_am_ok = alphamissense_agreement(cons, p.alphamissense)
    by_gene = {k[0]: g for k, g in cons.group_by("gene_symbol")}
    scores_dir = settings.data_dir("dms") / "mavedb" / "scores"

    frames, excl_rows, summaries, missing_scores = [], [], [], []
    for d in decisions.iter_rows(named=True):
        urn = d["score_set_urn"]
        summ = {"score_set_urn": urn, "gene_symbol": d["gene_symbol"], "decision": d["decision"],
                "raw_direction": d["raw_direction"], "direction_source": d["direction_source"],
                "direction_confidence": d["direction_confidence"], "assay_exclusion_reason": d["assay_exclusion_reason"]}
        path = scores_dir / f"{urn.replace(':', '_')}.csv.gz"
        if d["decision"] == "include" and not path.exists():
            missing_scores.append(urn)
            summ["assay_exclusion_reason"] = "scores_not_downloaded (re-run Stage 5 at this scope)"
        elif d["decision"] == "include" and d["gene_symbol"] not in by_gene:
            summ["assay_exclusion_reason"] = "gene_cds_not_translatable"
        elif d["decision"] == "include":
            mapped, excl, s = harmonize_assay(read_scores(path), urn, d["raw_direction"], by_gene[d["gene_symbol"]])
            mapped = mapped.with_columns(pl.lit(d["direction_source"]).alias("direction_source"))
            frames.append(mapped)
            excl_rows += [{"score_set_urn": urn, "gene_symbol": d["gene_symbol"], "reason": k, "n_measurements": v}
                          for k, v in sorted(excl.items())]
            miss = mapped.filter(pl.col("consequence") == "missense")
            summ.update(s, n_snvs_mapped=mapped.height, n_missense_snvs_mapped=miss.height)
            if "assay_numbering_mismatch" in excl:
                summ["assay_exclusion_reason"] = "assay_numbering_mismatch"
            elif mapped.height == 0:
                summ["assay_exclusion_reason"] = "no_variant_mapped"
        summaries.append(summ)

    dms = pl.concat(frames) if frames else pl.DataFrame(schema=OUT_SCHEMA)
    assays = pl.DataFrame(summaries, infer_schema_length=None)
    for c, t in (("numbering_offset", pl.Int64), ("n_rows", pl.Int64), ("n_scored", pl.Int64), ("mapping_method", pl.String), ("ref_match_rate", pl.Float64),
                 ("n_changes", pl.Int64), ("n_missense_changes", pl.Int64), ("n_snvs_mapped", pl.Int64),
                 ("n_missense_snvs_mapped", pl.Int64)):
        if c not in assays.columns:
            assays = assays.with_columns(pl.lit(None, dtype=t).alias(c))
    assays = assays.with_columns(
        pl.when(pl.col("assay_exclusion_reason") == "").then(pl.lit("harmonized")).otherwise(pl.lit("not_harmonized"))
        .alias("harmonization_status")).sort("gene_symbol", "score_set_urn")
    excl_df = pl.DataFrame(excl_rows, schema={"score_set_urn": pl.String, "gene_symbol": pl.String,
                                              "reason": pl.String, "n_measurements": pl.Int64})

    dms.write_parquet(p.dms)
    excl_df.write_parquet(p.dms_exclusions)
    t_assays = common.write_tsv(assays, p.tables / "dms_assay_harmonization.tsv")
    t_excl = common.write_tsv(excl_df, p.exclusions / "dms_harmonization_exclusions.tsv")
    common.record(settings, "stage6", "dms_harmonize", "config/dms_curation.tsv, dms/mavedb/scores/*, reference/*",
                  p.consequences, p.dms, p.dms_exclusions, t_assays, t_excl)

    dup = dms.height - dms.select("variant_id", "score_set_urn").n_unique()
    pct = dms["dms_percentile"].drop_nulls()
    harmonized = assays.filter(pl.col("harmonization_status") == "harmonized")
    included = assays.filter(pl.col("decision") == "include")
    checks = [
        status.Check("review_gate", True, f"mode={mode}; {decisions.height - len(unreviewed)}/{decisions.height} "
                     "in-scope assays reviewed" + ("" if not unreviewed else f"; {len(unreviewed)} use proposals (PROVISIONAL)")),
        status.Check("consequence_table", cons.height > 0, f"{cons.height:,} coding SNVs in {cons['gene_symbol'].n_unique()} "
                     "genes" + (f"; CDS not translatable: {skipped}" if skipped else "")),
        status.Check("alphamissense_cross_check", n_am > 0 and n_am_ok / n_am >= AM_MIN_AGREEMENT,
                     f"protein change agrees for {n_am_ok:,} of {n_am:,} AlphaMissense SNVs"),
        status.Check("scores_available", not missing_scores, "all included assays have score files" if not missing_scores
                     else f"missing: {missing_scores[:5]}"),
        status.Check("assays_harmonized", harmonized.height > 0, f"{harmonized.height} of {included.height} included assays"),
        status.Check("unique_variant_per_assay", dup == 0, f"{dup} duplicate (variant, assay) rows"),
        status.Check("percentile_range", pct.len() == 0 or (pct.min() > 0 and pct.max() <= 1), "within (0, 1]"),
    ]
    u = budget.usage(settings)
    checks.append(status.Check("disk_budget", u.ok, f"{u.used_mb:.1f} MB of {settings.local_budget_gb} GB"))
    write_report(settings, p, ds.scope, mode, unreviewed, assays, excl_df, dms, cons, skipped, n_am, n_am_ok)
    return status.write(settings, "stage6", checks, notes=notes + (
        f" PROVISIONAL: {len(unreviewed)} assays not yet reviewed." if unreviewed else ""))


def write_report(settings, p, scope, mode, unreviewed, assays, excl, dms, cons, skipped, n_am, n_am_ok):
    banner = common.provisional_banner(scope, [f"{len(unreviewed)} assay(s) have no reviewer decision and use the "
                                               "pipeline's proposed direction"] if unreviewed else [])
    cons_counts = cons.group_by("gene_symbol", "consequence").len().pivot(on="consequence", index="gene_symbol", values="len").sort("gene_symbol")
    a = assays.select("score_set_urn", "gene_symbol", "decision", "raw_direction", "direction_source", "mapping_method",
                      "numbering_offset", "n_scored", "ref_match_rate", "n_missense_changes", "n_missense_snvs_mapped", "assay_exclusion_reason")
    by_reason = excl.group_by("reason").agg(pl.col("n_measurements").sum(), pl.col("score_set_urn").n_unique().alias("assays")
                                            ).sort(["n_measurements", "reason"], descending=[True, False])
    miss = dms.filter(pl.col("consequence") == "missense")
    per_gene = (dms.group_by("gene_symbol").agg(
        pl.col("score_set_urn").n_unique().alias("assays"), pl.col("variant_id").n_unique().alias("snvs_with_a_measurement"),
        pl.col("variant_id").filter(pl.col("consequence") == "missense").n_unique().alias("missense_snvs_with_a_measurement"),
    ).sort("gene_symbol"))
    text = f"""# DMS harmonization report (Stage 6, protocol Step 13)

Generated by `ancal stage6`. Do not edit. Scope: `{scope}`. Review gate: `{mode}`.

{banner}## Summary

| | Count |
|---|---|
| In-scope candidate assays | {assays.height} |
| Included (reviewer or provisional proposal) | {assays.filter(pl.col('decision') == 'include').height} |
| Harmonized | {assays.filter(pl.col('harmonization_status') == 'harmonized').height} |
| Assays without a reviewer decision | {len(unreviewed)} |
| (SNV, assay) measurements on GRCh38 | {dms.height:,} |
| Distinct SNVs with at least one measurement | {dms['variant_id'].n_unique():,} |
| Distinct missense SNVs with at least one measurement | {miss['variant_id'].n_unique():,} |

## Method

- **Unit of measurement.** A protein change (`hgvs_pro`, or the submitter's protein column), or a
  coding SNV when at least {NT_MODE_MIN_FRACTION:.0%} of an assay's scored rows are single SNVs in `c.` notation
  whose reference bases match the MANE Select CDS.
- **Reference check.** The reference amino acid (or base) of each change is compared with MANE
  Select. Many assays number residues within the construct they tested. If fewer than {MIN_REF_MATCH:.0%} of an
  assay's residues agree as written, the pipeline looks for the single shift (`numbering_offset`)
  at which at least {MIN_REF_MATCH:.0%} of at least {MIN_OFFSET_RESIDUES} distinct residues agree, and applies it. If there is
  none, the assay is not harmonized (`assay_numbering_mismatch`); it is never partially mapped.
  **Every non-zero offset should be confirmed by the reviewer.**
- **Duplicates.** Several measurements of the same change are averaged; the count is kept
  (`n_measurements`). The score used is the assay-reported `score` column.
- **Direction.** `dms_score_harmonized` = the score, negated when the raw direction is
  `lower_is_more_damaging`. `dms_score_original` is kept unchanged.
- **Percentile.** `dms_percentile` = rank of the harmonized score among the assay's missense
  changes, divided by their number. Only missense changes have a percentile.
- **Placement on GRCh38.** A protein change is placed on every SNV that encodes it
  (`n_snvs_for_change`). Changes that need two or three base changes cannot be observed as one SNV
  and are counted as `not_reachable_by_one_snv`; this is expected, not an error.
- **Cross-check.** The pipeline's own protein consequences agree with AlphaMissense's
  `protein_variant` for {n_am_ok:,} of {n_am:,} missense SNVs in scope.

## Coding SNVs per gene (MANE Select)

{common.md_table(cons_counts)}
{f"{chr(10)}CDS not translatable (skipped): {', '.join(skipped)}{chr(10)}" if skipped else ""}
## Assays

{common.md_table(a)}

## Measurements with a GRCh38 placement, per gene

{common.md_table(per_gene)}

## Measurements not carried forward

Counts are measurements (score rows). Per assay: `reports/exclusions/dms_harmonization_exclusions.tsv`.

{common.md_table(by_reason)}

## Outputs

- `data_store/dms/{p.dms.name}`: one row per (SNV, assay).
- `data_store/intermediate/{p.consequences.name}`: every coding SNV of the scope genes with its consequence.
- `tables/dms_assay_harmonization.tsv`: one row per in-scope assay.
"""
    (p.qc / "dms_harmonization_report.md").write_text(text)
