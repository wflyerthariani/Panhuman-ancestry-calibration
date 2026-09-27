"""Stage 5: DMS candidate inventory (item 9, protocol Step 12) and the direction/quality review.

1. Snapshot every published human MaveDB score set (full records).
2. Resolve each to one human gene (MaveDB mapping > Ensembl ID > UniProt > target name).
3. Automatic pre-screen: removes only what the protocol rules out outright (noncoding
   target, multi-gene target, superseded, meta-analysis, < 20 variants, gene unresolved or
   not in MANE Select). Everything else is a candidate.
4. For candidates in scope, download scores and gather direction evidence:
     a. MaveDB calibrations (ranges labelled normal / abnormal)
     b. built-in controls: nonsense vs synonymous variants
     c. author text excerpts (shown to the reviewer, never parsed into a decision)
5. Write the review sheet config/dms_curation.tsv (reviewer columns are preserved across
   re-runs) and the reviewer guide reports/review/dms_direction_review.md.

Direction convention (Step 13): after harmonization, higher = more functionally damaging.
"""

from __future__ import annotations

import gzip
import io
import re
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from .. import budget, config, provenance, status
from ..settings import Settings
from ..sources import uniprot
from ..sources.mavedb import MaveDB

MIN_VARIANTS = 20
MIN_CONTROLS = 10          # nonsense and synonymous variants needed to call direction from controls
AUC_CLEAR = 0.25           # |AUC - 0.5| >= this -> clear separation
REVIEWER_COLUMNS = ["reviewer_decision", "reviewer_direction", "reviewer_notes", "reviewer_name", "reviewed_on"]
DIRECTIONS = ("higher_is_more_damaging", "lower_is_more_damaging")
TEXT_CUES = re.compile(
    r"[^.]*\b(higher|lower|more negative|more positive|increas\w*|decreas\w*|wild[- ]type[- ]like|"
    r"loss[- ]of[- ]function|gain[- ]of[- ]function|normali[sz]ed|null[- ]like|depleted|enriched)\b[^.]*\.",
    re.I,
)
KEYWORDS = {
    "assay_method": "Phenotypic Assay Method",
    "assay_mechanism": "Phenotypic Assay Mechanism",
    "model_system": "Phenotypic Assay Model System",
    "molecular_mechanism": "Molecular Mechanism Assessed",
    "library_method": "Variant Library Creation Method",
}


# --------------------------------------------------------------- gene symbols


def resolve_gene(rec: dict, mapped: dict, mane_by_name: dict, mane_by_id: dict, uni: dict) -> tuple[str | None, str, str]:
    """(gene symbol or None, method, note). Prefers MaveDB's own genomic mapping."""
    tg = rec["targetGenes"][0]
    candidates = []
    m = mapped.get(rec["urn"]) or []
    if len(m) == 1:
        candidates.append((m[0], "mavedb_mapping"))
    for e in tg.get("externalIdentifiers", []):
        ident = e["identifier"]
        if ident["dbName"] == "Ensembl":
            g = mane_by_id.get(ident["identifier"].split(".")[0])
            if g:
                candidates.append((g, "ensembl_id"))
        elif ident["dbName"] == "UniProt":
            g = uni.get(ident["identifier"].split("-")[0])
            if g:
                candidates.append((g, "uniprot"))
    token = re.split(r"[\s_(\-:/,]+", tg["name"].strip())[0].upper()
    if token in mane_by_name:
        candidates.append((token, "target_name"))
    if not candidates:
        return None, "unresolved", f"target name '{tg['name']}'"
    symbols = {c[0] for c in candidates}
    note = "" if len(symbols) == 1 else "conflict: " + "; ".join(f"{g} via {how}" for g, how in candidates)
    return candidates[0][0], candidates[0][1], note


# ------------------------------------------------------------------ pre-screen


def superseded_map(records: list[dict]) -> dict[str, str]:
    """old urn -> newer urn (MaveDB only records the link on the newer record)."""
    out = {}
    for r in records:
        old = r.get("supersededScoreSet")
        if old:
            out[old["urn"] if isinstance(old, dict) else old] = r["urn"]
    return out


def prescreen(rec: dict, gene: str | None, mane_by_name: dict, superseded: dict) -> str | None:
    tgs = rec["targetGenes"]
    if len(tgs) != 1:
        return f"multi_gene_target ({len(tgs)} targets)"
    if tgs[0].get("category") != "protein_coding":
        return f"noncoding_target ({tgs[0].get('category')}; Step 32 - separate benchmark)"
    if rec["urn"] in superseded:
        return f"superseded_by {superseded[rec['urn']]}"
    if rec.get("metaAnalyzesScoreSetUrns"):
        return "meta_analysis (derived from other score sets, not a direct measurement)"
    if (rec.get("numVariants") or 0) < MIN_VARIANTS:
        return f"fewer_than_{MIN_VARIANTS}_variants"
    if gene is None:
        return "gene_unresolved"
    if gene not in mane_by_name:
        return f"gene_not_in_mane_select ({gene})"
    return None


# --------------------------------------------------------- evidence: calibration


def calibration_evidence(rec: dict, score_median: float | None) -> tuple[str | None, str]:
    """Direction implied by MaveDB calibrations, plus a one-line summary."""
    normal, abnormal, titles = [], [], []
    for cal in rec.get("scoreCalibrations") or []:
        titles.append(cal.get("title", "untitled") + (" [research use only]" if cal.get("researchUseOnly") else ""))
        for fc in cal.get("functionalClassifications") or []:
            lo, hi = fc.get("range") or [None, None]
            lo = float("-inf") if lo is None else lo
            hi = float("inf") if hi is None else hi
            if fc.get("functionalClassification") == "normal":
                normal.append((lo, hi))
            elif fc.get("functionalClassification") == "abnormal":
                abnormal.append((lo, hi))
    if not titles:
        return None, "none"
    summary = "; ".join(titles)
    if abnormal and normal:
        ab_hi, nm_lo = max(h for _, h in abnormal), min(l for l, _ in normal)
        ab_lo, nm_hi = min(l for l, _ in abnormal), max(h for _, h in normal)
        if ab_hi <= nm_lo:
            return "lower_is_more_damaging", f"abnormal range(s) below normal range(s) — {summary}"
        if ab_lo >= nm_hi:
            return "higher_is_more_damaging", f"abnormal range(s) above normal range(s) — {summary}"
        return None, f"abnormal and normal ranges overlap — {summary}"
    if abnormal and score_median is not None:
        mids = [(max(l, -1e9) + min(h, 1e9)) / 2 for l, h in abnormal]
        if all(m > score_median for m in mids):
            return "higher_is_more_damaging", f"abnormal range(s) above the score median (no normal range given) — {summary}"
        if all(m < score_median for m in mids):
            return "lower_is_more_damaging", f"abnormal range(s) below the score median (no normal range given) — {summary}"
    return None, f"calibration present but direction not determinable — {summary}"


# ------------------------------------------------------------ evidence: controls

_SINGLE = re.compile(r"^p\.\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|=|\*)\)?$")


def classify_hgvs_pro(h: str | None) -> str:
    if not h or h in ("NA", "_wt", "_sy"):
        return "synonymous" if h == "_sy" else "other"
    if h in ("p.=", "p.(=)"):
        return "synonymous"
    m = _SINGLE.match(h)
    if not m:
        return "other"  # multi-mutants, indels, frameshifts, p.? ...
    alt = m.group(3)
    if alt == "=" or alt == m.group(1):
        return "synonymous"
    if alt in ("Ter", "*"):
        return "nonsense"
    return "missense"


def auc(pos: list[float], neg: list[float]) -> float:
    """P(score of a positive > score of a negative), ties counting 1/2 (Mann-Whitney)."""
    ranks = pl.Series([*pos, *neg]).rank("average").to_list()
    r_pos = sum(ranks[: len(pos)])
    return (r_pos - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


@dataclass
class ControlStats:
    n_scored: int
    n_missense: int
    n_synonymous: int
    n_nonsense: int
    median_missense: float | None
    median_synonymous: float | None
    median_nonsense: float | None
    median_all: float | None
    auc_nonsense_vs_synonymous: float | None
    direction: str | None
    summary: str


def control_stats(scores_csv: Path) -> ControlStats:
    raw = gzip.decompress(scores_csv.read_bytes())
    df = pl.read_csv(io.BytesIO(raw), infer_schema_length=0, null_values=["NA", ""])
    if "score" not in df.columns:
        return ControlStats(0, 0, 0, 0, None, None, None, None, None, None, "no 'score' column")
    df = df.with_columns(pl.col("score").cast(pl.Float64, strict=False)).filter(
        pl.col("score").is_not_null() & pl.col("score").is_finite())
    hgvs = df["hgvs_pro"] if "hgvs_pro" in df.columns else pl.Series([None] * df.height, dtype=pl.String)
    df = df.with_columns(pl.Series("cls", [classify_hgvs_pro(h) for h in hgvs.to_list()]))

    def med(c):
        s = df.filter(pl.col("cls") == c)["score"]
        return float(s.median()) if s.len() else None

    syn = df.filter(pl.col("cls") == "synonymous")["score"].to_list()
    non = df.filter(pl.col("cls") == "nonsense")["score"].to_list()
    n_mis = df.filter(pl.col("cls") == "missense").height
    a, direction = None, None
    if len(syn) >= MIN_CONTROLS and len(non) >= MIN_CONTROLS:
        a = auc(non, syn)
        if a >= 0.5 + AUC_CLEAR:
            direction = "higher_is_more_damaging"
        elif a <= 0.5 - AUC_CLEAR:
            direction = "lower_is_more_damaging"
        summary = (f"nonsense (n={len(non)}, median {med('nonsense'):.3g}) vs synonymous (n={len(syn)}, median "
                   f"{med('synonymous'):.3g}); AUC {a:.2f}")
    elif "hgvs_pro" not in df.columns or df.filter(pl.col("cls") != "other").height == 0:
        summary = "no protein-level HGVS; controls not classifiable"
    else:
        summary = f"too few controls (nonsense {len(non)}, synonymous {len(syn)}; need {MIN_CONTROLS} each)"
    return ControlStats(df.height, n_mis, len(syn), len(non), med("missense"), med("synonymous"), med("nonsense"),
                        float(df["score"].median()) if df.height else None, a, direction, summary)


# ----------------------------------------------------------------- author text


def author_excerpts(rec: dict, limit: int = 3) -> list[str]:
    text = " ".join(filter(None, [rec.get("shortDescription"), rec.get("methodText"), rec.get("abstractText")]))
    text = re.sub(r"\s+", " ", re.sub(r"[#*`]", "", text))
    hits = []
    for m in TEXT_CUES.finditer(text):
        s = m.group(0).strip()
        if "score" in s.lower() and s not in hits:
            hits.append(s[:300])
        if len(hits) == limit:
            break
    return hits


# -------------------------------------------------------------- proposal logic


def propose(cal_dir: str | None, ctl: ControlStats | None) -> tuple[str, str]:
    ctl_dir = ctl.direction if ctl else None
    if cal_dir and ctl_dir:
        return (cal_dir, "high") if cal_dir == ctl_dir else ("conflict", "none")
    if cal_dir or ctl_dir:
        return cal_dir or ctl_dir, "medium"
    return "unknown", "none"


def quality_flags(rec: dict, ctl: ControlStats | None, gene_note: str) -> list[str]:
    flags = []
    if ctl is None:
        flags.append("scores_not_downloaded_at_this_scope")
    else:
        if ctl.auc_nonsense_vs_synonymous is None:
            flags.append("no_nonsense_synonymous_controls")
        elif abs(ctl.auc_nonsense_vs_synonymous - 0.5) < AUC_CLEAR:
            flags.append("weak_control_separation")
        if ctl.n_missense < MIN_VARIANTS:
            flags.append(f"fewer_than_{MIN_VARIANTS}_scored_missense")
    ms = rec.get("mappingState")
    if ms != "complete":
        flags.append(f"genomic_mapping_{ms or 'missing'}")
    lic = rec["license"]["shortName"]
    if "NC" in lic:
        flags.append(f"licence_{lic.replace(' ', '_')}")
    if gene_note:
        flags.append("gene_resolution_conflict")
    return flags


def replicate_info(rec: dict) -> str:
    cols = (rec.get("datasetColumns") or {}).get("scoreColumns") or []
    rep = [c for c in cols if re.search(r"rep|replicate", c, re.I)]
    err = [c for c in cols if re.search(r"(^|_)(se|sd|sem|std|stderr|var|variance|ci)($|_)", c, re.I)]
    parts = []
    if rep:
        parts.append("replicate columns: " + ", ".join(rep[:6]))
    if err:
        parts.append("error columns: " + ", ".join(err[:6]))
    return "; ".join(parts) or "none reported in score columns"


def keywords(rec: dict) -> dict[str, str]:
    kws = (rec.get("experiment") or {}).get("keywords") or []
    got = {}
    for k in kws:
        key, label = k.get("keyword", {}).get("key"), k.get("keyword", {}).get("label")
        for col, want in KEYWORDS.items():
            if key == want and label:
                got[col] = label
    return got


def publication(rec: dict) -> str:
    pubs = rec.get("primaryPublicationIdentifiers") or []
    if pubs:
        p = pubs[0]
        ref = f"{p.get('dbName')}:{p.get('identifier')}"
        return f"{ref} ({p['doi']})" if p.get("doi") else ref
    dois = rec.get("doiIdentifiers") or []
    return f"doi:{dois[0]['identifier']}" if dois else ""


# ------------------------------------------------------------------------- run


def run(settings: Settings, refresh: bool = False) -> status.StageStatus:
    budget.check(settings)
    ds = config.load_data_sources(settings.root)
    mave = MaveDB(settings, ds.sources["mavedb"].uri, refresh=refresh)
    records = mave.human_score_sets()
    mapped = mave.mapped_genes()

    tx = pl.read_parquet(settings.data_dir("reference") / "mane_select_transcripts.parquet")
    mane_by_name = dict(zip(tx["gene_name"], tx["gene_id"]))
    mane_by_id = {gid.split(".")[0]: name for name, gid in mane_by_name.items()}
    accs = {e["identifier"]["identifier"] for r in records for t in r["targetGenes"]
            for e in t.get("externalIdentifiers", []) if e["identifier"]["dbName"] == "UniProt"}
    uni = uniprot.gene_names(settings, ds.sources["uniprot"].uri, accs)
    superseded = superseded_map(records)

    scope_genes = ds.scopes[ds.scope].genes
    in_scope = (lambda g: True) if not isinstance(scope_genes, list) else (lambda g: g in scope_genes)

    rows = []
    for rec in sorted(records, key=lambda r: r["urn"]):
        gene, method, note = (None, "multi_gene_target", "") if len(rec["targetGenes"]) != 1 else \
            resolve_gene(rec, mapped, mane_by_name, mane_by_id, uni)
        reason = prescreen(rec, gene, mane_by_name, superseded)
        scope_hit = gene is not None and in_scope(gene)
        ctl = None
        if reason is None and scope_hit:
            ctl = control_stats(mave.scores_csv(rec["urn"]))
        cal_dir, cal_summary = calibration_evidence(rec, ctl.median_all if ctl else None)
        prop, conf = propose(cal_dir, ctl) if reason is None else ("n/a", "n/a")
        flags = quality_flags(rec, ctl, note) if reason is None else []
        kw = keywords(rec)
        rows.append({
            "score_set_urn": rec["urn"],
            "gene_symbol": gene,
            "ensembl_gene_id": mane_by_name.get(gene) if gene else None,
            "gene_resolution": method + (f" ({note})" if note else ""),
            "title": rec["title"],
            "experimental_phenotype": (rec.get("shortDescription") or "").strip()[:300],
            "assay_method": kw.get("assay_method", ""),
            "assay_mechanism": kw.get("assay_mechanism", ""),
            "molecular_mechanism": kw.get("molecular_mechanism", ""),
            "model_system": kw.get("model_system", ""),
            "publication": publication(rec),
            "license": rec["license"]["shortName"],
            "published_date": rec.get("publishedDate"),
            "n_variants": rec.get("numVariants"),
            "replicate_availability": replicate_info(rec),
            "genomic_mapping_state": rec.get("mappingState"),
            "calibrations": cal_summary,
            "calibration_direction": cal_dir or "",
            "control_summary": ctl.summary if ctl else "not computed at this scope",
            "control_direction": (ctl.direction or "") if ctl else "",
            "auc_nonsense_vs_synonymous": round(ctl.auc_nonsense_vs_synonymous, 3)
            if ctl and ctl.auc_nonsense_vs_synonymous is not None else None,
            "n_scored": ctl.n_scored if ctl else None,
            "n_missense_scored": ctl.n_missense if ctl else None,
            "proposed_direction": prop,
            "direction_confidence": conf,
            "quality_flags": ",".join(flags),
            "in_scope": scope_hit,
            "prescreen_status": "excluded" if reason else "candidate",
            "exclusion_reason": reason or "",
            "_excerpts": author_excerpts(rec) if (reason is None and scope_hit) else [],
        })
    inv = pl.DataFrame(rows, infer_schema_length=None)

    tables = settings.root / "tables"
    inv.drop("_excerpts").write_csv(tables / "dms_candidate_inventory.tsv", separator="\t")
    curation_path = settings.root / "config" / "dms_curation.tsv"
    curation = merge_curation(inv, curation_path)
    curation.write_csv(curation_path, separator="\t")

    review = inv.filter(pl.col("in_scope") & (pl.col("prescreen_status") == "candidate"))
    write_review_guide(settings, review, curation, ds.scope)
    write_inventory_report(settings, inv, ds.scope)
    for p in (tables / "dms_candidate_inventory.tsv", curation_path):
        provenance.record(settings, p, stage="stage5", source_name="derived", command="ancal stage5",
                          software_name="ancal.stages.dms_inventory",
                          software_version=provenance.code_version(settings.root),
                          input_uri="dms/mavedb/score_sets.json.gz, mapped_genes.json.gz, scores/*")

    genes = scope_genes if isinstance(scope_genes, list) else []
    missing = [g for g in genes if review.filter(pl.col("gene_symbol") == g).height == 0]
    reviewed = curation.filter(pl.col("in_scope") & (pl.col("prescreen_status") == "candidate")
                               & (pl.col("reviewer_decision").fill_null("") != ""))
    checks = [
        status.Check("mavedb_snapshot", len(records) > 0, f"{len(records)} published human score sets"),
        status.Check("every_record_classified", inv.height == len(records),
                     f"{(inv['prescreen_status'] == 'candidate').sum()} candidates, "
                     f"{(inv['prescreen_status'] == 'excluded').sum()} excluded by pre-screen"),
        status.Check("scope_genes_have_candidates", not missing,
                     f"{review.height} in-scope candidate assays" + (f"; none for {missing}" if missing else "")),
        status.Check("disk_budget", budget.usage(settings).ok, f"{budget.usage(settings).used_mb:.1f} MB"),
    ]
    notes = (f"DMS inventory (item 9). REVIEW STOP: {reviewed.height}/{review.height} in-scope assays reviewed. "
             "Biologist review: reports/review/dms_direction_review.md -> config/dms_curation.tsv.")
    return status.write(settings, "stage5", checks, notes=notes)


# ------------------------------------------------------------ curation sheet


def merge_curation(inv: pl.DataFrame, path: Path) -> pl.DataFrame:
    """Regenerate the evidence columns; keep whatever reviewers have written."""
    cols = ["score_set_urn", "gene_symbol", "in_scope", "prescreen_status", "exclusion_reason", "title",
            "experimental_phenotype", "assay_method", "n_variants", "calibration_direction", "control_direction",
            "auc_nonsense_vs_synonymous", "proposed_direction", "direction_confidence", "quality_flags"]
    fresh = inv.select(cols)
    if path.exists():
        old = pl.read_csv(path, separator="\t", infer_schema_length=0).select(["score_set_urn", *REVIEWER_COLUMNS])
        fresh = fresh.join(old, on="score_set_urn", how="left")
    else:
        fresh = fresh.with_columns([pl.lit(None, dtype=pl.String).alias(c) for c in REVIEWER_COLUMNS])
    return fresh.sort(pl.col("in_scope").not_(), pl.col("prescreen_status"), "gene_symbol", "score_set_urn",
                      nulls_last=True)


# ----------------------------------------------------------------- reports


def _fmt(v, digits=3):
    return "—" if v is None else (f"{v:.{digits}g}" if isinstance(v, float) else str(v))


def write_review_guide(settings, review: pl.DataFrame, curation: pl.DataFrame, scope: str):
    by_urn = {r["score_set_urn"]: r for r in curation.iter_rows(named=True)}
    counts = review.group_by("proposed_direction").len().sort("proposed_direction")
    cards = []
    for i, r in enumerate(review.sort("gene_symbol", "score_set_urn").iter_rows(named=True), start=1):
        excerpts = "\n".join(f"   > {e}" for e in r["_excerpts"]) or "   > (no sentence about score meaning found — see the MaveDB page)"
        done = by_urn.get(r["score_set_urn"], {}).get("reviewer_decision") or "**not yet reviewed**"
        cards.append(f"""### {i}. {r['gene_symbol']} — {r['title']}

`{r['score_set_urn']}` · [MaveDB page](https://www.mavedb.org/score-sets/{r['score_set_urn']}) · {r['publication'] or 'no publication listed'} · licence {r['license']}

- **What was measured:** {r['experimental_phenotype'] or '—'}
- **Assay:** {r['assay_method'] or '—'}{(' · ' + r['assay_mechanism']) if r['assay_mechanism'] else ''}{(' · ' + r['model_system']) if r['model_system'] else ''}
- **Variants:** {r['n_variants']} in MaveDB, {_fmt(r['n_scored'])} with a numeric score, {_fmt(r['n_missense_scored'])} single missense
- **Replicates:** {r['replicate_availability']}
- **Genome mapping by MaveDB:** {r['genomic_mapping_state']}

| Evidence | What it shows | Implied direction |
|---|---|---|
| A. Calibration | {r['calibrations']} | {r['calibration_direction'] or '—'} |
| B. Controls | {r['control_summary']} | {r['control_direction'] or '—'} |

C. What the authors say about the score:
{excerpts}

**Proposed:** `{r['proposed_direction']}` (confidence: {r['direction_confidence']}) · **Flags:** {r['quality_flags'] or 'none'} · **Review status:** {done}
""")
    text = f"""# DMS assay review: score direction and assay quality

**Who this is for:** someone with a molecular biology or genetics background. No coding needed.
**Time needed:** about 5–10 minutes per assay. **Assays to review at this scope (`{scope}`):** {review.height}.

## Why this review is needed

The audit compares AVI predictions with deep mutational scanning (DMS) measurements. Each DMS
assay reports a score per variant, but assays disagree on what a *high* score means. In a
growth assay, a high score can mean the variant works (it grows). In a drug-resistance screen,
a high score can mean the variant is broken (it escapes the drug). Before scores can be pooled,
every assay must be oriented so that **a higher score always means more functional damage**
(protocol Step 13). The protocol also requires an explicit judgement of **assay quality**
(Step 12).

A wrong direction silently turns a good predictor into a bad one for that gene, so a person
has to sign off on each assay. The pipeline has gathered the evidence and made a proposal.
Your job is to confirm or correct it.

## The three kinds of evidence on each card

- **A. Calibration.** Some MaveDB score sets have been calibrated, often against ClinVar
  variants, into ranges labelled *normal* and *abnormal* function. If the abnormal range sits
  below the normal range, lower scores mean damage. Only a minority of assays have one.
- **B. Built-in controls.** Nonsense variants (premature stop, e.g. `p.Arg213Ter`) almost always
  destroy function, and synonymous variants (e.g. `p.Arg213=`) almost never do. The pipeline
  compares their scores:
  - **AUC** is the probability that a random nonsense variant scores higher than a random
    synonymous one.
  - **AUC near 1:** high score = damaging. **AUC near 0:** low score = damaging.
  - **AUC near 0.5:** the assay does not separate them, which is a warning sign for quality,
    or the assay measures something where truncation is not "damaging" in the assay's terms.
  - **Caveat:** nonsense variants near the protein's C-terminus may escape nonsense-mediated
    decay and behave like missense variants. Assays that start downstream of the protein may
    show little separation.
- **C. Author text.** Sentences from the MaveDB description that mention the score. They are
  shown verbatim; the pipeline does not interpret them. Please check them against A and B,
  and follow the MaveDB link or the paper when they are unclear.

## What to decide for each assay

Fill in these columns of **`config/dms_curation.tsv`**. It opens in Excel or Google Sheets.
The rows under review have `in_scope = true` and `prescreen_status = candidate`, and they are
listed first.

| Column | Allowed values | Meaning |
|---|---|---|
| `reviewer_decision` | `include` / `exclude` | Should this assay enter the DMS benchmark? |
| `reviewer_direction` | `higher_is_more_damaging` / `lower_is_more_damaging` | The direction of the **raw** MaveDB score. Required when `include`. The pipeline reverses the score if it is `lower_is_more_damaging`. |
| `reviewer_notes` | free text | **Required** if you change the proposed direction, exclude an assay, or have a caveat (e.g. "measures abundance only, not activity"). |
| `reviewer_name` | your name | |
| `reviewed_on` | YYYY-MM-DD | |

Please do **not** edit the other columns. They are regenerated each time the pipeline runs, and
your five columns are kept.

### Questions to ask yourself

1. **Does the assay measure a loss of the protein's normal function?** For example: growth that
   depends on the protein, reporter activity, abundance or stability, binding, or drug response
   through the protein.
   - Assays of gain of function, or of a neomorphic or dominant-negative effect, can still be
     included if "damage" is defined consistently.
   - Say so in the notes. Such assays deserve extra care over direction.
2. **Do A, B and C agree on the direction?** If they conflict, or the proposal is `unknown`,
   decide from the paper and explain in the notes.
3. **Is the quality acceptable?** Warning signs:
   - `weak_control_separation`, or an AUC between about 0.35 and 0.65
   - very few scored missense variants
   - no replicate or error information
   - an assay the authors themselves describe as noisy or preliminary
4. **Is it the right version?** When several score sets come from the same experiment (for
   example nucleotide-level and amino-acid-level scores, or per-condition and combined scores),
   prefer the **amino-acid-level, primary** score. Exclude the others with the note
   "duplicate of <urn>". The protocol pools variants within a gene, so counting the same
   measurement twice would bias results.
5. **Anything that makes this assay unsuitable for comparing populations?** For example: only a
   short domain was tested, or an unusual isoform or construct was used, such as a common
   polymorphism in the backbone (`TP53 (P72R)`). Note it; it matters for Stage 9.

### Flags you may see

| Flag | Meaning |
|---|---|
| `weak_control_separation` | Nonsense vs synonymous AUC is between 0.25 and 0.75 |
| `no_nonsense_synonymous_controls` | Fewer than {MIN_CONTROLS} of either class had a score, or no protein-level HGVS |
| `fewer_than_{MIN_VARIANTS}_scored_missense` | Too few single missense variants to be useful |
| `genomic_mapping_incomplete` / `_failed` | MaveDB could not map every variant to GRCh38. Stage 6 falls back to protein-level mapping |
| `licence_CC_BY-NC-SA_4.0` | Non-commercial licence (see limitation L-05) |
| `gene_resolution_conflict` | Different identifiers pointed to different genes. Please check the target |

## Proposal summary

{chr(10).join(f"- `{r[0]}`: {r[1]}" for r in counts.iter_rows())}

## Assays to review

{chr(10).join(cards)}

## When you are done

Tell Arman. The next pipeline stage (DMS harmonization) refuses to run until every in-scope
assay has a `reviewer_decision`, and every included assay has a `reviewer_direction`.
"""
    out = settings.root / "reports" / "review" / "dms_direction_review.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text)


def write_inventory_report(settings, inv: pl.DataFrame, scope: str):
    by_reason = (inv.filter(pl.col("prescreen_status") == "excluded")
                 .with_columns(pl.col("exclusion_reason").str.replace(r" \(.*$| urn:.*$", "").alias("reason"))
                 .group_by("reason").len().sort(["len", "reason"], descending=[True, False]))
    cand = inv.filter(pl.col("prescreen_status") == "candidate")
    res = inv.group_by(pl.col("gene_resolution").str.replace(r" \(.*$", "")).len().sort(["len", "gene_resolution"], descending=[True, False])
    lic = inv.group_by("license").len().sort(["len", "license"], descending=[True, False])
    mapping = cand.group_by("genomic_mapping_state").len().sort(["len", "genomic_mapping_state"], descending=[True, False], nulls_last=True)
    rows = lambda df: "\n".join("| " + " | ".join(_fmt(v) for v in r) + " |" for r in df.iter_rows())  # noqa: E731
    text = f"""# DMS candidate inventory report (Stage 5, protocol Step 12)

Generated by `ancal stage5`. Do not edit. Scope: `{scope}`. The MaveDB snapshot is recorded in
`reports/manifests/data_manifest.tsv` (re-runs reuse it; `ancal stage5 --refresh` takes a new one).

| | Count |
|---|---|
| Published human score sets in MaveDB | {inv.height} |
| Candidates after automatic pre-screen | {cand.height} |
| Distinct candidate genes | {cand['gene_symbol'].n_unique()} |
| Candidates in the current scope | {cand.filter(pl.col('in_scope')).height} |

## Automatic pre-screen exclusions

These follow directly from the protocol; nothing here involves judging quality. Every excluded
score set stays in `tables/dms_candidate_inventory.tsv` with its reason.

| Reason | Score sets |
|---|---|
{rows(by_reason)}

## How target genes were identified

`mavedb_mapping` means MaveDB's own GRCh38 mapping. Then, in order: an Ensembl gene ID matched to
MANE Select, the UniProt accession's primary gene name, and the first word of the target name.
Conflicts between methods are flagged for review.

| Method | Score sets |
|---|---|
{rows(res)}

## Genomic mapping of candidates (used in Stage 6)

| MaveDB mapping state | Candidates |
|---|---|
{rows(mapping)}

## Licences

| Licence | Score sets |
|---|---|
{rows(lic)}

## Not decided here

- **Direction and assay quality** need a biologist's review (`reports/review/dms_direction_review.md`).
- **"At least 20 eligible variants after harmonization"** is applied in Stage 9, once variants
  are joined to the population data.
- **Known data leakage** (assays used to train AVI's components) is assessed in Stage 9 and in
  the leakage log (Step 35).
"""
    (settings.root / "reports" / "qc" / "dms_inventory_report.md").write_text(text)
