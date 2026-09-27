import gzip

import polars as pl
import pytest

from ancal.stages import dms_inventory as di


@pytest.mark.parametrize("h,cls", [
    ("p.Arg175His", "missense"), ("p.(Arg175His)", "missense"), ("p.Arg175Ter", "nonsense"),
    ("p.Arg175*", "nonsense"), ("p.Arg175=", "synonymous"), ("p.=", "synonymous"), ("_sy", "synonymous"),
    ("p.Arg175Arg", "synonymous"), ("p.[Arg175His;Gly176Ser]", "other"), ("p.Arg175fs", "other"),
    ("NA", "other"), (None, "other"), ("_wt", "other"),
])
def test_classify_hgvs_pro(h, cls):
    assert di.classify_hgvs_pro(h) == cls


def test_auc():
    assert di.auc([3, 4, 5], [0, 1, 2]) == 1.0
    assert di.auc([0, 1], [5, 6]) == 0.0
    assert di.auc([1, 2], [1, 2]) == 0.5


def _cal(*classes):
    return {"scoreCalibrations": [{"title": "T", "functionalClassifications": [
        {"functionalClassification": c, "range": r} for c, r in classes]}]}


def test_calibration_direction():
    low_bad = _cal(("abnormal", [None, -0.7]), ("normal", [-0.5, None]))
    assert di.calibration_evidence(low_bad, None)[0] == "lower_is_more_damaging"
    high_bad = _cal(("normal", [None, 0.2]), ("abnormal", [0.8, 1.2]))
    assert di.calibration_evidence(high_bad, None)[0] == "higher_is_more_damaging"
    overlap = _cal(("normal", [0, 1]), ("abnormal", [0.5, 2]))
    assert di.calibration_evidence(overlap, None)[0] is None
    abnormal_only = _cal(("abnormal", [0.8, 1.2]))
    assert di.calibration_evidence(abnormal_only, 0.1)[0] == "higher_is_more_damaging"
    assert di.calibration_evidence({"scoreCalibrations": []}, 0.1) == (None, "none")


def _ctl(direction):
    return di.ControlStats(100, 80, 10, 10, 0, 0, 0, 0, 0.9, direction, "")


def test_propose():
    assert di.propose("lower_is_more_damaging", _ctl("lower_is_more_damaging")) == ("lower_is_more_damaging", "high")
    assert di.propose("lower_is_more_damaging", _ctl("higher_is_more_damaging")) == ("conflict", "none")
    assert di.propose(None, _ctl("higher_is_more_damaging")) == ("higher_is_more_damaging", "medium")
    assert di.propose(None, None) == ("unknown", "none")


def test_control_stats_reads_direction(tmp_path):
    rows = ["accession,hgvs_nt,hgvs_pro,score"]
    rows += [f"x#{i},NA,p.Arg{i}Ter,{-2 - i * 0.01}" for i in range(12)]
    rows += [f"y#{i},NA,p.Arg{i}=,{0 + i * 0.01}" for i in range(12)]
    rows += [f"z#{i},NA,p.Arg{i}His,{-1}" for i in range(30)] + ["w#1,NA,p.Arg1His,NA"]
    p = tmp_path / "s.csv.gz"
    p.write_bytes(gzip.compress("\n".join(rows).encode()))
    c = di.control_stats(p)
    assert c.direction == "lower_is_more_damaging" and c.auc_nonsense_vs_synonymous == 0.0
    assert (c.n_nonsense, c.n_synonymous, c.n_missense, c.n_scored) == (12, 12, 30, 54)


def _rec(urn, name="BRCA1", cat="protein_coding", ids=(), n=100, **kw):
    ext = [{"identifier": {"dbName": db, "identifier": i}} for db, i in ids]
    return {"urn": urn, "targetGenes": [{"name": name, "category": cat, "externalIdentifiers": ext}],
            "numVariants": n, **kw}


def test_resolve_gene_priority_and_conflict():
    mane = {"BRCA1": "ENSG00000012048.1", "TP53": "ENSG00000141510.1"}
    by_id = {"ENSG00000141510": "TP53", "ENSG00000012048": "BRCA1"}
    uni = {"P04637": "TP53"}
    assert di.resolve_gene(_rec("u1", "p53", ids=[("UniProt", "P04637")]), {}, mane, by_id, uni) == ("TP53", "uniprot", "")
    g, how, note = di.resolve_gene(_rec("u2", "BRCA1 RING", ids=[("UniProt", "P04637")]), {"u2": ["BRCA1"]}, mane, by_id, uni)
    assert (g, how) == ("BRCA1", "mavedb_mapping") and note.startswith("conflict")
    assert di.resolve_gene(_rec("u3", "alpha-synuclein"), {}, mane, by_id, uni)[0] is None


def test_prescreen_reasons():
    mane = {"BRCA1": "x"}
    assert di.prescreen(_rec("a"), "BRCA1", mane, {}) is None
    assert di.prescreen(_rec("a", cat="regulatory"), "BRCA1", mane, {}).startswith("noncoding_target")
    assert di.prescreen(_rec("a"), "BRCA1", mane, {"a": "b"}) == "superseded_by b"
    assert di.prescreen(_rec("a", metaAnalyzesScoreSetUrns=["x"]), "BRCA1", mane, {}).startswith("meta_analysis")
    assert di.prescreen(_rec("a", n=5), "BRCA1", mane, {}).startswith("fewer_than")
    assert di.prescreen(_rec("a"), None, mane, {}) == "gene_unresolved"
    assert di.prescreen(_rec("a"), "IGHG1", mane, {}).startswith("gene_not_in_mane")


def test_merge_curation_keeps_reviewer_columns(tmp_path):
    inv = pl.DataFrame({c: ["v"] for c in [
        "gene_symbol", "prescreen_status", "exclusion_reason", "title", "experimental_phenotype", "assay_method",
        "calibration_direction", "control_direction", "proposed_direction", "direction_confidence", "quality_flags"]}
    ).with_columns(pl.lit("urn:1").alias("score_set_urn"), pl.lit(True).alias("in_scope"),
                   pl.lit(10).alias("n_variants"), pl.lit(0.9).alias("auc_nonsense_vs_synonymous"))
    path = tmp_path / "cur.tsv"
    first = di.merge_curation(inv, path)
    assert first["reviewer_decision"].to_list() == [None]
    first.with_columns(pl.lit("include").alias("reviewer_decision"), pl.lit("Ana").alias("reviewer_name")) \
        .write_csv(path, separator="\t")
    again = di.merge_curation(inv.with_columns(pl.lit("changed").alias("title")), path)
    assert again["reviewer_decision"].to_list() == ["include"] and again["title"].to_list() == ["changed"]
