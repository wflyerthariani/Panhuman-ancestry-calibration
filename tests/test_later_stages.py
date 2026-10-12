"""Stages 7-11: table logic on small hand-made inputs (no network, no bcftools)."""

import polars as pl

from ancal import config
from ancal.settings import repo_root
from ancal.sources import avi as avi_src
from ancal.sources import remote_vcf
from ancal.stages import avi_join, clinvar_eligibility as ce, common, gene_universe, package
from ancal.stages import population_variants as pv

CFG = config.load_analysis(repo_root())


# ------------------------------------------------------------------- stage 7


def _samples():
    rows = [("s1", "PopA", "AFR", "primary"), ("s2", "PopA", "AFR", "sensitivity_only"), ("s3", "PopB", "EUR", "primary"),
            ("s4", "PopB", "EUR", "excluded")]
    return pl.DataFrame(rows, orient="row", schema=["sample_id", "original_population", "population_group", "inclusion_status"]
                        ).with_columns((pl.col("inclusion_status") != "excluded").alias("post_qc"), pl.lit(True).alias("in_vcf"))


def test_sample_groups_keep_related_samples_out_of_the_primary_set():
    groups, membership = pv.sample_groups(_samples())
    assert "s4" not in membership
    assert membership["s1"] == ["u_all", "u_gAFR", "u_p001", "a_all", "a_gAFR", "a_p001"]
    assert membership["s2"] == ["a_all", "a_gAFR", "a_p001"]
    n = {g.token: g.n_samples for g in groups}
    assert n == {"u_all": 2, "u_gAFR": 1, "u_gEUR": 1, "u_p001": 1, "u_p002": 1,
                 "a_all": 3, "a_gAFR": 2, "a_gEUR": 1, "a_p001": 2, "a_p002": 1}


def test_variant_type():
    assert [pv.variant_type(*x) for x in (("A", "G"), ("AT", "A"), ("A", "*"), ("AT", "GC"), ("A", "<DEL>"))] == \
        ["snv", "indel", "star_or_symbolic", "mnv", "star_or_symbolic"]


def test_norm_counts():
    assert pv.norm_counts("x\nLines   total/split/joined/realigned/removed/skipped:\t61/2/0/1/0/0\n") == \
        {"total": 61, "split": 2, "joined": 0, "realigned": 1, "removed": 0, "skipped": 0}


def _wide(groups, records):
    """records: (pos, ref, alt, filter, {token: (AN, AC, AC_Hom, AC_Het, NS)})"""
    rows = []
    for pos, ref, alt, flt, counts in records:
        row = {"chrom": "chr1", "pos": pos, "ref": ref, "alt": alt, "source_filter": flt, "old_rec": None}
        for g in groups:
            an, ac, hom, het, ns = counts.get(g.token, (2 * g.n_samples, 0, 0, 0, g.n_samples))
            row.update({f"AN_{g.token}": an, f"AC_{g.token}": ac, f"AC_Hom_{g.token}": hom, f"AC_Het_{g.token}": het,
                        f"AC_Hemi_{g.token}": 0, f"NS_{g.token}": ns})
        rows.append(row)
    return pl.DataFrame(rows).with_columns(pl.col("old_rec").cast(pl.String))


def test_build_tables_qc_observed_and_frequency():
    groups, _ = pv.sample_groups(_samples())
    seq = "ACGTACGTAC"
    het_a = {"u_all": (4, 1, 0, 1, 2), "u_gAFR": (2, 1, 0, 1, 1), "u_p001": (2, 1, 0, 1, 1),
             "a_all": (6, 1, 0, 1, 3), "a_gAFR": (4, 1, 0, 1, 2), "a_p001": (4, 1, 0, 1, 2)}
    wide = _wide(groups, [
        (2, "C", "T", "PASS", het_a),
        (2, "C", "T", "PASS", het_a),                    # exact duplicate record
        (3, "G", "A", "AS_VQSR", het_a),                 # fails the source filter
        (4, "A", "C", "PASS", het_a),                    # REF is not the GRCh38 base (T)
        (5, "AC", "A", "PASS", {"u_p002": (0, 0, 0, 0, 0)}),   # indel; nobody called in PopB
    ])
    region = pl.DataFrame({"chrom": ["chr1"] * 4, "pos": [2, 3, 4, 5]})
    cons = pl.DataFrame(schema={"chrom": pl.String, "pos": pl.Int64, "ref": pl.String, "alt": pl.String, "gene_symbol": pl.String,
                                "ensembl_gene_id": pl.String, "transcript_id": pl.String, "strand": pl.String, "cds_pos": pl.Int64,
                                "cds_ref": pl.String, "cds_alt": pl.String, "protein_pos": pl.Int64, "ref_aa": pl.String,
                                "alt_aa": pl.String, "consequence": pl.String})
    sites, long = pv.build_tables(wide, groups, seq, region, cons, CFG)
    s = {r["variant_id"]: r for r in sites.iter_rows(named=True)}
    assert sites.height == 4 and long.height == 4 * len(groups)
    assert s["chr1:2:C:T"]["n_duplicate_records_removed"] == 1 and s["chr1:2:C:T"]["primary_snv"]
    assert s["chr1:2:C:T"]["af_global"] == 0.25 and s["chr1:2:C:T"]["af_bin"] == "common"
    assert s["chr1:2:C:T"]["n_populations_observed"] == 1 and s["chr1:2:C:T"]["n_groups_observed"] == 1
    assert s["chr1:3:G:A"]["qc_exclusion_reason"] == "fails_source_filter" and not s["chr1:3:G:A"]["primary_snv"]
    assert s["chr1:4:A:C"]["qc_exclusion_reason"] == "reference_allele_mismatch"
    assert s["chr1:5:AC:A"]["variant_type"] == "indel" and not s["chr1:5:AC:A"]["primary_snv"]
    row = long.filter((pl.col("variant_id") == "chr1:2:C:T") & (pl.col("label") == "PopA") & (pl.col("sample_set") == "unrelated")).row(0, named=True)
    assert (row["ac"], row["an"], row["af"], row["n_het"], row["call_rate"], row["observed"]) == (1, 2, 0.5, 1, 1.0, True)
    an0 = long.filter((pl.col("variant_id") == "chr1:5:AC:A") & (pl.col("label") == "PopB") & (pl.col("sample_set") == "unrelated")).row(0, named=True)
    assert an0["qc_exclusion_reason"] == "allele_number_zero" and not an0["callable"] and not an0["observed"]
    assert not long.filter(pl.col("variant_id") == "chr1:3:G:A")["observed"].any()
    assert pv.consistency(long) == (0, 0)


def test_tbi_bins_and_bgzf_reader(tmp_path, monkeypatch):
    import gzip, struct, zlib

    def bgzf(data: bytes) -> bytes:
        c = zlib.compressobj(6, zlib.DEFLATED, -15)
        body = c.compress(data) + c.flush()
        size = 18 + len(body) + 8
        return (b"\x1f\x8b\x08\x04" + b"\0" * 4 + b"\0\xff" + struct.pack("<H", 6) + b"BC" + struct.pack("<HH", 2, size - 1)
                + body + struct.pack("<II", zlib.crc32(data), len(data)))

    b1 = b"##h\n#CHROM\tPOS\tID\tREF\tALT\n" + b"chr1\t100\t.\tA\tG\n"
    b2 = b"chr1\t200\t.\tAT\tA\nchr1\t300\t.\tC\tT\nchr2\t50\t.\tG\tC\n"
    blob = bgzf(b1) + bgzf(b2)
    off2 = len(bgzf(b1))

    class R:
        def __init__(self, content):
            self.content = content

    def fake(method, uri, headers):
        a, b = map(int, headers["Range"].removeprefix("bytes=").split("-"))
        return R(blob[a : b + 1])

    monkeypatch.setattr(remote_vcf, "_request", fake)
    assert remote_vcf.header("u") == "##h\n#CHROM\tPOS\tID\tREF\tALT\n"
    # index: chr1 has one chunk in the smallest bin (4681) starting at the first record
    first_rec = len(b"##h\n#CHROM\tPOS\tID\tREF\tALT\n")
    idx = (b"TBI\x01" + struct.pack("<8i", 1, 2, 1, 2, 0, ord("#"), 0, 5) + b"chr1\0"
           + struct.pack("<i", 1) + struct.pack("<Ii", 4681, 1) + struct.pack("<QQ", first_rec, (off2 << 16) | 40)
           + struct.pack("<i", 1) + struct.pack("<Q", first_rec))
    p = tmp_path / "x.tbi"
    p.write_bytes(gzip.compress(idx))
    tbi = remote_vcf.Tbi(p)
    assert tbi.start_offset("chr1", 150, 250) == first_rec and tbi.start_offset("chrZ", 0, 10) is None
    got = [l.split(b"\t")[1] for l in remote_vcf.region_lines("u", tbi, "chr1", 199, 300)]
    assert got == [b"200", b"300"]                       # 0-based half-open [199, 300); stops before chr2
    assert [l.split(b"\t")[1] for l in remote_vcf.region_lines("u", tbi, "chr1", 200, 201)] == [b"200"]  # deletion overlaps


# ------------------------------------------------------------------- stage 8


def test_avi_columns_and_join_classes():
    src = config.Source(name="AVI", uri="x", access=None, format="tsv", score_field="avi")
    assert avi_src.resolve_columns(["CHROM", "POS", "REF", "ALT", "avi"], src)["pos"] == "POS"
    try:
        avi_src.resolve_columns(["c", "POS", "REF", "ALT", "avi"], src)
        raise AssertionError("expected AviConfigError")
    except avi_src.AviConfigError as e:
        assert "avi.columns.chrom" in str(e)

    raw = pl.DataFrame([
        ("chr1", 10, "A", "G", 0.9, None), ("1", 11, "c", "t", 0.2, None),              # exact; corrected
        ("chr1", 12, "T", "A", 0.5, None),                                               # population has A>T: swapped
        ("chr1", 13, "G", "C", 0.1, "t1"), ("chr1", 13, "G", "C", 0.7, "t2"),          # two different scores
        ("chr1", 14, "C", "T", None, None),                                              # null score
        ("chr1", 15, "G", "A", 0.3, None),                                               # AVI REF G; GRCh38 and population have T
    ], orient="row", schema=avi_src.RAW_SCHEMA)
    scores = avi_join.normalize(raw)
    sites = pl.DataFrame([("chr1", 10, "A", "G"), ("chr1", 11, "C", "T"), ("chr1", 12, "A", "T"), ("chr1", 13, "G", "C"),
                          ("chr1", 14, "C", "T"), ("chr1", 15, "T", "C"), ("chr1", 16, "A", "C"), ("chr1", 17, "G", "T")],
                         orient="row", schema=["chrom", "pos", "ref", "alt"]).with_columns(
        pl.format("{}:{}:{}:{}", "chrom", "pos", "ref", "alt").alias("variant_id"),
        (pl.col("pos") != 17).alias("ref_matches_grch38"))
    j = avi_join.classify_joins(sites, scores)
    assert j["avi_join_class"].to_list() == [
        "exact_allele_match", "corrected_representation_match", "allele_swapped_match", "multiple_possible_matches",
        "no_avi_score", "unresolved_allele_conflict", "no_avi_score", "unresolved_allele_conflict"]
    assert j["avi_score"].to_list() == [0.9, 0.2, None, None, None, None, None, None]
    assert j["avi_valid"].sum() == 2


def test_coverage_is_scored_over_eligible_observed():
    sites = pl.DataFrame({"variant_id": ["a", "b", "c"], "primary_snv": [True, True, False], "af_bin": ["rare", "common", "rare"],
                          "coding": [True, False, True], "consequence": ["missense", None, "missense"], "gene_symbol": ["G", None, "G"]})
    join = pl.DataFrame({"variant_id": ["a", "b", "c"], "avi_valid": [True, False, True]})
    long = pl.DataFrame({"variant_id": ["a", "b", "c", "a"], "sample_set": ["unrelated"] * 4, "level": ["group"] * 3 + ["population"],
                         "label": ["AFR"] * 3 + ["PopA"], "population_group": ["AFR"] * 4, "observed": [True, True, True, False]})
    cov = avi_join.coverage(long, sites, join)
    overall = cov.filter(pl.col("stratifier") == "overall").row(0, named=True)
    assert (overall["label"], overall["n_eligible_observed"], overall["n_with_valid_avi_score"], overall["avi_coverage"]) == ("AFR", 2, 1, 0.5)
    assert overall["no_evidence_rate"] == 0.5 and cov.filter(pl.col("label") == "PopA").height == 0


# ------------------------------------------------------------------- stage 9


def test_gene_universe_criteria():
    inv = pl.DataFrame({"gene_symbol": ["A", "A", "B", "C", "D"], "ensembl_gene_id": ["EA", "EA", "EB", "EC", "ED"],
                        "prescreen_status": ["candidate"] * 5})
    assays = pl.DataFrame({"score_set_urn": ["a1", "a2", "b1", "c1"], "gene_symbol": ["A", "A", "B", "C"],
                           "decision": ["include", "exclude", "include", "exclude"],
                           "harmonization_status": ["harmonized", "not_harmonized", "harmonized", "not_harmonized"],
                           "direction_source": ["reviewer", "reviewer", "provisional_proposal", "reviewer"]})
    dms = pl.DataFrame({"gene_symbol": ["A"] * 3 + ["B"], "variant_id": ["v1", "v2", "v3", "w1"], "consequence": ["missense"] * 4})
    sites = pl.DataFrame({"variant_id": ["v1", "v2", "w1"], "primary_snv": [True, True, True], "ac_global": [3, 0, 1]})
    u = gene_universe.build(inv, assays, dms, sites, None, ["A", "B", "C"], 2, "harmonized_missense_snvs")
    r = {x["gene_symbol"]: x for x in u.iter_rows(named=True)}
    assert r["A"]["inclusion_status"] == "proposed_include" and r["A"]["n_population_observed_missense_snvs"] == 1
    assert r["B"]["exclusion_reason"] == "fewer_than_2_eligible_variants" and r["B"]["caveat"] == "provisional_unreviewed_assays"
    assert r["C"]["exclusion_reason"] == "no_assay_with_a_direction_and_acceptable_quality"
    assert r["D"]["inclusion_status"] == "not_assessed" and r["A"]["leakage_status"] == "overlap_unknown"
    assert r["A"]["n_observed_with_valid_avi_score"] is None
    by_obs = gene_universe.build(inv, assays, dms, sites, None, ["A", "B", "C"], 2, "population_observed_missense_snvs")
    assert by_obs.filter(pl.col("gene_symbol") == "A")["inclusion_status"][0] == "proposed_exclude"


# ------------------------------------------------------------------ stage 10


def test_clinvar_significance_and_eligibility():
    pos, neg = CFG.clinvar.positive_class, CFG.clinvar.negative_class
    f = lambda s: ce.significance_class(s, pos, neg)  # noqa: E731
    assert [f(x) for x in ("Pathogenic", "Likely_benign|other", "Pathogenic/Likely_pathogenic", "Uncertain_significance",
                           "Conflicting_classifications_of_pathogenicity", "drug_response", None)] == \
        ["positive", "negative", "positive", "vus", "conflicting", "other", "not_provided"]

    cons = pl.DataFrame([("chr1", 10, "A", "G", "G1", "E1", "T1", "+", 1, "A", "G", 2, "K", "R", "missense"),
                         ("chr1", 11, "C", "T", "G1", "E1", "T1", "+", 2, "C", "T", 2, "K", "K", "synonymous"),
                         ("chr1", 12, "G", "T", "G1", "E1", "T1", "+", 3, "G", "T", 2, "K", "N", "missense")],
                        orient="row", schema=list(__import__("ancal.consequence", fromlist=["SCHEMA"]).SCHEMA))
    cv = pl.DataFrame([("chr1", 10, "A", "G", "Pathogenic"), ("chr1", 11, "C", "T", "Benign"), ("chr1", 12, "G", "T", "Uncertain_significance"),
                       ("chr1", 12, "C", "T", "Benign"), ("chr1", 13, "AT", "A", "Pathogenic"), ("chr1", 99, "A", "C", "Likely_benign")],
                      orient="row", schema=["chrom", "pos", "ref", "alt", "clnsig"]).with_columns(pl.lit(None, dtype=pl.String).alias("mc"))
    out = ce.classify(cv, cons, CFG, None)
    assert out["exclusion_reason"].to_list() == ["", "not_coding_missense", "uncertain_significance", "no_valid_allele_match",
                                                 "not_an_snv", "not_coding_missense"]
    assert out["eligibility"].to_list() == ["eligible_pending_avi"] + ["excluded"] * 5
    avi = pl.DataFrame({"chrom": ["chr1"], "pos": [10], "ref": ["A"], "alt": ["G"], "avi_score": [None]}, schema_overrides={"avi_score": pl.Float64})
    with_avi = ce.classify(cv, cons, CFG, avi)
    assert with_avi["exclusion_reason"][0] == "no_avi_score" and with_avi["eligibility"][0] == "excluded"


# ------------------------------------------------------------------ stage 11


def test_attrition_flags_and_forecast():
    sites = pl.DataFrame({
        "variant_id": ["m1", "m2", "syn", "indel", "filt"], "gene_symbol": ["G1", "G2", "G1", None, "G1"],
        "consequence": ["missense", "missense", "synonymous", None, "missense"], "coding": [True, True, True, False, True],
        "primary_snv": [True, True, True, False, False], "variant_type": ["snv", "snv", "snv", "indel", "snv"],
        "qc_status": ["pass", "pass", "pass", "pass", "excluded"], "qc_exclusion_reason": ["", "", "", "", "fails_source_filter"],
        "in_scope_region": [True] * 5})
    dms = pl.DataFrame({"variant_id": ["m1", "m1", "m2"], "consequence": ["missense"] * 3, "score_set_urn": ["a", "b", "c"]})
    uni = pl.DataFrame({"gene_symbol": ["G1", "G2"], "inclusion_status": ["proposed_include", "proposed_exclude"]})
    flags = package.variant_flags(sites, dms, uni, None)
    assert flags["lost_at"].to_list() == ["", "gene_not_in_proposed_universe", "coding_synonymous", "not_an_snv", "fails_source_filter"]
    assert flags["evaluable_ignoring_avi"].to_list() == [True, False, False, False, False]
    avi = pl.DataFrame({"variant_id": ["m1"], "avi_join_class": ["no_avi_score"], "avi_valid": [False]})
    assert package.variant_flags(sites, dms, uni, avi)["lost_at"][0] == "avi_no_avi_score"

    long = pl.DataFrame({"variant_id": ["m1", "m2", "syn", "indel", "filt"], "sample_set": ["unrelated"] * 5, "level": ["group"] * 5,
                         "label": ["AFR"] * 5, "population_group": ["AFR"] * 5, "n_samples": [10] * 5, "ac": [1, 2, 1, 1, 1],
                         "observed": [True, True, True, True, False]})
    att = package.attrition(long, flags, have_avi=False)
    r = att.row(0, named=True)
    assert (r["n_alt_allele_count_ge_1"], r["n_observed_after_qc"], r["n_primary_snvs"], r["n_coding"], r["n_missense"]) == (5, 4, 3, 3, 2)
    assert (r["n_snvs_with_dms_measurement"], r["n_evaluable_ignoring_avi"], r["n_evaluable"], r["n_genes"], r["n_dms_assays"]) == (2, 1, None, 1, 2)
    counts = pl.DataFrame({"level": ["group"], "population_group": ["AFR"], "original_population": [None], "n_post_qc": ["12"], "n_unrelated": ["10"]})
    fc = package.forecast(att, counts, CFG, have_avi=False).row(0, named=True)
    assert fc["expectation"] == "expected_to_fail" and "200" in fc["criteria_failed"] and "5 genes" in fc["criteria_failed"]
    assert fc["n_unrelated_individuals"] == 10 and fc["valid_spearman"] == package.NOT_EVALUABLE


def test_decisions_tables_and_af_bins():
    text = "# T\n\n## 1. Open review items\n\n| ID | Item |\n|---|---|\n| R-01 | x |\n\n## 2. Decisions\n\ntext\n\n## 3. Known limitations\n\n| ID | L |\n|---|---|\n| L-01 | y |\n| L-02 | z |\n"
    t = package.decisions_tables(text)
    assert t["1. Open review items"] == [["R-01", "x"]] and len(t["3. Known limitations"]) == 2 and "2. Decisions" not in t
    df = pl.DataFrame({"af_global": [0.0, 0.0009, 0.001, 0.0499, 0.05, 1.0, None]}).with_columns(common.af_bin_expr(CFG))
    assert df["af_bin"].to_list() == ["very_rare", "very_rare", "rare", "low_frequency", "common", "common", None]
