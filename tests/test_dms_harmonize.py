import polars as pl

from ancal import consequence as cq
from ancal.stages import dms_harmonize as dh

SEQ = "ATGGCTTGGAAAGATCTGCCCGAATTTCACATTAGCACCGTTTATCAGCGTAACTGTGGCGGA"  # M A W K D L P E F H I S T V Y Q R N C G G


def _cons():
    rows = cq.transcript_snvs("chr1", "+", [(1, len(SEQ))], SEQ)
    df = pl.DataFrame(rows, orient="row", schema=["pos", "ref", "alt", "cds_pos", "cds_ref", "cds_alt", "protein_pos",
                                                 "ref_aa", "alt_aa", "consequence"])
    return df.with_columns(pl.lit("chr1").alias("chrom"), pl.lit("G").alias("gene_symbol"), pl.lit("E").alias("ensembl_gene_id"),
                           pl.lit("T").alias("transcript_id"), pl.lit("+").alias("strand")).select(list(cq.SCHEMA)).cast(cq.SCHEMA)


AA1_TO_3 = {v: k for k, v in dh.AA3.items() if k != "*"}
SYN = [(i + 1, AA1_TO_3[aa]) for i, aa in enumerate("MAWKDLPEFHISTVYQRNCGG") if i + 1 >= 6]


def _scores(rows):
    return pl.DataFrame(rows, orient="row", schema=["hgvs_nt", "hgvs_pro", "score"]).cast(pl.String)


def test_parse_protein():
    assert dh.parse_protein("p.Ala2Thr") == (2, "A", "T")
    assert dh.parse_protein("NP_1.2:p.(Trp3Ter)") == (3, "W", "*")
    assert dh.parse_protein("p.Lys4=") == (4, "K", "K")
    assert dh.parse_protein("p.[Ala2Thr;Trp3Cys]") == "multiple_changes"
    assert dh.parse_protein("p.Lys4fs") == "frameshift"
    assert dh.parse_protein(None) == "no_protein_change"


def test_numbering_offset():
    ref_at = {i + 1: aa for i, aa in enumerate("MAWKDLPEFHISTVYQRNCGGACDEFGHIKLMNPQRSTVWY")}
    shifted = {(i + 1, ref_at[i + 11]) for i in range(25)}        # construct numbering starts at residue 11
    assert dh.numbering_offset(shifted, ref_at) == (10, 1.0)
    assert dh.numbering_offset({(p, a) for p, a in ref_at.items()}, ref_at) == (0, 1.0)
    assert dh.numbering_offset({(i + 1, "X") for i in range(25)}, ref_at)[0] == 0   # nothing fits: caller rejects
    assert dh.numbering_offset({(1, ref_at[11]), (2, ref_at[12])}, ref_at)[0] == 0  # too few residues to infer


def test_harmonize_reverses_averages_and_places_on_every_snv():
    df = _scores([
        (None, "p.Ala2Thr", "1.0"), (None, "p.Ala2Thr", "3.0"),   # duplicate -> mean 2.0
        (None, "p.Trp3Cys", "0.0"),                               # two SNVs encode W>C (TGT, TGC)
        (None, "p.Ala2Trp", "5.0"),                               # needs two base changes
        (None, "p.Lys4Arg", "-1.0"), (None, "p.Asp5Glu", "NA"), (None, "p.[Ala2Thr;Trp3Cys]", "9"),
        (None, "p.Gly5Ala", "4.0"),                               # reference residue is D, not G
    ] + [(None, f"p.{aa}{pos}=", "0.5") for pos, aa in SYN])     # enough matching residues to pass the reference check
    out, excl, summ = dh.harmonize_assay(df, "urn:x", "lower_is_more_damaging", _cons())
    assert summ["mapping_method"] == "hgvs_pro"
    assert excl == {"no_score": 1, "multiple_changes": 1, "reference_mismatch": 1, "not_reachable_by_one_snv": 1}
    a2t = out.filter((pl.col("protein_pos") == 2) & (pl.col("alt_aa") == "T"))
    assert a2t["dms_score_original"].to_list() == [2.0] and a2t["dms_score_harmonized"].to_list() == [-2.0]
    assert a2t["n_measurements"].to_list() == [2]
    w3c = out.filter((pl.col("protein_pos") == 3) & (pl.col("consequence") == "missense"))
    assert w3c.height == 2 and w3c["n_snvs_for_change"].to_list() == [2, 2]
    # harmonized scores: A2T -2, W3C 0, A2W -5, K4R 1 -> K4R is the most damaging of four missense changes
    assert out.filter(pl.col("protein_pos") == 4)["dms_percentile"].to_list() == [1.0]
    assert a2t["dms_percentile"].to_list() == [0.5]
    assert out.select("variant_id", "score_set_urn").n_unique() == out.height


def test_assay_with_other_numbering_is_rejected_not_partly_mapped():
    df = _scores([(None, f"p.Gly{i}Ala", "1.0") for i in range(1, 8)])
    out, excl, _ = dh.harmonize_assay(df, "urn:x", "higher_is_more_damaging", _cons())
    assert out.height == 0 and excl == {"assay_numbering_mismatch": 7}


def test_nucleotide_level_assay_maps_by_cds_position():
    df = _scores([("c.4G>A", "p.Ala2Thr", "1.0"), ("c.5C>T", "p.Ala2Val", "2.0"), ("c.6T>C", "p.Ala2=", "0.1"),
                  ("c.6+1G>A", None, "3.0")])
    out, excl, summ = dh.harmonize_assay(df, "urn:x", "higher_is_more_damaging", _cons())
    assert summ["mapping_method"] == "hgvs_nt_cds" and excl == {"not_a_single_coding_snv": 1}
    assert out["variant_id"].to_list() == ["chr1:4:G:A", "chr1:5:C:T", "chr1:6:T:C"]
    assert out["consequence"].to_list() == ["missense", "missense", "synonymous"]
    assert out["dms_percentile"].to_list() == [0.5, 1.0, None]


def test_review_gate_modes():
    cur = pl.DataFrame({
        "score_set_urn": ["a", "b", "c", "d"], "gene_symbol": ["G"] * 4, "in_scope": ["true"] * 4,
        "prescreen_status": ["candidate"] * 4, "proposed_direction": ["lower_is_more_damaging", "unknown", "conflict", "higher_is_more_damaging"],
        "direction_confidence": ["high", "none", "none", "low"], "reviewer_decision": [None, None, "include", "exclude"],
        "reviewer_direction": [None, None, "higher_is_more_damaging", None]})
    prov, unreviewed = dh.assay_decisions(cur, "provisional")
    assert unreviewed == ["a", "b"]
    assert prov["decision"].to_list() == ["include", "exclude", "include", "exclude"]
    assert prov["direction_source"].to_list() == ["provisional_proposal", "provisional_proposal", "reviewer", "reviewer"]
    strict, _ = dh.assay_decisions(cur, "strict")
    assert strict["decision"].to_list() == ["pending", "pending", "include", "exclude"]
