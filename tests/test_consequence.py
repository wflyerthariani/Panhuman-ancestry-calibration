import polars as pl

from ancal import consequence as cq


def test_plus_strand_codons():
    #       pos: 1234567890
    seq = "NNATGGCTTGGNN"  # CDS 3..11 on +: ATG GCT TGG = M A W
    rows = cq.transcript_snvs("chr1", "+", [(3, 11)], seq)
    assert len(rows) == 27
    by = {(r[0], r[2]): r for r in rows}
    assert by[(3, "C")][6:] == (1, "M", "L", "start_lost")
    assert by[(8, "A")][6:] == (2, "A", "A", "synonymous")
    assert by[(6, "A")][6:] == (2, "A", "T", "missense")
    assert by[(10, "A")][6:] == (3, "W", "*", "nonsense")
    assert by[(10, "A")][3:6] == (8, "G", "A")


def test_minus_strand_split_cds():
    # Transcript on -: reading the complement right to left gives ATG GCT (M A).
    # Genomic 1-based: segment A = 2..4, segment B = 7..9 (an intron at 5..6).
    seq = "N" + "AGC" + "NN" + "CAT" + "N"
    rows = cq.transcript_snvs("chr1", "-", [(2, 4), (7, 9)], seq)
    by = {(r[0], r[2]): r for r in rows}
    assert by[(9, "A")][3:6] == (1, "A", "T")            # genomic T>A is CDS A>T at CDS position 1
    assert by[(3, "T")][6:] == (2, "A", "D", "missense")   # GCT -> GAT
    assert by[(2, "T")][6:] == (2, "A", "A", "synonymous")


def test_untranslatable_cds_is_skipped():
    assert cq.transcript_snvs("chr1", "+", [(1, 4)], "ATGG") is None
    assert cq.transcript_snvs("chr1", "+", [(1, 3)], "ANG") is None


def test_per_variant_prefers_the_more_severe_consequence():
    base = {"chrom": "chr1", "pos": 5, "ref": "A", "alt": "T", "ensembl_gene_id": "E", "transcript_id": "T", "strand": "+",
            "cds_pos": 1, "cds_ref": "A", "cds_alt": "T", "protein_pos": 2, "ref_aa": "K", "alt_aa": "K"}
    cons = pl.DataFrame([{**base, "gene_symbol": "B", "consequence": "synonymous"},
                         {**base, "gene_symbol": "A", "consequence": "missense", "alt_aa": "M"}]).cast(cq.SCHEMA)
    out = cq.per_variant(cons)
    assert out.height == 1 and out["consequence"][0] == "missense" and out["n_overlapping_genes"][0] == 2
