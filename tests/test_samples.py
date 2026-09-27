import gzip
import struct

import polars as pl
import pytest
import zstandard

from ancal.sources import hail_table
from ancal.stages import samples as st


def _meta(rows):
    cols = list(st.META)
    return pl.DataFrame([dict(zip(cols, r)) for r in rows], schema={c: pl.String for c in cols}).with_columns(
        pl.col("latitude").cast(pl.Float64), pl.col("longitude").cast(pl.Float64)
    )


def _row(sid, pop, group, project="HGDP", hard="false", hard_reason=None, gnomad_label=None):
    return (sid, project, pop, gnomad_label or pop, group, "region", "1.0", "2.0", "XX", None, hard, hard_reason)


def _pcs(ids_coords):
    data = {"sample_id": [i for i, _ in ids_coords]}
    for k in range(1, st.N_PCS + 1):
        data[f"pc{k}"] = [float(c[0] if k == 1 else c[1] if k == 2 else 0.0) for _, c in ids_coords]
    return pl.DataFrame(data)


@pytest.fixture
def table():
    meta = _meta([
        _row("A1", "Yoruba", "AFR"), _row("A2", "Yoruba", "AFR"), _row("A3", "ACB", "AFR", "1000 Genomes"),
        _row("E1", "French", "EUR"), _row("E2", "French", "EUR"),
        _row("X1", "French", "EUR", hard="true", hard_reason="[sex_aneuploidy]"),
        _row("O1", "French", "EUR"), _row("Q1", "Yoruba", "AFR"),
    ])
    pcs = _pcs([("A1", (0, 0)), ("A2", (0.1, 0)), ("A3", (5, 5)), ("E1", (5, 5)), ("E2", (5.1, 5))])
    vcf = ["A1", "A2", "A3", "E1", "E2", "X1", "O1", "Q1", "CHMI"]
    return st.build_sample_table(meta, vcf, pcs, related={"A2", "X1"}, outliers={"O1"})


def test_one_row_per_vcf_sample(table):
    assert table.height == 9
    assert table["sample_id"].n_unique() == 9


def test_inclusion_and_exclusion(table):
    r = {row["sample_id"]: row for row in table.iter_rows(named=True)}
    assert r["A1"]["inclusion_status"] == "primary" and r["A1"]["relatedness_status"] == "unrelated"
    assert r["A2"]["inclusion_status"] == "sensitivity_only" and r["A2"]["relatedness_status"] == "related"
    assert r["X1"]["exclusion_reason"] == "gnomad_hard_filter:sex_aneuploidy"
    assert r["X1"]["relatedness_status"] == "not_assessed"  # relatedness only assessed post-QC
    assert r["O1"]["exclusion_reason"] == "pca_outlier"
    assert r["Q1"]["exclusion_reason"] == "not_in_koenig_post_qc_release"
    assert r["CHMI"]["exclusion_reason"] == "not_in_sample_metadata"
    assert all(r[s]["exclusion_reason"] is None for s in ("A1", "A2", "A3", "E1", "E2"))


def test_pc_discordance_flags_but_keeps(table):
    r = {row["sample_id"]: row for row in table.iter_rows(named=True)}
    # A3 is labelled AFR but sits on the EUR centroid
    assert r["A3"]["pc_group_discordant"] is True and r["A3"]["pc_nearest_group"] == "EUR"
    assert r["A3"]["inclusion_status"] == "primary"
    assert r["A1"]["pc_group_discordant"] is False
    assert r["X1"]["pc_group_discordant"] is None  # not scored: excluded, no PCs


def test_admixed_flag_and_mapping(table):
    m = st.population_mapping(table)
    acb = m.filter(pl.col("original_population") == "ACB").row(0, named=True)
    assert acb["recently_admixed"] is True and acb["n_pc_group_discordant"] == 1
    assert not m["ambiguous_group_label"].any()


def test_counts_levels_and_koenig_diff(table):
    koenig = pl.DataFrame({
        "population_group": ["AFR", "AFR", "EUR"], "original_population": ["Yoruba", "ACB", "French"],
        "koenig_n_post_qc": [2, 1, 2], "koenig_n_unrelated": [1, 1, 1],
    })
    c = st.population_counts(table, koenig)
    total = c.filter(pl.col("level") == "total").row(0, named=True)
    assert total["n_post_qc"] == 5 and total["n_unrelated"] == 4 and total["n_related"] == 1
    assert total["diff_post_qc_vs_koenig"] == 0
    fr = c.filter(pl.col("original_population") == "French").row(0, named=True)
    assert fr["n_post_qc"] == 2 and fr["n_excluded"] == 2 and fr["diff_unrelated_vs_koenig_summary"] == 1
    assert c["level"].to_list()[0] == "total"


# ---------------------------------------------------------------- Hail decoder


def _hail_block(payload: bytes) -> bytes:
    frame = zstandard.ZstdCompressor().compress(payload)
    return struct.pack("<ii", 4 + len(frame), len(payload)) + frame


def _encode_rows(values):
    out = bytearray()
    for v in values:
        out += b"\x01\x00" + bytes([len(v)]) + v.encode()
    return bytes(out + b"\x00")


def test_hail_decoder_round_trip():
    raw = _hail_block(_encode_rows(["HG00116", "NA12878"]))
    assert hail_table.decode_string_rows(raw) == ["HG00116", "NA12878"]


def test_hail_decoder_multiple_blocks():
    payload = _encode_rows(["S1", "S2", "S3"])
    raw = _hail_block(payload[:5]) + _hail_block(payload[5:])
    assert hail_table.decode_string_rows(raw) == ["S1", "S2", "S3"]


def test_hail_decoder_rejects_garbage():
    with pytest.raises(hail_table.UnsupportedHailTable):
        hail_table.decode_string_rows(_hail_block(b"\x07junk"))


def test_parse_vcf_header(tmp_path):
    h = tmp_path / "h.txt"
    h.write_text("##fileformat=VCFv4.2\n##contig=<ID=chr17,length=83257441>\n"
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n")
    samples, contigs, _ = st.parse_vcf_header(h)
    assert samples == ["S1", "S2"] and contigs["chr17"] == 83257441
