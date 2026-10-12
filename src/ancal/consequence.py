"""Protein consequence of every possible SNV in a MANE Select CDS.

Built from the GRCh38 reference and the MANE Select CDS segments (Stage 4). One row per
(gene, genomic SNV): its CDS position, codon, amino-acid change and consequence class.
This table is the single definition of "coding", "missense", "synonymous" and "nonsense"
used by Stages 6-11, and it is how DMS protein changes are placed on GRCh38.

GENCODE CDS rows exclude the stop codon, so stop-lost variants are outside the table.
"""

from __future__ import annotations

from typing import Callable, Iterable

import polars as pl

_BASES = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {a + b + c: _AA[16 * i + 4 * j + k]
               for i, a in enumerate(_BASES) for j, b in enumerate(_BASES) for k, c in enumerate(_BASES)}
_COMP = str.maketrans("ACGT", "TGCA")
SEVERITY = {"nonsense": 0, "start_lost": 1, "missense": 2, "synonymous": 3}

SCHEMA = {
    "chrom": pl.String, "pos": pl.Int64, "ref": pl.String, "alt": pl.String, "gene_symbol": pl.String,
    "ensembl_gene_id": pl.String, "transcript_id": pl.String, "strand": pl.String, "cds_pos": pl.Int64,
    "cds_ref": pl.String, "cds_alt": pl.String, "protein_pos": pl.Int64, "ref_aa": pl.String,
    "alt_aa": pl.String, "consequence": pl.String,
}


def classify(ref_aa: str, alt_aa: str, protein_pos: int) -> str:
    if alt_aa == ref_aa:
        return "synonymous"
    if alt_aa == "*":
        return "nonsense"
    if protein_pos == 1:
        return "start_lost"
    return "missense"


def transcript_snvs(chrom: str, strand: str, segments: Iterable[tuple[int, int]], seq: str) -> list[tuple] | None:
    """All SNVs of one transcript's CDS. `segments` are 1-based inclusive; `seq` is the contig.

    Returns None if the CDS is not a whole number of codons or contains a non-ACGT base.
    Rows: (pos, ref, alt, cds_pos, cds_ref, cds_alt, protein_pos, ref_aa, alt_aa, consequence).
    """
    positions: list[int] = []
    for start, end in sorted(segments):
        positions.extend(range(start, end + 1))
    bases = "".join(seq[p - 1] for p in positions)
    if strand == "-":
        positions.reverse()
        bases = bases[::-1].translate(_COMP)
    if len(bases) % 3 or set(bases) - set("ACGT"):
        return None
    rows = []
    for c in range(0, len(bases), 3):
        codon = bases[c : c + 3]
        ref_aa = CODON_TABLE[codon]
        ppos = c // 3 + 1
        for k in range(3):
            cds_ref = codon[k]
            g_ref = cds_ref.translate(_COMP) if strand == "-" else cds_ref
            for cds_alt in "ACGT":
                if cds_alt == cds_ref:
                    continue
                alt_aa = CODON_TABLE[codon[:k] + cds_alt + codon[k + 1 :]]
                g_alt = cds_alt.translate(_COMP) if strand == "-" else cds_alt
                rows.append((positions[c + k], g_ref, g_alt, c + k + 1, cds_ref, cds_alt, ppos, ref_aa, alt_aa,
                             classify(ref_aa, alt_aa, ppos)))
    return rows


def build(tx: pl.DataFrame, cds: pl.DataFrame, genes: list[str],
          get_seq: Callable[[str], str]) -> tuple[pl.DataFrame, list[str]]:
    """Consequence table for `genes`, plus the genes skipped (CDS not translatable)."""
    tx = tx.filter(pl.col("gene_name").is_in(genes)).sort("chrom", "tx_start")
    segs = {k[0]: list(zip(g["start"], g["end"])) for k, g in
            cds.filter(pl.col("gene_name").is_in(genes)).group_by("transcript_id")}
    frames, skipped = [], []
    for gene, gene_id, tid, chrom, strand in tx.select("gene_name", "gene_id", "transcript_id", "chrom", "strand").iter_rows():
        rows = transcript_snvs(chrom, strand, segs.get(tid, []), get_seq(chrom)) if tid in segs else None
        if not rows:
            skipped.append(gene)
            continue
        df = pl.DataFrame(rows, orient="row", schema=["pos", "ref", "alt", "cds_pos", "cds_ref", "cds_alt",
                                                     "protein_pos", "ref_aa", "alt_aa", "consequence"])
        frames.append(df.with_columns(
            pl.lit(chrom).alias("chrom"), pl.lit(gene).alias("gene_symbol"), pl.lit(gene_id).alias("ensembl_gene_id"),
            pl.lit(tid).alias("transcript_id"), pl.lit(strand).alias("strand")).select(list(SCHEMA)))
    out = pl.concat(frames) if frames else pl.DataFrame(schema=SCHEMA)
    return out.cast(SCHEMA), sorted(skipped)


def per_variant(cons: pl.DataFrame) -> pl.DataFrame:
    """One row per genomic SNV. Where MANE CDS of two genes overlap, the more severe consequence wins."""
    return (
        cons.with_columns(pl.col("consequence").replace_strict(SEVERITY, return_dtype=pl.Int32).alias("_sev"))
        .sort("chrom", "pos", "ref", "alt", "_sev", "gene_symbol")
        .group_by("chrom", "pos", "ref", "alt", maintain_order=True)
        .agg(pl.col("gene_symbol").first(), pl.col("ensembl_gene_id").first(), pl.col("transcript_id").first(),
             pl.col("protein_pos").first(), pl.col("ref_aa").first(), pl.col("alt_aa").first(),
             pl.col("consequence").first(), pl.col("gene_symbol").n_unique().alias("n_overlapping_genes"))
    )
