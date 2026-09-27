import polars as pl

from ancal.stages import reference as ref


def _gtf(chrom, feature, tags='tag "MANE_Select";'):
    return f'{chrom}\tHAVANA\t{feature}\t100\t200\t.\t+\t0\tgene_id "G"; {tags}\n'


def test_gtf_keep():
    assert ref.gtf_keep(_gtf("chr17", "CDS"))
    assert ref.gtf_keep(_gtf("chrX", "transcript"))
    assert not ref.gtf_keep(_gtf("chr17", "exon"))
    assert not ref.gtf_keep(_gtf("chr17", "CDS", tags='tag "basic";'))
    assert not ref.gtf_keep(_gtf("chrY", "CDS"))  # chrY / PAR copies excluded (chromosomes = chr1-22, X)
    assert not ref.gtf_keep("##description: header\n")


def _cds(rows):
    return pl.DataFrame(rows, schema=["chrom", "start", "end", "gene_name"], orient="row")


def test_regions_bed_merges_and_reports_missing():
    cds = _cds([("chr17", 100, 200, "A"), ("chr17", 150, 300, "B"), ("chr17", 500, 600, "A"), ("chr2", 10, 20, "C")])
    bed, missing = ref.regions_bed(cds, ["A", "B", "C", "NOPE"])
    assert missing == ["NOPE"]
    assert bed.rows() == [("chr2", 9, 20, "C"), ("chr17", 99, 300, "A,B"), ("chr17", 499, 600, "A")]


def test_interval_index_is_half_open_zero_based():
    bed = pl.DataFrame([("chr17", 99, 200, "A")], schema=["chrom", "start", "end", "name"], orient="row")
    idx = ref.IntervalIndex(bed)
    assert not idx.contains("chr17", 99)
    assert idx.contains("chr17", 100) and idx.contains("chr17", 200)
    assert not idx.contains("chr17", 201)
    assert not idx.contains("chr1", 150)
