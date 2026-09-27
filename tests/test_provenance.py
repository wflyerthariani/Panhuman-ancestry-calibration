import gzip
import hashlib
import threading
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import pytest
import yaml

from ancal import budget, provenance
from ancal.settings import load_settings


class _QuietHandler(SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


@pytest.fixture
def server(tmp_path: Path):
    served = tmp_path / "served"
    served.mkdir()
    httpd = ThreadingHTTPServer(("127.0.0.1", 0), partial(_QuietHandler, directory=str(served)))
    t = threading.Thread(target=httpd.serve_forever, daemon=True)
    t.start()
    yield served, f"http://127.0.0.1:{httpd.server_address[1]}"
    httpd.shutdown()


def _settings(tmp_path: Path, budget_gb: float = 1.0):
    repo = tmp_path / "repo"
    (repo / "config").mkdir(parents=True)
    cfg = {"scope": "smoke", "storage": {"root": "./data_store", "local_budget_gb": budget_gb,
                                          "counted_output_dirs": ["reports"]}}
    (repo / "config" / "data_sources.yaml").write_text(yaml.safe_dump(cfg))
    (repo / "data_store").mkdir()
    return load_settings(repo)


def test_fetch_records_complete_row(tmp_path, server):
    served, base = server
    payload = b"sample\tpop\nA\tYRI\n"
    (served / "meta.tsv").write_bytes(payload)
    s = _settings(tmp_path)
    dest = s.data_dir("population_metadata") / "meta.tsv"

    row = provenance.fetch(s, f"{base}/meta.tsv", dest, stage="stageX", source_name="test meta")

    assert dest.read_bytes() == payload
    assert row.input_sha256 == row.output_sha256 == hashlib.sha256(payload).hexdigest()
    assert row.output_path == "population_metadata/meta.tsv"
    assert row.input_bytes == str(len(payload))
    assert provenance.incomplete_rows(s) == []


def test_fetch_is_idempotent(tmp_path, server):
    served, base = server
    (served / "a.txt").write_bytes(b"hello\n")
    s = _settings(tmp_path)
    dest = s.data_dir("x") / "a.txt"
    provenance.fetch(s, f"{base}/a.txt", dest, stage="s", source_name="a")
    provenance.fetch(s, f"{base}/a.txt", dest, stage="s", source_name="a")
    assert len(provenance.read_manifest(s)) == 1
    dest.write_bytes(b"tampered\n")  # no longer matches the manifest -> refetched
    provenance.fetch(s, f"{base}/a.txt", dest, stage="s", source_name="a")
    assert dest.read_bytes() == b"hello\n"
    assert len(provenance.read_manifest(s)) == 2


def test_fetch_refuses_over_budget(tmp_path, server):
    served, base = server
    (served / "big.bin").write_bytes(b"x" * 5000)
    s = _settings(tmp_path, budget_gb=1e-6)  # ~1 KB
    dest = s.data_dir("x") / "big.bin"
    with pytest.raises(budget.BudgetExceeded, match="would exceed"):
        provenance.fetch(s, f"{base}/big.bin", dest, stage="s", source_name="big")
    assert not dest.exists()
    assert not dest.with_name("big.bin.part").exists()
    assert provenance.read_manifest(s) == []


def test_stream_filter_multimember_gzip(tmp_path, server):
    served, base = server
    # bgzip-style: two concatenated gzip members
    raw = gzip.compress(b"#header\nchr17\t1\tkeep\n") + gzip.compress(b"chr1\t2\tdrop\nchr17\t3\tkeep\n")
    (served / "table.tsv.gz").write_bytes(raw)
    s = _settings(tmp_path)
    dest = s.data_dir("x") / "subset.tsv"

    row = provenance.stream_filter(
        s, f"{base}/table.tsv.gz", dest, lambda l: l.startswith(("#", "chr17\t")),
        stage="s", source_name="tbl", filter_description="chr17 rows",
    )

    assert dest.read_text() == "#header\nchr17\t1\tkeep\nchr17\t3\tkeep\n"
    assert row.input_sha256 == hashlib.sha256(raw).hexdigest()  # hash of the full source
    assert row.output_sha256 != row.input_sha256
    assert "chr17 rows" in row.command


def test_stream_filter_aborts_when_output_outgrows_budget(tmp_path, server, monkeypatch):
    served, base = server
    (served / "t.tsv").write_text("".join(f"row{i}\t{'x' * 50}\n" for i in range(200)))
    s = _settings(tmp_path, budget_gb=2e-6)  # ~2 KB
    monkeypatch.setattr(provenance, "BUDGET_CHECK_EVERY", 100)
    dest = s.data_dir("x") / "all.tsv"
    with pytest.raises(budget.BudgetExceeded):
        provenance.stream_filter(s, f"{base}/t.tsv", dest, lambda l: True, stage="s", source_name="t",
                                 filter_description="everything")
    assert not dest.exists()
    assert not dest.with_name("all.tsv.part").exists()


def test_record_remote_tabix_row(tmp_path):
    s = _settings(tmp_path)
    out = s.data_dir("variants") / "sites.vcf.gz"
    out.parent.mkdir(parents=True)
    out.write_bytes(b"x")
    row = provenance.record(s, out, stage="s", source_name="gnomad", command="bcftools view -R r.bed URL",
                            software_name="bcftools", software_version="1.21", input_uri="https://x/y.vcf.bgz",
                            access_method="remote_tabix")
    assert row.input_sha256 == provenance.NOT_HASHED
    assert provenance.incomplete_rows(s) == []


def test_incomplete_rows_detected(tmp_path):
    s = _settings(tmp_path)
    p = provenance.manifest_path(s)
    p.parent.mkdir(parents=True)
    cols = [f.name for f in provenance.fields(provenance.ManifestRow)]
    p.write_text("\t".join(cols) + "\n" + "\t".join("" for _ in cols) + "\n")
    [(line, missing)] = provenance.incomplete_rows(s)
    assert line == 2 and "output_sha256" in missing


def test_software_manifest(tmp_path):
    s = _settings(tmp_path)
    path = provenance.write_software_manifest(s)
    text = path.read_text()
    assert text.startswith("kind\tname\tversion\n")
    assert "\tbcftools\t" in text and "\tpolars\t" in text


def test_stream_filter_transform_bgzip(tmp_path, server):
    import subprocess

    served, base = server
    raw = gzip.compress(b">chr17\nacgtNN\nACgt\n")
    (served / "chr17.fa.gz").write_bytes(raw)
    s = _settings(tmp_path)
    dest = s.data_dir("reference") / "chr17.fa.gz"
    provenance.stream_filter(s, f"{base}/chr17.fa.gz", dest,
                             transform=lambda l: l if l.startswith(">") else l.upper(),
                             stage="s", source_name="fa", filter_description="upper", bgzip_output=True)
    assert gzip.decompress(dest.read_bytes()) == b">chr17\nACGTNN\nACGT\n"
    subprocess.run(["samtools", "faidx", str(dest)], check=True)  # only works on true bgzip
    assert Path(f"{dest}.fai").read_text().startswith("chr17\t10\t")


def test_stream_filter_needs_exactly_one_of_keep_transform(tmp_path):
    s = _settings(tmp_path)
    with pytest.raises(ValueError):
        provenance.stream_filter(s, "x", tmp_path / "o", stage="s", source_name="x", filter_description="x")
