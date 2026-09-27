"""Provenance: every fetched or derived data file gets a manifest row (protocol Step 1).

Three ways data enters the store, each recorded with the same columns:
  fetch()          - download a whole (small) file
  stream_filter()  - stream a large file, keep only matching lines; the full source
                     is hashed in flight but never written to disk
  record()         - a file produced by a command (e.g. a bcftools remote-tabix
                     region query), with the command line and input identifiers

Manifests live in reports/manifests/ (tracked in git), not in the data store.
"""

from __future__ import annotations

import csv
import getpass
import gzip
import hashlib
import importlib.metadata
import io
import os
import platform
import shutil
import subprocess
import sys
import zlib
from dataclasses import asdict, dataclass, fields
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterator

import requests

from . import __version__, budget
from .settings import Settings

CHUNK = 1 << 20  # 1 MiB
BUDGET_CHECK_EVERY = 64 * CHUNK
NOT_HASHED = "not_hashed:remote_region_query"


@dataclass
class ManifestRow:
    recorded_at: str
    stage: str
    scope: str
    source_name: str
    database_version: str
    reference_build: str
    access_method: str            # download | stream_filter | remote_tabix | derived | api
    input_uri: str
    input_etag: str
    input_last_modified: str
    input_bytes: str
    input_sha256: str
    output_path: str              # relative to storage_root (or repo, for repo outputs)
    output_bytes: int
    output_sha256: str
    software_name: str
    software_version: str
    command: str
    code_version: str             # ancal version + git commit
    storage_root: str
    operator: str


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(CHUNK), b""):
            h.update(block)
    return h.hexdigest()


def git_commit(root: Path) -> str:
    try:
        out = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "--short", "HEAD"], capture_output=True, text=True, check=True
        ).stdout.strip()
        dirty = subprocess.run(
            ["git", "-C", str(root), "status", "--porcelain", "--untracked-files=no", "--", "src", "workflow", "config"],
            capture_output=True, text=True,
        ).stdout.strip()
        return out + ("-dirty" if dirty else "")
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


def code_version(root: Path) -> str:
    return f"ancal {__version__} @ {git_commit(root)}"


def manifest_path(settings: Settings) -> Path:
    return settings.root / "reports" / "manifests" / "data_manifest.tsv"


def _append(settings: Settings, row: ManifestRow) -> None:
    path = manifest_path(settings)
    path.parent.mkdir(parents=True, exist_ok=True)
    new = not path.exists()
    with open(path, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[f.name for f in fields(ManifestRow)], delimiter="\t", lineterminator="\n")
        if new:
            w.writeheader()
        w.writerow(asdict(row))


def read_manifest(settings: Settings) -> list[dict]:
    path = manifest_path(settings)
    if not path.exists():
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _rel(settings: Settings, path: Path) -> str:
    path = path.resolve()
    for base in (settings.storage_root, settings.root):
        try:
            return str(path.relative_to(base))
        except ValueError:
            continue
    return str(path)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


@dataclass(frozen=True)
class RemoteInfo:
    bytes: int | None
    etag: str
    last_modified: str


def head(uri: str, timeout: float = 30) -> RemoteInfo:
    """Size and version identifiers of a remote file, without downloading it."""
    if uri.startswith("file://") or "://" not in uri:
        p = Path(uri.removeprefix("file://"))
        st = p.stat()
        return RemoteInfo(st.st_size, "", datetime.fromtimestamp(st.st_mtime, timezone.utc).isoformat())
    r = requests.head(uri, allow_redirects=True, timeout=timeout)
    r.raise_for_status()
    n = r.headers.get("Content-Length")
    return RemoteInfo(int(n) if n else None, r.headers.get("ETag", "").strip('"'), r.headers.get("Last-Modified", ""))


def _open_stream(uri: str, timeout: float = 60) -> Iterator[bytes]:
    if uri.startswith("file://") or "://" not in uri:
        with open(uri.removeprefix("file://"), "rb") as fh:
            yield from iter(lambda: fh.read(CHUNK), b"")
        return
    with requests.get(uri, stream=True, timeout=timeout) as r:
        r.raise_for_status()
        yield from r.iter_content(CHUNK)


def _row(settings, *, stage, source_name, database_version, reference_build, access_method, uri, info,
         input_sha256, out, software_name, software_version, command) -> ManifestRow:
    return ManifestRow(
        recorded_at=_now(), stage=stage, scope=settings.scope, source_name=source_name,
        database_version=database_version, reference_build=reference_build, access_method=access_method,
        input_uri=uri, input_etag=info.etag if info else "", input_last_modified=info.last_modified if info else "",
        input_bytes="" if not info or info.bytes is None else str(info.bytes), input_sha256=input_sha256,
        output_path=_rel(settings, out), output_bytes=out.stat().st_size, output_sha256=sha256_file(out),
        software_name=software_name, software_version=software_version, command=command,
        code_version=code_version(settings.root), storage_root=str(settings.storage_root),
        operator=getpass.getuser(),
    )


def recorded_row(settings: Settings, dest: Path, uri: str) -> dict | None:
    """The latest manifest row for `dest` from `uri`, if the file on disk still matches it."""
    if not dest.exists():
        return None
    rel = _rel(settings, dest)
    rows = [r for r in read_manifest(settings) if r["output_path"] == rel and r["input_uri"] == uri]
    if rows and rows[-1]["output_sha256"] == sha256_file(dest):
        return rows[-1]
    return None


def fetch(
    settings: Settings, uri: str, dest: Path, *, stage: str, source_name: str,
    database_version: str = "", reference_build: str = "", force: bool = False,
) -> ManifestRow | dict:
    """Download a whole file. Refuses if the file would push usage over the budget.

    Idempotent: if `dest` already matches its manifest row, nothing is downloaded or recorded.
    """
    if not force and (existing := recorded_row(settings, dest, uri)):
        return existing
    info = head(uri)
    if info.bytes is not None:
        budget.ensure_room(settings, info.bytes, f"{source_name} ({uri})")
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".part")
    h = hashlib.sha256()
    written = 0
    try:
        with open(tmp, "wb") as out:
            for block in _open_stream(uri):
                h.update(block)
                out.write(block)
                written += len(block)
                if info.bytes is None and written % BUDGET_CHECK_EVERY < CHUNK:
                    budget.check(settings)
        tmp.replace(dest)
    finally:
        tmp.unlink(missing_ok=True)
    row = _row(
        settings, stage=stage, source_name=source_name, database_version=database_version,
        reference_build=reference_build, access_method="download", uri=uri, info=info,
        input_sha256=h.hexdigest(), out=dest, software_name="ancal.provenance.fetch",
        software_version=__version__, command=f"fetch {uri} -> {dest.name}",
    )
    _append(settings, row)
    budget.check(settings)
    return row


def _gunzip_multi(blocks: Iterator[bytes], on_raw: Callable[[bytes], None]) -> Iterator[bytes]:
    """Decompress a gzip stream that may contain several members (e.g. bgzip)."""
    d = zlib.decompressobj(wbits=31)
    for block in blocks:
        on_raw(block)
        data = block
        while data:
            yield d.decompress(data)
            if d.eof:
                data = d.unused_data
                d = zlib.decompressobj(wbits=31)
            else:
                data = b""
    yield d.flush()


def _lines(chunks: Iterator[bytes]) -> Iterator[str]:
    buf = b""
    for chunk in chunks:
        buf += chunk
        *complete, buf = buf.split(b"\n")
        for line in complete:
            yield line.decode() + "\n"
    if buf:
        yield buf.decode() + "\n"


def stream_filter(
    settings: Settings, uri: str, dest: Path, keep: Callable[[str], bool], *, stage: str, source_name: str,
    database_version: str = "", reference_build: str = "", filter_description: str,
    gzipped_input: bool | None = None, gzip_output: bool | None = None,
) -> ManifestRow:
    """Stream `uri`, write only lines where keep(line) is true. The source is hashed, never stored.

    The output is budget-checked as it grows, so a filter that keeps too much aborts
    instead of filling the disk.
    """
    info = head(uri)
    if gzipped_input is None:
        gzipped_input = uri.endswith((".gz", ".bgz"))
    if gzip_output is None:
        gzip_output = dest.name.endswith(".gz")
    h = hashlib.sha256()

    raw = _open_stream(uri)
    if gzipped_input:
        chunks = _gunzip_multi(raw, h.update)
    else:
        def _hashed(blocks):
            for b in blocks:
                h.update(b)
                yield b
        chunks = _hashed(raw)

    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".part")
    kept = 0
    try:
        opener = (lambda p: gzip.open(p, "wt", compresslevel=6)) if gzip_output else (lambda p: open(p, "w"))
        with opener(tmp) as out:
            for line in _lines(chunks):
                if keep(line):
                    out.write(line)
                    kept += len(line)
                    if kept % BUDGET_CHECK_EVERY < len(line):
                        budget.check(settings)
        tmp.replace(dest)
    finally:
        tmp.unlink(missing_ok=True)
    row = _row(
        settings, stage=stage, source_name=source_name, database_version=database_version,
        reference_build=reference_build, access_method="stream_filter", uri=uri, info=info,
        input_sha256=h.hexdigest(), out=dest, software_name="ancal.provenance.stream_filter",
        software_version=__version__, command=f"stream_filter [{filter_description}] {uri} -> {dest.name}",
    )
    _append(settings, row)
    budget.check(settings)
    return row


def record(
    settings: Settings, out: Path, *, stage: str, source_name: str, command: str, software_name: str,
    software_version: str, input_uri: str = "", access_method: str = "derived", database_version: str = "",
    reference_build: str = "", input_info: RemoteInfo | None = None, input_sha256: str = "",
) -> ManifestRow | dict:
    """Record a file produced by a command (remote tabix query, derived table, ...).

    Idempotent: nothing is appended if the latest row for this output has the same
    checksum and command.
    """
    rel = _rel(settings, out)
    prior = [r for r in read_manifest(settings) if r["output_path"] == rel]
    if prior and prior[-1]["output_sha256"] == sha256_file(out) and prior[-1]["command"] == command:
        return prior[-1]
    row = _row(
        settings, stage=stage, source_name=source_name, database_version=database_version,
        reference_build=reference_build, access_method=access_method, uri=input_uri, info=input_info,
        input_sha256=input_sha256 or (NOT_HASHED if access_method == "remote_tabix" else ""),
        out=out, software_name=software_name, software_version=software_version, command=command,
    )
    _append(settings, row)
    budget.check(settings)
    return row


# ------------------------------------------------------------ software manifest

TOOLS = ("bcftools", "tabix", "bgzip", "samtools", "snakemake", "conda", "git")
PY_PACKAGES = ("polars", "pyarrow", "pandas", "pydantic", "PyYAML", "requests", "fsspec", "pytest", "snakemake")


def _tool_version(tool: str) -> str:
    exe = shutil.which(tool) or (os.environ.get("CONDA_EXE") if tool == "conda" else None)
    if not exe:
        return "not found"
    out = subprocess.run([exe, "--version"], capture_output=True, text=True)
    text = (out.stdout or out.stderr).strip()
    return text.splitlines()[0] if text else "?"


def write_software_manifest(settings: Settings) -> Path:
    rows = [("platform", "os", platform.platform()), ("platform", "machine", platform.machine()),
            ("platform", "python", sys.version.split()[0]), ("platform", "env_prefix", sys.prefix),
            ("code", "ancal", code_version(settings.root))]
    rows += [("tool", t, _tool_version(t)) for t in TOOLS]
    for pkg in PY_PACKAGES:
        try:
            rows.append(("python", pkg, importlib.metadata.version(pkg)))
        except importlib.metadata.PackageNotFoundError:
            rows.append(("python", pkg, "not installed"))
    path = settings.root / "reports" / "manifests" / "software_manifest.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    buf = io.StringIO()
    w = csv.writer(buf, delimiter="\t", lineterminator="\n")
    w.writerow(["kind", "name", "version"])
    w.writerows(rows)
    path.write_text(buf.getvalue())
    return path


REQUIRED_ROW_FIELDS = (
    "recorded_at", "stage", "source_name", "access_method", "input_uri", "output_path",
    "output_bytes", "output_sha256", "software_name", "software_version", "command", "code_version",
)


def incomplete_rows(settings: Settings) -> list[tuple[int, list[str]]]:
    """Manifest rows missing any required Step 1 field."""
    bad = []
    for i, r in enumerate(read_manifest(settings), start=2):
        missing = [f for f in REQUIRED_ROW_FIELDS if not r.get(f)]
        if r.get("access_method") in ("download", "stream_filter") and not r.get("input_sha256"):
            missing.append("input_sha256")
        if missing:
            bad.append((i, missing))
    return bad
