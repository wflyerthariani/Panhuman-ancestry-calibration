"""MaveDB API client with snapshotting.

MaveDB is a live database: records are added and revised. Every response used by the
pipeline is saved under data_store/dms/mavedb/ and recorded in the manifest, and later
runs reuse that snapshot unless refresh=True. The inventory is therefore reproducible, and
re-running Stage 5 does not silently pick up new MaveDB records.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import requests

from .. import provenance
from ..settings import Settings

PAGE = 100
WORKERS = 4
RETRIES = 5


def _request(method: str, url: str, **kw) -> requests.Response:
    for attempt in range(RETRIES):
        try:
            r = requests.request(method, url, timeout=120, **kw)
            if r.status_code in (429, 500, 502, 503, 504):
                raise requests.HTTPError(f"{r.status_code} from {url}")
            r.raise_for_status()
            return r
        except (requests.ConnectionError, requests.Timeout, requests.HTTPError):
            if attempt == RETRIES - 1:
                raise
            time.sleep(2 ** attempt)
    raise AssertionError("unreachable")


class MaveDB:
    def __init__(self, settings: Settings, base: str, refresh: bool = False):
        self.s = settings
        self.base = base.rstrip("/")
        self.dir = settings.data_dir("dms") / "mavedb"
        self.dir.mkdir(parents=True, exist_ok=True)
        self.refresh = refresh

    # ------------------------------------------------------------ snapshots

    def _save_json(self, name: str, obj, source_uri: str, source_name: str) -> Path:
        path = self.dir / name
        data = json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()
        with gzip.GzipFile(path, "wb", mtime=0) as fh:  # mtime=0 -> byte-identical for identical content
            fh.write(data)
        provenance.record(
            self.s, path, stage="stage5", source_name=source_name, command=f"GET/POST {source_uri}",
            software_name="ancal.sources.mavedb", software_version=provenance.code_version(self.s.root),
            input_uri=source_uri, access_method="api", database_version=f"MaveDB API {self.api_version()}",
            input_sha256=hashlib.sha256(data).hexdigest(),
        )
        return path

    def _load_json(self, name: str):
        path = self.dir / name
        if path.exists() and not self.refresh:
            with gzip.open(path, "rb") as fh:
                return json.loads(fh.read())
        return None

    _version: str | None = None

    def api_version(self) -> str:
        if MaveDB._version is None:
            MaveDB._version = _request("GET", f"{self.base}/api/version").json().get("version", "?")
        return MaveDB._version

    # ------------------------------------------------------------- queries

    def human_score_sets(self) -> list[dict]:
        """Full records of every published human score set (snapshot: score_sets.json.gz)."""
        cached = self._load_json("score_sets.json.gz")
        if cached is not None:
            return cached
        url = f"{self.base}/score-sets/search"
        body = {"targetOrganismNames": ["Homo sapiens"], "published": True, "limit": PAGE}
        short, offset, total = [], 0, None
        while total is None or offset < total:
            page = _request("POST", url, json={**body, "offset": offset}).json()
            total = page["numScoreSets"]
            short += page["scoreSets"]
            offset += PAGE
        urns = sorted({s["urn"] for s in short})
        with ThreadPoolExecutor(WORKERS) as pool:
            full = list(pool.map(lambda u: _request("GET", f"{self.base}/score-sets/{u}").json(), urns))
        if len(full) != total:
            raise RuntimeError(f"MaveDB search reported {total} score sets but {len(full)} were retrieved")
        self._save_json("score_sets.json.gz", full, f"{url} {json.dumps(body)} + GET /score-sets/{{urn}}",
                        "MaveDB human score sets (full records)")
        return full

    def mapped_genes(self) -> dict[str, list[str]]:
        cached = self._load_json("mapped_genes.json.gz")
        if cached is not None:
            return cached
        url = f"{self.base}/score-sets/mapped-genes"
        data = _request("GET", url).json()
        self._save_json("mapped_genes.json.gz", data, url, "MaveDB mapped HGNC genes per score set")
        return data

    def scores_csv(self, urn: str) -> Path:
        """Score CSV for one score set (snapshot: scores/<urn>.csv.gz)."""
        path = self.dir / "scores" / f"{urn.replace(':', '_')}.csv.gz"
        url = f"{self.base}/score-sets/{urn}/scores"
        if path.exists() and not self.refresh and provenance.recorded_row(self.s, path, url):
            return path
        path.parent.mkdir(parents=True, exist_ok=True)
        text = _request("GET", url).content
        with gzip.GzipFile(path, "wb", mtime=0) as fh:
            fh.write(text)
        provenance.record(
            self.s, path, stage="stage5", source_name=f"MaveDB scores {urn}", command=f"GET {url}",
            software_name="ancal.sources.mavedb", software_version=provenance.code_version(self.s.root),
            input_uri=url, access_method="api", database_version=f"MaveDB API {self.api_version()}",
            input_sha256=hashlib.sha256(text).hexdigest(),
        )
        return path
