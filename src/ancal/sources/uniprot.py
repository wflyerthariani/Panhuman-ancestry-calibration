"""UniProt accession -> primary gene name, snapshotted like the MaveDB client."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
from pathlib import Path

from .. import provenance
from ..settings import Settings
from .mavedb import _request

BATCH = 200


def gene_names(settings: Settings, base: str, accessions: set[str]) -> dict[str, str]:
    """Map each accession (isoform suffix removed) to its primary gene name. Cached in dms/uniprot_genes.json.gz."""
    path = settings.data_dir("dms") / "uniprot_genes.json.gz"
    cache: dict[str, str] = json.loads(gzip.open(path).read()) if path.exists() else {}
    wanted = sorted({a.split("-")[0].strip() for a in accessions if a} - set(cache))
    release = ""
    for i in range(0, len(wanted), BATCH):
        chunk = wanted[i : i + BATCH]
        r = _request("GET", base, params={"accessions": ",".join(chunk), "fields": "accession,gene_primary",
                                          "format": "tsv"})
        release = r.headers.get("X-UniProt-Release", release)
        for row in csv.DictReader(io.StringIO(r.text), delimiter="\t"):
            cache[row["Entry"]] = row["Gene Names (primary)"].split(";")[0].strip()
        for acc in chunk:
            cache.setdefault(acc, "")  # remember misses so they are not re-queried
    if wanted:
        data = json.dumps(cache, sort_keys=True, separators=(",", ":")).encode()
        with gzip.GzipFile(path, "wb", mtime=0) as fh:
            fh.write(data)
        provenance.record(
            settings, path, stage="stage5", source_name="UniProt accession -> gene", command=f"GET {base}",
            software_name="ancal.sources.uniprot", software_version=provenance.code_version(settings.root),
            input_uri=base, access_method="api", database_version=f"UniProt {release}",
            input_sha256=hashlib.sha256(data).hexdigest(),
        )
    return cache
