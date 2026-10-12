"""AVI (AlphaGenome Atlas, Aggregate Variant Impact) score loader.

The Atlas distribution format was not known when this was written, so the loader is driven
entirely by the `avi:` block of config/data_sources.yaml and handles the common cases:

  format: parquet   a file, directory or glob; filtered lazily to the scope positions
  format: tsv|csv   plain or gzipped; region access through tabix when <uri><index_suffix>
                    exists, otherwise the file is streamed once and filtered
  format: vcf       bgzipped and indexed; the score is an INFO tag

Only records inside the scope regions are stored. Nothing here interprets the score.
"""

from __future__ import annotations

import gzip
import subprocess
from pathlib import Path

import polars as pl

from .. import config, provenance
from ..settings import Settings

CANDIDATES = {
    "chrom": ("chrom", "#chrom", "chr", "chromosome", "contig"),
    "pos": ("pos", "position", "start", "pos_grch38"),
    "ref": ("ref", "reference", "ref_allele", "reference_allele"),
    "alt": ("alt", "alternate", "alt_allele", "alternate_allele"),
    "transcript": ("transcript", "transcript_id", "enst", "feature"),
}
RAW_SCHEMA = {"chrom_raw": pl.String, "pos": pl.Int64, "ref_raw": pl.String, "alt_raw": pl.String,
              "avi_score": pl.Float64, "transcript": pl.String}


class AviConfigError(ValueError):
    pass


def resolve_columns(header: list[str], src: config.Source) -> dict[str, str | None]:
    """Map chrom/pos/ref/alt/transcript/score to the file's column names."""
    lower = {h.lower(): h for h in header}
    out: dict[str, str | None] = {}
    for key, names in CANDIDATES.items():
        given = (src.columns or {}).get(key)
        found = given if given in header else next((lower[n] for n in names if n in lower), None)
        if given and given not in header:
            raise AviConfigError(f"avi.columns.{key} = {given!r} is not a column of the AVI file. Columns: {header}")
        if found is None and key != "transcript":
            raise AviConfigError(f"could not find the {key} column in the AVI file; set avi.columns.{key}. Columns: {header}")
        out[key] = found
    if src.score_field not in header:
        raise AviConfigError(f"avi.score_field = {src.score_field!r} is not a column of the AVI file. Columns: {header}")
    out["score"] = src.score_field
    return out


def _to_raw(df: pl.DataFrame | pl.LazyFrame, cols: dict) -> pl.DataFrame | pl.LazyFrame:
    return df.select(
        pl.col(cols["chrom"]).cast(pl.String).alias("chrom_raw"), pl.col(cols["pos"]).cast(pl.Int64).alias("pos"),
        pl.col(cols["ref"]).cast(pl.String).alias("ref_raw"), pl.col(cols["alt"]).cast(pl.String).alias("alt_raw"),
        pl.col(cols["score"]).cast(pl.Float64, strict=False).alias("avi_score"),
        (pl.col(cols["transcript"]).cast(pl.String) if cols.get("transcript") else pl.lit(None, dtype=pl.String)).alias("transcript"),
    )


def _both_styles(positions: pl.DataFrame) -> pl.DataFrame:
    """Scope positions under both chromosome spellings (chr17 and 17)."""
    return pl.concat([positions, positions.with_columns(pl.col("chrom").str.strip_prefix("chr"))]).rename({"chrom": "chrom_raw"})


def _is_local(uri: str) -> bool:
    return uri.startswith("file://") or "://" not in uri


def load(settings: Settings, src: config.Source, bed: pl.DataFrame, positions: pl.DataFrame, scope: str) -> tuple[pl.DataFrame, dict]:
    """Raw AVI records inside the scope regions, plus facts about the source for the report."""
    if not src.format or not src.score_field:
        raise AviConfigError("avi.format and avi.score_field must both be set in config/data_sources.yaml")
    d = settings.data_dir("avi")
    d.mkdir(parents=True, exist_ok=True)
    uri = src.uri.removeprefix("file://")
    want = _both_styles(positions)
    info = {"uri": src.uri, "format": src.format, "access": "", "columns": {}}

    if src.format == "parquet":
        scan = pl.scan_parquet(uri)
        cols = resolve_columns(scan.collect_schema().names(), src)
        raw = _to_raw(scan, cols).join(want.lazy(), on=["chrom_raw", "pos"], how="semi").collect()
        info.update(access="lazy parquet scan, filtered to scope positions", columns=cols)

    elif src.format == "vcf":
        tbi = d / f"{Path(uri).name}{src.index_suffix or '.tbi'}"
        provenance.fetch(settings, src.uri + (src.index_suffix or ".tbi"), tbi, stage="stage8", source_name="AVI index")
        regions = d / f"regions.{scope}.both.bed"
        pl.concat([bed, bed.with_columns(pl.col("chrom").str.strip_prefix("chr"))]).select("chrom", "start", "end").write_csv(
            regions, separator="\t", include_header=False)
        cmd = (f"set -o pipefail; bcftools view --no-version -R {regions} '{uri}##idx##{tbi}' -Ou | bcftools norm --no-version "
               f"-m -any -Ou | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/{src.score_field}\\n'")
        out = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True, check=True, cwd=d).stdout
        raw = pl.read_csv(out.encode(), separator="\t", has_header=False, null_values=["."], infer_schema_length=0,
                          new_columns=["chrom_raw", "pos", "ref_raw", "alt_raw", "avi_score"]) if out else pl.DataFrame(
            schema={k: pl.String for k in ("chrom_raw", "pos", "ref_raw", "alt_raw", "avi_score")})
        raw = raw.with_columns(pl.col("pos").cast(pl.Int64), pl.col("avi_score").cast(pl.Float64, strict=False),
                               pl.lit(None, dtype=pl.String).alias("transcript"))
        info.update(access="bcftools region query", columns={"score": f"INFO/{src.score_field}"})

    else:
        sep = "," if src.format == "csv" else "\t"
        subset = d / f"avi_subset.{scope}.{src.format}.gz"
        header = _header(uri, sep)
        cols = resolve_columns(header, src)
        index_uri = src.uri + (src.index_suffix or ".tbi")
        indexed = sep == "\t" and _exists(index_uri)
        if indexed:
            tbi = d / Path(index_uri.removeprefix("file://")).name
            provenance.fetch(settings, index_uri, tbi, stage="stage8", source_name="AVI index")
            regions = d / f"regions.{scope}.both.bed"
            pl.concat([bed, bed.with_columns(pl.col("chrom").str.strip_prefix("chr"))]).select("chrom", "start", "end").write_csv(
                regions, separator="\t", include_header=False)
            cmd = ["tabix", "-R", str(regions), f"{uri}##idx##{tbi}"]
            body = subprocess.run(cmd, capture_output=True, check=True, cwd=d).stdout
            with gzip.GzipFile(subset, "wb", mtime=0) as fh:
                fh.write((sep.join(header) + "\n").encode() + body)
            provenance.record(settings, subset, stage="stage8", source_name="AVI scores (region subset)", command=" ".join(cmd),
                              software_name="tabix", software_version="htslib", input_uri=src.uri, access_method="remote_tabix",
                              database_version=src.version or "", reference_build="GRCh38", input_info=provenance.head(src.uri))
            info["access"] = "tabix region query"
        else:
            ci, pi = header.index(cols["chrom"]), header.index(cols["pos"])
            keys = set(zip(want["chrom_raw"], want["pos"]))
            first = [True]

            def keep(line: str) -> bool:
                if first[0]:
                    first[0] = False
                    return True
                parts = line.rstrip("\n").split(sep)
                try:
                    return (parts[ci], int(parts[pi])) in keys
                except (IndexError, ValueError):
                    return False

            if not provenance.recorded_row(settings, subset, src.uri):
                provenance.stream_filter(settings, src.uri, subset, keep, stage="stage8", source_name="AVI scores",
                                         database_version=src.version or "", reference_build="GRCh38",
                                         filter_description="rows at positions inside the scope regions", gzip_output=True)
            info["access"] = "full stream, filtered to scope positions (no index found)"
        df = pl.read_csv(subset, separator=sep, infer_schema_length=0, null_values=["NA", "", ".", "nan", "NaN"])
        raw = _to_raw(df, cols).join(want, on=["chrom_raw", "pos"], how="semi")
        info["columns"] = cols

    head = provenance.head(src.uri) if _is_local(src.uri) or "*" not in src.uri else None
    info["last_modified"] = head.last_modified if head else ""
    return raw.select(list(RAW_SCHEMA)).cast(RAW_SCHEMA), info


def _exists(uri: str) -> bool:
    try:
        provenance.head(uri)
        return True
    except Exception:
        return False


def _header(uri: str, sep: str) -> list[str]:
    """Column names: the first line that is not a `##` comment."""
    for chunk in provenance._lines(provenance._gunzip_multi(provenance._open_stream(uri), lambda b: None)
                                   if uri.endswith((".gz", ".bgz")) else provenance._open_stream(uri)):
        line = chunk if isinstance(chunk, str) else chunk.decode()
        if not line.startswith("##"):
            return line.rstrip("\n").split(sep)
    raise AviConfigError("the AVI file is empty")
