"""Paths and small helpers shared by Stages 6-11."""

from __future__ import annotations

import subprocess
from pathlib import Path

import polars as pl

from .. import config, provenance
from ..settings import Settings


class Paths:
    """Every scope-dependent file written by Stages 6-11."""

    def __init__(self, s: Settings):
        sc = s.scope
        self.bed = s.root / "resources" / "regions" / f"{sc}.bed"
        self.tx = s.data_dir("reference") / "mane_select_transcripts.parquet"
        self.cds = s.data_dir("reference") / "mane_select_cds.parquet"
        self.alphamissense = s.data_dir("reference") / f"alphamissense_hg38.{sc}.tsv.gz"
        self.samples = s.data_dir("population_metadata") / "sample_metadata.parquet"
        self.consequences = s.data_dir("intermediate") / f"cds_snv_consequences.{sc}.parquet"
        self.dms = s.data_dir("dms") / f"dms_harmonized.{sc}.parquet"
        self.dms_exclusions = s.data_dir("dms") / f"dms_harmonization_exclusions.{sc}.parquet"
        self.sites = s.data_dir("variants") / f"population_sites.{sc}.parquet"
        self.long = s.data_dir("variants") / f"population_variants_long.{sc}.parquet"
        self.avi_scores = s.data_dir("avi") / f"avi_scores.{sc}.parquet"
        self.avi_join = s.data_dir("analysis_datasets") / f"avi_join.{sc}.parquet"
        self.clinvar = s.data_dir("clinvar") / f"clinvar_eligibility.{sc}.parquet"
        self.exclusion_log = s.data_dir("analysis_datasets") / f"variant_exclusion_log.{sc}.parquet"
        self.tables = s.root / "tables"
        self.qc = s.root / "reports" / "qc"
        self.exclusions = s.root / "reports" / "exclusions"
        self.d1 = s.root / "deliverables" / "d1"


def read_bed(path: Path) -> pl.DataFrame:
    return pl.read_csv(path, separator="\t", has_header=False, new_columns=["chrom", "start", "end", "name"],
                       schema_overrides={"start": pl.Int64, "end": pl.Int64})


def bed_positions(bed: pl.DataFrame) -> pl.DataFrame:
    """Every 1-based position inside a 0-based half-open BED, as (chrom, pos)."""
    return (bed.select("chrom", pl.int_ranges(pl.col("start") + 1, pl.col("end") + 1).alias("pos"))
            .explode("pos").unique().sort("chrom", "pos"))


def scope_genes(bed: pl.DataFrame) -> list[str]:
    return sorted({g for name in bed["name"] for g in name.split(",")})


def af_bin_expr(cfg: config.AnalysisConfig, col: str = "af_global") -> pl.Expr:
    """Global allele-frequency bin name (Step 16); null where the frequency is null."""
    bins = cfg.allele_frequency_bins.bins
    e = pl.when(pl.col(col).is_null()).then(None)
    for b in bins[:-1]:
        e = e.when(pl.col(col) < b.upper).then(pl.lit(b.name))
    return e.otherwise(pl.lit(bins[-1].name)).alias("af_bin")


def spearman(df: pl.DataFrame, a: str, b: str) -> float | None:
    d = df.select(a, b).drop_nulls()
    if d.height < 3 or d[a].n_unique() < 2 or d[b].n_unique() < 2:
        return None
    return d.select(pl.corr(pl.col(a).rank("average"), pl.col(b).rank("average"))).item()


def tool_version(tool: str) -> str:
    return subprocess.run([tool, "--version"], capture_output=True, text=True).stdout.splitlines()[0]


def fmt(v) -> str:
    if v is None:
        return "—"
    if isinstance(v, float):
        return f"{v:.4g}"
    if isinstance(v, int) and not isinstance(v, bool):
        return f"{v:,}"
    return str(v)


def md_table(df: pl.DataFrame, headers: list[str] | None = None) -> str:
    head = headers or df.columns
    lines = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    lines += ["| " + " | ".join(fmt(v) for v in r) + " |" for r in df.iter_rows()]
    return "\n".join(lines)


def write_tsv(df: pl.DataFrame, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(path, separator="\t", float_precision=6)
    return path


def record(s: Settings, stage: str, module: str, inputs: str, *paths: Path) -> None:
    for p in paths:
        provenance.record(s, p, stage=stage, source_name="derived", command=f"ancal {stage}",
                          software_name=f"ancal.stages.{module}", software_version=provenance.code_version(s.root),
                          input_uri=inputs, reference_build="GRCh38")


def provisional_banner(scope: str, extra: list[str] | None = None) -> str:
    notes = list(extra or [])
    if scope == "smoke":
        notes.insert(0, "produced at `smoke` scope (a few genes on one chromosome); rerun at `dms_genes` scope")
    if not notes:
        return ""
    return "> **PROVISIONAL.** " + "; ".join(notes) + ".\n\n"
