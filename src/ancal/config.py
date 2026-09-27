"""Schema, validation and freezing for config/analysis.yaml and config/data_sources.yaml."""

from __future__ import annotations

import hashlib
from datetime import date
from pathlib import Path
from typing import Any, Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

GROUP_LABELS = ["AFR", "AMR", "CSA", "EAS", "EUR", "MID", "OCE"]


class Strict(BaseModel):
    model_config = ConfigDict(extra="forbid")


class Deferred(Strict):
    """A value the protocol leaves open or that depends on data not yet seen."""

    status: Literal["pending", "proposed", "confirmed"]
    settled_in: str
    value: Any = None
    note: str = ""

    @model_validator(mode="after")
    def _confirmed_needs_value(self):
        if self.status == "confirmed" and self.value is None:
            raise ValueError("a confirmed deferred value cannot be null")
        return self


# ---------------------------------------------------------------- analysis.yaml


class Protocol(Strict):
    source: str
    author: str
    canvas_created: date
    config_drafted: date


class PopulationGroups(Strict):
    labels: list[str]
    mapping_file: str
    mapping: Deferred

    @field_validator("labels")
    @classmethod
    def _seven_groups(cls, v):
        if v != GROUP_LABELS:
            raise ValueError(f"population group labels must be exactly {GROUP_LABELS}")
        return v


class Samples(Strict):
    primary_set: Literal["unrelated"]
    sensitivity_set: Literal["all_eligible"]
    exclude_admixed: Literal[False]
    hprc: Literal["excluded_from_primary"]


class DmsAssays(Strict):
    eligibility_criteria: list[str]
    min_variants_per_gene: int = Field(ge=1)
    eligible_assays: Deferred


class DmsScore(Strict):
    harmonized_direction: Literal["higher_is_more_damaging"]
    per_assay_direction_file: str
    pooling: Literal["within_assay_percentile"]
    replicate_combination: Deferred


class AfBin(Strict):
    name: str
    lower: float = Field(ge=0, le=1)
    upper: float = Field(ge=0, le=1)


class AfBins(Strict):
    frequency_source: Literal["global_combined"]
    bins: list[AfBin]
    interval: Literal["left_closed"]

    @field_validator("bins")
    @classmethod
    def _contiguous(cls, bins):
        if not bins or bins[0].lower != 0.0 or bins[-1].upper != 1.0:
            raise ValueError("allele-frequency bins must span [0, 1]")
        for a, b in zip(bins, bins[1:]):
            if a.upper != b.lower:
                raise ValueError(f"gap or overlap between bins {a.name} and {b.name}")
        for b in bins:
            if b.lower >= b.upper:
                raise ValueError(f"bin {b.name} has lower >= upper")
        return bins


class PrecisionCriteria(Strict):
    min_evaluable_variants: int = Field(ge=1)
    min_represented_genes: int = Field(ge=1)
    require_valid_spearman: bool
    max_bootstrap_ci_width: float = Field(gt=0, le=2)


class Metrics(Strict):
    primary: dict[str, str]
    secondary: list[str]


class Bootstrap(Strict):
    procedure: Literal["hierarchical_gene_then_variant"]
    iterations: int = Field(ge=1)
    ci_level: float = Field(gt=0, lt=1)
    seed: Deferred


class AfStandardization(Strict):
    method: Literal["stratified_resampling"]
    target_distribution: Literal["pooled_eligible_variants"]


class ClassicalComparator(Strict):
    name: str
    score_field: str
    version: Deferred
    evaluation_set: Literal["intersection_of_variants_scored_by_both"]


class ClinVar(Strict):
    release: Deferred
    consequence: Literal["missense"]
    positive_class: list[str]
    negative_class: list[str]
    exclude: list[str]
    keep_fields: list[str]
    vus_as_benign: Literal[False]


class VariantRules(Strict):
    primary_variant_type: Literal["biallelic_snv"]
    retain_other_types_in_master: bool
    variant_id_format: Literal["{chrom}:{pos}:{ref}:{alt}"]
    join_keys: list[Literal["chrom", "pos", "ref", "alt"]]
    source_filter: str
    observed_if: str
    exclude_population_variant_if: list[str]
    callability_proxy: str

    @field_validator("join_keys")
    @classmethod
    def _full_allele_key(cls, v):
        if sorted(v) != ["alt", "chrom", "pos", "ref"]:
            raise ValueError("joins must use chrom, pos, ref and alt (never rsID)")
        return v


class AnalysisConfig(Strict):
    protocol: Protocol
    reference_build: Literal["GRCh38"]
    chromosomes: Deferred
    population_groups: PopulationGroups
    samples: Samples
    dms_gene_universe: Deferred
    dms_assays: DmsAssays
    dms_score: DmsScore
    avi_score: Deferred
    missing_score_definition: Deferred
    allele_frequency_bins: AfBins
    precision_criteria: PrecisionCriteria
    metrics: Metrics
    bootstrap: Bootstrap
    af_standardization: AfStandardization
    classical_comparator: ClassicalComparator
    clinvar: ClinVar
    coding_definition: Deferred
    noncoding_analysis: str
    variant_rules: VariantRules


# ------------------------------------------------------------ data_sources.yaml


class Storage(Strict):
    root: str
    local_budget_gb: float = Field(gt=0)
    counted_output_dirs: list[str]


class ScopeDef(Strict):
    genes: list[str] | str
    contigs: list[str] | str


class Source(Strict):
    name: str
    uri: str | None
    access: str | None
    version: str | None = None
    index_suffix: str | None = None
    qc_summary_uri: str | None = None
    pc_scores_uri: str | None = None
    pca_outliers_uri: str | None = None
    related_ids_uri: str | None = None
    unrelated_mt_uri: str | None = None
    related_mt_uri: str | None = None
    use: str | None = None
    used_in: str | None = None
    format: str | None = None
    score_field: str | None = None
    direction: str | None = None


REQUIRED_SOURCES = (
    "population_vcf", "sample_metadata", "reference_fasta", "gene_annotation",
    "alphamissense", "clinvar", "mavedb", "cadd", "avi",
)


class DataSources(Strict):
    scope: Literal["smoke", "dms_genes", "exome"]
    storage: Storage
    scopes: dict[Literal["smoke", "dms_genes", "exome"], ScopeDef]
    sources: dict[str, Source]

    @model_validator(mode="after")
    def _all_sources_present(self):
        missing = [s for s in REQUIRED_SOURCES if s not in self.sources]
        if missing:
            raise ValueError(f"data_sources.yaml is missing sources: {missing}")
        for name, src in self.sources.items():
            if src.uri is None and name != "avi":
                raise ValueError(f"source {name!r} has no uri")
        if self.scope not in self.scopes:
            raise ValueError(f"scope {self.scope!r} has no entry under scopes:")
        return self


# ----------------------------------------------------------------- operations


def _load_yaml(path: Path) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh) or {}


def load_analysis(root: Path) -> AnalysisConfig:
    return AnalysisConfig.model_validate(_load_yaml(root / "config" / "analysis.yaml"))


def load_data_sources(root: Path) -> DataSources:
    return DataSources.model_validate(_load_yaml(root / "config" / "data_sources.yaml"))


def deferred_items(model: BaseModel, prefix: str = "") -> list[tuple[str, Deferred]]:
    """Every Deferred block in the config, as (dotted.path, block)."""
    out = []
    for name in type(model).model_fields:
        value = getattr(model, name)
        path = f"{prefix}{name}"
        if isinstance(value, Deferred):
            out.append((path, value))
        elif isinstance(value, BaseModel):
            out.extend(deferred_items(value, prefix=f"{path}."))
    return out


class NotFreezable(RuntimeError):
    pass


def freeze(root: Path, today: date | None = None) -> Path:
    """Write config/frozen/analysis_<date>.yaml (+ .sha256). Refuses while anything is unconfirmed."""
    cfg = load_analysis(root)
    open_items = [(p, d) for p, d in deferred_items(cfg) if d.status != "confirmed"]
    if open_items:
        lines = "\n".join(f"  - {p}: {d.status} (settled in {d.settled_in})" for p, d in open_items)
        raise NotFreezable(f"Cannot freeze: {len(open_items)} item(s) not confirmed:\n{lines}")

    today = today or date.today()
    src = root / "config" / "analysis.yaml"
    frozen_dir = root / "config" / "frozen"
    frozen_dir.mkdir(parents=True, exist_ok=True)
    out = frozen_dir / f"analysis_{today.isoformat()}.yaml"
    if out.exists():
        raise NotFreezable(f"{out.name} already exists; a frozen config is never overwritten")
    body = src.read_bytes()
    header = f"# FROZEN {today.isoformat()} from config/analysis.yaml. Do not edit.\n".encode()
    out.write_bytes(header + body)
    digest = hashlib.sha256(out.read_bytes()).hexdigest()
    out.with_suffix(".yaml.sha256").write_text(f"{digest}  {out.name}\n")
    return out
