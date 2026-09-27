"""Load storage/scope settings from config/data_sources.yaml."""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

import yaml

VALID_SCOPES = ("smoke", "dms_genes", "exome")


def repo_root() -> Path:
    """Repository root: $ANCAL_ROOT if set, else the nearest parent holding config/data_sources.yaml."""
    env = os.environ.get("ANCAL_ROOT")
    if env:
        return Path(env).resolve()
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "config" / "data_sources.yaml").exists():
            return parent
    return Path.cwd().resolve()


@dataclass(frozen=True)
class Settings:
    root: Path
    scope: str
    storage_root: Path
    local_budget_gb: float
    counted_output_dirs: tuple[Path, ...] = field(default_factory=tuple)

    @property
    def budget_bytes(self) -> int:
        return int(self.local_budget_gb * 1024**3)

    def data_dir(self, name: str) -> Path:
        return self.storage_root / name


def load_settings(root: Path | None = None, scope_override: str | None = None) -> Settings:
    root = (root or repo_root()).resolve()
    with open(root / "config" / "data_sources.yaml") as fh:
        cfg = yaml.safe_load(fh) or {}

    scope = scope_override or os.environ.get("ANCAL_SCOPE") or cfg.get("scope", "smoke")
    if scope not in VALID_SCOPES:
        raise ValueError(f"scope must be one of {VALID_SCOPES}, got {scope!r}")

    storage = cfg.get("storage", {})
    storage_root = Path(storage.get("root", "./data_store"))
    if not storage_root.is_absolute():
        storage_root = root / storage_root

    counted = tuple(root / d for d in storage.get("counted_output_dirs", []))
    return Settings(
        root=root,
        scope=scope,
        storage_root=storage_root.resolve(),
        local_budget_gb=float(storage.get("local_budget_gb", 1.0)),
        counted_output_dirs=counted,
    )
