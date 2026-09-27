"""Local disk budget guard.

Every stage calls `check()` after it runs and `ensure_room()` before any download,
so that local data never exceeds `storage.local_budget_gb`.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .settings import Settings


class BudgetExceeded(RuntimeError):
    pass


def dir_size(path: Path) -> int:
    if not path.exists():
        return 0
    if path.is_file():
        return path.stat().st_size
    total = 0
    for p in path.rglob("*"):
        if p.is_file() and not p.is_symlink():
            total += p.stat().st_size
    return total


@dataclass(frozen=True)
class Usage:
    used_bytes: int
    budget_bytes: int
    by_dir: dict[str, int]

    @property
    def used_mb(self) -> float:
        return self.used_bytes / 1024**2

    @property
    def remaining_bytes(self) -> int:
        return self.budget_bytes - self.used_bytes

    @property
    def ok(self) -> bool:
        return self.used_bytes <= self.budget_bytes


def usage(settings: Settings) -> Usage:
    dirs = [settings.storage_root, *settings.counted_output_dirs]
    by_dir = {}
    for d in dirs:
        try:
            label = str(d.relative_to(settings.root))
        except ValueError:
            label = str(d)
        by_dir[label] = dir_size(d)
    return Usage(sum(by_dir.values()), settings.budget_bytes, by_dir)


def check(settings: Settings) -> Usage:
    """Raise BudgetExceeded if local usage is over the cap."""
    u = usage(settings)
    if not u.ok:
        raise BudgetExceeded(
            f"Local data uses {u.used_mb:.1f} MB, over the {settings.local_budget_gb} GB budget. "
            "See SCALING_UP.md to raise storage.local_budget_gb or move storage.root."
        )
    return u


def ensure_room(settings: Settings, incoming_bytes: int, what: str) -> None:
    """Refuse a download that would push usage over the cap."""
    u = usage(settings)
    if u.used_bytes + incoming_bytes > u.budget_bytes:
        need_mb = (u.used_bytes + incoming_bytes - u.budget_bytes) / 1024**2
        raise BudgetExceeded(
            f"Fetching {what} (~{incoming_bytes / 1024**2:.1f} MB) would exceed the "
            f"{settings.local_budget_gb} GB budget by {need_mb:.1f} MB (scope={settings.scope}). "
            "See SCALING_UP.md."
        )
