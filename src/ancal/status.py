"""Per-stage status files (reports/status/stageN.json).

status is one of:
  pass    - all gate checks passed
  fail    - a gate check failed; downstream stages must not run
  blocked - the stage cannot run yet (missing access, budget); recorded in item 15
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from . import budget
from .settings import Settings

VALID = ("pass", "fail", "blocked")


@dataclass
class Check:
    name: str
    ok: bool
    detail: str = ""


@dataclass
class StageStatus:
    stage: str
    status: str
    scope: str
    disk_used_mb: float
    checks: list[Check] = field(default_factory=list)
    notes: str = ""
    blocked_reason: str = ""
    written_at: str = ""


def status_path(settings: Settings, stage: str) -> Path:
    return settings.root / "reports" / "status" / f"{stage}.json"


def write(
    settings: Settings,
    stage: str,
    checks: list[Check],
    notes: str = "",
    blocked_reason: str = "",
) -> StageStatus:
    if blocked_reason:
        status = "blocked"
    elif all(c.ok for c in checks):
        status = "pass"
    else:
        status = "fail"
    st = StageStatus(
        stage=stage,
        status=status,
        scope=settings.scope,
        disk_used_mb=round(budget.usage(settings).used_mb, 2),
        checks=checks,
        notes=notes,
        blocked_reason=blocked_reason,
        written_at=datetime.now(timezone.utc).isoformat(timespec="seconds"),
    )
    path = status_path(settings, stage)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(asdict(st), indent=2) + "\n")
    return st


def read(settings: Settings, stage: str) -> dict | None:
    path = status_path(settings, stage)
    if not path.exists():
        return None
    return json.loads(path.read_text())


def read_all(settings: Settings) -> list[dict]:
    d = settings.root / "reports" / "status"
    return [json.loads(p.read_text()) for p in sorted(d.glob("stage*.json"))]
