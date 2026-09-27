from pathlib import Path

import pytest
import yaml

from ancal import budget, status
from ancal.settings import load_settings


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    (tmp_path / "config").mkdir()
    cfg = {
        "scope": "smoke",
        "storage": {"root": "./data_store", "local_budget_gb": 1e-6, "counted_output_dirs": ["tables"]},
    }
    (tmp_path / "config" / "data_sources.yaml").write_text(yaml.safe_dump(cfg))
    (tmp_path / "data_store").mkdir()
    (tmp_path / "tables").mkdir()
    return tmp_path


def test_usage_counts_storage_and_outputs(repo):
    (repo / "data_store" / "a.bin").write_bytes(b"x" * 300)
    (repo / "tables" / "t.tsv").write_bytes(b"y" * 200)
    u = budget.usage(load_settings(repo))
    assert u.used_bytes == 500
    assert u.by_dir == {"data_store": 300, "tables": 200}


def test_check_raises_over_budget(repo):
    s = load_settings(repo)  # budget ~1073 bytes
    (repo / "data_store" / "big.bin").write_bytes(b"x" * 2000)
    with pytest.raises(budget.BudgetExceeded):
        budget.check(s)


def test_ensure_room_refuses_projected_overflow(repo):
    s = load_settings(repo)
    budget.ensure_room(s, 100, "small file")
    with pytest.raises(budget.BudgetExceeded, match="would exceed"):
        budget.ensure_room(s, 10_000, "large file")


def test_invalid_scope_rejected(repo):
    with pytest.raises(ValueError):
        load_settings(repo, scope_override="genome")


def test_status_pass_fail_blocked(repo):
    s = load_settings(repo)
    ok = status.Check("a", True)
    bad = status.Check("b", False)
    assert status.write(s, "stage9", [ok]).status == "pass"
    assert status.write(s, "stage9", [ok, bad]).status == "fail"
    assert status.write(s, "stage9", [ok], blocked_reason="no access").status == "blocked"
    assert status.read(s, "stage9")["blocked_reason"] == "no access"
