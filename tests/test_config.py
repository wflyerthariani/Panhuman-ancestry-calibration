import shutil
from datetime import date
from pathlib import Path

import pytest
import yaml
from pydantic import ValidationError

from ancal import config

REPO = Path(__file__).resolve().parents[1]


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    (tmp_path / "config").mkdir()
    for f in ("analysis.yaml", "data_sources.yaml"):
        shutil.copy(REPO / "config" / f, tmp_path / "config" / f)
    return tmp_path


def _edit(repo: Path, fn):
    path = repo / "config" / "analysis.yaml"
    cfg = yaml.safe_load(path.read_text())
    fn(cfg)
    path.write_text(yaml.safe_dump(cfg))


def test_repo_configs_validate():
    config.load_analysis(REPO)
    config.load_data_sources(REPO)


def test_draft_has_pending_items():
    items = dict(config.deferred_items(config.load_analysis(REPO)))
    assert items["avi_score"].status == "pending"
    assert items["dms_gene_universe"].status == "pending"
    assert items["population_groups.mapping"].settled_in == "stage3"


def test_freeze_refuses_while_items_open(repo):
    with pytest.raises(config.NotFreezable, match="avi_score"):
        config.freeze(repo)
    assert not list((repo / "config").glob("frozen/*.yaml"))


def test_freeze_when_all_confirmed(repo):
    def confirm_all(cfg):
        def walk(node):
            if isinstance(node, dict):
                if {"status", "settled_in"} <= node.keys():
                    node["status"] = "confirmed"
                    if node.get("value") is None:
                        node["value"] = "x"
                for v in node.values():
                    walk(v)
        walk(cfg)

    _edit(repo, confirm_all)
    out = config.freeze(repo, today=date(2026, 10, 1))
    assert out.name == "analysis_2026-10-01.yaml"
    assert out.read_text().startswith("# FROZEN 2026-10-01")
    assert out.with_suffix(".yaml.sha256").exists()
    with pytest.raises(config.NotFreezable, match="never overwritten"):
        config.freeze(repo, today=date(2026, 10, 1))


def test_af_bins_must_be_contiguous(repo):
    _edit(repo, lambda c: c["allele_frequency_bins"]["bins"][1].update(lower=0.002))
    with pytest.raises(ValidationError, match="gap or overlap"):
        config.load_analysis(repo)


def test_rsid_join_rejected(repo):
    _edit(repo, lambda c: c["variant_rules"].update(join_keys=["chrom", "pos", "alt"]))
    with pytest.raises(ValidationError, match="never rsID"):
        config.load_analysis(repo)


def test_group_labels_fixed(repo):
    _edit(repo, lambda c: c["population_groups"].update(labels=["AFR", "EUR"]))
    with pytest.raises(ValidationError):
        config.load_analysis(repo)


def test_confirmed_value_cannot_be_null(repo):
    _edit(repo, lambda c: c["dms_gene_universe"].update(status="confirmed"))
    with pytest.raises(ValidationError, match="cannot be null"):
        config.load_analysis(repo)


def test_admixed_exclusion_forbidden(repo):
    _edit(repo, lambda c: c["samples"].update(exclude_admixed=True))
    with pytest.raises(ValidationError):
        config.load_analysis(repo)
