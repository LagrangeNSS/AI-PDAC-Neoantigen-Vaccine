"""Unit tests for tools.plddt — uses the live smoke outputs as fixtures so we
catch any regressions in the AF2 / ESMFold file schemas at the same time.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from tools.plddt import (
    PlddtTrack,
    find_colabfold_rank1,
    from_colabfold_scores,
    from_esmfold_metrics,
)

ROOT = Path("/mnt/d/AssignmentOPUS")
CF_OUT = ROOT / "colabfold" / "smoke_ubiquitin"
ESM_METRICS = ROOT / "esmfold" / "smoke_ubiquitin" / "smoke_ubiquitin_P0CG48_1-76_metrics.json"


def test_find_rank1_returns_single_file() -> None:
    path = find_colabfold_rank1(CF_OUT)
    assert "_scores_rank_001" in path.name
    assert path.exists()


def test_colabfold_rank1_smoke_passes_gate() -> None:
    track = from_colabfold_scores(find_colabfold_rank1(CF_OUT))
    assert track.source == "colabfold_af2"
    assert track.length == 76
    assert track.mean >= 85.0  # brief-locked hard gate
    assert track.min > 0 and track.max <= 100.0
    assert track.passes_gate(85.0) is True


def test_esmfold_smoke_passes_gate_heavy_atom() -> None:
    track = from_esmfold_metrics(ESM_METRICS, convention="heavy_atom_mean")
    assert track.source == "esmfold"
    assert track.convention == "esmfold_heavy_atom_mean"
    assert track.length == 76
    assert track.passes_gate(85.0) is True


def test_esmfold_smoke_ca_convention_higher_than_heavy_atom() -> None:
    ca = from_esmfold_metrics(ESM_METRICS, convention="ca_only").mean
    ha = from_esmfold_metrics(ESM_METRICS, convention="heavy_atom_mean").mean
    # Cα-only excludes side-chain atoms that drag the mean down on the flexible
    # C-terminal tail; it should be the higher of the two.
    assert ca > ha


def test_invalid_convention_rejected() -> None:
    with pytest.raises(ValueError):
        from_esmfold_metrics(ESM_METRICS, convention="banana")  # type: ignore[arg-type]


def test_fraction_below_threshold_matches_manual() -> None:
    arr = np.array([10.0, 50.0, 90.0, 99.0], dtype=np.float64)
    track = PlddtTrack(source="esmfold", convention="esmfold_ca_only", per_residue=arr)
    assert track.fraction_below(70.0) == pytest.approx(0.5)
    assert track.fraction_below(100.0) == pytest.approx(1.0)
    assert track.fraction_below(0.0) == pytest.approx(0.0)


def test_colabfold_scores_missing_plddt_raises(tmp_path: Path) -> None:
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps({"ptm": 0.8}))
    with pytest.raises(KeyError):
        from_colabfold_scores(bad)


def test_find_rank1_errors_on_empty_dir(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        find_colabfold_rank1(tmp_path)
