"""Unified per-residue pLDDT extractor for ColabFold AF2 and ESMFold outputs.

Both stacks emit per-residue pLDDT but in different files and conventions:
- ColabFold rank-1 score JSON: top-level "plddt" array (0–100 scale), one value
  per residue. The file is named '<seed>_scores_rank_<NNN>_*.json'.
- ESMFold (our `tools/esmfold_run.py`): a per-record metrics JSON with two
  per-residue arrays — 'per_residue_plddt' (heavy-atom-masked mean) and
  'per_residue_plddt_ca' (Cα atom). Both already on the 0–100 scale.

This helper smooths over both shapes so callers in Step 3+ can work in one
data structure (`PlddtTrack`).
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np

ESM_CONVENTIONS: tuple[str, ...] = ("heavy_atom_mean", "ca_only")


@dataclass(frozen=True)
class PlddtTrack:
    """A single per-residue pLDDT track on the 0–100 scale."""
    source: Literal["colabfold_af2", "esmfold"]
    convention: str  # 'colabfold_af2' | 'esmfold_heavy_atom_mean' | 'esmfold_ca_only'
    per_residue: np.ndarray  # shape [L], dtype float64

    @property
    def length(self) -> int:
        return int(self.per_residue.shape[0])

    @property
    def mean(self) -> float:
        return float(self.per_residue.mean()) if self.length else float("nan")

    @property
    def min(self) -> float:
        return float(self.per_residue.min()) if self.length else float("nan")

    @property
    def max(self) -> float:
        return float(self.per_residue.max()) if self.length else float("nan")

    def passes_gate(self, threshold: float) -> bool:
        return self.mean >= threshold

    def fraction_below(self, threshold: float) -> float:
        """Fraction of residues with pLDDT strictly below `threshold`.
        Useful for spotting flexible tails / disordered regions."""
        if not self.length:
            return float("nan")
        return float((self.per_residue < threshold).mean())


def from_colabfold_scores(path: str | Path) -> PlddtTrack:
    """Read a ColabFold '*_scores_rank_*.json' file.

    The schema has a top-level 'plddt' list of floats on the 0–100 scale.
    """
    path = Path(path)
    data = json.loads(path.read_text())
    if "plddt" not in data:
        raise KeyError(f"{path}: no top-level 'plddt' key (ColabFold scores schema)")
    arr = np.asarray(data["plddt"], dtype=np.float64)
    return PlddtTrack(source="colabfold_af2", convention="colabfold_af2", per_residue=arr)


def from_esmfold_metrics(
    path: str | Path,
    convention: Literal["heavy_atom_mean", "ca_only"] = "heavy_atom_mean",
) -> PlddtTrack:
    """Read an ESMFold per-record metrics JSON written by `tools/esmfold_run.py`.

    Both 'per_residue_plddt' (heavy-atom masked) and 'per_residue_plddt_ca' (Cα
    slot) are available; pick by `convention`.
    """
    if convention not in ESM_CONVENTIONS:
        raise ValueError(f"convention must be one of {ESM_CONVENTIONS}, got {convention!r}")
    path = Path(path)
    data = json.loads(path.read_text())
    key = "per_residue_plddt" if convention == "heavy_atom_mean" else "per_residue_plddt_ca"
    if key not in data:
        raise KeyError(f"{path}: no '{key}' key (ESMFold metrics schema)")
    arr = np.asarray(data[key], dtype=np.float64)
    return PlddtTrack(
        source="esmfold",
        convention=f"esmfold_{convention}",
        per_residue=arr,
    )


def find_colabfold_rank1(colabfold_outdir: str | Path) -> Path:
    """Return the rank-1 score JSON path for a ColabFold output directory.

    Assumes ColabFold's default file naming with '_scores_rank_001' in the name.
    """
    out = Path(colabfold_outdir)
    candidates = sorted(out.glob("*_scores_rank_001*.json"))
    if not candidates:
        raise FileNotFoundError(f"{out}: no '*_scores_rank_001*.json' file found")
    if len(candidates) > 1:
        raise RuntimeError(
            f"{out}: multiple rank-1 score JSONs found, expected one: "
            f"{[c.name for c in candidates]}"
        )
    return candidates[0]
