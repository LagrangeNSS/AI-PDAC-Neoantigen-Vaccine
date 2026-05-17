"""Cα RMSD between two PDB structures via Biopython.

Used in Step 3 to cross-validate ColabFold AF2 vs ESMFold predictions and to
compare WT vs KRAS-G12D structures. We deliberately keep this minimal:
- Cα atoms only (no heavy-atom RMSD).
- Sequence-aligned residue pairing by integer residue number on a single chain
  (default: chain 'A') — both predictors emit chain A for monomer inputs.
- Kabsch superposition via Biopython's SVDSuperimposer.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.SVDSuperimposer import SVDSuperimposer

_PARSER = PDBParser(QUIET=True)


@dataclass(frozen=True)
class RmsdResult:
    rmsd_angstrom: float
    n_atoms: int
    chain_id: str
    residue_range: tuple[int, int]  # 1-based inclusive (first_resnum, last_resnum) intersected


def _load_ca(path: str | Path, chain_id: str) -> dict[int, np.ndarray]:
    """Return {residue_number: Cα xyz} for the named chain of the first model."""
    structure = _PARSER.get_structure(Path(path).stem, str(path))
    model = next(iter(structure))  # first model only
    try:
        chain = model[chain_id]
    except KeyError:
        chains = [c.id for c in model]
        raise KeyError(
            f"{path}: chain {chain_id!r} not found; available chains: {chains}"
        )
    ca_by_resnum: dict[int, np.ndarray] = {}
    for residue in chain:
        # Skip HETATM and waters: standard residues have hetflag ' '.
        hetflag, resseq, _icode = residue.id
        if hetflag.strip():
            continue
        if "CA" not in residue:
            continue
        ca_by_resnum[int(resseq)] = np.asarray(residue["CA"].coord, dtype=np.float64)
    return ca_by_resnum


def ca_rmsd(
    pdb_a: str | Path,
    pdb_b: str | Path,
    *,
    chain_id: str = "A",
    residue_subset: range | tuple[int, int] | None = None,
) -> RmsdResult:
    """Cα RMSD between two PDB files after Kabsch superposition.

    Args:
        pdb_a, pdb_b: PDB paths.
        chain_id: single-letter chain id to compare in both structures.
        residue_subset: optional restriction. If a (lo, hi) tuple is given,
                       only residues with `lo <= resnum <= hi` are used.
                       If a range is given, only resnums in `range` are used.
                       If None, the intersection of both chains' resnums is used.

    Returns:
        RmsdResult with `rmsd_angstrom`, `n_atoms`, and the resnum range used.

    Raises:
        ValueError: if no common Cα atoms remain after intersection.
    """
    a = _load_ca(pdb_a, chain_id)
    b = _load_ca(pdb_b, chain_id)

    common = sorted(set(a) & set(b))
    if residue_subset is not None:
        if isinstance(residue_subset, range):
            allowed = set(residue_subset)
        else:
            lo, hi = residue_subset
            allowed = set(range(int(lo), int(hi) + 1))
        common = [r for r in common if r in allowed]
    if not common:
        raise ValueError(
            f"No common Cα residues between {pdb_a} and {pdb_b} "
            f"(chain {chain_id!r}, subset={residue_subset})"
        )

    coords_a = np.stack([a[r] for r in common], axis=0)
    coords_b = np.stack([b[r] for r in common], axis=0)

    sup = SVDSuperimposer()
    sup.set(coords_a, coords_b)
    sup.run()
    return RmsdResult(
        rmsd_angstrom=float(sup.get_rms()),
        n_atoms=len(common),
        chain_id=chain_id,
        residue_range=(min(common), max(common)),
    )
