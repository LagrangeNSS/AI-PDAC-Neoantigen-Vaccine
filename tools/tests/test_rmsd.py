"""Unit tests for tools.rmsd — uses the smoke PDBs as live fixtures."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from Bio.PDB import PDBIO, Atom, Chain, Model, Residue, Structure

from tools.rmsd import ca_rmsd

ROOT = Path("/mnt/d/AssignmentOPUS")
CF_RANK1_PDB = next((ROOT / "colabfold" / "smoke_ubiquitin").glob("*_unrelaxed_rank_001*.pdb"))
ESM_PDB = ROOT / "esmfold" / "smoke_ubiquitin" / "smoke_ubiquitin_P0CG48_1-76.pdb"


def _write_simple_chain(path: Path, coords: np.ndarray, chain_id: str = "A") -> None:
    """Write a synthetic PDB with Cα atoms at the given coords (1-indexed residues)."""
    structure = Structure.Structure("test")
    model = Model.Model(0)
    chain = Chain.Chain(chain_id)
    for i, xyz in enumerate(coords, start=1):
        residue = Residue.Residue((" ", i, " "), "GLY", " ")
        atom = Atom.Atom(
            name="CA",
            coord=np.asarray(xyz, dtype=np.float64),
            bfactor=0.0,
            occupancy=1.0,
            altloc=" ",
            fullname=" CA ",
            serial_number=i,
            element="C",
        )
        residue.add(atom)
        chain.add(residue)
    model.add(chain)
    structure.add(model)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(path))


def test_identical_coords_give_zero_rmsd(tmp_path: Path) -> None:
    coords = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]], dtype=np.float64)
    a = tmp_path / "a.pdb"
    b = tmp_path / "b.pdb"
    _write_simple_chain(a, coords)
    _write_simple_chain(b, coords)
    r = ca_rmsd(a, b)
    assert r.rmsd_angstrom == pytest.approx(0.0, abs=1e-6)
    assert r.n_atoms == 4
    assert r.chain_id == "A"
    assert r.residue_range == (1, 4)


def test_pure_translation_gives_zero_after_kabsch(tmp_path: Path) -> None:
    coords_a = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64)
    coords_b = coords_a + np.array([10.0, -5.0, 7.0])
    a = tmp_path / "a.pdb"
    b = tmp_path / "b.pdb"
    _write_simple_chain(a, coords_a)
    _write_simple_chain(b, coords_b)
    # Kabsch removes the translation; RMSD should be ~0.
    r = ca_rmsd(a, b)
    assert r.rmsd_angstrom == pytest.approx(0.0, abs=1e-6)


def test_residue_subset_tuple(tmp_path: Path) -> None:
    coords_a = np.array([[0, 0, 0], [10, 0, 0], [0, 0, 0], [0, 0, 0]], dtype=np.float64)
    coords_b = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]], dtype=np.float64)
    a = tmp_path / "a.pdb"
    b = tmp_path / "b.pdb"
    _write_simple_chain(a, coords_a)
    _write_simple_chain(b, coords_b)
    # Excluding the displaced residue 2 → identical coords → zero RMSD.
    r = ca_rmsd(a, b, residue_subset=(3, 4))
    assert r.rmsd_angstrom == pytest.approx(0.0, abs=1e-6)
    assert r.n_atoms == 2
    assert r.residue_range == (3, 4)


def test_smoke_cf_vs_esm_ubiquitin_topology_matches() -> None:
    """ColabFold AF2 vs ESMFold on the same ubiquitin sequence should give
    a globally similar fold. Cα RMSD on the structured core (residues 1–71,
    dropping the flexible Gly-Gly tail) should be < 3 Å."""
    r = ca_rmsd(CF_RANK1_PDB, ESM_PDB, residue_subset=(1, 71))
    assert r.n_atoms == 71
    # Empirically this is well under 2 Å for ubiquitin core. Use 3 Å as a
    # conservative gate so the test stays stable across runs.
    assert r.rmsd_angstrom < 3.0, (
        f"Cα RMSD on ubiquitin core (1-71) was {r.rmsd_angstrom:.2f} Å, "
        f"expected < 3.0 Å. This may indicate one stack diverged from the other."
    )


def test_missing_chain_raises(tmp_path: Path) -> None:
    coords = np.array([[0, 0, 0]], dtype=np.float64)
    a = tmp_path / "a.pdb"
    _write_simple_chain(a, coords, chain_id="A")
    with pytest.raises(KeyError):
        ca_rmsd(a, a, chain_id="Z")
