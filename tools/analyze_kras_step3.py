"""Step 3 analysis: WT vs G12D structural diff + per-residue pLDDT around mutation.

Reads ColabFold AF2 rank-1 + ESMFold outputs for both WT and G12D KRAS (1-169),
computes pairwise Cα RMSD (full G-domain, structured core, P-loop), and dumps
per-residue pLDDT tracks. All paths are absolute under WORKDIR.
"""
from __future__ import annotations

import json
from pathlib import Path

from tools.plddt import find_colabfold_rank1, from_colabfold_scores, from_esmfold_metrics
from tools.rmsd import ca_rmsd

WORKDIR = Path("/mnt/d/AssignmentOPUS")

CF_WT_DIR = WORKDIR / "colabfold" / "kras_wt_1-169"
CF_G12D_DIR = WORKDIR / "colabfold" / "kras_g12d_1-169"
ESM_WT_DIR = WORKDIR / "esmfold" / "kras_wt_1-169"
ESM_G12D_DIR = WORKDIR / "esmfold" / "kras_g12d_1-169"
OUT_DIR = WORKDIR / "analysis"
OUT_DIR.mkdir(exist_ok=True)


def _find_rank1_pdb(d: Path) -> Path:
    pdbs = list(d.glob("*unrelaxed_rank_001*.pdb"))
    if not pdbs:
        raise FileNotFoundError(f"No rank-1 PDB under {d}")
    return pdbs[0]


def _esmfold_pdb(d: Path) -> Path:
    pdbs = [p for p in d.glob("*.pdb") if "_metrics" not in p.name and "_summary" not in p.name]
    if len(pdbs) != 1:
        raise FileNotFoundError(f"Expected exactly 1 ESMFold PDB under {d}, got {len(pdbs)}")
    return pdbs[0]


def main() -> None:
    cf_wt_pdb = _find_rank1_pdb(CF_WT_DIR)
    cf_g12d_pdb = _find_rank1_pdb(CF_G12D_DIR)
    esm_wt_pdb = _esmfold_pdb(ESM_WT_DIR)
    esm_g12d_pdb = _esmfold_pdb(ESM_G12D_DIR)

    print(f"AF2 WT   : {cf_wt_pdb.name}")
    print(f"AF2 G12D : {cf_g12d_pdb.name}")
    print(f"ESM WT   : {esm_wt_pdb.name}")
    print(f"ESM G12D : {esm_g12d_pdb.name}")
    print()

    # KRAS structural regions (1-based, inclusive).
    P_LOOP = range(10, 18)            # residues 10-17
    SWITCH_I = range(30, 39)          # residues 30-38
    SWITCH_II = range(59, 73)         # residues 59-72
    STRUCTURED_CORE = range(1, 167)   # 1-166 (drops 167-169 linker → tight comparison)
    FULL = range(1, 170)              # 1-169

    pair_specs = [
        ("AF2_WT_vs_AF2_G12D", cf_wt_pdb, cf_g12d_pdb),
        ("ESM_WT_vs_ESM_G12D", esm_wt_pdb, esm_g12d_pdb),
        ("AF2_WT_vs_ESM_WT", cf_wt_pdb, esm_wt_pdb),
        ("AF2_G12D_vs_ESM_G12D", cf_g12d_pdb, esm_g12d_pdb),
    ]

    region_specs = [
        ("full_1-169", FULL),
        ("core_1-166", STRUCTURED_CORE),
        ("P-loop_10-17", P_LOOP),
        ("Switch_I_30-38", SWITCH_I),
        ("Switch_II_59-72", SWITCH_II),
    ]

    results = {}
    print("=== Pairwise Cα RMSD (Å) ===")
    print(f"{'pair':<24} | " + " | ".join(f"{r[0]:>14}" for r in region_specs))
    print("-" * (26 + len(region_specs) * 17))
    for pair_name, a, b in pair_specs:
        results[pair_name] = {}
        row = [f"{pair_name:<24}"]
        for region_name, residues in region_specs:
            r = ca_rmsd(str(a), str(b), chain_id="A", residue_subset=residues)
            results[pair_name][region_name] = {
                "rmsd_A": r.rmsd_angstrom,
                "n_ca": r.n_atoms,
            }
            row.append(f"{r.rmsd_angstrom:>10.3f} ({r.n_atoms:>2})")
        print(" | ".join(row))

    # Per-residue pLDDT tracks
    af2_wt_track = from_colabfold_scores(find_colabfold_rank1(CF_WT_DIR))
    af2_g12d_track = from_colabfold_scores(find_colabfold_rank1(CF_G12D_DIR))
    esm_wt_metrics_path = next(ESM_WT_DIR.glob("*_metrics.json"))
    esm_g12d_metrics_path = next(ESM_G12D_DIR.glob("*_metrics.json"))
    # heavy-atom-masked convention (matches the smoke fixture convention).
    esm_wt_track = from_esmfold_metrics(esm_wt_metrics_path, convention="heavy_atom_mean")
    esm_g12d_track = from_esmfold_metrics(esm_g12d_metrics_path, convention="heavy_atom_mean")

    means = {
        "AF2_WT_mean_plddt": af2_wt_track.mean,
        "AF2_G12D_mean_plddt": af2_g12d_track.mean,
        "ESM_WT_mean_plddt_heavy_atom": esm_wt_track.mean,
        "ESM_G12D_mean_plddt_heavy_atom": esm_g12d_track.mean,
    }
    print()
    print("=== Mean pLDDT (rank-1 / single ESMFold pass) ===")
    for k, v in means.items():
        print(f"{k:<35} = {v:.2f}   gate≥85: {'PASS' if v >= 85 else 'FAIL'}")

    # Per-residue pLDDT around the mutation site (window pos 4-20)
    print()
    print("=== Per-residue pLDDT, residues 4-20 (mutation site context) ===")
    print(f"{'pos':>3} | {'AF2 WT':>8} | {'AF2 G12D':>9} | {'ESM WT':>8} | {'ESM G12D':>9} | mark")
    print("----+----------+-----------+----------+-----------+-----")
    per_res = []
    for pos in range(4, 21):
        a_wt = af2_wt_track.per_residue[pos - 1]
        a_mu = af2_g12d_track.per_residue[pos - 1]
        e_wt = esm_wt_track.per_residue[pos - 1]
        e_mu = esm_g12d_track.per_residue[pos - 1]
        mark = " <-- G12D mutation site" if pos == 12 else ""
        print(f"{pos:>3} | {a_wt:>8.2f} | {a_mu:>9.2f} | {e_wt:>8.2f} | {e_mu:>9.2f} |{mark}")
        per_res.append({
            "pos_1based": pos,
            "af2_wt_plddt": float(a_wt),
            "af2_g12d_plddt": float(a_mu),
            "esm_wt_plddt": float(e_wt),
            "esm_g12d_plddt": float(e_mu),
            "is_mutation_site": pos == 12,
        })

    # Persist all results.
    out = {
        "inputs": {
            "af2_wt_pdb": str(cf_wt_pdb),
            "af2_g12d_pdb": str(cf_g12d_pdb),
            "esm_wt_pdb": str(esm_wt_pdb),
            "esm_g12d_pdb": str(esm_g12d_pdb),
            "af2_wt_scores": str(find_colabfold_rank1(CF_WT_DIR)),
            "af2_g12d_scores": str(find_colabfold_rank1(CF_G12D_DIR)),
            "esm_wt_metrics": str(esm_wt_metrics_path),
            "esm_g12d_metrics": str(esm_g12d_metrics_path),
        },
        "means": means,
        "rmsd_A": results,
        "per_residue_window_4-20": per_res,
        "full_per_residue_plddt": {
            "af2_wt": [float(x) for x in af2_wt_track.per_residue],
            "af2_g12d": [float(x) for x in af2_g12d_track.per_residue],
            "esm_wt_heavy_atom_mean": [float(x) for x in esm_wt_track.per_residue],
            "esm_g12d_heavy_atom_mean": [float(x) for x in esm_g12d_track.per_residue],
        },
    }
    out_path = OUT_DIR / "kras_wt_vs_g12d_structural.json"
    out_path.write_text(json.dumps(out, indent=2))
    print()
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
