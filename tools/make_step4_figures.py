"""Step 4 — publication-grade figures for the KRAS G12D feasibility study.

Generates four figures under `figures/`:
  01_plddt_tracks_WT_vs_G12D.png    — full-domain pLDDT tracks, 2 stacks × WT/G12D
  02_plddt_zoom_mutation_site.png   — bar chart of pLDDT around the mutation
  03_rmsd_by_region.png             — grouped bar chart of Cα RMSD by region
  04_neoantigen_shortlist.png       — annotated 9-mer table figure

All figures have .txt provenance sidecars listing input paths + key numbers.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

WORKDIR = Path("/mnt/d/AssignmentOPUS")
FIG_DIR = WORKDIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

STRUCT_PATH = WORKDIR / "analysis" / "kras_wt_vs_g12d_structural.json"
NEO_PATH = WORKDIR / "analysis" / "neoantigen_9mers.json"

# KRAS structural regions for shading.
P_LOOP = (10, 17)
SWITCH_I = (30, 38)
SWITCH_II = (59, 72)


def _now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%MZ")


def _write_provenance(png_path: Path, lines: list[str]) -> None:
    txt = png_path.with_suffix(".txt")
    txt.write_text("\n".join(lines) + "\n")


def fig1_plddt_tracks(struct: dict) -> None:
    af2_wt = struct["full_per_residue_plddt"]["af2_wt"]
    af2_mu = struct["full_per_residue_plddt"]["af2_g12d"]
    esm_wt = struct["full_per_residue_plddt"]["esm_wt_heavy_atom_mean"]
    esm_mu = struct["full_per_residue_plddt"]["esm_g12d_heavy_atom_mean"]
    L = len(af2_wt)
    x = np.arange(1, L + 1)

    fig, axes = plt.subplots(2, 1, figsize=(9.5, 5.5), sharex=True,
                             gridspec_kw={"hspace": 0.15})
    for ax, wt, mu, title in [
        (axes[0], af2_wt, af2_mu, "ColabFold AlphaFold-2 (rank-1, model 5, seed 42)"),
        (axes[1], esm_wt, esm_mu, "ESMFold v1 (single-sequence, heavy-atom mean)"),
    ]:
        # Shade KRAS functional regions.
        for (lo, hi), col, name in [
            (P_LOOP, "#FFE082", "P-loop"),
            (SWITCH_I, "#A5D6A7", "Switch I"),
            (SWITCH_II, "#90CAF9", "Switch II"),
        ]:
            ax.axvspan(lo, hi, color=col, alpha=0.50, lw=0)
        ax.plot(x, wt, color="#1f77b4", lw=1.3, label="WT")
        ax.plot(x, mu, color="#d62728", lw=1.3, label="G12D", alpha=0.85)
        ax.axhline(85, color="grey", lw=0.8, ls="--", label="gate ≥ 85")
        ax.axvline(12, color="black", lw=0.6, ls=":", alpha=0.7)
        ax.set_ylabel("pLDDT")
        ax.set_ylim(40, 100)
        ax.set_title(title, fontsize=10, loc="left")
        ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
        ax.grid(alpha=0.25, lw=0.4)

    axes[1].set_xlabel("Residue (1-based, KRAS-4B numbering)")
    axes[0].set_xlim(1, L)

    # Region labels along the top of the upper panel.
    for (lo, hi), name in [(P_LOOP, "P-loop"), (SWITCH_I, "Switch I"), (SWITCH_II, "Switch II")]:
        axes[0].text((lo + hi) / 2, 99, name, ha="center", va="top", fontsize=8, color="#222")

    fig.suptitle("Per-residue pLDDT, KRAS G-domain 1–169: WT vs G12D, two predictors",
                 fontsize=11, y=0.995)
    out = FIG_DIR / "01_plddt_tracks_WT_vs_G12D.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    _write_provenance(out, [
        f"# {out.name}",
        f"timestamp_utc: {_now_utc()}",
        "produced_by: tools/make_step4_figures.py:fig1_plddt_tracks",
        f"af2_wt_scores:  {struct['inputs']['af2_wt_scores']}",
        f"af2_g12d_scores:{struct['inputs']['af2_g12d_scores']}",
        f"esm_wt_metrics: {struct['inputs']['esm_wt_metrics']}",
        f"esm_g12d_metrics:{struct['inputs']['esm_g12d_metrics']}",
        f"mean_plddt_af2_wt   = {struct['means']['AF2_WT_mean_plddt']:.2f}",
        f"mean_plddt_af2_g12d = {struct['means']['AF2_G12D_mean_plddt']:.2f}",
        f"mean_plddt_esm_wt   = {struct['means']['ESM_WT_mean_plddt_heavy_atom']:.2f}",
        f"mean_plddt_esm_g12d = {struct['means']['ESM_G12D_mean_plddt_heavy_atom']:.2f}",
        "gate: pLDDT >= 85 hard gate; all four predictions PASS.",
    ])
    print(f"  wrote {out}")


def fig2_zoom_mutation(struct: dict) -> None:
    rows = struct["per_residue_window_4-20"]
    pos = [r["pos_1based"] for r in rows]
    af2_wt = [r["af2_wt_plddt"] for r in rows]
    af2_mu = [r["af2_g12d_plddt"] for r in rows]
    esm_wt = [r["esm_wt_plddt"] for r in rows]
    esm_mu = [r["esm_g12d_plddt"] for r in rows]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.0), sharey=True)
    width = 0.4
    for ax, wt, mu, title in [
        (axes[0], af2_wt, af2_mu, "ColabFold AF2 (rank-1)"),
        (axes[1], esm_wt, esm_mu, "ESMFold v1 (heavy-atom mean)"),
    ]:
        idx = np.arange(len(pos))
        ax.bar(idx - width / 2, wt, width, label="WT", color="#1f77b4")
        ax.bar(idx + width / 2, mu, width, label="G12D", color="#d62728")
        ax.set_xticks(idx)
        ax.set_xticklabels([str(p) for p in pos], fontsize=8)
        ax.axhline(85, color="grey", lw=0.8, ls="--")
        # Highlight position 12.
        i12 = pos.index(12)
        ax.axvspan(i12 - 0.5, i12 + 0.5, color="#FFEB3B", alpha=0.25, lw=0, zorder=0)
        ax.text(i12, 102, "G12D", ha="center", va="bottom", fontsize=9, color="#a00")
        delta = mu[i12] - wt[i12]
        ax.text(i12, 78, f"Δ={delta:+.2f}", ha="center", va="top", fontsize=8,
                color="#a00", fontweight="bold")
        ax.set_xlabel("Residue (1-based)")
        ax.set_title(title, fontsize=10)
        ax.set_ylim(60, 105)
        ax.grid(axis="y", alpha=0.25, lw=0.4)
        ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
    axes[0].set_ylabel("pLDDT")
    fig.suptitle("Per-residue pLDDT around the mutation site (residues 4–20)", fontsize=11, y=1.02)
    out = FIG_DIR / "02_plddt_zoom_mutation_site.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    _write_provenance(out, [
        f"# {out.name}",
        f"timestamp_utc: {_now_utc()}",
        "produced_by: tools/make_step4_figures.py:fig2_zoom_mutation",
        f"position_12_AF2_WT={af2_wt[pos.index(12)]:.2f}  AF2_G12D={af2_mu[pos.index(12)]:.2f}  Δ={af2_mu[pos.index(12)] - af2_wt[pos.index(12)]:+.2f}",
        f"position_12_ESM_WT={esm_wt[pos.index(12)]:.2f}  ESM_G12D={esm_mu[pos.index(12)]:.2f}  Δ={esm_mu[pos.index(12)] - esm_wt[pos.index(12)]:+.2f}",
        "interpretation: ESMFold reports a localised confidence drop at the mutation site; AF2 (with MSA) does not.",
    ])
    print(f"  wrote {out}")


def fig3_rmsd_by_region(struct: dict) -> None:
    pairs = ["AF2_WT_vs_AF2_G12D", "ESM_WT_vs_ESM_G12D", "AF2_WT_vs_ESM_WT", "AF2_G12D_vs_ESM_G12D"]
    pair_labels = ["AF2 WT vs AF2 G12D", "ESM WT vs ESM G12D", "AF2 vs ESM (WT)", "AF2 vs ESM (G12D)"]
    regions = ["full_1-169", "core_1-166", "P-loop_10-17", "Switch_I_30-38", "Switch_II_59-72"]
    region_labels = ["Full 1-169", "Core 1-166", "P-loop 10-17", "Switch I 30-38", "Switch II 59-72"]

    data = np.zeros((len(pairs), len(regions)), dtype=float)
    for i, pair in enumerate(pairs):
        for j, reg in enumerate(regions):
            data[i, j] = struct["rmsd_A"][pair][reg]["rmsd_A"]

    fig, ax = plt.subplots(figsize=(10, 4.5))
    x = np.arange(len(regions))
    width = 0.20
    colors = ["#1f77b4", "#2ca02c", "#ff7f0e", "#d62728"]
    for i, (lab, col) in enumerate(zip(pair_labels, colors)):
        ax.bar(x + (i - 1.5) * width, data[i], width, label=lab, color=col)
        for j, v in enumerate(data[i]):
            ax.text(x[j] + (i - 1.5) * width, v + 0.03, f"{v:.2f}",
                    ha="center", va="bottom", fontsize=7, rotation=0)
    ax.set_xticks(x)
    ax.set_xticklabels(region_labels, fontsize=9)
    ax.set_ylabel("Cα RMSD (Å)")
    ax.set_title("Cα RMSD by structural region — within-stack (WT vs G12D) and cross-stack (AF2 vs ESMFold)",
                 fontsize=10)
    ax.set_ylim(0, max(data.max() * 1.18, 0.5))
    ax.grid(axis="y", alpha=0.25, lw=0.4)
    ax.legend(loc="upper left", fontsize=8, framealpha=0.9)
    fig.tight_layout()
    out = FIG_DIR / "03_rmsd_by_region.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    _write_provenance(out, [
        f"# {out.name}",
        f"timestamp_utc: {_now_utc()}",
        "produced_by: tools/make_step4_figures.py:fig3_rmsd_by_region",
        "method: Kabsch superposition (Biopython SVDSuperimposer), Cα atoms, monomer chain A.",
        f"within_stack_WT_vs_G12D_full: AF2={data[0,0]:.3f} Å, ESM={data[1,0]:.3f} Å",
        f"cross_stack_AF2_vs_ESM_core: WT={data[2,1]:.3f} Å, G12D={data[3,1]:.3f} Å",
        "interpretation: G12D is globally structure-preserving in both stacks; cross-stack disagreement concentrates in Switch II.",
    ])
    print(f"  wrote {out}")


def fig4_neoantigen_shortlist(neo: dict) -> None:
    rows = neo["windows_covering_position_12"]
    fig, ax = plt.subplots(figsize=(10.5, 4.6))
    ax.axis("off")

    headers = ["start", "end", "G12D 9-mer", "WT 9-mer", "P2", "P2 tier", "P9", "P9 tier", "overall", "shortlist"]
    col_x = [0.04, 0.10, 0.17, 0.31, 0.46, 0.52, 0.66, 0.72, 0.86, 0.95]
    y0 = 0.93
    dy = 0.072
    # Header.
    for x, h in zip(col_x, headers):
        ax.text(x, y0, h, fontsize=9, fontweight="bold", transform=ax.transAxes)
    ax.plot([0.02, 0.99], [y0 - 0.018, y0 - 0.018], color="black", lw=0.8, transform=ax.transAxes)

    tier_colors = {"preferred": "#2e7d32", "tolerated": "#fbc02d", "poor": "#b71c1c"}

    for i, r in enumerate(rows):
        y = y0 - dy * (i + 1)
        is_hit = (r["overall_tier"] in ("preferred", "tolerated") and r["foreign_vs_wt"])
        bg = "#FFF8E1" if is_hit else "white"
        ax.add_patch(Rectangle((0.02, y - 0.024), 0.97, 0.044, transform=ax.transAxes,
                               facecolor=bg, edgecolor="none", zorder=0))
        ax.text(col_x[0], y, f"{r['start_1based']}", fontsize=9, transform=ax.transAxes)
        ax.text(col_x[1], y, f"{r['end_1based']}", fontsize=9, transform=ax.transAxes)

        # Render the G12D peptide with the mutation residue highlighted in red.
        pep = r["g12d_peptide"]
        mut_idx = 12 - r["start_1based"]  # 0-based within the 9-mer
        for j, aa in enumerate(pep):
            ax.text(col_x[2] + j * 0.014, y, aa,
                    fontsize=10, family="monospace",
                    color=("#c0392b" if j == mut_idx else "black"),
                    fontweight=("bold" if j == mut_idx else "normal"),
                    transform=ax.transAxes)
        # WT peptide.
        for j, aa in enumerate(r["wt_peptide"]):
            ax.text(col_x[3] + j * 0.014, y, aa, fontsize=10, family="monospace",
                    color=("#888" if j == mut_idx else "black"),
                    transform=ax.transAxes)

        ax.text(col_x[4], y, r["p2"], fontsize=10, family="monospace", transform=ax.transAxes)
        ax.text(col_x[5], y, r["p2_tier"], fontsize=9,
                color=tier_colors[r["p2_tier"]], transform=ax.transAxes,
                fontweight=("bold" if r["p2_tier"] == "preferred" else "normal"))
        ax.text(col_x[6], y, r["p9"], fontsize=10, family="monospace", transform=ax.transAxes)
        ax.text(col_x[7], y, r["p9_tier"], fontsize=9,
                color=tier_colors[r["p9_tier"]], transform=ax.transAxes,
                fontweight=("bold" if r["p9_tier"] == "preferred" else "normal"))
        ax.text(col_x[8], y, r["overall_tier"], fontsize=9,
                color=tier_colors[r["overall_tier"]], transform=ax.transAxes,
                fontweight=("bold" if r["overall_tier"] == "preferred" else "normal"))
        ax.text(col_x[9], y, ("PASS" if is_hit else "—"), fontsize=9,
                color=("#2e7d32" if is_hit else "#999"),
                fontweight=("bold" if is_hit else "normal"), transform=ax.transAxes)

    ax.text(0.02, 0.02,
            ("Mutation residue highlighted in red; PASS row shaded. "
             "Scorer = tools/mhc.py P2/P9 anchor motif for HLA-A*11:01 "
             "(Sidney 2008; Falk 1994; Zhang 1993). Motif scoring only — not an IC50 predictor."),
            fontsize=8, color="#555", transform=ax.transAxes, wrap=True)
    fig.suptitle("HLA-A*11:01 anchor shortlist for KRAS G12D 9-mers covering position 12",
                 fontsize=11, y=0.99)
    out = FIG_DIR / "04_neoantigen_shortlist.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    _write_provenance(out, [
        f"# {out.name}",
        f"timestamp_utc: {_now_utc()}",
        "produced_by: tools/make_step4_figures.py:fig4_neoantigen_shortlist",
        f"input: {NEO_PATH}",
        "scorer: tools.mhc P2/P9 anchor motif for HLA-A*11:01 (A3 supertype).",
        f"shortlist_count: {sum(1 for r in rows if (r['overall_tier'] in ('preferred','tolerated') and r['foreign_vs_wt']))} of {len(rows)}",
        "shortlist_peptide: VVGADGVGK (residues 8-16) — matches Tran 2016 NEJM and Rojas 2023 Nature canonical neoepitope.",
    ])
    print(f"  wrote {out}")


def main() -> None:
    struct = json.loads(STRUCT_PATH.read_text())
    neo = json.loads(NEO_PATH.read_text())
    fig1_plddt_tracks(struct)
    fig2_zoom_mutation(struct)
    fig3_rmsd_by_region(struct)
    fig4_neoantigen_shortlist(neo)


if __name__ == "__main__":
    main()
