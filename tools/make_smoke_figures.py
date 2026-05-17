"""Render the two smoke figures (per-residue pLDDT) for ColabFold AF2 and ESMFold
on the ubiquitin smoke target. Each figure gets a .txt sidecar with provenance.

Outputs:
    figures/00_smoke_colabfold.png  + .txt
    figures/00_smoke_esmfold.png    + .txt
"""
from __future__ import annotations

import datetime as dt
import json
import os
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path("/mnt/d/AssignmentOPUS")
FIG_DIR = ROOT / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

CF_SCORES = next((ROOT / "colabfold" / "smoke_ubiquitin").glob("*scores_rank_001*.json"))
CF_LOG = ROOT / "logs" / "colabfold_runs" / "smoke_ubiquitin.log"
ESM_METRICS = ROOT / "esmfold" / "smoke_ubiquitin" / "smoke_ubiquitin_P0CG48_1-76_metrics.json"
SEQ_FASTA = ROOT / "data" / "sequences" / "ubiquitin_smoke.fasta"


def read_fasta_single(path: Path) -> tuple[str, str]:
    rid, body = None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            rid = line[1:].split()[0]
        elif line.strip():
            body.append(line.strip())
    return rid, "".join(body)


def gpu_info() -> str:
    try:
        out = subprocess.check_output(
            ["nvidia-smi", "--query-gpu=name,driver_version,memory.total",
             "--format=csv,noheader"], text=True, timeout=10).strip()
        return out
    except Exception as e:
        return f"nvidia-smi unavailable: {e}"


def write_sidecar(png: Path, lines: list[str]) -> None:
    side = png.with_suffix(".txt")
    side.write_text("\n".join(lines) + "\n")


def plot_plddt(per_res: list[float], title: str, png: Path, color: str,
               gate: float = 85.0) -> tuple[float, float, float]:
    arr = np.asarray(per_res, dtype=float)
    L = len(arr)
    fig, ax = plt.subplots(figsize=(7.5, 3.6))
    ax.plot(np.arange(1, L + 1), arr, color=color, linewidth=1.8)
    ax.fill_between(np.arange(1, L + 1), arr, 0, color=color, alpha=0.12)
    ax.axhline(gate, color="grey", linestyle="--", linewidth=1.0,
               label=f"hard gate pLDDT={gate:.0f}")
    ax.axhline(arr.mean(), color="black", linestyle=":", linewidth=1.0,
               label=f"mean = {arr.mean():.1f}")
    ax.set_xlabel("residue index")
    ax.set_ylabel("pLDDT (0–100)")
    ax.set_title(title)
    ax.set_xlim(1, L)
    ax.set_ylim(0, 100)
    ax.grid(True, alpha=0.25)
    ax.legend(loc="lower left", fontsize=8)
    fig.tight_layout()
    fig.savefig(png, dpi=200)
    plt.close(fig)
    return float(arr.mean()), float(arr.min()), float(arr.max())


def main() -> int:
    rid, seq = read_fasta_single(SEQ_FASTA)
    gpu = gpu_info()
    ts = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    # ---- ColabFold figure ----
    cf = json.loads(CF_SCORES.read_text())
    cf_plddt = cf["plddt"]
    cf_ptm = cf.get("ptm")
    cf_png = FIG_DIR / "00_smoke_colabfold.png"
    cf_mean, cf_min, cf_max = plot_plddt(
        cf_plddt,
        title=f"ColabFold AF2 rank 1 — ubiquitin {rid}",
        png=cf_png,
        color="#1f77b4",
    )
    write_sidecar(cf_png, [
        f"# {cf_png.name} — provenance",
        f"generated_utc: {ts}",
        f"tool: ColabFold (sokrypton/ColabFold HEAD, alphafold-minus-jax weights via colabfold_batch)",
        "flags (brief-locked): --num-models 5 --num-recycle 3 --rank plddt --model-type alphafold2_ptm --random-seed 42 --msa-mode mmseqs2_uniref_env --disable-unified-memory",
        f"sequence: {rid} (len={len(seq)}); {seq}",
        f"rank-1 model: {CF_SCORES.name}",
        f"rank-1 mean pLDDT: {cf_mean:.2f}  min: {cf_min:.2f}  max: {cf_max:.2f}",
        f"rank-1 pTM: {cf_ptm if cf_ptm is not None else 'NA'}",
        f"hard gate (≥85 mean pLDDT): {'PASS' if cf_mean >= 85 else 'FAIL'}",
        f"GPU: {gpu}",
        f"source score json: {CF_SCORES}",
        f"source run log:    {CF_LOG}",
    ])
    print(f"[CF ] mean_pLDDT={cf_mean:.2f} -> {cf_png}")

    # ---- ESMFold figure ----
    esm = json.loads(ESM_METRICS.read_text())
    # Use the masked-heavy-atom mean (matches the smoke metric gated above).
    esm_plddt = esm["per_residue_plddt"]
    esm_png = FIG_DIR / "00_smoke_esmfold.png"
    esm_mean, esm_min, esm_max = plot_plddt(
        esm_plddt,
        title=f"ESMFold v1 — ubiquitin {esm['record_id']}",
        png=esm_png,
        color="#d62728",
    )
    write_sidecar(esm_png, [
        f"# {esm_png.name} — provenance",
        f"generated_utc: {ts}",
        "tool: ESMFold v1 (facebook/esmfold_v1 via transformers 5.8.0, single-sequence, no MSA)",
        "precision: ESM-2 stem fp16, folding trunk + heads + adapters fp32 on GPU "
        "(loaded fp16 to fit WSL2 7.4 GiB host RAM, then upcast non-esm parts to fp32 "
        "to avoid NaN in pTM softmax)",
        f"sequence: {esm['record_id']} (len={esm['length']}); {seq}",
        "per_residue_plddt convention: mean over masked heavy atoms (atom37_atom_exists), scaled x100",
        f"mean pLDDT (heavy-atom mask): {esm_mean:.2f}  min: {esm_min:.2f}  max: {esm_max:.2f}",
        f"mean pLDDT (Cα-only):         {esm['mean_plddt_ca']:.2f}",
        f"inference_seconds: {esm['inference_seconds']}",
        f"hard gate (≥85 mean pLDDT): {'PASS' if esm_mean >= 85 else 'FAIL'}",
        f"GPU: {gpu}",
        f"source metrics json: {ESM_METRICS}",
        "seed: torch.manual_seed(42)",
    ])
    print(f"[ESM] mean_pLDDT={esm_mean:.2f} -> {esm_png}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
