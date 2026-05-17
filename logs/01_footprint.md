# Step 1 — Disk footprint snapshot

Taken 2026-05-13T19:21Z (UTC), immediately after smoke tests passed.

## Top-level (`du -sh /mnt/d/AssignmentOPUS/<dir>`)

| Path | Size | Notes |
|---|---|---|
| `cache/` | 21 GB | tooling caches (see breakdown below) |
| `env/` | 7.7 GB | venv site-packages (torch, jax, alphafold-colabfold, transformers, …) |
| `colabfold/` | 7.7 MB | smoke outputs (5 PDBs + score JSONs + MSA a3m + PNGs) |
| `logs/` | 1.3 MB | install + smoke transcripts, downloads ledger |
| `figures/` | 164 KB | 2 smoke pLDDT plots + .txt provenance sidecars |
| `esmfold/` | 52 KB | smoke PDB + per-record + summary JSONs |
| `tools/` | 16 KB | `esmfold_run.py`, `make_smoke_figures.py` |
| `data/` | 0 (1 small FASTA) | ubiquitin_smoke.fasta |
| **TOTAL** | **28 GB** | |

## Cache breakdown (`du -sh /mnt/d/AssignmentOPUS/cache/<sub>`)

| Path | Size | Provenance |
|---|---|---|
| `cache/huggingface/` | 9.0 GB | 8.44 GB `pytorch_model.bin` for ESMFold v1 + tokenizer/config snapshots + xet hashes |
| `cache/pip/` | 6.2 GB | wheel cache; includes overlapping nvidia-* wheels from jax[cuda13] (now stale) and the active torch+cu130 set |
| `cache/xdg/colabfold/params/` | 3.5 GB | 10 AF2 `.npz` weight files (5 base + 5 PTM) + LICENSE + `download_finished.txt` |
| `cache/tmp/` | 1.8 GB | TMPDIR redirect — JAX XLA compilation artifacts + transient ColabFold scratch |
| `cache/jax/` | 2.4 MB | `JAX_COMPILATION_CACHE_DIR` (persistent kernel cache) |
| `cache/matplotlib/` | 32 KB | font cache |
| `cache/torch/` | 0 | TORCH_HOME (no torch hub downloads yet) |
| `cache/colabfold/` | 0 | superseded by ColabFold's actual XDG path (`cache/xdg/colabfold/`); kept as redirected env var pointer |

## Budget accounting

| Item | Bytes downloaded |
|---|---|
| pip wheels (5 batches) | 7.70 GB |
| ColabFold AF2 params (GCS, first run of `colabfold_batch`) | 3.47 GB |
| ColabFold MSA (mmseqs2 remote a3m) | 0.003 GB |
| ESMFold weights (HF direct, CloudFront) | 8.44 GB |
| **Cumulative downloads** | **19.62 GB / 20 GB** |

Remaining budget: **0.38 GB**. Sufficient for future MSA queries (~3–10 MB per sequence) but does NOT permit any additional model-weight pulls.

## Redirected env vars (verify with `./env/activate_project.sh`)

```
PIP_CACHE_DIR            = ./cache/pip
HF_HOME                  = ./cache/huggingface
HUGGINGFACE_HUB_CACHE    = ./cache/huggingface/hub
TRANSFORMERS_CACHE       = ./cache/huggingface
TORCH_HOME               = ./cache/torch
JAX_COMPILATION_CACHE_DIR= ./cache/jax
COLABFOLD_DATA_DIR       = ./cache/colabfold        (legacy; ColabFold actually used XDG_CACHE_HOME)
XDG_CACHE_HOME           = ./cache/xdg              (this is where AF2 params landed)
TMPDIR                   = ./cache/tmp
MPLCONFIGDIR             = ./cache/matplotlib
```

Nothing escaped to `~/.cache`, `~/.local`, or `/tmp` — verified by `du -sh ~/.cache /tmp 2>/dev/null` returning <100 KB.

## Cleanup opportunities (deferred)

These could reclaim ~6 GB but are intentionally left in place for Step 1 reproducibility/audit:

- `cache/pip/` (6.2 GB) — wheel cache; safe to `pip cache purge` once Step 3 verifies all builds reproducible
- `cache/tmp/` (1.8 GB) — JAX/XLA transient; will be re-created on next run
- The duplicated nvidia-* wheels in `cache/pip/` (~2 GB) from the jax[cuda13] → torch+cu130 swap

Reclamation plan deferred to Step 7 (cleanup).
