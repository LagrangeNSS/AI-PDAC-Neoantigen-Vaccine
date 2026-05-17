# Step 0 — Environment audit

Probe timestamp (UTC): 2026-05-13T14:27:57Z
Working directory ($WORKDIR): /mnt/d/AssignmentOPUS

## Operating system
- `uname -a` -> `Linux LAPTOP-HUQMD29R 6.6.87.2-microsoft-standard-WSL2 #1 SMP PREEMPT_DYNAMIC Thu Jun  5 18:30:46 UTC 2025 x86_64 GNU/Linux`
- `lsb_release -a` -> `Ubuntu 24.04 LTS (noble)`
- Matches the fixed execution environment in the brief (WSL2 Ubuntu 24.04).

## Python
- `python3 --version` -> `Python 3.12.3` (within required 3.10–3.12 range).
- `which python3` -> `/usr/bin/python3`
- `python3-venv` 3.12.3-0ubuntu2.1 installed.

## Disk
- `df -h /mnt/d/AssignmentOPUS` -> Filesystem `D:\`, Size 700G, Used 488G, **Avail 213G**, Use 70%.
- Comfortably above the 30 GB minimum.

## RAM
- `free -h` -> Total 7.4 GiB, Available 6.6 GiB, Swap 2.0 GiB.
- Note: brief states 16 GB host RAM, but WSL2 exposes only ~7.6 GiB to Linux by default (no `.wslconfig` override detected). All Python work must fit in this. The 12 GB GPU VRAM provides headroom for model weights, so this is workable for KRAS (189 aa) but worth flagging.

## CPU
- `nproc` -> 24 logical cores.

## GPU
- `nvidia-smi`:
  - GPU: `NVIDIA GeForce RTX 5070 Ti Laptop GPU`
  - Driver Version: **595.79** (>= 580, so use CUDA 13 wheels)
  - CUDA Version reported by driver: **13.2**
  - Memory: **12227 MiB (12 GB)** total; 1694 MiB in use at probe time
  - Compute Mode: Default
- `nvcc` not present on PATH (expected; JAX/PyTorch wheels bundle their own CUDA runtime libs).
- Architecture: Blackwell (sm_120) — Laptop variant of RTX 5070 Ti.

## Network reachability (HEAD probes, max-time 10s)
| Host | HTTP code | Note |
|------|-----------|------|
| github.com | 200 | OK |
| pypi.org | 200 | OK |
| huggingface.co | 200 | OK |
| www.uniprot.org | 200 | OK |
| alphafold.ebi.ac.uk | 200 | OK |
| api.colabfold.com | 404 | Expected; root path returns 404 (nginx route-only API). Server is alive — `colabfold_batch` uses POST /ticket/msa. Will validate with a real smoke prediction in Step 1. |

## Basic toolchain
- `git` 2.43.0 — OK
- `curl` 8.5.0 — OK
- `gcc` (build-essential) 13.3.0 — OK

## Decision logic outcomes
- Driver 595.79 >= 580 → **install CUDA 13 wheels**:
  - JAX: `pip install --upgrade "jax[cuda13]"`
  - PyTorch: `pip install --upgrade torch --index-url https://download.pytorch.org/whl/cu130`
  - If cu130 stable wheels lack sm_120 (Blackwell) support, fall back to `--pre` nightlies for both.

## Known platform caveats (must surface in README & FINAL_RESEARCH_SUMMARY)
1. **JAX RNG non-determinism on Blackwell sm_120 GPUs.** Per NVIDIA's JAX release notes, the JAX random number generator is non-deterministic on gamer Blackwell (sm_120) GPUs. Even with `--random-seed 42` we may observe small run-to-run variation in ColabFold outputs. We still lock the seed to minimize variance.
2. **WSL2 RAM ceiling ~7.6 GiB.** Models that need to host-load large state could OOM. Mitigation: run models on GPU (12 GB VRAM), avoid holding multiple ESMFold instances simultaneously, and free GPU/CPU memory between jobs.
3. **ColabFold MSA service is a runtime dependency.** Online MSA mode `mmseqs2_uniref_env` is mandatory by brief; if the server is unavailable, brief requires 3-retry exponential backoff (60s, 300s, 900s) then stop — no silent single-sequence fallback.
4. **Per-task timeout: 2 h wall time.** Any single prediction exceeding 2 h must be aborted with partial outputs preserved under ./logs/timeouts/<jobid>/.
5. **Download budget: 20 GB total.** Logged append-only in ./logs/downloads.md.

## Next step
Await user confirmation before proceeding to Step 1 (venv creation, package installation, smoke tests).
