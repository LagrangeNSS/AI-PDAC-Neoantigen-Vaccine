# Download ledger (append-only, ISO 8601 UTC)

Budget: 20 GB total cumulative.
Format: `timestamp | URL | destination | size | running_total`

| Timestamp | URL | Destination | Size | Running total |
|-----------|-----|-------------|------|---------------|
| 2026-05-13T14:50Z | git+https://github.com/sokrypton/ColabFold (via Tsinghua PyPI mirror for deps; ColabFold itself from github.com direct) | ./env + ./cache/pip | ~2.2 GB downloaded (pip wheels for colabfold + alphafold-minus-jax deps: numpy, scipy, biopython, tensorflow_cpu, dm-tree, dm-haiku, chex, ml_collections, jaxlib, tqdm, requests, absl-py, etc.) | 2.2 GB |
| 2026-05-14T00:00Z | `pip install jax[cuda13]` (BFSU mirror — Tsinghua/pypi.org both hung on fastly CDN) | ./env + ./cache/pip | ~2.7 GB downloaded; nvidia-cublas 13.4.1.1 405 MB, nvidia-cudnn-cu13 9.22.0.52 513 MB, nvidia-cufft 12.2.0.46 218 MB, nvidia-cusolver 12.2.0.1 227 MB, nvidia-cusparse 12.7.10.1 151 MB, nvidia-nccl-cu13 2.30.4 213 MB, nvidia-cuda-runtime 13.2.75, nvidia-cuda-nvrtc 13.2.78, jax-cuda13-pjrt 119 MB, jax-cuda13-plugin 6 MB, others ~70 MB | ~4.9 GB cumulative |
| 2026-05-14T00:46Z | `pip install torch --index-url https://download.pytorch.org/whl/cu130` (PyTorch wheel CDN; direct, ~3.5 MB/s) | ./env + ./cache/pip | ~2.75 GB downloaded; torch-2.12.0+cu130 532.7 MB, triton 3.7.0 201.5 MB, nvidia-cublas 13.1.1.3 423.1 MB (replaces 13.4.1.1), nvidia-cudnn-cu13 9.20.0.48 366.2 MB (replaces 9.22.0.52), nvidia-cusparselt-cu13 0.8.1 170.1 MB, nvidia-nccl-cu13 2.29.7 206.0 MB (replaces 2.30.4), nvidia-nvshmem-cu13 3.4.5 60.4 MB (replaces 3.6.5), nvidia-cufft 12.0.0.61 214.1 MB (replaces 12.2.0.46), nvidia-cusolver 12.0.4.66 200.9 MB (replaces 12.2.0.1), nvidia-cusparse 12.6.3.3 145.9 MB (replaces 12.7.10.1), nvidia-cuda-nvrtc 13.0.88 90.2 MB, nvidia-curand 10.4.0.35 59.5 MB, nvidia-nvjitlink 13.0.88 40.7 MB, cuda-bindings 13.0.3 12.1 MB, sympy 1.14.0 6.3 MB, networkx 3.6.1 2.1 MB, nvidia-cuda-runtime 13.0.96 2.2 MB, nvidia-cuda-cupti 13.0.85 10.7 MB, nvidia-cufile 1.15.1.6 1.2 MB, smaller deps ~5 MB | ~7.65 GB cumulative |
| 2026-05-14T00:55Z | `pip install transformers accelerate pytest` (BFSU mirror) | ./env + ./cache/pip | ~50 MB downloaded; tokenizers 0.22.2 ~5 MB, regex 2026.5.9 ~0.8 MB, safetensors 0.7.0 ~0.5 MB, hf-xet 1.5.0 4.5 MB, transformers 5.8.0, accelerate 1.13.0, pytest 9.0.3, huggingface-hub 1.14.0, typer/click/httpx/anyio/h11/etc. | ~7.7 GB cumulative |
| 2026-05-14T01:10Z | ColabFold AF2 params via `colabfold_batch` first invocation (downloaded from storage.googleapis.com/alphafold/) | ./cache/xdg/colabfold/params | 3.47 GB; params_model_{1..5}.npz + params_model_{1..5}_ptm.npz (5 base + 5 ptm = 10 files, ~347 MB each); LICENSE file; download_finished.txt sentinel | ~11.17 GB cumulative |
| 2026-05-14T01:25Z | ColabFold MSA via mmseqs2 (api.colabfold.com → uniref/env DBs) | ./colabfold/smoke_ubiquitin | ~3 MB a3m (2 999 793 bytes) for the 76 aa ubiquitin smoke; remote search, no DB local download | ~11.17 GB cumulative |
| 2026-05-14T17:34Z | ESMFold v1 weights from huggingface.co (direct, CloudFront JFK50-P9 PoP, ~2.5–2.8 MB/s) | ./cache/huggingface/hub/models--facebook--esmfold_v1 | 8.44 GB pytorch_model.bin (8 442 062 570 B) + ~10 MB tokenizer/config/snapshot bookkeeping | ~19.62 GB cumulative |

Footprint on disk after Step 1 smoke (du -sh):
- `./env` 7.7 GB
- `./cache/pip` 6.2 GB
- `./cache/huggingface` 9.0 GB (includes 8.44 GB ESMFold blob)
- `./cache/xdg/colabfold/params` 3.5 GB
- `./cache/tmp` 1.8 GB (JAX compilation + transient)
- `./cache/jax` 2.4 MB
- `./cache/matplotlib` 32 KB

Cumulative downloads vs budget: 19.62 GB / 20 GB. Headroom: 0.38 GB for the remaining pipeline.

Notes:
- jax's nvidia-* wheels (downloaded 2.7 GB) were later replaced by PyTorch's pinned versions (another 2.0 GB of overlapping nvidia-* wheels). The newer-version copies sit in `./cache/pip` but are not installed in `./env`. This is the cost of running jax[cuda13] and torch+cu130 side-by-side.
- Headroom is tight. Future steps must avoid model-weight downloads. MSA queries (a3m text streams) and per-residue feature JSONs are well within remaining 380 MB.
