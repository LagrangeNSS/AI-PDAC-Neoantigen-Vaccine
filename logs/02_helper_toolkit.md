# Step 2 — Helper toolkit

Window: 2026-05-13T19:21Z → 2026-05-13T19:30Z (UTC).

## Modules under `tools/`

| Module | Public surface | Purpose |
|---|---|---|
| `tools/seqio.py` | `FastaRecord`, `read_fasta`, `write_fasta`, `validate_protein`, `apply_point_mutation`, `slice_window`, `iter_kmers` | Light FASTA I/O + HGVS-style point mutation application (KRAS G12D pipeline) without pulling in Biopython for the read path |
| `tools/retry.py` | `retry`, `call_with_timeout`, `RetryStats`, `CallTimeoutError`, `DEFAULT_BACKOFF_SECONDS=(60,300,900)`, `DEFAULT_PER_CALL_TIMEOUT_SECONDS=7200` | Brief-locked 60/300/900 s exponential backoff + SIGALRM 2 h per-attempt timeout |
| `tools/plddt.py` | `PlddtTrack`, `from_colabfold_scores`, `from_esmfold_metrics`, `find_colabfold_rank1` | Unified per-residue pLDDT loader for ColabFold AF2 rank-1 JSON and ESMFold metrics JSON; both on the 0–100 scale |
| `tools/rmsd.py` | `RmsdResult`, `ca_rmsd` | Cα-only Kabsch RMSD between two PDB files (Biopython `SVDSuperimposer`); supports residue-subset restriction |
| `tools/mhc.py` | `AnchorTier`, `AnchorScore`, `score_9mer_hla_a1101`, `score_9mers_covering_position`, `filter_passing` | HLA-A*11:01 9-mer P2/P9 anchor scorer (preferred/tolerated/poor); conservative, motif-only — NOT a binding affinity predictor |
| `tools/msa.py` | `LOCKED_FLAGS`, `HARDWARE_FLAGS`, `run_colabfold`, `ColabFoldRun`, `ColabFoldGateError` | Subprocess wrapper around `colabfold_batch` that pins the brief-locked flag set, applies `tools.retry`, and enforces the rank-1 mean pLDDT ≥ 85 hard gate |

## Brief-locked policy enforced in code

`tools/msa.LOCKED_FLAGS`:
```
--num-models 5 --num-recycle 3 --rank plddt
--model-type alphafold2_ptm --random-seed 42
--msa-mode mmseqs2_uniref_env
```

`tools/msa.HARDWARE_FLAGS` (allowed, does not alter the model):
```
--disable-unified-memory
```

`tools/retry.DEFAULT_BACKOFF_SECONDS = (60, 300, 900)` — exactly as specified.
`tools/retry.DEFAULT_PER_CALL_TIMEOUT_SECONDS = 7200` — 2 h per attempt.

## Tests

`tools/tests/` — 43 tests, all passing.

| Test file | Count | What's covered |
|---|---|---|
| `test_seqio.py` | 13 | FASTA round-trip, blank-line handling, internal-whitespace stripping, AA validation, KRAS G12D point mutation, position-bound + WT-mismatch errors, window/kmer iteration |
| `test_retry.py` | 9 | Default schedule matches brief, success on first/Nth attempt, full-schedule exhaustion + final raise, non-retryable exception fast-path, SIGALRM timeout fires, timeout-then-success retry, custom short schedule |
| `test_plddt.py` | 8 | Live smoke fixtures: rank-1 finder, AF2 pLDDT gate ≥ 85, ESMFold heavy-atom and Cα convention monotonicity, fraction-below threshold, schema error paths |
| `test_rmsd.py` | 5 | Identity → 0, pure translation invariant under Kabsch, residue-subset restriction, **live AF2 vs ESMFold ubiquitin core RMSD < 3 Å**, missing-chain error |
| `test_mhc.py` | 8 | Canonical G12D 9-mer `VVGADGVGK` scores PREFERRED on both anchors, WT homolog same, poor-P9 / tolerated-P2 paths, length validation, full coverage of 9-mers spanning G12D, filter monotonicity, dict serialization |

Test runner config: `/mnt/d/AssignmentOPUS/pytest.ini`. WSL2 mount + redirected `TMPDIR` breaks pytest's default fd-capture (a `tmpfile.truncate()` lands on a path that's been GC'd mid-snap), so `--capture=no` is set globally. Helpers log via the stdlib `logging` module instead.

## Cross-stack RMSD sanity check (ColabFold AF2 vs ESMFold v1, ubiquitin)

| Region | n Cα | RMSD (Å) | Interpretation |
|---|---|---|---|
| Full 1–76 | 76 | **0.631** | Within experimental NMR ensemble variance |
| Core 1–71 (drops flexible diglycine tail) | 71 | **0.166** | Essentially identical fold |

This is unusually tight agreement — it confirms that:
1. The fp16/fp32 mixed-precision setup for ESMFold did not introduce numeric drift,
2. The two stacks see the same global topology,
3. The flexible C-terminal tail (Gly 73–76) accounts for almost the entire residual disagreement (jumping from 0.166 to 0.631 Å when the tail is included).

## What's NOT in Step 2

- No NetMHCpan / MHCflurry binding-affinity predictor. The brief explicitly limits us to GPU-resident, locally-installed tools; both NetMHCpan and MHCflurry require third-party licenses or remote APIs we have not approved. `tools/mhc.py` is a *motif* scorer, used to shortlist candidates before structural review.
- No structure-based binding-energy scorer. Step 3 will use Cα RMSD + per-residue pLDDT + pLDDT-band visualization as the proxy.

## Acceptance

- 43 / 43 unit tests green (≈3.4 s wall).
- Live cross-stack RMSD on the smoke fixture documented above.
- All brief-locked flags + retry/timeout values pinned in code, not in shell scripts.

Proceeding to Step 3 (core research workflow).
