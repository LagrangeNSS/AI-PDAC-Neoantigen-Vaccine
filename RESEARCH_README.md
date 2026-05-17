# KRAS G12D PDAC neoantigen feasibility study

A reproducible, two-stack structure-prediction + anchor-motif feasibility
study for an HLA-A\*11:01-restricted KRAS G12D neoantigen mRNA-vaccine
pipeline, modelled on Rojas et al. *Nature* 2023 (autogene cevumeran).

All work is contained under `/mnt/d/AssignmentOPUS/`; nothing escapes
WORKDIR to `~/.cache`, `/tmp`, etc.

## Quickstart

```bash
# Activate the venv with redirected caches.
source env/activate_project.sh

# Re-run unit tests (43 tests; ~3 s wall).
pytest

# Re-run the Step-3 core analysis from the existing predictions.
python -m tools.analyze_kras_step3

# Re-render Step-4 figures.
python -m tools.make_step4_figures
```

To re-run the structure predictions themselves (~5–7 min each, GPU):

```bash
# ColabFold AF2 (brief-locked flags pinned in tools/msa.LOCKED_FLAGS).
python -m tools.run_kras_colabfold \
    --fasta data/sequences/kras_wt_human_1-169.fasta \
    --out-dir colabfold/kras_wt_1-169 \
    --log-path logs/colabfold/kras_wt_1-169.log

# ESMFold (single-sequence; fp16 stem + fp32 trunk recipe).
python -m tools.esmfold_run \
    --fasta data/sequences/kras_wt_human_1-169.fasta \
    --out-dir esmfold/kras_wt_1-169
```

## Repository layout

| Path | Contents |
|---|---|
| `data/sequences/` | FASTA inputs — KRAS-4B WT and G12D, full 188 aa and G-domain 1–169 |
| `tools/` | Helper modules (`seqio`, `retry`, `plddt`, `rmsd`, `mhc`, `msa`, `esmfold_run`) + Step-3/4 drivers + 43 unit tests under `tools/tests/` |
| `colabfold/` | ColabFold AF2 outputs for both constructs (rank-1 PDB + score JSON + MSA a3m + PAE PNG) |
| `esmfold/` | ESMFold v1 outputs (PDB with pLDDT in B-factor + per-record metrics JSON) |
| `analysis/` | `kras_wt_vs_g12d_structural.json` (RMSD + per-residue pLDDT) and `neoantigen_9mers.json` (HLA-A\*11:01 anchor shortlist) |
| `figures/` | 4 publication figures with `.txt` provenance sidecars |
| `logs/` | Step-by-step transcripts (`00`–`06`) + download ledger + per-run tool logs |
| `refs/` | `citations.bib` (14 verified entries) + `notes.md` (Step-5 verification log) |
| `cache/` | Tooling caches (HF, ColabFold params, JAX kernel cache, matplotlib font cache); fully redirected via `env/activate_project.sh` |
| `env/` | Project Python 3.12 venv |

## Brief-locked policy enforced in code

Pinned in `tools/msa.LOCKED_FLAGS`, not editable from any driver:

```
--num-models 5 --num-recycle 3 --rank plddt
--model-type alphafold2_ptm --random-seed 42
--msa-mode mmseqs2_uniref_env
```

Hardware-tuning flag `--disable-unified-memory` is allowed (does not
alter the AF2 model). Retry schedule: 60 / 300 / 900 s exponential
backoff, 2 h SIGALRM per attempt (`tools/retry.py`). Mean-pLDDT
hard gate: ≥ 85.

## Results at a glance

| Construct | Stack | Rank-1 mean pLDDT | Gate ≥ 85 |
|---|---|---|---|
| KRAS WT 1–169 | ColabFold AF2 | 93.66 | ✅ |
| KRAS G12D 1–169 | ColabFold AF2 | 93.28 | ✅ |
| KRAS WT 1–169 | ESMFold v1 (heavy-atom mean) | 86.06 | ✅ |
| KRAS G12D 1–169 | ESMFold v1 (heavy-atom mean) | 86.08 | ✅ |

WT-vs-G12D Cα RMSD < 0.11 Å on the full G-domain in both stacks (fold
preserved). HLA-A\*11:01 anchor shortlist: 1 / 9 candidate 9-mers
passes (`VVGADGVGK`, residues 8–16) — the canonical Wang 2016 *CIR*
A\*11:01 / G12D neoepitope.

## Provenance

Every figure in `figures/` has a paired `.txt` sidecar listing input
paths and key numbers. Every citation in `refs/citations.bib` was
verified against the publisher / PubMed / PMC record in Step 5; the
verification log lives at `refs/notes.md`. Three citation errors
caught and corrected during Step 5 / Step 6 are documented inline in
`refs/notes.md`.

Read the step transcripts in order:
1. `logs/00_environment.md` — host audit + caveats
2. `logs/01_install.md` — venv setup + smoke tests + bug log
3. `logs/01_footprint.md` — disk-budget accounting
4. `logs/02_helper_toolkit.md` — module inventory + test results
5. `logs/03_core_workflow.md` — predictions + RMSD + shortlist
6. `logs/04_figures.md` — figure rationale
7. `logs/06_self_audit.md` — anti-hallucination checks + caveats
8. `FINAL_RESEARCH_SUMMARY.md` — Chinese-language abstract-level summary

## Hard caveats (see `logs/06_self_audit.md` §D for the full list)

- The HLA-A\*11:01 scorer in `tools/mhc.py` is a **P2/P9 anchor motif
  shortlister**, not an IC50 binding predictor.
- The productive HLA-A\*11:01 / KRAS G12D epitope is the **10-mer**
  VVVGADGVGK (Zhu et al. 2026 *Commun Biol* 9:26); our 9-mer scan
  identifies the correct window but cannot rank 9-mer vs 10-mer.
- Predictions use the **G-domain 1–169** only (HVR 167–188 is
  intrinsically disordered and would fail the gate without informing
  the G12D analysis).
- No OpenMM relaxation (would have exceeded the 20 GB download budget);
  `unrelaxed_rank_001_*.pdb` used throughout.
- JAX RNG is non-deterministic on consumer Blackwell GPUs; per-residue
  pLDDT may drift ~0.1 between reruns. Mean pLDDT margin is large
  enough that the gate is unaffected.
