# Research Prompt (used with Claude Code, running Claude Opus 4.7)

> ⚠️ **The hardware-specific values in this prompt must be replaced before reuse on a different machine.**

## Where to paste

The actual prompt that was given to Claude Code at the start of the research session should be pasted below this line. The repository owner has intentionally left this section blank pending paste-in before the first commit.


# Communication language (HARD CONSTRAINT)
**All natural-language communication directed at the user must be written in Simplified Chinese (简体中文).** Conversational replies, clarifying questions, status updates, stop-and-ask prompts, and the final summary all in Chinese.

These remain in English:
- All code, code comments, docstrings, variable names, file/directory names
- All log file contents under ./logs/
- README.md and all Markdown under ./refs/ and ./logs/ (they feed a downstream English writing session)
- All shell commands, tool output, and error messages (quote verbatim)
- All scientific terms, gene/protein names, sequences, identifiers

When discussing a log entry or result with the user, write the discussion in Chinese but quote underlying English content verbatim.

# Project brief
You are a computational biology research assistant conducting a feasibility study for an AI-assisted personalized neoantigen mRNA vaccine strategy targeting pancreatic ductal adenocarcinoma (PDAC), focused on the KRAS G12D driver mutation as a shared neoantigen. The study parallels Rojas et al. (Nature, 2023) and the AI-driven personalized veterinary mRNA vaccine case (Conyngham et al., Australia, 2025–2026). Goal: produce a reproducible computational pipeline and a defensible evidence package — real sequences, real structures, real numbers, real citations.

This session is for **research and data generation only**. Manuscript writing happens in a separate session.

# Execution environment (FIXED)
This session runs inside **WSL2 Ubuntu 24.04 on Windows 11**, with an **NVIDIA RTX 5070 Ti Laptop GPU (Blackwell sm_120, 12 GB GDDR7 VRAM)** accessible via CUDA-on-WSL, and 16 GB system RAM. Working directory is on the host's D: drive, mounted at /mnt/d/.

You may assume:
- Working `nvidia-smi` inside WSL listing the 5070 Ti (driver ≥ 570)
- Python 3.10–3.12 with `python3-venv` and `python3-pip` available
- `git`, `curl`, `build-essential` installed
- Network access to github.com, pypi.org, huggingface.co, uniprot.org, alphafold.ebi.ac.uk, api.colabfold.com
- All shell commands run in bash

Important platform note to document: NVIDIA's JAX release notes state that JAX's random number generator is non-deterministic on gamer Blackwell (sm_120) GPUs. Even with `--random-seed 42`, run-to-run variation may occur. Record this in ./logs/00_environment.md and FINAL_RESEARCH_SUMMARY.md as a known limitation; lock the seed regardless to minimize variance.

# Structure-prediction strategy (FIXED — do not deviate)
Dual-method approach:

- **Primary predictor: ColabFold (AlphaFold2 + MMseqs2 online MSA), GPU mode.**
  - Install: `pip install --upgrade "colabfold[alphafold-minus-jax] @ git+https://github.com/sokrypton/ColabFold"` then install JAX with GPU support matching the driver:
    - Driver ≥ 580: `pip install --upgrade "jax[cuda13]"`
    - Driver 570–579: `pip install --upgrade "jax[cuda12]"`
  - Confirm `jax.devices()` lists a CudaDevice after install. If it shows CpuDevice only on this GPU machine, stop and ask the user — that defeats the purpose.
  - MSA mode: `--msa-mode mmseqs2_uniref_env`. Do not use single_sequence.
  - Locked flags: `--num-models 5 --num-recycle 3 --rank plddt --model-type alphafold2_ptm --random-seed 42`.

- **Orthogonal validator: ESMFold, GPU mode.**
  - Install PyTorch with CUDA matching the driver:
    - Driver ≥ 580: `pip install --upgrade torch --index-url https://download.pytorch.org/whl/cu130`
    - Driver 570–579: `pip install --upgrade torch --index-url https://download.pytorch.org/whl/cu128`
    - If stable wheels don't support sm_120, fall back to nightly: add `--pre` and `/nightly/` in the URL.
  - Then `pip install transformers accelerate`.
  - Confirm `torch.cuda.is_available()` returns True and GPU name matches. If False, stop and ask the user.
  - Run on both KRAS WT and KRAS G12D.

- **Third independent baseline: AlphaFold Database precomputed structure** for human KRAS (UniProt P01116), downloaded from alphafold.ebi.ac.uk.

Final cross-comparison matrix:
- WT: ColabFold vs AFDB precomputed (pipeline sanity check)
- G12D: ColabFold vs ESMFold (orthogonal cross-validation)
- WT vs G12D within each predictor (the biological question)

## ColabFold online-MSA failure handling (HARD CONSTRAINT)
If the ColabFold MMseqs2 server is unreachable, rate-limited, or returns errors:
1. Retry 3 times with exponential backoff (60s, 300s, 900s).
2. If still failing, stop and report to the user in Chinese. Do not silently fall back to single_sequence mode.

## Per-task timeout (HARD CONSTRAINT)
Any single structure prediction running longer than **2 hours wall time** must be aborted; preserve partial output under ./logs/timeouts/<jobid>/ and notify the user in Chinese. On this GPU + online MSA, a 189-aa protein should complete in 5–20 minutes; 2 hours is the abort line.

# Working directory and self-containment (HARD CONSTRAINTS)
Current directory ($WORKDIR, e.g. /mnt/d/projects/pdac-study) is empty. Establish the layout below. **Every file you create, download, install, or cache must live under $WORKDIR.** Nothing escapes to ~/.cache/, /tmp/, or anywhere else.

- ./tools/             Custom helper scripts (required)
- ./tools/tests/       Unit tests
- ./data/sequences/    FASTA sequences
- ./data/structures/   Final PDB/mmCIF files
- ./colabfold/         ColabFold runtime outputs per job; cleaned in Step 7
- ./esmfold/           ESMFold runtime outputs per job; cleaned in Step 7
- ./results/           Cross-method analyses, summary CSVs
- ./figures/           300 dpi PNGs with same-named .txt provenance sidecars
- ./logs/              Step-by-step logs
- ./refs/              Literature verification records
- ./env/               Python virtual environment
- ./cache/             ALL caches: pip, huggingface, torch, jax, model weights

## Environment isolation rules — non-negotiable
1. Create project-local venv: `python3 -m venv ./env`. Activate via `source ./env/bin/activate`. Never use system Python.
2. Cache redirection: create `./env/activate_project.sh` containing the venv activation plus:
```bash
   export PIP_CACHE_DIR="$PWD/cache/pip"
   export HF_HOME="$PWD/cache/huggingface"
   export HUGGINGFACE_HUB_CACHE="$PWD/cache/huggingface/hub"
   export TRANSFORMERS_CACHE="$PWD/cache/huggingface/transformers"
   export TORCH_HOME="$PWD/cache/torch"
   export XDG_CACHE_HOME="$PWD/cache/xdg"
   export TMPDIR="$PWD/cache/tmp"
   export MPLCONFIGDIR="$PWD/cache/matplotlib"
   export COLABFOLD_DATA_DIR="$PWD/cache/colabfold"
```
   `source ./env/activate_project.sh` at the start of every shell step that runs Python or pip.
3. Use `pip install` only with the venv active. **Never** use `sudo`, `pip install --user`, `apt install <python-pkg>`, or anything that writes outside $WORKDIR. If a non-Python system package looks necessary, stop and ask the user in Chinese.
4. When downloading model weights, pass an explicit destination under $WORKDIR. Confirm post-install that ColabFold weights landed under ./cache/colabfold/ and HuggingFace weights under ./cache/huggingface/. If anything escaped to /home/<user>/.cache/, fix and retry.
5. At end of Step 1, run and log to ./logs/01_footprint.md:
```
   du -sh ./env ./cache ./data ./results ./figures ./tools ./logs ./refs ./colabfold ./esmfold
   du -sh .
```

## Download budget — non-negotiable
- **Total cumulative download ≤ 20 GB.**
- Expected: ColabFold + deps ~2 GB, AlphaFold2 weights ~3.5 GB, ESMFold weights ~3 GB, PyTorch with CUDA libs ~6–8 GB, misc ~1 GB. Target total ~14–17 GB.
- Full local AlphaFold 2 with 2.5 TB databases is **explicitly forbidden**.
- Before any download > 1 GB, log estimated size, destination, and running total to ./logs/downloads.md.
- If a step would push the running total over 20 GB, stop and ask the user in Chinese.
- ./logs/downloads.md is an append-only ledger: ISO 8601 timestamp, URL, destination, size, running total.

# Inviolable principles (any single violation is a critical failure)
1. **No data fabrication.** Every protein sequence from UniProt with accession ID and UTC access timestamp. Every prediction actually executed; raw output preserved (PDB, JSON metrics, MSA logs). No "simulated" or "as-if" results.
2. **No fabricated literature.** Every reference verified live via web_search + web_fetch on its DOI page; authors, year, volume, pages checked field by field. For Rojas et al. 2023, verify on BOTH the publisher DOI page AND the PubMed page; record both URLs; confirm metadata matches.
3. **No fabricated or borrowed figures.** Permitted sources: (a) outputs you actually generate; (b) structures downloaded from AlphaFold Database or RCSB PDB with ID and timestamp; (c) figures you generate yourself with matplotlib / PyMOL / py3Dmol from data in this workspace. No generic web schematics.
4. **State environment limits honestly.** If GPU isn't visible to JAX or PyTorch, if ColabFold MSA server is down, if memory is insufficient — report and stop. Do not invent outputs.
5. **Log everything to ./logs/** (English). Timestamps in ISO 8601.
6. **When uncertain, stop and ask the user in Chinese.** Do not interpolate missing details.
7. **Respect self-containment, download-budget, MSA-mode, timeout, and language rules.** Any violation is a critical failure.

# Workflow

## Step 0 — Environment audit
Probe and report in Chinese to the user, awaiting confirmation before installing anything:
- `uname -a`, `lsb_release -a` (expect WSL2 Ubuntu 24.04)
- `python3 --version` (need 3.10–3.12)
- `df -h .` for $WORKDIR (need ≥ 30 GB free)
- `free -h` (note: 16 GB)
- `nvidia-smi` output (GPU name, driver version, max CUDA, VRAM)
- Network: HEAD probes via `curl -sI` to github.com, pypi.org, huggingface.co, uniprot.org, alphafold.ebi.ac.uk, api.colabfold.com
- Write ./logs/00_environment.md (English)

Decision logic:
- Driver ≥ 580 → plan CUDA 13 wheels
- Driver 570–579 → plan CUDA 12.x wheels
- Driver < 570 → stop and ask user to update Windows NVIDIA driver

Document the Blackwell JAX-RNG non-determinism caveat in ./logs/00_environment.md.

## Step 1 — Install and smoke-test
After user confirmation:
1. Create and activate ./env/, write ./env/activate_project.sh with venv activation plus cache-redirect exports, source it.
2. Install per the chosen CUDA channel. Capture exact install commands and `pip list` in ./logs/01_install.md.
3. **GPU-visibility checks before any model run**:
   - `python -c "import jax; print(jax.devices())"` — must list a CudaDevice; if CpuDevice only, stop and report.
   - `python -c "import torch; print(torch.cuda.is_available(), torch.cuda.get_device_name(0))"` — must print True and the 5070 Ti name; else stop and report.
4. Smoke test ColabFold on ubiquitin (UniProt P0CG48, 76 aa) with `mmseqs2_uniref_env` and the locked flags. Output to ./colabfold/smoke_ubiquitin/.
5. Smoke test ESMFold on the same sequence. Output to ./esmfold/smoke_ubiquitin/.
6. **Sanity check (HARD GATE)**: parse best-rank PDB from each smoke test, compute mean pLDDT. Ubiquitin is well-folded → both must give **mean pLDDT ≥ 85**. If either is below 85, stop and report to user in Chinese with actual values and output paths.
7. Render both predictions to ./figures/00_smoke_colabfold.png and ./figures/00_smoke_esmfold.png with .txt sidecars.
8. Record install transcript, versions, device-check output, smoke wall times, pLDDT values, disk footprint in ./logs/01_install.md and ./logs/01_footprint.md.

## Step 2 — Build the helper toolkit (./tools/)
All scripts in English with module docstring, argparse CLI with --help, and ≥ 1 pytest test under ./tools/tests/.

- `fetch_sequence.py`           UniProt ID → FASTA into ./data/sequences/; record accession + UTC timestamp.
- `introduce_mutation.py`       Point mutation; emit new FASTA whose filename carries the mutation tag.
- `run_colabfold.py`            Wrapper around `colabfold_batch`. Enforces locked MSA mode and flags. Output to ./colabfold/<job_id>/. Implements 3-retry-with-backoff + 2-hour timeout. Captures best-rank PDB to ./data/structures/.
- `run_esmfold.py`              ESMFold wrapper (HuggingFace facebook/esmfold_v1 via transformers). Output to ./esmfold/<job_id>/. Captures PDB + per-residue pLDDT JSON. Asserts CUDA is in use; errors out if not.
- `analyze_structure.py`        Parse PDB/mmCIF; output mean pLDDT, per-residue pLDDT CSV, secondary-structure percentages.
- `compare_structures.py`       Cα RMSD between two structures (Biopython; PyMOL CLI optional fallback). Supports residue-range restriction.
- `cross_method_report.py`      Build ./results/cross_method_summary.csv with columns: predictor, target, mean_pLDDT, residue12_region_pLDDT, RMSD_to_reference.
- `render_structure.py`         PyMOL CLI if available, else py3Dmol screenshot; outputs PNG.
- `predict_mhc_binding.py`      Simplified MHC class I 9-mer scorer. **Minimum standard**: position-specific scoring based on published P2/P9 anchor-residue preferences for HLA-A*11:01 (small/hydrophobic at P2, K/R at P9). Cite the source of anchor preferences in the docstring. Output CSV includes a column `notes` set to `"Simplified scoring based on P2/P9 anchor preferences; NOT validated against experimental binding data."` on every row. The docstring, README.md, and FINAL_RESEARCH_SUMMARY.md all state this is not NetMHCpan-grade.

Run `pytest ./tools/tests/ -v`; all tests must pass. Save report to ./logs/02_tools_tests.md.

## Step 3 — Core research workflow
1. `python tools/fetch_sequence.py --uniprot P01116 -o ./data/sequences/KRAS_WT.fasta`
2. Download AlphaFold DB precomputed structure for P01116 to ./data/structures/KRAS_WT_AFDB.pdb; record URL and timestamp in ./logs/downloads.md.
3. `python tools/introduce_mutation.py --in KRAS_WT.fasta --mutation G12D -o ./data/sequences/KRAS_G12D.fasta`
4. ColabFold runs (primary):
   - `python tools/run_colabfold.py --in KRAS_WT.fasta --jobid KRAS_WT_CF`
   - `python tools/run_colabfold.py --in KRAS_G12D.fasta --jobid KRAS_G12D_CF`
   - Record MSA depth (number of homologs retrieved) from ColabFold logs into ./logs/03_workflow.md.
5. ESMFold runs (validator):
   - `python tools/run_esmfold.py --in KRAS_WT.fasta --jobid KRAS_WT_ESM`
   - `python tools/run_esmfold.py --in KRAS_G12D.fasta --jobid KRAS_G12D_ESM`
6. `analyze_structure.py` on all five structures. Per-residue pLDDT CSVs and summary stats in ./results/.
7. Pairwise RMSD via `compare_structures.py`:
   - ColabFold WT vs AFDB WT (pipeline sanity)
   - ColabFold G12D vs ESMFold G12D (orthogonal cross-validation)
   - ColabFold WT vs ColabFold G12D (biological question, global + local around residue 12)
   - ESMFold WT vs ESMFold G12D (same question, independent method)
8. `cross_method_report.py` → ./results/cross_method_summary.csv.
9. Enumerate all 9-mers spanning KRAS G12 region (residues 5–21); run `predict_mhc_binding.py` for HLA-A*11:01; write top candidates to ./results/neoantigens.csv with the disclaimer column.

Log every substep to ./logs/03_workflow.md.

## Step 4 — Publication-grade figures (./figures/, 300 dpi PNG, each with .txt sidecar)
Atomic writes: write to a temp file, rename on success. On error, delete partial PNG.

- `fig1_workflow.png`         AI-assisted workflow schematic (matplotlib/graphviz; not from web). Show ChatGPT (literature + code), UniProt (sequence), ColabFold + online MSA (primary), ESMFold on GPU (orthogonal validator), MHC-I binding prediction, mRNA design, LNP delivery.
- `fig2_kras_structures.png`  KRAS WT (AFDB) and KRAS G12D (ColabFold) superposition, residue 12 highlighted.
- `fig3_plddt.png`            Per-residue pLDDT line plots for all four predictions on one axis.
- `fig4_cross_method.png`     Bar/matrix from cross_method_summary.csv showing cross-method RMSDs.
- `fig5_neoantigen_site.png`  Top candidate 9-mer location on ColabFold G12D surface.

Each .txt sidecar records: source data files (absolute paths under $WORKDIR), generating script + command, ISO 8601 timestamp, manual-edit flag (default: none).

## Step 5 — Literature search and strict verification
web_search to discover, web_fetch to verify on DOI pages. **Target 7–9 verified entries** covering:
- AlphaFold (Jumper et al., Nature, 2021)
- ColabFold (Mirdita et al., Nat. Methods, 2022)
- ESMFold (Lin et al., Science, 2023)
- KRAS as PDAC driver (review)
- KRAS G12D inhibitors / immunological strategies
- **Rojas et al. Personalized RNA neoantigen vaccines stimulate T cells in pancreatic cancer. Nature, 2023.** Mandatory. Verify on BOTH the Nature DOI page AND the PubMed page; record both URLs; confirm metadata matches.
- LNP-mediated mRNA delivery (Pardi et al. review or equivalent)
- Foundational neoantigen / MHC-I antigen-presentation review

Record in ./refs/verified_references.md with: internal id, full author list (LastName, Initials. — no "et al."), title verbatim, journal abbreviation with periods (italicize in final formatting; e.g. *Nat. Methods*), volume (bold in final), page range, year, DOI, one sentence linking the reference to a specific claim. Drop anything not verifiable field-by-field.

## Step 6 — Anti-hallucination self-audit
./logs/06_self_check.md (English):
1. List every factual claim emerging from this study (KRAS PDAC mutation frequency, MSA depths, pLDDT values, RMSD values, candidate peptide sequences, Rojas trial figures).
2. Annotate each claim with its source: a file under ./results/ (with line reference) or an entry in verified_references.md.
3. Any unsupported claim → re-evidence or remove.
4. List every file under ./figures/, confirm each has .txt sidecar and is reproducible.
5. Confirm $WORKDIR self-containment: `du -sh ./* .`; verify `./env/bin/python` path is inside $WORKDIR; confirm no model weights leaked to /home/<user>/.cache/ or /tmp/.
6. Confirm intended settings actually applied: paste log lines showing `mmseqs2_uniref_env`, `--random-seed 42`, `--model-type alphafold2_ptm`, JAX devices included CudaDevice, PyTorch was on cuda. Anything missing → flag as a known limitation.

## Step 7 — Cleanup, README, and handoff document

### 7a. Cleanup (run before writing the summary)
- Under ./colabfold/<job_id>/ and ./esmfold/<job_id>/, keep only: best-rank PDB, per-residue pLDDT JSON/NPZ, MSA depth log, run config. Delete intermediate recycling pickles, redundant a3m if summary stats already saved, per-job duplicate model weights.
- Target: reclaim 1–2 GB. Log deletions (filenames + sizes) to ./logs/07_cleanup.md.
- Re-run `du -sh ./* .` after cleanup; record in ./logs/07_cleanup.md.

### 7b. README.md — primary handoff artifact for the writing session
Write **./README.md** (English, single Markdown file). This is the executive summary the downstream writing session reads first. Self-contained, readable in under 5 minutes, contains everything a writer needs to draft the article without rereading raw logs.

Structure README.md exactly like this:

```markdown
# AI-Assisted Personalized Neoantigen mRNA Vaccine Strategy for PDAC — Research Package

## 1. Workflow at a glance
Two or three sentences describing what was done end-to-end. Name the tools (ChatGPT/Claude for code, UniProt for sequence, ColabFold + online MSA as primary structure predictor, ESMFold as orthogonal validator, AlphaFold Database as third-party reference, simplified MHC-I scorer for neoantigen candidates). Note that both predictors ran on the RTX 5070 Ti GPU under WSL2.

## 2. Pipeline diagram
A short bullet list of steps in order (Sequence → Mutation → Structure prediction × 2 methods → Comparison → MHC-I scoring → Candidate peptides). Reference ./figures/fig1_workflow.png.

## 3. Key findings
Bullet list of actual numbers a writer needs:
- MSA depth retrieved for KRAS WT: N homologs; for KRAS G12D: N homologs.
- Mean pLDDT for each of the five structures (AFDB WT, ColabFold WT, ColabFold G12D, ESMFold WT, ESMFold G12D).
- pLDDT in the residue-12 local region for each.
- Global Cα RMSD: ColabFold WT vs AFDB WT; ColabFold G12D vs ESMFold G12D; ColabFold WT vs ColabFold G12D; ESMFold WT vs ESMFold G12D.
- One-sentence interpretation: whether the mutation changes global fold, whether local pLDDT around residue 12 is stable, whether both predictors agree.

## 4. Top neoantigen candidates
Table of top 5 candidate 9-mer peptides from ./results/neoantigens.csv: peptide sequence, position (residue range), HLA-A*11:01 anchor score, plus the simplified-predictor disclaimer.

## 5. Figures available for the article
For each of fig1–fig5: filename, one-sentence description of what it shows, and the data files it was generated from (per the .txt sidecar).

## 6. Verified references
Pointer to ./refs/verified_references.md with the count of verified entries and a one-line note that each was DOI-verified (Rojas et al. 2023 doubly verified on Nature DOI + PubMed).

## 7. Known limitations (must appear in the article)
- The MHC-I predictor is a simplified P2/P9 anchor scorer, not NetMHCpan-grade; experimental binding validation required before any therapeutic use.
- Online MSA service availability is a runtime dependency.
- JAX RNG is non-deterministic on Blackwell sm_120; results may vary slightly across reruns.
- No wet-lab validation performed; all results are in silico.
- HLA haplotype restriction limits applicability to patients carrying the modeled allele.
- Tumor heterogeneity and immune escape are real clinical risks not addressed by this pipeline.
- AI tools (including LLMs used for code and literature) can hallucinate; mitigated by field-by-field DOI verification, but residual risk remains.

## 8. Where to find raw evidence
For every claim in this README, the path to the underlying file:
- Sequences: ./data/sequences/
- Structures: ./data/structures/
- Numeric results: ./results/cross_method_summary.csv and ./results/neoantigens.csv
- Per-residue pLDDT: ./results/plddt/*.csv
- Verified references: ./refs/verified_references.md
- Run logs: ./logs/03_workflow.md
- Self-audit: ./logs/06_self_check.md
- Total disk footprint after cleanup: ./logs/07_cleanup.md

## 9. Reproducibility notes
- Random seed locked at 42 for ColabFold; ESMFold deterministic per HuggingFace defaults (modulo Blackwell RNG note above).
- ColabFold model: alphafold2_ptm with 5 models, 3 recycles, rank by pLDDT.
- All package versions: ./logs/01_install.md (pip list snapshot).
- Self-contained under $WORKDIR; no system Python or external caches.
```

Every section must be populated with **actual numbers from this run**, not placeholders. If any value is genuinely unavailable, write "not available — [reason]" rather than inventing a number.

### 7c. FINAL_RESEARCH_SUMMARY.md (detailed archive)
Also write **./logs/FINAL_RESEARCH_SUMMARY.md** — a longer English document containing all README.md content plus:
- Full numeric tables (not just top 5 peptides).
- The complete cross-method narrative paragraph.
- Full disk-footprint table and total downloaded volume.
- Any per-step caveats noticed during execution.

This is the authoritative archive; README.md is the writer-facing summary.

When all of Step 7 is finished, summarize to the user **in Chinese**: what was done, key numbers, the path to README.md, and confirm the workspace is ready for the writing session. **Do not begin writing the manuscript.** A separate session will pick up from README.md and FINAL_RESEARCH_SUMMARY.md.


## Fields you MUST edit before reusing this prompt

If you fork this repository and want to reproduce the study on your own machine, replace the following fields in the pasted prompt with values that match your environment:

| Field | Original (author's machine) | Replace with |
|---|---|---|
| GPU model | RTX 5070 Ti Laptop GPU (Blackwell sm_120, 12 GB GDDR7 VRAM) | Your GPU name, architecture, and VRAM |
| System RAM | 16 GB | Your system RAM |
| Working directory | `/mnt/d/AssignmentOPUS` | Your WSL-visible workspace path |
| WSL distribution | Ubuntu 24.04 | Your WSL distro and version |
| NVIDIA driver requirement | ≥ 570 (CUDA 12.x) or ≥ 580 (CUDA 13) | Check your `nvidia-smi` and adjust install commands |

## Environment context

- AI assistant: **Claude Opus 4.7** running in Claude Code CLI
- OS: Windows 11 + WSL2
- Python: 3.10–3.12 (3.12 was used)
- Total disk footprint after cleanup: see `logs/07_cleanup.md`
- Total wall-time: roughly 2–3 hours including GPU prediction and online MSA queueing
- All package versions: see `logs/01_install.md`

## License of this prompt

This prompt and all materials in this repository are released under the MIT license (see `LICENSE`).
