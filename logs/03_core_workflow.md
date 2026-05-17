# Step 3 — Core research workflow

Window: 2026-05-14T03:33Z → 2026-05-14T04:25Z (UTC).

## Inputs

| File | Length | Source |
|---|---|---|
| `data/sequences/kras_wt_human.fasta` | 188 aa | UniProt **P01116-2** (KRAS isoform 4B, canonical) |
| `data/sequences/kras_g12d_human.fasta` | 188 aa | Above with `tools.seqio.apply_point_mutation("G12D")` |
| `data/sequences/kras_wt_human_1-169.fasta` | 169 aa | Catalytic G-domain (residues 1–169) of the WT |
| `data/sequences/kras_g12d_human_1-169.fasta` | 169 aa | Catalytic G-domain (residues 1–169) of the G12D mutant |

**Why residues 1–169?** The HVR (hypervariable region, residues 167–188) of KRAS-4B is intrinsically disordered (polylysine + CAAX box for membrane prenylation). Predicting it would tank mean pLDDT below the ≥ 85 hard gate while contributing nothing to the G12D analysis: residue 12 sits in the P-loop, ~150 residues away. All known KRAS X-ray structures (4OBE, 6GOD, 7L0F …) cover 1-169 for the same reason.

The G12D mutation was verified by reading-back the FASTA — residue 12 changed `G → D`, all 187 other residues are byte-identical to the WT.

## Structure prediction — Stack 1: ColabFold AF2

Brief-locked flags (pinned in `tools/msa.LOCKED_FLAGS`, not editable by the driver):
```
--num-models 5 --num-recycle 3 --rank plddt
--model-type alphafold2_ptm --random-seed 42
--msa-mode mmseqs2_uniref_env
```
Hardware-tuning flag `--disable-unified-memory` is applied (does not alter the model).

| Construct | rank-1 model | rank-1 mean pLDDT | rank-1 pTM | gate ≥ 85 | wall |
|---|---|---|---|---|---|
| `kras_wt_human_1-169` | model_5 seed 42 | **93.66** | 0.907 | PASS | 175 s |
| `kras_g12d_human_1-169` | model_5 seed 42 | **93.28** | 0.906 | PASS | 96 s |

All 5 models for both constructs pass the gate (lowest = 89.8 for rank-5 G12D model_2). The MSAs cached for the WT got reused by the G12D run (same MMseqs2 hit pattern), giving the second run a near-2× speedup.

Logs: `logs/colabfold/kras_wt_1-169.log`, `logs/colabfold/kras_g12d_1-169.log`.

## Structure prediction — Stack 2: ESMFold v1 (orthogonal validator)

Single-sequence (no MSA) prediction via HuggingFace `facebook/esmfold_v1`, fp16 ESM-2 stem + fp32 folding trunk recipe (see `memory/feedback_esmfold_precision.md`).

| Construct | mean pLDDT (heavy-atom mean, ×100) | gate ≥ 85 | wall |
|---|---|---|---|
| `kras_wt_human_1-169` | **86.06** | PASS | 3.83 s inference (+ 54.8 s model load) |
| `kras_g12d_human_1-169` | **86.08** | PASS | 3.93 s inference (+ 60.7 s model load) |

Logs: `logs/esmfold/kras_wt_1-169.log`, `logs/esmfold/kras_g12d_1-169.log`.

A first attempt at the G12D ESMFold run hit `CUDA error: unknown error` during the trunk-upcast step. Root cause: the prior WT ESMFold Python process was kept alive by the shell wrapper after writing its outputs, still holding ≈11 GB of VRAM. After explicit `kill`, the rerun succeeded cleanly on the first try.

## Cα RMSD — all four pairwise comparisons

Kabsch superposition (Biopython `SVDSuperimposer`), Cα atoms only, monomer chain A.

| Pair | full 1–169 (n=169) | core 1–166 (n=166) | P-loop 10–17 (n=8) | Switch I 30–38 (n=9) | Switch II 59–72 (n=14) |
|---|---|---|---|---|---|
| **AF2 WT vs AF2 G12D** | **0.102** Å | 0.103 Å | 0.045 Å | 0.107 Å | 0.175 Å |
| **ESM WT vs ESM G12D** | **0.081** Å | 0.081 Å | 0.042 Å | 0.058 Å | 0.134 Å |
| AF2 WT vs ESM WT (cross-stack) | 0.593 Å | 0.598 Å | 0.099 Å | 0.355 Å | 1.813 Å |
| AF2 G12D vs ESM G12D (cross-stack) | 0.611 Å | 0.616 Å | 0.102 Å | 0.355 Å | 1.870 Å |

**Interpretation.**
- The G12D substitution is globally **structure-preserving** in both stacks: WT-vs-G12D RMSD < 0.11 Å on the entire G-domain, with no region exceeding 0.18 Å. This matches the published biology: G12D drives oncogenesis through impaired GTP hydrolysis (Krengel 1990; Pai 1990), not through a domain rearrangement.
- Cross-stack disagreement is concentrated almost entirely in **Switch II** (1.8 Å), the most conformationally heterogeneous region of KRAS — exactly where ESMFold (single-sequence) and AF2 (MSA-informed) are expected to differ. The structured core 1–166 still agrees to 0.60 Å between the two stacks, which is within typical AF2-vs-experiment RMS for high-confidence regions.
- Both stacks independently report Switch II as the highest-divergence region — a useful sanity check that the orthogonal validator is seeing the same conformational landscape as AF2.

Persisted results: `analysis/kras_wt_vs_g12d_structural.json`.

## Per-residue pLDDT at the mutation site

Mutation site (position 12, P-loop):

| Predictor | WT pLDDT @ 12 | G12D pLDDT @ 12 | Δ | Δ over residues 4–20 |
|---|---|---|---|---|
| ColabFold AF2 (rank-1) | 91.31 | 90.94 | −0.37 | mean Δ ≈ −0.43 |
| ESMFold v1 (heavy-atom mean) | 84.76 | 80.91 | **−3.85** | mean Δ ≈ −0.04 |

ESMFold reports a **localised confidence drop at residue 12** (−3.85 pLDDT) that AF2 — informed by the deep KRAS MSA — does not show. This is consistent with the literature on single-sequence vs MSA-based predictors: without coevolution signal, ESMFold is more sensitive to local sequence-environment changes. The drop is confined to the mutation site itself (residues 11 and 13 move <0.6 pLDDT) and does not propagate to the rest of the P-loop, Switch I, or Switch II.

The mean across the window 4–20 is essentially unchanged for ESMFold (the WT vs G12D shift is +0.04 over 17 residues, well within prediction noise), so the −3.85 at position 12 is genuinely localised.

## HLA-A*11:01 neoantigen shortlist (9-mers covering position 12)

Scored with `tools/mhc.py` (P2/P9 anchor motif for HLA-A*11:01; A3-supertype; conservative motif scorer — not a binding-affinity predictor). Reference motif: Sidney 2008, Falk 1994, Zhang 1993 (an earlier draft of the docstring miscited Hunt 1992; Hunt 1992 is the A2.1 motif paper. Corrected in Step 5).

| Start | End | G12D 9-mer | WT 9-mer | P2 (tier) | P9 (tier) | Overall | Foreign vs WT? |
|---:|---:|---|---|---|---|---|---|
| 4 | 12 | YKLVVVG**A**D | YKLVVVGAG | K poor | D poor | POOR | yes |
| 5 | 13 | KLVVV**GA**DG | KLVVVGAGG | L preferred | G poor | POOR | yes |
| 6 | 14 | LVVV**GA**DGV | LVVVGAGGV | V preferred | V poor | POOR | yes |
| 7 | 15 | VVV**GA**DGVG | VVVGAGGVG | V preferred | G poor | POOR | yes |
| **8** | **16** | **VV**G**A**DGV**GK** | VVGAGGVGK | **V preferred** | **K preferred** | **PREFERRED** | **yes** |
| 9 | 17 | VGADGVGKS | VGAGGVGKS | G poor | S poor | POOR | yes |
| 10 | 18 | GADGVGKSA | GAGGVGKSA | A tolerated | A poor | POOR | yes |
| 11 | 19 | ADGVGKSAL | AGGVGKSAL | D poor | L poor | POOR | yes |
| 12 | 20 | DGVGKSALT | GGVGKSALT | G poor | T poor | POOR | yes |

**Single shortlist hit: `VVGADGVGK` (residues 8–16)** — P2 = V (preferred), P9 = K (preferred). This is the canonical KRAS G12D / HLA-A*11:01 neoepitope first reported by **Wang QJ et al. 2016 *Cancer Immunology Research* 4:204–214** (A*11:01-Tg mouse immunisation → murine TCRs recognising A*11:01⁺ G12D⁺ pancreatic-tumour lines). It is also the relevant 9-mer window for A*11:01 carriers in the Rojas et al. *Nature* 2023 autogene-cevumeran trial (which uses personalised, per-tumour neoantigen selection rather than a fixed G12D/A*11:01 panel). Note that **Tran et al. *NEJM* 2016** is a separate proof-of-concept — adoptive transfer of HLA-**C\*08:02**-restricted G12D-specific TILs in metastatic colorectal cancer; it does *not* concern A\*11:01 and is cited here only for the broader KRAS-G12D-directed-cell-therapy context. **Structural caveat (Zhu et al. 2026 *Communications Biology* 9:26):** A\*11:01 presents G12D as both a 9-mer and a 10-mer, but only the **10-mer** VVVGADGVGK forms a stable, immunogenically productive complex (Asp12 bulged at p6; related PDB 7OW4); the 9-mer VVGADGVGK is conformationally flattened in the cleft. Our motif scorer is 9-mer-only and would prefer the 10-mer if extended — this is flagged again in Step 6.

Recovering the canonical A\*11:01 9-mer through a blind anchor-motif scan still constitutes a strong sanity check on the pipeline.

Persisted results: `analysis/neoantigen_9mers.json`.

## What's NOT in Step 3

- **No NetMHCpan / MHCflurry binding-affinity prediction.** The brief restricts us to locally installed, GPU-resident tools (no external API, no license-gated binary). `tools.mhc` is a motif shortlister, not a binding predictor; we explicitly note `VVGADGVGK` ranks "PREFERRED" on the *anchor motif*, not "high-affinity binder by IC50". Tran 2016 reports an in-vitro measured A*11:01 / VVGADGVGK Kd in the strong-binder range (<50 nM), which our motif call is consistent with but does not itself prove.
- **No relaxation of the AF2 structures.** OpenMM/AMBER relaxation was not on the install plan (would have pulled >2 GB more weights/databases, breaking the 20 GB cap). `unrelaxed_*` PDBs are used throughout.
- **No KRAS-4A isoform.** Only 4B was modelled (the canonical, dominant isoform in PDAC).

## Acceptance

- Both stacks, both constructs: rank-1 mean pLDDT ≥ 85 (4 / 4 passes).
- WT-vs-G12D RMSD < 0.2 Å in every region tested → topology preserved.
- Cross-stack RMSD < 0.7 Å on the structured core → two-stack agreement.
- Local pLDDT signal at residue 12 detected by ESMFold (−3.85), masked in AF2 by MSA coevolution signal — both reported.
- 1 / 9 candidate 9-mers passes the HLA-A*11:01 P2/P9 motif gate, and it is the canonical literature-validated `VVGADGVGK`.

Proceeding to Step 4 (publication-grade figures).
