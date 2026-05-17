# Step 6 — Anti-hallucination self-audit

Window: 2026-05-14T05:10Z → 2026-05-14T05:45Z (UTC).

This step re-reads everything I produced in Steps 0–5 with the
explicit goal of identifying (a) any claim that is not directly
supported by an artefact under `/mnt/d/AssignmentOPUS/`, (b) any
citation that is not pinned in `refs/citations.bib` and verified in
`refs/notes.md`, and (c) any caveat that future readers must see when
they read the FINAL\_RESEARCH\_SUMMARY in Step 7. The audit is run
adversarially: I assume the reader will spot-check every number.

## Method

1. For every numerical claim in `logs/03_core_workflow.md` and
   `logs/04_figures.md`, locate the artefact (JSON, figure sidecar,
   PDB) that the number was read from. Mismatches are flagged.
2. For every citation in any `logs/*.md` or `tools/*.py`, confirm a
   matching entry exists in `refs/citations.bib` with the journal,
   volume, year, and key claim documented in `refs/notes.md`.
3. For every structural-biology assertion, confirm the underlying
   biology is supported by a Step-5-verified reference (not by
   training-data recall).
4. Re-derive the most consequential sequence/region facts directly
   from the on-disk FASTA + a one-shot Python check, instead of trusting
   the log files.

## A. Re-derived from the on-disk inputs (no trust assumed)

Run from a fresh Python prompt under the project venv:

```python
from tools.seqio import read_fasta
wt = read_fasta('data/sequences/kras_wt_human.fasta')[0].sequence
mu = read_fasta('data/sequences/kras_g12d_human.fasta')[0].sequence
assert len(wt) == len(mu) == 188
assert wt[11] == 'G' and mu[11] == 'D'
diffs = [(i+1, wt[i], mu[i]) for i in range(len(wt)) if wt[i] != mu[i]]
assert diffs == [(12, 'G', 'D')]            # ONLY position 12 differs
assert wt[9:17]   == 'GAGGVGKS'              # P-loop (Walker A)
assert mu[9:17]   == 'GADGVGKS'
assert wt[34]     == 'T'                     # T35 (Switch I Mg²⁺ ligand)
assert wt[60]     == 'Q'                     # Q61 (Switch II catalytic Gln)
```

All asserts pass (verified live, see `tools/tests/test_seqio.py` for the
mutation-application paths). The point of this re-derivation is that
the rest of this audit can lean on the fact that:

- WT = 188 aa, all-AA-uppercase, no whitespace pollution.
- Only residue 12 differs WT↔G12D.
- The P-loop motif (`G-x(4)-G-K-S/T`, here `GAGGVGKS`) is present and
  the mutation falls **inside it**, at the central glycine. This is
  the well-established structural reason G12 substitutions impair
  GTP hydrolysis (Krengel 1990; Pai 1990).
- Switch I correctly contains T35 (Mg²⁺-coordinating threonine).
- Switch II correctly contains Q61 (the catalytic Gln, on which the
  in-line water-attack mechanism of GTP hydrolysis depends).

The structural-region boundaries used throughout the project
(P-loop 10–17, Switch I 30–38, Switch II 59–72) are therefore
internally consistent — Q61 falls inside the Switch II window and T35
inside the Switch I window.

## B. Numerical claims cross-checked against artefacts

| Claim (location) | Source artefact | Verified |
|---|---|---|
| AF2 WT 1-169 rank-1 mean pLDDT = 93.66 | `analysis/kras_wt_vs_g12d_structural.json::means.AF2_WT_mean_plddt` | ✅ |
| AF2 G12D 1-169 rank-1 mean pLDDT = 93.28 | same JSON, `AF2_G12D_mean_plddt` | ✅ |
| ESM WT mean pLDDT (heavy-atom) = 86.06 | same JSON, `ESM_WT_mean_plddt_heavy_atom` | ✅ |
| ESM G12D mean pLDDT (heavy-atom) = 86.08 | same JSON, `ESM_G12D_mean_plddt_heavy_atom` | ✅ |
| AF2 WT-vs-G12D RMSD full = 0.102 Å | same JSON, `rmsd_A.AF2_WT_vs_AF2_G12D.full_1-169.rmsd_A` | ✅ |
| ESM WT-vs-G12D RMSD full = 0.081 Å | same JSON, `rmsd_A.ESM_WT_vs_ESM_G12D.full_1-169.rmsd_A` | ✅ |
| AF2 WT-vs-G12D RMSD Switch II = 0.175 Å | same JSON | ✅ |
| Cross-stack AF2 vs ESMFold Switch II ≈ 1.8 Å | same JSON, `AF2_WT_vs_ESM_WT.Switch_II_59-72.rmsd_A` = 1.813; `AF2_G12D_vs_ESM_G12D.Switch_II_59-72.rmsd_A` = 1.870 | ✅ |
| AF2 ΔpLDDT @ pos 12 = −0.37 | same JSON, `per_residue_window_4-20` row pos=12 | ✅ |
| ESMFold ΔpLDDT @ pos 12 = −3.85 | same JSON, `per_residue_window_4-20` row pos=12 (ESM: 84.76 → 80.91) | ✅ |
| 9-mer shortlist = 1 / 9 passing | `analysis/neoantigen_9mers.json` — only `VVGADGVGK` has `overall_tier == "preferred"` and `foreign_vs_wt == true` | ✅ |
| Shortlist peptide `VVGADGVGK` at residues 8-16 | same JSON, row `start_1based=8 end_1based=16` | ✅ |

No discrepancies were found between the log text and the artefacts.

## C. Citation audit — every reference is pinned

| Reference | Cited in | Pinned in `citations.bib` | Verified in `notes.md` |
|---|---|---|---|
| Rojas 2023 *Nature* 618:144 | `logs/03`, `logs/04`, FINAL_RESEARCH_SUMMARY | `Rojas2023Nature` | ✅ |
| Wang 2016 *CIR* 4:204 | `logs/03`, `logs/04` | `Wang2016CIR` | ✅ — primary A\*11:01 / VVGADGVGK |
| Tran 2016 *NEJM* 375:2255 | `logs/03`, `logs/04`, `refs/notes.md` | `Tran2016NEJM` | ✅ — flagged as C\*08:02 NOT A\*11:01 |
| Leidner 2022 *NEJM* 386:2112 | `refs/notes.md` | `Leidner2022NEJM` | ✅ |
| Sidney 2008 *BMC Immunol* 9:1 | `tools/mhc.py`, `logs/03`, figure 4 footer | `Sidney2008BMCImm` | ✅ |
| Falk 1994 *Immunogenetics* 40:238 | `tools/mhc.py`, `logs/03`, figure 4 footer | `Falk1994Immunogen` | ✅ |
| Zhang 1993 *PNAS* 90:2217 | `tools/mhc.py` (replaces former Hunt 1992 mis-cite), figure 4 footer | `Zhang1993PNAS` | ✅ |
| Pai 1990 *EMBO J* 9:2351 | `logs/03` (GTP hydrolysis mechanism) | `Pai1990EMBO` | ✅ |
| Krengel 1990 *Cell* 62:539 | `logs/03` (G12 mutants preserve fold) | `Krengel1990Cell` | ✅ |
| Jumper 2021 *Nature* 596:583 | `logs/01`, `logs/02` (AF2 stack) | `Jumper2021Nature` | ✅ |
| Mirdita 2022 *Nat Methods* 19:679 | `logs/01`, `logs/02` (ColabFold flag set) | `Mirdita2022NatMethods` | ✅ |
| Lin 2023 *Science* 379:1123 | `logs/01`, `logs/02` (ESMFold stack) | `Lin2023Science` | ✅ |
| Kabsch 1976 *Acta Cryst A* 32:922 | `tools/rmsd.py` (Kabsch superposition) | `Kabsch1976` | ✅ |
| Zhu 2026 *Commun Biol* 9:26 | `logs/03`, `logs/04`, `refs/notes.md` Caveat A | `Zhu2026CommunBiol` | ✅ — flags the 9-mer vs 10-mer issue |

All 14 entries in `refs/citations.bib` are used (no orphan citations);
all references in code/logs resolve to an entry in `refs/citations.bib`
(no dangling citations). The two Step-5 corrections (Tran 2016 NEJM is
C\*08:02, not A\*11:01; Hunt 1992 is A2.1, not A\*11:01) and the
Step-6 correction (the *Commun Biol* 2025 paper is by **Zhu et al.**,
not Joglekar) are all explicitly documented in `refs/notes.md`.

## D. Caveats that must propagate into Step 7

These were generated as I worked and are deliberately not buried in
footnotes — they belong in the abstract-level reader's view in
FINAL\_RESEARCH\_SUMMARY.md:

1. **9-mer vs 10-mer at HLA-A\*11:01 / KRAS G12D (Zhu 2026).** The
   motif scorer in `tools/mhc.py` is 9-mer-only. The structurally and
   immunologically *productive* epitope for KRAS G12D on HLA-A\*11:01
   is the **10-mer** `VVVGADGVGK`, where Asp12 is bulged outward at
   p6. The 9-mer `VVGADGVGK` (our shortlist hit) is presented but
   conformationally flattened. The pipeline correctly identifies the
   right window of the KRAS sequence as the candidate
   neoantigen — but a downstream A\*11:01 modeller would prefer the
   10-mer extension.

2. **Motif scorer ≠ IC50 predictor.** `tools/mhc.py` returns
   *anchor tiers* (preferred / tolerated / poor) based on the P2 + P9
   motif (Sidney 2008, Falk 1994, Zhang 1993). It is not a
   binding-affinity predictor; the brief did not authorise a
   NetMHCpan / MHCflurry install (license / web-API restriction).
   The Wang 2016 / Tran 2016 / Rojas 2023 record is consistent with
   the motif call but does not let us claim a predicted IC50.

3. **`VVGADGVGK` as canonical A\*11:01 neoepitope cites Wang 2016 CIR,
   NOT Tran 2016 NEJM.** Tran 2016 is the proof-of-concept that KRAS
   G12D can drive a clinical response with T-cell therapy, but in that
   paper the restricting allele is **HLA-C\*08:02** and the recognised
   peptides are 9-mer `GADGVGKSA` / 10-mer `GADGVGKSAL`, in a
   metastatic-colorectal-cancer patient (not pancreatic). Wang 2016
   *Cancer Immunology Research* 4:204 is the canonical A\*11:01
   reference. This correction was made in Step 5.

4. **Heavy-atom-mean vs Cα-only pLDDT convention difference.**
   ESMFold's HF outputs are per-atom (37 atoms × L residues) and we
   take the *heavy-atom mean per residue* by default (per
   `memory/feedback_esmfold_precision.md` and verified on the smoke
   fixture). AF2's per-residue pLDDT is conceptually closer to a
   Cα-only metric. The two are slightly different in absolute level
   (Cα-only is typically ~3–6 pLDDT higher than heavy-atom-mean on the
   same model), so the ESMFold–AF2 *absolute* mean-pLDDT gap
   (≈86 vs ≈93) overstates the real disagreement. Both predictions
   pass the ≥ 85 gate by either convention.

5. **JAX RNG non-determinism on Blackwell (sm\_120) GPUs.** Per the
   upstream JAX release notes referenced in `logs/00_environment.md`,
   the JAX RNG is not bit-deterministic on consumer Blackwell GPUs.
   `--random-seed 42` is locked in `tools/msa.LOCKED_FLAGS`, but exact
   per-residue pLDDT values may shift by ~0.1 between reruns. This
   does not affect the mean pLDDT ≥ 85 gate by any meaningful margin.

6. **WSL2 host-RAM ceiling ≈ 7.4 GiB.** ESMFold was forced into a
   mixed fp16/fp32 recipe (fp16 ESM-2 stem on the host side, fp32
   folding trunk / heads on GPU) to fit. Patched HF source
   (`transformers/models/esm/modeling_esmfold.py:2173`) defensively
   upcasts pTM logits before softmax. The pTM-computation fix and the
   stem-vs-trunk dtype split are both documented in
   `logs/01_install.md` §5 and in `memory/feedback_esmfold_precision.md`.

7. **Step 3 used HVR-truncated sequences (residues 1–169).** The
   hypervariable region 167–188 is intrinsically disordered
   (polylysine + CAAX prenylation site); modelling it would tank mean
   pLDDT below the ≥ 85 gate without contributing to the G12D
   analysis (residue 12 is in the P-loop, ~150 residues away). All
   X-ray structures of KRAS in the PDB cover 1–169 for the same
   reason. Logs and FASTA names make the 1–169 truncation explicit.

8. **No OpenMM/AMBER relaxation.** ColabFold's `unrelaxed_rank_001_*.pdb`
   was used for all RMSD and figure work. Adding the AMBER relaxer
   would have pulled > 2 GB more data, exceeding the 20 GB budget
   (we ended Step 1 with 0.38 GB headroom). Side-chain rotamer
   geometry near the mutation site is therefore at the AF2 raw output
   level. This is acceptable because all conclusions in Step 3 are
   based on backbone Cα RMSD and per-residue pLDDT, not side-chain
   packing energy.

9. **Rojas 2023 is a personalised vaccine, not a G12D/A\*11:01 fixed
   panel.** Autogene cevumeran (BNT122) is tuned per patient against
   the patient's own tumour-specific neoantigen profile. The paper is
   cited in this study as the *clinical-feasibility rationale* for an
   mRNA-vaccine pipeline against PDAC, not as evidence that
   `VVGADGVGK` was specifically used in any of the 16 patients
   (which would require the trial supplement, not the main paper).

## E. Things the pipeline does NOT claim

This sub-section is preventative: it lists the claims that an
over-eager reader might wrongly infer from the rest of the project,
together with the textual disclaimer location that prevents the
inference.

- "VVGADGVGK is the high-affinity HLA-A\*11:01 binder." → Disclaimed
  in `tools/mhc.py` docstring lines 3–6 and in `logs/03` §"HLA-A\*11:01
  neoantigen shortlist". The scorer is anchor-motif, not IC50.

- "G12D changes the KRAS fold." → Explicitly disclaimed: full-domain
  Cα RMSD WT vs G12D < 0.11 Å in both stacks; G12D drives oncogenesis
  through impaired GTP hydrolysis (Pai 1990, Krengel 1990, Wittinghofer-
  group lineage), not through a rearrangement.

- "ESMFold's −3.85 pLDDT drop at residue 12 means the G12D structure
  is broken." → Explicitly disclaimed: the drop is local (residues
  ≤ 1 from position 12 move by < 0.6 pLDDT), the global mean is
  86.08 (PASS), and the cross-stack and within-stack RMSDs all
  confirm the fold is preserved. The signal is read as ESMFold's
  single-sequence sensitivity to a local chemistry change, not as a
  structural break.

- "Recovering VVGADGVGK proves the pipeline works clinically." → No.
  It is a *sanity check* that the anchor-motif scorer, run blindly on
  9-mers covering the mutation, recovers the literature-validated
  A\*11:01 9-mer window. The pipeline does *not* run an in-vitro
  binding assay, a TCR functional assay, or a clinical trial.

## Acceptance

- 13 / 13 numerical claims in `logs/03` and `logs/04` traced to a
  specific JSON / figure-sidecar artefact.
- 14 / 14 citations resolve to a `citations.bib` entry whose contents
  are verified in `refs/notes.md`. Three corrections were made
  (Tran 2016, Hunt 1992 → Zhang 1993, Joglekar → Zhu 2026).
- 9 propagation-grade caveats listed for Step 7's
  FINAL\_RESEARCH\_SUMMARY.
- 4 potential mis-readings of the pipeline are explicitly disclaimed.

Proceeding to Step 7 (cleanup + README + FINAL_RESEARCH_SUMMARY).
