# refs/notes.md — Step 5 literature verification log

Window: 2026-05-14T04:40Z → 2026-05-14T05:10Z (UTC).

Every claim in `refs/citations.bib` was cross-checked against the
publisher / PubMed / PMC record via WebSearch before the BibTeX entry was
finalised. This file records which claims were verified, which had to be
corrected, and the residual caveats that propagate into Step 6.

## Method

For each candidate reference cited in code or logs:

1. WebSearch for `<first author> <year> <journal> <key claim>`.
2. Pull the publisher landing page (or the PMC mirror) to confirm
   journal, volume, issue, page range, DOI, PMID, and the key claim that
   the citation is being used to support.
3. If the canonical record contradicts the way I had cited it, write
   down the discrepancy in the corresponding `note` field of
   `citations.bib` and propagate the correction back into the code/log
   that originally cited it.

No reference is cited from training-data memory alone; every entry has a
direct link in its `url` field.

## What was verified cleanly (no correction needed)

| Citation key | Used for | Verified? |
|---|---|---|
| `Rojas2023Nature` | Personalised mRNA-vaccine context (autogene cevumeran in resected PDAC) | ✅ Nature 618:144-150, DOI 10.1038/s41586-023-06063-y |
| `Leidner2022NEJM` | TCR-T case report in metastatic PDAC (HLA-C\*08:02 G12D) | ✅ NEJM 386:2112-2119, DOI 10.1056/NEJMoa2119662 |
| `Sidney2008BMCImm` | HLA-A\*11:01 ↦ A3 supertype, basic C-term anchor | ✅ BMC Immunol 9:1, DOI 10.1186/1471-2172-9-1 |
| `Falk1994Immunogen` | Primary A\*11:01 pool-sequencing motif (P2 hydrophobic, P9=K) | ✅ Immunogenetics 40:238-241, DOI 10.1007/BF00167086 |
| `Pai1990EMBO` | H-ras p21 GTP-bound 1.35 Å structure, in-line water hydrolysis | ✅ EMBO J 9:2351-2359, DOI 10.1002/j.1460-2075.1990.tb07409.x |
| `Krengel1990Cell` | H-ras p21 mutant X-ray structures show preserved fold | ✅ Cell 62:539-548, DOI 10.1016/0092-8674(90)90018-A |
| `Jumper2021Nature` | AlphaFold-2 method paper | ✅ Nature 596:583-589, DOI 10.1038/s41586-021-03819-2 |
| `Mirdita2022NatMethods` | ColabFold method paper + flag set | ✅ Nat Methods 19:679-682, DOI 10.1038/s41592-022-01488-1 |
| `Lin2023Science` | ESMFold method paper | ✅ Science 379:1123-1130, DOI 10.1126/science.ade2574 |
| `Kabsch1976` | Kabsch superposition (Biopython SVDSuperimposer) | ✅ Acta Cryst A 32:922-923, DOI 10.1107/S0567739476001873 |

## What had to be corrected

### Correction 1 — Tran 2016 NEJM is **NOT** the A\*11:01 reference

**Original drafting error.** Early drafts of `logs/03_core_workflow.md`
and `logs/04_figures.md` named **Tran et al. *NEJM* 2016** as the source
of the canonical KRAS G12D / **HLA-A\*11:01** / `VVGADGVGK` triad. That
attribution is wrong.

What Tran 2016 actually reports (verified at
[nejm.org/doi/full/10.1056/NEJMoa1609279](https://www.nejm.org/doi/full/10.1056/NEJMoa1609279)):

- Patient with metastatic **colorectal** cancer (NOT pancreatic).
- TIL transfer of **HLA-C\*08:02-restricted** (NOT A\*11:01) TCR clonotypes.
- Recognised peptides were the G12D 9-mer **GADGVGKSA** and 10-mer
  **GADGVGKSAL** (NOT VVGADGVGK).
- Tumour escape at 9 months by chr-6 haplotype loss eliminating
  HLA-C\*08:02.

**Correct primary A\*11:01 reference:** Wang QJ, Yu Z, Griffith K,
Hanada K, Restifo NP, Yang JC. *Identification of T-cell Receptors
Targeting KRAS-Mutated Human Tumors.* **Cancer Immunology Research**
**4**(3):204–214 (2016). DOI: 10.1158/2326-6066.CIR-15-0188 — entered
as `Wang2016CIR` in `citations.bib`.

Wang 2016 used HLA-A\*11:01-transgenic mice immunised with the G12V and
G12D 10-mer peptides, isolated murine TCRs (TRAV4-4\*01 / TRBV12-2\*01
for G12D), and showed recognition of HLA-A\*11:01⁺ G12D⁺ pancreatic
tumour lines. This is the canonical primary reference for the
VVGADGVGK / HLA-A\*11:01 axis used in this study.

**Files fixed:**
- `logs/03_core_workflow.md` — shortlist paragraph now cites Wang 2016
  as the A\*11:01 primary and explicitly notes that Tran 2016 is the
  C\*08:02 case in colorectal cancer.
- `logs/04_figures.md` — Figure 4 caption rewritten accordingly.
- `refs/citations.bib` — `note` fields on both `Tran2016NEJM` and
  `Wang2016CIR` spell out the correction.

### Correction 2 — Hunt 1992 Science is **NOT** the A\*11:01 motif paper

**Original drafting error.** Early drafts of `tools/mhc.py` cited
**Hunt et al. 1992 Science 255:1261** as a cross-check on the HLA-A\*11:01
P2/P9 motif.

What Hunt 1992 actually reports: pool sequencing of self-peptides
eluted from **HLA-A2.1** (= HLA-A\*02:01), NOT HLA-A\*11:01.

**Correct replacement:** Zhang QJ, Gavioli R, Klein G, Masucci MG.
*An HLA-A11-specific motif in nonamer peptides derived from viral and
cellular proteins.* **PNAS** **90**(6):2217–2221 (1993). PMID 7681587.
URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC46057/ — entered as
`Zhang1993PNAS` in `citations.bib`. This paper confirms the A\*11:01
nonamer motif (hydrophobic P2, basic P9 = K) using viral and cellular
peptides, which is exactly what `tools/mhc.py` needs to back its
P2/P9 anchor tiers.

**Files fixed:**
- `tools/mhc.py` — docstring citation block now lists Zhang 1993 with
  an explicit comment noting the prior mis-cite.
- `tools/make_step4_figures.py` — figure-4 footer string now reads
  `(Sidney 2008; Falk 1994; Zhang 1993)`.
- `logs/03_core_workflow.md` — shortlist preamble updated.
- `refs/citations.bib` — `note` on `Zhang1993PNAS` calls out the
  Hunt-1992 / A2.1 mix-up explicitly.

Figure 4 was re-rendered after the fix; the PNG and `.txt` sidecar both
carry the corrected citation block.

### Correction 3 — Authorship of the A\*11:01 / G12D 10-mer structural paper

**Original drafting error.** I initially attributed the
*Communications Biology* 2025 paper "Structure guided analysis of KRAS
G12 mutants in HLA-A\*11:01 reveals a length encoded immunogenic
advantage in G12D" (DOI 10.1038/s42003-025-09285-0) to **Alok V.
Joglekar** as first author.

What the canonical record actually shows (verified via WebSearch on the
DOI; see also the [ResearchGate listing](https://www.researchgate.net/publication/398294614)):

- First author: **Jiali Zhu**.
- Co-authors: Zhifeng Chen, Xi Xu, Yuxuan Wang, Pei Liu, Maorong Wen,
  Quanmeng Wang, Yuchen He, Hongwei Jin, Hongjuan Xue, Shuqing Wang,
  Ke Xu, Linlin Zhao.
- Volume 9, article 26; the journal lists the formal publication
  year as **2026** (received 21 May 2025, accepted 17 Nov 2025).

The paper's *content* — that the A\*11:01-presented 10-mer
VVVGADGVGK is the productive TCR-engaging epitope while the 9-mer is
conformationally flattened — was reported correctly in early drafts.
Only the author attribution was wrong.

**Files fixed:**
- `refs/citations.bib` — entry renamed `Joglekar2025CommunBiol` →
  `Zhu2026CommunBiol`, author list, year, volume, and pages corrected;
  the correction is documented inside the `note` field.
- `logs/03_core_workflow.md` — Step 3 shortlist paragraph now reads
  "Zhu et al. 2026 Communications Biology 9:26".
- `logs/04_figures.md` — Figure 4 caption updated.

## Residual caveats

### Caveat A — 9-mer vs 10-mer for HLA-A\*11:01 / G12D

**Zhu et al. 2026 *Communications Biology* (vol 9, p 26)** (
[s42003-025-09285-0](https://www.nature.com/articles/s42003-025-09285-0))
performed a structure-guided analysis of KRAS-G12 mutants in
HLA-A\*11:01 and found:

- The **9-mer** `VVGADGVGK` is conformationally **flattened** in the
  A\*11:01 cleft (Asp12 lies along the floor) and is less
  immunogenically productive than the 10-mer in their assays.
- The **10-mer** `VVVGADGVGK` is bulged at p6 (with Asp12 protruding
  toward the TCR) and is the productive TCR-engaging epitope.

`tools/mhc.py` and `analysis/neoantigen_9mers.json` are 9-mer-only by
construction (the scorer enumerates 9-mers covering the mutation; the
P2/P9 motif is a 9-mer-anchor scheme). Both the 9-mer and 10-mer share
the P2=V / P9=K anchors that the scorer evaluates, so a 10-mer scan
would also flag the same window. But the rank order between the 9-mer
and the 10-mer would only be settled by an explicit structural
follow-up at A\*11:01, which is outside the brief of this feasibility
study.

This caveat is referenced in `refs/citations.bib::Zhu2026CommunBiol`
and is propagated into:
- `logs/03_core_workflow.md` (shortlist paragraph)
- `logs/04_figures.md` (Figure 4 caption)
- Will appear again in `logs/06_self_audit.md` (Step 6).

### Caveat B — `tools/mhc.py` is a motif scorer, not a binding predictor

The motif scorer returns P2/P9 anchor *tiers* (preferred / tolerated /
poor), not predicted IC50 in nM. It was chosen because the brief
restricts the install to locally resident GPU tools (no NetMHCpan
license, no MHCflurry web call). Wang 2016 reports in-vitro measured
A\*11:01 / `VVGADGVGK` binding in the strong-binder range, which is
consistent with the motif call here but does not make this an IC50
prediction. This is restated in `logs/03_core_workflow.md` ("What's
NOT in Step 3") and will be flagged again in Step 6.

### Caveat C — Rojas 2023 cohort does not target G12D / A\*11:01 specifically

Autogene cevumeran (BNT122) is a **personalised** RNA-neoantigen vaccine:
each patient's vaccine is tailored to their own tumour's neoantigen
landscape. The Nature 2023 paper does not single out KRAS G12D /
HLA-A\*11:01 as a fixed panel across the cohort. `refs/citations.bib`
records this fact in the `Rojas2023Nature` note; the logs cite Rojas
2023 only as the clinical-feasibility rationale for an mRNA-vaccine
pipeline, not as a per-allele neoantigen reference.

## Acceptance

- 14 / 14 BibTeX entries verified against their canonical record.
- 2 mis-citations identified and corrected (Tran 2016 NEJM; Hunt 1992
  Science); all downstream files updated.
- 1 important structural caveat carried forward into Steps 6/7
  (9-mer vs 10-mer at A\*11:01, Zhu 2026).
- 2 scope caveats restated (motif scorer ≠ IC50; Rojas 2023 = personalised).

Proceeding to Step 6 (anti-hallucination self-audit).
