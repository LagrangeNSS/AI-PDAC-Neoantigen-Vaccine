# Step 4 — Publication-grade figures

Window: 2026-05-14T04:25Z → 2026-05-14T04:40Z (UTC).

All figures rendered by `tools/make_step4_figures.py` from
`analysis/kras_wt_vs_g12d_structural.json` + `analysis/neoantigen_9mers.json`.
Each `.png` is paired with a `.txt` provenance sidecar listing inputs and key numbers.

## Figure 1 — `figures/01_plddt_tracks_WT_vs_G12D.png`

Two stacked panels (AF2 top, ESMFold bottom). Each panel:
- X axis: residue 1–169 (KRAS-4B numbering).
- Y axis: pLDDT (0–100).
- WT (blue) and G12D (red) overlaid; near-perfect overlap everywhere except a small dip at residue 12 in ESMFold (visible as the red dip slightly below blue).
- Gold/green/blue shading marks the **P-loop**, **Switch I**, **Switch II** regions.
- Dashed grey line at pLDDT = 85 (hard gate).

**Reading the figure.** AF2 (with the deep KRAS MSA) is highly confident across the entire core (≥90 pLDDT) and only dips into the 50–80 range inside the two switch regions — exactly where solution-state ensembles show real conformational heterogeneity. ESMFold (no MSA) is uniformly less confident (mean ≈ 86) and dips harder in Switch II. **WT and G12D tracks are nearly indistinguishable globally**, supporting the "fold preserved" interpretation in Step 3.

## Figure 2 — `figures/02_plddt_zoom_mutation_site.png`

Two side-by-side bar panels (AF2 left, ESMFold right) showing pLDDT per residue from 4–20 with WT/G12D side-by-side bars. Residue 12 (G12D) is highlighted with a yellow vertical band and a Δ annotation.

- AF2: Δ at residue 12 = **−0.37** (essentially noise; MSA buffers the local sequence change).
- ESMFold: Δ at residue 12 = **−3.85** (clearly localised; from 84.76 to 80.91, dropping below the gate at that single residue).

This is the single most informative figure in the package: it shows that the two stacks see the mutation differently — AF2 sees no problem at position 12 because the KRAS MSA "knows" position 12 tolerates many residues, while ESMFold (single-sequence) reads only the local chemistry and flags the substitution. The fact that the drop is local (not propagated to neighbours) makes it a credible signal, not a numerical artefact.

## Figure 3 — `figures/03_rmsd_by_region.png`

Grouped bar chart, 5 regions × 4 pair-comparisons:
- Within-stack WT-vs-G12D (AF2 blue, ESM green): everywhere **< 0.20 Å** — fold preserved by both stacks.
- Cross-stack AF2-vs-ESMFold (WT orange, G12D red): ~0.6 Å on the structured core, **rising to 1.8 Å at Switch II** — the expected divergence region between an MSA-based and a single-sequence predictor.

Numeric labels printed above each bar. The visual punch: the within-stack bars are vanishingly small (G12D preserves topology) while the Switch II cross-stack bars are an order of magnitude larger (Switch II is genuinely flexible, both predictors place it differently).

## Figure 4 — `figures/04_neoantigen_shortlist.png`

Annotated 9-mer table (9 rows × 10 columns). The G12D 9-mer column shows the peptide with the substituted Asp residue rendered in red+bold; the WT 9-mer column shows the same window with the original Gly greyed out. Tier labels are coloured by tier (preferred=green-bold, tolerated=amber, poor=red).

The single shortlist hit `VVGADGVGK` (residues 8–16) is highlighted with a pale yellow row and a green "PASS" badge. This is the canonical KRAS G12D / HLA-A*11:01 9-mer first reported in:
- **Wang QJ et al. 2016 *Cancer Immunology Research* 4:204–214** — A\*11:01-Tg mouse immunisation; murine TCRs recognising A\*11:01⁺ G12D⁺ pancreatic-tumour lines. This is the primary A\*11:01 / VVGADGVGK reference for this study.
- **Rojas et al. 2023 *Nature*** — autogene cevumeran mRNA vaccine in resected PDAC; vaccines were personalised per tumour, so `VVGADGVGK` is the relevant A\*11:01 window only for A\*11:01-positive G12D-bearing patients (it is not used as a fixed panel across the cohort).
- **Tran et al. 2016 *NEJM*** — TCR-T adoptive transfer in metastatic colorectal cancer. **HLA-C\*08:02-restricted (NOT A\*11:01)**; cited here only as KRAS-G12D-directed-cell-therapy proof of concept, *not* as the A\*11:01 reference.

**Structural caveat (Zhu et al. 2026 *Communications Biology* 9:26):** A\*11:01 presents G12D as both a 9-mer and a 10-mer, but only the **10-mer** VVVGADGVGK adopts a stable, immunogenically productive conformation in the cleft (Asp12 bulged at p6); the 9-mer VVGADGVGK is conformationally flattened. The motif scorer here is 9-mer-only by construction. Recovering the canonical A\*11:01 9-mer through a blind anchor-motif scan is still a non-trivial sanity check on the pipeline.

## Acceptance

- 4 / 4 figures rendered without warnings.
- Each figure has a `.txt` provenance sidecar listing input paths + key numbers + interpretation.
- All figures use absolute pLDDT (0–100 scale) and Å units explicitly.
- No emojis, no decorative styling — straight matplotlib defaults with curated colour choices.

Proceeding to Step 5 (literature verification).
