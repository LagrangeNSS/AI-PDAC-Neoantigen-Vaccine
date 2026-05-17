"""MHC-I 9-mer P2/P9 anchor scorer for HLA-A*11:01.

This is a deliberately conservative, motif-based scorer — NOT a binding predictor
(no NetMHCpan / MHCflurry call is performed). The role of this helper is to
shortlist candidate 9-mers around a point mutation by anchor compatibility
before downstream structural/visual review.

Anchor preferences for HLA-A*11:01 (P2 + P9):
- P2  preferred:  V, L, I, M, T
       tolerated: F, A, S
       poor:      every other AA
- P9  preferred:  K, R          (the canonical A3-supertype basic C-terminus)
       tolerated: Y, F
       poor:      every other AA

Sources / cross-checks (cited in refs/citations.bib; verified in refs/notes.md):
- Sidney et al. 2008, BMC Immunology 9:1 — HLA-A*11:01 motif (A3 supertype).
- Falk et al. 1994, Immunogenetics 40:238 — initial A*11:01 peptide motif (K/R at C-term).
- Zhang QJ et al. 1993, PNAS 90:2217 — A*11:01 nonamer motif from viral/cellular peptides.
  (Note: a previous draft erroneously cited Hunt et al. 1992 Science 255:1261 here;
  that paper is the HLA-A2.1 motif paper, not A*11:01. Corrected in Step 5.)

The scorer returns a tier (preferred / tolerated / poor) per anchor and an
overall combined tier (the worse of the two).
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum
from typing import Iterable

# Per-anchor preference dictionaries for HLA-A*11:01.
A1101_P2_PREFERRED = frozenset("VLIMT")
A1101_P2_TOLERATED = frozenset("FAS")
A1101_P9_PREFERRED = frozenset("KR")
A1101_P9_TOLERATED = frozenset("YF")


class AnchorTier(IntEnum):
    """Ordered tier; bigger is better."""
    POOR = 0
    TOLERATED = 1
    PREFERRED = 2

    @property
    def label(self) -> str:
        return self.name.lower()


@dataclass(frozen=True)
class AnchorScore:
    """Per-9mer anchor result for a single HLA-I allele."""
    peptide: str
    allele: str
    p2: str
    p9: str
    p2_tier: AnchorTier
    p9_tier: AnchorTier

    @property
    def overall_tier(self) -> AnchorTier:
        return AnchorTier(min(int(self.p2_tier), int(self.p9_tier)))

    def to_dict(self) -> dict:
        return {
            "peptide": self.peptide,
            "allele": self.allele,
            "p2": self.p2,
            "p9": self.p9,
            "p2_tier": self.p2_tier.label,
            "p9_tier": self.p9_tier.label,
            "overall_tier": self.overall_tier.label,
        }


def _classify(residue: str, preferred: frozenset, tolerated: frozenset) -> AnchorTier:
    if residue in preferred:
        return AnchorTier.PREFERRED
    if residue in tolerated:
        return AnchorTier.TOLERATED
    return AnchorTier.POOR


def score_9mer_hla_a1101(peptide: str) -> AnchorScore:
    """Score a single 9-mer peptide against the HLA-A*11:01 P2/P9 motif."""
    if len(peptide) != 9:
        raise ValueError(f"peptide must be 9 aa, got len={len(peptide)}: {peptide!r}")
    peptide = peptide.upper()
    p2 = peptide[1]
    p9 = peptide[8]
    return AnchorScore(
        peptide=peptide,
        allele="HLA-A*11:01",
        p2=p2,
        p9=p9,
        p2_tier=_classify(p2, A1101_P2_PREFERRED, A1101_P2_TOLERATED),
        p9_tier=_classify(p9, A1101_P9_PREFERRED, A1101_P9_TOLERATED),
    )


def score_9mers_covering_position(
    sequence: str,
    one_based_position: int,
) -> list[AnchorScore]:
    """Return AnchorScore for every 9-mer that covers `one_based_position`.

    There are at most 9 such 9-mers (positions [pos-8..pos] in 1-based start
    coordinates), clipped by sequence bounds. The mutation must lie within
    the 9-mer for it to be considered a candidate neoantigen.
    """
    if not (1 <= one_based_position <= len(sequence)):
        raise ValueError(f"position {one_based_position} out of bounds for L={len(sequence)}")
    lo = max(1, one_based_position - 8)
    hi = min(len(sequence) - 8, one_based_position)  # last valid 9-mer start
    out: list[AnchorScore] = []
    for start in range(lo, hi + 1):
        pep = sequence[start - 1: start + 8]
        if len(pep) == 9:
            out.append(score_9mer_hla_a1101(pep))
    return out


def filter_passing(
    scores: Iterable[AnchorScore],
    min_tier: AnchorTier = AnchorTier.TOLERATED,
) -> list[AnchorScore]:
    """Keep only AnchorScore entries with `overall_tier >= min_tier`.

    Default = TOLERATED, i.e. drop strict POOR-on-either-anchor entries.
    """
    return [s for s in scores if s.overall_tier >= min_tier]
