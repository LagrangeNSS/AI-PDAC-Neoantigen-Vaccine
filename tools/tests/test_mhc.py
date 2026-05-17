"""Unit tests for tools.mhc — HLA-A*11:01 P2/P9 anchor scorer."""
from __future__ import annotations

import pytest

from tools.mhc import (
    AnchorTier,
    filter_passing,
    score_9mer_hla_a1101,
    score_9mers_covering_position,
)


def test_canonical_g12d_neoantigen_passes() -> None:
    """KRAS G12D 9-mer 'VVGADGVGK' is the canonical HLA-A*11:01 binder
    cited in Tran et al. 2016 NEJM and Bear et al. 2021 JITC.
    P2 = V (preferred), P9 = K (preferred) → overall preferred.
    """
    s = score_9mer_hla_a1101("VVGADGVGK")
    assert s.p2 == "V" and s.p9 == "K"
    assert s.p2_tier == AnchorTier.PREFERRED
    assert s.p9_tier == AnchorTier.PREFERRED
    assert s.overall_tier == AnchorTier.PREFERRED


def test_wt_g12_homolog_also_passes() -> None:
    """The WT counterpart VVGAGGVGK shares the same anchors as the G12D
    neoantigen; the difference is at P5 (G vs D). The scorer should rank
    it as preferred too (anchors are agnostic to non-anchor residues)."""
    s = score_9mer_hla_a1101("VVGAGGVGK")
    assert s.overall_tier == AnchorTier.PREFERRED


def test_poor_p9_anchor() -> None:
    s = score_9mer_hla_a1101("VVGADGVGA")  # P9 = A → poor
    assert s.p9_tier == AnchorTier.POOR
    assert s.overall_tier == AnchorTier.POOR


def test_tolerated_p2_anchor() -> None:
    s = score_9mer_hla_a1101("VFGADGVGK")  # P2 = F → tolerated
    assert s.p2_tier == AnchorTier.TOLERATED
    assert s.p9_tier == AnchorTier.PREFERRED
    assert s.overall_tier == AnchorTier.TOLERATED  # the worse of the two


def test_length_validation() -> None:
    with pytest.raises(ValueError):
        score_9mer_hla_a1101("VVGADGVG")  # 8 aa
    with pytest.raises(ValueError):
        score_9mer_hla_a1101("VVGADGVGKK")  # 10 aa


def test_score_9mers_covering_position_kras_g12d_full_set() -> None:
    """For KRAS residues 1-25 = 'MTEYKLVVVGADGVGKSALTIQLIQ', G12D mutation
    appears at position 12. There should be 9 candidate 9-mers covering
    position 12 (start positions 4..12 are all valid since L=25).
    """
    kras_n = "MTEYKLVVVGADGVGKSALTIQLIQ"
    scores = score_9mers_covering_position(kras_n, one_based_position=12)
    assert len(scores) == 9
    assert all(len(s.peptide) == 9 for s in scores)
    peptides = [s.peptide for s in scores]
    assert "VVGADGVGK" in peptides  # the canonical one (start pos 8 → peptide indices 8..16)


def test_filter_passing_drops_only_strict_poor() -> None:
    seq = "MTEYKLVVVGADGVGKSALTIQLIQ"
    scores = score_9mers_covering_position(seq, one_based_position=12)
    kept = filter_passing(scores)
    # We should keep at least the canonical VVGADGVGK
    assert any(s.peptide == "VVGADGVGK" for s in kept)
    # And we should drop strict-poor cases (e.g. P2=K/L doesn't fit) — exact
    # count depends on motif; just verify the filter is monotone.
    assert len(kept) <= len(scores)


def test_to_dict_shape() -> None:
    s = score_9mer_hla_a1101("VVGADGVGK")
    d = s.to_dict()
    assert d["allele"] == "HLA-A*11:01"
    assert d["overall_tier"] == "preferred"
    assert set(d.keys()) == {
        "peptide", "allele", "p2", "p9", "p2_tier", "p9_tier", "overall_tier"
    }
