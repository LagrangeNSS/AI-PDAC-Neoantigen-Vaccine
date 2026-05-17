"""Unit tests for tools.seqio."""
from __future__ import annotations

from pathlib import Path

import pytest

from tools.seqio import (
    FastaRecord,
    apply_point_mutation,
    iter_kmers,
    read_fasta,
    slice_window,
    validate_protein,
    write_fasta,
)


def test_read_fasta_single_record(tmp_path: Path) -> None:
    p = tmp_path / "x.fa"
    p.write_text(">my_id some description\nMQIF\nVKTL\n")
    [rec] = read_fasta(p)
    assert rec.id == "my_id"
    assert rec.description == "some description"
    assert rec.sequence == "MQIFVKTL"
    assert rec.length == 8


def test_read_fasta_multi_records_skips_blank_lines(tmp_path: Path) -> None:
    p = tmp_path / "x.fa"
    p.write_text(">a\nMQI\n\n>b extra info\nFV\n  KTL\n\n")
    a, b = read_fasta(p)
    assert (a.id, a.sequence) == ("a", "MQI")
    assert (b.id, b.sequence, b.description) == ("b", "FVKTL", "extra info")


def test_read_fasta_strips_whitespace_inside_sequence(tmp_path: Path) -> None:
    p = tmp_path / "x.fa"
    p.write_text(">a\nMQ I\tF\n")
    [rec] = read_fasta(p)
    assert rec.sequence == "MQIF"


def test_write_then_read_roundtrip(tmp_path: Path) -> None:
    p = tmp_path / "out.fa"
    rec = FastaRecord(id="kras_g12d", description="KRAS 1-189 G12D", sequence="MTEYKLVVVGADGVG")
    write_fasta([rec], p)
    [got] = read_fasta(p)
    assert got == rec


def test_validate_protein_rejects_non_canonical() -> None:
    with pytest.raises(ValueError):
        validate_protein("MQIFVKTLXXXZ")


def test_validate_protein_accepts_canonical() -> None:
    validate_protein("ACDEFGHIKLMNPQRSTVWY")  # should not raise


def test_apply_point_mutation_g12d() -> None:
    # KRAS positions 1-15 (WT): MTEYKLVVVGAGGVG
    wt = "MTEYKLVVVGAGGVG"
    mut = apply_point_mutation(wt, "G12D")
    assert mut == "MTEYKLVVVGADGVG"
    assert mut[11] == "D"


def test_apply_point_mutation_wrong_wt_rejected() -> None:
    with pytest.raises(ValueError, match="WT mismatch"):
        apply_point_mutation("MTEYKLVVVGAGGVG", "A12D")


def test_apply_point_mutation_position_out_of_bounds() -> None:
    with pytest.raises(ValueError, match="out of bounds"):
        apply_point_mutation("MQIF", "M99D")


def test_apply_point_mutation_malformed() -> None:
    with pytest.raises(ValueError):
        apply_point_mutation("MQIF", "Q?")


def test_slice_window_clips_at_left_edge() -> None:
    seq, lo, hi = slice_window("MQIFVKTL", center_pos=2, flank=4)
    assert seq == "MQIFVK"
    assert (lo, hi) == (1, 6)


def test_slice_window_clips_at_right_edge() -> None:
    # MQIFVKTL (1-indexed: M1 Q2 I3 F4 V5 K6 T7 L8)
    # center=7, flank=4 → start=max(1,3)=3, end=min(8,11)=8 → "IFVKTL"
    seq, lo, hi = slice_window("MQIFVKTL", center_pos=7, flank=4)
    assert seq == "IFVKTL"
    assert (lo, hi) == (3, 8)


def test_iter_kmers_basic() -> None:
    out = list(iter_kmers("ABCDE", 3))
    assert out == [(1, "ABC"), (2, "BCD"), (3, "CDE")]
