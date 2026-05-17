"""Minimal FASTA I/O + point-mutation helpers for the KRAS-G12D study.

We deliberately avoid Biopython here to keep the helper light and dependency-cheap
for the parts of the pipeline that only need to read/write FASTA.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

AA_ONE = set("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True)
class FastaRecord:
    """A single FASTA record: identifier, optional description, and sequence."""
    id: str
    description: str  # the part of the header line after the first whitespace
    sequence: str

    @property
    def length(self) -> int:
        return len(self.sequence)

    def to_fasta(self, line_width: int = 80) -> str:
        header = f">{self.id}" + (f" {self.description}" if self.description else "")
        chunks = [self.sequence[i:i + line_width] for i in range(0, len(self.sequence), line_width)]
        return header + "\n" + "\n".join(chunks) + "\n"


def read_fasta(path: str | Path) -> list[FastaRecord]:
    """Parse a FASTA file into FastaRecord objects.

    - Empty lines are skipped.
    - The first whitespace-separated token after '>' is the id.
    - The remainder of the header is the description.
    - Sequences are upper-cased; whitespace within sequence lines is removed.
    - We do NOT validate the alphabet here (let downstream helpers decide).
    """
    path = Path(path)
    records: list[FastaRecord] = []
    rid: str | None = None
    desc = ""
    seq_chunks: list[str] = []
    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n").rstrip("\r")
            if not line.strip():
                continue
            if line.startswith(">"):
                if rid is not None:
                    records.append(FastaRecord(rid, desc, "".join(seq_chunks).upper()))
                parts = line[1:].split(None, 1)
                rid = parts[0] if parts else ""
                desc = parts[1] if len(parts) > 1 else ""
                seq_chunks = []
            else:
                seq_chunks.append(line.replace(" ", "").replace("\t", ""))
        if rid is not None:
            records.append(FastaRecord(rid, desc, "".join(seq_chunks).upper()))
    return records


def write_fasta(records: Iterable[FastaRecord], path: str | Path, line_width: int = 80) -> None:
    """Write one or more FastaRecords to a FASTA file (UNIX line endings)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="\n") as fh:
        for rec in records:
            fh.write(rec.to_fasta(line_width=line_width))


def validate_protein(sequence: str) -> None:
    """Raise ValueError if `sequence` contains a character outside the 20 canonical AAs."""
    bad = sorted({c for c in sequence if c not in AA_ONE})
    if bad:
        raise ValueError(f"Non-canonical amino acid(s) found: {bad}")


def apply_point_mutation(sequence: str, mutation: str) -> str:
    """Apply a single-residue point mutation in HGVS-like 1-based notation.

    Args:
        sequence: WT protein sequence (e.g. KRAS 1–189).
        mutation: 'X<pos><Y>' where X is the WT AA, pos is 1-based, Y is the mutant AA.
                  Example: 'G12D' for KRAS G12D.

    Returns:
        The mutated sequence.

    Raises:
        ValueError: if pos is out of bounds, X doesn't match the WT residue, or
                    the format is malformed.
    """
    if len(mutation) < 3:
        raise ValueError(f"mutation string too short: {mutation!r}")
    wt_aa = mutation[0]
    mut_aa = mutation[-1]
    try:
        pos = int(mutation[1:-1])
    except ValueError as e:
        raise ValueError(f"could not parse position from {mutation!r}") from e
    if wt_aa not in AA_ONE or mut_aa not in AA_ONE:
        raise ValueError(f"both flanks of {mutation!r} must be canonical AAs")
    if not (1 <= pos <= len(sequence)):
        raise ValueError(
            f"position {pos} out of bounds for sequence of length {len(sequence)}"
        )
    actual = sequence[pos - 1]
    if actual != wt_aa:
        raise ValueError(
            f"WT mismatch: mutation {mutation} expected {wt_aa} at pos {pos}, "
            f"but sequence has {actual}"
        )
    return sequence[: pos - 1] + mut_aa + sequence[pos:]


def slice_window(sequence: str, center_pos: int, flank: int) -> tuple[str, int, int]:
    """Return a window of `2*flank+1` residues centered on `center_pos` (1-based),
    clipped to sequence bounds. Also returns the (start, end) 1-based inclusive
    coordinates of the window in the original sequence.

    Useful for extracting neoantigen peptide context around a point mutation.
    """
    if not (1 <= center_pos <= len(sequence)):
        raise ValueError(f"center_pos {center_pos} out of bounds")
    start = max(1, center_pos - flank)
    end = min(len(sequence), center_pos + flank)
    return sequence[start - 1: end], start, end


def iter_kmers(sequence: str, k: int) -> Iterator[tuple[int, str]]:
    """Yield (1-based start position, k-mer) for every k-mer in `sequence`."""
    if k <= 0:
        raise ValueError("k must be positive")
    for i in range(len(sequence) - k + 1):
        yield i + 1, sequence[i: i + k]
