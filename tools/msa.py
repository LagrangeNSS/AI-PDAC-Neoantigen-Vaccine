"""Thin wrapper for ColabFold's MSA pipeline + AF2 inference (`colabfold_batch`),
with brief-locked flags, exponential-backoff retry, and a per-prediction timeout.

Why a wrapper? Step 3 of the pipeline launches many predictions (WT, several
G12D peptide contexts, etc.); we want one canonical place that:
  - pins the locked flag set,
  - applies the 60/300/900 s exponential backoff schedule,
  - enforces a 2 h per-prediction timeout (SIGALRM via tools.retry),
  - and verifies the AF2 hard gate (rank-1 mean pLDDT >= 85) immediately
    after the run completes, raising ColabFoldGateError otherwise.
"""
from __future__ import annotations

import logging
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from tools.plddt import PlddtTrack, find_colabfold_rank1, from_colabfold_scores
from tools.retry import DEFAULT_BACKOFF_SECONDS, DEFAULT_PER_CALL_TIMEOUT_SECONDS, retry

log = logging.getLogger("tools.msa")

# Brief-locked, do not modify without explicit user approval.
LOCKED_FLAGS: tuple[str, ...] = (
    "--num-models", "5",
    "--num-recycle", "3",
    "--rank", "plddt",
    "--model-type", "alphafold2_ptm",
    "--random-seed", "42",
    "--msa-mode", "mmseqs2_uniref_env",
)

# Hardware-tuning flag, allowed (does NOT alter the model). Reason: WSL2 host RAM
# is 7.4 GiB; JAX's unified-memory allocator probes 50 GB → OOM-killer fires.
HARDWARE_FLAGS: tuple[str, ...] = ("--disable-unified-memory",)


class ColabFoldGateError(RuntimeError):
    """Raised when a ColabFold run completes but rank-1 mean pLDDT < gate."""


@dataclass(frozen=True)
class ColabFoldRun:
    """Summary of a single ColabFold prediction directory."""
    out_dir: Path
    rank1_scores_path: Path
    rank1_plddt: PlddtTrack

    @property
    def rank1_mean_plddt(self) -> float:
        return self.rank1_plddt.mean


def _build_cmd(fasta: Path, out_dir: Path, extra_flags: tuple[str, ...]) -> list[str]:
    binary = shutil.which("colabfold_batch")
    if binary is None:
        raise RuntimeError("colabfold_batch not found on PATH (did you source ./env/activate_project.sh?)")
    cmd = [binary, *LOCKED_FLAGS, *HARDWARE_FLAGS, *extra_flags, str(fasta), str(out_dir)]
    return cmd


def _run_colabfold_once(fasta: Path, out_dir: Path, extra_flags: tuple[str, ...],
                        log_path: Path | None) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = _build_cmd(fasta, out_dir, extra_flags)
    log.info("colabfold cmd: %s", " ".join(cmd))
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(log_path, "ab") as fh:
            fh.write(("$ " + " ".join(cmd) + "\n").encode())
            proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, check=False)
    else:
        proc = subprocess.run(cmd, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            f"colabfold_batch exited with code {proc.returncode}. See {log_path}"
        )


def run_colabfold(
    fasta: str | Path,
    out_dir: str | Path,
    *,
    plddt_gate: float = 85.0,
    extra_flags: tuple[str, ...] = (),
    log_path: str | Path | None = None,
    per_call_timeout_seconds: int = DEFAULT_PER_CALL_TIMEOUT_SECONDS,
    backoff_seconds: tuple[int, ...] = DEFAULT_BACKOFF_SECONDS,
    skip_if_done: bool = True,
) -> ColabFoldRun:
    """Run ColabFold + AF2 with brief-locked flags, retry, timeout, and gate check.

    Args:
        fasta: input FASTA path.
        out_dir: ColabFold output directory.
        plddt_gate: minimum rank-1 mean pLDDT; raise ColabFoldGateError below.
        extra_flags: additional flags appended after LOCKED_FLAGS + HARDWARE_FLAGS.
                     Do NOT pass flags that contradict LOCKED_FLAGS.
        log_path: append-only log for the subprocess stdout/stderr.
        per_call_timeout_seconds: SIGALRM timeout per attempt (default 2 h).
        backoff_seconds: retry schedule (default 60/300/900).
        skip_if_done: if True and `<out_dir>/*_scores_rank_001*.json` already
                      exists, do not re-run; just load and gate-check.

    Returns:
        ColabFoldRun with rank-1 pLDDT loaded.
    """
    fasta = Path(fasta)
    out_dir = Path(out_dir)
    log_path = Path(log_path) if log_path is not None else None

    if skip_if_done and out_dir.exists():
        existing = list(out_dir.glob("*_scores_rank_001*.json"))
        if existing:
            log.info("skip_if_done: rank-1 already present at %s", existing[0].name)
            track = from_colabfold_scores(existing[0])
            run = ColabFoldRun(out_dir=out_dir, rank1_scores_path=existing[0], rank1_plddt=track)
            _gate_or_raise(run, plddt_gate)
            return run

    def _fn() -> None:
        _run_colabfold_once(fasta, out_dir, extra_flags, log_path)

    retry(
        _fn,
        backoff_seconds=backoff_seconds,
        per_call_timeout_seconds=per_call_timeout_seconds,
        label=f"colabfold[{fasta.name}]",
    )

    rank1_path = find_colabfold_rank1(out_dir)
    track = from_colabfold_scores(rank1_path)
    run = ColabFoldRun(out_dir=out_dir, rank1_scores_path=rank1_path, rank1_plddt=track)
    _gate_or_raise(run, plddt_gate)
    return run


def _gate_or_raise(run: ColabFoldRun, threshold: float) -> None:
    mean = run.rank1_mean_plddt
    if mean < threshold:
        raise ColabFoldGateError(
            f"ColabFold rank-1 mean pLDDT {mean:.2f} < gate {threshold} "
            f"(out_dir={run.out_dir})"
        )
    log.info("colabfold gate: rank-1 mean pLDDT %.2f >= %.1f (PASS)", mean, threshold)
