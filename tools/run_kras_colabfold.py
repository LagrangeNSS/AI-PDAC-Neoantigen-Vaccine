"""Driver: ColabFold AF2 on KRAS WT/G12D G-domain (residues 1-169).

Uses tools.msa.run_colabfold (brief-locked flags + 60/300/900 retry +
2 h SIGALRM + ≥85 pLDDT hard gate).
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path

from tools.msa import run_colabfold, ColabFoldGateError


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--log-path", required=True, type=Path)
    parser.add_argument("--plddt-gate", type=float, default=85.0)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s %(message)s")

    t0 = time.time()
    try:
        run = run_colabfold(
            args.fasta,
            args.out_dir,
            plddt_gate=args.plddt_gate,
            log_path=args.log_path,
        )
    except ColabFoldGateError as e:
        print(f"GATE FAIL: {e}", file=sys.stderr)
        return 1
    elapsed = time.time() - t0
    print(json.dumps({
        "fasta": str(args.fasta),
        "out_dir": str(args.out_dir),
        "rank1_scores_path": str(run.rank1_scores_path),
        "rank1_mean_plddt": run.rank1_mean_plddt,
        "wall_seconds": round(elapsed, 1),
        "gate_pass": True,
    }, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
