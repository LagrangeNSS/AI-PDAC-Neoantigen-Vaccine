"""Pytest plumbing: ensure WORKDIR is on sys.path so `from tools.X import Y` resolves."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path("/mnt/d/AssignmentOPUS")
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
