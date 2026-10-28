"""Import the locked Part 1 package (``src`` of external/part1_locked) and verify the lock.

Every definition shared with Part 1 (AP with sklearn tie semantics, HT weights, Part 1's paired cluster bootstrap,
macro per-cluster AP, metric_suite, cluster construction, UniProt -> DepMap mapping, common-essential flags, the
scorer itself) is used from here, unchanged, so that the test evaluation uses Part 1's definitions by construction.
"""
from __future__ import annotations

import subprocess
import sys

from testeval.paths import LOCKED_COMMIT, PART1

if str(PART1) not in sys.path:
    sys.path.insert(0, str(PART1))

import src as _p1pkg  # noqa: E402

if not str(_p1pkg.__file__).startswith(str(PART1)):
    raise ImportError(f"'src' resolved to {_p1pkg.__file__}, not the locked Part 1 clone {PART1}")

from src import clusters, data, metrics, universe  # noqa: E402,F401
from src.config import SEEDS  # noqa: E402,F401


def lock_state() -> dict:
    """Commit of the clone and any modified tracked files (the untracked `datasets` symlink is expected)."""
    def g(*a):
        return subprocess.run(["git", "-C", str(PART1), *a], capture_output=True, text=True).stdout.strip()
    head = g("rev-parse", "HEAD")
    dirty = [l for l in g("status", "--porcelain").splitlines() if l.strip() and l[3:] != "datasets"]
    return dict(path=str(PART1), commit=head, expected=LOCKED_COMMIT, matches=head == LOCKED_COMMIT,
                modified=dirty, clean=not dirty)


def assert_locked() -> dict:
    st = lock_state()
    if not st["matches"] or not st["clean"]:
        raise RuntimeError(f"Part 1 clone is not the clean locked commit: {st}")
    return st
