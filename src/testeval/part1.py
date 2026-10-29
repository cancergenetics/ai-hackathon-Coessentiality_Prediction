"""Import the locked Part 1 package (``src`` of external/part1_locked) and verify the lock.

Every definition shared with Part 1 (AP with sklearn tie semantics, HT weights, Part 1's paired cluster bootstrap,
macro per-cluster AP, metric_suite, cluster construction, UniProt -> DepMap mapping, common-essential flags, the
scorer itself) is used from here, unchanged, so that the test evaluation uses Part 1's definitions by construction.
"""
from __future__ import annotations

import subprocess
import sys

from testeval.paths import LOCKED_COMMIT, LOCKED_OBJECTS, PART1, PUBLIC_MIRROR_COMMIT

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
    objects_match = all(g("rev-parse", f"HEAD:{p}") == oid for p, oid in LOCKED_OBJECTS.items())
    via = "locked commit" if head == LOCKED_COMMIT else \
        "public mirror, content-identical" if head == PUBLIC_MIRROR_COMMIT and objects_match else None
    return dict(path=str(PART1), commit=head, expected=LOCKED_COMMIT, matches=via is not None, via=via,
                modified=dirty, clean=not dirty)


def assert_locked() -> dict:
    st = lock_state()
    if not st["matches"] or not st["clean"]:
        raise RuntimeError(f"Part 1 clone is not the clean locked commit: {st}")
    return st


_FINAL_EVAL = None


def final_eval():
    """Part 1's scripts/final_eval.py as a module (for its locked verdict function ``decide`` and constants)."""
    global _FINAL_EVAL
    if _FINAL_EVAL is None:
        import importlib.util
        spec = importlib.util.spec_from_file_location("part1_final_eval", PART1 / "scripts" / "final_eval.py")
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        _FINAL_EVAL = mod
    return _FINAL_EVAL
