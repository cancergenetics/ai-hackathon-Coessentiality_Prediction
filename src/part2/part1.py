"""Read-only imports of the locked Part 1 code (never modified; settings never changed).

* ``coessential.method``  (external/part1_locked/mvp/src): ``Config``, ``CoessentialityModel`` -> zavg_BEF.
* ``coessential.all_pairs``: ``protein_genes`` (the gene table behind Part 1 genes.tsv).
* ``src.baselines``       (external/part1_locked): ``sample_covariance``, ``whitening`` -> GLS (dev variant gls_signed).
"""
from __future__ import annotations

import subprocess
import sys

from part2.paths import LOCKED_COMMIT, PART1, PART1_MVP_SRC

for _p in (str(PART1_MVP_SRC), str(PART1)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import coessential  # noqa: E402
import src as _p1src  # noqa: E402

for _mod, _root in ((coessential, PART1_MVP_SRC), (_p1src, PART1)):
    if not str(_mod.__file__).startswith(str(_root)):
        raise ImportError(f"{_mod.__name__} resolved to {_mod.__file__}, not the locked Part 1 clone {_root}")

from coessential import all_pairs as p1_all_pairs  # noqa: E402,F401
from coessential import data as p1_data  # noqa: E402,F401
from coessential.method import Config, CoessentialityModel  # noqa: E402,F401
from src import baselines as p1_baselines  # noqa: E402,F401

LOCKED_CONFIG = Config()      # the published zavg_BEF configuration (k=20, tau=0.1, q=0.25, d=512, ...)


def lock_state() -> dict:
    def g(*a):
        return subprocess.run(["git", "-C", str(PART1), *a], capture_output=True, text=True).stdout.strip()
    head = g("rev-parse", "HEAD")
    dirty = [l for l in g("status", "--porcelain").splitlines() if l.strip() and l[3:] != "datasets"]
    return dict(commit=head, expected=LOCKED_COMMIT, matches=head == LOCKED_COMMIT, modified=dirty,
                clean=not dirty)


def assert_locked() -> dict:
    st = lock_state()
    if not st["matches"] or not st["clean"]:
        raise RuntimeError(f"Part 1 clone is not the clean locked commit: {st}")
    if Config() != LOCKED_CONFIG:
        raise RuntimeError("locked Config changed")
    return st
