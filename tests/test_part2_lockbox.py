"""FZ-LB lockbox draw: iterative stratification on toy labels (blocks kept whole, rare-label targets met, totals,
determinism) and, when the data are present, a fresh-subprocess re-draw (different PYTHONHASHSEED) that must
reproduce the committed lists and manifest hashes byte for byte."""
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import lockbox, paths  # noqa: E402


def _toy(seed=0, n=600):
    rng = np.random.default_rng(seed)
    ids = [f"M{i:04d}" for i in range(n)]
    lin = rng.choice(["a", "b", "c", "d"], n, p=[0.4, 0.3, 0.2, 0.1])
    lab = pd.get_dummies(pd.Series(lin, index=ids), prefix="lineage", prefix_sep=":").astype(int)
    for g, p in (("g1", 0.03), ("g2", 0.05), ("g3", 0.2)):
        lab[g] = (rng.random(n) < p).astype(int)
    # ~10% of models share a patient with the next model (same lineage)
    blocks = pd.Series([f"P{i}" for i in range(n)], index=ids)
    for i in range(0, n - 1, 10):
        blocks.iloc[i + 1] = blocks.iloc[i]
        lab.iloc[i + 1, :4] = lab.iloc[i, :4].values
    return lab, blocks


def test_iterative_stratification_toy():
    lab, blocks = _toy()
    g = ["g1", "g2", "g3"]
    f1 = lockbox.iterative_stratification(lab, blocks, 0.25, 3000, genotype_cols=g)
    f2 = lockbox.iterative_stratification(lab.sample(frac=1, random_state=1), blocks, 0.25, 3000, genotype_cols=g)
    assert f1.sort_index().equals(f2.sort_index())                      # input order does not matter
    assert (f1.groupby(blocks).nunique() == 1).all()                     # patients never straddle
    nL = int((f1 == "L").sum())
    assert abs(nL - 150) <= 3
    t = lockbox.stratum_table(lab, f1, g)
    assert t.dev_L.abs().max() <= 1, t
    f3 = lockbox.iterative_stratification(lab, blocks, 0.25, 3001, genotype_cols=g)
    assert not f1.equals(f3)                                             # the seed matters


def test_target_rule():
    assert lockbox.target_L(16, True) == 4          # n_U < 17: plain 25 %
    assert lockbox.target_L(17, True) == 5          # >= 5 L while leaving >= 12 in A
    assert lockbox.target_L(19, True) == 5
    assert lockbox.target_L(18, False) == 5         # round half up of 4.5 (lineage: plain rule)
    assert lockbox.target_L(14, False) == 4


@pytest.mark.skipif(not paths.LOCKBOX_MANIFEST.exists() or not (paths.CONTEXT_CACHE / "context_table.parquet").exists(),
                    reason="lockbox or context table not built")
def test_lockbox_redraw_is_identical_in_fresh_process(tmp_path):
    env = dict(os.environ, PYTHONHASHSEED="12345")
    r = subprocess.run([sys.executable, str(ROOT / "scripts" / "33_make_lockbox.py"), "--verify"], env=env,
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout + r.stderr
    assert json.loads(r.stdout.strip().splitlines()[-1]) == dict(lockbox_identical=True, A_identical=True,
                                                                 manifest_match=True)
    A = paths.A_MODELS.read_text().split()
    L = paths.LOCKBOX_MODELS.read_text().split()
    assert not set(A) & set(L) and A == sorted(A) and L == sorted(L)
    m = json.loads(paths.LOCKBOX_MANIFEST.read_text())
    assert len(A) + len(L) == m["n_U"] == 1345 and m["n_patients_straddling"] == 0
