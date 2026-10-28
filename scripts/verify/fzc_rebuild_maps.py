#!/usr/bin/env python
"""Independent FZ-C determinism check: rebuild the zavg_BEF and GLS maps in a fresh process from committed code and
the committed selection (z0, k), writing ONLY to a scratch directory, then compare to the committed map files.

Reuses scripts/36_build_maps.py:stage_build (that is the code under review) but redirects its output directories;
nothing under results/ is written. Inputs are the cached A-only score blocks (data_cache/part2/blocks).

  PYTHONHASHSEED=0 .venv/bin/python scripts/verify/fzc_rebuild_maps.py --out <scratch dir> [--n-proc 3]
"""
import argparse
import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-proc", type=int, default=3)
    a = ap.parse_args()
    out = Path(a.out).resolve()
    if ROOT / "results" in out.parents or out == ROOT / "results":
        raise SystemExit("refusing to write under results/")
    spec = importlib.util.spec_from_file_location("build_maps", ROOT / "scripts" / "36_build_maps.py")
    bm = importlib.util.module_from_spec(spec)
    sys.modules["build_maps"] = bm                     # needed for spawn-pickling of module-level functions
    spec.loader.exec_module(bm)
    bm.C.ensure_hashseed()
    committed = ROOT / "results" / "part2" / "construction"
    bm.OUT = out
    bm.CACHE = out / "cache"
    t0 = time.perf_counter()
    gctx = bm.gene_context()
    sel = json.loads((committed / "selection" / "selection.json").read_text())
    res = bm.stage_build(gctx, sel, a.n_proc)
    wall = time.perf_counter() - t0
    cmp = {}
    for d in ("zavg_map", "gls_map"):
        man = {l.split()[1]: l.split()[0] for l in (committed / d / "MANIFEST.sha256").read_text().splitlines()}
        cmp[d] = {}
        for f, h in man.items():
            new = out / d / f
            cmp[d][f] = dict(identical=new.exists() and sha(new) == h)
        # parameters.json differs only by code_commit / timings if at all: compare without them
        p_old = json.loads((committed / d / "parameters.json").read_text())
        p_new = json.loads((out / d / "parameters.json").read_text())
        for p in (p_old, p_new):
            p.pop("code_commit", None)
        cmp[d]["parameters.json"]["identical_except_code_commit"] = p_old == p_new
    rep = dict(wall_s=round(wall, 1), n_proc=a.n_proc, results=res, compare=cmp)
    (out / "rebuild_report.json").write_text(json.dumps(rep, indent=1, default=str))
    print(json.dumps(rep, indent=1, default=str))


if __name__ == "__main__":
    main()
