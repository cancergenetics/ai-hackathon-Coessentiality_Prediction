#!/usr/bin/env python
"""Full dress rehearsal of the lockbox driver (review MUST-2c): the whole of scripts/51_lockbox_confirm.py, BENCH and
BRIDGE-L included, on the REAL A data in a THROW-AWAY `git clone` (not a worktree), with a dummy outcome loader that
returns random values for every non-A request. The dummy sees the lockbox model IDs and the gene keys only; it never
reads a lockbox value. A requests (A effects, DEMETER2 on A) go to the project loader.

Guards (refuses otherwise): the repository root is a plain clone (`.git` is a directory, an `origin` remote exists) and
carries the untracked marker file `.lockbox_rehearsal`; the marker is never created by any project script.

Setup (from a scratch directory; data_cache/part2 is a real directory so persisted reads and the GLS discovery tables
stay inside the clone):
    git clone <worktree> <clone>; touch <clone>/.lockbox_rehearsal
    ln -s <main>/datasets <clone>/datasets; ln -s <main>/external <clone>/external
    mkdir -p <clone>/data_cache/part2; symlink every child of <main>/data_cache/part2 except assoc, copy assoc
    cd <clone>; PART1_DIR=... PYTHONHASHSEED=0 caffeinate -i <venv>/bin/python scripts/verify/lockbox_rehearsal.py <out.json>
Delete the clone afterwards.
"""
from __future__ import annotations

import importlib.util
import json
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from part2 import assoc  # noqa: E402


def guard() -> None:
    if not (ROOT / ".lockbox_rehearsal").is_file():
        raise SystemExit("refused: no .lockbox_rehearsal marker (run only in a throw-away clone)")
    if not (ROOT / ".git").is_dir():
        raise SystemExit("refused: not a plain clone (.git is not a directory; worktrees are refused)")
    r = subprocess.run(["git", "-C", str(ROOT), "remote", "get-url", "origin"], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit("refused: no origin remote (not a clone)")


def dummy_loader(A: set, seed: int = 20261029):
    rng = np.random.default_rng(seed)
    calls = []

    def loader(model_ids, genes, source):
        ids = list(model_ids)
        if set(ids) <= A:
            return assoc._default_outcome_loader(ids, genes, source)
        cols = list(genes) if genes is not None else []
        calls.append(dict(source=source, n_models=len(ids), n_genes=len(cols)))
        vals = rng.normal(-0.3, 0.5, (len(ids), len(cols))).astype(np.float32).astype(float)
        return pd.DataFrame(vals, index=ids, columns=cols)

    return loader, calls


def main():
    guard()
    out = Path(sys.argv[1])
    spec = importlib.util.spec_from_file_location("script_51_lockbox_confirm", ROOT / "scripts" / "51_lockbox_confirm.py")
    D = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = D
    spec.loader.exec_module(D)
    t0 = time.time()
    A = set(assoc.a_models(assoc.Project(root=ROOT)))
    loader, calls = dummy_loader(A)
    world = D.real_world(outcome_loader=loader)
    t_world = time.time() - t0
    pf = D.preflight(world)
    t_pre = time.time() - t0
    if pf["gaps"]:
        raise SystemExit(f"preflight gaps: {pf['gaps']}")
    res = D.run_all(world)
    led = [json.loads(l) for l in world.proj.ledger.read_text().splitlines()]
    stage = {}
    for aid in D.confirm_order(world.proj) + [D.BENCH_ID, D.BRIDGE_ID]:
        p = world.out(aid) / "summary.json"
        stage[aid] = json.loads(p.read_text()).get("runtime_s") if p.exists() else None
    stage["BENCH_estimate"] = res["bench"].get("runtime_s")
    stage["BRIDGE_estimate"] = res["bridge"].get("runtime_s")
    rec = dict(note="dress rehearsal on real A data with a dummy L loader (random values; no lockbox value read)",
               head=subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True,
                                   text=True).stdout.strip(),
               world_s=round(t_world, 1), preflight_s=round(t_pre - t_world, 1), run_s=res["runtime_s"],
               total_s=round(time.time() - t0, 1),
               peak_rss_gb=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9, 2),
               n_ledger_lines=len(led), ledger_order=[r["analysis_id"] for r in led], dummy_loader_calls=calls,
               stage_runtime_s=stage, closed=all(D.confirm.is_closed(world.proj, r["analysis_id"]) for r in led),
               bench_primary_keys=sorted(res["bench"]["primary"]), bridge_primary_keys=sorted(res["bridge"]["primary"]),
               log=world.log)
    out.write_text(json.dumps(rec, indent=1, default=str) + "\n")
    print(json.dumps({k: v for k, v in rec.items() if k != "log"}, indent=1, default=str))


if __name__ == "__main__":
    main()
