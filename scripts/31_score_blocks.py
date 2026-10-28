#!/usr/bin/env python
"""Run the lean locked scorer (zavg_BEF) or Part 1's GLS on one FZ-C construction run set and cache the S x S block.

Run sets (src/part2/runsets.py; all are subsets of the construction set A, verified against the FZ-LB manifest):
  A, boot1000..boot1019 (PatientID-blocked 80%), half2000a/b..half2003a/b (PatientID-blocked, lineage-stratified).
There is deliberately no all-model / U / lockbox option: Part 2 construction sees A only (DESIGN v3.1 §0 F3/F5);
``assert_A_only`` re-checks every model ID against the manifest-verified A list before scoring (hard guard).
(The Part 1 bit-exact reproduction on all models lives in tests/test_part2_scores.py, which calls the scorer directly.)
Gene set S: dependency-bearing filter computed on A. The A run of zavg_BEF also stores the (non-S genes) x S
rectangle for the attachment rule (--rect, default on for set A).
Writes data_cache/part2/blocks/<method>_<hash>_S<hash>.npz/.json and appends a line to
results/part2/construction/runs/timings.jsonl (runtime, stage timings, peak RSS of this process = this run).
No hu.MAP file, label, context annotation or lockbox outcome is read.

  PYTHONHASHSEED=0 .venv/bin/python scripts/31_score_blocks.py --set boot1000 --method zavg_BEF
"""
import argparse
import json
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import assoc, consensus, part1, paths, runsets, scores  # noqa: E402


class NotAModelsError(SystemExit):
    pass


def assert_A_only(ids) -> list:
    """Hard guard (FZ-X): refuse to score any model outside the construction set A, read from A_models.txt verified
    against the committed FZ-LB manifest (assoc.a_models), whatever the run-set code returned."""
    A = set(assoc.a_models())
    ids = [str(i) for i in ids]
    bad = sorted(set(ids) - A)
    if bad or not ids:
        raise NotAModelsError(f"31_score_blocks: {len(bad)} model(s) outside the construction set A "
                              f"(e.g. {bad[:3]}); only A run sets may be scored" if bad else "empty model set")
    return ids


def main():
    consensus.ensure_hashseed()
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", default="A", choices=runsets.run_names())
    ap.add_argument("--method", default="zavg_BEF", choices=list(scores.METHODS))
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    lock = part1.assert_locked()
    ids = assert_A_only(runsets.models(a.set))
    S = runsets.S_mask()
    rect = a.set == "A" and a.method == "zavg_BEF"
    blk = scores.score_block(ids, a.method, S, force=a.force, rect=rect)
    m = blk.meta
    rec = dict(set=a.set, method=a.method, n_models=m["n_models"], n_S=m["n_S"], model_set_hash=m["model_set_hash"],
               S_hash=m["S_hash"], runtime_s=m["runtime_s"], stages_s=m["stages_s"], peak_rss_gb=m["peak_rss_gb"],
               null_moments=m.get("null_moments"), rect=blk.rect is not None, part1_lock=lock["commit"],
               cache=str(scores.cache_path(a.method, ids, S).relative_to(paths.PROJECT)))
    out = paths.RESULTS / "construction" / "runs" / "timings.jsonl"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "a") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
