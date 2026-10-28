#!/usr/bin/env python
"""Run the lean locked scorer (zavg_BEF) or Part 1's GLS on one model set and cache the S x S block.

Model sets:  all   = all CRISPRGeneEffect models (Part 1's data; for the bit-exact reproduction check)
             U     = universe (minus engineered / drug-adapted models, Model.csv EngineeredModel)
             sub   = PatientID-blocked subsample of U, --frac (0.8) with --seed
Gene set S:  dependency-bearing filter computed on --S-from (all | U); the construction set A is not frozen yet, so
             U stands in for A in the timing runs.
Writes data_cache/part2/blocks/<method>_<hash>_S<hash>.npz/.json and appends a line to results/part2/infra/timings.jsonl.
No hu.MAP file, label or context annotation is read.

  PYTHONHASHSEED=0 .venv/bin/python scripts/31_score_blocks.py --set U --method zavg_BEF
"""
import argparse
import json
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import consensus, data, genefilter, part1, paths, scores  # noqa: E402


def model_set(name, frac, seed):
    ge = data.load_gene_effect()
    if name == "all":
        return list(ge.models)
    U = data.universe_models()
    if name == "U":
        return U
    if name == "sub":
        return data.blocked_subsample(U, frac, seed)
    raise ValueError(name)


def main():
    consensus.ensure_hashseed()
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", default="U", choices=["all", "U", "sub"])
    ap.add_argument("--method", default="zavg_BEF", choices=list(scores.METHODS))
    ap.add_argument("--frac", type=float, default=0.8)
    ap.add_argument("--seed", type=int, default=1000)
    ap.add_argument("--S-from", dest="s_from", default="U", choices=["all", "U"])
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    lock = part1.assert_locked()
    ids = model_set(a.set, a.frac, a.seed)
    S = genefilter.dependency_bearing_for_models(model_set(a.s_from, a.frac, a.seed))
    blk = scores.score_block(ids, a.method, S, force=a.force)
    m = blk.meta
    rec = dict(set=a.set, method=a.method, frac=a.frac if a.set == "sub" else 1.0,
               seed=a.seed if a.set == "sub" else None, S_from=a.s_from, n_models=m["n_models"], n_S=m["n_S"],
               model_set_hash=m["model_set_hash"], runtime_s=m["runtime_s"], stages_s=m["stages_s"],
               peak_rss_gb=m["peak_rss_gb"], null_moments=m.get("null_moments"), part1_lock=lock["commit"],
               cache=str(scores.cache_path(a.method, ids, S).relative_to(paths.PROJECT)))
    out = paths.RESULTS / "infra" / "timings.jsonl"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "a") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
