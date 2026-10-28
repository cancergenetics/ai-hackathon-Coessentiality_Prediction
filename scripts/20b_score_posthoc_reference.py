"""Score a pair-ID file with Part 1's post hoc PC-removed Pearson reference (spec_k20_tinf_std); see
src/testeval/posthoc.py. Usage: python scripts/20b_score_posthoc_reference.py --pairs IDS --out TSV --record JSON"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from testeval import posthoc  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--pairs", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--record")
a = ap.parse_args()
rec = posthoc.run(a.pairs, a.out, a.record)
print(json.dumps(dict(method=rec["method"], counts=rec["counts"], runtime_s=rec["runtime"]["total_s"]), indent=1))
