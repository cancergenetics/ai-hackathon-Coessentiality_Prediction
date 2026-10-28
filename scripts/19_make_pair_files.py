"""Build the predictor's ID-only input files (and, for the dry run, the development label files).

  --set train   all hu.MAP TRAIN pairs (Part 1 development labels, read-only Part 1 folder)
                -> data_cache/test_eval/ids/train.tsv
  --set dryrun  Part 1's dev / IFT split of the train pairs (canonical Part 1 protein assignment, hash-checked):
                both proteins on the dev side / both on the IFT side (cross pairs dropped, as Part 1)
                -> data_cache/test_eval/dryrun/{dev,ift}_{positive,negative}_pairs.txt + ids/{dryrun_dev,dryrun_ift}.tsv
  --set test    hu.MAP TEST pairs (only after docs/TEST_PROTOCOL.md is committed) -> data_cache/test_eval/ids/test.tsv
ID files hold unique canonical non-self pairs of pos ∪ neg, sorted, no labels.
"""
import argparse
import hashlib
import io
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from testeval import pairs as P  # noqa: E402
from testeval import paths  # noqa: E402
from testeval.protocol import sha256_file  # noqa: E402


def load_assignment() -> pd.DataFrame:
    b = paths.P1_ASSIGNMENT.read_bytes()
    want = paths.P1_ASSIGNMENT_SHA.read_text().split()[0]
    got = hashlib.sha256(b).hexdigest()
    if got != want:
        raise RuntimeError(f"Part 1 protein assignment hash mismatch {got} != {want}")
    return pd.read_csv(io.BytesIO(b))


def write_labels(t: pd.DataFrame, pos_path: Path, neg_path: Path):
    pos_path.parent.mkdir(parents=True, exist_ok=True)
    t[t.label == 1][["a", "b"]].to_csv(pos_path, sep="\t", header=False, index=False, lineterminator="\n")
    t[t.label == 0][["a", "b"]].to_csv(neg_path, sep="\t", header=False, index=False, lineterminator="\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", required=True, choices=["train", "dryrun", "test"])
    a = ap.parse_args()
    ids = paths.CACHE / "ids"
    rec = dict(set=a.set)
    if a.set == "train":
        t, audit = P.labelled_table(paths.TRAIN_POS, paths.TRAIN_NEG)
        rec.update(inputs={str(p): sha256_file(p) for p in (paths.TRAIN_POS, paths.TRAIN_NEG)}, audit=audit,
                   ids=P.write_id_file([t], ids / "train.tsv"))
    elif a.set == "dryrun":
        t, audit = P.labelled_table(paths.TRAIN_POS, paths.TRAIN_NEG)
        if audit["self_pairs"] or audit["conflicts"]:
            raise RuntimeError("train files changed: self pairs / conflicts present")
        asg = load_assignment()
        side = dict(zip(asg.protein, asg.side))
        sa, sb = t.a.map(side).values, t.b.map(side).values
        t["side"] = np.where(sa == sb, sa, "cross")
        d = paths.CACHE / "dryrun"
        rec.update(inputs={str(p): sha256_file(p) for p in (paths.TRAIN_POS, paths.TRAIN_NEG, paths.P1_ASSIGNMENT)},
                   counts={s: dict(pos=int(((t.side == s) & (t.label == 1)).sum()),
                                   neg=int(((t.side == s) & (t.label == 0)).sum())) for s in ("dev", "ift", "cross")})
        for s in ("dev", "ift"):
            sub = t[t.side == s]
            write_labels(sub, d / f"{s}_positive_pairs.txt", d / f"{s}_negative_pairs.txt")
            rec[f"ids_{s}"] = P.write_id_file([sub], ids / f"dryrun_{s}.tsv")
    else:
        from testeval.protocol import require_frozen
        rec["git"] = require_frozen()
        pos, neg = paths.test_files()
        t, audit = P.labelled_table(pos, neg)                 # label-agnostic union is written; audit = counts only
        rec.update(inputs={str(p): sha256_file(p) for p in (pos, neg)}, ids=P.write_id_file([t], ids / "test.tsv"))
    out = paths.CACHE / f"ids_{a.set}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rec, indent=1, default=str))
    print(json.dumps({k: v for k, v in rec.items() if k != "inputs"}, indent=1, default=str))


if __name__ == "__main__":
    main()
