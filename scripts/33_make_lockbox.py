#!/usr/bin/env python
"""FZ-LB: draw the lockbox L (25% of U, seed 3000, PatientID-blocked, multilabel iterative stratification by
OncotreeLineage + the DESIGN §2.1 genotype list; algorithm in src/part2/lockbox.py). Model labels only: reads the
context table (verified against results/part2/contexts/checksums.json) and Model.csv PatientID; no outcome is read.

Outputs (results/part2/lockbox/, committed): lockbox_models.txt, A_models.txt (sorted ModelIDs, one per line),
strata.tsv (per-stratum n_U / n_A / n_L / target), manifest.json (SHA-256 of both lists and of the inputs, seed,
algorithm, counts). --out DIR writes elsewhere (determinism check); --verify re-draws in memory and compares with
the committed lists.

  PYTHONHASHSEED=0 .venv/bin/python scripts/33_make_lockbox.py
"""
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from part2 import contexts, data, lockbox, paths  # noqa: E402


def sha_text(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def draw():
    ck = json.loads((paths.CONTEXT_RESULTS / "checksums.json").read_text())
    tab_p = paths.CONTEXT_CACHE / "context_table.parquet"
    rel = str(tab_p.relative_to(paths.PROJECT))
    h = data.sha256_file(tab_p)
    if ck["files"].get(rel) != h:
        raise RuntimeError(f"{rel} SHA-256 {h} != checksums.json {ck['files'].get(rel)}")
    ctx = contexts.load_table(tab_p)
    U = data.universe_model_ids()
    if sorted(U) != sorted(ctx.index):
        raise RuntimeError("context table index != U")
    labels = lockbox.label_matrix(ctx)
    blocks = data.load_models().PatientID.reindex(labels.index)
    gcols = list(lockbox.GENOTYPE_LABELS)
    fold = lockbox.iterative_stratification(labels, blocks, lockbox.FRAC_L, lockbox.SEED, genotype_cols=gcols)
    strata = lockbox.stratum_table(labels, fold, gcols)
    return fold, strata, dict(context_table=rel, context_table_sha256=h, context_content_sha256=ck["content_sha256"]
                              ["context_table.parquet"])


def size_class(n_A: int) -> str:
    return "confirmatory" if n_A >= 12 else ("exploratory" if n_A >= 8 else "descriptive")


def context_counts(fold):
    """Per binary context of the context table: positives (n_*) and reference (ref_*) models in U / A / L, the
    DESIGN §2.2 size class in A and whether n_L >= 5 (lockbox-testable size). Labels only."""
    ctx = contexts.load_table()
    summ = {r.context: r.kind for r in pd.read_csv(paths.CONTEXT_RESULTS / "context_summary.tsv", sep="\t").itertuples()}
    rows = []
    for c, kind in summ.items():
        if kind != "binary":
            continue
        v = ctx[c].reindex(fold.index).values.astype(float)
        rec = dict(context=c)
        for lab, m in (("U", np.ones(len(v), bool)), ("A", (fold == "A").values), ("L", (fold == "L").values)):
            rec[f"n_{lab}"] = int(np.sum((v == 1) & m))
            rec[f"ref_{lab}"] = int(np.sum((v == 0) & m))
        rec["size_class_A"] = size_class(rec["n_A"])
        rec["nL_ge5"] = rec["n_L"] >= 5
        rows.append(rec)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(paths.LOCKBOX_DIR))
    ap.add_argument("--verify", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    fold, strata, inp = draw()
    L = sorted(fold.index[fold == "L"])
    A = sorted(fold.index[fold == "A"])
    txtL, txtA = "\n".join(L) + "\n", "\n".join(A) + "\n"
    if a.verify:
        okL = paths.LOCKBOX_MODELS.read_text() == txtL
        okA = paths.A_MODELS.read_text() == txtA
        m = json.loads(paths.LOCKBOX_MANIFEST.read_text())
        okM = m["files"]["lockbox_models.txt"] == sha_text(txtL) and m["files"]["A_models.txt"] == sha_text(txtA)
        print(json.dumps(dict(lockbox_identical=okL, A_identical=okA, manifest_match=okM)))
        sys.exit(0 if okL and okA and okM else 1)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    counts = context_counts(fold)
    counts.to_csv(out / "context_counts.tsv", sep="\t", index=False)
    (out / "lockbox_models.txt").write_text(txtL)
    (out / "A_models.txt").write_text(txtA)
    strata.to_csv(out / "strata.tsv", sep="\t", index=False)
    blocks = data.load_models().PatientID.reindex(fold.index)
    straddle = int(blocks.groupby(blocks).apply(lambda s: fold[s.index].nunique()).gt(1).sum())
    g = strata[strata.kind == "genotype"]
    man = dict(
        freeze="FZ-LB", design="DESIGN.md v3.1 §2.1", seed=lockbox.SEED, frac_L=lockbox.FRAC_L,
        algorithm="iterative stratification (Sechidis et al. 2011) on PatientID blocks; see src/part2/lockbox.py",
        genotype_labels=lockbox.GENOTYPE_LABELS, n_U=len(fold), n_L=len(L), n_A=len(A),
        n_patients_straddling=straddle, universe_hash=data.model_set_hash(fold.index),
        files={"lockbox_models.txt": sha_text(txtL), "A_models.txt": sha_text(txtA),
               "strata.tsv": data.sha256_file(out / "strata.tsv"),
               "context_counts.tsv": data.sha256_file(out / "context_counts.tsv")},
        inputs=dict(inp, model_csv_sha256=data.sha256_file(paths.MODEL_CSV)),
        balance=dict(max_abs_dev_L_genotype=int(g.dev_L.abs().max()),
                     max_abs_dev_L_lineage=int(strata[strata.kind == "lineage"].dev_L.abs().max()),
                     genotypes_with_nL_lt5=g.loc[g.n_L < 5, "stratum"].tolist()),
        runtime_s=round(time.time() - t0, 2))
    (out / "manifest.json").write_text(json.dumps(man, indent=1) + "\n")
    print(json.dumps({k: man[k] for k in ("n_U", "n_L", "n_A", "n_patients_straddling", "balance", "files")}, indent=1))


if __name__ == "__main__":
    main()
