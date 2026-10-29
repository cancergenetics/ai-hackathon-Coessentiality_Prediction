"""Score spot check (label-free): a 1% random sample (numpy default_rng(7), without replacement) of the scored test
pairs is re-scored from the raw CRISPRGeneEffect.csv through entry points other than the harness's
external/part1_locked/scripts/predict_pairs.py run on the full file, and compared with the score TSVs.

  stage 'sample'  : write the sample pair file (IDs only, from the score TSV; no label file is read)
  stage 'compare' : (1) gls/pearson/spearman via Part 1 src/baselines.py (dev-time estimators, not the predictor's
                        BaselineEngine) and via Part 1 src.predict.BaselineEngine called directly on the sample;
                    (2) zavg_BEF via the genome-wide all-pairs matrix (mvp/scripts/predict_all_pairs.py, float32);
                    (3) any extra TSVs given as name=path (mvp/scripts/predict_pairs.py, locked predict_pairs on the
                        sample file), compared on the mapped pairs.
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
P1 = ROOT / "external" / "part1_locked"
SC = ROOT / "data_cache" / "test_eval" / "scores" / "test"
ALLP = Path(os.environ.get("PART1_DIR", ROOT.parent / "coessentiality_protein complexes")) / "outputs" / "all_pairs"
VER = ROOT / "results" / "test_eval" / "verify"
stage, work = sys.argv[1], Path(sys.argv[2])
work.mkdir(parents=True, exist_ok=True)


def load(m):
    return pd.read_csv(SC / f"{m}.tsv", sep="\t", dtype={"uniprot_a": str, "uniprot_b": str})


if stage == "sample":
    t = load("zavg_BEF")[["uniprot_a", "uniprot_b"]]
    n = len(t)
    k = int(round(0.01 * n))
    idx = np.sort(np.random.default_rng(7).choice(n, size=k, replace=False))
    t.iloc[idx].to_csv(work / "sample_pairs.txt", sep="\t", header=False, index=False)
    print(json.dumps(dict(n_total=n, n_sample=k)))
    sys.exit()

sys.path.insert(0, str(P1))
from src import baselines, data, predict  # noqa: E402

samp = pd.read_csv(work / "sample_pairs.txt", sep="\t", header=None, names=["uniprot_a", "uniprot_b"], dtype=str)
key = samp.uniprot_a + "|" + samp.uniprot_b
ref = {}
for m in ("zavg_BEF", "B", "gls_signed", "spearman_signed", "pearson_signed", "spec_k20_tinf_std"):
    t = load(m)
    t.index = t.uniprot_a + "|" + t.uniprot_b
    ref[m] = t.loc[key]
mapped = ref["zavg_BEF"].mapped.values == 1
res = dict(n_sample=len(samp), n_sample_mapped=int(mapped.sum()), comparisons={})


def cmp(name, new, refm, mask):
    a, b = np.asarray(new, float)[mask], ref[refm].score.values[mask]
    d = np.abs(a - b)
    rel = d / np.maximum(np.abs(b), 1e-12)
    res["comparisons"][name] = dict(reference_tsv=refm, n=int(mask.sum()), n_nonfinite_new=int((~np.isfinite(a)).sum()),
                                    max_abs_diff=float(np.nanmax(d)), max_rel_diff=float(np.nanmax(rel)),
                                    n_exact=int((d == 0).sum()), spearman_rank_corr=float(pd.Series(a).corr(pd.Series(b), method="spearman")))
    print(name, res["comparisons"][name], flush=True)


# raw CSV, Part 1 parse + hash check (predict.load_matrix with cache=False reads the raw CSV)
X, genes, prov = predict.load_matrix(cache=False)
res["matrix"] = prov
umap = data.build_uniprot_map(pd.unique(samp.values.ravel()), genes)
ia = umap.gene_idx.reindex(samp.uniprot_a).values.astype(np.int64)
ib = umap.gene_idx.reindex(samp.uniprot_b).values.astype(np.int64)
ok = (ia >= 0) & (ib >= 0) & (ia != ib)
res["mapping_agrees_with_tsv_mapped_flag"] = bool((ok == mapped).all())

# (1) baselines
eng = predict.BaselineEngine(X)
for meth, refm in (("gls", "gls_signed"), ("pearson", "pearson_signed"), ("spearman", "spearman_signed")):
    s = np.full(len(samp), np.nan)
    s[ok] = baselines.score(meth, X, ia[ok], ib[ok])["signed"]
    cmp(f"{meth}: Part1 src/baselines.py", s, refm, mapped)
    s2 = np.full(len(samp), np.nan)
    s2[ok] = eng.compute(meth, ia[ok], ib[ok])
    cmp(f"{meth}: Part1 BaselineEngine on sample", s2, refm, mapped)

# (2) zavg_BEF from the genome-wide matrix
g = pd.read_csv(ALLP / "genes.tsv", sep="\t", dtype={"uniprot_ids": str})
u2i = {}
for i, ids in zip(g["index"], g.uniprot_ids.fillna("")):
    for u in ids.split("|"):
        if u:
            u2i.setdefault(u, []).append(i)
Z = np.load(ALLP / "zavg_BEF.npy", mmap_mode="r")
ga = np.array([u2i[u][0] if len(u2i.get(u, [])) == 1 else -1 for u in samp.uniprot_a])
gb = np.array([u2i[u][0] if len(u2i.get(u, [])) == 1 else -1 for u in samp.uniprot_b])
okz = (ga >= 0) & (gb >= 0) & (ga != gb)
sz = np.full(len(samp), np.nan)
sz[okz] = [float(Z[min(x, y), max(x, y)]) for x, y in zip(ga[okz], gb[okz])]
res["allpairs_mapped_and_tsv_mapped"] = int((okz & mapped).sum())
res["allpairs_unmapped_but_tsv_mapped"] = int((~okz & mapped).sum())
cmp("zavg_BEF: all-pairs matrix (mvp, float32)", sz, "zavg_BEF", okz & mapped)

# (3) extra TSVs
for arg in sys.argv[3:]:
    name, path = arg.split("=", 1)
    refm = name.split(":")[0]
    t = pd.read_csv(path, sep="\t", dtype={"uniprot_a": str, "uniprot_b": str})
    t.index = t.uniprot_a + "|" + t.uniprot_b
    t = t.loc[key]
    col = "score" if "score" in t.columns else t.columns[2]
    mk = mapped & (t.mapped.values == 1)
    res[f"{name}: mapped flag agrees"] = bool(((t.mapped.values == 1) == mapped).all())
    cmp(name, t[col].values, refm, mk)

(VER / "spotcheck_scores.json").write_text(json.dumps(res, indent=1))
