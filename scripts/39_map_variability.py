#!/usr/bin/env python
"""Map-to-map variability of the FZ-C complex maps over the 20 PatientID-blocked bootstrap subsamples (DESIGN §2.3,
"reported descriptively"). Purely descriptive; nothing here feeds a decision.

Inputs: the cached S x S score blocks of the construction runs (boot1000-1019, A; scripts/34_run_A_blocks.sh), the
frozen maps results/part2/construction/{zavg_map,gls_map}/clusters.tsv (read only) and
selection/selection.json (z0, k). A model set only; no validation data, external complex resources or lockbox.

Bootstrap map b = the frozen consensus + persistence + nesting pipeline (construct.build_map, ClusterConfig of
scripts/36_build_maps.py, cohesion reference = the full-A block as in the frozen build) applied to the single graph of
subsample b with the same 5 Leiden seeds, at the selected (z0, k). GLS: per-bootstrap graph edge-count-matched to the
zavg graph of the same subsample, cohesion on the GLS A block rank-mapped to the zavg A block (D4).

  stage maps    : build the 20 + 20 bootstrap maps (one process per subsample), write level-1 labels
  stage analyse : (a) best-match Jaccard of every frozen level-1 complex in every bootstrap map (by size bin) plus the
                  size-matched random-gene-set chance reference; (b)/(c) pairwise ARI/NMI of bootstrap maps
  stage verify  : (optional) rebuild the frozen zavg map from all 20 graphs and compare with the frozen clusters.tsv

  PYTHONHASHSEED=0 .venv/bin/python scripts/39_map_variability.py --stage all --n-proc 4
"""
import argparse
import json
import multiprocessing as mp
import resource
import sys
import time
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score  # noqa: E402

from part2 import construct, paths, runsets, scores  # noqa: E402
from part2 import consensus as C  # noqa: E402
from part2.graph import build_graph  # noqa: E402

CON = paths.RESULTS / "construction"
OUT = CON / "variability"
CCFG = C.ClusterConfig(support_mode="containment", placement="nested")      # as scripts/36_build_maps.py
BINS = [(3, 4), (5, 8), (9, 16), (17, 60)]
BIN_NAMES = ["3-4", "5-8", "9-16", "17-60"]
N_RANDOM = 100              # random gene sets per frozen complex
RANDOM_SEED = 7000
METHODS = ("zavg_BEF", "gls")
MAPDIR = {"zavg_BEF": "zavg_map", "gls": "gls_map"}


def bin_of(size):
    for name, (lo, hi) in zip(BIN_NAMES, BINS):
        if lo <= size <= hi:
            return name
    return None


def peak_gb_children():
    r = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return r / 1024 ** 3 if sys.platform == "darwin" else r / 1024 ** 2


# ------------------------------------------------------------------------------------------ stage maps
def _boot_job(args):
    """Both maps of one subsample. Returns level-1 labels (S-local, -1 = uncovered) and counts."""
    seed, z0, k, S = args
    name = f"boot{seed}"
    t0 = time.perf_counter()
    zb = construct.load(name, "zavg_BEF", S)
    zg = build_graph(zb.vals, k, z0)
    nz = int(zg.m)
    del zb
    refA = construct.load("A", "zavg_BEF", S).vals
    out = dict(seed=seed)
    # zavg
    cm, _, tim = construct.build_map([zg], refA, z0, CCFG, 1)
    out["zavg_BEF"] = dict(labels=cm.labels(zg.n, 1), n1=len(cm.level(1)), n2=len(cm.level(2)), edges=nz, **tim)
    # GLS (edge-matched to this subsample's zavg graph)
    (gg,), (gm,) = construct.gls_graphs([name], k, [zg], S)
    glsA = construct.load("A", "gls", S).vals
    ref = construct.rank_map(glsA, refA)
    del glsA, refA
    cm, _, tim = construct.build_map([gg], ref, z0, CCFG, 1)
    out["gls"] = dict(labels=cm.labels(zg.n, 1), n1=len(cm.level(1)), n2=len(cm.level(2)), edges=int(gg.m),
                      attainable=gm["attainable"], **tim)
    out["wall_s"] = round(time.perf_counter() - t0, 1)
    out["peak_rss_gb"] = round(scores.peak_rss_gb(), 2)
    print(f"boot{seed}: zavg L1={out['zavg_BEF']['n1']} gls L1={out['gls']['n1']} {out['wall_s']}s "
          f"peak_rss={out['peak_rss_gb']}GB", flush=True)
    return out


def stage_maps(S, z0, k, n_proc):
    t0 = time.perf_counter()
    jobs = [(s, z0, k, S) for s in runsets.BOOT_SEEDS]
    with mp.get_context("spawn").Pool(n_proc) as pool:
        res = pool.map(_boot_job, jobs, chunksize=1)
    OUT.mkdir(parents=True, exist_ok=True)
    sidx = np.flatnonzero(S)
    for m in METHODS:
        df = pd.DataFrame({f"boot{r['seed']}": r[m]["labels"] for r in res})
        df.insert(0, "s_index", np.arange(len(sidx)))
        df.insert(1, "gene_index", sidx)
        df.to_csv(OUT / f"{m}_bootstrap_level1_labels.tsv.gz", sep="\t", index=False)
    meta = dict(
        z0=z0, k=k, n_boot=len(res), boot_seeds=list(runsets.BOOT_SEEDS), leiden_seeds=list(CCFG.seeds),
        per_boot=[dict(seed=r["seed"], wall_s=r["wall_s"], peak_rss_gb=r["peak_rss_gb"],
                       **{m: {kk: vv for kk, vv in r[m].items() if kk != "labels"} for m in METHODS}) for r in res])
    meta["maps_wall_s"] = round(time.perf_counter() - t0, 1)
    meta["maps_peak_rss_gb_worker_max"] = max(r["peak_rss_gb"] for r in res)
    meta["maps_peak_rss_gb_children_max"] = round(peak_gb_children(), 2)
    meta["n_proc"] = n_proc
    (OUT / "maps_runtime.json").write_text(json.dumps(meta, indent=1, default=float))
    print("maps stage", meta["maps_wall_s"], "s; worker peak", meta["maps_peak_rss_gb_worker_max"], "GB", flush=True)


# ------------------------------------------------------------------------------------------ stage analyse
def load_frozen(method, S):
    """Frozen level-1 complexes as S-local index arrays."""
    d = CON / MAPDIR[method]
    cl = pd.read_csv(d / "clusters.tsv", sep="\t")
    mem = pd.read_csv(d / "membership.tsv", sep="\t")
    sidx = np.flatnonzero(S)
    loc = np.full(int(S.shape[0]), -1)
    loc[sidx] = np.arange(len(sidx))
    l1 = cl[cl.level == 1].reset_index(drop=True)
    m1 = mem[mem.level == 1]
    by = {c: loc[g.values] for c, g in m1.groupby("cluster_id").gene_index}
    assert all((v >= 0).all() for v in by.values())
    members = [np.sort(by[c]) for c in l1.cluster_id]
    assert [len(m) for m in members] == l1["size"].tolist()
    return l1, members


def q(x, p):
    return float(np.percentile(x, p)) if len(x) else float("nan")


def jaccard_matrix(members, labels_mat):
    """len(members) x n_maps best-match Jaccard against each map's level-1 partition."""
    out = np.zeros((len(members), labels_mat.shape[1]))
    for b in range(labels_mat.shape[1]):
        lab = labels_mat[:, b]
        sizes = np.bincount(lab[lab >= 0])
        for i, m in enumerate(members):
            out[i, b] = C.best_jaccard(m, lab, sizes)
    return out


def summarise(x):
    x = np.asarray(x).ravel()
    return dict(n=int(len(x)), median=q(x, 50), q25=q(x, 25), q75=q(x, 75), frac_ge_0_5=float((x >= 0.5).mean()),
                mean=float(x.mean()) if len(x) else float("nan"))


def analyse_method(method, S, rng):
    l1, members = load_frozen(method, S)
    lab_df = pd.read_csv(OUT / f"{method}_bootstrap_level1_labels.tsv.gz", sep="\t")
    cols = [c for c in lab_df.columns if c.startswith("boot")]
    L = lab_df[cols].values.astype(int)
    nS = L.shape[0]
    sizes = np.array([len(m) for m in members])
    binlab = np.array([bin_of(s) for s in sizes])
    Jobs = jaccard_matrix(members, L)                                    # complexes x boots
    # chance reference: size-matched random sets from S, and from the genes the frozen map covers
    covered = np.unique(np.concatenate(members))
    refs = {"S": np.arange(nS), "frozen_covered": covered}
    Jrand = {}
    for rname, pool in refs.items():
        Jr = np.zeros((len(members), N_RANDOM, L.shape[1]))
        for i, s in enumerate(sizes):
            sets = [rng.choice(pool, size=s, replace=False) for _ in range(N_RANDOM)]
            for b in range(L.shape[1]):
                lab = L[:, b]
                sz = np.bincount(lab[lab >= 0])
                for r_, ms in enumerate(sets):
                    Jr[i, r_, b] = C.best_jaccard(ms, lab, sz)
        Jrand[rname] = Jr
    # per-complex table
    per = l1[["cluster_id", "size", "locus", "stability", "support"]].copy()
    per["size_bin"] = binlab
    per["jaccard_median"] = np.median(Jobs, 1)
    per["jaccard_q25"] = np.percentile(Jobs, 25, axis=1)
    per["jaccard_q75"] = np.percentile(Jobs, 75, axis=1)
    per["jaccard_min"] = Jobs.min(1)
    per["frac_ge_0_5"] = (Jobs >= 0.5).mean(1)
    for rname, Jr in Jrand.items():
        thr = np.array([np.percentile(Jr[i], 95) for i in range(len(members))])
        per[f"random_{rname}_p95"] = thr
        per[f"median_above_random_{rname}_p95"] = per.jaccard_median.values > thr
    per.to_csv(OUT / f"{method}_per_complex_jaccard.tsv", sep="\t", index=False, float_format="%.4f")
    # per-bin summaries
    bins = {}
    for name in BIN_NAMES + ["all"]:
        sel = np.ones(len(members), bool) if name == "all" else binlab == name
        if not sel.any():
            continue
        d = dict(n_complexes=int(sel.sum()), observed_pooled=summarise(Jobs[sel]),
                 observed_per_complex_median=summarise(np.median(Jobs[sel], 1)),
                 observed_per_complex_frac_ge_0_5_median=float(np.median((Jobs[sel] >= 0.5).mean(1))))
        for rname, Jr in Jrand.items():
            d[f"random_{rname}_pooled"] = summarise(Jr[sel])
            d[f"random_{rname}_p95"] = q(Jr[sel].ravel(), 95)
            d[f"frac_complexes_median_above_random_{rname}_p95"] = float(per[f"median_above_random_{rname}_p95"].values[sel].mean())
        bins[name] = d
    # (b)/(c) pairwise ARI / NMI between bootstrap maps (level 1)
    pair_rows = []
    cov = L >= 0
    for a, b in combinations(range(L.shape[1]), 2):
        both = cov[:, a] & cov[:, b]
        either = cov[:, a] | cov[:, b]
        la, lb = L[both, a], L[both, b]
        ua = np.where(cov[either, a], L[either, a], -1 - np.arange(either.sum()))      # uncovered = own singleton
        ub = np.where(cov[either, b], L[either, b], -1 - np.arange(either.sum()))
        pair_rows.append(dict(map_a=cols[a], map_b=cols[b], n_both=int(both.sum()), n_either=int(either.sum()),
                              ari_both=adjusted_rand_score(la, lb), nmi_both=normalized_mutual_info_score(la, lb),
                              ari_either_singletons=adjusted_rand_score(ua, ub),
                              nmi_either_singletons=normalized_mutual_info_score(ua, ub)))
    pairs = pd.DataFrame(pair_rows)
    pairs.to_csv(OUT / f"{method}_pairwise_ari_nmi.tsv", sep="\t", index=False, float_format="%.4f")
    # each bootstrap map vs the frozen map (covered in both)
    flab = np.full(nS, -1)
    for i, m in enumerate(members):
        flab[m] = i
    vs = []
    for b in range(L.shape[1]):
        both = (L[:, b] >= 0) & (flab >= 0)
        vs.append(dict(map=cols[b], n_both=int(both.sum()), ari=adjusted_rand_score(flab[both], L[both, b]),
                       nmi=normalized_mutual_info_score(flab[both], L[both, b]),
                       n_clusters=int(len(np.unique(L[L[:, b] >= 0, b]))), n_genes_covered=int(cov[:, b].sum())))
    vs = pd.DataFrame(vs)
    vs.to_csv(OUT / f"{method}_bootstrap_vs_frozen.tsv", sep="\t", index=False, float_format="%.4f")
    pw = {c: summarise(pairs[c]) for c in ["ari_both", "nmi_both", "ari_either_singletons", "nmi_either_singletons"]}
    return dict(n_frozen_complexes=len(members), n_frozen_genes=int(len(covered)), size_bins=bins, pairwise=pw,
                n_pairs=len(pairs), bootstrap_vs_frozen=dict(ari=summarise(vs.ari), nmi=summarise(vs.nmi)),
                bootstrap_maps=dict(n_clusters=summarise(vs.n_clusters), n_genes_covered=summarise(vs.n_genes_covered)))


def stage_analyse(S):
    t0 = time.perf_counter()
    rng = np.random.default_rng(RANDOM_SEED)
    res = {m: analyse_method(m, S, rng) for m in METHODS}
    res["settings"] = dict(size_bins=BIN_NAMES, n_random_sets_per_complex=N_RANDOM, random_seed=RANDOM_SEED,
                           random_pools=dict(S="uniform over S (6,022 genes)",
                                             frozen_covered="uniform over genes in the frozen map's level-1 complexes"))
    res["analyse_wall_s"] = round(time.perf_counter() - t0, 1)
    res["analyse_peak_rss_gb"] = round(scores.peak_rss_gb(), 2)
    (OUT / "variability.json").write_text(json.dumps(res, indent=1, default=float))
    print("analyse", res["analyse_wall_s"], "s", flush=True)
    return res


# ------------------------------------------------------------------------------------------ stage verify
def stage_verify(S, z0, k, n_proc):
    """Pipeline check: rebuilding the zavg map from all 20 graphs must give the frozen level-1 member sets."""
    t0 = time.perf_counter()
    names = [f"boot{s}" for s in runsets.BOOT_SEEDS]
    zg, _ = construct.zavg_graphs(names, k, z0, S)
    refA = construct.load("A", "zavg_BEF", S).vals
    cm, _, _ = construct.build_map(zg, refA, z0, CCFG, n_proc)
    mine = {tuple(sorted(c["members"].tolist())) for c in cm.level(1)}
    _, fm = load_frozen("zavg_BEF", S)
    frozen = {tuple(m.tolist()) for m in fm}
    out = dict(n_rebuilt=len(mine), n_frozen=len(frozen), identical=mine == frozen, n_common=len(mine & frozen),
               wall_s=round(time.perf_counter() - t0, 1))
    (OUT / "verify_frozen_rebuild.json").write_text(json.dumps(out, indent=1))
    print("verify", out, flush=True)


def main():
    C.ensure_hashseed()
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["maps", "analyse", "verify", "all"], default="all")
    ap.add_argument("--n-proc", type=int, default=4)
    a = ap.parse_args()
    sel = json.loads((CON / "selection" / "selection.json").read_text())["zavg_BEF"]["selected"]
    z0, k = sel["z0"], sel["k"]
    assert (z0, k) == (6.0, 10), sel
    S = runsets.S_mask()
    OUT.mkdir(parents=True, exist_ok=True)
    if a.stage in ("maps", "all"):
        stage_maps(S, z0, k, a.n_proc)
    if a.stage in ("analyse", "all"):
        stage_analyse(S)
    if a.stage == "verify":
        stage_verify(S, z0, k, a.n_proc)


if __name__ == "__main__":
    main()
