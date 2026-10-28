"""FZ-C construction glue on DepMap blocks (DESIGN v3.1 §2.3; decisions in results/part2/construction/DECISIONS.md).

Inputs are only the cached S x S score blocks of the construction runs (src/part2/runsets.py: A, 20 PatientID-
blocked 80% subsamples, 4 split halves), the gene filter on A, HGNC location and Ensembl coordinates (locus flag).
No physical-interaction or protein-abundance data and no context label is read here.

* zavg_BEF graphs: score >= z0 AND mutual top-k, weight min(s, z0+3) - z0 + 1.
* GLS graphs (bridge comparator, round-2 algorithms MUST 5): same S and k; per run the GLS graph keeps the top
  mutual-top-k pairs up to the edge count of the zavg graph of the *same* model set; weights quantile-matched to
  that zavg graph's weights. Cohesion for GLS uses the GLS block rank-mapped onto the zavg block of the same model
  set (decision D4), so 'mean score >= z0 (+1)' means the same quantile in both maps.
"""
from __future__ import annotations

import multiprocessing as mp
import time

import numpy as np

from part2 import consensus as C
from part2 import runsets, scores
from part2.graph import Graph, build_graph, matched_graph


def load(name: str, method: str, S=None) -> scores.Block:
    """Cached block of a construction run (never computed here: run scripts/34_run_A_blocks.sh first)."""
    S = runsets.S_mask() if S is None else S
    ids = runsets.models(name)
    p = scores.cache_path(method, ids, S)
    if not p.exists():
        raise FileNotFoundError(f"missing block {p.name} for {name}/{method}; run scripts/34_run_A_blocks.sh")
    return scores.load_block(p)


def rank_map(src: np.ndarray, ref: np.ndarray) -> np.ndarray:
    """Symmetric block whose upper-triangle finite values are src's, replaced rank-for-rank by the quantiles of ref's
    upper-triangle finite values (monotone; NaN kept)."""
    n = len(src)
    iu = np.triu_indices(n, 1)
    v = src[iu].astype(np.float64)
    ok = np.isfinite(v)
    r = np.sort(ref[iu][np.isfinite(ref[iu])].astype(np.float64))
    o = np.argsort(v[ok], kind="stable")
    q = (np.arange(len(o)) + 0.5) / len(o)
    mapped = np.empty(len(o))
    mapped[o] = np.interp(q, (np.arange(len(r)) + 0.5) / len(r), r)
    out = np.full(len(v), np.nan)
    out[np.flatnonzero(ok)] = mapped
    B = np.full((n, n), np.nan, np.float32)
    B[iu] = out
    B[(iu[1], iu[0])] = out
    return B


def gls_matched(gls_vals: np.ndarray, zg: Graph, k: int) -> tuple[Graph, dict]:
    g, thr, att = matched_graph(gls_vals, k, zg.m, zg.w)
    return g, dict(target_edges=int(zg.m), edges=int(g.m), gls_threshold=thr, attainable=bool(att))


# ------------------------------------------------------------------------------------------ parallel stage 1
def _level_job(args):
    graphs, gamma, cfg = args
    M = C.stage1(graphs, gamma, cfg)
    r, c, f = C.coclustering(M)
    return C.Level(gamma, M, r, c, f)


def build_levels_parallel(graphs, cfg: C.ClusterConfig, n_proc: int = 4) -> list[C.Level]:
    """Same result as consensus.build_levels (Leiden seeds fixed per run; one gamma per process)."""
    if n_proc <= 1:
        return C.build_levels(graphs, cfg)
    with mp.get_context("spawn").Pool(n_proc) as pool:
        return pool.map(_level_job, [(graphs, g, cfg) for g in cfg.gammas])


# ------------------------------------------------------------------------------------------ maps
def build_map(graphs, ref, z0, cfg: C.ClusterConfig, n_proc: int = 4) -> tuple[C.ClusterMap, list, dict]:
    t = time.perf_counter()
    levels = build_levels_parallel(graphs, cfg, n_proc)
    t1 = time.perf_counter()
    cm = C.build_map(graphs, cfg, ref, z0, levels=levels)
    t2 = time.perf_counter()
    return cm, levels, dict(stage1_s=round(t1 - t, 1), persistence_nesting_s=round(t2 - t1, 1))


def zavg_graphs(names, k, z0, S=None):
    out, meta = [], []
    for n in names:
        b = load(n, "zavg_BEF", S)
        out.append(build_graph(b.vals, k, z0))
        meta.append(dict(run=n, model_set_hash=b.meta["model_set_hash"], n_models=b.meta["n_models"],
                         edges=int(out[-1].m), null_moments=b.meta.get("null_moments")))
        del b
    return out, meta


def gls_graphs(names, k, zgraphs, S=None):
    out, meta = [], []
    for n, zg in zip(names, zgraphs):
        b = load(n, "gls", S)
        g, info = gls_matched(b.vals, zg, k)
        out.append(g)
        meta.append(dict(run=n, model_set_hash=b.meta["model_set_hash"], n_models=b.meta["n_models"], **info))
        del b
    return out, meta
