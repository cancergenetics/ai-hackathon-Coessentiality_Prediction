"""Graph from an S x S score block: edge iff score >= z0 AND mutual top-k (within S); weight capped.

weight = min(score, z0 + cap) - z0 + 1  in [1, 1 + cap]  (cap = 3 -> [1, 4]; paralog edges cannot dominate CPM).
GLS bridge (round-2 algorithms MUST 5): same S and k as the selected zavg graph, threshold chosen to reproduce the
zavg edge count (or the largest attainable count if the mutual top-k cap binds), weights mapped to [1, 1+cap] by
quantile matching to the zavg weight distribution.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from part2.scores import select_topk


@dataclass
class Graph:
    n: int                 # number of nodes (= |S|; local indices)
    r: np.ndarray          # int64, r < c
    c: np.ndarray
    v: np.ndarray          # float64 score
    w: np.ndarray          # float64 weight

    @property
    def m(self) -> int:
        return len(self.r)

    def keys(self) -> np.ndarray:
        return self.r * self.n + self.c


def mutual_knn(vals: np.ndarray, k: int, topk=None):
    """All mutual top-k pairs (r < c) of a symmetric score block with their scores (no threshold)."""
    n = len(vals)
    I, _ = select_topk(vals, k) if topk is None else (topk[0][:, :k], None)
    r = np.repeat(np.arange(n), I.shape[1])
    c = I.ravel().astype(np.int64)
    ok = c >= 0
    r, c = r[ok], c[ok]
    key = r * n + c
    rkey = c * n + r
    ks = np.sort(key)
    pos = np.searchsorted(ks, rkey)
    pos[pos == len(ks)] = 0
    mutual = (ks[pos] == rkey) & (r < c)
    r, c = r[mutual], c[mutual]
    return r, c, vals[r, c].astype(np.float64)


def cap_weight(v, z0: float, cap: float = 3.0):
    return np.minimum(v, z0 + cap) - z0 + 1.0


def build_graph(vals: np.ndarray, k: int, z0: float, cap: float = 3.0, topk=None) -> Graph:
    r, c, v = mutual_knn(vals, k, topk)
    ok = v >= z0
    return Graph(len(vals), r[ok], c[ok], v[ok], cap_weight(v[ok], z0, cap))


def matched_graph(vals: np.ndarray, k: int, target_edges: int, ref_weights: np.ndarray, topk=None):
    """Edge-count-matched graph (for GLS / Pearson): top `target_edges` mutual top-k pairs by score; weights by
    quantile matching to ref_weights. Returns (Graph, threshold, attainable)."""
    r, c, v = mutual_knn(vals, k, topk)
    o = np.lexsort((c, r, -v))
    attainable = len(o) >= target_edges
    o = o[:target_edges]
    r, c, v = r[o], c[o], v[o]
    thr = float(v.min()) if len(v) else np.inf
    q = (np.argsort(np.argsort(v, kind="stable"), kind="stable") + 0.5) / max(len(v), 1)
    w = np.quantile(np.asarray(ref_weights, float), q) if len(v) else v
    s = np.lexsort((c, r))
    return Graph(len(vals), r[s], c[s], v[s], np.asarray(w)[s]), thr, attainable


def degree_preserving_rewire(r, c, n, n_swaps_per_edge: int = 10, seed: int = 0):
    """Double-edge swaps on a simple undirected graph; edges keep their index, so per-edge attributes (e.g. per-run
    presence and weights) travel with them. Returns new (r, c) with r < c."""
    rng = np.random.default_rng(seed)
    r, c = np.asarray(r).copy(), np.asarray(c).copy()
    m = len(r)
    if m < 2:
        return r, c
    present = set((r * n + c).tolist())
    tries = n_swaps_per_edge * m
    E1, E2 = rng.integers(0, m, tries), rng.integers(0, m, tries)
    flip = rng.random(tries) < 0.5
    for e1, e2, f in zip(E1.tolist(), E2.tolist(), flip.tolist()):
        if e1 == e2:
            continue
        a, b = int(r[e1]), int(c[e1])
        x, y = (int(c[e2]), int(r[e2])) if f else (int(r[e2]), int(c[e2]))
        # (a,b),(x,y) -> (a,y),(x,b)
        if a == y or x == b:
            continue
        k1 = min(a, y) * n + max(a, y)
        k2 = min(x, b) * n + max(x, b)
        if k1 in present or k2 in present or k1 == k2:
            continue
        present.discard(a * n + b if a < b else b * n + a)
        present.discard(min(x, y) * n + max(x, y))
        present.add(k1)
        present.add(k2)
        r[e1], c[e1] = min(a, y), max(a, y)
        r[e2], c[e2] = min(x, b), max(x, b)
    return r, c


def rewire_runs(graphs: list[Graph], seed: int = 0, n_swaps_per_edge: int = 10) -> list[Graph]:
    """Round-2 SHOULD 9: rewire the union graph once and carry each edge's per-run presence / weight with it."""
    n = graphs[0].n
    allk = np.unique(np.concatenate([g.keys() for g in graphs]))
    ur, uc = allk // n, allk % n
    nr, nc = degree_preserving_rewire(ur, uc, n, n_swaps_per_edge, seed)
    out = []
    for g in graphs:
        pos = np.searchsorted(allk, g.keys())
        out.append(Graph(n, nr[pos], nc[pos], g.v.copy(), g.w.copy()))
    for g in out:                                           # canonical order
        s = np.lexsort((g.c, g.r))
        g.r, g.c, g.v, g.w = g.r[s], g.c[s], g.v[s], g.w[s]
    return out
