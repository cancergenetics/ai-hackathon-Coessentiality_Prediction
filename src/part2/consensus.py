"""Consensus + persistence clustering (DESIGN v2 §2.3 with round-2 algorithms MUST 3 specifics).

Stage 1  For every run graph (bootstrap / subsample) and seed in SEEDS: Leiden-CPM on edge weights in [1, 4] at each
         gamma_j of a geometric grid gamma_j = 0.05 * 2^(j/2), j = 0..10.
Stage 2  Per gamma_j: pair set = union of pairs co-clustered in >= 1 stage-1 run; frequency = co-clustered runs /
         all runs (a pair absent from a run counts 0). Leiden-CPM on the frequency graph with gamma_2 = 0.5.
Persist. A stage-2 cluster c at gamma_j (size >= min_size) recurs at gamma_{j+1} if its best-match Jaccard there is
         >= 0.75; support(c) = fraction of stage-1 runs at gamma_j holding a >= 0.75-Jaccard match. A chain of
         recurrences is persistent if it holds >= 3 adjacent gamma values with support >= 0.6 each; its
         representative is the member of highest support (ties: coarser gamma).
Levels   Top-down. Level 1: persistent clusters on all nodes, coarse first, one cluster per gene (a candidate
         overlapping already-placed genes is trimmed and kept only if >= 80% of it remains). A candidate violating a
         hard rule (size > max_size, or cohesion = mean score over ALL within-cluster pairs from the reference
         block < z0 + 0 for coarse clusters, < z0 + 1 for sub-modules; DESIGN v3.1) is re-clustered: stage 2 re-run on the same frequencies induced on its genes at finer
         gamma, persistence again, candidates re-enter. Level 2 (sub-modules): the same restricted re-clustering
         inside each accepted level-1 cluster at gamma finer than its representative; proper (Jaccard < 0.75 to
         the parent) sub-clusters, one per gene. Children are subsets of parents by construction (containment 1).
Counts of candidates removed by each rule are returned.

Determinism: Leiden seeds fixed, igraph graphs built from sorted edge arrays, no hash-ordered iteration; run under
PYTHONHASHSEED=0 (``ensure_hashseed``) and check in a subprocess (tests/test_part2_consensus.py).
"""
from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field

import igraph as ig
import leidenalg as la
import numpy as np

from part2.graph import Graph

SEEDS = (0, 1, 2, 3, 4)


def gamma_grid(g0: float = 0.05, n: int = 11, step: float = 2 ** 0.5) -> tuple:
    return tuple(float(g0 * step ** j) for j in range(n))


@dataclass(frozen=True)
class ClusterConfig:
    gammas: tuple = gamma_grid()
    seeds: tuple = SEEDS
    gamma2: float = 0.5
    stage2_seed: int = 0
    n_iterations: int = 2                 # leidenalg default; -1 = until stable
    jaccard_recur: float = 0.75
    support_min: float = 0.6
    support_mode: str = "jaccard"         # 'jaccard' (round-2 review: a run holds a >= jaccard_recur match) |
                                          # 'containment' (a run puts >= jaccard_recur of the members in one cluster)
    min_adjacent: int = 3
    min_size: int = 3
    max_size: int = 60
    trim_keep: float = 0.8                # one cluster per gene: keep a trimmed candidate if >= 80% remains
    cohesion_delta: float | None = 1.0    # sub-modules: mean all-pairs score >= z0 + delta (None: no rule)
    cohesion_delta_level1: float | None = 0.0   # coarse clusters: mean all-pairs score >= z0 + delta (DESIGN v3.1)
    max_depth: int = 4                    # recursion limit for re-clustering oversize / non-cohesive candidates
    placement: str = "greedy"             # 'greedy' (coarse-first, trim, drop if < trim_keep remains) |
                                          # 'nested' (FZ-C decision D2: persistence-first; containment = nesting;
                                          #  partial overlap > overlap_conflict resolved to the more persistent)
    nest_contain: float = 0.8             # 'nested': |a & b| / |b| >= this -> b nested in a
    overlap_conflict: float = 0.5         # 'nested': partial overlap > this fraction of the candidate -> dropped


def ensure_hashseed():
    """Re-exec the current script with PYTHONHASHSEED=0 if it is not set (call at the top of scripts)."""
    if os.environ.get("PYTHONHASHSEED") != "0":
        os.environ["PYTHONHASHSEED"] = "0"
        os.execv(sys.executable, [sys.executable] + sys.argv)


# ------------------------------------------------------------------------------------------ Leiden
def leiden(n, r, c, w, gamma, seed, n_iterations=2) -> np.ndarray:
    g = ig.Graph(n=n, edges=np.column_stack([r, c]).tolist(), directed=False)
    g.es["weight"] = np.asarray(w, float).tolist()
    p = la.find_partition(g, la.CPMVertexPartition, weights="weight", resolution_parameter=gamma, seed=int(seed),
                          n_iterations=n_iterations)
    return np.asarray(p.membership, dtype=np.int64)


def stage1(graphs: list[Graph], gamma: float, cfg: ClusterConfig) -> np.ndarray:
    """(n_runs = len(graphs) * len(seeds)) x n memberships."""
    return np.array([leiden(g.n, g.r, g.c, g.w, gamma, s, cfg.n_iterations) for g in graphs for s in cfg.seeds])


def coclustering(memberships: np.ndarray, merge_every: int = 10):
    """Union of co-clustered pairs over runs and their frequency (denominator = all runs). Returns (r, c, f)."""
    R, n = memberships.shape
    acc_k = np.empty(0, np.int64)
    acc_c = np.empty(0, np.int64)
    buf = []
    for t in range(R):
        m = memberships[t]
        o = np.argsort(m, kind="stable")
        ms = m[o]
        bounds = np.flatnonzero(np.r_[True, ms[1:] != ms[:-1], True])
        keys = []
        for a, b in zip(bounds[:-1], bounds[1:]):
            if b - a < 2:
                continue
            nodes = np.sort(o[a:b])
            i, j = np.triu_indices(len(nodes), 1)
            keys.append(nodes[i] * n + nodes[j])
        if keys:
            buf.append(np.concatenate(keys))
        if len(buf) >= merge_every or t == R - 1:
            if buf:
                allk = np.concatenate([acc_k] + buf)
                allc = np.concatenate([acc_c] + [np.ones(len(b), np.int64) for b in buf])
                acc_k, inv = np.unique(allk, return_inverse=True)
                acc_c = np.bincount(inv, weights=allc).astype(np.int64)
                buf = []
    return acc_k // n, acc_k % n, acc_c / R


# ------------------------------------------------------------------------------------------ matching
def best_jaccard(members: np.ndarray, labels: np.ndarray, sizes: np.ndarray | None = None) -> float:
    """Best-match Jaccard of a member set against a partition (labels over the same nodes; -1 = unassigned)."""
    lab = labels[members]
    lab = lab[lab >= 0]
    if not len(lab):
        return 0.0
    sizes = np.bincount(labels[labels >= 0]) if sizes is None else sizes
    u, cnt = np.unique(lab, return_counts=True)
    return float(np.max(cnt / (len(members) + sizes[u] - cnt)))


def _restrict_labels(labels: np.ndarray, nodes_mask: np.ndarray | None) -> np.ndarray:
    if nodes_mask is None:
        return labels
    out = labels.copy()
    out[~nodes_mask] = -1
    return out


@dataclass
class Level:
    """Per-gamma consensus state (stage-1 memberships, stage-2 frequencies) shared by all re-clusterings."""
    gamma: float
    runs: np.ndarray               # n_runs x n stage-1 memberships
    fr: np.ndarray
    fc: np.ndarray
    ff: np.ndarray


def build_levels(graphs: list[Graph], cfg: ClusterConfig) -> list[Level]:
    out = []
    for gm in cfg.gammas:
        M = stage1(graphs, gm, cfg)
        r, c, f = coclustering(M)
        out.append(Level(gm, M, r, c, f))
    return out


def stage2(level: Level, n: int, nodes: np.ndarray | None, cfg: ClusterConfig) -> list[np.ndarray]:
    """Stage-2 clusters (size >= min_size, sorted member arrays) on the frequency graph induced on `nodes`."""
    r, c, f = level.fr, level.fc, level.ff
    if nodes is not None:                                   # induced subgraph, relabelled (small Leiden problem)
        nodes = np.sort(np.asarray(nodes))
        mask = np.zeros(n, bool)
        mask[nodes] = True
        k = mask[r] & mask[c]
        loc = np.full(n, -1)
        loc[nodes] = np.arange(len(nodes))
        sub = leiden(len(nodes), loc[r[k]], loc[c[k]], f[k], cfg.gamma2, cfg.stage2_seed, cfg.n_iterations)
        memb = np.full(n, -1)
        memb[nodes] = sub
    else:
        memb = leiden(n, r, c, f, cfg.gamma2, cfg.stage2_seed, cfg.n_iterations)
    o = np.argsort(memb, kind="stable")
    ms = memb[o]
    bounds = np.flatnonzero(np.r_[True, ms[1:] != ms[:-1], True])
    cl = [np.sort(o[a:b]) for a, b in zip(bounds[:-1], bounds[1:]) if ms[a] >= 0 and b - a >= cfg.min_size]
    cl.sort(key=lambda x: (-len(x), int(x[0])))
    return cl


def containment(members: np.ndarray, labels: np.ndarray) -> float:
    """Largest fraction of the members that one cluster of the partition holds."""
    lab = labels[members]
    lab = lab[lab >= 0]
    if not len(lab):
        return 0.0
    return float(np.unique(lab, return_counts=True)[1].max() / len(members))


def support(members, level: Level, nodes_mask, cfg: ClusterConfig) -> tuple[float, float]:
    """(support = fraction of stage-1 runs in which the cluster recurs [support_mode], mean best-match Jaccard
    over runs = stability)."""
    js, cs = [], []
    for lab in level.runs:
        lab = _restrict_labels(lab, nodes_mask)
        js.append(best_jaccard(members, lab))
        if cfg.support_mode == "containment":
            cs.append(containment(members, lab))
    js = np.array(js)
    rec = np.array(cs) if cfg.support_mode == "containment" else js
    if cfg.support_mode not in ("jaccard", "containment"):
        raise ValueError(cfg.support_mode)
    return float((rec >= cfg.jaccard_recur).mean()), float(js.mean())


def persistent(levels: list[Level], n: int, nodes, j_indices, cfg: ClusterConfig) -> list[dict]:
    """Persistent clusters over adjacent gamma indices j_indices (ascending, contiguous) on node subset `nodes`
    (None = all nodes)."""
    j_indices = list(j_indices)
    if not j_indices:
        return []
    nodes_mask = None
    if nodes is not None:
        nodes_mask = np.zeros(n, bool)
        nodes_mask[nodes] = True
    P = {j: stage2(levels[j], n, nodes, cfg) for j in j_indices}
    sup = {j: [support(c, levels[j], nodes_mask, cfg) for c in P[j]] for j in j_indices}
    labels = {}
    for j in j_indices:
        lab = np.full(n, -1)
        for i, c in enumerate(P[j]):
            lab[c] = i
        labels[j] = lab
    nxt = {}
    for a, b in zip(j_indices[:-1], j_indices[1:]):
        sizes = np.array([len(c) for c in P[b]]) if P[b] else np.zeros(0, int)
        for i, c in enumerate(P[a]):
            lab = labels[b][c]
            lab = lab[lab >= 0]
            if not len(lab):
                continue
            u, cnt = np.unique(lab, return_counts=True)
            jac = cnt / (len(c) + sizes[u] - cnt)
            t = int(np.argmax(jac))
            if jac[t] >= cfg.jaccard_recur:
                nxt[(a, i)] = (b, int(u[t]))
    has_prev = set(nxt.values())
    out = []
    for j in j_indices:
        for i in range(len(P[j])):
            if (j, i) in has_prev:
                continue
            chain = [(j, i)]
            while chain[-1] in nxt:
                chain.append(nxt[chain[-1]])
            # longest run of adjacent members with support >= support_min
            best, cur = [], []
            for node in chain:
                if sup[node[0]][node[1]][0] >= cfg.support_min:
                    cur.append(node)
                    if len(cur) > len(best):
                        best = list(cur)
                else:
                    cur = []
            if len(best) < cfg.min_adjacent:
                continue
            rep = max(best, key=lambda x: (sup[x[0]][x[1]][0], -x[0]))
            out.append(dict(members=P[rep[0]][rep[1]], j=rep[0], gamma=levels[rep[0]].gamma,
                            j_start=best[0][0], j_end=best[-1][0], n_adjacent=len(best),
                            support=sup[rep[0]][rep[1]][0], stability=sup[rep[0]][rep[1]][1]))
    out.sort(key=lambda d: (d["j_start"], -d["support"], -len(d["members"]), int(d["members"][0])))
    return out


# ------------------------------------------------------------------------------------------ hierarchy
@dataclass
class ClusterMap:
    clusters: list = field(default_factory=list)      # dicts: id, level, parent, members (local idx), gamma, ...
    counts: dict = field(default_factory=dict)

    def level(self, lev: int) -> list:
        return [c for c in self.clusters if c["level"] == lev]

    def labels(self, n: int, lev: int = 1) -> np.ndarray:
        lab = np.full(n, -1)
        for c in self.level(lev):
            lab[c["members"]] = c["id"]
        return lab


def cohesion(members, ref_scores: np.ndarray | None) -> float:
    if ref_scores is None:
        return np.nan
    sub = ref_scores[np.ix_(members, members)].astype(np.float64)
    iu = np.triu_indices(len(members), 1)
    v = sub[iu]
    return float(np.nanmean(v)) if np.isfinite(v).any() else -np.inf


def _place(cands, levels, n, cfg, ref_scores, z0, assigned, counts, depth, within=None, delta=None):
    """Greedy one-per-gene placement with re-clustering of rule violators. Returns accepted cluster dicts."""
    acc = []
    queue = list(cands)
    while queue:
        d = queue.pop(0)
        m = d["members"]
        if within is not None and len(m) >= len(within) and best_jaccard_sets(m, within) >= cfg.jaccard_recur:
            continue                                            # not a proper sub-cluster
        free = m[~assigned[m]]
        if len(free) < cfg.min_size or len(free) < cfg.trim_keep * len(m):
            counts["overlap_dropped"] += 1
            continue
        coh = cohesion(free, ref_scores)
        too_big = len(free) > cfg.max_size
        not_coh = delta is not None and ref_scores is not None and coh < z0 + delta
        if too_big or not_coh:
            counts["oversize" if too_big else "cohesion_fail"] += 1
            dd = d.get("depth", depth)
            if dd < cfg.max_depth and d["j"] + 1 < len(levels):
                sub = persistent(levels, n, free, range(d["j"] + 1, len(levels)), cfg)
                sub = [s for s in sub if best_jaccard_sets(s["members"], free) < cfg.jaccard_recur or len(s["members"]) < len(free)]
                counts["reclustered"] += 1
                queue = _merge_sorted(queue, [dict(s, depth=dd + 1) for s in sub])
            continue
        assigned[free] = True
        acc.append(dict(d, members=free, size=len(free), cohesion=coh, trimmed=len(m) - len(free)))
    return acc


def _merge_sorted(queue, new):
    allq = queue + new
    allq.sort(key=lambda d: (d["j_start"], -d["support"], -len(d["members"]), int(d["members"][0])))
    return allq


def best_jaccard_sets(a, b) -> float:
    inter = len(np.intersect1d(a, b, assume_unique=True))
    return inter / (len(a) + len(b) - inter)


def _persist_key(d):
    """'nested' placement order: most persistent first (adjacent gamma steps, then support), then larger, then
    coarser, then lowest member index (deterministic)."""
    return (-d["n_adjacent"], -d["support"], -len(d["members"]), d["j_start"], int(d["members"][0]))


def _place_nested(cands, levels, n, cfg, ref_scores, z0, counts, delta1, delta2, within=None):
    """FZ-C decision D2 (results/part2/construction/DECISIONS.md). Candidates (persistent chains) are visited most
    persistent first. Against the clusters already accepted (one partition per level):
      * duplicate (Jaccard >= jaccard_recur with an accepted cluster)             -> dropped ('duplicate');
      * nested child (>= nest_contain of it inside one accepted level-1 cluster)   -> sub-module candidate of that
        cluster (restricted to the parent), placed at level 2 after level 1 is complete;
      * parent (it holds >= nest_contain of one or more accepted level-1 clusters)  -> becomes the level-1 cluster
        (members = candidate U those clusters), the contained clusters become its level-2 sub-modules, provided
        the merged parent passes size and level-1 cohesion (else the candidate is dropped, 'parent_rule_fail');
      * partial overlap with other accepted clusters: if the overlapping genes are > overlap_conflict of the
        candidate it is dropped ('overlap_conflict', resolved to the more persistent cluster), otherwise they are
        trimmed from it.
    Rule violators (size > max_size, cohesion < z0 + delta) are re-clustered on the same frequencies as before.
    Returns (level-1 dicts, {level-1 index: [level-2 candidate dicts]})."""
    acc, kids = [], {}
    owner = np.full(n, -1)                       # level-1 index owning each gene
    queue = sorted(cands, key=_persist_key)
    while queue:
        d = queue.pop(0)
        m = np.asarray(d["members"])
        if within is not None:
            m = np.intersect1d(m, within)
        if len(m) < cfg.min_size:
            counts["overlap_dropped"] += 1
            continue
        own = owner[m]
        placed = np.unique(own[own >= 0])
        dup = parent_of = None
        children, partial = [], []
        for a in placed.tolist():
            am = acc[a]["members"]
            inter = int((own == a).sum())
            if inter / (len(m) + len(am) - inter) >= cfg.jaccard_recur:
                dup = a
                break
            if inter / len(m) >= cfg.nest_contain:
                parent_of = a
            elif inter / len(am) >= cfg.nest_contain:
                children.append(a)
            else:
                partial.append(a)
        if dup is not None:
            counts["duplicate"] += 1
            continue
        if parent_of is not None:
            kids.setdefault(parent_of, []).append(dict(d, members=np.intersect1d(m, acc[parent_of]["members"])))
            counts["nested_child"] += 1
            continue
        pmask = np.isin(own, partial) if partial else np.zeros(len(m), bool)
        if pmask.sum() > cfg.overlap_conflict * len(m):
            counts["overlap_conflict"] += 1
            continue
        free = m[~pmask]
        if partial:
            counts["trimmed"] += 1
        if children:
            merged = np.union1d(free, np.concatenate([acc[a]["members"] for a in children]))
            coh = cohesion(merged, ref_scores)
            if len(merged) > cfg.max_size or (delta1 is not None and ref_scores is not None and coh < z0 + delta1):
                counts["parent_rule_fail"] += 1
                continue
            # the merged parent replaces its children at level 1; they (and their own kids) become its sub-modules
            new_kids = []
            for a in children:
                new_kids.append(dict(acc[a], members=acc[a]["members"]))
                new_kids.extend(kids.pop(a, []))
                acc[a] = None
                owner[owner == a] = -1
            idx = len(acc)
            acc.append(dict(d, members=merged, size=len(merged), cohesion=coh, trimmed=len(m) - len(free)))
            owner[merged] = idx
            kids[idx] = new_kids
            counts["parent_absorbed"] += len(children)
            continue
        if len(free) < cfg.min_size:
            counts["overlap_dropped"] += 1
            continue
        coh = cohesion(free, ref_scores)
        too_big = len(free) > cfg.max_size
        not_coh = delta1 is not None and ref_scores is not None and coh < z0 + delta1
        if too_big or not_coh:
            counts["oversize" if too_big else "cohesion_fail"] += 1
            dd = d.get("depth", 0)
            if dd < cfg.max_depth and d["j"] + 1 < len(levels):
                sub = persistent(levels, n, free, range(d["j"] + 1, len(levels)), cfg)
                sub = [s for s in sub if best_jaccard_sets(s["members"], free) < cfg.jaccard_recur or len(s["members"]) < len(free)]
                counts["reclustered"] += 1
                queue = sorted(queue + [dict(s, depth=dd + 1) for s in sub], key=_persist_key)
            continue
        idx = len(acc)
        acc.append(dict(d, members=free, size=len(free), cohesion=coh, trimmed=len(m) - len(free)))
        owner[free] = idx
    out, out_kids = [], {}
    for i, a in enumerate(acc):
        if a is None:
            continue
        out_kids[len(out)] = kids.get(i, [])
        out.append(a)
    return out, out_kids


def _build_map_nested(graphs, cfg, ref_scores, z0, levels, sublevels):
    n = graphs[0].n
    counts = dict(oversize=0, cohesion_fail=0, reclustered=0, overlap_dropped=0, duplicate=0, nested_child=0,
                  overlap_conflict=0, trimmed=0, parent_rule_fail=0, parent_absorbed=0)
    top = persistent(levels, n, None, range(len(levels)), cfg)
    counts["persistent_level1_candidates"] = len(top)
    l1, nested = _place_nested(top, levels, n, cfg, ref_scores, z0, counts, cfg.cohesion_delta_level1, None)
    # deterministic level-1 order: coarse first as before
    order = sorted(range(len(l1)), key=lambda i: (l1[i]["j_start"], -len(l1[i]["members"]), int(l1[i]["members"][0])))
    cm = ClusterMap(counts=counts)
    for new_id, i in enumerate(order):
        cm.clusters.append(dict(l1[i], id=new_id, level=1, parent=-1))
    if sublevels:
        nid = len(l1)
        sc = dict(oversize=0, cohesion_fail=0, reclustered=0, overlap_dropped=0, duplicate=0, nested_child=0,
                  overlap_conflict=0, trimmed=0, parent_rule_fail=0, parent_absorbed=0, from_nesting=0,
                  from_reclustering=0)
        for new_id, i in enumerate(order):
            p = l1[i]
            pm = p["members"]
            if len(pm) < 2 * cfg.min_size:
                continue
            cands = [dict(c, origin="nesting") for c in nested.get(i, [])]
            if p["j"] + 1 < len(levels):
                cands += [dict(c, origin="reclustering")
                          for c in persistent(levels, n, pm, range(p["j"] + 1, len(levels)), cfg)]
            cands = [c for c in cands if len(c["members"]) < len(pm)
                     and best_jaccard_sets(np.intersect1d(c["members"], pm), pm) < cfg.jaccard_recur]
            # level 2: one partition per parent; no deeper nesting (depth = max), sub-module cohesion z0 + delta2
            sub, _ = _place_nested([dict(c, depth=cfg.max_depth) for c in cands], levels, n, cfg, ref_scores, z0,
                                   sc, cfg.cohesion_delta, None, within=pm)
            for d in sub:
                sc["from_" + d.get("origin", "reclustering")] += 1
                cm.clusters.append(dict(d, id=nid, level=2, parent=new_id))
                nid += 1
        counts["level2"] = sc
    counts["n_level1"] = len(cm.level(1))
    counts["n_level2"] = len(cm.level(2))
    return cm


def build_map(graphs: list[Graph], cfg: ClusterConfig = ClusterConfig(), ref_scores: np.ndarray | None = None,
              z0: float = 0.0, levels: list[Level] | None = None, sublevels: bool = True) -> ClusterMap:
    """Full consensus + persistence + top-down nesting. ref_scores: S x S reference block (full construction set)
    for the cohesion rule; z0: the graph threshold (cohesion >= z0 + cohesion_delta)."""
    n = graphs[0].n
    levels = build_levels(graphs, cfg) if levels is None else levels
    if cfg.placement == "nested":
        return _build_map_nested(graphs, cfg, ref_scores, z0, levels, sublevels)
    if cfg.placement != "greedy":
        raise ValueError(cfg.placement)
    counts = dict(oversize=0, cohesion_fail=0, reclustered=0, overlap_dropped=0)
    top = persistent(levels, n, None, range(len(levels)), cfg)
    counts["persistent_level1_candidates"] = len(top)
    assigned = np.zeros(n, bool)
    l1 = _place(top, levels, n, cfg, ref_scores, z0, assigned, counts, depth=0, delta=cfg.cohesion_delta_level1)
    cm = ClusterMap(counts=counts)
    l1 = [dict(d, id=i, level=1, parent=-1) for i, d in enumerate(l1)]
    cm.clusters.extend(l1)
    if sublevels:
        nid = len(l1)
        sub_counts = dict(oversize=0, cohesion_fail=0, reclustered=0, overlap_dropped=0)
        for p in l1:
            if len(p["members"]) < 2 * cfg.min_size or p["j"] + 1 >= len(levels):
                continue
            cands = persistent(levels, n, p["members"], range(p["j"] + 1, len(levels)), cfg)
            cands = [c for c in cands if len(c["members"]) < len(p["members"])
                     and best_jaccard_sets(c["members"], p["members"]) < cfg.jaccard_recur]
            asg = np.zeros(n, bool)
            for d in _place(cands, levels, n, cfg, ref_scores, z0, asg, sub_counts, depth=cfg.max_depth,
                            within=p["members"], delta=cfg.cohesion_delta):
                cm.clusters.append(dict(d, id=nid, level=2, parent=p["id"]))
                nid += 1
        counts["level2"] = sub_counts
    counts["n_level1"] = len(cm.level(1))
    counts["n_level2"] = len(cm.level(2))
    return cm
