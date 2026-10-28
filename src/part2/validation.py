"""Validation of the frozen FZ-C maps (DESIGN v3.1 §2.3 tiers + bridge). Validation inputs only: BioPlex 3.0
(293T and HCT116 networks, Entrez-mapped) and CCLE proteomics restricted to A models. Never a construction or
selection input: callers must pass ``require_frozen`` (the map manifests are committed in git and unchanged).

T1  within-cluster BioPlex edges among BioPlex-observed members (>= 3 observed, >= 2 edges) above a node-restricted,
    degree-preserving null: every observed member replaced by a random gene of the same BioPlex-degree bin
    (log2) and A-essentiality class from the universe (BioPlex nodes in S), 10^4 draws; p = (1 + #null >= obs) /
    (1 + B); BH over clusters (per level) q < 0.05.
T2  (not T1) mean pairwise protein co-abundance (pairwise-complete Pearson over A models with proteomics) among
    measured members (>= 3) above a size-matched null (random measured S genes), 10^4 draws; BH (per level)
    q < 0.05.
T3-neg tested (BioPlex >= 3 observed members or >= 3 measured members) but not T1/T2; T3-nc otherwise.
Bridge: size-standardised within-cluster BioPlex precision (pairs with both proteins observed; bins 3-4, 5-8, 9-16,
    17-60 by cluster size; bin weights = pooled zavg + GLS cluster-size distribution), paired difference with a
    protein-node bootstrap (B = 2000, seed 5000, clusters fixed, pair weight = product of node multiplicities).
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from part2 import paths

BIOPLEX_FILES = (paths.DATASETS / "BioPlex" / "BioPlex_293T_Network_10K_Dec_2019.tsv",
                 paths.DATASETS / "BioPlex" / "BioPlex_HCT116_Network_5.5K_Dec_2019.tsv")
SIZE_BINS = ((3, 4), (5, 8), (9, 16), (17, 60))


class NotFrozenError(RuntimeError):
    pass


def require_frozen(map_dirs) -> dict:
    """Refuse to run unless every map directory's MANIFEST.sha256 is committed at HEAD, the working-tree files match
    it, and the manifest itself is unchanged since HEAD (FZ-C: validation inputs only after the map commit)."""
    import hashlib
    out = {}
    for d in map_dirs:
        d = Path(d)
        man = d / "MANIFEST.sha256"
        rel = man.relative_to(paths.PROJECT)
        r = subprocess.run(["git", "-C", str(paths.PROJECT), "show", f"HEAD:{rel}"], capture_output=True, text=True)
        if r.returncode != 0 or r.stdout != man.read_text():
            raise NotFrozenError(f"{rel} is not committed at HEAD (or differs): commit the FZ-C maps first")
        for line in r.stdout.splitlines():
            h, name = line.split("  ", 1)
            if hashlib.sha256((d / name).read_bytes()).hexdigest() != h:
                raise NotFrozenError(f"{d.name}/{name} differs from the committed manifest")
        log = subprocess.run(["git", "-C", str(paths.PROJECT), "log", "-1", "--format=%H", "--", str(rel)],
                             capture_output=True, text=True).stdout.strip()
        out[d.name] = log
    return out


# ------------------------------------------------------------------------------------------ BioPlex
def load_bioplex(files=BIOPLEX_FILES) -> tuple[set, pd.Series]:
    """Union of undirected edges (Entrez pairs, min < max) of the given networks and node degree in the union."""
    keys = set()
    for f in files:
        df = pd.read_csv(f, sep="\t", usecols=["GeneA", "GeneB"])
        a = pd.to_numeric(df.GeneA, errors="coerce")
        b = pd.to_numeric(df.GeneB, errors="coerce")
        ok = a.notna() & b.notna() & (a != b)
        a, b = a[ok].astype(np.int64).values, b[ok].astype(np.int64).values
        keys.update(zip(np.minimum(a, b).tolist(), np.maximum(a, b).tolist()))
    e = np.array(sorted(keys), dtype=np.int64)
    deg = pd.Series(np.r_[e[:, 0], e[:, 1]]).value_counts().sort_index()
    return keys, deg


class Network:
    """BioPlex restricted to a gene universe (S, local indices): dense adjacency over observed S genes."""

    def __init__(self, entrez: np.ndarray, edges: set, deg: pd.Series):
        self.entrez = np.asarray(entrez, np.int64)
        self.observed = np.isin(self.entrez, deg.index.values)
        self.obs_idx = np.flatnonzero(self.observed)
        self.pos = np.full(len(self.entrez), -1)
        self.pos[self.obs_idx] = np.arange(len(self.obs_idx))
        n = len(self.obs_idx)
        A = np.zeros((n, n), bool)
        ez = self.entrez[self.obs_idx]
        loc = pd.Series(np.arange(n), index=ez)
        loc = loc[~loc.index.duplicated()]
        for a, b in edges:
            if a in loc.index and b in loc.index:
                i, j = loc[a], loc[b]
                A[i, j] = A[j, i] = True
        self.adj = A
        self.degree = deg.reindex(ez).fillna(0).astype(int).values     # degree in the full union network
        self.n_edges_universe = int(A.sum() // 2)

    def edges_within(self, members) -> tuple[int, int]:
        """(number of observed members, number of BioPlex edges among them)."""
        p = self.pos[np.asarray(members)]
        p = p[p >= 0]
        if len(p) < 2:
            return len(p), 0
        return len(p), int(np.triu(self.adj[np.ix_(p, p)], 1).sum())


def strata(degree: np.ndarray, ess: np.ndarray, min_n: int = 10) -> np.ndarray:
    """Stratum label per observed gene: (log2 degree bin, essentiality class); sparse degree bins are merged into
    the next lower bin within a class until every stratum has >= min_n genes (then classes are merged if needed)."""
    b = np.floor(np.log2(np.maximum(degree, 1))).astype(int)
    out = np.empty(len(b), object)
    for cl in np.unique(ess):
        m = ess == cl
        bins = b[m].copy()
        while True:
            u, c = np.unique(bins, return_counts=True)
            small = u[c < min_n]
            if not len(small) or len(u) == 1:
                break
            s = small[-1]
            lower = u[u < s]
            tgt = lower[-1] if len(lower) else u[u > s][0]
            bins[bins == s] = tgt
        out[m] = [f"{cl}:{x}" for x in bins]
    u, c = np.unique(out, return_counts=True)
    for s in u[c < min_n]:                                   # tiny class: pool with all genes of its degree bin
        out[out == s] = "pooled:" + s.split(":")[1]
    return out


def t1_test(net: Network, clusters: list[np.ndarray], ess: np.ndarray, B: int = 10_000, seed: int = 4100,
            min_obs: int = 3) -> pd.DataFrame:
    """Per cluster: observed members, edges, null mean, p (degree x essentiality preserving node replacement)."""
    rng = np.random.default_rng(seed)
    st = strata(net.degree, ess[net.obs_idx])
    pools = {s: np.flatnonzero(st == s) for s in np.unique(st)}
    rows = []
    for ci, m in enumerate(clusters):
        p = net.pos[np.asarray(m)]
        p = np.sort(p[p >= 0])
        n_obs = len(p)
        if n_obs < min_obs:
            rows.append(dict(i=ci, n_obs=n_obs, edges=np.nan, n_pairs=n_obs * (n_obs - 1) // 2, null_mean=np.nan,
                             p=np.nan))
            continue
        obs = int(np.triu(net.adj[np.ix_(p, p)], 1).sum())
        draw = np.empty((B, n_obs), np.int64)
        for k, g in enumerate(p):
            pool = pools[st[g]]
            draw[:, k] = pool[rng.integers(0, len(pool), B)]
        for _ in range(20):                                  # resample within-draw duplicates
            s = np.sort(draw, 1)
            dup = (s[:, 1:] == s[:, :-1]).any(1)
            if not dup.any():
                break
            for k, g in enumerate(p):
                pool = pools[st[g]]
                draw[dup, k] = pool[rng.integers(0, len(pool), dup.sum())]
        ia, ib = np.triu_indices(n_obs, 1)
        null = np.zeros(B, np.int64)
        for s0 in range(0, len(ia), 4096):
            null += net.adj[draw[:, ia[s0:s0 + 4096]], draw[:, ib[s0:s0 + 4096]]].sum(1)
        rows.append(dict(i=ci, n_obs=n_obs, edges=obs, n_pairs=len(ia), null_mean=float(null.mean()),
                         p=(1 + int((null >= obs).sum())) / (1 + B)))
    return pd.DataFrame(rows)


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    q = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return q
    pv = p[ok]
    o = np.argsort(pv, kind="stable")
    r = pv[o] * len(pv) / np.arange(1, len(pv) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    qq = np.empty(len(pv))
    qq[o] = np.minimum(r, 1)
    q[ok] = qq
    return q


# ------------------------------------------------------------------------------------------ proteomics
def t2_test(R: np.ndarray, measured: np.ndarray, clusters: list[np.ndarray], B: int = 10_000, seed: int = 4200,
            min_obs: int = 3) -> pd.DataFrame:
    """R: protein correlation matrix over the measured S genes (local order = np.flatnonzero(measured)).
    Statistic: mean pairwise r among measured members; null: random measured S genes of the same count."""
    rng = np.random.default_rng(seed)
    pos = np.full(len(measured), -1)
    mi = np.flatnonzero(measured)
    pos[mi] = np.arange(len(mi))
    rows = []
    for ci, m in enumerate(clusters):
        p = pos[np.asarray(m)]
        p = np.sort(p[p >= 0])
        n = len(p)
        if n < min_obs:
            rows.append(dict(i=ci, n_meas=n, mean_r=np.nan, null_mean=np.nan, p=np.nan))
            continue
        ia, ib = np.triu_indices(n, 1)
        obs = float(np.nanmean(R[p[ia], p[ib]]))
        draw = rng.integers(0, len(mi), (B, n))
        for _ in range(50):                                  # resample within-draw duplicates
            s = np.sort(draw, 1)
            dup = (s[:, 1:] == s[:, :-1]).any(1)
            if not dup.any():
                break
            draw[dup] = rng.integers(0, len(mi), (int(dup.sum()), n))
        null = np.nanmean(R[draw[:, ia], draw[:, ib]], axis=1)
        rows.append(dict(i=ci, n_meas=n, mean_r=obs, null_mean=float(np.nanmean(null)),
                         p=(1 + int((null >= obs).sum())) / (1 + B)))
    return pd.DataFrame(rows)


def assign_tiers(t1: pd.DataFrame, t2: pd.DataFrame, alpha: float = 0.05, min_edges: int = 2) -> pd.DataFrame:
    df = t1.rename(columns={"p": "p_bioplex", "null_mean": "null_mean_edges"}).merge(
        t2.rename(columns={"p": "p_prot", "null_mean": "null_mean_r"}), on="i")
    df["q_bioplex"] = bh(df.p_bioplex.values)
    df["q_prot"] = bh(df.p_prot.values)
    t1ok = (df.q_bioplex < alpha) & (df.edges >= min_edges)
    t2ok = ~t1ok & (df.q_prot < alpha)
    tested = (df.n_obs >= 3) | (df.n_meas >= 3)
    df["tier"] = np.where(t1ok, "T1", np.where(t2ok, "T2", np.where(tested, "T3-neg", "T3-nc")))
    return df


# ------------------------------------------------------------------------------------------ bridge
def size_bin(n: int) -> int:
    for b, (lo, hi) in enumerate(SIZE_BINS):
        if lo <= n <= hi:
            return b
    return -1


def pair_table(net: Network, clusters: list[np.ndarray]) -> pd.DataFrame:
    """Within-cluster pairs with both proteins BioPlex-observed: node positions (observed order), edge, bin."""
    out = []
    for ci, m in enumerate(clusters):
        b = size_bin(len(m))
        p = net.pos[np.asarray(m)]
        p = np.sort(p[p >= 0])
        if len(p) < 2:
            continue
        ia, ib = np.triu_indices(len(p), 1)
        out.append(pd.DataFrame(dict(cluster=ci, bin=b, a=p[ia], b=p[ib], edge=net.adj[p[ia], p[ib]])))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame(columns=["cluster", "bin", "a", "b", "edge"])


def bin_weights(sizes_a, sizes_b) -> np.ndarray:
    allb = np.array([size_bin(s) for s in list(sizes_a) + list(sizes_b)])
    return np.array([(allb == b).mean() for b in range(len(SIZE_BINS))])


def standardised(pairs: pd.DataFrame, w: np.ndarray, value: str = "edge", mult: np.ndarray | None = None):
    """Size-standardised precision (pair-pooled within bin, bins weighted by w; bins without pairs in this map are
    dropped and the remaining weights renormalised - the caller passes weights over bins defined in both maps)."""
    v = pairs[value].values.astype(float)
    ww = np.ones(len(v)) if mult is None else mult[pairs.a.values] * mult[pairs.b.values]
    per = np.full(len(SIZE_BINS), np.nan)
    for b in range(len(SIZE_BINS)):
        k = pairs.bin.values == b
        s = ww[k].sum()
        if s > 0:
            per[b] = (ww[k] * v[k]).sum() / s
    ok = np.isfinite(per) & (w > 0)
    return float((per[ok] * w[ok]).sum() / w[ok].sum()), per


def bridge_bootstrap(pz: pd.DataFrame, pg: pd.DataFrame, n_nodes: int, w: np.ndarray, B: int = 2000,
                     seed: int = 5000, value: str = "edge") -> dict:
    """Paired protein-node bootstrap of the difference (zavg - GLS) of size-standardised precision. Nodes = all
    observed proteins of the universe positions used by either map (multinomial resample of n_nodes)."""
    used = np.unique(np.r_[pz.a.values, pz.b.values, pg.a.values, pg.b.values])
    rng = np.random.default_rng(seed)
    dz, per_z = standardised(pz, w, value)
    dg, per_g = standardised(pg, w, value)
    diffs = np.empty(B)
    vz, vg = np.empty(B), np.empty(B)
    for t in range(B):
        mult = np.zeros(n_nodes)
        mult[used] = rng.multinomial(len(used), np.full(len(used), 1 / len(used)))
        vz[t] = standardised(pz, w, value, mult)[0]
        vg[t] = standardised(pg, w, value, mult)[0]
        diffs[t] = vz[t] - vg[t]
    lo, hi = np.nanpercentile(diffs, [2.5, 97.5])
    return dict(zavg=dz, gls=dg, diff=dz - dg, ci95=[float(lo), float(hi)], excludes_0=bool(lo > 0 or hi < 0),
                zavg_ci95=np.nanpercentile(vz, [2.5, 97.5]).tolist(), gls_ci95=np.nanpercentile(vg, [2.5, 97.5]).tolist(),
                per_bin_zavg=per_z.tolist(), per_bin_gls=per_g.tolist(), bin_weights=w.tolist(), n_nodes=int(len(used)),
                n_pairs_zavg=int(len(pz)), n_pairs_gls=int(len(pg)), B=B, seed=seed)


def match_coverage(tab_a: pd.DataFrame, tab_b: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Decision D5: level-1 gene coverage matching. The map with more clustered genes drops its least stable
    level-1 clusters (ascending stability, then ascending cohesion, then descending cluster id) until its coverage
    is <= the other map's."""
    ca, cb = int(tab_a["size"].sum()), int(tab_b["size"].sum())
    info = dict(genes_before=dict(a=ca, b=cb))

    def trim(t, target):
        o = t.sort_values(["stability", "cohesion", "cluster_id"], ascending=[True, True, False])
        keep = np.ones(len(o), bool)
        tot = int(o["size"].sum())
        for i in range(len(o)):
            if tot <= target:
                break
            keep[i] = False
            tot -= int(o["size"].iloc[i])
        return o[keep].sort_values("cluster_id"), int((~keep).sum())

    if ca > cb:
        tab_a, n = trim(tab_a, cb)
        info["trimmed"] = dict(map="a", n_clusters_dropped=n)
    elif cb > ca:
        tab_b, n = trim(tab_b, ca)
        info["trimmed"] = dict(map="b", n_clusters_dropped=n)
    info["genes_after"] = dict(a=int(tab_a["size"].sum()), b=int(tab_b["size"].sum()))
    return tab_a, tab_b, info
