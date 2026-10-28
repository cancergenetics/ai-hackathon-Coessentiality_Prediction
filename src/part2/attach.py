"""Attachment of genes outside S to frozen clusters (DESIGN v3.1 §2.3; FZ-C decision D3).

Gene g (not in S) is attached to cluster C if
  (i)  >= 50% of C's members have score >= z0 to g  (x_gC >= ceil(|C| / 2)), and
  (ii) a size-matched null gives BH FDR < 0.01 over all gene x cluster tests.
Null (round-2 algorithms SHOULD 9, closed form): for a random member set of size |C| drawn from S, the number of
members with score >= z0 to g is hypergeometric(N_g, K_g, |C|), with N_g = S genes with a finite score to g and
K_g = those with score >= z0; p_gC = P(X >= x_gC). BH runs over every (gene, cluster) test (pairs with x = 0 have
p = 1 and stay in the denominator). Attached genes are not members; they can define H1 events (DESIGN §2.6).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import hypergeom


def bh_partial(p: np.ndarray, m_total: int) -> np.ndarray:
    """BH q-values for the p-values given, treating m_total - len(p) further tests as p = 1 (ranked last)."""
    p = np.asarray(p, float)
    if not len(p):
        return p
    o = np.argsort(p, kind="stable")
    ranks = np.arange(1, len(p) + 1)
    q = p[o] * m_total / ranks
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[o] = np.minimum(q, 1.0)
    return out


def attach(rect: np.ndarray, clusters: list[np.ndarray], z0: float, frac: float = 0.5, fdr: float = 0.01,
           row_ids=None, cluster_ids=None) -> tuple[pd.DataFrame, dict]:
    """rect: genes (outside S) x S scores (NaN allowed); clusters: member arrays of S-local indices.
    Returns (table of attached (gene, cluster) pairs with x, |C|, p, q; summary counts)."""
    R = np.asarray(rect, np.float32)
    fin = np.isfinite(R)
    hit = fin & (R >= z0)
    N = fin.sum(1)
    K = hit.sum(1)
    row_ids = np.arange(len(R)) if row_ids is None else np.asarray(row_ids)
    cluster_ids = np.arange(len(clusters)) if cluster_ids is None else np.asarray(cluster_ids)
    tested_rows = N > 0
    sizes = np.array([len(c) for c in clusters])
    M = np.zeros((R.shape[1], len(clusters)), np.float32)
    for j, c in enumerate(clusters):
        M[c, j] = 1.0
    X = (hit.astype(np.float32) @ M).round().astype(np.int64)               # genes x clusters member hits
    m_total = int(tested_rows.sum()) * len(clusters)
    gi, cj = np.nonzero((X > 0) & tested_rows[:, None])
    x = X[gi, cj]
    p = hypergeom.sf(x - 1, N[gi], K[gi], sizes[cj])
    q = bh_partial(p, m_total)
    need = np.ceil(frac * sizes[cj]).astype(int)
    ok = (x >= need) & (q < fdr)
    df = pd.DataFrame(dict(gene=row_ids[gi[ok]], cluster=cluster_ids[cj[ok]], x=x[ok], size=sizes[cj[ok]],
                           frac_members=x[ok] / sizes[cj[ok]], K_g=K[gi[ok]], N_g=N[gi[ok]], p=p[ok], q=q[ok]))
    df = df.sort_values(["cluster", "q", "gene"]).reset_index(drop=True)
    summ = dict(n_genes_tested=int(tested_rows.sum()), n_clusters=len(clusters), n_tests=m_total,
                n_majority=int((x >= need).sum()), n_attached_pairs=int(ok.sum()),
                n_attached_genes=int(df.gene.nunique()), expected_false_at_fdr=fdr * int(ok.sum()))
    return df, summ
