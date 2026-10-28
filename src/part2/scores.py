"""Lean locked score recompute for a model subset (no n x n matrix over all genes is ever stored).

Per run (model set) the scorer writes the **S x S score block** over a fixed gene set S (the dependency-bearing
filter set, computed once on the construction set A; round-2 algorithms review MUST 2) as float32, plus, for
zavg_BEF, the per-gene top-K lists over *all* protein genes for the S rows (free by-product of the block scan).
Top-k within S, mutual edges and all-pairs cohesion are then exact.

Cache:  data_cache/part2/blocks/<method>_<model-set hash>_S<S hash>.npz + .json sidecar (model-set hash, model IDs,
        runtime, stage timings, peak RSS, code commit, null moments, ...).

zavg_BEF: ``CoessentialityModel(X, Config())`` from the locked Part 1 package, unchanged. Its genome-wide state is
recomputed on the given model subset by the locked procedure: same Config (k=20, tau=0.1, q=0.25, d=512,
max_nan=0.1, n_null=200,000, seed=0 -> same null-pair draw procedure and seeded ARPACK v0), null z-moments recomputed
on the subset's own data (never Part 1's stored moments) and recorded in the sidecar. Rows are scored with the
model's fixed gene blocks and the arithmetic of Part 1's ``score_all_pairs`` (z_block -> (B+E+F)/3), the upper
triangle (row gene index < column gene index, genes.tsv order) is the computed value and the lower triangle its
mirror, exactly as Part 1's stored matrix; on all models the block equals Part 1's zavg_BEF.npy restricted to S
(tests/test_part2_scores.py). Genes of S that are not background genes in the run (NaN fraction >= 10% in that
subset) get NaN rows/columns in that run only (``nan_non_bg``).

GLS: Part 1's dev-selected variant ``gls_signed`` (src/predict.py ``BaselineEngine._gls`` with src/baselines.py):
signed Wainberg GLS t; sample covariance = np.cov over all DepMap genes with no NaN in the subset
(``sample_covariance``); exact pairwise-complete marginalisation per missingness-pattern pair (whitening of
Sigma[S', S'], rank deficiency 'drop'); t = r sqrt(df / (1 - r^2)), r = cosine of whitened, intercept-projected
profiles. Only the batching of the matrix products differs from Part 1 (agreement ~1e-10, tested).
"""
from __future__ import annotations

import hashlib
import json
import logging
import platform
import resource
import subprocess
import sys
import time

import numpy as np

from part2 import paths

log = logging.getLogger(__name__)
METHODS = ("zavg_BEF", "gls", "pearson")


def peak_rss_gb() -> float:
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 1024 ** 3 if sys.platform == "darwin" else r / 1024 ** 2


def code_commit() -> dict:
    def g(*a):
        return subprocess.run(["git", "-C", str(paths.PROJECT), *a], capture_output=True, text=True).stdout.strip()
    return dict(commit=g("rev-parse", "HEAD"),
                dirty_part2=[l[3:] for l in g("status", "--porcelain", "--", "src/part2").splitlines()])


def mask_hash(mask) -> str:
    return hashlib.sha256(np.packbits(np.asarray(mask, bool)).tobytes()).hexdigest()[:12]


# ------------------------------------------------------------------------------------------ top-K helpers
def lexsort_rows(v, ids):
    """Per-row order: value descending, ties by id ascending."""
    o = np.argsort(ids, axis=1, kind="stable")
    o2 = np.argsort(-np.take_along_axis(v, o, 1), axis=1, kind="stable")
    return np.take_along_axis(o, o2, 1)


def select_topk(S: np.ndarray, K: int, col_ids: np.ndarray | None = None):
    """Row-wise top-K of S (NaN / -inf never selected). Returns (ids int32, vals float32) ordered by value desc,
    ties by column id asc; missing slots are (-1, NaN). col_ids maps S columns to ids (default 0..m-1)."""
    S = np.where(np.isnan(S), -np.inf, np.asarray(S, np.float64))
    m = S.shape[1]
    col_ids = np.arange(m) if col_ids is None else np.asarray(col_ids)
    k = min(K, m)
    if k < m:
        # include every column tied with the k-th value so tie-breaking by id is exact
        part = np.argpartition(-S, k - 1, axis=1)[:, :k]
        kth = np.take_along_axis(S, part, 1).min(1)
        extra = (S >= kth[:, None]).sum(1).max()
        if extra > k:
            part = np.argpartition(-S, extra - 1, axis=1)[:, :extra]
    else:
        part = np.tile(np.arange(m), (len(S), 1))
    pv = np.take_along_axis(S, part, 1)
    pid = col_ids[part]
    o = lexsort_rows(pv, pid)[:, :k]
    pv, pid = np.take_along_axis(pv, o, 1), np.take_along_axis(pid, o, 1)
    I = np.full((len(S), K), -1, np.int32)
    V = np.full((len(S), K), np.nan, np.float32)
    ok = np.isfinite(pv)
    I[:, :k] = np.where(ok, pid, -1)
    V[:, :k] = np.where(ok, pv, np.nan)
    return I, V


def mirror_upper(B: np.ndarray) -> np.ndarray:
    """Symmetrise in place: lower triangle := upper triangle (row < column), diagonal NaN (Part 1 convention)."""
    iu = np.triu_indices(len(B), 1)
    B[(iu[1], iu[0])] = B[iu]
    np.fill_diagonal(B, np.nan)
    return B


# ------------------------------------------------------------------------------------------ zavg_BEF
def zavg_model(X, config=None):
    """Locked CoessentialityModel with its genome-wide state computed (timed stages); null moments on X."""
    from part2.part1 import LOCKED_CONFIG, CoessentialityModel
    cfg = LOCKED_CONFIG if config is None else config
    t, tim = time.perf_counter(), {}
    m = CoessentialityModel(X, cfg)
    for name, fn in (("init", None), ("eig+profiles", m.profiles), ("f_row_stats", m.f_row_stats),
                     ("neighbours", m.neighbours), ("so_row_stats", m.so_row_stats), ("embedding_arpack", m.embedding)):
        if fn is not None:
            fn()
        tim[name] = round(time.perf_counter() - t, 2)
        t = time.perf_counter()
    for c in "BEF":
        m.null_moments(c)
    tim["null_moments"] = round(time.perf_counter() - t, 2)
    return m, tim


def zavg_rows(m, gidx, rows_wanted=None):
    """Iterate the model's fixed gene blocks: (positions into gidx of the block's wanted protein rows, zavg_BEF
    float64 rows restricted to the protein columns gidx). Arithmetic identical to Part 1 score_all_pairs."""
    pos = np.full(m.G, -1)
    pos[gidx] = np.arange(len(gidx))
    want = None if rows_wanted is None else np.isin(np.arange(len(gidx)), rows_wanted)
    for rows in m._blocks():
        loc = np.flatnonzero(pos[rows] >= 0)
        if want is not None:
            loc = loc[want[pos[rows[loc]]]]
        if not len(loc):
            continue
        Zb = m.z_block(rows)                                # fixed genome-wide block
        Z = {c: Zb[c][loc][:, gidx] for c in "BEF"}
        del Zb
        yield pos[rows[loc]], (Z["B"] + Z["E"] + Z["F"]) / 3.0


def zavg_block(X, gidx, S_mask, K_all=50, config=None, nan_non_bg=True, rect=False):
    """S x S zavg_BEF block (+ top-K over all protein genes for S rows). rect=True additionally returns, in
    info['_rect'] = (vals, pos), the (non-S protein genes) x S rectangle used by the attachment rule (DESIGN v3.1
    §2.3): entry (g, s) is the value Part 1's upper-triangle convention would store, i.e. computed in row min(g, s)
    (genes.tsv order). The S x S block is unchanged by rect (same rows, same arithmetic)."""
    m, tim = zavg_model(X, config)
    s_pos = np.flatnonzero(S_mask)                          # positions into gidx (ascending)
    ns = len(s_pos)
    loc_of = np.full(len(gidx), -1)
    loc_of[s_pos] = np.arange(ns)
    B = np.empty((ns, ns), np.float32)
    TI = np.full((ns, K_all), -1, np.int32)
    TV = np.full((ns, K_all), np.nan, np.float32)
    if rect:
        r_pos = np.flatnonzero(~np.asarray(S_mask, bool))  # non-S protein genes (positions into gidx)
        rloc = np.full(len(gidx), -1)
        rloc[r_pos] = np.arange(len(r_pos))
        R_row = np.full((len(r_pos), ns), np.nan, np.float32)   # computed in the non-S row
        R_col = np.full((ns, len(r_pos)), np.nan, np.float32)   # computed in the S row
    t = time.perf_counter()
    for ip, za in zavg_rows(m, gidx, None if rect else s_pos):
        z32 = za.astype(np.float32)                          # stored precision (as Part 1's matrix)
        ins = loc_of[ip] >= 0
        if ins.any():
            zs = z32[ins]
            B[loc_of[ip[ins]]] = zs[:, s_pos]
            if rect:
                R_col[loc_of[ip[ins]]] = zs[:, r_pos]
            if K_all:
                zz = zs.astype(np.float64)
                zz[np.arange(len(zz)), ip[ins]] = -np.inf
                TI[loc_of[ip[ins]]], TV[loc_of[ip[ins]]] = select_topk(zz, K_all)
        if rect and (~ins).any():
            R_row[rloc[ip[~ins]]] = z32[~ins][:, s_pos]
    tim["block_scan"] = round(time.perf_counter() - t, 2)
    mirror_upper(B)
    non_bg = ~m.bg[gidx[s_pos]]
    if nan_non_bg and non_bg.any():
        B[non_bg] = np.nan
        B[:, non_bg] = np.nan
    info = dict(null_moments={c: dict(mean=m._moments[c][0], sd=m._moments[c][1]) for c in "BEF"},
                n_bg=int(m.n_bg), knn_K=int(m.K), n_S_non_bg_nan=int(non_bg.sum()) if nan_non_bg else 0,
                config={k: getattr(m.cfg, k) for k in m.cfg.__dataclass_fields__})
    if rect:
        upper = r_pos[:, None] < s_pos[None, :]             # row gene index < column gene index
        Rv = np.where(upper, R_row, R_col.T).astype(np.float32)
        del R_row, R_col
        if nan_non_bg:
            Rv[~m.bg[gidx[r_pos]]] = np.nan
            Rv[:, non_bg] = np.nan
        info["_rect"] = (Rv, r_pos)
    return B, (TI, TV), tim, info


# ------------------------------------------------------------------------------------------ GLS
def gls_block(X, gidx, S_mask):
    """Signed GLS t (Part 1 gls_signed) over the gene set S (positions into gidx): symmetric float32 block."""
    from part2.part1 import p1_baselines as PB
    tim, t = {}, time.perf_counter()
    X = np.asarray(X)
    Sigma, ng = PB.sample_covariance(X)                    # all DepMap genes with no NaN, as Part 1
    tim["covariance"] = round(time.perf_counter() - t, 2)
    s_pos = np.flatnonzero(S_mask)
    XT = np.ascontiguousarray(np.asarray(X[:, gidx[s_pos]], dtype=np.float64).T)    # S genes x models
    valid = ~np.isnan(XT)
    keys = [valid[j].tobytes() for j in range(len(s_pos))]
    uniq = {k: i for i, k in enumerate(dict.fromkeys(keys))}
    pid = np.array([uniq[k] for k in keys])
    pmask = [np.frombuffer(k, dtype=bool) for k in uniq]
    members = [np.flatnonzero(pid == p) for p in range(len(pmask))]
    B = np.full((len(s_pos), len(s_pos)), np.nan, np.float32)
    white, n_white = {}, 0
    t = time.perf_counter()
    for P in range(len(pmask)):
        for Q in range(P, len(pmask)):
            Sm = pmask[P] & pmask[Q]
            key = Sm.tobytes()
            if key not in white:
                if Sm.sum() < 4:
                    white[key] = None
                else:
                    A, keep = PB.whitening(Sigma[np.ix_(Sm, Sm)], "drop")
                    w0 = A.sum(1)
                    white[key] = (A, np.flatnonzero(Sm)[keep], w0 / np.linalg.norm(w0))
                    n_white += 1
            if white[key] is None:
                continue
            A, Sk, w0 = white[key]
            df = A.shape[0] - 2
            if df < 1:
                continue

            def resid(g):
                W = XT[g][:, Sk] @ A.T
                R = W - np.outer(W @ w0, w0)                # project out the warped intercept
                with np.errstate(invalid="ignore", divide="ignore"):
                    return R / np.linalg.norm(R, axis=1, keepdims=True)
            gP, gQ = members[P], members[Q]
            RP = resid(gP)
            RQ = RP if P == Q else resid(gQ)
            for s in range(0, len(gP), 2048):
                r = np.clip(RP[s:s + 2048] @ RQ.T, -1.0, 1.0)
                with np.errstate(invalid="ignore", divide="ignore"):
                    T = (r * np.sqrt(df / (1.0 - r * r))).astype(np.float32)
                B[np.ix_(gP[s:s + 2048], gQ)] = T
                if P != Q:
                    B[np.ix_(gQ, gP[s:s + 2048])] = T.T
    mirror_upper(B)
    tim["pattern_pairs"] = round(time.perf_counter() - t, 2)
    info = dict(n_genes_for_covariance=ng, n_patterns=len(pmask), n_whitenings=n_white)
    return B, None, tim, info


def pairwise_pearson(Y: np.ndarray) -> np.ndarray:
    """Pairwise-complete Pearson r between the columns of Y (models x genes, NaN = missing); exact per pair
    (sums over the models where both genes are observed), equal to Part 1's pearson_signed."""
    Y = np.asarray(Y, np.float64)
    M = (~np.isnan(Y)).astype(np.float64)
    X0 = np.where(M > 0, Y, 0.0)
    n = M.T @ M
    Sx = X0.T @ M                                           # sum of x_i over models where i and j observed
    Sxx = (X0 * X0).T @ M
    Sxy = X0.T @ X0
    with np.errstate(invalid="ignore", divide="ignore"):
        num = n * Sxy - Sx * Sx.T
        den = np.sqrt((n * Sxx - Sx * Sx) * (n * Sxx.T - Sx.T * Sx.T))
        r = num / den
    r[n < 3] = np.nan
    return r


def pearson_block(X, gidx, S_mask):
    """Pearson (Part 1 baseline pearson_signed, pairwise complete) over S: symmetric float32 block (bridge
    secondary comparator)."""
    t = time.perf_counter()
    s_pos = np.flatnonzero(S_mask)
    B = pairwise_pearson(np.asarray(X)[:, gidx[s_pos]]).astype(np.float32)
    mirror_upper(B)
    return B, None, dict(pearson=round(time.perf_counter() - t, 2)), {}


def block_from_matrix(X, gidx, method: str, S_mask, K_all=50, config=None, rect=False):
    if method == "zavg_BEF":
        return zavg_block(X, gidx, S_mask, K_all, config, rect=rect)
    if method == "gls":
        return gls_block(X, gidx, S_mask)
    if method == "pearson":
        return pearson_block(X, gidx, S_mask)
    raise ValueError(f"method must be one of {METHODS}")


# ------------------------------------------------------------------------------------------ DepMap wrapper
class Block:
    """S x S score block of one run. ``vals`` float32 symmetric (NaN diagonal); ``s_index`` = genes.tsv indices of
    S (ascending); ``topk_all`` = (idx, val) over all protein genes for S rows (zavg only) or None."""

    def __init__(self, vals, s_index, meta, topk_all=None, rect=None):
        self.vals, self.s_index, self.meta, self.topk_all = vals, np.asarray(s_index), meta, topk_all
        self.rect = rect      # (vals: non-S genes x S float32, genes.tsv indices of the non-S rows) or None

    def topk(self, K: int = 50):
        """Exact top-K within S (local indices into S, values), self excluded."""
        return select_topk(self.vals, K)


def cache_path(method, model_ids, S_mask):
    from part2.data import model_set_hash
    return paths.CACHE / "blocks" / f"{method}_{model_set_hash(model_ids)}_S{mask_hash(S_mask)}.npz"


def load_block(path) -> Block:
    d = np.load(path)
    meta = json.loads(path.with_suffix(".json").read_text())
    tk = (d["topk_idx"], d["topk_val"]) if "topk_idx" in d.files else None
    rc = (d["rect_vals"], d["rect_index"]) if "rect_vals" in d.files else None
    return Block(d["vals"], d["s_index"], meta, tk, rc)


def score_block(model_ids, method: str, S_mask, K_all: int = 50, use_cache: bool = True,
                force: bool = False, rect: bool = False) -> Block:
    """Run the locked scorer (zavg_BEF) or Part 1's GLS on the given DepMap model set and return the S x S block
    over the gene set S (boolean mask over Part 1 genes.tsv). Serial; cached; sidecar records runtime / RSS."""
    from part2 import data
    model_ids = sorted(map(str, model_ids))
    S_mask = np.asarray(S_mask, bool)
    out = cache_path(method, model_ids, S_mask)
    if use_cache and out.exists() and not force:
        blk = load_block(out)
        if not rect or blk.rect is not None:
            return blk
    t0 = time.perf_counter()
    ge = data.load_gene_effect()
    genes = data.protein_genes()
    if len(S_mask) != len(genes):
        raise ValueError("S_mask must be over genes.tsv")
    if not (ge.columns[genes.depmap_index.values] == genes.depmap_column.values).all():
        raise RuntimeError("genes.tsv depmap_index does not match CRISPRGeneEffect columns")
    X = ge.subset(model_ids)
    t_load = time.perf_counter() - t0
    gidx = genes.depmap_index.values.astype(np.int64)
    vals, tk, tim, info = block_from_matrix(X, gidx, method, S_mask, K_all, rect=rect)
    del X
    rc = info.pop("_rect", None)
    if rc is not None:
        rc = (rc[0], np.asarray(rc[1], np.int64))          # positions into gidx == genes.tsv index
    s_index = np.flatnonzero(S_mask)
    meta = dict(method=method, model_set_hash=data.model_set_hash(model_ids), n_models=len(model_ids),
                model_ids=model_ids, S_hash=mask_hash(S_mask), n_S=int(S_mask.sum()), n_genes=len(genes),
                gene_index="Part 1 genes.tsv `index`", gene_effect_sha256=ge.csv_sha256,
                runtime_s=round(time.perf_counter() - t0, 1), load_s=round(t_load, 1), stages_s=tim,
                peak_rss_gb=round(peak_rss_gb(), 2), code=code_commit(), python=platform.python_version(),
                numpy=np.__version__, created=time.strftime("%Y-%m-%dT%H:%M:%S"), **info)
    if use_cache:
        out.parent.mkdir(parents=True, exist_ok=True)
        arrs = dict(vals=vals, s_index=s_index)
        if tk is not None:
            arrs.update(topk_idx=tk[0], topk_val=tk[1])
        if rc is not None:
            arrs.update(rect_vals=rc[0], rect_index=rc[1])
        np.savez(out, **arrs)
        out.with_suffix(".json").write_text(json.dumps(meta, indent=1, default=str))
    log.info("score_block %s n_models=%d |S|=%d: %.1f s, peak RSS %.2f GB", method, len(model_ids),
             len(s_index), meta["runtime_s"], peak_rss_gb())
    return Block(vals, s_index, meta, tk, rc)


def scores_topk(model_ids, method: str = "zavg_BEF", K: int = 50, S_mask=None, **kw):
    """Top-K neighbour lists (genes.tsv indices, float32 scores) for the S rows, exact within S.
    Derived from the cached S x S block (S_mask default: dependency-bearing filter on the given models)."""
    from part2 import genefilter
    if S_mask is None:
        S_mask = genefilter.dependency_bearing_for_models(model_ids)
    blk = score_block(model_ids, method, S_mask, **kw)
    I, V = blk.topk(K)
    return np.where(I >= 0, blk.s_index[np.clip(I, 0, None)], -1).astype(np.int32), V, blk.meta
