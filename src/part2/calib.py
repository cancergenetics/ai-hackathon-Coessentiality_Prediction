"""Registered negative-control calibrations not covered by robust.pseudo_context_calibration (registry `calibration`
fields of H3a_lineage_atlas, H3b_driver_family and H2_submodule_divergence; PI ruling R4, disclosure D8).

* matched_random_sets   random gene sets matched to each frozen coarse unit on size, essentiality (mean clipped A effect,
                        +/- 0.1) and coherence (mean pairwise zavg_BEF on A, +/- 1.0 z), drawn from the members of another
                        frozen coarse unit (one unit with >= k members, else the union of two), tolerances doubled once
                        if fewer than R sets are accepted -- the registered H1 matched-null rule, with R sets per unit.
* dispersed_never_expressed  sets of genes never expressed in A (log1p TPM < 1 in every RNA-profiled A model), no two
                        members within 10 Mb (Ensembl coordinates), sizes drawn from the coarse-unit size distribution.
* segment_pseudo_contexts  CN segment-matched null (H3b CN contexts): segments of the same event type at other loci
                        (results/part2/contexts/segments_ge8.tsv; > 10 Mb from the context gene or on another
                        chromosome) whose A count is within +/- 25% of the context's n_A; members within 10 Mb of the
                        segment span are removed for that pseudo-context.
* random_splits         H2: random member subsets of the parent with the sub-module's size, residualised on their
                        complement exactly like h2.residualised.
Each set or pseudo-context is scored like the real units (clip +/- 3, 50% rule, >= 2 members after removals) and tested
with the context's covariates by the same HC3 / Satterthwaite machinery (h1.fit_columns on assoc._design).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from part2 import stats
from part2.cmap import FrozenMap, symbol
from part2.cscore import CLIP, WINDOW_BP, Annotation, Removal, mean_score
from part2.h1 import coherence, fit_columns


# ------------------------------------------------------------------------------------------------ testing
def test_columns(Y: pd.DataFrame, ctx: pd.Series, cov: pd.DataFrame | None, alternative: str = "two-sided",
                 min_n: int = 12) -> pd.DataFrame:
    """HC3 fit of ctx for every column of Y (models x sets); rows with < min_n context or reference models after
    complete case are NaN (registered evaluability). Returns beta, se, t, df, p per column (index = Y.columns)."""
    from part2.assoc import _design
    X, _, ok = _design(ctx.reindex(Y.index), None if cov is None else cov.reindex(Y.index))
    V = Y.values[ok].astype(float)
    f = fit_columns(X, V, alternative)
    f.index = Y.columns
    obs = ~np.isnan(V)
    n_ctx = (obs * X[:, [1]]).sum(0)
    n_ref = obs.sum(0) - n_ctx
    f.loc[(n_ctx < min_n) | (n_ref < min_n)] = np.nan
    return f


def score_sets(E: pd.DataFrame, sets: dict, removal: Removal | None = None, annot: Annotation | None = None,
               min_members: int = 2) -> pd.DataFrame:
    """{name: gene keys} -> models x names primary scores after `removal` (NaN column if < min_members remain)."""
    cols = {}
    for k, genes in sets.items():
        rem = set(removal.removed(genes, annot or Annotation())) if removal is not None else set()
        use = [g for g in genes if g not in rem and g in E.columns]
        cols[k] = mean_score(E, use) if len(use) >= min_members else pd.Series(np.nan, index=E.index)
    return pd.DataFrame(cols, index=E.index)


# ------------------------------------------------------------------------------------------------ matched random sets
def matched_random_sets(E_A: pd.DataFrame, fmap: FrozenMap, units, essentiality: pd.Series,
                        pair_scores: pd.DataFrame | None, R: int, seed: int, tol_ess: float = 0.1,
                        tol_coh: float = 1.0, max_attempts: int = 20_000) -> tuple[dict, pd.DataFrame]:
    """{f'{unit}#r{k}': gene keys} and a per-unit status table. Seed: numpy default_rng([seed, unit index])."""
    coarse = fmap.coarse()
    allg = sorted({g for u in coarse for g in fmap.members[u] if g in E_A.columns})
    gi = {g: i for i, g in enumerate(allg)}
    ess_v = essentiality.reindex(allg).values.astype(float)
    PS = None if pair_scores is None else pair_scores.reindex(index=allg, columns=allg).values.astype(float)
    pools = {u: np.array([gi[g] for g in fmap.members[u] if g in gi]) for u in coarse}
    out, rows = {}, []
    for ui, u in enumerate(units):
        mem = [g for g in fmap.members[u] if g in gi]
        k = len(mem)
        idx = np.array([gi[g] for g in mem])
        t_ess = float(np.nanmean(ess_v[idx]))
        t_coh = coherence(mem, pair_scores) if PS is not None else np.nan
        iu = np.triu_indices(k, 1)
        others = [pools[v] for v in coarse if v != u and len(pools[v])]
        big = [p for p in others if len(p) >= k]
        rng = np.random.default_rng([seed, ui])
        status, sets = "unavailable", []
        for mult in (1.0, 2.0):
            seen, sets = set(), []
            for _ in range(max_attempts):
                if big:
                    cand = big[rng.integers(len(big))]
                elif len(others) >= 2:
                    i, j = rng.choice(len(others), 2, replace=False)
                    cand = np.concatenate([others[i], others[j]])
                    if len(cand) < k:
                        continue
                else:
                    break
                s = np.sort(rng.choice(cand, k, replace=False))
                key = s.tobytes()
                if key in seen:
                    continue
                if abs(float(np.nanmean(ess_v[s])) - t_ess) > tol_ess * mult:
                    continue
                if PS is not None and np.isfinite(t_coh) and \
                        not abs(float(np.nanmean(PS[np.ix_(s, s)][iu])) - t_coh) <= tol_coh * mult:
                    continue
                seen.add(key)
                sets.append(s)
                if len(sets) >= R:
                    break
            if len(sets) >= R:
                status = "ok" if mult == 1.0 else "tolerances doubled"
                break
        if status == "unavailable":
            sets = []
        for r, s in enumerate(sets):
            out[f"{u}#r{r}"] = [allg[i] for i in s]
        rows.append(dict(unit=u, size=k, target_ess=t_ess, target_coh=t_coh, status=status, n_sets=len(sets)))
    return out, pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------------ never expressed
def never_expressed_pool(expr_A: pd.DataFrame, gene_keys, annot: Annotation) -> list:
    """Gene keys (screened) with log1p TPM < 1 in every RNA-profiled A model and known coordinates."""
    ex = expr_A.dropna(how="all")
    mx = ex.max(axis=0, skipna=True)
    never = set(mx.index[(mx < 1.0) & ex.notna().any(axis=0)])
    return sorted(g for g in gene_keys if g in never and annot.pos.get(symbol(g), -1) >= 0
                  and annot.chrom.get(symbol(g), ""))


def dispersed_sets(pool, sizes, annot: Annotation, rng: np.random.Generator, window: int = WINDOW_BP) -> dict:
    """One set per size; members drawn in random order, a gene within `window` of a chosen member is skipped."""
    out = {}
    for i, k in enumerate(sizes):
        chosen = []
        for j in rng.permutation(len(pool)):
            g = pool[j]
            if any(annot.near(symbol(g), symbol(h), window) for h in chosen):
                continue
            chosen.append(g)
            if len(chosen) == k:
                break
        out[f"never{i}_k{k}"] = chosen
    return out


# ------------------------------------------------------------------------------------------------ CN segment null
def segment_pseudo_contexts(seg_tab: pd.DataFrame, seg_A: pd.DataFrame, context_gene: str, event: str, n_target: int,
                            annot: Annotation, window: int = WINDOW_BP, tol: float = 0.25) -> tuple[pd.DataFrame, dict]:
    """(models x qualifying segments 1/0/NaN, {segment: symbols within window of its span}) for one CN context."""
    cg_chr, cg_pos = annot.chrom.get(context_gene, ""), annot.pos.get(context_gene, -1)
    keep = []
    for r in seg_tab[seg_tab.event == event].itertuples(index=False):
        if r.segment not in seg_A.columns:
            continue
        if str(r.chrom) == cg_chr and r.start_pos - window <= cg_pos <= r.end_pos + window:
            continue                                              # same locus as the real context
        n = int((seg_A[r.segment] == 1).sum())
        if abs(n - n_target) <= tol * n_target:
            keep.append(r)
    near = {}
    for r in keep:
        near[r.segment] = sorted(s for s, c in annot.chrom.items() if c == str(r.chrom) and
                                 r.start_pos - window <= annot.pos.get(s, -1) <= r.end_pos + window)
    return seg_A[[r.segment for r in keep]].astype(float), near


# ------------------------------------------------------------------------------------------------ H2 random splits
def random_split_scores(E_A: pd.DataFrame, fmap: FrozenMap, sub: str, removal: Removal | None, annot: Annotation,
                        R: int, rng: np.random.Generator) -> pd.DataFrame:
    """R random subsets of the parent's members (after removal) with the sub-module's size, residualised on their
    complement in the parent like h2.residualised (OLS over A)."""
    from part2.h2 import residualised
    p = fmap.parent[sub]
    rm = set(removal.removed(fmap.members[p], annot)) if removal else set()
    pm = [g for g in fmap.members[p] if g not in rm and g in E_A.columns]
    k = len([g for g in fmap.members[sub] if g in pm])
    cols = {}
    for r in range(R):
        s = sorted(rng.choice(pm, k, replace=False).tolist()) if 2 <= k < len(pm) - 1 else []
        if not s:
            cols[f"{sub}#split{r}"] = pd.Series(np.nan, index=E_A.index)
            continue
        tmp = FrozenMap(members={p: pm, f"{p}.R": s}, level={p: 1, f"{p}.R": 2}, parent={p: "", f"{p}.R": p})
        y, _, _ = residualised(E_A, tmp, f"{p}.R")
        cols[f"{sub}#split{r}"] = y
    return pd.DataFrame(cols, index=E_A.index)


def bb_rejections(p: np.ndarray, fam: np.ndarray, q: float) -> dict:
    fin = np.isfinite(p)
    if not fin.any():
        return dict(rejections=0, families_selected=0)
    bb = stats.benjamini_bogomolov(p[fin], np.asarray(fam)[fin], q)
    return dict(rejections=int(bb["reject"].sum()), families_selected=int(bb["R"]))


def clip_essentiality(E_A: pd.DataFrame) -> pd.Series:
    return E_A.clip(-CLIP, CLIP).mean()
