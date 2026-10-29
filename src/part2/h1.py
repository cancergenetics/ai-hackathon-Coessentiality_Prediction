"""H1 within-complex genetic interactions (DESIGN v3.1 §2.6; docs/registry/H1_within_complex_GI.yaml).

Pipeline
  loss_calls        rules 1-2 of the enumeration rule: loss(g, m) = 1 / 0 / NaN and the loss kind per model.
  enumerate_events  rules 3-8 (labels, frozen map, attachments, paralogs, coordinates only; no dependency data).
  component_scores  member-resolved scores on A (y_P, aligned y_N, S_{C minus g}, raw N mean, baselines b_j) and the
                    per-context score frames + (context, unit) pairs for assoc.run_association (BH pooled over all
                    counted tests: kappa_N for every event with an evaluable N axis, beta_P for every |P| >= 1).
  classify          class calls per axis with TOST (+/-0.1 Chronos), P-vs-N contrast, member breadth and the matched
                    random-set null (K = 1000, seed 7000 + event index), plus the collapse index.

Alignment sign. This module implements
    y_N(m) = mean_j sign(b_j) * effect_j(m)      (ALIGN = +1)
i.e. the signed dependency magnitude in the baseline direction, so kappa_N > 0 when dependency grows (registry
stated meaning and class rules). The registry formula text read '-sign(b_j)' (a sign typo); corrected by erratum E1
(PI ruling R1, docs/registry/AMENDMENTS.md), semantics unchanged.

Collapse index = -kappa_N / mean_{j in N'} |b_j|: the fraction of the baseline dependency magnitude lost in loss
lines (1 = members' dependency fully gone; < 0 = sensitisation). Reported, not a class condition.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import stats as sst

from part2 import stats
from part2.cmap import FrozenMap, symbol
from part2.cscore import CLIP, Annotation, mean_score

ALIGN = +1.0
TOST_MARGIN = 0.1
MIN_BASELINE = 0.1


# ------------------------------------------------------------------------------------------------ rules 1-2
def loss_calls(lof: pd.DataFrame, homdel: pd.DataFrame, damaging: pd.DataFrame, expr: pd.DataFrame,
               lineage: pd.Series, lof_pass: set | None = None, silencing_min_lineage: int = 5):
    """All inputs models x gene symbols (aligned); lof = primary LoF call (MSI-H lines NaN), homdel = hom-del call,
    damaging = damaging-mutation count, expr = log1p TPM (NaN = not measured). lof_pass: genes passing the LoF
    functional check (None = all). Returns (loss 1/0/NaN, kind 'homdel' / 'LoF' / 'silencing' / '')."""
    genes = list(lof.columns)
    idx = lof.index
    lof, homdel, damaging, expr = (d.reindex(index=idx, columns=genes).astype(float)
                                   for d in (lof, homdel, damaging, expr))
    passes = np.array([lof_pass is None or g in lof_pass for g in genes])
    lin = lineage.reindex(idx)
    med_u = expr.median()
    med = pd.DataFrame(np.tile(med_u.values, (len(idx), 1)), index=idx, columns=genes)
    for l, rows in lin.groupby(lin).groups.items():
        sub = expr.loc[rows]
        ok = sub.notna().sum() >= silencing_min_lineage
        med.loc[rows, ok[ok].index] = np.tile(sub.median()[ok].values, (len(rows), 1))
    a = (lof == 1).values & passes[None, :]
    b = (homdel == 1).values
    c = (expr < 1).values & (med >= 3).values & ~b
    one = a | b | c
    zero = (damaging == 0).values & (homdel == 0).values & ~(expr < 1).values
    loss = np.where(one, 1.0, np.where(zero, 0.0, np.nan))
    kind = np.where(b, "homdel", np.where(a, "LoF", np.where(c, "silencing", "")))
    kind = np.where(one, kind, "")
    return pd.DataFrame(loss, index=idx, columns=genes), pd.DataFrame(kind, index=idx, columns=genes)


def event_type(kinds: pd.Series) -> str:
    """Rule 6 on the kinds of the A loss models: CN-type if >= 50% hom-del (ties with LoF go to CN), else LoF-type
    if >= 50% LoF (mutation-based LoF without hom-del), else silencing-type."""
    k = kinds[kinds != ""]
    if not len(k):
        return "silencing"
    sh = (k == "homdel").mean()
    sl = (k == "LoF").mean()
    if sh >= 0.5:
        return "CN"
    if sl >= 0.5:
        return "LoF"
    return "silencing"


COVARIATES = {"LoF": ["lineage_FE", "library", "ScreenType3DO", "ModelType", "log_burden", "MSIScore"],
              "CN": ["lineage_FE", "library", "ScreenType3DO", "ModelType", "Ploidy", "WGD", "TP53_altered"],
              "silencing": ["lineage_FE", "library", "ScreenType3DO", "ModelType"]}


# ------------------------------------------------------------------------------------------------ rules 3-8
@dataclass
class Event:
    index: int
    id: str
    unit: str
    gene: str
    link: str
    type: str
    P: list
    N: list
    n_A: int
    ref_A: int
    n_L: int = 0
    ref_L: int = 0
    testable: bool = False
    locus: bool = False
    covariates: list = field(default_factory=list)

    @property
    def context(self) -> str:
        return f"H1loss:{self.gene}"

    @property
    def axes(self) -> list:
        return (["N"] if len(self.N) >= 2 else []) + (["P"] if len(self.P) >= 1 else [])


def enumerate_events(fmap: FrozenMap, loss: pd.DataFrame, kind: pd.DataFrame, annot: Annotation, fold: pd.Series,
                     min_A: int = 12, min_L: int = 5, window: int = 10_000_000) -> list[Event]:
    """fold: ModelID -> 'A' / 'L' (labels only). Events sorted by (unit, gene); index = list position."""
    isA = fold.reindex(loss.index).eq("A").values
    isL = fold.reindex(loss.index).eq("L").values
    out = []
    for g in sorted(loss.columns):
        v = loss[g].values
        nA, rA = int(((v == 1) & isA).sum()), int(((v == 0) & isA).sum())
        if nA < min_A or rA < min_A:
            continue
        typ = event_type(kind[g][(v == 1) & isA])
        nL, rL = int(((v == 1) & isL).sum()), int(((v == 0) & isL).sum())
        links = [(u, "member") for u in fmap.units_of_symbol(g, 1)]
        links += [(u, "attached") for u, att in sorted(fmap.attached.items())
                  if g in {symbol(x) for x in att} and (u, "member") not in links]
        for u, how in links:
            rest = [x for x in fmap.members[u] if symbol(x) != g]
            P = [x for x in rest if annot.is_paralog(symbol(x), g)]
            N = [x for x in rest if x not in P and not (typ == "CN" and annot.near(symbol(x), g, window))]
            if len(N) < 2 and len(P) < 1:
                continue
            out.append(Event(index=-1, id=f"{u}|{g}", unit=u, gene=g, link=how, type=typ, P=P, N=N, n_A=nA, ref_A=rA,
                             n_L=nL, ref_L=rL, testable=bool(nL >= min_L and rL >= min_L),
                             locus=bool(fmap.locus.get(u, False)), covariates=list(COVARIATES[typ])))
    out.sort(key=lambda e: (e.unit, e.gene))
    for i, e in enumerate(out):
        e.index = i
    return out


def events_table(events) -> pd.DataFrame:
    return pd.DataFrame([dict(index=e.index, event=e.id, unit=e.unit, gene=e.gene, link=e.link, type=e.type,
                              n_P=len(e.P), n_N=len(e.N), P=";".join(e.P), N=";".join(e.N), n_A=e.n_A, ref_A=e.ref_A,
                              n_L=e.n_L, ref_L=e.ref_L, testable=e.testable, locus=e.locus,
                              covariates=";".join(e.covariates)) for e in events])


def units_of(events) -> list:
    return sorted([f"{e.id}|{a}" for e in events for a in e.axes])


# ------------------------------------------------------------------------------------------------ scores
def baselines(E: pd.DataFrame, genes, loss_col: pd.Series) -> pd.Series:
    """b_j = median clipped effect over the reference (loss = 0) rows of E."""
    ref = loss_col.reindex(E.index) == 0
    return E.loc[ref.values, list(genes)].astype(float).clip(-CLIP, CLIP).median()


def aligned_score(E: pd.DataFrame, genes, b: pd.Series, min_frac: float = 0.5) -> pd.Series:
    V = E.reindex(columns=list(genes)).astype(float).clip(-CLIP, CLIP) * (ALIGN * np.sign(b.reindex(genes).values))
    s = V.mean(axis=1)
    s[(V.notna().sum(axis=1) < min_frac * len(genes)).values] = np.nan
    return s


@dataclass
class Components:
    event: Event
    b: pd.Series            # baselines of N (all of N, incl. |b| < 0.1)
    Nprime: list            # N members with |b_j| >= 0.1
    y: dict                 # 'N' aligned, 'P', 'all' (S_{C minus g}), 'Nraw', 'contrast' -> Series over E rows


def component_scores(E: pd.DataFrame, event: Event, loss_col: pd.Series, b: pd.Series | None = None) -> Components:
    """Scores over E's rows. b (baselines) defaults to A reference medians on E; pass the A baselines for L rows."""
    if b is None:
        b = baselines(E, event.N, loss_col) if event.N else pd.Series(dtype=float)
    Np = [g for g in event.N if abs(b.get(g, 0.0)) >= MIN_BASELINE]
    y = {}
    if len(Np) >= 2:
        y["N"] = aligned_score(E, Np, b)
    if event.P:
        y["P"] = mean_score(E, event.P)
    if event.N:
        y["Nraw"] = mean_score(E, event.N)
    y["all"] = mean_score(E, event.P + event.N)
    if "P" in y and "Nraw" in y:
        y["contrast"] = y["P"] - y["Nraw"]
    return Components(event=event, b=b, Nprime=Np, y=y)


def assoc_inputs(E_A: pd.DataFrame, events, loss_A: pd.DataFrame):
    """(Y frames per context, context frame, pairs, components) for run_association on the H1 entry. Units: every
    '<event>|N' with |N'| >= 2 and every '<event>|P'; the unit set hash is part2.registry.set_hash(units_of)."""
    comps, cols = {}, {}
    for e in events:
        c = component_scores(E_A, e, loss_A[e.gene])
        comps[e.id] = c
        if "N" in c.y:
            cols[f"{e.id}|N"] = (e.context, c.y["N"])
        if "P" in c.y:
            cols[f"{e.id}|P"] = (e.context, c.y["P"])
    all_units = units_of(events)
    ctxs = sorted({e.context for e in events})
    Y = {}
    for cx in ctxs:
        df = pd.DataFrame(np.nan, index=E_A.index, columns=all_units)
        for u, (c2, s) in cols.items():
            if c2 == cx:
                df[u] = s
        Y[cx] = df
    contexts = pd.DataFrame({f"H1loss:{g}": loss_A[g] for g in sorted({e.gene for e in events})}).reindex(E_A.index)
    pairs = [(c2, u) for u, (c2, _) in cols.items()]
    return Y, contexts, pairs, comps


# ------------------------------------------------------------------------------------------------ fitting helpers
def _xy(ctx: pd.Series, cov: pd.DataFrame | None, idx):
    from part2.assoc import _design
    X, names, ok = _design(ctx.reindex(idx), None if cov is None else cov.reindex(idx))
    return X, ok


def fit_columns(X: np.ndarray, Ymat: np.ndarray, alternative: str = "two-sided") -> pd.DataFrame:
    """HC3 fit of column 1 of X for every column of Ymat (n x K; NaN-pattern groups)."""
    out = np.full((Ymat.shape[1], 5), np.nan)
    nanpat = np.isnan(Ymat)
    keys = [nanpat[:, k].tobytes() for k in range(Ymat.shape[1])]
    for key, cols in pd.Series(range(Ymat.shape[1])).groupby(keys):
        cols = cols.values
        r = ~nanpat[:, cols[0]]
        if r.sum() < X.shape[1] + 2 or X[r, 1].sum() < 2:
            continue
        try:
            d = stats.HC3Design(X[r], 1)
        except (np.linalg.LinAlgError, ValueError):
            continue
        b, se, t = d.fit(Ymat[r][:, cols])
        out[cols] = np.column_stack([b, se, t, np.full(len(cols), d.df), stats.pvalue_t(t, d.df, alternative)])
    return pd.DataFrame(out, columns=["beta", "se", "t", "df", "p"])


def ci(beta, se, df, level: float) -> tuple:
    h = sst.t.ppf(0.5 + level / 2, df) * se
    return beta - h, beta + h


# ------------------------------------------------------------------------------------------------ matched null
def coherence(genes, pair_scores: pd.DataFrame | None) -> float:
    if pair_scores is None or len(genes) < 2:
        return np.nan
    M = pair_scores.reindex(index=list(genes), columns=list(genes)).values.astype(float)
    iu = np.triu_indices(len(genes), 1)
    return float(np.nanmean(M[iu]))


def matched_null(E_A: pd.DataFrame, fmap: FrozenMap, event: Event, axis: str, comp: Components, ctx: pd.Series,
                 cov: pd.DataFrame | None, annot: Annotation, essentiality: pd.Series,
                 pair_scores: pd.DataFrame | None, seed0: int = 7000, K: int = 1000, max_attempts: int = 100_000,
                 min_sets: int = 200, tol_ess: float = 0.1, tol_coh: float = 1.0, window: int = 10_000_000) -> dict:
    """Null distribution of the axis statistic over random member sets matched on size, essentiality (+/- tol_ess)
    and coherence (+/- tol_coh; not used for single-gene sets). Seed seed0 + event index."""
    genes = event.P if axis == "P" else comp.Nprime
    k = len(genes)
    t_ess = float(essentiality.reindex(genes).mean())
    t_coh = coherence(genes, pair_scores)
    lossg = event.gene
    excl_units = {event.unit} | {u for u in fmap.coarse()
                                 if any(symbol(x) == lossg or annot.is_paralog(symbol(x), lossg)
                                        for x in fmap.members[u])}
    pool_units = []
    for u in fmap.coarse():
        if u in excl_units:
            continue
        mem = [x for x in fmap.members[u] if x in E_A.columns
               and not (event.type == "CN" and annot.near(symbol(x), lossg, window))]
        if mem:
            pool_units.append(mem)
    ref = ctx.reindex(E_A.index) == 0
    allg = sorted({x for m in pool_units for x in m})
    bpool = E_A.loc[ref.values, allg].astype(float).clip(-CLIP, CLIP).median()
    if axis == "N":                                       # members need a baseline direction to be aligned
        pool_units = [[x for x in m if abs(bpool[x]) >= MIN_BASELINE] for m in pool_units]
        pool_units = [m for m in pool_units if m]
    rng = np.random.default_rng(seed0 + event.index)
    gi = {g: i for i, g in enumerate(allg)}
    pool_idx = [np.array([gi[x] for x in m]) for m in pool_units]
    ess_v = essentiality.reindex(allg).values.astype(float)
    PS = None
    if pair_scores is not None and k >= 2 and np.isfinite(t_coh):
        PS = pair_scores.reindex(index=allg, columns=allg).values.astype(float)
        iu = np.triu_indices(k, 1)
    sizes = np.array([len(m) for m in pool_idx])
    big = np.flatnonzero(sizes >= k)
    sets, status = [], "ok"
    for tol_mult in (1.0, 2.0):
        seen, sets = set(), []
        te, tc = tol_ess * tol_mult, tol_coh * tol_mult
        for _ in range(max_attempts):
            if len(big):
                cand = pool_idx[big[rng.integers(len(big))]]
            elif len(pool_idx) >= 2:
                i, j = rng.choice(len(pool_idx), 2, replace=False)
                cand = np.concatenate([pool_idx[i], pool_idx[j]])
                if len(cand) < k:
                    continue
            else:
                break
            s = np.sort(rng.choice(cand, k, replace=False))
            key = s.tobytes()
            if key in seen:
                continue
            if abs(float(np.nanmean(ess_v[s])) - t_ess) > te:
                continue
            if PS is not None and not abs(float(np.nanmean(PS[np.ix_(s, s)][iu])) - t_coh) <= tc:
                continue
            seen.add(key)
            sets.append(s)
            if len(sets) >= K:
                break
        if len(sets) >= min_sets:
            status = "ok" if tol_mult == 1.0 else "tolerances doubled"
            break
    else:
        status = "null unavailable"
    if status == "null unavailable" or not sets:
        return dict(status="null unavailable", n_sets=len(sets), q025=np.nan, q975=np.nan, null=np.array([]))
    V = E_A[allg].astype(float).clip(-CLIP, CLIP).values
    if axis == "N":
        V = V * (ALIGN * np.sign(bpool.reindex(allg).values))[None, :]
    Ymat = np.empty((V.shape[0], len(sets)))
    for j, s in enumerate(sets):
        sub = V[:, s]
        nobs = np.isfinite(sub).sum(1)
        with np.errstate(invalid="ignore", divide="ignore"):
            Ymat[:, j] = np.where(nobs >= 0.5 * k, np.nansum(sub, 1) / np.maximum(nobs, 1), np.nan)
    X, ok = _xy(ctx, cov, E_A.index)
    f = fit_columns(X, Ymat[ok])
    nb = f.beta.dropna().values
    return dict(status=status, n_sets=len(sets), q025=float(np.quantile(nb, 0.025)),
                q975=float(np.quantile(nb, 0.975)), null=nb, target_ess=t_ess, target_coh=t_coh)


# ------------------------------------------------------------------------------------------------ classes
def member_profile(E_A, genes, ctx, cov, b: pd.Series | None = None) -> pd.DataFrame:
    X, ok = _xy(ctx, cov, E_A.index)
    Ymat = E_A.reindex(columns=list(genes)).astype(float).clip(-CLIP, CLIP).values[ok]
    f = fit_columns(X, Ymat)
    f.index = list(genes)
    lo, hi = ci(f.beta.values, f.se.values, f.df.values, 0.95)
    f["ci95_lo"], f["ci95_hi"] = lo, hi
    if b is not None:
        f["baseline"] = b.reindex(genes).values
        f["aligned_beta"] = ALIGN * np.sign(f.baseline) * f.beta
    return f


def _axis_fit(y: pd.Series, ctx, cov, alternative="two-sided") -> dict:
    X, ok = _xy(ctx, cov, y.index)
    f = fit_columns(X, y.values[ok][:, None], alternative).iloc[0].to_dict()
    lo90, hi90 = ci(f["beta"], f["se"], f["df"], 0.90)
    f["ci90"] = (lo90, hi90)
    f["tost_neutral"] = bool(np.isfinite(lo90) and lo90 > -TOST_MARGIN and hi90 < TOST_MARGIN)
    return f


def classify(events, comps: dict, res: pd.DataFrame, E_A: pd.DataFrame, loss_A: pd.DataFrame, covariates: dict,
             fmap: FrozenMap, annot: Annotation, essentiality: pd.Series, pair_scores: pd.DataFrame | None,
             q: float = 0.05, K: int = 1000, seed0: int = 7000, max_attempts: int = 100_000) -> pd.DataFrame:
    """One row per event: both axes' statistics, q, null quantiles, TOST, breadth, contrast, classes."""
    qcol = "q_adj_perm" if "q_adj_perm" in res else "q_adj"
    rmap = {r.unit: r for r in res.itertuples(index=False)}
    rows = []
    for e in events:
        c = comps[e.id]
        ctx = loss_A[e.gene].rename(e.context)
        cov = covariates.get(e.context)
        row = dict(event=e.id, unit=e.unit, gene=e.gene, type=e.type, link=e.link, n_A=e.n_A, ref_A=e.ref_A,
                   testable=e.testable, locus=e.locus, n_P=len(e.P), n_Nprime=len(c.Nprime),
                   N_no_baseline=";".join(g for g in e.N if g not in c.Nprime))
        f_all = _axis_fit(c.y["all"], ctx, cov)
        row["beta_all"] = f_all["beta"]
        if "Nraw" in c.y:
            row["beta_Nraw"] = _axis_fit(c.y["Nraw"], ctx, cov)["beta"]
        # ---- complex-wide axis
        cls_n = "not evaluable"
        if "N" in c.y:
            fN = _axis_fit(c.y["N"], ctx, cov)
            r = rmap.get(f"{e.id}|N")
            qn = float(getattr(r, qcol)) if r is not None and pd.notna(getattr(r, qcol)) else np.nan
            disc = bool(r.discovery) if r is not None else False
            prof = member_profile(E_A, c.Nprime, ctx, cov, c.b)
            frac_pos = float((prof.aligned_beta > 0).mean())
            frac_neg = float((prof.aligned_beta < 0).mean())
            nl = matched_null(E_A, fmap, e, "N", c, ctx, cov, annot, essentiality, pair_scores, seed0=seed0, K=K,
                              max_attempts=max_attempts)
            kap = fN["beta"]
            if nl["status"] == "null unavailable":
                cls_n = "INCONCLUSIVE (null unavailable)"
            elif kap > 0 and disc and qn < q and kap > nl["q975"] and frac_pos >= 0.5:
                cls_n = "RESIDUAL SENSITISATION"
            elif kap < 0 and disc and qn < q and kap < nl["q025"] and frac_neg >= 0.5:
                cls_n = "COLLAPSE"
            elif fN["tost_neutral"]:
                cls_n = "NEUTRAL"
            else:
                cls_n = "INCONCLUSIVE"
            row.update(kappa_N=kap, kappa_N_se=fN["se"], kappa_N_p=fN["p"], kappa_N_q=qn,
                       kappa_N_ci90_lo=fN["ci90"][0], kappa_N_ci90_hi=fN["ci90"][1], N_tost=fN["tost_neutral"],
                       N_null_q025=nl["q025"], N_null_q975=nl["q975"], N_null_status=nl["status"],
                       N_null_sets=nl["n_sets"], N_frac_aligned_pos=frac_pos, N_frac_aligned_neg=frac_neg,
                       collapse_index=float(-kap / np.mean(np.abs(c.b.reindex(c.Nprime).values))),
                       N_most_shifted=str(prof.beta.abs().idxmax()))
        # ---- paralog axis
        cls_p = "not evaluable"
        if "P" in c.y:
            fP = _axis_fit(c.y["P"], ctx, cov)
            r = rmap.get(f"{e.id}|P")
            qp = float(getattr(r, qcol)) if r is not None and pd.notna(getattr(r, qcol)) else np.nan
            disc = bool(r.discovery) if r is not None else False
            con = _axis_fit(c.y["contrast"], ctx, cov, "less") if "contrast" in c.y else None
            con_ok = bool(con is not None and con["beta"] < 0 and con["p"] < 0.05)
            nl = matched_null(E_A, fmap, e, "P", c, ctx, cov, annot, essentiality, pair_scores, seed0=seed0, K=K,
                              max_attempts=max_attempts)
            bp = fP["beta"]
            if nl["status"] == "null unavailable":
                cls_p = "INCONCLUSIVE (null unavailable)"
            elif bp < 0 and disc and qp < q and bp < nl["q025"] and con_ok:
                cls_p = "PARALOG COMPENSATION"
            elif fP["tost_neutral"]:
                cls_p = "NEUTRAL"
            elif bp > 0 and disc and qp < q:
                cls_p = "INCONCLUSIVE (paralog relief, exploratory)"
            else:
                cls_p = "INCONCLUSIVE"
            row.update(beta_P=bp, beta_P_se=fP["se"], beta_P_p=fP["p"], beta_P_q=qp, beta_P_ci90_lo=fP["ci90"][0],
                       beta_P_ci90_hi=fP["ci90"][1], P_tost=fP["tost_neutral"], P_null_q025=nl["q025"],
                       P_null_q975=nl["q975"], P_null_status=nl["status"], P_null_sets=nl["n_sets"],
                       contrast_beta=None if con is None else con["beta"],
                       contrast_p=None if con is None else con["p"])
        row["class_complex_wide"], row["class_paralog"] = cls_n, cls_p
        if cls_n == "COLLAPSE":
            row["class_note"] = "collapse: exploratory unless lockbox-confirmed"
        rows.append(row)
    return pd.DataFrame(rows)


def score_L_fn(events, comps: dict, loss_L: pd.DataFrame):
    """score_L callback for confirm.confirm_entry: same members, A baselines for the alignment."""
    ev = {e.id: e for e in events}

    def f(E_L, context, unit):
        eid, axis = unit.rsplit("|", 1)
        e, c = ev[eid], comps[eid]
        if axis == "P":
            return mean_score(E_L, e.P)
        return aligned_score(E_L, c.Nprime, c.b)
    return f
