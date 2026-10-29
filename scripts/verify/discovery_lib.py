"""Independent re-implementation for the discovery-run review (H3a / H3b / H1). Does NOT import src/part2.

Built from: registry definitions (docs/registry/*.yaml), the frozen units (results/part2/fzx/*), the frozen context
table, the A model list (results/part2/lockbox/A_models.txt) and CRISPRGeneEffect restricted to A (the float32 cache
of the CSV, rows filtered to A on load). Covariate sources: Model.csv (OncotreeLineage, ModelType, PatientID) and
CRISPRConfounders.csv (Library*, ScreenType3DO), as named in the registry covariate_rule.

Statistics written from the definitions:
  OLS fit + HC3 SE via statsmodels; Satterthwaite df (Bell-McCaffrey, normal working model) via numpy:
  A = M diag(a) M with a_i = (g_ij / (1 - h_i))^2, g_j = row j of (X'X)^-1 X'; df = tr(A)^2 / tr(A^2).
  BH step-up, Simes, Benjamini-Bogomolov (families by BH over Simes p at q; within-family BH at q R / M).
Design convention mirrored (from the registry covariate_rule plus the engine docstring, which the registry does not
spell out): one-hot FE drops the first sorted level as reference; indicator columns with < 2 models on either side
are merged into the reference; constant columns dropped; collinear columns dropped greedily in order.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats as ss

ROOT = Path(__file__).resolve().parents[2]
CLIP = 3.0


# ------------------------------------------------------------------------------------------------ data
def a_models():
    return [l for l in (ROOT / "results/part2/lockbox/A_models.txt").read_text().split() if l]


def l_models():
    return [l for l in (ROOT / "results/part2/lockbox/lockbox_models.txt").read_text().split() if l]


def effects_A():
    side = next((ROOT / "data_cache/part2").glob("CRISPRGeneEffect.*.f32.json"))
    meta = json.loads(side.read_text())
    X = np.load(side.with_suffix(".npy"), mmap_mode="r")
    A = a_models()
    pos = {m: i for i, m in enumerate(meta["models"])}
    rows = [pos[m] for m in A if m in pos]
    idx = [m for m in A if m in pos]
    E = pd.DataFrame(np.asarray(X[rows, :], dtype=np.float64), index=idx, columns=meta["columns"])
    assert set(E.index) <= set(A) and not (set(E.index) & set(l_models()))
    return E


def context_table():
    return pd.read_parquet(ROOT / "data_cache/part2/contexts/context_table.parquet")


def cov_source(models):
    m = pd.read_csv(ROOT / "datasets/DepMAP/Model.csv", low_memory=False).set_index("ModelID").reindex(models)
    cf = pd.read_csv(ROOT / "datasets/DepMAP/CRISPRConfounders.csv", index_col=0).reindex(models)
    ct = context_table().reindex(models)
    out = pd.DataFrame(index=models)
    out["OncotreeLineage"] = m.OncotreeLineage
    out["PatientID"] = m.PatientID
    out["ModelType"] = m.ModelType
    for c in cf.columns:
        if c.startswith("Library") or c == "ScreenType3DO":
            out[c] = cf[c]
    for c in ("MSIScore", "log_burden", "Ploidy", "WGD", "TP53_altered"):
        out[c] = ct[c]
    return out


def membership():
    mem = pd.read_csv(ROOT / "results/part2/fzx/zavg_engine_membership.tsv", sep="\t")
    units = [l for l in (ROOT / "results/part2/fzx/units_zavg_coarse.txt").read_text().split() if l]
    m1 = mem[(mem.level == 1) & (mem.role == "member") & mem.unit_id.isin(units)]
    members = {u: list(g.gene) for u, g in m1.groupby("unit_id")}
    locus = mem[mem.level == 1].groupby("unit_id").locus.first().to_dict()
    return units, members, locus


# ------------------------------------------------------------------------------------------------ scores
def score(E, genes):
    """S_C: mean of members' effect clipped to [-3, 3]; NaN members skipped; NaN if < 50% observed."""
    genes = [g for g in genes if g in E.columns]
    V = E[genes].clip(-CLIP, CLIP)
    obs = V.notna().sum(axis=1)
    s = V.mean(axis=1, skipna=True)
    s[obs < 0.5 * len(genes)] = np.nan
    return s


# ------------------------------------------------------------------------------------------------ design
FE = {"lineage_FE": "OncotreeLineage", "ModelType": "ModelType"}


def cov_frame(names, src):
    parts = []
    for n in names:
        if n in FE:
            col = src[FE[n]]
            levels = sorted(col.dropna().unique())
            d = pd.DataFrame({f"{n}:{lv}": (col == lv).astype(float) for lv in levels[1:]}, index=src.index)
            d[col.isna()] = np.nan
            parts.append(d)
        elif n == "library":
            parts.append(src[[c for c in sorted(src.columns) if c.startswith("Library")]].astype(float)
                         .add_prefix("library:"))
        else:
            parts.append(src[[n]].astype(float))
    return pd.concat(parts, axis=1) if parts else pd.DataFrame(index=src.index)


def design(ctx, y, Z):
    df = pd.concat([pd.Series(1.0, index=ctx.index, name="intercept"), ctx.astype(float).rename("context"), Z], axis=1)
    ok = df.notna().all(axis=1) & y.notna()
    D = df[ok]
    keep, dropped = [], []
    for c in D.columns:
        v = D[c].values
        if c not in ("intercept", "context"):
            if np.all(v == v[0]):
                dropped.append(c); continue
            if c.split(":")[0] in ("lineage_FE", "ModelType", "library") and (min((v != 0).sum(), (v == 0).sum()) < 2):
                dropped.append(c); continue
        trial = keep + [c]
        if np.linalg.matrix_rank(D[trial].values) < len(trial):
            dropped.append(c); continue
        keep = trial
    return D[keep].values, y[ok].values.astype(float), ok, dropped, keep


def hc3_satt(X, y, j=1):
    fit = sm.OLS(y, X).fit(cov_type="HC3")
    beta, se = float(fit.params[j]), float(fit.bse[j])
    XtXi = np.linalg.inv(X.T @ X)
    G = XtXi @ X.T                       # rows = coefficient weight vectors g_j
    H = X @ G
    h = np.diag(H)
    M = np.eye(len(y)) - H
    a = (G[j] / (1 - h)) ** 2
    trA = float(np.sum(a * (1 - h)))
    trA2 = float(a @ (M ** 2) @ a)
    dfree = trA ** 2 / trA2
    t = beta / se
    p = 2 * ss.t.sf(abs(t), dfree)
    return dict(beta=beta, se=se, t=t, df=dfree, p=p, n=len(y), n_context=int((X[:, j] != 0).sum()))


def fit_test(ctx, y, Z):
    X, yy, ok, dropped, keep = design(ctx, y, Z)
    r = hc3_satt(X, yy, keep.index("context"))
    r["dropped"] = dropped
    return r


# ------------------------------------------------------------------------------------------------ FDR
def bh_reject(p, q):
    p = np.asarray(p, float)
    m = len(p)
    o = np.argsort(p)
    thr = q * np.arange(1, m + 1) / m
    ok = p[o] <= thr
    rej = np.zeros(m, bool)
    if ok.any():
        k = np.max(np.nonzero(ok)[0])
        rej[o[: k + 1]] = True
    return rej


def bh_adjust(p):
    p = np.asarray(p, float)
    m = len(p)
    o = np.argsort(p)
    adj = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    out = np.empty(m)
    out[o] = np.minimum(adj, 1)
    return out


def simes(p):
    p = np.sort(np.asarray(p, float))
    return float(np.min(p * len(p) / np.arange(1, len(p) + 1)))


def bb(p, fam, q=0.05):
    p = np.asarray(p, float)
    fam = np.asarray(fam, object)
    fams = list(dict.fromkeys(fam))
    M = len(fams)
    sp = np.array([simes(p[fam == f]) for f in fams])
    sel = bh_reject(sp, q)
    R = int(sel.sum())
    rej = np.zeros(len(p), bool)
    for f, s in zip(fams, sel):
        if s:
            ix = np.nonzero(fam == f)[0]
            rej[ix] = bh_reject(p[ix], q * R / M)
    return rej, R, dict(zip(fams, sp))
