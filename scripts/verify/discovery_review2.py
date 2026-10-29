"""Independent review, part 2: H1S, H2, registered calibrations (D8), GPD sensitivity. Does NOT import src/part2.

    .venv/bin/python scripts/verify/discovery_review2.py --out <scratch>/dr2          # ~5-10 min

1. H1S / H2: recompute 5 + 5 analytic tests from the registry definitions (own HC3 + Bell-McCaffrey df);
   |N'| for every H1S event (evaluability: 338 registered -> 321 built -> 293 evaluable).
2. H2 stage-wise BB from the committed p_final (and analytic p) with own code; GPD p_cons for H1S / H2.
3. H3a matched random sets: (a) RNG-mirrored re-draw for a random subsample of units, re-tested with own code and
   compared to the committed calibration t / p; (b) an independent re-draw (own RNG stream, same rule) on the same
   units; (c) the 'discoveries exceeding all matched sets' count from the committed tables (full) and from own t
   (subsample).
4. Never-expressed pool from the raw expression CSV (own read); RNG-mirrored and independent sets, H3a rate.
5. H2 random member splits: own implementation (mirrored and independent RNG), stage-wise rejections per replicate.
"""
from __future__ import annotations

import argparse
import io
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import discovery_lib as L  # noqa: E402
import discovery_review as DR  # noqa: E402

SEED = 20261030
D = L.ROOT / "results/part2/directions"
H3A_COV = DR.H3A_COV
PART1 = Path(__import__("os").environ.get("PART1_DIR", L.ROOT.parent / "coessentiality_protein complexes"))


def tests(aid):
    return pd.read_csv(D / aid / "tests.tsv.gz", sep="\t")


# ------------------------------------------------------------------------------------------------ fast HC3 (many y)
def hc3_many(X, Y, j=1):
    """HC3 t / p for column j of X for every column of Y (no NaN in Y); Bell-McCaffrey df (depends on X only)."""
    from scipy import stats as ss
    XtXi = np.linalg.inv(X.T @ X)
    G = XtXi @ X.T
    h = np.einsum("ij,ji->i", X, G)
    B = G @ Y
    E = Y - X @ B
    w = (G[j] / (1 - h)) ** 2
    se = np.sqrt(w @ (E ** 2))
    M = np.eye(len(h)) - X @ G
    trA = float(np.sum(w * (1 - h)))
    df = trA ** 2 / float(w @ (M ** 2) @ w)
    t = B[j] / se
    return t, 2 * ss.t.sf(np.abs(t), df), df


def fit_many(ctx, Ys: pd.DataFrame, Z):
    """Group columns of Ys by missingness pattern; design rules of L.design per pattern. Returns DataFrame t, p."""
    out = {}
    miss = Ys.isna()
    keys = miss.T.apply(lambda r: hash(r.values.tobytes()), axis=1)
    for _, cols in keys.groupby(keys):
        c0 = cols.index[0]
        X, yy, ok, _, keep = (L.design(ctx, Ys[c0], Z) if Ys[c0].notna().any()
                              else (None, None, None, None, []))
        if "context" not in keep:
            for c in cols.index:
                out[c] = (np.nan, np.nan)
            continue
        n_ctx = int((X[:, keep.index("context")] != 0).sum())
        if n_ctx < 12 or len(yy) - n_ctx < 12:
            for c in cols.index:
                out[c] = (np.nan, np.nan)
            continue
        Y = Ys.loc[ok[ok].index, list(cols.index)].values.astype(float)
        t, p, _ = hc3_many(X, Y, keep.index("context"))
        for c, tt, pp in zip(cols.index, t, p):
            out[c] = (tt, pp)
    return pd.DataFrame(out, index=["t", "p"]).T


# ------------------------------------------------------------------------------------------------ stage-wise BB
def bh_level(p, level):
    return L.bh_reject(np.where(np.isfinite(p), p, 1.0), level)


def stagewise(t2, pcol, parents):
    """Registry fdr_rule: families (c, P) with P discovered for c in the parent direction (H3a lineage/melanoma,
    H3b driver) and >= 2 sub-modules; BH inside each family at q R_x / m_x of the parent direction."""
    lev = {}
    for x, tp in parents.items():
        R, m = int(tp.discovery.sum()), int(tp.evaluable.sum())
        for c in tp.context.unique():
            lev[c] = 0.05 * R / m
    disc = {x: set(zip(tp.context[tp.discovery], tp.unit[tp.discovery])) for x, tp in parents.items()}
    rej = np.zeros(len(t2), bool)
    nfam = 0
    for (c, pa), g in t2.groupby(["context", "parent"]):
        if not any((c, pa) in s for s in disc.values()):
            continue
        nfam += 1
        rej[t2.index.get_indexer(g.index)] = bh_level(g[pcol].values, lev[c])
    return rej, nfam, lev


# ------------------------------------------------------------------------------------------------ data helpers
def membership_all():
    mem = pd.read_csv(L.ROOT / "results/part2/fzx/zavg_engine_membership.tsv", sep="\t")
    m = mem[mem.role == "member"]
    members = {u: list(g.gene) for u, g in m.groupby("unit_id", sort=True)}
    level = mem.groupby("unit_id")["level"].first().to_dict()
    parent = mem.groupby("unit_id")["parent"].first().fillna("").astype(str).to_dict()
    return members, level, parent


def pair_scores():
    d = np.load(L.ROOT / "data_cache/part2/blocks/zavg_BEF_bd28f78d612bef9f_S8c264ca513f8.npz")
    genes = pd.read_csv(PART1 / "outputs/all_pairs/genes.tsv", sep="\t")
    keys = genes.depmap_column.values[d["s_index"]]
    return pd.DataFrame(d["vals"], index=keys, columns=keys)


def cov_source_plus(models):
    src = L.cov_source(models)
    cf = pd.read_csv(L.ROOT / "datasets/DepMAP/CRISPRConfounders.csv", index_col=0).reindex(models)
    src["NNMD"] = cf.ScreenNNMD
    src["CasActivity"] = cf.CasActivity
    return src


# ------------------------------------------------------------------------------------------------ matched sets (mirror)
def matched_sets_mirror(E, members, coarse, PSdf, R, seed, unit_idx, tol_ess=0.1, tol_coh=1.0, max_att=20000,
                        rng_override=None):
    allg = sorted({g for u in coarse for g in members[u] if g in E.columns})
    gi = {g: i for i, g in enumerate(allg)}
    ess = E.clip(-3, 3).mean().reindex(allg).values
    PS = PSdf.reindex(index=allg, columns=allg).values.astype(float)
    pools = {u: np.array([gi[g] for g in members[u] if g in gi]) for u in coarse}
    out, status = {}, {}
    for ui in unit_idx:
        u = coarse[ui]
        mem = [g for g in members[u] if g in gi]
        k = len(mem)
        idx = np.array([gi[g] for g in mem])
        t_ess = float(np.nanmean(ess[idx]))
        iu = np.triu_indices(k, 1)
        t_coh = float(np.nanmean(PS[np.ix_(idx, idx)][iu])) if k >= 2 else np.nan
        others = [pools[v] for v in coarse if v != u and len(pools[v])]
        big = [p for p in others if len(p) >= k]
        rng = np.random.default_rng([seed, ui]) if rng_override is None else rng_override
        st, sets = "unavailable", []
        for mult in (1.0, 2.0):
            seen, sets = set(), []
            for _ in range(max_att):
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
                if abs(float(np.nanmean(ess[s])) - t_ess) > tol_ess * mult:
                    continue
                if np.isfinite(t_coh) and not abs(float(np.nanmean(PS[np.ix_(s, s)][iu])) - t_coh) <= tol_coh * mult:
                    continue
                seen.add(key)
                sets.append(s)
                if len(sets) >= R:
                    break
            if len(sets) >= R:
                st = "ok" if mult == 1.0 else "tolerances doubled"
                break
        status[u] = st
        if st != "unavailable":
            for r, s in enumerate(sets):
                out[f"{u}#r{r}"] = [allg[i] for i in s]
    return out, status


def independent_sets(E, members, coarse, PSdf, R, rng, unit_idx, tol_ess=0.1, tol_coh=1.0, max_att=20000):
    """Own sampler, same acceptance rule, different stream and draw order: pick a random gene pool from ANY other
    unit with >= k members (sizes weighted equally), draw k genes, accept on |d ess| and |d coh| (doubled once)."""
    allg = sorted({g for u in coarse for g in members[u] if g in E.columns})
    ess = E.clip(-3, 3).mean().reindex(allg)
    PS = PSdf.reindex(index=allg, columns=allg)
    out = {}
    for ui in unit_idx:
        u = coarse[ui]
        mem = [g for g in members[u] if g in ess.index]
        k = len(mem)
        t_ess = ess[mem].mean()
        sub = PS.loc[mem, mem].values
        t_coh = np.nanmean(sub[np.triu_indices(k, 1)])
        donors = [v for v in coarse if v != u and len([g for g in members[v] if g in ess.index]) >= k]
        if not donors:
            continue
        got = []
        for mult in (1, 2):
            got = []
            for _ in range(max_att):
                v = donors[rng.integers(len(donors))]
                pool = [g for g in members[v] if g in ess.index]
                s = sorted(rng.choice(pool, k, replace=False).tolist())
                if abs(ess[s].mean() - t_ess) > 0.1 * mult:
                    continue
                ss_ = PS.loc[s, s].values
                if abs(np.nanmean(ss_[np.triu_indices(k, 1)]) - t_coh) > 1.0 * mult:
                    continue
                if s not in got:
                    got.append(s)
                if len(got) >= R:
                    break
            if len(got) >= R:
                break
        if len(got) >= R:
            for r, s in enumerate(got):
                out[f"{u}#i{r}"] = s
    return out


def score_cols(E, sets, removed=frozenset()):
    cols = {}
    for k, g in sets.items():
        use = [x for x in g if x not in removed and x in E.columns]
        cols[k] = L.score(E, use) if len(use) >= 2 else pd.Series(np.nan, index=E.index)
    return pd.DataFrame(cols, index=E.index)


# ------------------------------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-units", type=int, default=30)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    rep = {}

    E = L.effects_A()
    models = list(E.index)
    src = cov_source_plus(models)
    ct = L.context_table().reindex(models)
    an = DR.Annot()
    sym2ent = {DR.sym(c): DR.ent(c) for c in E.columns}
    members, level, parent = membership_all()
    coarse = sorted(u for u, lv in level.items() if lv == 1)
    h3b_cov = DR.parse_registry_contexts(L.ROOT / "docs/registry/H3b_driver_family.yaml")
    h2_cov = DR.parse_registry_contexts(L.ROOT / "docs/registry/H2_submodule_divergence.yaml")
    T = {k: tests(v) for k, v in dict(H3a="H3a_lineage_atlas", H3b="H3b_driver_family", H1S="H1S_silencing_GI",
                                       H2="H2_submodule_divergence").items()}
    rep["n_models"] = len(models)

    # ================================================================ 1a. H1S: |N'| and 5 tests
    ev = pd.read_csv(L.ROOT / "results/part2/fzx/h1_events.tsv", sep="\t").set_index("event")
    split = pd.read_csv(L.ROOT / "results/part2/fzx/h1_event_split_A2.tsv", sep="\t").set_index("event")
    h1s_events = split.index[split.entry == "H1S_silencing_GI"]
    loss = pd.read_parquet(L.ROOT / "data_cache/part2/fzx/h1_loss_table.parquet").reindex(models)
    cls = pd.read_csv(D / "H1S_silencing_GI/classes.tsv", sep="\t").set_index("event")

    def h1_y(event, axis):
        e = ev.loc[event]
        lc = loss[f"H1loss:{e.gene}"]
        if axis == "P":
            return L.score(E, e.P.split(";")), lc, None
        N = [g for g in str(e.N).split(";") if g in E.columns]
        b = E.loc[lc == 0, N].median(skipna=True)
        Np = [g for g in N if abs(b[g]) >= 0.1]
        if len(Np) == 0:
            return None, lc, 0
        V = E[Np].clip(-3, 3) * np.sign(b[Np])
        y = V.mean(axis=1, skipna=True)
        y[V.notna().sum(axis=1) < 0.5 * len(Np)] = np.nan
        return y, lc, len(Np)

    def h1_fit(event, axis):
        e = ev.loc[event]
        y, lc, nprime = h1_y(event, axis)
        covs = [c for c in e.covariates.split(";") if c and not (e.gene == "TP53" and c == "TP53_altered")]
        f = L.fit_test(lc, y, L.cov_frame(covs, src))
        f["n_Nprime"] = nprime
        return f

    nN = {}
    for event in h1s_events:
        e = ev.loc[event]
        if int(e.n_N) >= 2:
            _, _, npr = h1_y(event, "N")
            nN[event] = npr
    nN = pd.Series(nN)
    committed_np = cls.n_Nprime.reindex(nN.index)
    t1s = T["H1S"]
    rep["H1S_evaluability"] = dict(
        events=int(len(h1s_events)), N_axes_registered=int(len(nN)),
        P_axes_registered=int((ev.loc[h1s_events].n_P >= 1).sum()),
        N_axes_Nprime_lt2_mine=int((nN < 2).sum()),
        Nprime_equal_committed=int((nN.values == committed_np.values).sum()),
        Nprime_compared=int(committed_np.notna().sum()),
        tests_built=int(len(t1s)), evaluable=int(t1s.evaluable.sum()),
        not_evaluable_reasons=t1s[~t1s.evaluable].reason.value_counts().to_dict(),
        units_missing_are_Nprime_lt2=bool(set(f"{e}|N" for e in nN.index[nN < 2]) ==
                                         set(f"{e}|N" for e in nN.index) - set(t1s.unit)))
    d = t1s[t1s.evaluable & t1s.discovery]
    nd = t1s[t1s.evaluable & ~t1s.discovery]
    pick = list(rng.choice(d.index, 3, replace=False)) + list(rng.choice(nd.index, 2, replace=False))
    # make sure one P axis is included if any discovery is a P axis
    rows = []
    for i in pick:
        r = t1s.loc[i]
        event, axis = r.unit.rsplit("|", 1)
        f = h1_fit(event, axis)
        rows.append(dict(direction="H1S", unit=r.unit, discovery=bool(r.discovery), n=r.n, n_mine=f["n"],
                         beta=r.beta, beta_mine=f["beta"], se=r.se, se_mine=f["se"], df=r.df, df_mine=f["df"],
                         p=r.p, p_mine=f["p"]))

    # ================================================================ 1b. H2: 5 tests
    t2 = T["H2"]
    subs_all = sorted(u for u, lv in level.items() if lv == 2)
    kids = pd.Series({s: parent[s] for s in subs_all}).value_counts()
    subs = sorted(s for s in subs_all if kids[parent[s]] >= 2)
    rep["H2_units"] = dict(n_submodules_eligible_parents=len(subs), units_in_tests=int(t2.unit.nunique()))

    def h2_removed(c, pmem):
        if c.startswith("lineage:") or c == "melanoma":
            return set()
        return DR.h3b_removed(c, pmem, an, sym2ent)

    def h2_y(c, s, Eu=E):
        p = parent[s]
        rem = h2_removed(c, members[p])
        sm = [g for g in members[s] if g not in rem and g in Eu.columns]
        sib = [g for g in members[p] if g not in set(members[s]) and g not in rem and g in Eu.columns]
        if len(sm) < 2 or len(sib) < 2:
            return pd.Series(np.nan, index=Eu.index)
        Ss, Sb = L.score(Eu, sm), L.score(Eu, sib)
        ok = Ss.notna() & Sb.notna()
        X = np.c_[np.ones(ok.sum()), Sb[ok].values]
        a_, b_ = np.linalg.lstsq(X, Ss[ok].values, rcond=None)[0]
        return Ss - (a_ + b_ * Sb)

    def h2_Z(c):
        return L.cov_frame(h2_cov[c], src)

    d2 = t2[t2.discovery]
    nd2 = t2[~t2.discovery]
    pick2 = list(rng.choice(d2.index, 3, replace=False)) + list(rng.choice(nd2.index, 2, replace=False))
    for i in pick2:
        r = t2.loc[i]
        f = L.fit_test(ct[r.context], h2_y(r.context, r.unit), h2_Z(r.context))
        rows.append(dict(direction="H2", unit=f"{r.context}|{r.unit}", discovery=bool(r.discovery), n=r.n,
                         n_mine=f["n"], beta=r.beta, beta_mine=f["beta"], se=r.se, se_mine=f["se"], df=r.df,
                         df_mine=f["df"], p=r.p, p_mine=f["p"]))
    rc = pd.DataFrame(rows)
    for c in ("beta", "se", "df", "p"):
        rc[f"rel_{c}"] = (rc[f"{c}_mine"] - rc[c]).abs() / rc[c].abs().clip(lower=1e-300)
    rc["match"] = (rc.n == rc.n_mine) & (rc[["rel_beta", "rel_se", "rel_df"]].max(axis=1) < 1e-4) & (rc.rel_p < 1e-3)
    rc.to_csv(out / "recompute_h1s_h2.tsv", sep="\t", index=False, float_format="%.6g")
    rep["recompute"] = dict(n=len(rc), n_match=int(rc.match.sum()),
                            max_rel={c: float(rc[f"rel_{c}"].max()) for c in ("beta", "se", "df", "p")})
    print(rc[["direction", "unit", "discovery", "n", "n_mine", "beta", "beta_mine", "p", "p_mine", "match"]]
          .to_string())

    # ================================================================ 2. H2 stage-wise; GPD p_cons H1S / H2
    parents = {"H3a": T["H3a"], "H3b": T["H3b"]}
    for tp in parents.values():
        tp["evaluable"] = tp.evaluable.astype(bool)
    rej, nfam, lev = stagewise(t2, "p_final", parents)
    rej_a, _, _ = stagewise(t2, "p", parents)
    t2["p_cons"] = np.fmax(t2.p_final, t2.p)
    rej_c, _, _ = stagewise(t2, "p_cons", parents)
    rep["H2_stagewise"] = dict(families_with_tests=nfam, levels=sorted(set(round(v, 7) for v in lev.values())),
                               discoveries_mine=int(rej.sum()), equal_committed=bool((rej == t2.discovery).all()),
                               analytic_rejections_mine=int(rej_a.sum()),
                               analytic_equal_committed=bool((rej_a == t2.reject).all()),
                               p_cons_discoveries=int(rej_c.sum()), p_cons_lost=int((t2.discovery & ~rej_c).sum()))
    e1 = t1s[t1s.evaluable]
    r1 = L.bh_reject(e1.p_final.values, 0.05)
    r1c = L.bh_reject(np.fmax(e1.p_final.values, e1.p.values), 0.05)
    r1pad = L.bh_reject(np.r_[e1.p_final.values, np.ones(338 - len(e1))], 0.05)[:len(e1)]
    rep["H1S_bh"] = dict(m=len(e1), discoveries_mine=int(r1.sum()), equal_committed=bool((r1 == e1.discovery).all()),
                         p_cons_discoveries=int(r1c.sum()), lost_p_cons=list(e1.unit[e1.discovery.values & ~r1c]),
                         at_registered_m338=int(r1pad.sum()))

    # ================================================================ 3. H3a matched random sets
    PSdf = pair_scores()
    cal = pd.read_csv(D / "H3a_lineage_atlas/calibration_D8_tests.tsv.gz", sep="\t")
    mrs = cal[cal.calibration == "matched_random_set"].copy()
    h3a = T["H3a"]
    # (c) full count from committed tables
    tmax = mrs.assign(at=mrs.t.abs()).groupby(["context", "unit"]).at.max()
    dd = h3a[h3a.discovery].copy()
    dd["tmax"] = [tmax.get((c, u), np.nan) for c, u in zip(dd.context, dd.unit)]
    rep["H3a_exceeds_all_committed_tables"] = dict(n=int((dd.t.abs() > dd.tmax).sum()), of=len(dd),
                                                    without_sets=int(dd.tmax.isna().sum()),
                                                    committed_rate_p005=float((mrs.p < 0.05).mean()),
                                                    committed_n=int(mrs.p.notna().sum()))
    # (a) mirrored subsample
    uidx = sorted(rng.choice(len(coarse), a.n_units, replace=False).tolist())
    msets, mstat = matched_sets_mirror(E, members, coarse, PSdf, 10, 4600, uidx)
    ctxs = list(h3a.context.unique())
    Z = L.cov_frame(H3A_COV, src)
    Ys = score_cols(E, msets)
    mine = []
    for c in ctxs:
        f = fit_many(ct[c], Ys, Z)
        f["context"] = c
        mine.append(f.rename_axis("set").reset_index())
    mine = pd.concat(mine)
    cm = mrs.merge(mine, on=["context", "set"], how="inner", suffixes=("", "_mine"))
    sub_units = [coarse[i] for i in uidx]
    exp_sets = mrs[mrs.unit.isin(sub_units)]
    rep["H3a_matched_mirror"] = dict(
        units=len(uidx), status_mine=pd.Series(mstat).value_counts().to_dict(), sets_mine=len(msets),
        committed_sets_for_units=int(exp_sets.set.nunique()), tests_matched=len(cm), committed_tests=len(exp_sets),
        max_abs_dt=float((cm.t - cm.t_mine).abs().max()) if len(cm) else None,
        max_rel_p=float(((cm.p - cm.p_mine).abs() / cm.p.clip(lower=1e-300)).max()) if len(cm) else None,
        rate_p005_mine=float((mine.p < 0.05).sum() / mine.p.notna().sum()),
        rate_p005_committed_same_units=float((exp_sets.p < 0.05).sum() / exp_sets.p.notna().sum()),
        rate_p005_real_same_units=float((h3a[h3a.unit.isin(sub_units) & h3a.evaluable].p < 0.05).mean()))
    tm = mine.assign(unit=mine.set.str.split("#").str[0], at=mine.t.abs()).groupby(["context", "unit"]).at.max()
    ds = dd[dd.unit.isin(sub_units)]
    ex_mine = [abs(r.t) > tm.get((r.context, r.unit), np.inf) for r in ds.itertuples()]
    ex_comm = [abs(r.t) > r.tmax if np.isfinite(r.tmax) else False for r in ds.itertuples()]
    rep["H3a_matched_mirror"].update(discoveries_in_units=len(ds), exceed_all_mine=int(sum(ex_mine)),
                                     exceed_all_committed=int(sum(ex_comm)))
    # (b) independent re-draw on the same units
    isets = independent_sets(E, members, coarse, PSdf, 10, np.random.default_rng(SEED + 1), uidx)
    Yi = score_cols(E, isets)
    ind = []
    for c in ctxs:
        f = fit_many(ct[c], Yi, Z)
        f["context"] = c
        ind.append(f.rename_axis("set").reset_index())
    ind = pd.concat(ind)
    ti = ind.assign(unit=ind.set.str.split("#").str[0], at=ind.t.abs()).groupby(["context", "unit"]).at.max()
    ex_ind = [abs(r.t) > ti.get((r.context, r.unit), np.inf) for r in ds.itertuples()]
    rep["H3a_matched_independent"] = dict(sets=len(isets), rate_p005=float((ind.p < 0.05).sum() / ind.p.notna().sum()),
                                          exceed_all=int(sum(ex_ind)), discoveries_in_units=len(ds))
    print(json.dumps({k: rep[k] for k in rep if k.startswith("H3a") or k.startswith("H2_") or k.startswith("H1S")},
                     indent=1, default=str))

    # ================================================================ 4. never-expressed
    ex = pd.read_csv(L.ROOT / "datasets/DepMAP/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv", index_col=0,
                     low_memory=False)
    ex = ex[ex.IsDefaultEntryForModel == "Yes"].set_index("ModelID")
    exA = ex.reindex([m for m in models if m in ex.index])
    vals = [c for c in exA.columns if "(" in c]
    exA = exA[vals].astype(float).dropna(how="all")
    mx = exA.max(axis=0, skipna=True)
    never = set(mx.index[mx < 1.0])
    pool = sorted(g for g in E.columns if g in never and an.chrom.get(DR.ent(g)) and np.isfinite(an.mid.get(DR.ent(g), np.nan)))
    rep["never_expressed"] = dict(rna_profiled_A=int(len(exA)), pool_mine=len(pool))
    # mirrored draw
    r_ = np.random.default_rng([4600, 1_000_000])
    sizes = r_.choice([len(members[u]) for u in coarse], 100)

    def near(g, h):
        eg, eh = DR.ent(g), DR.ent(h)
        return an.chrom.get(eg) == an.chrom.get(eh) and abs(an.mid[eg] - an.mid[eh]) <= 10e6

    def dispersed(pool, sizes, r_):
        sets = {}
        for i, k in enumerate(sizes):
            ch = []
            for j in r_.permutation(len(pool)):
                g = pool[j]
                if any(near(g, h) for h in ch):
                    continue
                ch.append(g)
                if len(ch) == k:
                    break
            sets[f"never{i}_k{k}"] = ch
        return sets
    nsets = dispersed(pool, sizes, r_)
    nsets_ind = dispersed(pool, np.random.default_rng(SEED + 2).choice([len(members[u]) for u in coarse], 100),
                          np.random.default_rng(SEED + 3))
    calN = cal[cal.calibration == "never_expressed"]
    res = {}
    for name, ss_ in (("mirror", nsets), ("independent", nsets_ind)):
        Yn = score_cols(E, ss_)
        allp, per = [], {}
        for c in ctxs:
            f = fit_many(ct[c], Yn, Z)
            allp.append(f.p.values)
            per[c] = float(np.nanmean(f.p.values < 0.05))
            if name == "mirror":
                f["context"] = c
                res.setdefault("mirror_tab", []).append(f.rename_axis("set").reset_index())
        p = np.concatenate(allp)
        p = p[np.isfinite(p)]
        rb, Rb, _ = L.bb(p, np.repeat(ctxs, [len(x[np.isfinite(x)]) for x in allp]), 0.05)
        res[name] = dict(n_tests=len(p), rate_p005=float((p < 0.05).mean()), bb_rejections=int(rb.sum()),
                         top_contexts=dict(sorted(per.items(), key=lambda kv: -kv[1])[:6]))
    mt = pd.concat(res.pop("mirror_tab"))
    cmN = calN.merge(mt, on=["context", "set"], how="inner", suffixes=("", "_mine"))
    rep["never_expressed"].update(res, sizes_equal=None, committed_rate=float((calN.p < 0.05).mean()),
                                  committed_n=len(calN), mirror_tests_matched=len(cmN),
                                  mirror_max_abs_dt=float((cmN.t - cmN.t_mine).abs().max()) if len(cmN) else None,
                                  committed_per_context_top=calN.groupby("context").p.apply(
                                      lambda s: float((s < 0.05).mean())).sort_values(ascending=False).head(6).to_dict())
    print(json.dumps(rep["never_expressed"], indent=1, default=str))

    # ================================================================ 5. H2 random splits
    ctx2 = [c for c in h2_cov]                                  # registry order
    calH2 = pd.read_csv(D / "H2_submodule_divergence/calibration_D8_tests.tsv.gz", sep="\t")
    rs = calH2[calH2.calibration == "random_split"]
    sel = {}
    for x, tp in parents.items():
        for c, u in zip(tp.context[tp.discovery], tp.unit[tp.discovery]):
            sel[(c, u)] = True
    split_res = {}
    for mode in ("mirror", "independent"):
        recs = []
        for r in t2.itertuples(index=False):
            p_ = parent[r.unit]
            rem = h2_removed(r.context, members[p_])
            pm = [g for g in members[p_] if g not in rem and g in E.columns]
            k = len([g for g in members[r.unit] if g in pm])
            rg = (np.random.default_rng([4800, subs.index(r.unit), ctx2.index(r.context)]) if mode == "mirror"
                  else np.random.default_rng([SEED + 4, subs.index(r.unit), ctx2.index(r.context)]))
            cols = {}
            for j in range(20):
                s = sorted(rg.choice(pm, k, replace=False).tolist()) if 2 <= k < len(pm) - 1 else []
                if not s:
                    cols[f"{r.unit}#split{j}"] = pd.Series(np.nan, index=E.index)
                    continue
                sib = [g for g in pm if g not in set(s)]
                Ss, Sb = L.score(E, s), L.score(E, sib)
                ok = Ss.notna() & Sb.notna()
                X = np.c_[np.ones(ok.sum()), Sb[ok].values]
                a_, b_ = np.linalg.lstsq(X, Ss[ok].values, rcond=None)[0]
                cols[f"{r.unit}#split{j}"] = Ss - (a_ + b_ * Sb)
            f = fit_many(ct[r.context], pd.DataFrame(cols, index=E.index), h2_Z(r.context))
            for j, (nm, v) in enumerate(f.iterrows()):
                recs.append(dict(rep=j, context=r.context, parent=r.parent, set=nm, unit=r.unit, t=v.t, p=v.p))
        d_ = pd.DataFrame(recs)
        nrej = []
        for j, g in d_.groupby("rep"):
            g = g.reset_index(drop=True)
            rj, _, _ = stagewise(g.assign(discovery=False), "p", parents)
            nrej.append(int(rj.sum()))
        split_res[mode] = dict(rate_p005=float((d_.p < 0.05).sum() / d_.p.notna().sum()), n_tests=int(d_.p.notna().sum()),
                               stagewise_mean=float(np.mean(nrej)), stagewise_max=int(max(nrej)),
                               stagewise_min=int(min(nrej)))
        if mode == "mirror":
            cmS = rs.merge(d_, on=["context", "set"], how="inner", suffixes=("", "_mine"))
            split_res["mirror"].update(tests_matched=len(cmS), committed=len(rs),
                                       max_abs_dt=float((cmS.t - cmS.t_mine).abs().max()))
    rep["H2_random_splits"] = split_res
    rep["H2_random_splits"]["committed_rate_p005"] = float((rs.p < 0.05).mean())
    print(json.dumps(rep["H2_random_splits"], indent=1))
    (out / "review2.json").write_text(json.dumps(rep, indent=1, default=lambda o: o.item() if hasattr(o, "item")
                                                 else str(o)) + "\n")


if __name__ == "__main__":
    main()
