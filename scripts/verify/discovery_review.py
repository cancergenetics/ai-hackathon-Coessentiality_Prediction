"""Independent review of the A-side discovery runs H3a / H3b / H1 (single look). Does NOT import src/part2.

    .venv/bin/python scripts/verify/discovery_review.py --out <scratch>/discovery_review      # ~3-5 min

1. Recompute analytic beta / SE / df / p for 30 tests: 10 random H3a discoveries, 10 random H3b discoveries,
   5 random non-discoveries (H3a + H3b pooled), and the 5 H1 discovered tests (kappa_N or beta_P per score_rule,
   ALIGN = +1: y_N = mean over N' of sign(b_j) * effect_j).
2. BB (Simes family selection; within-family BH at q R / M) from the committed p_final; exact discovery sets.
   H1: pooled BH at m = 24 (and at the registered m = 28, padding).
3. TOST neutral calls for 3 random NEUTRAL H1 events.
4. GPD tail: p = max(p_perm, p_analytic) and p = max(p_perm, 1 / 10,001) for screened tests; BB re-run.
5. Sanity: H3a discoveries per lineage vs n, effect sizes, locus share; leave-top-k-models-out for all H3a/H3b
   discoveries (k = 1, 2, 3, 5; most supportive models by exact dfbeta), detail for the 5 top H3a hits.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as ss

sys.path.insert(0, str(Path(__file__).resolve().parent))
import discovery_lib as L  # noqa: E402

SEED = 20261029
DIRS = {"H3a": "H3a_lineage_atlas", "H3b": "H3b_driver_family", "H1": "H1_within_complex_GI"}
H3A_COV = ["library", "ScreenType3DO", "ModelType"]


def tests(d):
    return pd.read_csv(L.ROOT / f"results/part2/directions/{d}/tests.tsv.gz", sep="\t")


def parse_registry_contexts(path):
    """Minimal parser: context id -> covariate list (H3b) from the registry YAML (no yaml module in the venv)."""
    out, cur = {}, None
    for line in Path(path).read_text().splitlines():
        s = line.strip()
        if s.startswith("- id: "):
            cur = s[6:].strip().strip('"')
        elif s.startswith("covariates: [") and cur is not None and line.startswith("      "):
            out[cur] = [x.strip() for x in s[len("covariates: ["):-1].split(",")]
    return out


# ------------------------------------------------------------------------------------------------ annotation
class Annot:
    """Entrez-keyed Ensembl 116 paralogs (registry paralog_definition) and coordinates."""

    def __init__(self):
        co = pd.read_csv(L.ROOT / "datasets/Paralogs/ensembl116_gene_coordinates.tsv", sep="\t", dtype=str)
        co.columns = ["ens", "sym", "entrez", "chrom", "start", "end"]
        self.ens_chrom = dict(zip(co.ens, co.chrom))
        # one location per Entrez: primary-assembly chromosomes before patches; X before Y for PAR genes
        rank = {str(i): i for i in range(1, 23)} | {"X": 23, "Y": 24, "MT": 25}
        c2 = co.dropna(subset=["entrez"]).assign(_r=lambda d: d.chrom.map(rank).fillna(99))
        c2 = c2.sort_values("_r").drop_duplicates("entrez")
        self.chrom = dict(zip(c2.entrez, c2.chrom))
        self.mid = dict(zip(c2.entrez, (c2.start.astype(float) + c2.end.astype(float)) / 2))
        pt = pd.read_csv(L.ROOT / "datasets/Paralogs/ensembl116_human_paralogues.tsv", sep="\t", dtype=str)
        pt.columns = ["g", "gs", "h", "hs", "idt", "idq", "type", "lca"]
        pt = pt[pt.type.isin(["within_species_paralog", "other_paralog"])]
        pt = pt[np.minimum(pt.idt.astype(float), pt.idq.astype(float)) >= 20]
        self.par_ens = {}
        for g, h in zip(pt.g, pt.h):
            self.par_ens.setdefault(g, set()).add(h)
            self.par_ens.setdefault(h, set()).add(g)
        self.ent2ens = {}                                     # one Ensembl ID can carry several Entrez IDs
        for e, n in zip(co.ens, co.entrez):
            if isinstance(n, str):
                self.ent2ens.setdefault(n, set()).add(e)

    def paralogs(self, ent):
        out = set()
        for e in self.ent2ens.get(ent, ()):
            out |= self.par_ens.get(e, set())
        return out                                           # Ensembl IDs

    def is_paralog(self, a, b):
        pa = self.paralogs(a)
        return any(e in pa for e in self.ent2ens.get(b, ()))

    def has_chrY_paralog(self, ent):
        return any(self.ens_chrom.get(e) == "Y" for e in self.paralogs(ent))


def ent(key):
    return key.rsplit("(", 1)[1].rstrip(")").strip()


def sym(key):
    return key.rsplit("(", 1)[0].strip()


H3B_GENES = {"TP53_altered": ["TP53"], "MSI_H": [], "WGD": [], "9p21_loss": ["MTAP", "CDKN2A", "CDKN2B"],
             "MYC_amp": ["MYC"], "MYCN_amp": ["MYCN"], "NRF2_active": ["KEAP1", "CUL3", "NFE2L2"],
             "fusion:EWSR1::FLI1|ERG": ["EWSR1", "FLI1", "ERG"], "LoY": [], "sex_female": []}
H3B_WINDOW = {"9p21_loss": "MTAP", "MYC_amp": "MYC", "MYCN_amp": "MYCN"}


def h3b_removed(ctx, members, an, sym2ent):
    """Registry removed_members for an H3b context, restricted to the unit's members."""
    genes = H3B_GENES.get(ctx, [ctx.split(":", 1)[1]] if ":" in ctx else [])
    gent = [sym2ent[g] for g in genes if g in sym2ent]
    rem = set()
    for m in members:
        e = ent(m)
        if e in gent or any(an.is_paralog(g, e) for g in gent):
            rem.add(m)
        if ctx in H3B_WINDOW:
            w = sym2ent[H3B_WINDOW[ctx]]
            if an.chrom.get(e) == an.chrom.get(w) and abs(an.mid.get(e, np.inf) - an.mid[w]) <= 10e6:
                rem.add(m)
        if ctx in ("LoY", "sex_female"):
            if an.chrom.get(e) == "Y" or (an.chrom.get(e) == "X" and an.has_chrY_paralog(e)):
                rem.add(m)
    return rem


# ------------------------------------------------------------------------------------------------ influence
def leave_top_k(ctx, y, Z, ks=(1, 2, 3, 5)):
    X, yy, ok, _, keep = L.design(ctx, y, Z)
    j = keep.index("context")
    XtXi = np.linalg.inv(X.T @ X)
    G = XtXi @ X.T
    h = np.einsum("ij,ji->i", X, G)
    e = yy - X @ (G @ yy)
    beta = G[j] @ yy
    dfb = G[j] * e / (1 - h)                                  # beta - beta_(-i)
    order = np.argsort(-np.sign(beta) * dfb)                  # most supportive first
    ids = np.asarray(ok[ok].index)
    out = {}
    for k in ks:
        yk = y.copy()
        yk[ids[order[:k]]] = np.nan                           # design rules re-applied on the reduced rows
        r = L.fit_test(ctx, yk, Z)
        out[k] = dict(beta=r["beta"], p=r["p"], removed=list(ids[order[:k]]),
                      removed_in_context=int((X[order[:k], j] != 0).sum()))
    return out


# ------------------------------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    rep = {}

    E = L.effects_A()
    A, Lm = L.a_models(), L.l_models()
    rep["A_n"], rep["L_n"], rep["A_and_L"] = len(A), len(Lm), len(set(A) & set(Lm))
    rep["E_rows_equal_A"] = sorted(E.index) == sorted(A)
    units, mem, locus = L.membership()
    src = L.cov_source(list(E.index))
    ct = L.context_table().reindex(E.index)
    an = Annot()
    sym2ent = {sym(c): ent(c) for c in E.columns}
    h3b_cov = parse_registry_contexts(L.ROOT / "docs/registry/H3b_driver_family.yaml")
    T = {k: tests(v) for k, v in DIRS.items()}

    # n_used never exceeds A complete-case count; all n <= |A|
    rep["max_n_used"] = {k: int(t.n.max()) for k, t in T.items()}

    # ---------------------------------------------------------------- 1. recompute sample
    def h3_fit(kind, ctx, unit):
        members = mem[unit]
        if kind == "H3a":
            covs, rem = H3A_COV, set()
        else:
            covs, rem = h3b_cov[ctx], h3b_removed(ctx, members, an, sym2ent)
        kept = [m for m in members if m not in rem]
        y = L.score(E, kept)
        return L.fit_test(ct[ctx], y, L.cov_frame(covs, src)), y, covs, sorted(rem)

    rows = []
    picks = []
    for k in ("H3a", "H3b"):
        d = T[k][T[k].discovery]
        picks += [(k, "discovery", i) for i in rng.choice(d.index.values, 10, replace=False)]
    pool = [(k, i) for k in ("H3a", "H3b") for i in T[k].index[~T[k].discovery]]
    for j in rng.choice(len(pool), 5, replace=False):
        picks.append((pool[j][0], "non-discovery", pool[j][1]))
    for k, lab, i in picks:
        r = T[k].loc[i]
        f, _, covs, rem = h3_fit(k, r.context, r.unit)
        rows.append(dict(direction=k, sample=lab, context=r.context, unit=r.unit, n=r.n, n_mine=f["n"],
                         beta=r.beta, beta_mine=f["beta"], se=r.se, se_mine=f["se"], df=r.df, df_mine=f["df"],
                         p=r.p, p_mine=f["p"], removed_mine=";".join(rem)))

    # H1
    ev = pd.read_csv(L.ROOT / "results/part2/fzx/h1_events.tsv", sep="\t").set_index("event")
    loss = pd.read_parquet(L.ROOT / "data_cache/part2/fzx/h1_loss_table.parquet").reindex(E.index)
    cls = pd.read_csv(L.ROOT / "results/part2/directions/H1_within_complex_GI/classes.tsv", sep="\t").set_index("event")

    def h1_fit(event, axis):
        e = ev.loc[event]
        lc = loss[f"H1loss:{e.gene}"]
        covs = [c for c in e.covariates.split(";") if c and not (e.gene == "TP53" and c == "TP53_altered")]
        if axis == "P":
            y = L.score(E, e.P.split(";"))
            nprime = None
        else:
            N = [g for g in e.N.split(";") if g in E.columns]
            ref = E.loc[lc == 0, N]
            b = ref.median(skipna=True)
            Np = [g for g in N if abs(b[g]) >= 0.1]
            V = E[Np].clip(-L.CLIP, L.CLIP) * np.sign(b[Np])        # ALIGN = +1
            y = V.mean(axis=1, skipna=True)
            y[V.notna().sum(axis=1) < 0.5 * len(Np)] = np.nan
            nprime = len(Np)
        f = L.fit_test(lc, y, L.cov_frame(covs, src))
        f["n_Nprime"] = nprime
        return f

    h1t = T["H1"]
    for i, r in h1t[h1t.discovery].iterrows():
        event, axis = r.unit.rsplit("|", 1)
        f = h1_fit(event, axis)
        rows.append(dict(direction="H1", sample=f"discovery ({'kappa_N' if axis == 'N' else 'beta_P'})",
                         context=r.context, unit=r.unit, n=r.n, n_mine=f["n"], beta=r.beta, beta_mine=f["beta"],
                         se=r.se, se_mine=f["se"], df=r.df, df_mine=f["df"], p=r.p, p_mine=f["p"],
                         removed_mine=f"n_Nprime={f['n_Nprime']} (committed {cls.loc[event, 'n_Nprime']})"
                         if axis == "N" else ""))
    rc = pd.DataFrame(rows)
    for c in ("beta", "se", "df", "p"):
        rc[f"rel_{c}"] = (rc[f"{c}_mine"] - rc[c]).abs() / rc[c].abs().clip(lower=1e-300)
    rc["match"] = (rc.n == rc.n_mine) & (rc[["rel_beta", "rel_se", "rel_df"]].max(axis=1) < 1e-4) & (rc.rel_p < 1e-3)
    rc.to_csv(out / "recompute_30.tsv", sep="\t", index=False, float_format="%.6g")
    rep["recompute"] = dict(n=len(rc), n_match=int(rc.match.sum()), max_rel_beta=float(rc.rel_beta.max()),
                            max_rel_se=float(rc.rel_se.max()), max_rel_df=float(rc.rel_df.max()),
                            max_rel_p=float(rc.rel_p.max()))

    # ---------------------------------------------------------------- 2. FDR from committed p_final
    fdr = {}
    for k in ("H3a", "H3b"):
        t = T[k]
        rej_a, Ra, _ = L.bb(t.p.values, t.context.values, 0.05)
        scr, _, _ = L.bb(t.p.values, t.context.values, 0.1)
        rej, R, sp = L.bb(t.p_final.values, t.context.values, 0.05)
        fdr[k] = dict(M=int(t.context.nunique()), m=len(t), R_analytic=Ra, reject_analytic=int(rej_a.sum()),
                      reject_analytic_equal=bool((rej_a == t.reject.values).all()),
                      screen_bb_q01=int(scr.sum()), screen_equal=bool((scr == t.p_perm.notna().values).all()),
                      R_final=R, discoveries=int(rej.sum()), discoveries_equal=bool((rej == t.discovery.values).all()),
                      families_not_selected=[f for f in sp if f not in set(t.context[rej])])
    t = T["H1"]
    ev_t = t[t.evaluable]
    rej = L.bh_reject(ev_t.p_final.values, 0.05)
    rej28 = L.bh_reject(np.r_[ev_t.p_final.values, np.ones(28 - len(ev_t))], 0.05)[: len(ev_t)]
    rej_an = L.bh_reject(ev_t.p.values, 0.05)
    fdr["H1"] = dict(m=len(ev_t), discoveries=int(rej.sum()), discoveries_equal=bool((rej == ev_t.discovery.values).all()),
                     discoveries_at_registered_m28=int(rej28.sum()), analytic_reject=int(rej_an.sum()),
                     disc_units=list(ev_t.unit[rej]),
                     disc_p_final=dict(zip(ev_t.unit[rej], ev_t.p_final[rej])),
                     disc_p_analytic=dict(zip(ev_t.unit[rej], ev_t.p[rej])))
    rep["fdr"] = fdr

    # ---------------------------------------------------------------- 3. TOST for 3 neutral H1 events
    neutral = cls.index[cls.class_complex_wide == "NEUTRAL"].tolist()
    tost = []
    for event in rng.choice(neutral, 3, replace=False):
        f = h1_fit(event, "N")
        tc = ss.t.ppf(0.95, f["df"])
        lo, hi = f["beta"] - tc * f["se"], f["beta"] + tc * f["se"]
        p_lo = ss.t.sf((f["beta"] + 0.1) / f["se"], f["df"])
        p_hi = ss.t.cdf((f["beta"] - 0.1) / f["se"], f["df"])
        c = cls.loc[event]
        tost.append(dict(event=event, kappa=f["beta"], kappa_committed=c.kappa_N, ci90=[lo, hi],
                         ci90_committed=[c.kappa_N_ci90_lo, c.kappa_N_ci90_hi], tost_p=max(p_lo, p_hi),
                         neutral_mine=bool(lo > -0.1 and hi < 0.1), neutral_committed=bool(c.N_tost),
                         sig=bool(f["p"] < 0.05)))
    rep["tost"] = tost

    # ---------------------------------------------------------------- 4. GPD sensitivity
    gpd = {}
    for k in ("H3a", "H3b"):
        t = T[k]
        s = t.p_perm.notna()
        base = t.discovery.values
        res = {}
        for name, pf in (("max_perm_analytic", np.where(s, np.maximum(t.p_perm, t.p), t.p)),
                         ("floor_1_over_10001", np.where(s, np.maximum(t.p_perm, 1 / 10001), t.p))):
            rej, R, _ = L.bb(pf, t.context.values, 0.05)
            res[name] = dict(discoveries=int(rej.sum()), lost=int((base & ~rej).sum()), gained=int((~base & rej).sum()),
                             R=R)
        res["n_screened"] = int(s.sum())
        res["n_pperm_lt_analytic"] = int((t.p_perm < t.p).sum())
        res["n_pperm_lt_1e-4"] = int((t.p_perm < 1e-4).sum())
        res["n_pperm_below_analytic_by_10x"] = int((t.p_perm < t.p / 10).sum())
        res["n_pperm_above_analytic"] = int((t.p_perm > t.p).sum())
        gpd[k] = res
    t = T["H1"][T["H1"].evaluable]
    for name, pf in (("max_perm_analytic", np.where(t.p_perm.notna(), np.maximum(t.p_perm, t.p), t.p)),
                     ("floor_1_over_10001", np.where(t.p_perm.notna(), np.maximum(t.p_perm, 1 / 10001), t.p))):
        rej = L.bh_reject(pf, 0.05)
        gpd.setdefault("H1", {})[name] = dict(discoveries=int(rej.sum()), lost=list(t.unit[t.discovery.values & ~rej]))
    rep["gpd"] = gpd

    # ---------------------------------------------------------------- 5. sanity
    h3a = T["H3a"]
    d = h3a[h3a.discovery]
    nA = h3a.groupby("context").n_context.first()
    lin = pd.DataFrame(dict(n_context=nA, discoveries=d.groupby("context").size(),
                            pos=d[d.beta > 0].groupby("context").size(), neg=d[d.beta < 0].groupby("context").size(),
                            med_abs_beta=d.groupby("context").beta.apply(lambda b: b.abs().median()),
                            min_abs_beta=d.groupby("context").beta.apply(lambda b: b.abs().min()),
                            locus=d[d.locus].groupby("context").size())).fillna(0)
    lin["frac_of_367"] = lin.discoveries / 367
    ess = pd.Series({u: L.score(E, mem[u]).mean() for u in units})   # mean S_C over A (scale check)
    lin["rho_beta_vs_meanSC"] = [ss.spearmanr(g.beta.values, ess[g.unit].values).statistic
                                 for c in lin.index for g in [h3a[h3a.context == c]]]
    lin["disc_lt005_abs_beta"] = [int((d[d.context == c].beta.abs() < 0.05).sum()) for c in lin.index]
    lin = lin.sort_values("discoveries", ascending=False)
    lin.to_csv(out / "h3a_per_lineage.tsv", sep="\t", float_format="%.4g")
    rho = ss.spearmanr(lin.n_context, lin.discoveries)
    rep["h3a_lineage"] = dict(spearman_n_vs_disc=[float(rho.statistic), float(rho.pvalue)],
                              n_disc=int(len(d)), n_locus=int(d.locus.sum()), frac_locus=float(d.locus.mean()),
                              n_units_with_any_disc=int(d.unit.nunique()),
                              units_locus_total=int(sum(bool(locus.get(u)) for u in units)))
    h3b = T["H3b"]
    db = h3b[h3b.discovery]
    rep["h3b_locus"] = dict(n_disc=int(len(db)), n_locus=int(db.locus.sum()), frac_locus=float(db.locus.mean()))

    infl = []
    for k, t in (("H3a", h3a), ("H3b", h3b)):
        for i, r in t[t.discovery].iterrows():
            members = mem[r.unit]
            covs = H3A_COV if k == "H3a" else h3b_cov[r.context]
            rem = set() if k == "H3a" else h3b_removed(r.context, members, an, sym2ent)
            y = L.score(E, [m for m in members if m not in rem])
            lt = leave_top_k(ct[r.context], y, L.cov_frame(covs, src))
            row = dict(direction=k, context=r.context, unit=r.unit, n_context=r.n_context, beta=r.beta, p=r.p,
                       p_final=r.p_final, locus=r.locus)
            for kk, v in lt.items():
                row[f"p_drop{kk}"] = v["p"]
                row[f"beta_drop{kk}"] = v["beta"]
                row[f"ctx_in_drop{kk}"] = v["removed_in_context"]
            row["drop5_models"] = ";".join(lt[5]["removed"])
            infl.append(row)
    infl = pd.DataFrame(infl)
    infl.to_csv(out / "influence_all_discoveries.tsv", sep="\t", index=False, float_format="%.4g")
    summ = {}
    for k, g in infl.groupby("direction"):
        summ[k] = {f"drop{kk}": dict(p_ge_005=int((g[f"p_drop{kk}"] >= 0.05).sum()),
                                     sign_flip=int((np.sign(g[f"beta_drop{kk}"]) != np.sign(g.beta)).sum()),
                                     p_ge_1e3=int((g[f"p_drop{kk}"] >= 1e-3).sum()))
                   for kk in (1, 2, 3, 5)}
        summ[k]["n"] = len(g)
        fragile = g[g.p_drop2 >= 0.05].sort_values("p_drop2", ascending=False)
        summ[k]["fragile_drop2"] = fragile[["context", "unit", "n_context", "beta", "p", "p_drop2"]].head(15)\
            .to_dict("records")
    rep["influence"] = summ
    h1inf = []
    for i, r in h1t[h1t.discovery].iterrows():
        event, axis = r.unit.rsplit("|", 1)
        e = ev.loc[event]
        lc = loss[f"H1loss:{e.gene}"]
        covs = [c for c in e.covariates.split(";") if c and not (e.gene == "TP53" and c == "TP53_altered")]
        if axis == "P":
            y = L.score(E, e.P.split(";"))
        else:
            N = [g for g in e.N.split(";") if g in E.columns]
            b = E.loc[lc == 0, N].median(skipna=True)
            Np = [g for g in N if abs(b[g]) >= 0.1]
            V = E[Np].clip(-L.CLIP, L.CLIP) * np.sign(b[Np])
            y = V.mean(axis=1, skipna=True)
            y[V.notna().sum(axis=1) < 0.5 * len(Np)] = np.nan
        lt = leave_top_k(lc, y, L.cov_frame(covs, src), ks=(1, 2, 3))
        h1inf.append(dict(unit=r.unit, beta=r.beta, p=r.p, **{f"p_drop{k}": v["p"] for k, v in lt.items()},
                          **{f"beta_drop{k}": v["beta"] for k, v in lt.items()},
                          **{f"loss_models_in_drop{k}": v["removed_in_context"] for k, v in lt.items()}))
    rep["influence_H1"] = h1inf
    top = d.assign(_t=-d.t.abs()).sort_values(["p_final", "_t"]).head(5)
    rep["top5_h3a"] = infl[infl.direction == "H3a"].set_index(["context", "unit"]).loc[
        list(zip(top.context, top.unit))].reset_index()[["context", "unit", "n_context", "beta", "p", "p_final",
                                                           "beta_drop1", "p_drop1", "beta_drop2", "p_drop2",
                                                           "beta_drop5", "p_drop5", "drop5_models"]].to_dict("records")
    (out / "review.json").write_text(json.dumps(rep, indent=1, default=lambda o: o.item() if hasattr(o, "item")
                                                else str(o)) + "\n")
    print(rc[["direction", "sample", "context", "unit", "n", "n_mine", "beta", "beta_mine", "se", "se_mine", "p",
              "p_mine", "match"]].to_string())
    print(json.dumps({k: v for k, v in rep.items() if k not in ("influence", "top5_h3a")}, indent=1, default=str))
    print(json.dumps(rep["influence"], indent=1, default=str))
    print(pd.DataFrame(rep["top5_h3a"]).to_string())
    print(lin.to_string())


if __name__ == "__main__":
    main()
