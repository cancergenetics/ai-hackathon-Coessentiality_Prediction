#!/usr/bin/env python
"""Independent review of the lockbox results: rates, labels, BRIDGE-L, BENCH and the power aid, from committed per-hit
files only (no L access; nothing from scripts/55_lockbox_r7_report.py or src/part2/bench.py is imported).
Output: <out>/rates_independent.tsv, defect_sensitive_independent.tsv, power_aid.tsv, tables_review.json.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sst

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "part2" / "lockbox_confirm"
DIRS = ROOT / "results" / "part2" / "directions"
FZX = ROOT / "results" / "part2" / "fzx"
H = ("H3a_lineage_atlas", "H3b_driver_family", "H1_within_complex_GI", "H1S_silencing_GI", "H2_submodule_divergence")
BR = "BRIDGE-L_confirmation_rate"
CAL = {"H3a_lineage_atlas": "matched_random_set", "H3b_driver_family": "matched_random_set",
       "H2_submodule_divergence": "random_split"}


def tb(s):
    return s.astype(str).str.lower() == "true"


def load():
    rows = []
    for p in sorted(OUT.glob("*/confirmation.tsv")):
        t = pd.read_csv(p, sep="\t")
        t["entry"], t["hit_index"] = p.parent.name, range(len(t))
        rows.append(t)
    t = pd.concat(rows, ignore_index=True)
    for c in ("eligible", "lockbox_testable", "confirmed_final", "testable_L", "same_sign"):
        t[c] = tb(t[c])
    c = pd.read_csv(OUT / "corrected_refit_nonconfirmatory.tsv", sep="\t")
    t = t.merge(c[["entry", "hit_index", "beta_L_relook", "se_L_relook", "p_L_relook", "testable_L_relook",
                   "same_sign_relook", "confirmed_relook", "affected_cause"]], on=["entry", "hit_index"], how="left")
    t["affected"] = t.affected_cause.notna()
    lb = t.eligible & t.lockbox_testable
    # primary (R7): as executed; an eligible lockbox fit without an L estimate is NOT EVALUABLE
    t["P"] = np.where(~lb, "", np.where(t.beta_L.isna(), "NE", np.where(t.confirmed_final, "C", "N")))
    # corrected: my own application of the rule to the relook statistics
    rc = tb(t.testable_L_relook) & tb(t.same_sign_relook) & (t.p_L_relook.astype(float) < 0.05)
    corr = np.where(t.beta_L_relook.isna(), "NE", np.where(rc, "C", "N"))
    t["Cc"] = np.where(~lb, "", np.where(t.affected, corr, t.P))
    t["relook_rule_mismatch"] = t.affected & (rc != tb(t.confirmed_relook))
    # primary rule re-applied to the committed L statistics
    pr = t.testable_L & t.same_sign & (t.p_L < 0.05) & t.beta_L.notna()
    t["primary_rule_mismatch"] = lb & (pr != t.confirmed_final)
    return t


def flags(t):
    out = pd.Series(np.nan, index=t.index, dtype=object)
    for aid, cal in CAL.items():
        ct = pd.read_csv(DIRS / aid / "calibration_D8_tests.tsv.gz", sep="\t")
        ct = ct[ct.calibration == cal]
        mx = ct.assign(a=ct.t.abs()).groupby(["context", "unit"]).a.max().to_dict()
        m = t.entry == aid
        out[m] = [np.nan if (c, u) not in mx else bool(abs(tt) > mx[(c, u)])
                  for c, u, tt in zip(t.context[m], t.unit[m], t.t[m])]
    return out


def rate(df, col):
    lb = df[df.eligible & df.lockbox_testable]
    ne, c = int((lb[col] == "NE").sum()), int((lb[col] == "C").sum())
    return dict(n_lockbox=len(lb), n_ne=ne, n_conf=c, rate_ev=c / (len(lb) - ne) if len(lb) > ne else np.nan,
                rate_all=c / len(lb) if len(lb) else np.nan)


def rates(t):
    zm = pd.read_csv(FZX / "zavg_engine_membership.tsv", sep="\t")
    loc = zm.groupby("unit_id").locus.first().to_dict()
    t = t.assign(locus=[bool(loc.get(str(u).split("|")[0], False)) for u in t.unit])
    groups = {k: t[t.entry == k] for k in H}
    groups["controls (23)"] = t[t.entry.str.startswith("CTRL_")]
    groups[BR + " (GLS hits)"] = t[t.entry == BR]
    rows = []
    for g, df in groups.items():
        subs = {"overall": df}
        if "BRIDGE" not in g:
            subs["non_locus"] = df[~df.locus]
        if g in CAL:
            subs["beats_all_matched_sets"] = df[df.beats_all == True]  # noqa: E712
            subs["not_beats_all"] = df[df.beats_all == False]  # noqa: E712
            subs["no_sets"] = df[df.beats_all.isna()]
        for s, d in subs.items():
            p, c = rate(d, "P"), rate(d, "Cc")
            lb = d[d.eligible & d.lockbox_testable]
            rows.append(dict(direction=g, subset=s, **{f"primary_{k}": v for k, v in p.items()},
                             **{f"corrected_{k}": v for k, v in c.items()},
                             n_defect_sensitive=int((lb.P != lb.Cc).sum())))
    return pd.DataFrame(rows)


SIZE_BINS = ((3, 4), (5, 8), (9, 16), (17, 60))


def sbin(n):
    for lo, hi in SIZE_BINS:
        if lo <= n <= hi:
            return f"{lo}-{hi}"
    return None


def bridge(t, col, B=2000, seed=777, failed_as_not=False):
    sizes = {}
    for m, f in (("zavg", "zavg_engine_membership.tsv"), ("gls", "gls_engine_membership.tsv")):
        mm = pd.read_csv(FZX / f, sep="\t")
        sizes[m] = mm[mm.role == "member"].groupby("unit_id").gene.nunique().to_dict() if "role" in mm else {}
        sizes[m + "_all"] = mm.groupby("unit_id").gene.nunique().to_dict()
    h = t[t.entry.isin(["H3a_lineage_atlas", "H3b_driver_family", BR]) & t.eligible & t.lockbox_testable].copy()
    if not failed_as_not:
        h = h[h[col] != "NE"]
    h["map"] = np.where(h.entry == BR, "gls", "zavg")
    h["y"] = (h[col] == "C").astype(float)
    out = {}
    for skey in ("all", "member"):
        sz = [sizes[m + "_all"].get(u) if skey == "all" else sizes[m].get(u) for m, u in zip(h["map"], h.unit)]
        hh = h.assign(bin=[sbin(int(s)) if s is not None else None for s in sz]).dropna(subset=["bin"])
        w = hh.bin.value_counts(normalize=True)

        def r(d):
            g = d.groupby("bin").y.mean()
            ww = w.reindex(g.index)
            return float((g * ww).sum() / ww.sum())
        obs = r(hh[hh["map"] == "zavg"]) - r(hh[hh["map"] == "gls"])
        rng = np.random.default_rng(seed)
        cells = [(m, g.y.values, b) for (m, b), g in hh.groupby(["map", "bin"])]
        diffs = []
        for _ in range(B):
            acc = {"zavg": {}, "gls": {}}
            for m, y, b in cells:
                acc[m][b] = y[rng.integers(len(y), size=len(y))].mean()
            rz = sum(acc["zavg"][b] * w[b] for b in acc["zavg"]) / sum(w[b] for b in acc["zavg"])
            rg = sum(acc["gls"][b] * w[b] for b in acc["gls"]) / sum(w[b] for b in acc["gls"])
            diffs.append(rz - rg)
        out[skey] = dict(rate_zavg=r(hh[hh["map"] == "zavg"]), rate_gls=r(hh[hh["map"] == "gls"]), diff=obs,
                         ci95=np.quantile(diffs, [0.025, 0.975]).tolist(), n=hh["map"].value_counts().to_dict(),
                         weights=w.round(4).to_dict())
    return out


# ------------------------------------------------------------------------------------------------ power aid
def tweedie(z, grid_n=120, deg=7):
    """Efron's Lindsey-method estimate of log f(z) (Poisson GLM on a degree-7 polynomial of histogram counts); returns
    E[mu|z] = z + (log f)'(z) and Var[mu|z] = 1 + (log f)''(z) at each z."""
    import statsmodels.api as sm
    lo, hi = np.quantile(z, [0.0005, 0.9995])
    edges = np.linspace(lo - 0.5, hi + 0.5, grid_n + 1)
    cnt, _ = np.histogram(z, edges)
    mid = 0.5 * (edges[1:] + edges[:-1])
    sc = (mid - mid.mean()) / mid.std()
    X = np.column_stack([sc ** k for k in range(deg + 1)])
    fit = sm.GLM(cnt, X, family=sm.families.Poisson()).fit()
    b = fit.params
    s = mid.std()
    zz = np.clip(z, mid[0], mid[-1])
    u = (zz - mid.mean()) / s
    d1 = sum(k * b[k] * u ** (k - 1) for k in range(1, deg + 1)) / s
    d2 = sum(k * (k - 1) * b[k] * u ** (k - 2) for k in range(2, deg + 1)) / s ** 2
    return z + d1, np.clip(1 + d2, 0.0, None)


def power_aid(t):
    out, summ = [], {}
    for aid in ("H3a_lineage_atlas", "H3b_driver_family"):
        d = pd.read_parquet(ROOT / "data_cache" / "part2" / "assoc" / aid / "discovery.parquet")
        d = d[d.evaluable.astype(bool)]
        z_all = sst.norm.ppf(sst.t.cdf(d.t.values, d.df.values)).clip(-12, 12)
        disc = d[d.discovery.astype(bool)].copy()
        disc["z"] = sst.norm.ppf(sst.t.cdf(disc.t.values, disc.df.values)).clip(-12, 12)
        Em, Vm = tweedie(z_all)
        # evaluate posterior at the discoveries' z by interpolation over the sorted grid of all z
        o = np.argsort(z_all)
        em = np.interp(disc.z, z_all[o], Em[o])
        vm = np.interp(disc.z, z_all[o], Vm[o])
        disc["shrink_eb"] = np.where(np.abs(disc.z) > 0, em / disc.z, np.nan)
        disc["post_var"] = vm
        h = t[(t.entry == aid) & t.eligible & t.lockbox_testable].merge(
            disc[["context", "unit", "n", "n_context", "se", "z", "shrink_eb", "post_var"]], on=["context", "unit"],
            how="left")
        # use corrected statistics where available (failed primary fits have no L se)
        se_obs = np.where(h.affected & h.se_L_relook.notna(), h.se_L_relook, h.se_L).astype(float)
        nc_L = h.n_context_L.astype(float).values
        nr_L = (h.n_L - h.n_context_L).astype(float).values
        nc_A, nr_A = h.n_context.astype(float).values, (h.n - h.n_context).astype(float).values
        se_pred = h.se.values * np.sqrt((1 / nc_L + 1 / nr_L) / (1 / nc_A + 1 / nr_A))
        df_L = np.where(np.isfinite(h.df_L.astype(float)), h.df_L.astype(float), nc_L - 1).clip(2, None)
        crit = sst.t.ppf(0.95, df_L)
        absb = np.abs(h.beta.values)
        rows = {}
        for se_name, se in (("se_pred_from_A", se_pred), ("se_observed_L", se_obs)):
            k = h.se.values / se                                   # A-to-L noncentrality scale
            for name, eff in (("naive (no shrinkage)", absb), ("50% winner's-curse factor", 0.5 * absb),
                              ("EB posterior mean (Tweedie)", absb * np.clip(h.shrink_eb.values, 0, None))):
                pw = sst.nct.sf(crit, df_L, eff / se)
                rows[(se_name, name)] = pw
            # EB with posterior uncertainty: Z_L ~ N(k m, 1 + k^2 v), m = |E[mu|z]|
            m = np.abs(h.z.values) * np.clip(h.shrink_eb.values, 0, None)
            rows[(se_name, "EB posterior (mean + variance)")] = sst.norm.sf((crit - k * m) / np.sqrt(1 + k ** 2 * h.post_var.values))
        obs_conf = (h.Cc == "C").values if aid == "H3b_driver_family" else (h.P == "C").values
        obs_prim_ev = h[h.P != "NE"]
        summ[aid] = dict(n=len(h), observed_primary_rate_evaluable=float((obs_prim_ev.P == "C").mean()),
                         observed_corrected_rate=float((h.Cc == "C").mean()),
                         median_shrink_eb=float(np.nanmedian(h.shrink_eb)),
                         median_nc_L=float(np.nanmedian(nc_L)), median_se_ratio_obs_over_pred=float(
                             np.nanmedian(se_obs / se_pred)))
        for (sn, nm), pw in rows.items():
            ok = np.isfinite(pw)
            summ[aid][f"{sn} | {nm}"] = dict(expected_rate=float(np.mean(pw[ok])),
                                             sd_rate=float(np.sqrt(np.sum(pw[ok] * (1 - pw[ok]))) / ok.sum()),
                                             n=int(ok.sum()))
        h = h.assign(se_pred=se_pred, se_obs=se_obs, obs_conf=obs_conf,
                     **{f"pw::{sn}::{nm}": pw for (sn, nm), pw in rows.items()})
        out.append(h[["entry", "context", "unit", "beta", "z", "shrink_eb", "post_var", "n_context", "n_context_L",
                      "se", "se_pred", "se_obs", "obs_conf"] + [c for c in h.columns if c.startswith("pw::")]])
    return pd.concat(out), summ


# ------------------------------------------------------------------------------------------------ BENCH
def bench():
    b = pd.read_csv(OUT / "BENCH_complex_vs_gene" / "bench_pairs.tsv", sep="\t")
    js = json.loads((OUT / "BENCH_complex_vs_gene" / "bench.json").read_text())
    r = b.breadth.rank(method="first")
    n = len(b)
    # terciles fixed on the primary pairs' A breadth (equal thirds, ties by pair order)
    terc = np.where(r <= n / 3, "low", np.where(r <= 2 * n / 3, "mid", "high"))
    out = dict(n_pairs=n, n_delta_nan=int(b.delta.isna().sum()), mean=float(b.delta.mean()),
               tercile_labels_match=bool((terc == b.tercile.values).all()),
               n_tercile_mismatch=int((terc != b.tercile.values).sum()),
               **{f"mean_{k}": float(b.delta[b.tercile == k].mean()) for k in ("low", "mid", "high")},
               **{f"n_{k}": int((b.tercile == k).sum()) for k in ("low", "mid", "high")},
               spearman=float(sst.spearmanr(b.delta, b.breadth).statistic),
               reported=js["primary"], reported_spearman=js["trend_spearman_delta_breadth"])
    # how the per-pair spread compares with the bootstrap CI (descriptive)
    out["naive_se_mean_iid_pairs"] = float(b.delta.std() / np.sqrt(n))
    hi = b[b.tercile == "high"]
    out["naive_se_high_iid_pairs"] = float(hi.delta.std() / np.sqrt(len(hi)))
    out["n_distinct_units"] = int(b.unit.nunique())
    out["mean_delta_by_complex_then_mean"] = float(b.groupby("unit").delta.mean().mean())
    out["high_by_complex_then_mean"] = float(hi.groupby("unit").delta.mean().mean())
    out["frac_pairs_delta_pos"] = float((b.delta > 0).mean())
    out["frac_high_delta_pos"] = float((hi.delta > 0).mean())
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    t = load()
    t["beats_all"] = flags(t)
    res = {}
    res["primary_rule_mismatch"] = int(t.primary_rule_mismatch.sum())
    res["relook_rule_mismatch"] = int(t.relook_rule_mismatch.sum())
    rt = rates(t)
    rt.to_csv(a.out / "rates_independent.tsv", sep="\t", index=False, float_format="%.4f")
    prod = pd.read_csv(OUT / "comparison_primary_vs_corrected_R7.tsv", sep="\t")
    mm = rt.merge(prod, on=["direction", "subset"], how="left")
    cmp = []
    for my, pc in (("primary_n_lockbox", "primary_n_lockbox"), ("primary_n_ne", "primary_n_not_evaluable"),
                   ("primary_n_conf", "primary_n_confirmed"), ("primary_rate_ev", "primary_rate_evaluable"),
                   ("primary_rate_all", "primary_rate_all_in"), ("corrected_n_ne", "corrected_n_not_evaluable"),
                   ("corrected_n_conf", "corrected_n_confirmed"), ("corrected_rate_ev", "corrected_rate_evaluable"),
                   ("corrected_rate_all", "corrected_rate_all_in"), ("n_defect_sensitive_x", "n_defect_sensitive_y")):
        if my in mm and pc in mm:
            x, y = mm[my].astype(float), mm[pc].astype(float)
            bad = mm[(~np.isclose(x, y, atol=5e-5)) & y.notna()]
            cmp.append(dict(field=pc, n_rows=int(y.notna().sum()), n_mismatch=int(len(bad)),
                            rows=bad[["direction", "subset"]].values.tolist()))
    res["rate_table_vs_producer"] = cmp
    lb = t[t.eligible & t.lockbox_testable]
    ds = lb[lb.P != lb.Cc]
    ds[["entry", "context", "unit", "beta", "beta_L", "p_L", "P", "Cc", "affected_cause", "beta_L_relook",
        "p_L_relook"]].to_csv(a.out / "defect_sensitive_independent.tsv", sep="\t", index=False, float_format="%.6g")
    pds = pd.read_csv(OUT / "defect_sensitive_calls_R7.tsv", sep="\t")
    mine = set(zip(ds.entry, ds.context, ds.unit))
    theirs = set(zip(pds.entry, pds.context, pds.unit))
    res["defect_sensitive"] = dict(n_mine=len(mine), n_producer=len(theirs), only_mine=sorted(mine - theirs),
                                   only_producer=sorted(theirs - mine),
                                   by_transition=ds.groupby(["entry", "P", "Cc"]).size().reset_index()
                                   .values.tolist(),
                                   flips_among_fitted=ds[ds.P != "NE"][["entry", "context", "unit", "p_L",
                                                                       "p_L_relook"]].values.tolist())
    # E2 / promotion
    e2 = []
    for aid in ("H1_within_complex_GI", "H1S_silencing_GI"):
        c = t[t.entry == aid]
        lab_diff = c.label_final.astype(str) != c.class_engine.astype(str)
        not_he = c.headline_eligible.astype(str) == "False"
        bad = c[(lab_diff | not_he) & tb(c.promotable)]
        e2.append(dict(entry=aid, n_e2_or_not_headline=int((lab_diff | not_he).sum()),
                       n_of_those_confirmed=int(((lab_diff | not_he) & c.confirmed_final).sum()),
                       n_promoted=int(len(bad)), n_promotable_total=int(tb(c.promotable).sum()),
                       promotable_units=c[tb(c.promotable)].unit.tolist()))
    ev = pd.read_csv(OUT / "h1_events_primary_vs_corrected_R7.tsv", sep="\t")
    e2.append(dict(file="h1_events_primary_vs_corrected_R7.tsv",
                   promotable_primary_true=ev[tb(ev.promotable_primary)][["event", "axis"]].values.tolist(),
                   corrected_confirmed=ev[ev.corrected_call == "confirmed"][["event", "axis", "label_final"]]
                   .values.tolist()))
    res["e2"] = e2
    res["bridge"] = {nm: bridge(t, col, failed_as_not=fa) for nm, col, fa in
                     (("as_executed", "P", True), ("primary_evaluable", "P", False),
                      ("corrected", "Cc", False))}
    pw, ps = power_aid(t)
    pw.to_csv(a.out / "power_aid.tsv", sep="\t", index=False, float_format="%.5g")
    res["power"] = ps
    # expected vs observed within the 'beats all matched sets' split (is the split explained by power alone?)
    fl = t[["entry", "context", "unit", "beats_all", "P", "Cc"]]
    pm = pw.merge(fl, on=["entry", "context", "unit"], how="left")
    key = "pw::se_observed_L::EB posterior (mean + variance)"
    key2 = "pw::se_pred_from_A::EB posterior (mean + variance)"
    sp = {}
    for (aid, ba), g in pm.groupby(["entry", pm.beats_all.astype(str)]):
        sp[f"{aid}|beats_all={ba}"] = dict(n=len(g), obs_primary=float((g.P == "C").sum() / max(1, (g.P != "NE").sum())),
                                           obs_corrected=float((g.Cc == "C").mean()),
                                           exp_eb_seL=float(g[key].mean()), exp_eb_seA=float(g[key2].mean()),
                                           sd_exp=float(np.sqrt((g[key] * (1 - g[key])).sum()) / len(g)))
    res["power_by_beats_all"] = sp
    res["bench"] = bench()
    (a.out / "tables_review.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
    print(json.dumps(res, indent=1, default=str)[:12000])


if __name__ == "__main__":
    main()
