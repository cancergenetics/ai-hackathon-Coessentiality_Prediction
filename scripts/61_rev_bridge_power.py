#!/usr/bin/env python
"""R-4 (referee M3.3): power-expected BRIDGE-L confirmation rates for the zavg-map and GLS-map discoveries.
POST HOC, DESCRIPTIVE (results/part2/revision/INTENT.md).

The lockbox results review's power aid (scripts/verify/lockbox_results_review_tables.py: `tweedie`, empirical-Bayes
posterior from all evaluable A tests of the same map and family; one-sided alpha = 0.05 L test) with A estimates and
committed L SAMPLE SIZES only: the L SE is predicted from the A SE and the L / A context and reference counts, the L
df from the A df scaled by the context counts. No L estimate, SE, df or p is read. The observed rates use the
committed R7 primary call labels (primary-evaluable BRIDGE-L, -0.065).
Outputs: results/part2/revision/r4_*.  Usage: .venv/bin/python scripts/61_rev_bridge_power.py
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import bench, revision as rv  # noqa: E402
from part2.cmap import FrozenMap  # noqa: E402

OUT = ROOT / "results" / "part2" / "revision"
LC = ROOT / "results" / "part2" / "lockbox_confirm"
BR = "BRIDGE-L_confirmation_rate"
ENTRIES = ("H3a_lineage_atlas", "H3b_driver_family", BR)
A_TABLES = {("zavg", "H3a"): ("data_cache/part2/assoc/H3a_lineage_atlas/discovery.parquet",
                              "results/part2/assoc/H3a_lineage_atlas/discovery_summary.json"),
            ("zavg", "H3b"): ("data_cache/part2/assoc/H3b_driver_family/discovery.parquet",
                              "results/part2/assoc/H3b_driver_family/discovery_summary.json"),
            ("gls", "H3a"): (f"data_cache/part2/assoc/{BR}/discovery_H3a.parquet",
                             f"results/part2/assoc/{BR}/discovery_H3a_summary.json"),
            ("gls", "H3b"): (f"data_cache/part2/assoc/{BR}/discovery_H3b.parquet",
                             f"results/part2/assoc/{BR}/discovery_H3b_summary.json")}
MAPS = {"zavg": "results/part2/fzx/zavg_engine_membership.tsv", "gls": "results/part2/fzx/gls_engine_membership.tsv"}


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def a_table(key) -> pd.DataFrame:
    p, s = (ROOT / x for x in A_TABLES[key])
    if sha(p) != json.loads(s.read_text())["full_table_sha256"]:
        raise SystemExit(f"{p} differs from its committed summary")
    d = pd.read_parquet(p)
    return d[d.evaluable.astype(bool)].copy()


def hits_table() -> pd.DataFrame:
    """Eligible, lockbox-testable, evaluable BRIDGE-L hits with committed L sample sizes and R7 primary calls."""
    lab = pd.read_csv(LC / "primary_relabelled_R7.tsv.gz", sep="\t",
                      usecols=["entry", "hit_index", "context", "unit", "eligible", "primary_call"])
    rows = []
    for e in ENTRIES:
        c = pd.read_csv(LC / e / "confirmation.tsv", sep="\t",
                        usecols=["context", "unit", "eligible", "lockbox_testable", "n_L", "n_context_L"])
        c["entry"], c["hit_index"] = e, range(len(c))
        rows.append(c)
    c = pd.concat(rows, ignore_index=True)
    h = lab[lab.entry.isin(ENTRIES)].merge(c, on=["entry", "hit_index"], suffixes=("", "_c"), validate="one_to_one")
    assert (h.context == h.context_c).all() and (h.unit == h.unit_c).all()
    tb = lambda s: s.astype(str).str.lower().eq("true")  # noqa: E731
    h = h[tb(h.eligible) & tb(h.lockbox_testable)]
    n_all = len(h)
    h = h[~h.primary_call.astype(str).str.startswith("NOT EVALUABLE")].copy()
    h["map"] = np.where(h.entry == BR, "gls", "zavg")
    h["confirmed"] = (h.primary_call == "confirmed").astype(float)
    h.attrs["n_lockbox_testable"] = n_all
    return h.drop(columns=["context_c", "unit_c", "eligible_c"])


def main():
    t0 = time.time()
    spec = importlib.util.spec_from_file_location("lrr", ROOT / "scripts" / "verify" / "lockbox_results_review_tables.py")
    lrr = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(lrr)
    h = hits_table()
    parts = []
    for (m, fam), _ in A_TABLES.items():
        d = a_table((m, fam))
        d["z"] = rv.t_to_z(d.t, d.df)
        sub = h[(h["map"] == m) & h.context.isin(set(d.context))].copy()
        sub["family"] = fam
        em, vm = rv.tweedie_posterior(d.z.values, d.z.values, lrr.tweedie)
        d["shrink_eb"] = np.where(np.abs(d.z) > 0, em / d.z, np.nan)
        d["post_var"] = vm
        sub = sub.merge(d[["context", "unit", "beta", "se", "t", "df", "n", "n_context", "z", "shrink_eb",
                           "post_var"]].rename(columns={"beta": "beta_A"}), on=["context", "unit"], how="left",
                        validate="one_to_one")
        parts.append(sub)
    x = pd.concat(parts, ignore_index=True)
    assert len(x) == len(h) and x.z.notna().all(), "every hit must map to one A test"
    x["beta"] = x.beta_A
    ncA, nrA = x.n_context.astype(float), (x.n - x.n_context).astype(float)
    ncL, nrL = x.n_context_L.astype(float), (x.n_L - x.n_context_L).astype(float)
    x["se_pred_L"] = x.se * np.sqrt((1 / ncL + 1 / nrL) / (1 / ncA + 1 / nrA))
    x["df_pred_L"] = (x.df * (ncL - 1) / (ncA - 1)).clip(lower=2)
    x["df_ncL_minus_1"] = (ncL - 1).clip(lower=2)
    for dfn in ("df_pred_L", "df_ncL_minus_1"):
        pw = rv.power_eb(x.z.values, x.shrink_eb.values, x.post_var.values, x.se.values, x.se_pred_L.values,
                         x[dfn].values, np.abs(x.beta.values))
        for k, v in pw.items():
            x[f"pw_{k}|{dfn}"] = v
    fm = {m: FrozenMap.from_membership(pd.read_csv(ROOT / p, sep="\t")) for m, p in MAPS.items()}
    x["size"] = [len(fm[m].members[u]) for m, u in zip(x["map"], x.unit)]
    x["bin"] = [bench.size_bin(int(s)) for s in x["size"]]
    x = x.dropna(subset=["bin"])
    w = x.bin.value_counts(normalize=True)
    res = dict(label="POST HOC, DESCRIPTIVE (revision R-4); A estimates and committed L sample sizes only",
               n_hits={m: int((x["map"] == m).sum()) for m in ("zavg", "gls")},
               n_lockbox_testable_incl_not_evaluable=int(h.attrs["n_lockbox_testable"]),
               bin_weights=w.round(4).to_dict())
    obs = {m: rv.standardised_rate(x[x["map"] == m].confirmed.values, x[x["map"] == m].bin.values, w)[0]
           for m in ("zavg", "gls")}
    res["observed"] = dict(rate_zavg=obs["zavg"], rate_gls=obs["gls"], diff=obs["zavg"] - obs["gls"],
                           committed_R7_primary_evaluable=dict(diff=-0.0648, ci95=[-0.1069, -0.0246]))
    rows = []
    for col in [c for c in x.columns if c.startswith("pw_")]:
        ver, dfn = col[3:].split("|")
        e = {m: rv.standardised_rate(x[x["map"] == m][col].values, x[x["map"] == m].bin.values, w) for m in ("zavg", "gls")}
        diff = e["zavg"][0] - e["gls"][0]
        sd = float(np.hypot(e["zavg"][1], e["gls"][1]))
        od = obs["zavg"] - obs["gls"]
        rows.append(dict(power=ver, df=dfn, exp_zavg=e["zavg"][0], sd_zavg=e["zavg"][1], exp_gls=e["gls"][0],
                         sd_gls=e["gls"][1], exp_diff=diff, sd_diff=sd, obs_zavg=obs["zavg"], obs_gls=obs["gls"],
                         obs_diff=od, obs_minus_exp=od - diff, z_obs_minus_exp=(od - diff) / sd,
                         exp_zavg_unweighted=float(x[x["map"] == "zavg"][col].mean()),
                         exp_gls_unweighted=float(x[x["map"] == "gls"][col].mean()),
                         obs_zavg_unweighted=float(x[x["map"] == "zavg"].confirmed.mean()),
                         obs_gls_unweighted=float(x[x["map"] == "gls"].confirmed.mean())))
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "r4_bridge_power.tsv", sep="\t", index=False, float_format="%.4f")
    head = "pw_eb_mean_var|df_pred_L"
    desc = []
    for (m, fam), g in x.groupby(["map", "family"]):
        desc.append(dict(map=m, family=fam, n=len(g), median_abs_z_A=float(g.z.abs().median()),
                         median_shrink_eb=float(g.shrink_eb.median()), median_n_context_A=float(g.n_context.median()),
                         median_n_context_L=float(g.n_context_L.median()), median_size=float(g["size"].median()),
                         median_abs_beta_A=float(g.beta.abs().median()), expected_rate=float(g[head].mean()),
                         observed_rate=float(g.confirmed.mean())))
    for m, g in x.groupby("map"):
        desc.append(dict(map=m, family="all", n=len(g), median_abs_z_A=float(g.z.abs().median()),
                         median_shrink_eb=float(g.shrink_eb.median()), median_n_context_A=float(g.n_context.median()),
                         median_n_context_L=float(g.n_context_L.median()), median_size=float(g["size"].median()),
                         median_abs_beta_A=float(g.beta.abs().median()), expected_rate=float(g[head].mean()),
                         observed_rate=float(g.confirmed.mean())))
    pd.DataFrame(desc).to_csv(OUT / "r4_by_map_family.tsv", sep="\t", index=False, float_format="%.4f")
    # expected vs observed per size bin (headline power)
    pb = x.groupby(["map", "bin"]).agg(n=("confirmed", "size"), observed=("confirmed", "mean"),
                                        expected=(head, "mean")).reset_index()
    pb.to_csv(OUT / "r4_by_size_bin.tsv", sep="\t", index=False, float_format="%.4f")
    x[["map", "family", "entry", "context", "unit", "size", "bin", "beta", "se", "df", "n_context", "n",
       "n_context_L", "n_L", "z", "shrink_eb", "post_var", "se_pred_L", "df_pred_L", "confirmed"]
      + [c for c in x.columns if c.startswith("pw_")]].to_csv(OUT / "r4_per_hit.tsv.gz", sep="\t", index=False,
                                                               float_format="%.5g",
                                                               compression=dict(method="gzip", mtime=0))
    res["headline"] = tab[(tab.power == "eb_mean_var") & (tab.df == "df_pred_L")].iloc[0].to_dict()
    res["runtime_s"] = round(time.time() - t0, 1)
    res["peak_rss_bytes"] = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (OUT / "r4_summary.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    print(json.dumps(res, indent=1, default=float))
    print(tab.to_string(index=False))
    print(pd.DataFrame(desc).to_string(index=False))
    print(pb.to_string(index=False))


if __name__ == "__main__":
    main()
