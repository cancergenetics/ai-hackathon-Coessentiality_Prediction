#!/usr/bin/env python
"""PI ruling R7 reporting layer (results/part2/lockbox_confirm/README.md). It reads committed tables only: no lockbox
outcome and no persisted read.

PRIMARY: the single read as executed (5bc9444). Eligible lockbox-route hits with no L estimate are relabelled NOT
  EVALUABLE, as either "driver defect" or "genuine": the cause comes from the corrected relook, where "genuine" means
  the fixed design still cannot be fitted. Rates are shown on the evaluable denominator and all-in, where NOT EVALUABLE
  counts as not confirmed.
CORRECTED (disclosed, NON-CONFIRMATORY): the 560 affected fits from corrected_refit_nonconfirmatory.tsv replace the
  primary fits; every other fit is identical (reproduction check 68f0d3f).
DEFECT-SENSITIVE: an eligible call whose status (confirmed / not confirmed / NOT EVALUABLE) differs between primary and
  corrected. It cannot support a claim.
Outputs (results/part2/lockbox_confirm/): primary_relabelled_R7.tsv.gz, rates_primary_R7.tsv,
  comparison_primary_vs_corrected_R7.tsv, defect_sensitive_calls_R7.tsv, h1_events_primary_vs_corrected_R7.tsv,
  bridge_primary_vs_corrected_R7.json.
"""
from __future__ import annotations

import gzip
import io
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import bench, registry, validation  # noqa: E402
from part2.cmap import FrozenMap  # noqa: E402

OUT = ROOT / "results" / "part2" / "lockbox_confirm"
DIRS = ROOT / "results" / "part2" / "directions"
CORR = OUT / "corrected_refit_nonconfirmatory.tsv"
ORDER = ("H3a_lineage_atlas", "H3b_driver_family", "H1_within_complex_GI", "H1S_silencing_GI",
         "H2_submodule_divergence")
BRIDGE = "BRIDGE-L_confirmation_rate"
CAL = {"H3a_lineage_atlas": "matched_random_set", "H3b_driver_family": "matched_random_set",
       "H2_submodule_divergence": "random_split"}


def committed(p: Path) -> Path:
    registry.committed_blob(ROOT, p)
    return p


def _bool(s: pd.Series) -> pd.Series:
    return s.astype(str).str.lower() == "true"


def load_primary() -> pd.DataFrame:
    rows = []
    for p in sorted(OUT.glob("*/confirmation.tsv")):
        t = pd.read_csv(committed(p), sep="\t")
        t.insert(0, "entry", p.parent.name)
        t.insert(1, "hit_index", range(len(t)))
        rows.append(t)
    t = pd.concat(rows, ignore_index=True)
    for c in ("eligible", "lockbox_testable", "confirmed_final"):
        t[c] = _bool(t[c])
    return t


def flags(t: pd.DataFrame) -> pd.Series:
    out = pd.Series(np.nan, index=t.index, dtype=object)
    for aid, cal in CAL.items():
        m = (t.entry == aid).values
        ct = pd.read_csv(committed(DIRS / aid / "calibration_D8_tests.tsv.gz"), sep="\t")
        ct = ct[ct.calibration == cal]
        mx = ct.assign(at=ct.t.abs()).groupby(["context", "unit"]).at.max()
        sub = t[m]
        v = [mx.get((c, u), np.nan) for c, u in zip(sub.context, sub.unit)]
        out[m] = [np.nan if not np.isfinite(x) else bool(abs(tt) > x) for tt, x in zip(sub.t, v)]
    return out


def statuses(t: pd.DataFrame, corr: pd.DataFrame) -> pd.DataFrame:
    key = corr.set_index(["entry", "hit_index"])
    t = t.copy()
    lb = t.eligible & t.lockbox_testable
    fail = lb & t.beta_L.isna()
    cause = [key.affected_cause.get((e, i), "") for e, i in zip(t.entry, t.hit_index)]
    t["affected_cause"] = cause
    t["primary_call"] = np.where(~lb, np.where(t.eligible, t.status.astype(str), ""),
                                 np.where(t.confirmed_final, "confirmed", "not confirmed"))
    t.loc[fail, "primary_call"] = ["NOT EVALUABLE (genuine)" if c == "failed fit (genuine)" else
                                   "NOT EVALUABLE (driver defect)" for c in t.affected_cause[fail]]
    if (fail & ~t.affected_cause.str.startswith("failed")).any():
        raise RuntimeError("a failed primary fit is missing from the corrected relook")
    cc = t.primary_call.copy()
    for i in t.index[lb & (t.affected_cause != "")]:
        r = key.loc[(t.entry[i], t.hit_index[i])]
        if not np.isfinite(float(r.beta_L_relook)):
            cc[i] = "NOT EVALUABLE (genuine)"
        else:
            cc[i] = "confirmed" if str(r.confirmed_relook).lower() == "true" else "not confirmed"
    t["corrected_call"] = cc
    t["defect_sensitive"] = lb & (t.primary_call != t.corrected_call)
    return t


def _rate(df: pd.DataFrame, col: str) -> dict:
    lb = df[df.eligible & df.lockbox_testable]
    ne = lb[col].str.startswith("NOT EVALUABLE")
    conf = int((lb[col] == "confirmed").sum())
    ev = int((~ne).sum())
    return dict(n_lockbox=int(len(lb)), n_not_evaluable=int(ne.sum()), n_confirmed=conf,
                rate_evaluable=conf / ev if ev else None, rate_all_in=conf / len(lb) if len(lb) else None)


def rates(t: pd.DataFrame, zlocus: dict) -> pd.DataFrame:
    t = t.assign(locus=[bool(zlocus.get(str(u).split("|")[0], False)) for u in t.unit])
    groups = {k: t[t.entry == k] for k in ORDER}
    groups["controls (23)"] = t[t.entry.str.startswith("CTRL_")]
    groups[BRIDGE + " (GLS hits)"] = t[t.entry == BRIDGE]
    rows = []
    for g, df in groups.items():
        subsets = {"overall": df, "non_locus": df[~df.locus] if "BRIDGE" not in g else None}
        if g in CAL:
            subsets["beats_all_matched_sets"] = df[df.beats_all == True]  # noqa: E712
            subsets["not_beats_all"] = df[df.beats_all == False]  # noqa: E712
        for s, d in subsets.items():
            if d is None:
                continue
            p, c = _rate(d, "primary_call"), _rate(d, "corrected_call")
            rows.append(dict(direction=g, subset=s, **{f"primary_{k}": v for k, v in p.items()},
                             **{f"corrected_{k}": v for k, v in c.items()},
                             n_defect_sensitive=int(d.defect_sensitive.sum())))
    return pd.DataFrame(rows)


def bridge(t: pd.DataFrame) -> dict:
    zm = FrozenMap.from_membership(pd.read_csv(committed(ROOT / "results/part2/fzx/zavg_engine_membership.tsv"),
                                               sep="\t"))
    gm = FrozenMap.from_membership(pd.read_csv(committed(ROOT / "results/part2/fzx/gls_engine_membership.tsv"),
                                               sep="\t"))
    cl = {k: pd.read_csv(committed(ROOT / f"results/part2/construction/{k}_map/clusters.tsv"), sep="\t")
          for k in ("zavg", "gls")}
    mz, mg, _ = validation.match_coverage(cl["zavg"][cl["zavg"].level == 1], cl["gls"][cl["gls"].level == 1])
    keep = {("zavg", u) for u in mz.cluster_id} | {("gls", u) for u in mg.cluster_id}
    fm = {"zavg": zm, "gls": gm}
    executed = json.loads(committed(OUT / BRIDGE / "bridge.json").read_text())
    out = dict(as_executed=dict(primary=executed["primary"], multiverse=executed["multiverse"]))
    for name, col in (("primary_evaluable", "primary_call"), ("corrected_nonconfirmatory", "corrected_call")):
        h = t[t.entry.isin(["H3a_lineage_atlas", "H3b_driver_family", BRIDGE]) & t.lockbox_testable & t.eligible]
        h = h[~h[col].str.startswith("NOT EVALUABLE")]
        hits = pd.DataFrame({"map": np.where(h.entry == BRIDGE, "gls", "zavg"), "unit": h.unit.values,
                             "confirmed": (h[col] == "confirmed").values})
        hits["size"] = [len(fm[m].members[u]) for m, u in zip(hits["map"], hits.unit)]
        hits["tier"] = [fm[m].tier.get(u, "") for m, u in zip(hits["map"], hits.unit)]
        cov = hits[[(m, u) in keep for m, u in zip(hits["map"], hits.unit)]]
        out[name] = dict(n_hits={m: int((hits["map"] == m).sum()) for m in ("zavg", "gls")},
                         primary=bench.bridge_rate(hits, B=2000, seed=5200),
                         multiverse=dict(unweighted=bench.bridge_unweighted(hits, B=2000, seed=5200),
                                         matched_gene_coverage=bench.bridge_rate(cov, B=2000, seed=5200),
                                         T1_T2_clusters=bench.bridge_rate(hits[hits.tier.isin(["T1", "T2"])],
                                                                          B=2000, seed=5200)))
    return out


def h1_table(t: pd.DataFrame, corr: pd.DataFrame) -> pd.DataFrame:
    key = corr.set_index(["entry", "hit_index"])
    h = t[t.entry == "H1_within_complex_GI"].copy()
    rows = []
    for r in h.itertuples(index=False):
        c = key.loc[("H1_within_complex_GI", r.hit_index)] if ("H1_within_complex_GI", r.hit_index) in key.index \
            else None
        rows.append(dict(event=r.event, axis=r.axis, class_engine=r.class_engine, label_final=r.label_final,
                         eligible=r.eligible, route=r.route, beta_A=r.beta,
                         primary_beta_L=r.beta_L, primary_ci95=[r.ci95_L_lo, r.ci95_L_hi], primary_p_L=r.p_L,
                         primary_call=r.primary_call or "(not eligible)",
                         corrected_beta_L=r.beta_L if c is None else c.beta_L_relook,
                         corrected_ci95=[r.ci95_L_lo, r.ci95_L_hi] if c is None else
                         [float(c.ci95_L_lo_relook), float(c.ci95_L_hi_relook)],
                         corrected_p_L=r.p_L if c is None else c.p_L_relook,
                         corrected_call=r.corrected_call or "(not eligible)", affected=c is not None,
                         defect_sensitive=r.defect_sensitive, promotable_primary=r.promotable,
                         rare_route_pass=r.rare_route_pass))
    return pd.DataFrame(rows)


def main():
    corr = pd.read_csv(committed(CORR), sep="\t")
    t = load_primary()
    t["beats_all"] = flags(t)
    t = statuses(t, corr)
    zm = pd.read_csv(committed(ROOT / "results/part2/fzx/zavg_engine_membership.tsv"), sep="\t")
    zlocus = zm.groupby("unit_id").locus.first().to_dict()
    rt = rates(t, zlocus)
    cols = ["entry", "hit_index", "context", "unit", "eligible", "route", "beta", "beta_L", "p_L", "confirmed_final",
            "status", "primary_call", "corrected_call", "affected_cause", "defect_sensitive", "beats_all"]
    buf = io.StringIO()
    t[[c for c in cols if c in t]].to_csv(buf, sep="\t", index=False, float_format="%.6g")
    with open(OUT / "primary_relabelled_R7.tsv.gz", "wb") as fh, gzip.GzipFile(fileobj=fh, mode="wb", mtime=0) as gz:
        gz.write(buf.getvalue().encode())
    rt[[c for c in rt.columns if not c.startswith("corrected_")]].to_csv(OUT / "rates_primary_R7.tsv", sep="\t",
                                                                          index=False, float_format="%.4f")
    rt.to_csv(OUT / "comparison_primary_vs_corrected_R7.tsv", sep="\t", index=False, float_format="%.4f")
    ds = t[t.defect_sensitive][["entry", "context", "unit", "beta", "beta_L", "p_L", "primary_call",
                                "corrected_call", "affected_cause"]]
    ds = ds.merge(corr[["entry", "unit", "context", "beta_L_relook", "p_L_relook"]], on=["entry", "unit", "context"],
                  how="left")
    ds.to_csv(OUT / "defect_sensitive_calls_R7.tsv", sep="\t", index=False, float_format="%.6g")
    h1_table(t, corr).to_csv(OUT / "h1_events_primary_vs_corrected_R7.tsv", sep="\t", index=False,
                             float_format="%.6g")
    b = bridge(t)
    (OUT / "bridge_primary_vs_corrected_R7.json").write_text(json.dumps(b, indent=1, default=float) + "\n")
    print(rt.to_string())
    print(ds.groupby(["entry", "primary_call", "corrected_call"]).size().to_string())
    print(json.dumps({k: v["primary"] for k, v in b.items()}, indent=1, default=float))


if __name__ == "__main__":
    main()
