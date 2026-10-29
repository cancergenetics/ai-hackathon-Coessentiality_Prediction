#!/usr/bin/env python
"""Referee round 2: post hoc computations A1-A3 (intent: results/part2/revision/INTENT_round2.md, committed first).

POST HOC and DESCRIPTIVE. Reads committed result tables only: no dataset, no A/L dependency matrix, no lockbox value,
no persisted read, no test label. Nothing is refitted; every value is a filter, count or re-sort of committed columns.
No registered result, label, call, rank or ledger is changed.

A1  H1 neutral events (label_final NEUTRAL, registered margin +/-0.1) whose committed 90% CI lies inside +/-0.05.
A2  Strongest-findings counts and order under the alternative readings of the two ambiguous registry items
    (H1 LOSO unit; "A statistic x tier"), using src/part2/strongest.py (apply_rule logic, TIER_CLASS, LOSO_UNIT).
A3  Per-atlas range of the R-1(b) "beats all others" rate among random-atlas discoveries.

Outputs (results/part2/revision/): r5_h1_margin_005.tsv, r5_ranking_alternatives.tsv, r5_atlas_flag_range.tsv,
r5_summary.json.

Usage: PYTHONHASHSEED=0 .venv/bin/python scripts/62_rev_round2.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import strongest  # noqa: E402

REV = ROOT / "results" / "part2" / "revision"
ROB = ROOT / "results" / "part2" / "robustness"
H1C = ROOT / "results" / "part2" / "directions" / "H1_within_complex_GI" / "classes.tsv"
MARGIN, HALF = 0.1, 0.05
HIGHLIGHT = [("lineage:Ovary/Fallopian Tube", "Z0035"), ("lineage:Bone", "Z0303"), ("lineage:Head and Neck", "Z0051"),
             ("lineage:Myeloid", "Z0254"), ("LoF:RB1", "Z0069")]
ALT_WEIGHT = {1: 3.0, 2: 2.0, 3: 1.0}       # R-alt2 tier weights (INTENT_round2.md; not registered)


def _tsv(p: Path) -> pd.DataFrame:
    return pd.read_csv(p, sep="\t")


# ------------------------------------------------------------------------------------------------------------- A1
def a1() -> dict:
    c = _tsv(H1C)
    n = c[c.label_final.eq("NEUTRAL")].copy()
    assert n.N_tost.astype(str).eq("True").all(), "a NEUTRAL event without TOST at the registered margin"
    assert ((n.kappa_N_ci90_lo > -MARGIN) & (n.kappa_N_ci90_hi < MARGIN)).all()
    n["inside_005_strict"] = (n.kappa_N_ci90_lo > -HALF) & (n.kappa_N_ci90_hi < HALF)
    n["inside_005_closed"] = (n.kappa_N_ci90_lo >= -HALF) & (n.kappa_N_ci90_hi <= HALF)
    n["min_margin"] = np.maximum(n.kappa_N_ci90_lo.abs(), n.kappa_N_ci90_hi.abs())
    out = n[["event", "gene", "n_A", "kappa_N", "kappa_N_ci90_lo", "kappa_N_ci90_hi", "min_margin",
             "inside_005_strict", "inside_005_closed", "label_final"]].sort_values("min_margin")
    out.to_csv(REV / "r5_h1_margin_005.tsv", sep="\t", index=False, float_format="%.6g")
    return dict(n_neutral=int(len(n)), n_inside_005=int(n.inside_005_strict.sum()),
                n_inside_005_closed=int(n.inside_005_closed.sum()),
                min_margin_range=[round(float(n.min_margin.min()), 4), round(float(n.min_margin.max()), 4)],
                min_margin_median=round(float(n.min_margin.median()), 4),
                n_events_total=int(len(c)))


# ------------------------------------------------------------------------------------------------------------- A2
def _rank(s: pd.DataFrame, key: pd.Series) -> pd.Series:
    """Rank 1..n by key descending; ties broken by context, unit (as in strongest.apply_rule)."""
    x = s.assign(_k=key.values)
    x = x.sort_values(["_k", "context", "unit"], ascending=[False, True, True], kind="mergesort")
    return pd.Series(np.arange(1, len(x) + 1), index=x.index).reindex(s.index)


def a2() -> dict:
    s = _tsv(ROB / "strongest_findings.tsv")
    cand = _tsv(ROB / "candidates.tsv")
    per = _tsv(ROB / "robustness_per_hit.tsv")
    # registered order reproduces from the committed columns with the apply_rule logic
    tc = s.tier.map(strongest.TIER_CLASS).fillna(9).astype(int)
    assert tc.equals(s.tier_class.astype(int))
    reg = (s.assign(_tc=tc, _abs=s.t_A.abs())
           .sort_values(["_tc", "_abs", "context", "unit"], ascending=[True, False, True, True], kind="mergesort"))
    assert list(reg["rank"]) == list(range(1, len(s) + 1)), "registered rank does not reproduce"

    # (i) H1 LOSO unit
    assert set(per[per.direction == "H3a"].loso_unit) == {strongest.LOSO_UNIT["H3a"]}
    assert set(per[per.direction == "H3b"].loso_unit) == {strongest.LOSO_UNIT["H3b"]}
    h1 = cand[cand.direction == "H1"]
    h1c = h1[h1.candidate.astype(str).eq("True")]
    why = {}
    for r in h1[~h1.candidate.astype(str).eq("True")].itertuples(index=False):
        reasons = []
        if str(r.defect_sensitive) == "True":
            reasons.append("defect-sensitive")
        if str(r.lockbox_pass) != "True":
            reasons.append(f"lockbox call: {r.primary_call}")
        if str(r.headline_eligible) != "True":
            reasons.append("not headline-eligible (E2)")
        why[r.unit] = "; ".join(reasons)
    h1s = s[s.direction == "H1"][["unit", "rank", "route", "loso_unit", "loso_folds", "loso_pass", "loso_skipped"]]
    counts = {d: int((s.direction == d).sum()) for d in ("H3a", "H3b", "H1")}
    k = len(h1c)
    loso = dict(registered_unit=strongest.LOSO_UNIT["H1"], alternative_unit="OncotreePrimaryDisease",
                counts_registered=counts, total_registered=int(len(s)),
                H3a_H3b_unchanged=True, n_h1_candidates=k, h1_candidates=sorted(h1c.unit),
                h1_not_candidates=why,
                total_alternative_bound=[int(len(s)) - k, int(len(s))],
                h1_alternative_bound=[counts["H1"] - k, counts["H1"]],
                h1_status_under_alternative="not determined: per-fold H1 fits by primary disease are not committed "
                                            "and no refit is made in this round",
                h1_lineage_folds=h1s.to_dict(orient="records"),
                rank_effect="if an H1 call dropped, every finding ranked below it moves up by one; order otherwise "
                            "unchanged")

    # (ii) ranking alternatives
    s = s.assign(abs_t=s.t_A.abs())
    alts = {"registered": s["rank"],
            "alt1_abs_t": _rank(s, s.abs_t),
            "alt2_abs_t_x_tierweight": _rank(s, s.abs_t * tc.map(ALT_WEIGHT))}
    tab = s[["rank", "direction", "context", "unit", "tier", "t_A"]].rename(columns={"rank": "rank_registered"})
    tab["rank_alt1_abs_t"] = alts["alt1_abs_t"].values
    tab["rank_alt2_abs_t_x_tierweight"] = alts["alt2_abs_t_x_tierweight"].values
    tab.sort_values("rank_registered").to_csv(REV / "r5_ranking_alternatives.tsv", sep="\t", index=False,
                                              float_format="%.6g")
    from scipy.stats import spearmanr
    rk = {}
    reg_r = s["rank"]
    for name in ("alt1_abs_t", "alt2_abs_t_x_tierweight"):
        a = alts[name]
        top = {}
        for n in (10, 30):
            top[f"top{n}_overlap"] = int(len(set(s.index[reg_r <= n]) & set(s.index[a <= n])))
        hl = {}
        for ctx, u in HIGHLIGHT + [(c, u) for c, u in zip(s[s.direction == "H1"].context, s[s.direction == "H1"].unit)]:
            m = (s.context == ctx) & (s.unit == u)
            if m.sum() == 1:
                i = s.index[m][0]
                hl[f"{ctx} {u}"] = [int(reg_r[i]), int(a[i])]
        rk[name] = dict(n_rows_rank_changed=int((a != reg_r).sum()),
                        spearman_vs_registered=round(float(spearmanr(reg_r, a)[0]), 4),
                        max_abs_rank_shift=int((a - reg_r).abs().max()),
                        median_abs_rank_shift=float((a - reg_r).abs().median()),
                        highlighted_registered_vs_alt=hl, **top)
    return dict(loso_h1=loso, ranking=dict(membership_changed=False, n=int(len(s)), weights_alt2=ALT_WEIGHT,
                                           **rk))


# ------------------------------------------------------------------------------------------------------------- A3
def a3() -> dict:
    d = pd.read_csv(REV / "r1b_random_discoveries.tsv.gz", sep="\t")
    d = d[d.discovery.astype(str).eq("True")]
    summ = json.loads((REV / "r1_summary.json").read_text())
    rows, out = [], {}
    for aid, g in d.groupby("direction"):
        b = summ["directions"][aid]["b"]
        flag = g.exceeds_all_others.astype(str).eq("True")
        assert int(flag.sum()) == b["summary"]["random_symmetric_flag"], aid
        assert len(g) == b["summary"]["n_random_discoveries_total"], aid
        per = g.assign(flag=flag).groupby("rep").agg(n_discovery=("flag", "size"), n_flag=("flag", "sum"))
        exp = {r["rep"]: r["n_discovery"] for r in b["replicates"]}
        assert all(int(per.loc[k, "n_discovery"]) == v for k, v in exp.items()), aid
        per["rate"] = per.n_flag / per.n_discovery
        for rep, r in per.iterrows():
            rows.append(dict(direction=aid, rep=int(rep), n_discovery=int(r.n_discovery), n_flag=int(r.n_flag),
                             rate=float(r.rate)))
        out[aid] = dict(per_atlas_min=round(float(per.rate.min()), 4), per_atlas_max=round(float(per.rate.max()), 4),
                        per_atlas_median=round(float(per.rate.median()), 4),
                        pooled=round(float(flag.mean()), 4), n_atlases=int(len(per)),
                        real_rate=round(float(b["summary"]["real_beats_all_rate"]), 4),
                        n_atlases_at_or_above_real=int((per.rate >= b["summary"]["real_beats_all_rate"]).sum()))
    pd.DataFrame(rows).to_csv(REV / "r5_atlas_flag_range.tsv", sep="\t", index=False, float_format="%.6g")
    return out


def main():
    res = {"standing": "POST HOC, descriptive; committed tables only; no refit, no lockbox or test-label access "
                       "(INTENT_round2.md)",
           "A1_h1_margin_005": a1(), "A2_alternative_readings": a2(), "A3_atlas_flag_range": a3()}
    (REV / "r5_summary.json").write_text(json.dumps(res, indent=1, sort_keys=True, default=str) + "\n")
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
