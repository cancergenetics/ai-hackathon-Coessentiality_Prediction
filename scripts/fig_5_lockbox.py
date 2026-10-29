"""Figure 5: lockbox confirmation (single read, R7 primary with the disclosed corrected relook), BENCH and BRIDGE-L.

Inputs (committed only; no dataset, data_cache or persisted lockbox read is opened):
  results/part2/lockbox_confirm/comparison_primary_vs_corrected_R7.tsv        rates per direction and subset
  results/part2/verify/lockbox_results_review/tables_review.json              post hoc power aid (expected rates)
  results/part2/lockbox_confirm/h1_events_primary_vs_corrected_R7.tsv         H1 per event axis, primary / corrected
  results/part2/lockbox_confirm/H1_within_complex_GI_events.tsv                H1 L context-line counts
  results/part2/lockbox_confirm/BENCH_complex_vs_gene/{bench.json, bench_pairs.tsv}
  results/part2/verify/lockbox_results_review/bench_recompute.json            BENCH CIs with replicate-assembled design
  results/part2/lockbox_confirm/BRIDGE-L_confirmation_rate/bridge_hits.tsv   per-hit confirmation by map and size
  results/part2/lockbox_confirm/primary_relabelled_R7.tsv.gz                  R7 calls (NOT EVALUABLE flags)
  results/part2/lockbox_confirm/bridge_primary_vs_corrected_R7.json           BRIDGE-L estimates and CIs

Panels
 a  confirmation rate per direction and subset: primary evaluable, primary all-in, corrected (non-confirmatory);
    grey band = post hoc power-expected rate (H3a, H3b; empirical-Bayes-shrunk A effects, L power)
 b  H1 lockbox forest: L estimate (95% CI) for the five eligible calls and the 11 A-neutral events, with the A
    estimate; FOCAD shows only the corrected relook (defect-sensitive, no claim)
 c  BENCH: mean Delta = t_L(complex) - t_L(best gene) overall and per A-breadth tercile, 95% bootstrap CI as executed
    and with the covariate design assembled within each replicate
 d  BRIDGE-L: confirmation rate per cluster-size bin for zavg-map and GLS-map discoveries (R7 primary, NOT
    EVALUABLE excluded), and the size-standardised difference in every version with 95% CI

Usage: python scripts/fig_5_lockbox.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import fig_style as S

LB = "results/part2/lockbox_confirm"
REV = "results/part2/verify/lockbox_results_review"
C_EVAL, C_ALL, C_COR = S.OI["blue"], S.OI["vermillion"], S.OI["green"]


def _tsv(rel: str) -> pd.DataFrame:
    return pd.read_csv(S.ROOT / rel, sep="\t")


def panel_a(ax) -> None:
    comp = _tsv(f"{LB}/comparison_primary_vs_corrected_R7.tsv").set_index(["direction", "subset"])
    rv = S.load_json(f"{REV}/tables_review.json")
    pw, pwb = rv["power"], rv["power_by_beats_all"]

    def eb(ID):
        v = [pw[ID][k]["expected_rate"] for k in pw[ID] if isinstance(pw[ID][k], dict) and "EB" in k]
        return min(v), max(v)

    def ebb(ID, flag):
        e = pwb[f"{ID}|beats_all={flag}"]
        return min(e["exp_eb_seL"], e["exp_eb_seA"]), max(e["exp_eb_seL"], e["exp_eb_seA"])

    rows = [
        ("H3a lineage: all", ("H3a_lineage_atlas", "overall"), eb("H3a_lineage_atlas")),
        ("   beats all matched sets", ("H3a_lineage_atlas", "beats_all_matched_sets"), ebb("H3a_lineage_atlas", "1.0")),
        ("   not beats all", ("H3a_lineage_atlas", "not_beats_all"), ebb("H3a_lineage_atlas", "0.0")),
        ("H3b driver/state: all", ("H3b_driver_family", "overall"), eb("H3b_driver_family")),
        ("   beats all matched sets", ("H3b_driver_family", "beats_all_matched_sets"), ebb("H3b_driver_family", "1.0")),
        ("   not beats all", ("H3b_driver_family", "not_beats_all"), ebb("H3b_driver_family", "0.0")),
        ("H1 eligible calls", ("H1_within_complex_GI", "overall"), None),
        ("H1S (exploratory)", ("H1S_silencing_GI", "overall"), None),
        ("H2 (exploratory)", ("H2_submodule_divergence", "overall"), None),
        ("Controls (units)", ("controls (23)", "overall"), None),
        ("BRIDGE-L GLS-map hits", ("BRIDGE-L_confirmation_rate (GLS hits)", "overall"), None),
    ]
    y = np.arange(len(rows))[::-1]
    for yy, (lab, key, band) in zip(y, rows):
        r = comp.loc[key]
        if band is not None:
            ax.add_patch(plt.Rectangle((band[0], yy - 0.32), band[1] - band[0], 0.64, color="#CFCFCF", lw=0, zorder=0))
        ax.plot([r.primary_rate_all_in, r.corrected_rate_evaluable, r.primary_rate_evaluable], [yy] * 3,
                color="#888888", lw=0.5, zorder=1)
        ax.plot(r.primary_rate_evaluable, yy, "o", color=C_EVAL, ms=4, zorder=3)
        ax.plot(r.primary_rate_all_in, yy, "s", mfc="white", mec=C_ALL, mew=0.9, ms=3.6, zorder=3)
        ax.plot(r.corrected_rate_evaluable, yy, "D", mfc="white", mec=C_COR, mew=0.9, ms=3.2, zorder=4)
        ne = int(r.primary_n_not_evaluable)
        txt = f"{int(r.primary_n_confirmed)}/{int(r.primary_n_lockbox) - ne}"
        if ne:
            txt += f" +{ne} NE\ncorr. {int(r.corrected_n_confirmed)}/{int(r.corrected_n_lockbox) - int(r.corrected_n_not_evaluable)}"
        ax.text(1.01, yy, txt, va="center", fontsize=4.4, transform=ax.get_yaxis_transform(), linespacing=0.95)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=5.4)
    ax.set_xlim(0.25, 1.0)
    ax.set_xlabel("Lockbox confirmation rate")
    ax.set_title("Confirmation in L (338 lines)", loc="left")
    hs = [plt.Line2D([], [], marker="o", ls="", color=C_EVAL, ms=4, label="primary, evaluable"),
          plt.Line2D([], [], marker="s", ls="", mfc="white", mec=C_ALL, ms=3.6, label="primary, all-in"),
          plt.Line2D([], [], marker="D", ls="", mfc="white", mec=C_COR, ms=3.2, label="corrected (non-confirmatory)"),
          plt.Rectangle((0, 0), 1, 1, color="#CFCFCF", label="expected from L power (post hoc)")]
    ax.legend(handles=hs, frameon=False, fontsize=4.8, loc="upper left", bbox_to_anchor=(-0.05, -0.13), ncol=2)


def panel_b(ax) -> None:
    ev = _tsv(f"{LB}/h1_events_primary_vs_corrected_R7.tsv")
    nl = _tsv(f"{LB}/H1_within_complex_GI_events.tsv")[["event", "axis", "n_context_L"]]
    ev = ev.merge(nl, on=["event", "axis"], how="left")
    for col, src in (("p", "primary_ci95"), ("c", "corrected_ci95")):
        ci = ev[src].str.strip("[]").str.split(",", expand=True).apply(pd.to_numeric, errors="coerce")
        ev[f"{col}_lo"], ev[f"{col}_hi"] = ci[0], ci[1]
    elig = ev[ev.eligible].copy()
    order = ["Z0109|ARID1A", "Z0186|EP300", "Z0286|FOCAD", "Z0158|TP53", "Z0130|PTEN"]
    elig["o"] = elig.event.map({e: i for i, e in enumerate(order)})
    elig = elig.sort_values("o")
    neu = ev[(ev.label_final == "NEUTRAL") & (ev.axis == "N")].sort_values("primary_beta_L")
    rows = list(elig.itertuples()) + list(neu.itertuples())
    y = np.arange(len(rows))[::-1].astype(float)
    y[len(elig):] -= 0.8
    XL = (-1.0, 0.55)
    ax.axvspan(-0.1, 0.1, color="#EEEEEE", zorder=0)
    ax.axvline(0, color="black", lw=0.5)
    labels = []
    for yy, r in zip(y, rows):
        gene = r.event.split("|")[1]
        defect = r.primary_call.startswith("NOT EVALUABLE")
        if defect:
            b, lo, hi = r.corrected_beta_L, r.c_lo, r.c_hi
        else:
            b, lo, hi = r.primary_beta_L, r.p_lo, r.p_hi
        call = r.primary_call
        if defect:
            col, mfc, tag = "#9E9E9E", "white", "NE (defect); corrected confirms — no claim"
        elif call == "confirmed" and r.promotable_primary:
            col, mfc, tag = S.OI["blue"], S.OI["blue"], "confirmed"
        elif call == "confirmed":
            col, mfc, tag = S.OI["orange"], S.OI["orange"], "confirmed, within margin (E2): not promotable"
        elif "rare-context route passed" in call:
            col, mfc, tag = S.OI["purple"], "white", "rare-context route passed (not a replication)"
        elif call == "not confirmed":
            col, mfc, tag = S.OI["vermillion"], S.OI["vermillion"], "not confirmed"
        else:
            col, mfc, tag = S.OI["green"], S.OI["green"], ""
        if pd.notna(b):
            ax.plot([max(lo, XL[0]), min(hi, XL[1])], [yy, yy], color=col, lw=0.9, zorder=2)
            if lo < XL[0]:
                ax.plot(XL[0] + 0.02, yy, marker="<", color=col, ms=2.6, zorder=3)
            if hi > XL[1]:
                ax.plot(XL[1] - 0.02, yy, marker=">", color=col, ms=2.6, zorder=3)
            ax.plot(b, yy, "o", color=col, mfc=mfc, ms=3.2, zorder=4)
        ax.plot(r.beta_A, yy, marker="x", color="black", ms=3.2, mew=0.8, zorder=5)
        nL = "" if pd.isna(r.n_context_L) else f", L n={int(r.n_context_L)}"
        labels.append(f"{gene} {r.axis}{nL}")
        if tag:
            ax.text(XL[1] + 0.02, yy, tag, fontsize=4.4, va="center", color=col)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=4.8)
    ax.set_xlim(*XL)
    ax.set_ylim(y.min() - 0.7, y.max() + 0.7)
    ax.text(XL[0] + 0.02, y[len(elig)] + 0.55, "A-neutral events (TOST ±0.1 in A; not eligible)", fontsize=4.6,
            color="#555555", va="bottom")
    ax.set_xlabel("β_P / κ_N in L (Chronos; 95% CI)    × = A estimate")
    ax.set_title("H1 in the lockbox", loc="left")


def panel_c(ax) -> None:
    b = S.load_json(f"{LB}/BENCH_complex_vs_gene/bench.json")
    br = S.load_json(f"{REV}/bench_recompute.json")
    pairs = _tsv(f"{LB}/BENCH_complex_vs_gene/bench_pairs.tsv")
    pr = b["primary"]
    cats = [("", "all pairs", None), ("_low", "low", "low"), ("_mid", "mid", "mid"), ("_high", "high\n(registered)", "high")]
    ax.axhline(0, color="black", lw=0.5)
    for i, (k, lab, terc) in enumerate(cats):
        m = pr[f"mean{k}"]
        lo, hi = pr[f"ci95{k}"]
        alo, ahi = br["bootstrap_assembled"][f"ci95{k}"]
        ax.errorbar(i - 0.12, m, yerr=[[m - lo], [hi - m]], fmt="o", color=S.OI["blue"], ms=3.6, lw=0.9, capsize=1.5)
        ax.errorbar(i + 0.12, m, yerr=[[m - alo], [ahi - m]], fmt="o", color=S.OI["orange"], mfc="white", ms=3.4,
                    lw=0.9, capsize=1.5)
        sub = pairs if terc is None else pairs[pairs.tercile == terc]
        ax.text(i, -0.66, f"n={len(sub)}\n{(sub.delta > 0).mean():.0%} Δ>0", ha="center", va="bottom", fontsize=4.6,
                color="#555555")
    ax.set_xticks(range(len(cats)))
    ax.set_xticklabels([c[1] for c in cats], fontsize=5.6)
    ax.set_xlabel("A breadth tercile (fraction of members with one-sided p < 0.05 in A)")
    ax.set_ylabel("Δ = t_L(complex) − t_L(best gene)")
    ax.set_ylim(-0.68, 0.62)
    ax.axvspan(2.6, 3.4, color="#F2F2F2", zorder=0)
    hs = [plt.Line2D([], [], marker="o", color=S.OI["blue"], ms=3.4, label="95% CI as executed"),
          plt.Line2D([], [], marker="o", color=S.OI["orange"], mfc="white", ms=3.4,
                     label="95% CI, covariates assembled per replicate")]
    ax.legend(handles=hs, frameon=False, fontsize=4.8, loc="upper left")
    ax.text(0.02, 0.85, f"{br['bootstrap_executed']['frac_pair_rep_nan']:.1%} of pair × replicate fits fail\n"
                        f"(dropped); {b['distinct']['n_pairs']:,} pairs in {b['distinct']['n_distinct_complexes_in_pairs']} complexes",
            transform=ax.transAxes, ha="left", va="top", fontsize=4.6, color="#555555")
    ax.set_title("BENCH: complex vs best gene", loc="left")


def panel_d(ax_bins, ax_diff) -> None:
    h = _tsv(f"{LB}/BRIDGE-L_confirmation_rate/bridge_hits.tsv")
    rel = _tsv(f"{LB}/primary_relabelled_R7.tsv.gz")
    h["entry"] = np.where(h["map"].eq("gls"), "BRIDGE-L_confirmation_rate", h["source"])
    x = h.merge(rel[["entry", "context", "unit", "primary_call"]], on=["entry", "context", "unit"], how="left")
    assert x.primary_call.notna().all()
    x = x[~x.primary_call.str.startswith("NOT EVALUABLE")].copy()
    x["bin"] = pd.cut(x["size"], [2, 4, 8, 16, 60], labels=["3–4", "5–8", "9–16", "17–60"])
    r = x.groupby(["map", "bin"], observed=True).confirmed.agg(["mean", "size"])
    bins = ["3–4", "5–8", "9–16", "17–60"]
    for j, (mp, col, lab) in enumerate((("zavg", S.COLOR["zavg"], "zavg_BEF map"), ("gls", S.COLOR["gls"], "GLS map"))):
        xs = np.arange(len(bins)) + (j - 0.5) * 0.38
        vals = [r.loc[(mp, b), "mean"] for b in bins]
        ax_bins.bar(xs, vals, width=0.36, color=col, label=lab)
        for xx, b, v in zip(xs, bins, vals):
            ax_bins.text(xx, v + 0.01, str(int(r.loc[(mp, b), "size"])), ha="center", va="bottom", fontsize=4.2,
                         rotation=90)
    ax_bins.set_xticks(range(len(bins)))
    ax_bins.set_xticklabels(bins, fontsize=5.0, rotation=45)
    ax_bins.set_xlabel("Cluster size (members)")
    ax_bins.set_ylabel("Confirmation rate (R7 primary)")
    ax_bins.set_ylim(0, 0.9)
    ax_bins.legend(frameon=False, fontsize=4.6, loc="upper left", ncol=1)
    ax_bins.set_title("BRIDGE-L: zavg vs GLS", loc="left")

    bj = S.load_json(f"{LB}/bridge_primary_vs_corrected_R7.json")
    vers = [("R7 primary", bj["primary_evaluable"]["primary"], True),
            ("as executed", bj["as_executed"]["primary"], False),
            ("corrected", bj["corrected_nonconfirmatory"]["primary"], False),
            ("unweighted", bj["primary_evaluable"]["multiverse"]["unweighted"], False),
            ("matched cov.", bj["primary_evaluable"]["multiverse"]["matched_gene_coverage"], False),
            ("T1/T2 only", bj["primary_evaluable"]["multiverse"]["T1_T2_clusters"], False)]
    y = np.arange(len(vers))[::-1]
    ax_diff.axvline(0, color="black", lw=0.5)
    for yy, (lab, v, head) in zip(y, vers):
        ax_diff.errorbar(v["diff"], yy, xerr=[[v["diff"] - v["ci95"][0]], [v["ci95"][1] - v["diff"]]], fmt="o",
                         color="black" if head else "#777777", ms=3.4 if head else 2.8, lw=1.0 if head else 0.7,
                         capsize=1.2)
    ax_diff.set_yticks(y)
    ax_diff.set_yticklabels([v[0] for v in vers], fontsize=5.0)
    ax_diff.set_xlim(-0.14, 0.03)
    ax_diff.set_xticks([-0.1, 0.0])
    ax_diff.set_xlabel("Δ rate zavg − GLS (95% CI)")
    ax_diff.set_title("size-standardised", loc="left", fontsize=6.5)


def main() -> None:
    S.setup()
    fig = plt.figure(figsize=(7.2, 7.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1.0], width_ratios=[1.0, 1.0], hspace=0.5, wspace=0.95)
    ax = fig.add_subplot(gs[0, 0])
    S.panel_label(ax, "a")
    panel_a(ax)
    ax = fig.add_subplot(gs[0, 1])
    S.panel_label(ax, "b")
    panel_b(ax)
    ax = fig.add_subplot(gs[1, 0])
    S.panel_label(ax, "c")
    panel_c(ax)
    sub = gs[1, 1].subgridspec(1, 2, width_ratios=[1.0, 0.75], wspace=1.6)
    a1, a2 = fig.add_subplot(sub[0, 0]), fig.add_subplot(sub[0, 1])
    S.panel_label(a1, "d")
    panel_d(a1, a2)
    S.save(fig, "fig5_lockbox")


if __name__ == "__main__":
    main()
