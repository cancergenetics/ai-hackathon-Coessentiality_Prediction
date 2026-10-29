"""Figure 1: the single confirmatory look of the locked method on the withheld hu.MAP test pairs.

Input: results/test_eval/test_results.json only (committed in 0861e59; written by scripts/21_evaluate_test.py).
No label or score file is read. Full precision-recall curves would need the test labels joined to the score TSVs,
neither of which is committed, so panel a shows the protocol's pre-specified operating points instead.

Panels
 a  precision at recall 0.10 / 0.25 / 0.50 (ITT) and the rule (b) operating point with cluster-bootstrap CIs
 b  paired Delta AP with 95% cluster-bootstrap CIs (B = 2,000, seed 0), ITT and S1
 c  isotonic reliability (train-fitted calibrator, log-spaced bins, ITT)
 d  pooled vs cluster-structured AP (ITT): the macro-AP caveat
 e  AP by protein-novelty stratum C1/C2/C3 with Delta AP vs GLS

Usage: python scripts/fig_1_test_evaluation.py
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

import fig_style as S


def main() -> None:
    S.setup()
    d = S.load_json(S.TEST_RESULTS)
    R = d["results"]
    itt, s1 = R["itt"], R["mapped_only"]
    shown = ["zavg_BEF", "B", "gls_signed", "spearman_signed", "pearson_signed", "spec_k20_tinf_std"]

    fig = plt.figure(figsize=(7.2, 7.4))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 0.9], hspace=0.55, wspace=0.42)

    # ---------- a: operating points
    ax = fig.add_subplot(gs[0, 0])
    S.panel_label(ax, "a")
    stat = itt["bootstrap"]["cluster"]["stat"]
    for m in shown:
        mm = itt["metrics"][m]
        rec = [0.10, 0.25, 0.50]
        prec = [mm["P@R0.1"], mm["P@R0.25"], mm["P@R0.5"]]
        ax.plot(rec, prec, ls=":", lw=0.8, color=S.COLOR[m], marker=S.MARKER[m], ms=3, mfc="white")
        p, r = stat[f"{m}|precision_rule_b"], stat[f"{m}|recall_rule_b"]
        ax.errorbar(r["point"], p["point"], xerr=[[r["point"] - r["ci_low"]], [r["ci_high"] - r["point"]]],
                    yerr=[[p["point"] - p["ci_low"]], [p["ci_high"] - p["point"]]], fmt=S.MARKER[m],
                    color=S.COLOR[m], ms=4.5, lw=0.7, capsize=1.5, label=S.SHORT[m])
    prev = itt["prevalence_design"]
    ax.axhline(prev, color=S.OI["grey"], lw=0.6, ls="--")
    ax.text(0.01, prev + 0.015, f"prevalence {prev:.4f}", fontsize=5.5, color=S.OI["grey"], ha="left")
    ax.set_xlim(0, 0.55)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("Operating points (ITT)", loc="left")
    ax.legend(frameon=False, loc="upper right", handletextpad=0.3, borderaxespad=0.2)

    # ---------- b: Delta AP forest
    ax = fig.add_subplot(gs[0, 1])
    S.panel_label(ax, "b")
    rows = []
    for cand in ["zavg_BEF", "B"]:
        for b in S.BASELINES:
            rows.append((f"{S.SHORT[cand]} − {S.SHORT[b]}", f"{cand}|ap - {b}|ap", cand))
    rows.append(("B − zavg_BEF", "B|ap - zavg_BEF|ap", "B"))
    rows.append(("PC-rm − zavg_BEF", "spec_k20_tinf_std|ap - zavg_BEF|ap", "spec_k20_tinf_std"))
    y = np.arange(len(rows))[::-1].astype(float)
    for (lab, key, cand), yy in zip(rows, y):
        for sset, off, mk, name in ((itt, 0.17, "o", "ITT"), (s1, -0.17, "s", "S1 (mapped-only)")):
            v = sset["bootstrap"]["cluster"]["diff"][key]
            ax.errorbar(v["diff"], yy + off, xerr=[[v["diff"] - v["ci_low"]], [v["ci_high"] - v["diff"]]],
                        fmt=mk, color=S.COLOR[cand], mfc=S.COLOR[cand] if name == "ITT" else "white",
                        ms=3.5, lw=0.8, capsize=1.5)
    ax.axvline(0, color="black", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows])
    ax.set_xlabel("Δ AP (95% cluster-bootstrap CI)")
    ax.set_title("Paired differences in AP", loc="left")
    h1 = plt.Line2D([], [], marker="o", color="k", ls="", ms=3.5, label="ITT (primary metric)")
    h2 = plt.Line2D([], [], marker="s", color="k", mfc="white", ls="", ms=3.5, label="S1 (locked verdict set)")
    ax.legend(handles=[h1, h2], frameon=False, loc="lower right", fontsize=5.5)

    # ---------- c: reliability
    ax = fig.add_subplot(gs[1, 0])
    S.panel_label(ax, "c")
    lo, hi = 1e-4, 1.0
    ax.plot([lo, hi], [lo, hi], color=S.OI["grey"], lw=0.6, ls="--")
    n_zero = 0
    for m in ["zavg_BEF", "gls_signed", "spearman_signed", "pearson_signed"]:
        rel = itt["metrics"][m]["calibration"]["isotonic"]["reliability_log16"]
        xs = np.array([b["mean_pred"] for b in rel])
        ys = np.array([b["obs_rate"] for b in rel])
        ns = np.array([b["n"] for b in rel])
        keep = ys > 0
        n_zero += int((~keep).sum())
        ax.plot(xs[keep], ys[keep], marker=S.MARKER[m], color=S.COLOR[m], lw=0.8, ms=3, label=S.SHORT[m])
        del ns
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e-3, 1)
    ax.set_ylim(1e-3, 1)
    ax.set_xlabel("Mean predicted probability (isotonic, fitted on train)")
    ax.set_ylabel("Observed co-complex rate")
    ax.set_title("Calibration (ITT, log-spaced bins)", loc="left")
    iso = itt["metrics"]["zavg_BEF"]["calibration"]["isotonic"]
    txt = (f"zavg_BEF: Brier skill {iso['bss_vs_fit_prevalence']:.3f}\n"
           f"calibration-in-the-large {iso['calibration_in_the_large']:.3f}\n"
           f"slope {iso['calibration_slope']:.3f}")
    ax.text(0.03, 0.97, txt, transform=ax.transAxes, fontsize=5.5, va="top")
    if n_zero:
        ax.text(0.98, 0.02, f"{n_zero} bins with no positives not shown", transform=ax.transAxes,
                fontsize=5, ha="right", va="bottom", color=S.OI["grey"])
    ax.legend(frameon=False, loc="lower right", bbox_to_anchor=(1.0, 0.07))

    # ---------- d: pooled vs cluster-structured AP
    ax = fig.add_subplot(gs[1, 1])
    S.panel_label(ax, "d")
    groups = [("AP\n(pooled)", "ap"), ("AP without\nlargest comp.", "ap_no_largest_component"),
              ("AP without\n10 largest", "ap_no_top10_clusters"), ("macro\nper-cluster AP", "macro_cluster_ap")]
    meths = ["zavg_BEF", "gls_signed", "spearman_signed", "pearson_signed"]
    w = 0.19
    x0 = np.arange(len(groups))
    for i, m in enumerate(meths):
        for j, (_, key) in enumerate(groups):
            val = itt["metrics"][m][key]
            st = stat.get(f"{m}|{key}")
            xx = x0[j] + (i - 1.5) * w
            ax.bar(xx, val, width=w, color=S.COLOR[m], label=S.SHORT[m] if j == 0 else None)
            if st is not None:
                ax.errorbar(xx, st["point"], yerr=[[st["point"] - st["ci_low"]], [st["ci_high"] - st["point"]]],
                            color="black", lw=0.6, capsize=1)
    ax.set_xticks(x0)
    ax.set_xticklabels([g[0] for g in groups], fontsize=6)
    ax.set_ylabel("AP (ITT)")
    ax.set_title("Pooled vs cluster-structured AP (ITT)", loc="left")
    dm = itt["bootstrap"]["cluster"]["diff"]["zavg_BEF|macro_cluster_ap - gls_signed|macro_cluster_ap"]
    ax.text(0.99, 0.86, "zavg_BEF − GLS, macro-AP:\n" + S.fmt_ci(dm["diff"], dm["ci_low"], dm["ci_high"], 3, True),
            transform=ax.transAxes, ha="right", va="top", fontsize=5.5)
    ax.legend(frameon=False, loc="upper left", fontsize=5.2, ncol=4, columnspacing=0.8, handlelength=1.2)
    ax.set_ylim(0, 0.6)

    # ---------- e: C1/C2/C3 strata
    ax = fig.add_subplot(gs[2, :])
    S.panel_label(ax, "e")
    strata = [("C1: both proteins\nin train positives", "stratum_C1"),
              ("C2: one protein", "stratum_C2"), ("C3: neither protein", "stratum_C3"),
              ("all (ITT)", "itt")]
    meths = ["zavg_BEF", "B", "gls_signed", "spearman_signed", "pearson_signed"]
    w = 0.15
    x0 = np.arange(len(strata))
    for i, m in enumerate(meths):
        for j, (_, key) in enumerate(strata):
            st = R[key]["bootstrap"]["cluster"]["stat"][f"{m}|ap"]
            xx = x0[j] + (i - 2) * w
            ax.bar(xx, st["point"], width=w, color=S.COLOR[m], label=S.SHORT[m] if j == 0 else None)
            ax.errorbar(xx, st["point"], yerr=[[st["point"] - st["ci_low"]], [st["ci_high"] - st["point"]]],
                        color="black", lw=0.6, capsize=1)
    ticks = []
    for j, (name, key) in enumerate(strata):
        r = R[key]
        dv = r["bootstrap"]["cluster"]["diff"]["zavg_BEF|ap - gls_signed|ap"]
        ticks.append(f"{name.replace(chr(10), ' ')}\npositives = {r['n_pos']:,}\nΔ AP vs GLS:\n"
                     + S.fmt_ci(dv["diff"], dv["ci_low"], dv["ci_high"], 3, True))
    ax.set_xticks(x0)
    ax.set_xticklabels(ticks, fontsize=5.8)
    ax.set_ylabel("AP (95% cluster-bootstrap CI)")
    ax.set_title("Generalisation to proteins not seen among train positives", loc="left")
    ax.legend(frameon=False, ncol=5, loc="upper left", fontsize=5.5)
    ax.set_ylim(0, max(ax.get_ylim()[1], 0.6))

    S.save(fig, "fig1_test_evaluation")


if __name__ == "__main__":
    main()
