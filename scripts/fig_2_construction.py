"""Figure 2: complex construction on the construction set A, tiers, map-to-map variability and the bridge endpoint.

Inputs (committed only): results/part2/construction/{summary/summary.json, selection/selection.json,
tiers/summary.json, variability/variability.json, variability/*_per_complex_jaccard.tsv, bridge/bridge.json}.
This is the pre-written NULL version of Fig 2 (DESIGN v3.1 §4): the pre-registered bridge endpoint's CI includes 0.

Panels
 a  level-1 cluster size distribution, zavg_BEF vs GLS map
 b  coverage of S genes by A-essentiality class
 c  split-half selection grid (z0, k): the criterion is flat across z0 at k = 10
 d  tiers (T1 BioPlex, T2 protein co-abundance, T3-neg, T3-nc)
 e  map-to-map variability: per-complex median best-match Jaccard in 20 subsample (80%, PatientID-blocked) maps vs random size-matched sets
 f  bridge forest plot: primary endpoint (null), sensitivities, secondary endpoints

Usage: python scripts/fig_2_construction.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import fig_style as S

BINS = ["3-4", "5-8", "9-16", "17-60"]


def main() -> None:
    S.setup()
    summ = S.load_json(S.SUMMARY)
    sel = S.load_json(S.SELECTION)
    tiers = S.load_json(S.TIERS)
    var = S.load_json(S.VARIAB)
    br = S.load_json(S.BRIDGE)
    jz = pd.read_csv(S.ROOT / S.JAC_ZAVG, sep="\t")
    jg = pd.read_csv(S.ROOT / S.JAC_GLS, sep="\t")

    fig = plt.figure(figsize=(7.2, 7.6))
    gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 0.9], hspace=0.75, wspace=0.75)
    cz, cg = S.COLOR["zavg"], S.COLOR["gls"]

    # ---------- a: size distribution
    ax = fig.add_subplot(gs[0, 0])
    S.panel_label(ax, "a")
    x = np.arange(len(BINS))
    w = 0.38
    for off, key, col, lab in ((-w / 2, "zavg", cz, "zavg_BEF"), (w / 2, "gls", cg, "GLS")):
        sb = summ[key]["size_bins_level1"]
        ax.bar(x + off, [sb[b] for b in BINS], width=w, color=col,
               label=f"{lab}: {summ[key]['n_level1']} complexes, {summ[key]['genes_level1']:,} genes")
    ax.set_xticks(x)
    ax.set_xticklabels(BINS)
    ax.set_xlabel("Level-1 cluster size (genes)")
    ax.set_ylabel("Clusters")
    ax.set_title("Map overview", loc="left")
    ax.legend(frameon=False, fontsize=5, loc="upper left", bbox_to_anchor=(0, 1.02))
    ax.set_ylim(0, max(summ["gls"]["size_bins_level1"].values()) * 1.45)

    # ---------- b: coverage by class
    ax = fig.add_subplot(gs[0, 1])
    S.panel_label(ax, "b")
    classes = [("pan", "pan-\nessential"), ("selective", "selective"), ("other", "other")]
    x = np.arange(len(classes))
    for off, key, col in ((-w / 2, "zavg", cz), (w / 2, "gls", cg)):
        cov = summ[key]["coverage_by_class"]
        ax.bar(x + off, [cov[c]["frac"] for c, _ in classes], width=w, color=col)
        for xi, (c, _) in zip(x, classes):
            ax.text(xi + off, cov[c]["frac"] + 0.02, f"{cov[c]['clustered']:,}", ha="center", fontsize=4.8, rotation=90)
    ax.set_xticks(x)
    ax.set_xticklabels([c[1] for c in classes], fontsize=5.8)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel(f"Fraction of S genes clustered\n(S = {summ['zavg']['genes_S']:,} genes)")
    ax.set_title("Coverage", loc="left")

    # ---------- c: selection grid
    ax = fig.add_subplot(gs[0, 2])
    S.panel_label(ax, "c")
    for key, col, mk, lab in (("zavg_BEF", cz, "o", "zavg_BEF"), ("gls", cg, "D", "GLS (edge-matched)")):
        grid = sel[key]["grid"]
        for k, ls in ((10, "-"), (20, ":")):
            pts = sorted([g for g in grid if g["k"] == k], key=lambda g: g["z0"])
            ax.errorbar([g["z0"] for g in pts], [g["stat"] for g in pts], yerr=[g["stat_sd"] for g in pts],
                        color=col, marker=mk, ms=3, ls=ls, lw=0.8, capsize=1.5, label=f"{lab}, k = {k}")
        rew = max(g["rewired_stat"] for g in grid)
    ax.axhline(rew, color=S.OI["grey"], lw=0.6, ls="--")
    ax.text(6.0, rew + 0.01, f"rewired baseline = {rew:g}", fontsize=5, color=S.OI["grey"], ha="right")
    s0 = sel["zavg_BEF"]["selected"]
    ax.plot([s0["z0"]], [s0["stat"]], marker="*", ms=8, color="black", ls="", zorder=5)
    ax.annotate(f"selected z0 = {s0['z0']:g}, k = {s0['k']}", (s0["z0"], s0["stat"]), xytext=(4.0, 0.42),
                fontsize=5.2, arrowprops=dict(arrowstyle="-", lw=0.5))
    ax.set_xticks([4, 5, 6])
    ax.set_xlabel("Edge threshold z0")
    ax.set_ylabel("Split-half recovery statistic\n(mean ± SD over 4 splits)")
    ax.set_ylim(0, 0.5)
    ax.set_title("Selection", loc="left")
    ax.legend(frameon=False, fontsize=4.5, loc="lower left", bbox_to_anchor=(0.0, 0.04), handlelength=1.5)

    # ---------- d: tiers
    ax = fig.add_subplot(gs[1, 0:2])
    S.panel_label(ax, "d")
    rows = [("zavg_BEF L1", tiers["tiers"]["zavg"]["level1"]),
            ("zavg_BEF L1 non-locus", summ["zavg"]["tiers_nonlocus_level1"]),
            ("zavg_BEF L2", tiers["tiers"]["zavg"]["level2"]),
            ("GLS L1", tiers["tiers"]["gls"]["level1"]),
            ("GLS L1 non-locus", summ["gls"]["tiers_nonlocus_level1"]),
            ("GLS L2", tiers["tiers"]["gls"]["level2"])]
    y = np.arange(len(rows))[::-1]
    for (lab, t), yy in zip(rows, y):
        left = 0
        tot = sum(t[k] for k in ["T1", "T2", "T3-neg", "T3-nc"])
        for k in ["T1", "T2", "T3-neg", "T3-nc"]:
            ax.barh(yy, t[k] / tot, left=left, color=S.TIER_COLOR[k], height=0.7,
                    label=k if yy == y[0] else None, edgecolor="white", lw=0.3)
            if t[k] / tot > 0.06:
                ax.text(left + t[k] / tot / 2, yy, str(t[k]), ha="center", va="center", fontsize=5.5,
                        color="white" if k in ("T1",) else "black")
            left += t[k] / tot
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r[0]} (n = {sum(r[1][k] for k in ['T1', 'T2', 'T3-neg', 'T3-nc'])})" for r in rows])
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of clusters")
    ax.set_title("Tiers", loc="left")
    ax.legend(frameon=False, ncol=4, loc="lower left", bbox_to_anchor=(0, -0.42), fontsize=5.5)

    # ---------- e: variability
    ax = fig.add_subplot(gs[1, 2])
    S.panel_label(ax, "e")
    pos = np.arange(len(BINS))
    for off, df, key, col in ((-0.2, jz, "zavg_BEF", cz), (0.2, jg, "gls", cg)):
        data = [df.loc[df["size_bin"] == b, "jaccard_median"].to_numpy() for b in BINS]
        bp = ax.boxplot(data, positions=pos + off, widths=0.32, patch_artist=True, showfliers=False,
                        medianprops=dict(color="black", lw=0.8), whiskerprops=dict(lw=0.6),
                        capprops=dict(lw=0.6), boxprops=dict(lw=0.6))
        for b in bp["boxes"]:
            b.set_facecolor(col)
            b.set_alpha(0.85)
        p95 = [var[key]["size_bins"][b]["random_S_p95"] for b in BINS]
        ax.plot(pos + off, p95, ls="", marker="_", ms=7, mew=1.2, color="black")
    ax.set_xticks(pos)
    ax.set_xticklabels(BINS)
    ax.set_ylabel("Median best-match Jaccard\nin 20 subsample maps")
    ax.set_ylim(0, 1.02)
    ax.set_title("Variability", loc="left")
    az = var["zavg_BEF"]["pairwise"]["ari_both"]["median"]
    ag = var["gls"]["pairwise"]["ari_both"]["median"]
    ax.set_xlabel(f"Complex size bin\npairwise ARI (median, 190 map pairs):\nzavg_BEF {az:.3f}, GLS {ag:.3f}", fontsize=6)

    # ---------- f: bridge forest
    ax = fig.add_subplot(gs[2, :])
    S.panel_label(ax, "f")
    pr = br["primary"]
    sens = br["sensitivity"]
    sec = br["secondary"]["protein_coabundance"]
    rows = [
        ("PRIMARY: BioPlex precision", pr, "black"),
        ("no coverage matching", sens["no_coverage_matching"], S.OI["grey"]),
        ("locus clusters excluded", sens["locus_excluded"], S.OI["grey"]),
        ("secondary: co-abundance r", sec, S.OI["vermillion"]),
    ]
    y = np.arange(len(rows))[::-1]
    for (lab, v, col), yy in zip(rows, y):
        lo, hi = v["ci95"]
        ax.errorbar(v["diff"], yy, xerr=[[v["diff"] - lo], [hi - v["diff"]]], fmt="o" if col == "black" else "s",
                    color=col, ms=4 if col == "black" else 3.2, lw=1.0 if col == "black" else 0.8, capsize=2)
        ax.text(0.085, yy, f"zavg {v['zavg']:.3f} vs GLS {v['gls']:.3f};  Δ " + S.fmt_ci(v["diff"], lo, hi, 3, True),
                va="center", fontsize=5.5)
    ax.axvline(0, color="black", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=6)
    ax.set_xlim(-0.13, 0.26)
    ax.set_xticks([-0.1, -0.05, 0, 0.05, 0.1])
    ax.set_xlabel("Difference zavg_BEF map − GLS map (95% protein-node bootstrap CI; ← favours GLS | favours zavg_BEF →)")
    ax.set_title(f"Bridge endpoint: {br['decision']} (the 95% CI of the primary endpoint includes 0)", loc="left")
    pa = br["secondary"]["pair_level_AP_bioplex"]
    ax.text(-0.125, y[-1] - 0.85,
            f"pair level, S × S ({pa['n_pairs']:,} pairs): AP vs BioPlex zavg_BEF {pa['zavg_BEF']['AP']:.4f}, "
            f"GLS {pa['gls']['AP']:.4f}, Pearson {pa['pearson']['AP']:.4f} (prevalence {pa['prevalence']:.4f})",
            fontsize=5.3, va="center", color=S.OI["blue"])
    ax.set_ylim(y[-1] - 1.3, y[0] + 0.6)

    S.save(fig, "fig2_construction")


if __name__ == "__main__":
    main()
