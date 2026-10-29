"""Figure 3: context atlas of the A-side H3a (lineage) and H3b (driver / genome-state) discoveries, read against the
registered calibrations run after discovery (D8), plus the exploratory H2 rank-against-splits panel.

Inputs (committed only): results/part2/directions/{H3a_lineage_atlas,H3b_driver_family}/{discoveries.tsv,
discoveries_report.tsv, calibration_D8.json, calibration_D8_tests.tsv.gz, gpd_sensitivity.json,
sensitivity_screen_quality/{discoveries.tsv, comparison_per_context.tsv, scale_effect_rho.tsv}},
results/part2/directions/H2_submodule_divergence/{discoveries.tsv, calibration_D8_tests.tsv.gz},
results/part2/fzx/zavg_engine_membership.tsv (member symbols for row labels; gene content only).

Panels
 a  fraction of tests with p < 0.05, real complexes vs their matched random gene sets, per context
 b  lineage 'scale effect': Spearman rho between unit essentiality and the lineage beta, primary vs sensitivity
 c  H2 (exploratory): rank of each discovery's |t| among its own 20 random member splits
 d  heatmap of complex-specific non-locus discoveries (beta) for the units with the strongest such hits
    (coloured: beats all matched sets and not locus; grey dot: other discovery of that unit in that context)
 e  discoveries per context: complex-specific non-locus / other non-locus / locus, with the never-expressed-set
    rate per lineage context (H3a)
 f  survival of discoveries under the screen-quality and GPD-tail sensitivities, per context

Usage: python scripts/fig_3_context_atlas.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

import fig_style as S
import part2_report as P

C_FOCUS, C_OTHER, C_LOCUS = S.OI["blue"], S.OI["skyblue"], S.OI["grey"]


def main() -> None:
    S.setup()
    mem = P.members()
    da, db = P.discoveries(P.H3A), P.discoveries(P.H3B)
    ma, mb = P.matched_rates(P.H3A), P.matched_rates(P.H3B)
    nea = P.never_expressed(P.H3A)
    cal_a = S.load_json(f"{P.DIR}/{P.H3A}/calibration_D8.json")["matched_random_sets"]
    cal_b = S.load_json(f"{P.DIR}/{P.H3B}/calibration_D8.json")["matched_random_sets"]

    fig = plt.figure(figsize=(7.2, 10.2))
    gs = fig.add_gridspec(3, 3, height_ratios=[1.0, 2.3, 1.05], hspace=0.5, wspace=0.95)

    # ---------- a: real vs matched random rate
    ax = fig.add_subplot(gs[0, 0])
    S.panel_label(ax, "a")
    for m, col, mk, lab, cal in ((ma, S.OI["blue"], "o", "H3a lineage", cal_a),
                                 (mb, S.OI["vermillion"], "s", "H3b driver / state", cal_b)):
        ax.scatter(m.random, m.real, s=9, color=col, marker=mk, lw=0, alpha=0.85,
                   label=f"{lab.split()[0]}: {cal['real_pooled']['frac_p_lt_005']:.3f} vs "
                         f"{cal['pooled']['frac_p_lt_005']:.3f}")
    ax.plot([0, 0.65], [0, 0.65], color="black", lw=0.6, ls="--")
    ax.set_xlim(0, 0.65)
    ax.set_ylim(0, 0.65)
    ax.set_xlabel("Matched random gene sets\n(fraction p < 0.05)")
    ax.set_ylabel("Real complexes (fraction p < 0.05)")
    ax.set_title("Real vs matched sets", loc="left")
    ax.legend(frameon=False, fontsize=4.8, loc="lower right", title="pooled: real vs random", title_fontsize=4.8)

    # ---------- e: lineage scale effect (rho)
    ax = fig.add_subplot(gs[0, 1])
    S.panel_label(ax, "b")
    rho = pd.read_csv(S.ROOT / f"{P.DIR}/{P.H3A}/sensitivity_screen_quality/scale_effect_rho.tsv", sep="\t")
    rho = rho.dropna(subset=["rho_primary"]).sort_values("rho_primary")
    y = np.arange(len(rho))
    ax.scatter(rho.rho_primary, y, s=8, color=S.OI["blue"], label="primary", zorder=3)
    ax.scatter(rho.rho_sens, y, s=8, facecolor="white", edgecolor=S.OI["orange"], lw=0.7,
               label="+ NNMD, CasActivity", zorder=3)
    for yy, a, b in zip(y, rho.rho_primary, rho.rho_sens):
        ax.plot([a, b], [yy, yy], color=S.OI["grey"], lw=0.5)
    ax.axvline(0, color="black", lw=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels([P.short_ctx(c) for c in rho.context], fontsize=4.6)
    ax.set_xlabel("Spearman ρ (unit essentiality,\nlineage β)")
    ax.set_title("Lineage scale effect", loc="left")
    ax.legend(frameon=False, fontsize=4.4, loc="upper left", handletextpad=0.2)

    # ---------- f: H2 rank against own splits
    ax = fig.add_subplot(gs[0, 2])
    S.panel_label(ax, "c")
    h2 = P.h2_ranks()
    nr = int(h2.n_splits.iloc[0]) + 1
    counts = h2.rank_vs_splits.value_counts().reindex(range(1, nr + 1), fill_value=0)
    ax.bar(counts.index, counts.values, color=S.OI["purple"], width=0.8)
    ax.axhline(len(h2) / nr, color="black", lw=0.6, ls="--")
    ax.text(nr, len(h2) / nr + 1.5, f"exchangeable: {len(h2) / nr:.1f} per rank", ha="right", fontsize=4.8)
    ax.set_xlabel(f"Rank of |t| among own {nr - 1} splits (1 = beats all)")
    ax.set_ylabel("H2 discoveries")
    ax.set_xlim(0.3, nr + 0.7)
    ax.set_title(f"H2 (exploratory), n = {len(h2)}", loc="left")
    ax.text(0.97, 0.80, f"rank 1: {int((h2.rank_vs_splits == 1).sum())}\nrank ≤ 2: "
            f"{int((h2.rank_vs_splits <= 2).sum())}", transform=ax.transAxes, ha="right", fontsize=5.2)

    # ---------- heatmap of complex-specific non-locus discoveries
    ax = fig.add_subplot(gs[1, 0:2])
    S.panel_label(ax, "d")
    both = pd.concat([da.assign(direction="H3a"), db.assign(direction="H3b")])
    foc = both[both.focus]
    k = 34
    units = (foc.assign(a=foc.t.abs()).groupby("unit").a.max().sort_values(ascending=False).head(k).index.tolist())
    sub = both[both.unit.isin(units)]
    ctx_order = ([c for c in ma.index if c in set(sub[sub.focus].context)]
                 + [c for c in mb.index if c in set(sub[sub.focus].context)])
    mat = np.full((len(units), len(ctx_order)), np.nan)
    other = np.zeros_like(mat, dtype=bool)
    for _, r in sub.iterrows():
        if r.context not in ctx_order:
            continue
        i, j = units.index(r.unit), ctx_order.index(r.context)
        if r.focus:
            mat[i, j] = r.beta
        else:
            other[i, j] = True
    lim = np.nanmax(np.abs(mat))
    im = ax.imshow(mat, cmap="RdBu", norm=TwoSlopeNorm(0, -lim, lim), aspect="auto", interpolation="none")
    oi, oj = np.where(other)
    ax.scatter(oj, oi, s=3, color=S.OI["grey"], marker="o", lw=0)
    n_a = sum(c in ma.index for c in ctx_order)
    ax.axvline(n_a - 0.5, color="black", lw=0.8)
    ax.set_xticks(range(len(ctx_order)))
    ax.set_xticklabels([P.short_ctx(c) for c in ctx_order], rotation=90, fontsize=4.8)
    ax.set_yticks(range(len(units)))
    ax.set_yticklabels([f"{u}  {P.gene_label(u, mem, 4)}" for u in units], fontsize=4.6)
    ax.set_xticks(np.arange(-0.5, len(ctx_order)), minor=True)
    ax.set_yticks(np.arange(-0.5, len(units)), minor=True)
    ax.grid(which="minor", color="#EEEEEE", lw=0.3)
    ax.tick_params(which="minor", length=0)
    ax.text(n_a / 2 - 0.5, -1.2, "H3a lineage", ha="center", fontsize=5.5)
    ax.text(n_a + (len(ctx_order) - n_a) / 2 - 0.5, -1.2, "H3b driver / state", ha="center", fontsize=5.5)
    ax.set_title(f"Complex-specific, non-locus discoveries ({len(units)} units with the largest |t|)", loc="left",
                 pad=12)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    cb.set_label("β (Chronos)", fontsize=5)
    cb.ax.tick_params(labelsize=5)

    # ---------- per-context discovery composition + never-expressed rate
    ax = fig.add_subplot(gs[1, 2])
    S.panel_label(ax, "e")
    rows = []
    for d, mref in ((da, ma), (db, mb)):
        for c in mref.index:
            g = d[d.context == c]
            if len(g) == 0:
                continue
            rows.append((c, int(g.focus.sum()), int((~g.focus & ~g.locus).sum()), int(g.locus.sum()),
                         nea.get(c, np.nan) if d is da else np.nan, d is da))
    comp = pd.DataFrame(rows, columns=["context", "focus", "other", "locus", "never", "h3a"])
    comp["total"] = comp.focus + comp.other + comp.locus
    comp = pd.concat([comp[comp.h3a].sort_values("total", ascending=False),
                      comp[~comp.h3a].sort_values("total", ascending=False)]).reset_index(drop=True)
    y = np.arange(len(comp))[::-1]
    ax.barh(y, comp.focus, color=C_FOCUS, height=0.75, label="complex-specific")
    ax.barh(y, comp.other, left=comp.focus, color=C_OTHER, height=0.75, label="other, non-locus")
    ax.barh(y, comp.locus, left=comp.focus + comp.other, color=C_LOCUS, height=0.75, label="locus")
    ax.set_yticks(y)
    ax.set_yticklabels([P.short_ctx(c) for c in comp.context], fontsize=4.8)
    ax.set_xlabel("Discoveries (A)")
    ax.set_xlim(0, 215)
    sep = int((~comp.h3a).sum()) - 0.5
    ax.axhline(sep, color="black", lw=0.6)
    ax.text(210, sep + 0.6, "H3a", fontsize=5, ha="right", va="bottom")
    ax.text(210, sep - 0.6, "H3b", fontsize=5, ha="right", va="top")
    ax2 = ax.twiny()
    ax2.scatter(comp.never, y, s=7, color=S.OI["vermillion"], marker="D", lw=0, zorder=4)
    ax2.axvline(0.05, color=S.OI["vermillion"], lw=0.5, ls=":")
    ax2.set_xlim(0, 0.4)
    ax2.set_xlabel("Never-expressed sets, fraction p < 0.05 (◆; H3a)", fontsize=5, color=S.OI["vermillion"])
    ax2.tick_params(axis="x", labelsize=5, colors=S.OI["vermillion"])
    ax2.spines["top"].set_visible(True)
    ax.legend(frameon=False, fontsize=4.8, loc="lower right", bbox_to_anchor=(1.0, 0.04))

    # ---------- d: sensitivity survival per context
    ax = fig.add_subplot(gs[2, :])
    S.panel_label(ax, "f")
    rows = []
    for d, mref, lab in ((da, ma, "H3a"), (db, mb, "H3b")):
        for c in mref.index:
            g = d[d.context == c]
            if len(g) < 3:
                continue
            f = g[g.focus]
            rows.append((c, lab, len(g), g.screen.mean(), g.cons.mean(),
                         f.screen.mean() if len(f) else np.nan, len(f)))
    sv = pd.DataFrame(rows, columns=["context", "dir", "n", "screen", "gpd", "screen_focus", "n_focus"])
    sv = pd.concat([sv[sv.dir == "H3a"].sort_values("n", ascending=False),
                    sv[sv.dir == "H3b"].sort_values("n", ascending=False)])
    x = np.arange(len(sv))
    ax.bar(x - 0.2, sv.screen, width=0.38, color=S.OI["orange"], label="screen-quality sensitivity (all discoveries)")
    ax.bar(x + 0.2, sv.gpd, width=0.38, color=S.OI["green"], label="GPD-tail sensitivity, p = max(p_perm, p)")
    ax.scatter(x - 0.2, sv.screen_focus, s=9, color="black", marker="_", lw=1.0, zorder=4,
               label="screen-quality, complex-specific non-locus subset")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{P.short_ctx(c)}\n({n})" for c, n in zip(sv.context, sv.n)], rotation=90, fontsize=4.8)
    nsep = int((sv.dir == "H3a").sum()) - 0.5
    ax.axvline(nsep, color="black", lw=0.6)
    ax.set_ylim(0, 1.3)
    ax.set_ylabel("Fraction of primary\ndiscoveries kept")
    ax.set_title("Sensitivity survival per context (contexts with ≥ 3 discoveries; n in brackets)", loc="left")
    ax.legend(frameon=False, ncol=3, fontsize=5, loc="upper left")

    S.save(fig, "fig3_context_atlas")


if __name__ == "__main__":
    main()
