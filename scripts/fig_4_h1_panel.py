"""Figure 4: H1 within-complex genetic interactions on construction set A (headline mechanistic panel).

Inputs (committed only): results/part2/directions/H1_within_complex_GI/{classes.tsv, discoveries.tsv,
sensitivity_SMARCB1_CN.tsv}, results/part2/directions/h1_label_final.json (E2 labels),
results/part2/fzx/h1_events.tsv (member composition), results/part2/controls/controls_summary.tsv (gene-level
control readouts quoted in panel c), results/part2/directions/H1S_silencing_GI/summary.json.

Panels
 a  complex-wide axis: kappa_N with 90% CI for every event with an evaluable N axis, TOST margin +/-0.1 (shaded),
    matched random-member null 2.5-97.5% (grey bar), reported label (colour), BH discovery (star); SMARCB1 with the
    pre-declared copy-number covariates as an extra row
 b  paralog axis: beta_P with 90% CI and its matched null, for the four events with an in-complex paralog
 c  member composition of the three headline-eligible events (ARID1A, EP300, FOCAD): event gene, paralogs (P),
    counted non-paralog members (N'), members without a baseline direction, most shifted member; the per-member
    coefficients (registered report item) are not in the committed run outputs and are not drawn
 d  reported labels of the 24 H1 events (complex-wide axis; paralog compensation counted on the paralog axis)
 e  reported labels of the 317 H1S events (exploratory)

Usage: python scripts/fig_4_h1_panel.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import fig_style as S
import part2_report as P

LAB_COL = {
    "PARALOG COMPENSATION": S.OI["vermillion"],
    "RESIDUAL SENSITISATION": S.OI["blue"],
    "COLLAPSE": S.OI["purple"],
    "DETECTABLE, WITHIN EQUIVALENCE MARGIN (negligible)": S.OI["orange"],
    "NEUTRAL": S.OI["green"],
    "INCONCLUSIVE": S.OI["grey"],
    "INCONCLUSIVE (null unavailable)": S.OI["grey"],
    "NOT EVALUABLE": "#BBBBBB",
    "not evaluable": "#BBBBBB",
}
SHORT_LAB = {
    "DETECTABLE, WITHIN EQUIVALENCE MARGIN (negligible)": "detectable, within margin (E2)",
    "INCONCLUSIVE (null unavailable)": "inconclusive (null unavailable)",
    "NOT EVALUABLE": "not evaluable (counted test)",
    "not evaluable": "not evaluable (|N'| < 2)",
}


def _syms(s) -> list[str]:
    return [] if not isinstance(s, str) or not s else [x.split(" ")[0] for x in s.split(";")]


def main() -> None:
    S.setup()
    c = P.h1_events()
    sens = pd.read_csv(S.ROOT / f"{P.DIR}/{P.H1}/sensitivity_SMARCB1_CN.tsv", sep="\t")
    ctl = pd.read_csv(S.ROOT / S.CONTROLS, sep="\t")
    lab = S.load_json(f"{P.DIR}/h1_label_final.json")

    fig = plt.figure(figsize=(7.2, 8.4))
    gs = fig.add_gridspec(3, 2, height_ratios=[2.1, 1.35, 0.85], width_ratios=[1.35, 1], hspace=0.48, wspace=0.5)

    # ---------- a: kappa_N forest
    ax = fig.add_subplot(gs[0, 0])
    S.panel_label(ax, "a")
    n = c[c.kappa_N.notna()].copy().sort_values("kappa_N")
    rows = []
    for _, r in n.iterrows():
        rows.append(dict(label=f"{r.gene}  ({r.unit}, {r.type}, n={r.n_A})", k=r.kappa_N, lo=r.kappa_N_ci90_lo,
                         hi=r.kappa_N_ci90_hi, nlo=r.N_null_q025, nhi=r.N_null_q975, lab=r.label_final,
                         disc=bool(r.disc_N)))
        if r.gene == "SMARCB1":
            cn = sens[sens.covariate_set == "CN-covariate sensitivity"].iloc[0]
            rows.append(dict(label=f"  SMARCB1 + CN covariates (n={int(cn.n_loss_complete)} of {int(cn.n_complete)})",
                             k=cn.kappa_N, lo=cn.kappa_N_ci90_lo, hi=cn.kappa_N_ci90_hi, nlo=cn.N_null_q025,
                             nhi=cn.N_null_q975, lab=cn.class_complex_wide, disc=False, sens=True))
    y = np.arange(len(rows))[::-1]
    ax.axvspan(-0.1, 0.1, color="#EEEEEE", zorder=0)
    ax.axvline(0, color="black", lw=0.5)
    for yy, r in zip(y, rows):
        if pd.notna(r["nlo"]):
            ax.plot([r["nlo"], r["nhi"]], [yy, yy], color="#9E9E9E", lw=3.2, alpha=0.55, solid_capstyle="butt",
                    zorder=1)
        col = LAB_COL.get(r["lab"], "black")
        ax.errorbar(r["k"], yy, xerr=[[r["k"] - r["lo"]], [r["hi"] - r["k"]]], fmt="o", color=col,
                    mfc="white" if r.get("sens") else col, ms=3.2, lw=0.8, capsize=1.2, zorder=3)
        if r["disc"]:
            ax.plot(r["hi"] + 0.015, yy, marker="*", color="black", ms=4, ls="")
    ax.set_yticks(y)
    ax.set_yticklabels([r["label"] for r in rows], fontsize=5.0)
    ax.set_xlabel("κ_N (Chronos; 90% CI)  < 0 collapse · > 0 residual sensitisation")
    ax.set_xlim(-0.18, 0.43)
    ax.set_title("Complex-wide axis (non-paralog members N')", loc="left")
    handles = []
    for k in ["RESIDUAL SENSITISATION", "DETECTABLE, WITHIN EQUIVALENCE MARGIN (negligible)", "NEUTRAL",
              "INCONCLUSIVE"]:
        handles.append(plt.Line2D([], [], marker="o", ls="", color=LAB_COL[k], ms=3.2, label=SHORT_LAB.get(k, k.lower())))
    handles.append(plt.Line2D([], [], color="#9E9E9E", lw=3.2, alpha=0.55, label="matched null 2.5–97.5%"))
    handles.append(plt.Line2D([], [], marker="*", ls="", color="black", ms=4, label="BH discovery (q < 0.05)"))
    handles.append(Rectangle((0, 0), 1, 1, color="#EEEEEE", label="TOST margin ±0.1"))
    ax.legend(handles=handles, frameon=False, fontsize=4.8, loc="upper right")

    # ---------- b: beta_P forest
    ax = fig.add_subplot(gs[0, 1])
    S.panel_label(ax, "b")
    p = c[c.beta_P.notna()].copy().sort_values("beta_P", ascending=False)
    y = np.arange(len(p))[::-1]
    ax.axvspan(-0.1, 0.1, color="#EEEEEE", zorder=0)
    ax.axvline(0, color="black", lw=0.5)
    for yy, (_, r) in zip(y, p.iterrows()):
        ax.plot([r.P_null_q025, r.P_null_q975], [yy, yy], color="#9E9E9E", lw=3.2, alpha=0.55,
                solid_capstyle="butt", zorder=1)
        col = LAB_COL.get(r.label_final_paralog, "black")
        ax.errorbar(r.beta_P, yy, xerr=[[r.beta_P - r.beta_P_ci90_lo], [r.beta_P_ci90_hi - r.beta_P]], fmt="o",
                    color=col, ms=3.4, lw=0.8, capsize=1.2, zorder=3)
        if r.disc_P:
            ax.plot(r.beta_P_ci90_lo - 0.06, yy, marker="*", color="black", ms=4, ls="")
        ax.text(0.35, yy, f"{SHORT_LAB.get(r.label_final_paralog, r.label_final_paralog.lower())}",
                fontsize=4.8, va="center", color=col)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.gene} → {', '.join(_syms(r.P))}\n({r.unit}, n={r.n_A})" for _, r in p.iterrows()],
                       fontsize=5.2)
    ax.set_xlim(-1.55, 0.95)
    ax.set_ylim(-0.7, len(p) - 0.3)
    ax.set_xlabel("β_P (Chronos; 90% CI)  < 0 = paralog more essential")
    ax.set_title("Paralog axis (in-complex paralogs P)", loc="left")

    # ---------- c: member composition of headline-eligible events
    ax = fig.add_subplot(gs[1, :])
    S.panel_label(ax, "c")
    ev = [e for e in lab[P.H1]["headline_eligible"]]
    role_col = {"event": "black", "P": S.OI["vermillion"], "N'": S.OI["blue"], "nobase": "#CCCCCC"}
    ctl_g = ctl[ctl.level == "G"]
    for i, e in enumerate(ev):
        r = c[c.event == e].iloc[0]
        yy = len(ev) - 1 - i
        nob = set(_syms(r.N_no_baseline))
        genes = ([(r.gene, "event")] + [(g, "P") for g in _syms(r.P)]
                 + [(g, "N'") for g in _syms(r.N) if g not in nob] + [(g, "nobase") for g in _syms(r.N) if g in nob])
        shifted = r.N_most_shifted.split(" ")[0] if isinstance(r.N_most_shifted, str) else None
        for j, (g, role) in enumerate(genes):
            ax.add_patch(Rectangle((j * 1.0, yy - 0.4), 0.92, 0.8, color=role_col[role], lw=0))
            ax.text(j + 0.46, yy, g, ha="center", va="center", fontsize=4.6,
                    color="white" if role in ("event", "P", "N'") else "black", rotation=90)
            if g == shifted:
                ax.plot(j + 0.46, yy + 0.5, marker="v", color="black", ms=3)
        if pd.notna(r.beta_P) and r.label_final_paralog == "PARALOG COMPENSATION":
            stat = (f"β_P {r.beta_P:+.2f} (95% CI {r.ci95_lo_P:+.2f} to {r.ci95_hi_P:+.2f}), null "
                    f"{r.P_null_q025:+.2f}…{r.P_null_q975:+.2f}; P−N contrast {r.contrast_beta:+.2f} "
                    f"(p {r.contrast_p:.1g}); {int(r.n_context_P)} loss lines")
        else:
            stat = (f"κ_N {r.kappa_N:+.3f} (95% CI {r.ci95_lo_N:+.3f} to {r.ci95_hi_N:+.3f}), null "
                    f"{r.N_null_q025:+.3f}…{r.N_null_q975:+.3f}; {r.N_frac_aligned_pos:.0%} of N' aligned > 0; "
                    f"{int(r.n_context_N)} loss lines")
        if pd.notna(r.kappa_N) and r.label_final_paralog == "PARALOG COMPENSATION":
            stat += f"\ncomplex-wide κ_N {r.kappa_N:+.3f} (90% CI {r.kappa_N_ci90_lo:+.3f} to {r.kappa_N_ci90_hi:+.3f}): " \
                    f"{r.label_final.lower()}"
        elif pd.isna(r.kappa_N):
            stat += f"\ncomplex-wide axis not evaluable (|N'| = {int(r.n_Nprime)})"
        ax.text(len(genes) + 0.4, yy, stat, va="center", fontsize=5.0)
        ax.text(-0.3, yy, f"{r.gene} loss\n{r.unit} ({r.type})", ha="right", va="center", fontsize=5.6,
                fontweight="bold")
    # gene-level control readouts on A, quoted with their own contexts
    notes = []
    for ent, gene in (("CTRL_PC6", "ARID1B"), ("CTRL_PC10", "PELO"), ("CTRL_PC10", "HBS1L")):
        g = ctl_g[(ctl_g.entry == ent) & ctl_g.unit.str.startswith(gene + " ")]
        if len(g):
            g = g.iloc[0]
            notes.append(f"{ent[5:]} {gene} in {g.context}: β {g.beta:+.2f} (one-sided p {g.p:.1g})")
    ax.text(-3.0, -0.95, "Gene-level control readouts on A (own contexts, Table S1): " + "; ".join(notes),
            fontsize=4.8, color="#555555", va="center")
    ax.set_xlim(-3.2, 30)
    ax.set_ylim(-1.2, len(ev) - 0.3)
    ax.axis("off")
    hs = [Rectangle((0, 0), 1, 1, color=role_col[k]) for k in role_col]
    ax.legend(hs, ["lost gene", "in-complex paralog (P)", "counted member (N')", "no baseline direction (|b| < 0.1)"],
              frameon=False, fontsize=4.8, ncol=4, loc="upper left", bbox_to_anchor=(0.0, 1.12))
    ax.set_title("Member composition of the headline-eligible events (▼ most shifted N member)", loc="left", pad=14)

    # ---------- d: label counts H1 vs H1S
    order = ["PARALOG COMPENSATION", "RESIDUAL SENSITISATION", "COLLAPSE",
             "DETECTABLE, WITHIN EQUIVALENCE MARGIN (negligible)", "NEUTRAL", "INCONCLUSIVE",
             "INCONCLUSIVE (null unavailable)", "NOT EVALUABLE", "not evaluable"]
    h1s = pd.read_csv(S.ROOT / f"{P.DIR}/{P.H1S}/classes.tsv", sep="\t")
    for k, (name, cnt) in enumerate((("H1 (24 events)", c.label_final.value_counts()),
                                     ("H1S (317 events, exploratory)", h1s.label_final.value_counts()))):
        ax = fig.add_subplot(gs[2, k])
        S.panel_label(ax, "d" if k == 0 else "e")
        vals = [int(cnt.get(o, 0)) for o in order]
        if k == 0:
            vals[0] = int((c.label_final_paralog == "PARALOG COMPENSATION").sum())
        yy = np.arange(len(order))[::-1]
        ax.barh(yy, vals, color=[LAB_COL[o] for o in order], height=0.7)
        for v, y0 in zip(vals, yy):
            ax.text(v + max(vals) * 0.02, y0, str(v), va="center", fontsize=5)
        ax.set_yticks(yy)
        ax.set_yticklabels([SHORT_LAB.get(o, o.lower()) + (" (paralog axis)" if (o == order[0]) else "")
                            for o in order], fontsize=5)
        ax.set_xlim(0, max(vals) * 1.15)
        ax.set_xlabel("Events")
        ax.set_title(name, loc="left")

    S.save(fig, "fig4_h1_panel")


if __name__ == "__main__":
    main()
