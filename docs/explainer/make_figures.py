"""Figures for docs/explainer/METHOD.md (matplotlib PNG, <= 150 dpi, colour-blind-safe palette).

  worked_figures(J, ex)  -> figures/worked_heatmap.png, worked_scores.png, worked_null.png (called by worked_example.py;
                            inputs: data/worked_example.json and data/worked_example_*.tsv / .npy, all construction set A)
  results_figure()       -> figures/results.png, from committed result files only (paths cited in the code)
figures/pipeline.svg is hand-written and not produced here.

Usage (project root):  .venv/bin/python docs/explainer/make_figures.py     # redraws every PNG from data/ + results/
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FIG = HERE / "figures"
DPI = 150
# reference categorical palette (validated light-mode slots 1-3) + neutrals
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED, GRID, LIGHT = "#14161a", "#5a5f69", "#d9dce2", "#a9adb5"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlelocation": "left", "figure.facecolor": "white",
    "axes.facecolor": "white", "savefig.facecolor": "white", "legend.frameon": False})


def _save(fig, name):
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / name, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------------------------------------ worked example
def worked_figures(J: dict, ex: pd.DataFrame):
    heatmap(J)
    scores_fig(J, ex)
    null_fig(J)


def heatmap(J):
    s2 = J["step2"]
    mem = s2["members"]
    M = np.array([[np.nan if v is None else v for v in r] for r in s2["member_score_matrix"]], float)
    near = s2["near_miss"][:3]
    # near-miss rows: their score to each member (from the JSON pair table is members only; recompute from A block
    # is done in worked_example.py and stored per gene as max/mean only), so show members x members plus a side bar
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(9.6, 4.0), gridspec_kw=dict(width_ratios=[1.15, 1], wspace=0.75))
    vmax = np.nanmax(M)
    im = ax.imshow(M, cmap="Blues", vmin=0, vmax=vmax)
    edges = {(p["a"], p["b"]): p for p in s2["pair_scores"]}
    for i in range(len(mem)):
        for j in range(len(mem)):
            if i == j:
                ax.text(j, i, "–", ha="center", va="center", color=MUTED)
                continue
            p = edges.get((mem[i], mem[j])) or edges.get((mem[j], mem[i]))
            col = "white" if M[i, j] > 0.55 * vmax else INK
            ax.text(j, i - 0.12, f"{M[i, j]:.1f}", ha="center", va="center", color=col, fontsize=9, fontweight="bold")
            ax.text(j, i + 0.22, f"{p['edge_in_subsample_graphs']}/20", ha="center", va="center", color=col, fontsize=7)
    ax.set_xticks(range(len(mem)), mem, rotation=40, ha="right")
    ax.set_yticks(range(len(mem)), mem)
    ax.spines[:].set_visible(False)
    ax.set_title("a  Member pair scores (zavg_BEF on A)")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("zavg_BEF (z units)")
    cb.ax.axhline(s2["z0"], color=ORANGE, lw=2)
    ax.annotate("small text: subsample graphs (of 20)\nin which the pair is an edge", xy=(0, -0.42),
                xycoords="axes fraction", fontsize=7.5, color=MUTED)
    # side panel: mean score to the five members, members vs near-miss genes
    rows = [(m_, float(np.nanmean(np.delete(M[i], i)))) for i, m_ in enumerate(mem)]
    rows += [(f"{n['gene']} ({n['cluster'] or 'no cluster'})", n["mean_score_to_members"]) for n in near]
    y = np.arange(len(rows))[::-1]
    cols = [BLUE] * len(mem) + [LIGHT] * len(near)
    bx.barh(y, [r[1] for r in rows], color=cols, height=0.6)
    bx.axvline(s2["z0"], color=ORANGE, lw=1.5, ls="--")
    bx.text(s2["z0"] + 0.2, -0.75, "z0 = 6 (edge threshold)", color=INK, fontsize=7.5, va="center")
    bx.set_yticks(y, [r[0] for r in rows], fontsize=8.5)
    bx.set_ylim(-1.1, len(rows) - 0.4)
    for yy, (_, v) in zip(y, rows):
        bx.text(v + 0.3, yy, f"{v:.1f}", va="center", fontsize=8, color=INK)
    bx.set_xlabel("mean zavg_BEF to the 5 SNAPC members")
    bx.set_title("b  Members vs nearest outsiders")
    bx.grid(axis="x", color=GRID, lw=0.6)
    bx.set_axisbelow(True)
    _save(fig, "worked_heatmap.png")


def scores_fig(J, ex):
    fig, ax = plt.subplots(figsize=(8.6, 3.9))
    rng = np.random.default_rng(0)
    out = ex[ex.bone == 0].score.dropna()
    inb = ex[ex.bone == 1].dropna(subset=["score"])
    groups = [("All other A lines", out.values, LIGHT)]
    for d, g in inb.groupby("disease"):
        groups.append((f"{d}", g.score.values, BLUE))
    groups = groups[:1] + sorted(groups[1:], key=lambda t: -len(t[1]))
    for i, (lab, v, col) in enumerate(groups):
        x = i + rng.uniform(-0.18, 0.18, len(v))
        ax.scatter(x, v, s=10 if i == 0 else 22, color=col, alpha=0.55 if i == 0 else 0.9, lw=0)
        ax.hlines(np.median(v), i - 0.3, i + 0.3, color=INK, lw=2)
        ax.text(i, ax.get_ylim()[0] if False else np.min(v) - 0.12, f"n={len(v)}", ha="center", va="top", fontsize=8,
                color=MUTED)
    s3 = J["step3"]
    ax.set_xticks(range(len(groups)), [g[0].replace(" ", "\n", 1) if len(g[0]) > 14 else g[0] for g in groups],
                  fontsize=8.5)
    ax.set_ylabel("SNAPC complex score\n(mean Chronos effect; lower = more dependent)")
    ax.axhline(np.median(out.values), color=LIGHT, lw=1, ls=":")
    ax.set_title(f"Complex score in Bone lines vs all other lines (A only; Bone median {s3['dist_in']['median']:.2f}, "
                 f"others {s3['dist_out']['median']:.2f})", fontsize=10)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_ylim(np.min(ex.score) - 0.35, np.max(ex.score) + 0.15)
    _save(fig, "worked_scores.png")


def null_fig(J):
    s4, s5 = J["step4"], J["step5"]
    t_obs = s4["refit"]["t"]
    tnull = np.load(HERE / "data/worked_example_perm_tnull.npy")
    allm = pd.read_csv(HERE / "data/worked_example_matched_t_bone_all_units.tsv", sep="\t").t.dropna().values
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(9.2, 3.6), gridspec_kw=dict(wspace=0.3))
    ax.hist(tnull, bins=60, color=LIGHT)
    ax.axvline(t_obs, color=BLUE, lw=2)
    ax.text(t_obs, ax.get_ylim()[1] * 0.9, f" observed t = {t_obs:.2f}", color=INK, fontsize=8.5)
    ax.set_xlabel("HC3 t under permutation (10,000 Freedman–Lane draws)")
    ax.set_ylabel("count")
    ax.set_title("a  Permutation null for SNAPC × Bone")
    ax.set_xlim(min(t_obs, tnull.min()) - 0.5, max(4.5, tnull.max() + 0.5))
    bins = np.linspace(-12, 12, 73)
    bx.hist(allm, bins=bins, color=LIGHT, density=True, label="all matched sets, Bone (3,500)")
    ms = [m["t_committed"] for m in s5["matched_sets"]]
    ne = [v for v in s5["never_expressed"]["t_values"] if v is not None]
    top = bx.get_ylim()[1]
    bx.vlines(ms, 0, top * 0.35, color=ORANGE, lw=1.6, label="SNAPC's 10 matched sets")
    bx.vlines(ne, top * 0.40, top * 0.58, color=AQUA, lw=1.0, label="never-expressed sets (100)")
    bx.axvline(t_obs, color=BLUE, lw=2, label="SNAPC observed")
    bx.set_xlabel("HC3 t of the Bone coefficient")
    bx.set_ylabel("density")
    bx.set_title("b  Observed t against reference gene sets")
    bx.legend(fontsize=7.5, loc="upper right")
    bx.set_xlim(-12, 10)
    _save(fig, "worked_null.png")


# ------------------------------------------------------------------------------------------------ results figure
def results_figure():
    tr = json.loads((ROOT / "results/test_eval/test_results.json").read_text())       # Part 1 single look (aggregate)
    st = tr["results"]["itt"]["bootstrap"]["cluster"]["stat"]
    rates = pd.read_csv(ROOT / "results/part2/lockbox_confirm/rates_primary_R7.tsv", sep="\t")
    r1a = pd.read_csv(ROOT / "results/part2/revision/r1a_rates.tsv", sep="\t")
    c3a = json.loads((ROOT / "results/part2/directions/H3a_lineage_atlas/calibration_D8.json").read_text())
    p3a = json.loads((ROOT / "results/part2/directions/H3a_lineage_atlas/calibration.json").read_text())
    br = json.loads((ROOT / "results/part2/construction/bridge/bridge.json").read_text())["primary"]
    bl = json.loads((ROOT / "results/part2/lockbox_confirm/bridge_primary_vs_corrected_R7.json").read_text())
    bl = bl["primary_evaluable"]["primary"]

    fig, axs = plt.subplots(2, 2, figsize=(10.5, 7.2), gridspec_kw=dict(hspace=0.62, wspace=0.62))
    # a: Part 1 test AP
    ax = axs[0, 0]
    meth = [("zavg_BEF (locked)", "zavg_BEF"), ("GLS", "gls_signed"), ("Spearman", "spearman_signed"),
            ("Pearson", "pearson_signed")]
    for i, (lab, k) in enumerate(meth):
        s = st[f"{k}|ap"]
        col = BLUE if i == 0 else MUTED
        ax.errorbar(s["point"], i, xerr=[[s["point"] - s["ci_low"]], [s["ci_high"] - s["point"]]], fmt="o",
                    color=col, ms=7, capsize=0, lw=2)
        ax.text(s["ci_high"] + 0.01, i, f"{s['point']:.3f}", va="center", fontsize=8.5, color=INK)
    ax.set_yticks(range(len(meth)), [m[0] for m in meth])
    ax.invert_yaxis()
    ax.axvline(0.0055, color=LIGHT, ls=":", lw=1)
    ax.set_xlim(0, 0.5)
    ax.set_xlabel("average precision on held-out hu.MAP 3.0 test pairs\n(95% cluster-bootstrap CI; prevalence 0.0055)")
    ax.set_title("a  Input score: one held-out test")
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)

    # b: lockbox confirmation rates
    ax = axs[0, 1]
    rows = [("Lineage atlas", "H3a_lineage_atlas"), ("Driver / genome state", "H3b_driver_family"),
            ("Sub-module (exploratory)", "H2_submodule_divergence"), ("Silencing (exploratory)", "H1S_silencing_GI"),
            ("Positive controls", "controls (23)")]
    for i, (lab, k) in enumerate(rows):
        r = rates[(rates.direction == k) & (rates.subset == "overall")].iloc[0]
        n_ev = int(r.primary_n_lockbox - r.primary_n_not_evaluable)
        ax.barh(i, r.primary_rate_evaluable, color=BLUE if i < 2 else LIGHT, height=0.6)
        ax.text(r.primary_rate_evaluable + 0.02, i, f"{int(r.primary_n_confirmed)}/{n_ev}", va="center", fontsize=8.5)
    ax.axvspan(0.47, 0.52, color=ORANGE, alpha=0.18, lw=0)
    ax.text(0.54, -0.85, "orange band: power-expected rate, lineage atlas (0.47–0.52)", ha="left", fontsize=7.5, color=INK)
    ax.set_yticks(range(len(rows)), [r[0] for r in rows])
    ax.invert_yaxis()
    ax.set_xlim(0, 1.12)
    ax.set_ylim(len(rows) - 0.5, -1.3)
    ax.set_xlabel("fraction of discoveries confirmed in lockbox L\n(evaluable fits; same sign, one-sided p < 0.05)")
    ax.set_title("b  Discoveries that replicate in held-out lines")
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)

    # c: matched-set parity (H3a, fraction of tests with p < 0.05)
    ax = axs[1, 0]
    g = lambda scope, typ: float(r1a[(r1a.direction == "H3a_lineage_atlas") & (r1a.scope == scope) &  # noqa: E731
                                     (r1a.set_type == typ) & (r1a.context == "pooled")].frac_p_lt_005.iloc[0])
    bars = [("Real complexes (367)", g("all_units", "real"), BLUE),
            ("Sets from other complexes,\nmatched size+ess.+coherence", g("all_units", "member_sets"), MUTED),
            ("Non-member sets, matched\nsize+essentiality only*", g("all_units", "nonmember_ess"), MUTED),
            ("Never-expressed gene sets", c3a["never_expressed_sets"]["pooled"]["frac_p_lt_005"], MUTED),
            ("Pseudo-contexts (random lines)", p3a["frac_p_lt_005"], MUTED)]
    for i, (lab, v, col) in enumerate(bars):
        ax.barh(i, v, color=col if i == 0 else LIGHT, height=0.6)
        ax.text(v + 0.006, i, f"{v:.3f}", va="center", fontsize=8.5)
    ax.axvline(0.05, color=ORANGE, ls="--", lw=1.4)
    ax.text(0.054, -0.62, "5% expected under no effect", fontsize=7.5, color=INK, va="center")
    ax.set_yticks(range(len(bars)), [b[0] for b in bars], fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlim(0, 0.36)
    ax.set_ylim(len(bars) - 0.5, -0.9)
    ax.set_xlabel("fraction of lineage tests with p < 0.05 (construction set A)\n* post hoc revision analysis")
    ax.set_title("c  Complexes vs reference gene sets")
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)

    # d: bridges and BENCH (differences with 95% CI)
    ax = axs[1, 1]
    items = [("BioPlex bridge (primary)\nzavg − GLS map precision", br["diff"], br["ci95"]),
             ("BRIDGE-L (secondary)\nzavg − GLS confirmation rate", bl["diff"], bl["ci95"])]
    for i, (lab, d, ci) in enumerate(items):
        ax.errorbar(d, i, xerr=[[d - ci[0]], [ci[1] - d]], fmt="o", color=BLUE, ms=7, lw=2, capsize=0)
        ax.text(d, i - 0.28, f"{d:+.3f} [{ci[0]:+.3f}, {ci[1]:+.3f}]", ha="center", fontsize=8, color=INK)
    ax.axvline(0, color=MUTED, lw=1)
    ax.set_xticks([-0.10, -0.05, 0.0, 0.025], ["−0.10", "−0.05", "0", "+0.025"])
    ax.set_yticks(range(len(items)), [it[0] for it in items], fontsize=8.5)
    ax.invert_yaxis()
    ax.set_ylim(len(items) - 0.5, -0.8)
    ax.set_xlabel("zavg map − GLS map (95% CI); below 0 favours GLS")
    ax.set_title("d  Does the zavg map beat a GLS map?")
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    _save(fig, "results.png")


def main():
    results_figure()
    J = json.loads((HERE / "data/worked_example.json").read_text())
    ex = pd.read_csv(HERE / "data/worked_example_scores_A.tsv", sep="\t", index_col=0)
    worked_figures(J, ex)


if __name__ == "__main__":
    main()
