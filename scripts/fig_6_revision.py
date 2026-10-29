"""Figure 6 (revision round 1): post hoc, descriptive reference analyses R-1, R-3 and R-4 (A only).

Every panel is POST HOC and DESCRIPTIVE (results/part2/revision/README.md; intent 4157f1e). Nothing here is a
registered test, and no lockbox value is shown except the committed R7 primary BRIDGE-L rates in panel d.

Inputs (committed only; no dataset, data_cache or persisted lockbox read is opened):
  results/part2/revision/r1a_rates.tsv            R-1(a) fraction p < 0.05: real units, complex-member matched sets,
                                                  non-member coherent sets (N-nbr), essentiality-only sets (N-ess)
  results/part2/revision/r1b_per_context.tsv      R-1(b) 'beats all matched sets' among real discoveries against the
                                                  same flag among random-atlas discoveries (selection-aware reference)
  results/part2/revision/r1_summary.json          R-1(b) pooled rates
  results/part2/revision/r3_demeter2_per_hit.tsv  R-3 DEMETER2 beta on A models for the 207 H3a strongest findings
  results/part2/revision/r3_summary.json          R-3 Spearman correlations
  results/part2/revision/r4_by_size_bin.tsv       R-4 BRIDGE-L observed and power-expected rate per size bin and map
  results/part2/revision/r4_summary.json          R-4 headline (EB mean + variance)

Panels
 a  R-1(a): fraction of tests with p < 0.05 for real complexes and three kinds of matched gene sets, each on the units
    where that kind of set could be matched (common-unit comparison)
 b  R-1(b): per context, rate at which real discoveries beat all their matched sets (y) against the rate at which
    random-atlas discoveries beat all others in their group (x); identity line; pooled values as large symbols
 c  R-3: A beta (Chronos) against DEMETER2 beta on A models for the H3a strongest findings; revised leads highlighted,
    discordant leads open
 d  R-4: BRIDGE-L confirmation rate per size bin, observed (R7 primary) and expected from L's power (EB mean +
    variance), for zavg-map and GLS-map discoveries

Usage: python scripts/fig_6_revision.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import fig_style as S

REV = "results/part2/revision"
DIRS = (("H3a_lineage_atlas", "H3a lineage"), ("H3b_driver_family", "H3b driver/state"))


def _tsv(name: str) -> pd.DataFrame:
    return pd.read_csv(S.ROOT / REV / name, sep="\t")


def panel_a(ax) -> None:
    r = _tsv("r1a_rates.tsv")
    r = r[r.context.eq("pooled")].set_index(["direction", "scope", "set_type"]).frac_p_lt_005
    groups = [
        ("units with\ncoherent sets", "common_units_nbr",
         [("real", "Real complexes", S.OI["blue"]), ("member_sets", "Sets from other complexes", S.OI["skyblue"]),
          ("nonmember_nbr", "Coherent sets, no complex", S.OI["green"])]),
        ("units with\ness.-only sets", "common_units_ess",
         [("real", "Real complexes", S.OI["blue"]), ("member_sets", "Sets from other complexes", S.OI["skyblue"]),
          ("nonmember_ess", "Essentiality-only sets", S.OI["orange"])]),
    ]
    w = 0.13
    xt, xl = [], []
    x0 = 0.0
    seen = set()
    for did, dlab in DIRS:
        for glab, scope, bars in groups:
            for i, (st, lab, col) in enumerate(bars):
                v = r.get((did, scope, st), np.nan)
                ax.bar(x0 + (i - 1) * w, v, w * 0.92, color=col, label=None if lab in seen else lab)
                seen.add(lab)
            xt.append(x0)
            xl.append(f"{dlab.split()[0]}\n{glab}")
            x0 += 0.55
        x0 += 0.15
    ax.axhline(0.05, color="#555555", lw=0.6, ls=":")
    ax.text(x0 - 0.25, 0.052, "0.05", fontsize=5.5, color="#555555", va="bottom", ha="right")
    ax.set_xticks(xt)
    ax.set_xticklabels(xl, fontsize=5.3)
    ax.set_ylabel("Fraction of tests with p < 0.05")
    ax.set_title("Context association: complexes vs matched sets", fontsize=7.5)
    ax.legend(frameon=False, fontsize=5.3, loc="upper right")
    ax.set_ylim(0, 0.36)


def panel_b(ax) -> None:
    p = _tsv("r1b_per_context.tsv")
    summ = S.load_json(f"{REV}/r1_summary.json")["directions"]
    cols = {"H3a_lineage_atlas": S.OI["blue"], "H3b_driver_family": S.OI["vermillion"]}
    mk = {"H3a_lineage_atlas": "o", "H3b_driver_family": "s"}
    for did, dlab in DIRS:
        q = p[p.direction.eq(did) & p.n_real_disc.ge(10)]
        ax.scatter(q.rand_symmetric_flag_rate, q.real_beats_all_rate, s=np.clip(q.n_real_disc, 10, 200) / 4,
                   facecolor="none", edgecolor=cols[did], lw=0.6, marker=mk[did], alpha=0.8)
        sm = summ[did]["b"]["summary"]
        ax.scatter([sm["random_symmetric_flag_rate"]], [sm["real_beats_all_rate"]], s=40, color=cols[did],
                   marker=mk[did], zorder=5, label=f"{dlab}, pooled")
    ax.axvline(1 / 11, color="#999999", lw=0.6, ls="--")
    ax.text(1 / 11 + 0.01, 0.97, "1/11 (no selection)", fontsize=5.3, color="#666666", va="top")
    ax.plot([0, 1], [0, 1], color="#555555", lw=0.6, ls=":")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Random-atlas discoveries that beat all\nothers in their group (selection-aware reference)")
    ax.set_ylabel("Real discoveries that beat\nall 10 matched sets")
    ax.set_title("'Beats all matched sets', conditional on discovery", fontsize=7.5)
    ax.legend(frameon=False, fontsize=5.3, loc="lower right")


def panel_c(ax) -> None:
    d = _tsv("r3_demeter2_per_hit.tsv").dropna(subset=["beta_D2"])
    s3 = S.load_json(f"{REV}/r3_summary.json")
    other = d[~d.lead_corrected.astype(bool)]
    lead = d[d.lead_corrected.astype(bool)]
    ax.scatter(other.beta_A, other.beta_D2, s=6, color=S.OI["grey"], lw=0, alpha=0.7,
               label="Other H3a strongest findings")
    conc = lead[lead.same_sign.astype(bool)]
    disc = lead[~lead.same_sign.astype(bool)]
    ax.scatter(conc.beta_A, conc.beta_D2, s=14, color=S.OI["blue"], lw=0, label="Revised leads, same sign")
    ax.scatter(disc.beta_A, disc.beta_D2, s=14, facecolor="none", edgecolor=S.OI["vermillion"], lw=0.8,
               label="Revised leads, opposite sign")
    ax.axhline(0, color="#555555", lw=0.5)
    ax.axvline(0, color="#555555", lw=0.5)
    ax.set_xlabel("β in A, CRISPR (Chronos)")
    ax.set_ylabel("β in A, DEMETER2 (RNAi)")
    lr, ar = s3["leads"], s3["all_207"]
    ax.set_title("Orthogonal RNAi check of H3a strongest findings", fontsize=7.5)
    ax.text(0.02, 0.98, f"Leads: {lr['same_sign']}/{lr['n_evaluable']} same sign, "
                        f"Spearman {lr['spearman_beta_A_vs_D2']:.2f}\n"
                        f"All: {ar['same_sign']}/{ar['n_evaluable']} same sign, "
                        f"Spearman {ar['spearman_beta_A_vs_D2']:.2f}",
            transform=ax.transAxes, fontsize=5.3, va="top", bbox=dict(fc="white", ec="none", pad=0.5))
    ax.legend(frameon=False, fontsize=5.3, loc="lower right")


def panel_d(ax) -> None:
    b = _tsv("r4_by_size_bin.tsv")
    h = S.load_json(f"{REV}/r4_summary.json")["headline"]
    bins = ["3-4", "5-8", "9-16", "17-60"]
    x = np.arange(len(bins))
    for off, mp, col, lab in ((-0.12, "zavg", S.COLOR["zavg"], "zavg map"), (0.12, "gls", S.COLOR["gls"], "GLS map")):
        q = b[b["map"].eq(mp)].set_index("bin").loc[bins]
        ax.plot(x + off, q.observed, "o", color=col, ms=4, label=f"{lab}: observed")
        ax.plot(x + off, q.expected, "_", color=col, ms=9, mew=1.4, label=f"{lab}: expected from power")
        for xx, o, e in zip(x + off, q.observed, q.expected):
            ax.plot([xx, xx], [e, o], color=col, lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(bins)
    ax.set_xlabel("Cluster size bin (genes)")
    ax.set_ylabel("BRIDGE-L confirmation rate")
    ax.set_ylim(0.22, 0.66)
    ax.set_title("BRIDGE-L: observed vs power-expected", fontsize=7.5)
    ax.text(0.02, 0.98, f"Size-standardised Δ (zavg − GLS):\nobserved {h['obs_diff']:+.3f}, "
                        f"expected {h['exp_diff']:+.3f} (SD {h['sd_diff']:.3f})",
            transform=ax.transAxes, fontsize=5.3, va="top")
    ax.legend(frameon=False, fontsize=5.0, loc="lower right", ncol=2)


def main() -> None:
    S.setup()
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 6.2))
    for f, ax, lab in ((panel_a, axs[0, 0], "a"), (panel_b, axs[0, 1], "b"), (panel_c, axs[1, 0], "c"),
                       (panel_d, axs[1, 1], "d")):
        f(ax)
        S.panel_label(ax, lab)
    fig.tight_layout(w_pad=2.0, h_pad=2.2)
    for p in S.save(fig, "fig6_revision_references"):
        print(p.relative_to(S.ROOT))


if __name__ == "__main__":
    main()
