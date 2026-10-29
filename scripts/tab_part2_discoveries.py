"""Part 2 discovery tables (A side), built only from committed result files.

 part2_derived.json   derived counts cited in the manuscript (re-tabulations of committed outputs, not new tests)
 Table S3  H3a / H3b discoveries with the D8 'beats all matched sets' flag, the post hoc 95% filter, locus flag,
           GPD and screen-quality sensitivity survival; Markdown shows the top complex-specific non-locus hits
 Table S4  H1 panel: 24 genetic-loss events, engine class and reported label (E2), kappa_N / beta_P with CIs,
           matched-null bands, member composition
 Table S5  H2 (exploratory): each discovery's rank against its own 20 random member splits

Usage: python scripts/tab_part2_discoveries.py
"""
from __future__ import annotations

import pandas as pd

import fig_style as S
import part2_report as P
from tab_manuscript import write


def h3_tables(derived: dict) -> None:
    mem = P.members()
    allrows = []
    for ID, short in ((P.H3A, "H3a"), (P.H3B, "H3b")):
        d = P.discoveries(ID)
        ne = P.never_expressed(ID)
        mr = P.matched_rates(ID)
        per = {}
        for ctx, g in d.groupby("context"):
            per[ctx] = {
                "discoveries": len(g), "beats_all": int(g.beats_all.sum()), "with_sets": int(g.has_sets.sum()),
                "focus": int(g.focus.sum()), "locus": int(g.locus.sum()), "screen_kept": int(g.screen.sum()),
                "gpd_kept": int(g.cons.sum()), "outside_beta": int(g.outside_beta.sum()),
            }
        for ctx in mr.index:
            per.setdefault(ctx, {"discoveries": 0})
            per[ctx]["real_frac_p_lt_005"] = float(mr.loc[ctx, "real"])
            per[ctx]["matched_frac_p_lt_005"] = float(mr.loc[ctx, "random"])
            per[ctx]["never_expressed_frac_p_lt_005"] = float(ne.get(ctx, float("nan")))
        f = d[d.focus]
        derived[short] = {
            "n_discovery": len(d), "n_locus": int(d.locus.sum()), "n_nonlocus": int((~d.locus).sum()),
            "locus_frac": float(d.locus.mean()), "n_beats_all": int(d.beats_all.sum()),
            "n_without_sets": int((~d.has_sets).sum()), "n_outside_beta": int(d.outside_beta.sum()),
            "n_focus_beats_all_nonlocus": len(f), "n_focus_units": int(f.unit.nunique()),
            "n_units_any_discovery": int(d.unit.nunique()),
            "n_focus_screen_kept": int(f.screen.sum()), "n_focus_gpd_kept": int(f.cons.sum()),
            "n_focus_screen_and_gpd_kept": int((f.screen & f.cons).sum()),
            "screen_kept_beats_all": [int((d.beats_all & d.screen).sum()), int(d.beats_all.sum())],
            "screen_kept_not_beats_all": [int((~d.beats_all & d.screen).sum()), int((~d.beats_all).sum())],
            "n_abs_beta_lt_0.05": int((d.beta.abs() < 0.05).sum()),
            "n_focus_abs_beta_lt_0.05": int((f.beta.abs() < 0.05).sum()),
            "n_families_with_discoveries": int(d.context.nunique()),
            "per_context": per,
        }
        d = d.assign(direction=short, members=[P.gene_label(u, mem, 8) for u in d.unit])
        allrows.append(d)
    t = pd.concat(allrows)
    t = t.assign(abs_t=t.t.abs()).sort_values(["direction", "abs_t"], ascending=[True, False])
    cols = ["direction", "context", "unit", "members", "n", "n_context", "beta", "ci95_lo", "ci95_hi", "t", "p",
            "p_perm_report", "q_final", "tier", "locus", "beats_all", "has_sets", "outside_beta", "cons", "screen",
            "control"]
    t[cols].to_csv(S.TAB_DIR / "tableS3_context_discoveries_full.tsv", sep="\t", index=False)
    top = []
    for short, k in (("H3a", 25), ("H3b", 15)):
        x = t[(t.direction == short) & t.focus].head(k).copy()
        x["β [95% CI]"] = [S.fmt_ci(b, lo, hi, 3, True) for b, lo, hi in zip(x.beta, x.ci95_lo, x.ci95_hi)]
        x["abs t"] = x.abs_t.map(lambda v: f"{v:.1f}")
        x["context"] = x.context.map(P.short_ctx)
        x["kept: GPD / screen"] = [f"{'yes' if a else 'no'} / {'yes' if b else 'no'}" for a, b in zip(x.cons, x.screen)]
        top.append(x[["direction", "context", "unit", "members", "n_context", "β [95% CI]", "abs t", "tier",
                      "kept: GPD / screen", "control"]].rename(columns={"n_context": "context lines (A)"}))
    md = pd.concat(top)
    write("tableS3_context_discoveries", "Supplementary Table S3: complex-specific context discoveries (A side)",
          ["Rows: the strongest H3a / H3b discoveries (ranked by |t|, never by the GPD p) that are not in a locus "
           "cluster and whose |t| exceeds that of every one of the unit's 10 matched random gene sets in the same "
           "context (D8 'beats all matched sets' flag). tableS3_context_discoveries_full.tsv lists all 1,306 discoveries with every flag.",
           "Members: first 8 member symbols of the frozen unit (gene content only; no curated or literature "
           "label). β in Chronos units: negative = more dependent in the context.",
           "'control readout': the same context and unit as a registered positive control (Table S1); 'control unit, other context': the unit is a control readout in a different (sometimes overlapping) context. Such rows are "
           "not independent discoveries.",
           "Lockbox confirmation: [[LOCKBOX: per-discovery L effect, 95% CI and confirmation flag]].",
           "Sources: results/part2/directions/{H3a_lineage_atlas,H3b_driver_family}/{discoveries,"
           "discoveries_report}.tsv, calibration_D8_tests.tsv.gz, sensitivity_screen_quality/discoveries.tsv; "
           "specificity_discoveries.tsv.gz; script scripts/tab_part2_discoveries.py."], md)


def h1_table(derived: dict) -> None:
    c = P.h1_events()
    lab = S.load_json(f"{P.DIR}/h1_label_final.json")
    derived["H1"] = {
        "label_counts": lab[P.H1]["label_counts"], "headline_eligible": lab[P.H1]["headline_eligible"],
        "n_events": len(c), "n_neutral": int((c.label_final == "NEUTRAL").sum()),
        "neutral_events": sorted(c.loc[c.label_final == "NEUTRAL", "gene"]),
        "inconclusive_events": sorted(c.loc[c.label_final.str.startswith("INCONCLUSIVE"), "gene"]),
    }
    derived["H1S"] = {"label_counts": lab[P.H1S]["label_counts"]}

    def f(v, nd=3):
        return "" if pd.isna(v) else f"{v:+.{nd}f}"

    rows = []
    for _, r in c.iterrows():
        rows.append({
            "event": r.event.replace("|", " · "), "type": r.type, "link": r.link, "loss lines (A)": r.n_A,
            "N' / P": f"{r.n_Nprime} / {r.n_P}",
            "κ_N [90% CI]": "" if pd.isna(r.kappa_N) else f"{f(r.kappa_N)} [{f(r.kappa_N_ci90_lo)}, {f(r.kappa_N_ci90_hi)}]",
            "κ_N null 2.5–97.5%": "" if pd.isna(r.N_null_q025) else f"{f(r.N_null_q025)} to {f(r.N_null_q975)}",
            "β_P [90% CI]": "" if pd.isna(r.beta_P) else f"{f(r.beta_P)} [{f(r.beta_P_ci90_lo)}, {f(r.beta_P_ci90_hi)}]",
            "β_P null 2.5–97.5%": "" if pd.isna(r.P_null_q025) else f"{f(r.P_null_q025)} to {f(r.P_null_q975)}",
            "q (N / P)": " / ".join("" if pd.isna(q) else f"{q:.2g}" for q in (r.q_N, r.q_P)),
            "engine class (complex-wide / paralog)": f"{r.class_complex_wide} / {r.class_paralog}",
            "reported label (complex-wide / paralog)": f"{r.label_final} / {r.label_final_paralog}",
            "most shifted N member": r.N_most_shifted if isinstance(r.N_most_shifted, str) else "",
            "P (paralogs)": r.P if isinstance(r.P, str) else "",
        })
    df = pd.DataFrame(rows)
    write("tableS4_h1_panel", "Supplementary Table S4: H1 within-complex genetic interactions (A side)",
          ["κ_N: loss coefficient on the aligned score of non-paralog members with a baseline direction (N'); "
           "> 0 = their dependency grows in its baseline direction (residual sensitisation), < 0 = it shrinks "
           "toward 0 (collapse). β_P: loss coefficient on the mean in-complex paralog effect; < 0 = paralogs more "
           "essential. Chronos units. 90% CIs (the TOST interval; margin ±0.1).",
           "q: BH over the 24 evaluable counted tests (registered m 28; PTEN's q is 0.051 at m = 28).",
           "Reported labels apply E2 (post hoc, conservative): a non-neutral class whose 90% CI also lies inside "
           "±0.1 is 'DETECTABLE, WITHIN EQUIVALENCE MARGIN (negligible)'. MAP2K4 (11 complete-case loss lines) "
           "is NOT EVALUABLE (R2).",
           "Per-member coefficients (β_j with 95% CI) are a registered report item but are not in the committed "
           "run outputs; see Fig. 4 legend.",
           "Lockbox: [[LOCKBOX: per-event L confirmation of κ_N / β_P]].",
           "Sources: results/part2/directions/H1_within_complex_GI/{classes,discoveries}.tsv, "
           "results/part2/fzx/h1_events.tsv; script scripts/tab_part2_discoveries.py."], df)
    c.to_csv(S.TAB_DIR / "tableS4_h1_panel_full.tsv", sep="\t", index=False)


def h2_table(derived: dict) -> None:
    mem = P.members()
    d = P.h2_ranks()
    vc = d.rank_vs_splits.value_counts().sort_index()
    derived["H2"] = {
        "n_discovery": len(d), "n_submodules": int(d.unit.nunique()), "n_parents": int(d.parent.nunique()),
        "n_locus": int(d.locus.sum()), "rank_counts": {int(k): int(v) for k, v in vc.items()},
        "n_rank1": int((d.rank_vs_splits == 1).sum()), "n_rank_le2": int((d.rank_vs_splits <= 2).sum()),
        "n_rank_le5": int((d.rank_vs_splits <= 5).sum()),
        "expected_rank1_if_exchangeable": len(d) / 21, "median_rank": float(d.rank_vs_splits.median()),
        "n_splits_per_discovery": sorted(set(d.n_splits)),
    }
    d = d.assign(abs_t=d.t.abs()).sort_values(["rank_vs_splits", "abs_t"], ascending=[True, False])
    d["members"] = [P.gene_label(u, mem, 8) for u in d.unit]
    d.drop(columns=["abs_t"]).to_csv(S.TAB_DIR / "tableS5_h2_split_ranks_full.tsv", sep="\t", index=False)
    x = d.copy()
    x["context"] = x.context.map(P.short_ctx)
    x["β [95% CI]"] = [S.fmt_ci(b, lo, hi, 3, True) for b, lo, hi in zip(x.beta, x.ci95_lo, x.ci95_hi)]
    x["t"] = x.t.map(lambda v: f"{v:+.1f}")
    x["rank vs 20 splits"] = x.rank_vs_splits
    md = x[["context", "unit", "members", "β [95% CI]", "t", "rank vs 20 splits", "locus", "tier"]]
    write("tableS5_h2_split_ranks", "Supplementary Table S5: H2 sub-module divergence (exploratory), rank against "
          "own random splits",
          ["All 109 stage-wise discoveries. Rank = 1 + number of the discovery's own 20 random member splits "
           "(same parent, same size, residualised the same way, seed 4800) whose |t| is at least the real |t|; "
           "rank 1 = beats all 20. Under exchangeability rank 1 has probability 1/21.",
           "β: coefficient of the context on the sub-module score residualised on its sibling complement "
           "(Chronos units). H2 is exploratory (no evaluable control).",
           "Sources: results/part2/directions/H2_submodule_divergence/{discoveries.tsv,"
           "calibration_D8_tests.tsv.gz}; script scripts/tab_part2_discoveries.py."], md)


def main() -> None:
    S.setup()
    derived: dict = {"label": "re-tabulation of committed Part 2 A-side outputs by scripts/tab_part2_discoveries.py; "
                              "not a trial, not a test"}
    h3_tables(derived)
    h1_table(derived)
    h2_table(derived)
    P.jdump(derived, "manuscript/tables/part2_derived.json")


if __name__ == "__main__":
    main()
