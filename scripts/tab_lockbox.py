"""Lockbox tables (Supplementary Tables S6-S8) and derived counts, built only from committed result files.

Inputs (all committed; no dataset, data_cache, persisted lockbox read or hu.MAP file is opened):
  results/part2/lockbox_confirm/comparison_primary_vs_corrected_R7.tsv, rates.tsv   rates per direction (R7)
  results/part2/lockbox_confirm/primary_relabelled_R7.tsv.gz                         per-hit primary / corrected calls
  results/part2/lockbox_confirm/corrected_refit_nonconfirmatory.tsv                  relook statistics (560 fits)
  results/part2/lockbox_confirm/<ID>/confirmation.tsv                                L statistics per hit
  results/part2/lockbox_confirm/h1_events_primary_vs_corrected_R7.tsv, H1_within_complex_GI_events.tsv
  results/part2/lockbox_confirm/BENCH_complex_vs_gene/{bench.json, bench_pairs.tsv}
  results/part2/lockbox_confirm/BRIDGE-L_confirmation_rate/bridge_hits.tsv, bridge_primary_vs_corrected_R7.json
  results/part2/verify/lockbox_results_review/{tables_review.json, bench_recompute.json}   power aid (post hoc), BENCH
  plus the A-side discovery files read through scripts/part2_report.py.

Outputs (manuscript/tables/):
  lockbox_derived.json            counts cited in the manuscript (re-tabulation; not a trial, not a test)
  tableS6_lockbox_rates.{md,tsv}  confirmation rate per direction and subset: primary evaluable, primary all-in,
                                  corrected (non-confirmatory), defect-sensitive count, post hoc power expectation
  tableS7_lockbox_h3.{md,tsv}, tableS7_lockbox_h3_full.tsv
                                  every H3a / H3b discovery with its L effect and calls; the Markdown lists the
                                  registered strongest findings (results/part2/robustness/strongest_findings.tsv,
                                  = Table S9) in registered rank order. Revision round 1 (2026-10-29): the Markdown
                                  previously listed an earlier, unregistered 'strongest specific' subset, kept in the
                                  full TSV as the column earlier_unregistered_set.
  tableS8_lockbox_h1_h2_bench_bridge.{md,tsv}
                                  H1 per event axis (primary and corrected), H2 by split rank, BENCH, BRIDGE-L

Definitions used here (fixed by D13/D14, not chosen after seeing these tables):
  defect-sensitive = the primary (as executed) and corrected relook calls differ (53 calls); no claim.
  beats all matched sets = D8 reporting flag (|t| above every matched random set of the unit in that context).

Usage: python scripts/tab_lockbox.py
"""
from __future__ import annotations

import pandas as pd

import fig_style as S
import part2_report as P
from tab_manuscript import write

LB = "results/part2/lockbox_confirm"
REV = "results/part2/verify/lockbox_results_review"


def _tsv(rel: str) -> pd.DataFrame:
    return pd.read_csv(S.ROOT / rel, sep="\t")


def relabelled() -> pd.DataFrame:
    r = _tsv(f"{LB}/primary_relabelled_R7.tsv.gz")
    r["defect_sensitive"] = r.defect_sensitive.astype(str).eq("True")
    return r


def h3_lockbox(ID: str, rel: pd.DataFrame) -> pd.DataFrame:
    """A-side discoveries (part2_report.discoveries) joined with their L statistics and R7 calls."""
    d = P.discoveries(ID)
    c = _tsv(f"{LB}/{ID}/confirmation.tsv")[
        ["context", "unit", "route", "n_context_L", "beta_L", "ci95_L_lo", "ci95_L_hi", "p_L"]]
    r = rel[rel.entry == ID][["context", "unit", "primary_call", "corrected_call", "defect_sensitive"]]
    x = d.merge(c, on=["context", "unit"], how="left").merge(r, on=["context", "unit"], how="left")
    rf = _tsv(f"{LB}/corrected_refit_nonconfirmatory.tsv")
    rf = rf[rf.entry == ID][["context", "unit", "beta_L_relook", "p_L_relook", "ci95_L_lo_relook",
                             "ci95_L_hi_relook"]]
    x = x.merge(rf, on=["context", "unit"], how="left")
    assert len(x) == len(d)
    x["confirmed_primary"] = x.primary_call.eq("confirmed")
    x["rare_pass"] = x.primary_call.eq("not lockbox-testable: rare-context route passed")
    x["strongest"] = (x.confirmed_primary & x.beats_all & ~x.locus & ~x.defect_sensitive
                      & x.control.ne("control readout"))
    x["scale_caveat"] = x.context.eq("lineage:Lymphoid") & (x.beta > 0)
    return x


def rates_table(derived: dict) -> pd.DataFrame:
    comp = _tsv(f"{LB}/comparison_primary_vs_corrected_R7.tsv")
    rv = S.load_json(f"{REV}/tables_review.json")
    pw, pwb = rv["power"], rv["power_by_beats_all"]

    def band(ID):
        v = [pw[ID][k]["expected_rate"] for k in pw[ID] if isinstance(pw[ID][k], dict) and "EB" in k]
        return min(v), max(v)

    exp = {("H3a_lineage_atlas", "overall"): band("H3a_lineage_atlas"),
           ("H3b_driver_family", "overall"): band("H3b_driver_family")}
    for ID in ("H3a_lineage_atlas", "H3b_driver_family"):
        for flag, sub in (("1.0", "beats_all_matched_sets"), ("0.0", "not_beats_all")):
            e = pwb[f"{ID}|beats_all={flag}"]
            exp[(ID, sub)] = (min(e["exp_eb_seL"], e["exp_eb_seA"]), max(e["exp_eb_seL"], e["exp_eb_seA"]))
    rows = []
    for _, r in comp.iterrows():
        lo, hi = exp.get((r.direction, r.subset), (float("nan"), float("nan")))
        rows.append({
            "direction": r.direction, "subset": r.subset, "hits (lockbox route)": int(r.primary_n_lockbox),
            "NOT EVALUABLE (primary)": int(r.primary_n_not_evaluable),
            "confirmed (primary)": int(r.primary_n_confirmed),
            "rate, primary evaluable": round(r.primary_rate_evaluable, 3),
            "rate, primary all-in": round(r.primary_rate_all_in, 3),
            "confirmed (corrected)": int(r.corrected_n_confirmed),
            "rate, corrected (non-confirmatory)": round(r.corrected_rate_evaluable, 3),
            "rate, corrected all-in": round(r.corrected_rate_all_in, 3),
            "defect-sensitive calls": int(r.n_defect_sensitive),
            "expected (post hoc power aid)": "" if pd.isna(lo) else f"{lo:.2f}–{hi:.2f}",
        })
    df = pd.DataFrame(rows)
    derived["power_band"] = {f"{k[0]}|{k[1]}": [round(v[0], 4), round(v[1], 4)] for k, v in sorted(exp.items())}
    derived["power_shrink_median"] = {ID: round(pw[ID]["median_shrink_eb"], 3) for ID in pw}
    rr = _tsv(f"{LB}/rates.tsv")
    derived["no_sets"] = {f"{x.entry}": [int(x.n_confirmed), int(x.n_lockbox)]
                          for x in rr[rr.subset == "no_sets"].itertuples() if x.n_lockbox > 0}
    write("tableS6_lockbox_rates", "Supplementary Table S6: lockbox confirmation rates per direction",
          ["Primary = the single read as executed (5bc9444) under ruling R7. Fits that failed because of the "
           "driver defect (the L covariate design was assembled over all models and then subset to L) are NOT "
           "EVALUABLE; 'evaluable' excludes them from the denominator, 'all-in' counts them as not confirmed. "
           "Three further fits failed for a genuine reason (CTRL_PC21 ×2, H1S HNF1B).",
           "Corrected = the disclosed, NON-CONFIRMATORY relook of all 560 affected fits (510 whose design changed, "
           "50 that had failed) through the registered resume path, after the fix and an intent commit (R6, R7). "
           "A call that differs between primary and corrected is defect-sensitive and carries no claim.",
           "H3b is always reported as evaluable / all-in / corrected together (D14): 26 of its 28 NOT EVALUABLE "
           "fits are in whole-genome doubling, which is effectively not evaluable in the primary.",
           "Expected (post hoc, descriptive; not registered, not a test): mean power of the one-sided L test over "
           "hits, with the A effect shrunk by an empirical-Bayes (Tweedie) posterior; the range spans the "
           "posterior mean / mean-plus-variance variants and the two SE choices (independent review §4). The gap "
           "between 'beats all matched sets' and the rest is about what this expectation predicts and is not "
           "evidence of complex-specificity.",
           "For H2 the subset 'beats_all_matched_sets' means: beats all 20 of its own random member splits.",
           "Rates exclude rare-context hits (H3b 11, H1 1, H1S 6), which are reported separately (Table S8). The "
           "'no sets' subset (H3a 30/62, H3b 3/6) is in results/part2/lockbox_confirm/rates.tsv.",
           "Sources: results/part2/lockbox_confirm/comparison_primary_vs_corrected_R7.tsv, rates.tsv; "
           "results/part2/verify/lockbox_results_review/tables_review.json (power, power_by_beats_all); script "
           "scripts/tab_lockbox.py."], df)
    return df


def h3_table(derived: dict, rel: pd.DataFrame) -> pd.DataFrame:
    mem = P.members()
    allx = []
    for ID, short in ((P.H3A, "H3a"), (P.H3B, "H3b")):
        x = h3_lockbox(ID, rel)
        s = x[x.strongest]
        lb = x[x.route == "lockbox"]
        derived[short] = {
            "n_discovery": len(x), "n_lockbox_route": len(lb),
            "n_confirmed_primary": int(x.confirmed_primary.sum()),
            "n_defect_sensitive": int(x.defect_sensitive.sum()),
            "n_rare": int((x.route != "lockbox").sum()), "n_rare_pass": int(x.rare_pass.sum()),
            "n_strongest": len(s), "n_strongest_units": int(s.unit.nunique()),
            "n_strongest_screen_and_gpd_kept": int((s.screen & s.cons).sum()),
            "n_strongest_control_unit_other_context": int(s.control.eq("control unit, other context").sum()),
            "n_strongest_lymphoid": int(s.context.eq("lineage:Lymphoid").sum()),
            "n_strongest_lymphoid_positive_beta": int(s.scale_caveat.sum()),
            "n_focus_beats_all_nonlocus": int(x.focus.sum()),
            "n_focus_confirmed_primary": int((x.focus & x.confirmed_primary).sum()),
            "n_focus_defect_sensitive": int((x.focus & x.defect_sensitive).sum()),
            "strongest_per_context": {k: int(v) for k, v in s.context.value_counts().sort_index().items()},
            "confirmed_per_context": {
                k: [int(g.confirmed_primary.sum()), int(len(g)), int(g.defect_sensitive.sum())]
                for k, g in lb.groupby("context")},
        }
        x = x.assign(direction=short, members=[P.gene_label(u, mem, 8) for u in x.unit])
        allx.append(x)
    t = pd.concat(allx)
    reg = _tsv("results/part2/robustness/strongest_findings.tsv")[["direction", "context", "unit", "rank"]]
    reg = reg.rename(columns={"rank": "registered_rank_S9"})
    t = t.merge(reg, on=["direction", "context", "unit"], how="left")
    t["registered_strongest_S9"] = t.registered_rank_S9.notna()
    t = t.rename(columns={"strongest": "earlier_unregistered_set"})
    t = t.assign(abs_t=t.t.abs()).sort_values(["direction", "registered_strongest_S9", "registered_rank_S9", "abs_t"],
                                              ascending=[True, False, True, False])
    cols = ["direction", "context", "unit", "members", "n_context", "beta", "ci95_lo", "ci95_hi", "t", "q_final",
            "locus", "beats_all", "control", "cons", "screen", "route", "n_context_L", "beta_L", "ci95_L_lo",
            "ci95_L_hi", "p_L", "primary_call", "corrected_call", "defect_sensitive", "beta_L_relook", "p_L_relook",
            "registered_strongest_S9", "registered_rank_S9", "earlier_unregistered_set", "scale_caveat"]
    t[cols].to_csv(S.TAB_DIR / "tableS7_lockbox_h3_full.tsv", sep="\t", index=False)
    top = []
    for short, k in (("H3a", 30), ("H3b", 45)):
        x = t[(t.direction == short) & t.registered_strongest_S9].head(k).copy()
        x["context"] = x.context.map(P.short_ctx)
        x["β_A [95% CI]"] = [S.fmt_ci(b, lo, hi, 3, True) for b, lo, hi in zip(x.beta, x.ci95_lo, x.ci95_hi)]
        x["β_L [95% CI]"] = [S.fmt_ci(b, lo, hi, 3, True) for b, lo, hi in zip(x.beta_L, x.ci95_L_lo, x.ci95_L_hi)]
        x["p_L"] = x.p_L.map(lambda v: f"{v:.1e}")
        x["lines A / L"] = [f"{a} / {int(b)}" for a, b in zip(x.n_context, x.n_context_L)]
        x["kept: GPD / screen"] = [f"{'yes' if a else 'no'} / {'yes' if b else 'no'}" for a, b in zip(x.cons, x.screen)]
        x["notes"] = [("; ".join(z for z in (c, "lymphoid β > 0: scale caveat" if sc else "") if z))
                      for c, sc in zip(x.control, x.scale_caveat)]
        top.append(x[["direction", "context", "unit", "members", "lines A / L", "β_A [95% CI]", "β_L [95% CI]",
                      "p_L", "kept: GPD / screen", "notes"]])
    md = pd.concat(top)
    write("tableS7_lockbox_h3", "Supplementary Table S7: lockbox results for H3a / H3b discoveries; rows shown are "
          "the registered strongest findings (Table S9)",
          ["Rows: the registered strongest findings (results/part2/robustness/strongest_findings.tsv; Supplementary "
           "Table S9), in registered rank order (tier, then A |t|): the first 30 of 207 H3a findings and all 45 H3b "
           "findings (rare-context-route rows have no L estimate). tableS7_lockbox_h3_full.tsv lists all 1,306 "
           "discoveries with every flag and both calls.",
           "Revision round 1: an earlier version of this view listed an unregistered 'strongest specific' subset "
           "(beats all matched sets; not a control readout). That subset is kept only as the column "
           "earlier_unregistered_set of the full TSV.",
           "β in Chronos units; negative = more dependent in the context. β_L is the same model refitted in the "
           "338 lockbox lines; its 95% CI is the HC3 interval, so for contexts with < 20 L lines (decided by the "
           "permutation p) it can include 0 while the call is 'confirmed'.",
           "'control unit, other context': the unit is a registered control readout in another context; not an "
           "independent finding. 'Lymphoid β > 0: scale caveat': lower lymphoid dependency on an essential "
           "unit is partly the lymphoid scale effect (main text, 'Two lineage-level artefacts').",
           "'kept: GPD / screen' = kept under the two A-side sensitivities (p = max(p_perm, analytic p); "
           "ScreenNNMD + CasActivity covariates).",
           "Members: first 8 member symbols of the frozen unit; any functional reading in the text is post hoc.",
           "Sources: results/part2/lockbox_confirm/{H3a_lineage_atlas,H3b_driver_family}/confirmation.tsv, "
           "primary_relabelled_R7.tsv.gz, corrected_refit_nonconfirmatory.tsv; A-side flags from "
           "results/part2/directions (via scripts/part2_report.py); script scripts/tab_lockbox.py."], md)
    return t


def h1_h2_bench_bridge(derived: dict, rel: pd.DataFrame) -> None:
    ev = _tsv(f"{LB}/h1_events_primary_vs_corrected_R7.tsv")
    e90 = _tsv(f"{LB}/H1_within_complex_GI_events.tsv")[["event", "axis", "ci90_L_lo", "ci90_L_hi",
                                                          "n_context_L"]]
    ev = ev.merge(e90, on=["event", "axis"], how="left")
    ci = ev.primary_ci95.str.strip("[]").str.split(",", expand=True).apply(pd.to_numeric, errors="coerce")
    ev["ci95_lo_num"], ev["ci95_hi_num"] = ci[0], ci[1]
    neu = ev[ev.label_final.eq("NEUTRAL") & ev.axis.eq("N")]
    derived["H1_lockbox"] = {
        "n_neutral_events": len(neu),
        "neutral_L_ci90_inside_margin": int(((neu.ci90_L_lo > -0.1) & (neu.ci90_L_hi < 0.1)).sum()),
        "neutral_L_ci95_excludes_0": int(((neu.ci95_lo_num > 0) | (neu.ci95_hi_num < 0)).sum()),
        "neutral_L_p_lt_005": int((pd.to_numeric(neu.primary_p_L, errors="coerce") < 0.05).sum()),
        "neutral_L_beta_range": [round(float(pd.to_numeric(neu.primary_beta_L).min()), 3),
                                 round(float(pd.to_numeric(neu.primary_beta_L).max()), 3)],
        "neutral_L_median_context_lines": float(neu.n_context_L.median()),
        "neutral_events": sorted(neu.event.str.split("|").str[1]),
    }
    rows = []
    for _, r in ev.sort_values(["eligible", "event"], ascending=[False, True]).iterrows():
        rows.append({
            "part": "H1", "item": f"{r.event.split('|')[1]} ({r.event.split('|')[0]}) {r.axis}",
            "A label": r.label_final, "route": r.route,
            "β/κ A": f"{r.beta_A:+.3f}",
            "L estimate [95% CI]": "" if pd.isna(r.primary_beta_L) else
            f"{r.primary_beta_L:+.3f} [{r.ci95_lo_num:+.3f}, {r.ci95_hi_num:+.3f}]",
            "L context lines": "" if pd.isna(r.n_context_L) else str(int(r.n_context_L)),
            "p_L": "" if pd.isna(r.primary_p_L) else f"{r.primary_p_L:.2g}",
            "primary call": r.primary_call, "corrected call (non-confirmatory)": r.corrected_call,
            "defect-sensitive": "yes" if bool(r.defect_sensitive) else "no",
            "promotable": "yes" if bool(r.promotable_primary) else "no",
        })
    # H2 by rank against own random splits
    h2 = P.h2_ranks()
    c2 = _tsv(f"{LB}/H2_submodule_divergence/confirmation.tsv")[["context", "unit", "beta_L", "p_L"]]
    r2 = rel[rel.entry == P.H2][["context", "unit", "primary_call", "defect_sensitive"]]
    h2 = h2.merge(c2, on=["context", "unit"]).merge(r2, on=["context", "unit"])
    for lab, g in (("rank 1 (beats all 20 splits)", h2[h2.rank_vs_splits == 1]),
                   ("rank ≥ 2", h2[h2.rank_vs_splits > 1])):
        ok = g[~g.primary_call.str.startswith("NOT EVALUABLE")]
        rows.append({"part": "H2 (exploratory)", "item": lab, "A label": f"{len(g)} discoveries",
                     "route": "lockbox", "primary call": f"{int(g.primary_call.eq('confirmed').sum())} confirmed "
                     f"of {len(ok)} evaluable", "defect-sensitive": str(int(g.defect_sensitive.sum()))})
    derived["H2_lockbox"] = {
        "rank1": [int(h2[h2.rank_vs_splits == 1].primary_call.eq("confirmed").sum()),
                  int((h2.rank_vs_splits == 1).sum())],
        "rank_ge2": [int(h2[h2.rank_vs_splits > 1].primary_call.eq("confirmed").sum()),
                     int((h2.rank_vs_splits > 1).sum())],
    }
    # BENCH
    b = S.load_json(f"{LB}/BENCH_complex_vs_gene/bench.json")
    br = S.load_json(f"{REV}/bench_recompute.json")
    pr = b["primary"]
    for k, lab in (("", "all pairs"), ("_low", "low breadth tercile"), ("_mid", "mid tercile"),
                   ("_high", "high tercile (registered criterion)")):
        ex = pr[f"ci95{k}"]
        asb = br["bootstrap_assembled"][f"ci95{k}"]
        rows.append({"part": "BENCH", "item": lab, "A label": "", "route": "",
                     "L estimate [95% CI]": f"Δ {pr['mean' + k]:+.3f} [{ex[0]:+.3f}, {ex[1]:+.3f}]",
                     "corrected call (non-confirmatory)": f"replicate-assembled covariates: [{asb[0]:+.3f}, "
                                                          f"{asb[1]:+.3f}]"})
    derived["BENCH"] = {
        "n_pairs": b["distinct"]["n_pairs"], "n_distinct_complexes": b["distinct"]["n_distinct_complexes_in_pairs"],
        "n_distinct_gstar": b["distinct"]["n_distinct_gstar"], "mean": pr["mean"], "ci95": pr["ci95"],
        "high": pr["mean_high"], "ci95_high": pr["ci95_high"], "ci95_high_assembled": br["bootstrap_assembled"]["ci95_high"],
        "low": pr["mean_low"], "ci95_low": pr["ci95_low"], "mid": pr["mean_mid"], "ci95_mid": pr["ci95_mid"],
        "rho_breadth": b["trend_spearman_delta_breadth"]["rho"], "frac_pair_rep_nan": br["bootstrap_executed"]["frac_pair_rep_nan"],
        "n_rep_with_any_nan": br["bootstrap_executed"]["n_rep_with_any_nan_pair"],
        "top_tercile_wins": pr["top_tercile_wins"],
        "multiverse": {k: {kk: v[kk] for kk in ("mean", "mean_high", "n_pairs")}
                       for k, v in b["multiverse"].items() if isinstance(v, dict)},
    }
    # BRIDGE-L
    bj = S.load_json(f"{LB}/bridge_primary_vs_corrected_R7.json")
    for key, lab in (("primary_evaluable", "primary (R7; NOT EVALUABLE excluded) — headline"),
                     ("as_executed", "as executed (failed fits counted not confirmed) — disclosure"),
                     ("corrected_nonconfirmatory", "corrected relook (non-confirmatory)")):
        p = bj[key]["primary"]
        rows.append({"part": "BRIDGE-L", "item": lab, "A label": "", "route": "",
                     "L estimate [95% CI]": f"zavg {p['rate_zavg']:.3f} vs GLS {p['rate_gls']:.3f}; Δ "
                                            f"{p['diff']:+.3f} [{p['ci95'][0]:+.3f}, {p['ci95'][1]:+.3f}]"})
    for mk, lab in (("unweighted", "multiverse: unweighted"), ("matched_gene_coverage", "multiverse: matched gene coverage"),
                    ("T1_T2_clusters", "multiverse: T1/T2 clusters only")):
        p = bj["primary_evaluable"]["multiverse"][mk]
        rows.append({"part": "BRIDGE-L", "item": lab, "A label": "", "route": "",
                     "L estimate [95% CI]": f"zavg {p['rate_zavg']:.3f} vs GLS {p['rate_gls']:.3f}; Δ "
                                            f"{p['diff']:+.3f} [{p['ci95'][0]:+.3f}, {p['ci95'][1]:+.3f}]"})
    derived["BRIDGE_L"] = {k: bj[k]["primary"] for k in bj}
    derived["BRIDGE_L"]["n_hits_primary_evaluable"] = bj["primary_evaluable"]["n_hits"]
    df = pd.DataFrame(rows).astype(object).fillna("")
    write("tableS8_lockbox_h1_h2_bench_bridge", "Supplementary Table S8: lockbox results for H1, H2, BENCH and "
          "BRIDGE-L",
          ["H1: every event axis with its A estimate (κ_N for N, β_P for P; Chronos units) and its L estimate with "
           "95% CI. Only the five non-neutral calls are eligible for confirmation; TP53 and PTEN carry the E2 "
           "label 'detectable, within equivalence margin' and cannot be promoted by a confirmation. EP300 P has "
           "4 L loss lines and takes the rare-context route (leave-one-lineage-out, leave-one-line-out and "
           "DEMETER2 on A, all same sign); this is not an independent replication. FOCAD N failed in the primary "
           "(driver defect); its corrected confirmation is defect-sensitive and carries no claim.",
           "H2 (exploratory): confirmation by rank of the discovery's |t| among its own 20 random member splits "
           "(Table S5). Confirmation in L shows that a divergence replicates, not that it is sub-module-specific.",
           "BENCH: Δ = t_L(complex score) − t_L(best member g*, chosen in A), signed toward the A direction; "
           "1,295 primary pairs in 353 complexes. Registered criterion: the complex wins if the high-breadth tercile "
           "CI excludes 0; it does not, so complex scores did not out-replicate the best gene. 16.5% of "
           "pair × replicate bootstrap fits fail and are dropped; with the covariate design assembled within each "
           "replicate, the high-tercile lower bound is −0.122 instead of −0.033 (independent review §5). The "
           "Spearman trend with breadth is reported without a p-value (pairs are not independent).",
           "BRIDGE-L: size-standardised lockbox confirmation rate of H3 discoveries made on the zavg map versus "
           "the GLS map (bootstrap B = 2,000, seed 5200). The R7 primary is the headline; the as-executed value is "
           "a disclosure. Every version favours GLS.",
           "Sources: results/part2/lockbox_confirm/{h1_events_primary_vs_corrected_R7.tsv, "
           "H1_within_complex_GI_events.tsv, H2_submodule_divergence/confirmation.tsv, BENCH_complex_vs_gene/"
           "bench.json, bridge_primary_vs_corrected_R7.json}; results/part2/verify/lockbox_results_review/"
           "bench_recompute.json; script scripts/tab_lockbox.py."], df)


def main() -> None:
    S.setup()
    derived: dict = {"label": "re-tabulation of committed lockbox outputs by scripts/tab_lockbox.py; not a trial, "
                              "not a test, no lockbox read"}
    rel = relabelled()
    derived["defect_sensitive_total"] = int(rel.defect_sensitive.sum())
    derived["defect_sensitive_by_entry"] = {k: int(v) for k, v in
                                            rel[rel.defect_sensitive].entry.value_counts().sort_index().items()}
    rates_table(derived)
    h3_table(derived, rel)
    h1_h2_bench_bridge(derived, rel)
    P.jdump(derived, "manuscript/tables/lockbox_derived.json")


if __name__ == "__main__":
    main()
