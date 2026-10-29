#!/usr/bin/env python
"""R-2 (referee M2): lineage scale effect with a non-linear trend. POST HOC, DESCRIPTIVE (results/part2/revision/
INTENT.md). Committed A-side files only: H3a tests.tsv.gz, robustness/unit_essentiality_A.tsv, the non-L columns of
robustness/strongest_findings.tsv and the never-expressed rows of calibration_D8_tests.tsv.gz.

Per H3a context: natural cubic spline (df = 4; ns() knots) OLS of beta on unit essentiality over all evaluable tested
units; residual z for the 207 H3a strongest findings; uniform rule scale_artefact_possible = concordant & |z| < 2 (any
context); corrected lead = T1, not a control readout, |z_nl| >= 2, not scale_artefact_possible. Never-expressed
p < 0.05 rate per context as a reliability column.
Outputs: results/part2/revision/r2_*.  Usage: .venv/bin/python scripts/59_rev_scale_effect.py
"""
from __future__ import annotations

import json
import resource
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import revision as rv  # noqa: E402

OUT = ROOT / "results" / "part2" / "revision"
ROB = ROOT / "results" / "part2" / "robustness"
H3A = ROOT / "results" / "part2" / "directions" / "H3a_lineage_atlas"
SF_COLS = ["rank", "direction", "context", "unit", "members", "tier", "tier_class", "beta_A", "t_A", "q_final",
           "control_recovery", "beats_all_matched_sets", "unit_essentiality_A", "scale_effect_fitted_beta",
           "scale_effect_concordant", "scale_effect_residual_z"]          # no L column is read
D16_ARTEFACT_CONTEXTS = ("lineage:Lymphoid", "lineage:Myeloid")


def main():
    t0 = time.time()
    tests = pd.read_csv(H3A / "tests.tsv.gz", sep="\t", usecols=["context", "unit", "beta", "evaluable"])
    ess = pd.read_csv(ROB / "unit_essentiality_A.tsv", sep="\t").set_index("unit").essentiality
    by_ctx, per_test = rv.spline_trend(tests, ess, df=4)
    sf = pd.read_csv(ROB / "strongest_findings.tsv", sep="\t", usecols=SF_COLS)
    h = sf[sf.direction == "H3a"].copy()
    h = h.merge(per_test[["context", "unit", "beta", "essentiality", "fitted_nl", "residual_beta_nl", "residual_z_nl",
                          "concordant_nl", "fitted_lin", "residual_z_lin", "concordant_lin"]],
                on=["context", "unit"], how="left", validate="one_to_one")
    assert np.allclose(h.beta, h.beta_A, atol=1e-5), "tests.tsv beta differs from the strongest-findings beta"
    assert np.allclose(h.residual_z_lin, h.scale_effect_residual_z, atol=1e-3), "linear residual z does not reproduce"
    readout = h.control_recovery.fillna("").astype(str).str.startswith("readout")
    conc_lin = h.scale_effect_concordant.astype(str).str.lower().eq("true")
    h["lead_D16"] = (h.tier == "T1") & ~readout & (h.scale_effect_residual_z.abs() >= 2)
    h["artefact_D16"] = h.context.isin(D16_ARTEFACT_CONTEXTS) & conc_lin & (h.scale_effect_residual_z.abs() < 2)
    lab = rv.scale_labels(h.residual_z_nl, h.concordant_nl, h.tier, h.control_recovery, h.direction)
    h["scale_artefact_possible"] = lab.scale_artefact_possible.values
    h["lead_corrected"] = lab.lead_corrected.values
    h["abs_residual_z_nl"] = h.residual_z_nl.abs()

    # never-expressed reliability column (committed calibration rows; H3a)
    ct = pd.read_csv(H3A / "calibration_D8_tests.tsv.gz", sep="\t")
    ne = ct[ct.calibration == "never_expressed"].groupby("context").p.agg(
        never_expressed_n=lambda p: int(np.isfinite(p).sum()), never_expressed_frac_p_lt_005=lambda p: float(
            np.mean(p[np.isfinite(p)] < 0.05)))
    cnt = h.groupby("context").agg(n_strongest=("unit", "size"),
                                   n_abs_z_nl_ge_2=("abs_residual_z_nl", lambda z: int((z >= 2).sum())),
                                   n_abs_z_lin_ge_2=("residual_z_lin", lambda z: int((z.abs() >= 2).sum())),
                                   n_scale_artefact_possible=("scale_artefact_possible", "sum"),
                                   n_artefact_D16=("artefact_D16", "sum"), n_lead_corrected=("lead_corrected", "sum"),
                                   n_lead_D16=("lead_D16", "sum"))
    tab = by_ctx.set_index("context").join(cnt, how="left").join(ne, how="left").reset_index()
    for c in cnt.columns:
        tab[c] = tab[c].fillna(0).astype(int)
    tab = tab.sort_values(["n_strongest", "context"], ascending=[False, True])
    tab.to_csv(OUT / "r2_trend_by_context.tsv", sep="\t", index=False, float_format="%.4f")

    cols = ["rank", "context", "unit", "tier", "beta_A", "t_A", "essentiality", "fitted_lin", "residual_z_lin",
            "concordant_lin", "fitted_nl", "residual_beta_nl", "residual_z_nl", "concordant_nl",
            "scale_artefact_possible", "lead_corrected", "lead_D16", "artefact_D16", "control_recovery",
            "beats_all_matched_sets", "members"]
    h = h.sort_values("rank")
    h[cols].to_csv(OUT / "r2_strongest_h3a.tsv", sep="\t", index=False, float_format="%.5g")
    leads = h[h.lead_corrected].sort_values("abs_residual_z_nl", ascending=False)
    leads[cols].to_csv(OUT / "r2_leads_corrected.tsv", sep="\t", index=False, float_format="%.5g")

    def xt(a, b):
        return {f"{x}->{y}": int(((h[a] == x) & (h[b] == y)).sum()) for x in (True, False) for y in (True, False)}
    summ = dict(
        label="POST HOC, DESCRIPTIVE (revision R-2); committed A-side files only",
        trend="natural cubic spline, df = 4 (boundary knots = range, interior knots at 25/50/75% quantiles), OLS, "
              "residual SD ddof = 5",
        n_h3a_strongest=int(len(h)), abs_residual_z_nl_ge_2=int((h.abs_residual_z_nl >= 2).sum()),
        abs_residual_z_lin_ge_2=int((h.residual_z_lin.abs() >= 2).sum()),
        concordant_nl=int(h.concordant_nl.sum()), concordant_lin=int(h.concordant_lin.sum()),
        scale_artefact_possible=int(h.scale_artefact_possible.sum()),
        lead_corrected=int(h.lead_corrected.sum()), lead_D16=int(h.lead_D16.sum()),
        artefact_D16=int(h.artefact_D16.sum()),
        lead_corrected_concordant_nl=int((h.lead_corrected & h.concordant_nl).sum()),
        lead_corrected_beats_all_matched_sets=int((h.lead_corrected & h.beats_all_matched_sets.astype(str)
                                                   .str.lower().eq("true")).sum()),
        transitions_lead=xt("lead_D16", "lead_corrected"),
        transitions_artefact=xt("artefact_D16", "scale_artefact_possible"),
        D16_leads_lost=h[h.lead_D16 & ~h.lead_corrected][["context", "unit", "residual_z_lin", "residual_z_nl"]]
        .round(3).values.tolist(),
        new_leads=h[~h.lead_D16 & h.lead_corrected][["context", "unit", "residual_z_lin", "residual_z_nl"]]
        .round(3).values.tolist(),
        per_context={r.context: dict(n=int(r.n_strongest), abs_z_nl_ge_2=int(r.n_abs_z_nl_ge_2),
                                     scale_artefact_possible=int(r.n_scale_artefact_possible),
                                     lead_corrected=int(r.n_lead_corrected),
                                     never_expressed_frac_p_lt_005=float(r.never_expressed_frac_p_lt_005))
                     for r in tab.itertuples() if r.n_strongest > 0},
        runtime_s=round(time.time() - t0, 2), peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    (OUT / "r2_summary.json").write_text(json.dumps(summ, indent=1, default=float) + "\n")
    print(json.dumps({k: v for k, v in summ.items() if k != "per_context"}, indent=1, default=float))
    print(tab[["context", "n_units", "spearman_rho", "r2_lin", "r2_nl", "n_strongest", "n_abs_z_nl_ge_2",
               "n_scale_artefact_possible", "n_lead_corrected", "n_lead_D16", "never_expressed_frac_p_lt_005"]]
          .to_string(index=False))
    print(leads[["rank", "context", "unit", "beta_A", "essentiality", "residual_z_lin", "residual_z_nl",
                 "concordant_nl"]].head(20).to_string(index=False))


if __name__ == "__main__":
    main()
