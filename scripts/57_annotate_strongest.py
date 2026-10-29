#!/usr/bin/env python
"""Post hoc DESCRIPTIVE annotations of the strongest-findings list (PI ruling after the robustness run; not filters,
not tests). Adds to results/part2/robustness/strongest_findings.tsv:

  beats_all_matched_sets   H3: |t| above every matched random set of the unit in that context (calibration_D8_tests,
                           via scripts/part2_report.discoveries); H1: NA (H1 uses its own matched null; classes.tsv)
  control_recovery         'readout (PCx)': the (context, complex unit) pair is a registered complex-level control
                           readout, counting by-definition equivalent contexts (strongest.CONTROL_CONTEXT_OVERLAP, e.g.
                           H1loss:ARID1A = LoF:ARID1A); 'same unit, other control context (...)'; '' otherwise
  loso_partially_testable  >= 1 LOSO fold skipped (< 2 context or reference models left); loso_folds_tested / _skipped
  scale_effect_*           H3a only. Unit essentiality = mean over members of the gene's mean clipped A effect (as in
                           the calibration matching); per H3a context, OLS of beta on essentiality over ALL evaluable
                           tested units (tests.tsv.gz); concordant = sign(beta) == sign(fitted); residual beta and
                           residual z (residual / residual SD of that context)
  locus_flag, tier         for completeness

Unit essentiality needs the A gene-effect matrix (assoc.load_A_effects; A models only). `--essentiality` writes the
per-unit table results/part2/robustness/unit_essentiality_A.tsv once; the annotation step (default) then reads
committed files only.

Usage (PART1_DIR set, PYTHONHASHSEED=0):
    .venv/bin/python scripts/57_annotate_strongest.py --essentiality
    .venv/bin/python scripts/57_annotate_strongest.py
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from part2 import strongest  # noqa: E402

OUT = ROOT / "results" / "part2" / "robustness"
ESS = OUT / "unit_essentiality_A.tsv"
MAP = ROOT / "results" / "part2" / "fzx" / "zavg_engine_membership.tsv"
TESTS = ROOT / "results" / "part2" / "directions" / "H3a_lineage_atlas" / "tests.tsv.gz"
CONTROLS = ROOT / "results" / "part2" / "controls" / "controls_summary.tsv"


def write_essentiality():
    from part2 import assoc, calib
    from part2.cmap import FrozenMap
    fm = FrozenMap.from_membership(pd.read_csv(MAP, sep="\t"))
    E_A = assoc.load_A_effects()
    gene_ess = calib.clip_essentiality(E_A)
    units = fm.coarse()
    mem = {u: [g for g in fm.members[u] if g in E_A.columns] for u in units}
    ess = strongest.unit_essentiality(gene_ess, mem)
    pd.DataFrame(dict(unit=units, n_members_scored=[len(mem[u]) for u in units],
                      essentiality=ess.reindex(units).values)).to_csv(ESS, sep="\t", index=False, float_format="%.6g")
    print(f"{len(units)} units, A models {len(E_A)}")


def annotate():
    s = pd.read_csv(OUT / "strongest_findings.tsv", sep="\t")
    base_cols = [c for c in s.columns if c not in ANNOT_COLS]
    s = s[base_cols]
    per = pd.read_csv(OUT / "robustness_per_hit.tsv", sep="\t")[
        ["direction", "context", "unit", "locus", "beats_all"]]
    s = s.merge(per, on=["direction", "context", "unit"], how="left", validate="one_to_one")
    s["beats_all_matched_sets"] = np.where(s.direction == "H1", None, s.beats_all.astype(str).eq("True"))
    ctl = pd.read_csv(CONTROLS, sep="\t")
    ctl = ctl[ctl.level == "C"][["design_id", "context", "unit"]]
    s["control_recovery"] = [strongest.control_recovery(c, u, ctl) for c, u in zip(s.context, s.unit)]
    s["loso_folds_tested"] = s.loso_folds.astype(int)
    s["loso_folds_skipped"] = s.loso_skipped.astype(int)
    s["loso_partially_testable"] = s.loso_folds_skipped > 0
    s["locus_flag"] = s.locus.astype(str).eq("True")
    ess = pd.read_csv(ESS, sep="\t").set_index("unit").essentiality
    by_ctx, per_test = strongest.scale_effect(pd.read_csv(TESTS, sep="\t"), ess)
    pt = per_test.rename(columns={"essentiality": "unit_essentiality_A", "fitted": "scale_effect_fitted_beta",
                                  "residual_beta": "scale_effect_residual_beta",
                                  "residual_z": "scale_effect_residual_z"}).drop(columns="beta")
    pt["direction"] = "H3a"
    s = s.merge(pt, on=["direction", "context", "unit"], how="left", validate="one_to_one")
    s = s.merge(by_ctx[["context", "slope", "spearman_rho"]].rename(
        columns={"slope": "scale_effect_context_slope", "spearman_rho": "scale_effect_context_rho"}),
        on="context", how="left")
    s.loc[s.direction != "H3a", "scale_effect_context_slope"] = np.nan
    s.loc[s.direction != "H3a", "scale_effect_context_rho"] = np.nan
    s["scale_effect_concordant"] = s.scale_effect_concordant.astype(object).where(s.direction == "H3a", None)
    s = s.drop(columns=["locus", "beats_all"])
    s = s[base_cols + [c for c in ANNOT_COLS if c in s.columns]]
    s.to_csv(OUT / "strongest_findings.tsv", sep="\t", index=False, float_format="%.6g")
    by_ctx.to_csv(OUT / "scale_effect_by_context.tsv", sep="\t", index=False, float_format="%.6g")

    h = s[s.direction == "H3a"]
    conc = h.scale_effect_concordant.astype(bool)
    big = h.scale_effect_residual_z.abs() >= 2
    summ = dict(
        note="post hoc, descriptive; not filters, not tests",
        n=int(len(s)), locus_flag_true=int(s.locus_flag.sum()),
        tier={k: int(v) for k, v in s.tier.value_counts().sort_index().items()},
        beats_all_matched_sets={d: int(g.beats_all_matched_sets.astype(object).eq(True).sum())
                                for d, g in s.groupby("direction")},
        control_recovery={k or "none": int(v) for k, v in
                          s.control_recovery.str.split(" \\(").str[0].value_counts().items()},
        loso_partially_testable={d: int(g.loso_partially_testable.sum()) for d, g in s.groupby("direction")},
        scale_effect_H3a=dict(
            n=int(len(h)), concordant=int(conc.sum()), abs_residual_z_ge_2=int(big.sum()),
            concordant_and_abs_residual_z_ge_2=int((conc & big).sum()),
            per_context={c: dict(n=int(len(g)), concordant=int(g.scale_effect_concordant.astype(bool).sum()),
                                 abs_residual_z_ge_2=int((g.scale_effect_residual_z.abs() >= 2).sum()),
                                 slope=round(float(g.scale_effect_context_slope.iloc[0]), 4),
                                 rho=round(float(g.scale_effect_context_rho.iloc[0]), 3))
                         for c, g in h.groupby("context")}))
    (OUT / "annotations_summary.json").write_text(json.dumps(summ, indent=1, sort_keys=True) + "\n")
    print(json.dumps(summ, indent=1))
    v = h.assign(az=h.scale_effect_residual_z.abs()).sort_values(["tier_class", "az"], ascending=[True, False])
    print(v.head(15)[["rank", "context", "unit", "tier", "beta_A", "beta_L", "scale_effect_residual_z",
                      "scale_effect_concordant", "control_recovery", "members"]].to_string())


ANNOT_COLS = ["beats_all_matched_sets", "control_recovery", "loso_folds_tested", "loso_folds_skipped",
              "loso_partially_testable", "unit_essentiality_A", "scale_effect_context_slope",
              "scale_effect_context_rho", "scale_effect_fitted_beta", "scale_effect_concordant",
              "scale_effect_residual_beta", "scale_effect_residual_z", "locus_flag"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--essentiality", action="store_true")
    a = ap.parse_args()
    write_essentiality() if a.essentiality else annotate()


if __name__ == "__main__":
    main()
