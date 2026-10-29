#!/usr/bin/env python3
"""Supplementary Table S9: the 254 registered strongest findings with their post hoc annotations.

Input (committed): results/part2/robustness/strongest_findings.tsv (registered rule, D16; annotation columns from
scripts/57_annotate_strongest.py). Nothing is selected or re-ranked here: the 254 rows and their registered rank are
copied unchanged. Two descriptive reporting labels are added, both post hoc (RESEARCH_TREE.md D16 reporting rule):

- lead_example: H3a, tier T1, no control_recovery annotation, |scale_effect_residual_z| >= 2. These are the hits the
  manuscript names first.
- possible_lineage_scale_artefact: H3a in Lymphoid or Myeloid, scale_effect_concordant, |scale_effect_residual_z| < 2,
  i.e. the hit's sign and size are in line with the lineage-wide shift of beta with unit essentiality.

Outputs (deterministic): manuscript/tables/tableS9_strongest_findings.tsv (all 254 rows),
manuscript/tables/tableS9_strongest_findings.md (summary and the lead rows) and
manuscript/tables/strongest_derived.json (counts cited in the manuscript). Reads no dataset and no lockbox file.
"""
from __future__ import annotations

import json
import os
import re

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "results/part2/robustness/strongest_findings.tsv")
UNITS = os.path.join(ROOT, "results/part2/robustness/unit_essentiality_A.tsv")
OUT = os.path.join(ROOT, "manuscript/tables")
SCALE_LINEAGES = ("lineage:Lymphoid", "lineage:Myeloid")

COLS = ["rank", "direction", "context", "unit", "members", "tier", "beta_A", "t_A", "q_final", "route", "beta_L",
        "ci95_L_lo", "ci95_L_hi", "p_L", "n_context_L", "loso_unit", "loso_folds_tested", "loso_folds_skipped",
        "loso_partially_testable", "top2_robust", "breadth", "beats_all_matched_sets", "control_recovery",
        "unit_essentiality_A", "scale_effect_context_rho", "scale_effect_fitted_beta", "scale_effect_concordant",
        "scale_effect_residual_z", "lead_example", "possible_lineage_scale_artefact"]


def symbols(s: str) -> str:
    return ";".join(re.sub(r"\s*\(\d+\)$", "", x.strip()) for x in str(s).split(";"))


def fmt(x, nd=3):
    return "" if pd.isna(x) else f"{x:+.{nd}f}"


def main():
    d = pd.read_csv(SRC, sep="\t")
    assert len(d) == 254, len(d)
    d["members"] = d["members"].map(symbols)
    h3a = d.direction.eq("H3a")
    rz = d.scale_effect_residual_z.abs()
    noctl = d.control_recovery.isna() | d.control_recovery.eq("")
    d["lead_example"] = h3a & d.tier.eq("T1") & noctl & rz.ge(2)
    d["possible_lineage_scale_artefact"] = (h3a & d.context.isin(SCALE_LINEAGES)
                                            & d.scale_effect_concordant.eq(True) & rz.lt(2))
    t = d[COLS].sort_values("rank")
    t.to_csv(os.path.join(OUT, "tableS9_strongest_findings.tsv"), sep="\t", index=False, float_format="%.6g",
             lineterminator="\n")

    units = pd.read_csv(UNITS, sep="\t")
    lead = t[t.lead_example]
    art = t[t.possible_lineage_scale_artefact]
    lm = t[h3a & t.context.isin(SCALE_LINEAGES)]
    counts = {
        "n": int(len(t)),
        "per_direction": {k: int(v) for k, v in t.direction.value_counts().sort_index().items()},
        "units_per_direction": {k: int(v) for k, v in t.groupby("direction").unit.nunique().sort_index().items()},
        "tier": {k: int(v) for k, v in t.tier.value_counts().sort_index().items()},
        "lead_example": int(len(lead)),
        "lead_example_lymphoid_myeloid": int(lead.context.isin(SCALE_LINEAGES).sum()),
        "lead_example_per_context": {k: int(v) for k, v in lead.context.value_counts().sort_index().items()},
        "h3a_lymphoid_myeloid": int(len(lm)),
        "h3a_lymphoid_myeloid_concordant": int(lm.scale_effect_concordant.eq(True).sum()),
        "possible_lineage_scale_artefact": int(len(art)),
        "possible_lineage_scale_artefact_per_context": {k: int(v) for k, v in art.context.value_counts().sort_index().items()},
        "possible_lineage_scale_artefact_positive_beta": int((art.beta_A > 0).sum()),
        "median_unit_essentiality_artefact": round(float(art.unit_essentiality_A.median()), 3),
        "median_unit_essentiality_h3a_strongest": round(float(t[h3a].unit_essentiality_A.median()), 3),
        "median_unit_essentiality_all_units": round(float(units.essentiality.median()), 3),
        "h3b_T1_noncontrol": int((t.direction.eq("H3b") & t.tier.eq("T1") & noctl).sum()),
        "note": "post hoc, descriptive reporting labels (D16); they select and rank nothing",
    }
    with open(os.path.join(OUT, "strongest_derived.json"), "w") as fh:
        json.dump(counts, fh, indent=1, sort_keys=True)
        fh.write("\n")

    lines = [
        "# Supplementary Table S9: the 254 registered strongest findings, with post hoc annotations",
        "",
        "Registered rule (DESIGN §2.5; results/part2/robustness/README.md): discovery q < 0.05; primary lockbox call "
        "confirmed, or rare-context route passed; not defect-sensitive (R7); not a locus cluster; robust to "
        "leave-one-out by primary disease (H3a) or lineage (H3b, H1); for H1 also headline-eligible (E2). Ranked by "
        "tier class, then |t_A|. Computed on A after the lockbox read (disclosed).",
        "",
        f"Counts: {counts['per_direction'].get('H3a', 0)} H3a, {counts['per_direction'].get('H3b', 0)} H3b, "
        f"{counts['per_direction'].get('H1', 0)} H1; tiers T1 {counts['tier'].get('T1', 0)}, "
        f"T2 {counts['tier'].get('T2', 0)}, T3-neg {counts['tier'].get('T3-neg', 0)}.",
        "",
        "Annotation columns are post hoc and descriptive (none is a filter or a test): beats_all_matched_sets, "
        "control_recovery, loso_partially_testable, unit_essentiality_A and the scale_effect_* columns "
        "(scripts/57_annotate_strongest.py). The two reporting labels below are also post hoc (D16):",
        f"- lead_example ({counts['lead_example']} rows): H3a, T1, no control annotation, |residual z| >= 2.",
        f"- possible_lineage_scale_artefact ({counts['possible_lineage_scale_artefact']} rows): H3a in Lymphoid or "
        "Myeloid, scale-effect-concordant, |residual z| < 2 (beta in line with the lineage-wide essentiality trend; "
        f"median unit essentiality {counts['median_unit_essentiality_artefact']:+.2f} against "
        f"{counts['median_unit_essentiality_all_units']:+.2f} over all 367 units).",
        "",
        "Lead rows (full table: tableS9_strongest_findings.tsv). Members are listed in full in the TSV; beta in "
        "Chronos units, negative = more dependent in the context.",
        "",
        "| Rank | Context | Unit | Members (first 8) | A β | L β [95% CI] | n L | Residual z | Scale-concordant |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for _, r in lead.sort_values("rank").iterrows():
        mem = r.members.split(";")
        ms = ", ".join(mem[:8]) + ("…" if len(mem) > 8 else "")
        lines.append(f"| {r['rank']} | {r.context.replace('lineage:', '')} | {r.unit} | {ms} | {fmt(r.beta_A)} | "
                     f"{fmt(r.beta_L)} [{fmt(r.ci95_L_lo)}, {fmt(r.ci95_L_hi)}] | {int(r.n_context_L)} | "
                     f"{r.scale_effect_residual_z:+.2f} | {'yes' if r.scale_effect_concordant else 'no'} |")
    lines.append("")
    with open(os.path.join(OUT, "tableS9_strongest_findings.md"), "w") as fh:
        fh.write("\n".join(lines))


if __name__ == "__main__":
    main()
