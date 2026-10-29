"""Reporting layer for H1 / H1S class calls (PI ruling R5, erratum E2; post hoc, conservative). Adds to classes.tsv:

  label_final           complex-wide axis: "DETECTABLE, WITHIN EQUIVALENCE MARGIN (negligible)" when a non-neutral
                        class (RESIDUAL SENSITISATION / COLLAPSE) AND the TOST-neutral criterion (90% CI within +/-0.1)
                        both hold; otherwise the engine's class. "NOT EVALUABLE" where the counted test is not evaluable
                        (ruling R2, applied in the reporting layer only for runs made before the fix: H1 is not re-run).
  label_final_paralog   the same rule for the paralog axis (PARALOG COMPENSATION and P_tost).
  headline_eligible     False for DETECTABLE-WITHIN-MARGIN labels (never headline; lockbox confirmation still runs on the
                        registered statistic but cannot promote them), for NOT EVALUABLE and for exploratory entries.
The engine's class columns, q-values and trial lines are unchanged. Also reports the H1 BH sensitivity of each
discovery to m (evaluable m vs the registered family_m.value).
Usage: .venv/bin/python scripts/47_h1_label_final.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import assoc, registry, stats  # noqa: E402

LABEL = "DETECTABLE, WITHIN EQUIVALENCE MARGIN (negligible)"
NONNEUTRAL = {"N": ("RESIDUAL SENSITISATION", "COLLAPSE"), "P": ("PARALOG COMPENSATION",)}


def label(cls, tost, evaluable, ax):
    if not evaluable:
        return "NOT EVALUABLE"
    if cls in NONNEUTRAL[ax] and bool(tost):
        return LABEL
    return cls


def main():
    proj = assoc.Project()
    out = {}
    for aid in ("H1_within_complex_GI", "H1S_silencing_GI"):
        od = proj.results_dir / "directions" / aid
        cls = pd.read_csv(od / "classes.tsv", sep="\t")
        cls = cls.drop(columns=[c for c in ("label_final", "label_final_paralog", "headline_eligible") if c in cls])
        tests = pd.read_csv(od / "tests.tsv.gz", sep="\t").set_index("unit")
        entry, _ = registry.load_entry(proj.registry_dir, aid, root=proj.root)
        expl = entry.get("standing") == "exploratory"
        ev = lambda u: bool(u in tests.index and tests.loc[u, "evaluable"])  # noqa: E731
        lab_n, lab_p = [], []
        for r in cls.itertuples(index=False):
            lab_n.append(r.class_complex_wide if r.class_complex_wide in ("not evaluable", "NOT EVALUABLE") else
                         label(r.class_complex_wide, getattr(r, "N_tost", False), ev(f"{r.event}|N"), "N"))
            lab_p.append(r.class_paralog if r.class_paralog in ("not evaluable", "NOT EVALUABLE") else
                         label(r.class_paralog, getattr(r, "P_tost", False), ev(f"{r.event}|P"), "P"))
        cls["label_final"], cls["label_final_paralog"] = lab_n, lab_p
        strong = {"RESIDUAL SENSITISATION", "COLLAPSE", "PARALOG COMPENSATION"}
        cls["headline_eligible"] = [(not expl) and (a in strong or b in strong) for a, b in zip(lab_n, lab_p)]
        cls.to_csv(od / "classes.tsv", sep="\t", index=False, float_format="%.6g")
        changed = cls[(cls.label_final != cls.class_complex_wide) | (cls.label_final_paralog != cls.class_paralog)]
        out[aid] = dict(changed=changed[["event", "class_complex_wide", "label_final", "class_paralog",
                                         "label_final_paralog"]].to_dict("records"),
                        label_counts=cls.label_final.value_counts().to_dict(),
                        headline_eligible=cls.event[cls.headline_eligible].tolist())
        # BH sensitivity to m: evaluable m (run) vs registered m (padding with p = 1 for missing tests)
        t = tests[tests.evaluable]
        m_reg = int(entry["family_m"]["value"])
        pf = t.p_final.values.astype(float)
        pad = np.concatenate([pf, np.ones(max(m_reg - len(pf), 0))])
        rej_reg, adj_reg = stats.bh(pad, float(entry["fdr_family"]["q"]))
        d = t.assign(q_at_registered_m=adj_reg[:len(pf)], discovery_at_registered_m=rej_reg[:len(pf)])
        d = d[d.discovery.astype(str) == "True"]
        out[aid]["m_sensitivity"] = dict(m_evaluable=int(len(pf)), m_registered=m_reg,
                                         discoveries=[dict(unit=u, q=float(r.q_final), q_at_registered_m=float(
                                             r.q_at_registered_m), kept=bool(r.discovery_at_registered_m))
                                                      for u, r in d.iterrows()])
    (proj.results_dir / "directions" / "h1_label_final.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
