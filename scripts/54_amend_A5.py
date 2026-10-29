#!/usr/bin/env python
"""Amendment A5 (PI rulings on the pre-read code review, results/part2/verify/lockbox_code_review.md): before any
lockbox read; no L values read.

- Erratum to A4 G1's reason (BENCH_complex_vs_gene): the random-set comparison draws from frozen-complex members
  (calib.matched_random_sets), so it would not need L values of non-member genes. Under-specification remains the
  reason the A-L random-set correlation is NOT RUN.
- A4 reversed for the distinct-complex count among g*: it is A-only and is reported. `estimand_secondary_not_run`
  keeps only the A-L random-set correlation.
- G3 (H1S_silencing_GI): the six rare-route N axes use the A CRISPR alignment parameters, as the registry's
  "alignment from A baselines" prescribes (`lockbox.rare_route_fitted_parameters: A CRISPR parameters`, listed in
  `lockbox.rare_route_fitted_hits`); LOSO and leave-one-model-out run for them too. The A4 "NOT EVALUABLE" list is
  removed.
- BRIDGE-L_confirmation_rate: `multiverse.ci: bootstrap as primary` (costs < 1 min); the BENCH multiverse stays point
  estimates only.
Readings recorded in docs/registry/AMENDMENTS.md (A5): lineage-like controls (melanoma, heme, heme_vs_breast) get an
unrestricted PatientID-blocked permutation (decides PC17); H1/H1S permutations are PatientID-blocked within lineage;
breadth over all members (removed ones included) kept literal and disclosed; tie-breaking (N11) disclosed.

Usage: .venv/bin/python scripts/54_amend_A5.py [--check]
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import registry  # noqa: E402

REG = ROOT / "docs" / "registry"
TIMESTAMP = "2026-10-29T03:25:03Z"                   # UTC time of the PI rulings record (docs/registry/AMENDMENTS.md)
AMENDS = "b64109790478261569cdbc74d2eda6c31154cd4c"   # committed state amended (A4 + review + main merged)
STATEMENT = "before any lockbox read; no L values read"
H1S_FITTED = ["Z0018|IGSF11|N", "Z0018|SOX10|N", "Z0035|PAX8|N", "Z0084|NFATC2|N", "Z0149|SPDEF|N", "Z0316|IGF1R|N"]
BENCH_REASON_FIELD = ('estimand_secondary_not_run_reason: "PI rulings A4 (G1) and A5: under-specified (correlation '
                      'statistic, number of sets, matching tolerances, pair set and summary are not registered); '
                      'specifying it now would add analyst choices. Erratum (A5): the matched sets are drawn from '
                      'frozen-complex members, so the comparison would not need L values of non-member genes. '
                      'Reported as a limitation."')
REASONS = {
    "BENCH_complex_vs_gene": (
        "PI rulings A5 on the pre-read code review, " + STATEMENT + ". Erratum to A4 G1's reason: the random-set "
        "comparison draws from frozen-complex members (calib.matched_random_sets), so it would not need L values of "
        "non-member genes; under-specification remains the reason the A-L random-set correlation is NOT RUN. A4 is "
        "reversed for the number of distinct complexes among the g*, which is A-only and is reported. The BENCH "
        "multiverse stays point estimates only. Record: docs/registry/AMENDMENTS.md."),
    "BRIDGE-L_confirmation_rate": (
        "PI ruling A5 on the pre-read code review, " + STATEMENT + ". The BRIDGE-L multiverse reports bootstrap CIs "
        "(the primary resampling, B 2000, seed 5200; under 1 min), reversing A4 G2 for this entry. Record: "
        "docs/registry/AMENDMENTS.md."),
    "H1S_silencing_GI": (
        "PI ruling A5 (G3) on the pre-read code review, " + STATEMENT + ". The rare-route N axes of the six "
        "non-lockbox-testable H1S calls use the A CRISPR alignment parameters (baselines b_j and N'), as the "
        "registry's 'alignment from A baselines' prescribes; LOSO and leave-one-model-out run for them too. The A4 "
        "NOT EVALUABLE ruling is withdrawn. Record: docs/registry/AMENDMENTS.md."),
}


def git(*a) -> str:
    return subprocess.run(["git", "-C", str(ROOT), *a], capture_output=True, text=True, check=True).stdout


def _item(reason: str, paths: list) -> list:
    out = ["  - id: A5", f"    timestamp_utc: {TIMESTAMP}", f"    amends_commit: {AMENDS}", "    reason: >"]
    cur = "     "
    for w in reason.split():
        if len(cur) + len(w) + 1 > 118:
            out.append(cur)
            cur = "     "
        cur += " " + w
    out.append(cur)
    out.append("    changed_paths: [" + ", ".join(paths) + "]")
    return out


def amend_text(k: str, text: str) -> tuple[str, list]:
    lines = text.rstrip("\n").split("\n")
    if k == "BENCH_complex_vs_gene":
        i = lines.index("estimand_secondary_not_run:")
        j = next(n for n, l in enumerate(lines) if l.startswith("estimand_secondary_not_run_reason:"))
        lines = lines[:i] + ["estimand_secondary_not_run:",
                             '  - "A-L correlation of context effects vs matched random gene sets (seed 7100)"',
                             BENCH_REASON_FIELD] + lines[j + 1:]
        paths = ["estimand_secondary_not_run", "estimand_secondary_not_run_reason"]
    elif k == "BRIDGE-L_confirmation_rate":
        i = lines.index("  ci: point estimate only")
        lines[i] = "  ci: bootstrap as primary"
        paths = ["multiverse.ci"]
    else:
        i = lines.index("  rare_route_fitted_parameters: not evaluable")
        j = lines.index("  rare_route_not_evaluable:")
        end = j + 1
        while lines[end].startswith("    - "):
            end += 1
        lines = lines[:i] + ["  rare_route_fitted_parameters: A CRISPR parameters", "  rare_route_fitted_hits:"] + \
            [f'    - "{u}"' for u in H1S_FITTED] + lines[end:]
        paths = ["lockbox.rare_route_fitted_parameters", "lockbox.rare_route_not_evaluable",
                 "lockbox.rare_route_fitted_hits"]
    st = lines.index("status: open")
    lines = lines[:st] + _item(REASONS[k], paths) + lines[st:]
    return "\n".join(lines) + "\n", paths


def plan() -> dict:
    out = {}
    for k in REASONS:
        old_txt = git("show", f"{AMENDS}:docs/registry/{k}.yaml")
        old = registry.load_yaml(old_txt)
        new_txt, paths = amend_text(k, old_txt)
        new = registry.load_yaml(new_txt)
        errs = registry.validate(new, k)
        if errs:
            raise SystemExit(f"{k}: {errs}")
        if registry.amendment_changes(old, new):
            raise SystemExit(f"{k}: undeclared changes {registry.amendment_changes(old, new)}")
        if registry.changed_top_paths(old, new, paths) != set(paths):
            raise SystemExit(f"{k}: declared paths that did not change")
        out[k] = (new_txt, paths, old, new)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    p = plan()
    if a.check:
        bad = [k for k, (t, *_) in p.items() if (REG / f"{k}.yaml").read_text() != t]
        print("A5 applied" if not bad else f"differs: {bad}")
        raise SystemExit(1 if bad else 0)
    for k, (t, paths, _, _) in p.items():
        (REG / f"{k}.yaml").write_text(t)
        print(k, paths)


if __name__ == "__main__":
    main()
