#!/usr/bin/env python
"""Amendment A4 (PI rulings on results/part2/lockbox_confirm/GAPS.md): before any lockbox read; no L values read.

G1  BENCH_complex_vs_gene: both secondary components (A-L correlation vs matched random gene sets; number of distinct
    complexes among the g*) are declared NOT RUN (`estimand_secondary_not_run` + reason): under-specified, specifying
    them now adds analyst choices, and the random-set comparison needs L values of non-member genes (would widen the
    registered member-only BENCH read; same reasoning as A3(b)). Reported as limitations.
G2  BENCH_complex_vs_gene and BRIDGE-L_confirmation_rate: `multiverse.ci: point estimate only` (registry silent; the
    multiverse is descriptive). The primary estimands keep their registered CIs.
G3  H1S_silencing_GI: `lockbox.rare_route_fitted_parameters: not evaluable` -- the N-axis calls on non-lockbox-testable
    events, whose DEMETER2 rare route would need fitted baselines, are NOT EVALUABLE in L (no CRISPR or DEMETER2
    refit); listed in `lockbox.rare_route_not_evaluable`. H1S is exploratory.
The other literal readings of GAPS.md section 2 are accepted as written (record: docs/registry/AMENDMENTS.md, A4).

Usage: .venv/bin/python scripts/50_amend_A4.py [--check]
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
TIMESTAMP = "2026-10-29T02:42:37Z"                   # UTC time of the PI rulings record (docs/registry/AMENDMENTS.md)
AMENDS = "b349e98868ec8fb5d75b9f0f10e162ade9b6efe2"   # committed state amended (lockbox machinery + A3)
STATEMENT = "before any lockbox read; no L values read"
H1S_NOT_EVALUABLE = ["Z0018|IGSF11|N", "Z0018|SOX10|N", "Z0035|PAX8|N", "Z0084|NFATC2|N", "Z0149|SPDEF|N",
                     "Z0316|IGF1R|N"]
SECONDARY = ["A-L correlation of context effects vs matched random gene sets (seed 7100)",
             "number of distinct complexes among the g*"]
REASONS = {
    "BENCH_complex_vs_gene": (
        "PI rulings A4 (G1, G2) on results/part2/lockbox_confirm/GAPS.md, " + STATEMENT + ". G1: both BENCH secondary "
        "components (A-L correlation vs matched random sets; number of distinct complexes among the g*) are NOT RUN: "
        "under-specified, specifying them now would add analyst choices, and the random-set comparison needs L values of "
        "non-member genes, which would widen the registered member-only BENCH read (same reasoning as A3(b)); reported "
        "as limitations. G2: the multiverse reports point estimates only (registry silent; descriptive); the primary "
        "estimand keeps its registered CIs. Record: docs/registry/AMENDMENTS.md."),
    "BRIDGE-L_confirmation_rate": (
        "PI ruling A4 (G2) on results/part2/lockbox_confirm/GAPS.md, " + STATEMENT + ". The multiverse reports point "
        "estimates only (registry silent; descriptive); the primary estimand keeps its registered CI. Record: "
        "docs/registry/AMENDMENTS.md."),
    "H1S_silencing_GI": (
        "PI ruling A4 (G3) on results/part2/lockbox_confirm/GAPS.md, " + STATEMENT + ". The N-axis calls on "
        "non-lockbox-testable events would need fitted baselines on DEMETER2 for the rare-context route (registry "
        "silent); they are NOT EVALUABLE in L, with no CRISPR or DEMETER2 refit (H1S is exploratory). Record: "
        "docs/registry/AMENDMENTS.md."),
}


def git(*a) -> str:
    return subprocess.run(["git", "-C", str(ROOT), *a], capture_output=True, text=True, check=True).stdout


def _block_end(lines: list, header: str) -> int:
    i = lines.index(header)
    j = i + 1
    while j < len(lines) and (lines[j].startswith("  ") or not lines[j].strip()):
        j += 1
    return j


def _item(reason: str, paths: list) -> list:
    out = ["  - id: A4", f"    timestamp_utc: {TIMESTAMP}", f"    amends_commit: {AMENDS}", "    reason: >"]
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
        j = _block_end(lines, "multiverse:")
        lines = lines[:j] + ["  ci: point estimate only"] + lines[j:]
        a = lines.index("amendments:")
        lines = lines[:a] + ["estimand_secondary_not_run:"] + [f'  - "{s}"' for s in SECONDARY] + \
            ['estimand_secondary_not_run_reason: "PI ruling A4 (G1): under-specified; specifying it now would add '
             'analyst choices, and the random-set comparison needs L values of non-member genes (would widen the '
             'member-only BENCH read). Reported as a limitation."'] + lines[a:]
        paths = ["estimand_secondary_not_run", "estimand_secondary_not_run_reason", "multiverse.ci"]
    elif k == "BRIDGE-L_confirmation_rate":
        j = _block_end(lines, "multiverse:")
        lines = lines[:j] + ["  ci: point estimate only"] + lines[j:]
        paths = ["multiverse.ci"]
    else:
        j = _block_end(lines, "lockbox:")
        lines = lines[:j] + ["  rare_route_fitted_parameters: not evaluable", "  rare_route_not_evaluable:"] + \
            [f'    - "{u}"' for u in H1S_NOT_EVALUABLE] + lines[j:]
        paths = ["lockbox.rare_route_fitted_parameters", "lockbox.rare_route_not_evaluable"]
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
        print("A4 applied" if not bad else f"differs: {bad}")
        raise SystemExit(1 if bad else 0)
    for k, (t, paths, _, _) in p.items():
        (REG / f"{k}.yaml").write_text(t)
        print(k, paths)


if __name__ == "__main__":
    main()
