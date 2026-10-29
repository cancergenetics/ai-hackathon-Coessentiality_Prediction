#!/usr/bin/env python
"""Amendment A3 (PI rulings on results/part2/lockbox_confirm/BLOCKED.md): before any lockbox read; no L values read.

(a) Lockbox permutation seeds (the exact PatientID-blocked Freedman-Lane permutation of a confirmation in L, used by
    confirm.confirm_entry as seed + hit index, hits in (context, unit) order): H3a 6100, H3b 6200, H1 6300, H1S 6400,
    H2 6500, BRIDGE-L 6600, CTRL entries 6700 + registry index (the index of the registered control permutation seed
    4400 + index: CTRL_PC1 -> 0 ... CTRL_PC22 -> 22). Written as `seeds.lockbox_permutation`. Checked against every
    registered seed and every seed named in DESIGN.md / RESEARCH_TREE.md / docs/registry: no collision, so no block
    was shifted (+50 rule not needed).
(b) BENCH_complex_vs_gene multiverse score `within_line_rank` is declared NOT RUN (`multiverse.not_run`,
    `multiverse.not_run_reason`): its registered computation needs genome-wide within-line ranks in L, while the
    registered BENCH read covers complex members only; the PI chose not to widen the read (reported as a limitation).

Each touched entry gets an `amendments` item (id A3, UTC timestamp, amends_commit = the committed state amended,
reason, changed_paths); registry.amendment_changes(committed, new) must be [] and every declared path must change.
Text edits only (the registry YAML subset has no writer); the result is re-parsed and validated before writing.

Usage: .venv/bin/python scripts/52_amend_A3.py [--check]
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import registry  # noqa: E402

REG = ROOT / "docs" / "registry"
TIMESTAMP = "2026-10-29T02:03:27Z"                   # UTC time of the PI rulings record (docs/registry/AMENDMENTS.md)
AMENDS = "21fe1be1cdfe9d9b5b68c02c1b5a46b2b5b66a69"   # committed state amended (dry-check stop note; registry = 75bf5e6)
STATEMENT = "before any lockbox read; no L values read"
LOCKBOX_SEEDS = {"H3a_lineage_atlas": 6100, "H3b_driver_family": 6200, "H1_within_complex_GI": 6300,
                 "H1S_silencing_GI": 6400, "H2_submodule_divergence": 6500, "BRIDGE-L_confirmation_rate": 6600}
CTRL_BASE, CTRL_REGISTERED_BASE = 6700, 4400
SHIFT = 50
REASON_SEEDS = ("PI ruling A3(a) on the lockbox dry check (results/part2/lockbox_confirm/BLOCKED.md), " + STATEMENT +
                ". The lockbox confirmation permutation seed was not registered: lockbox_permutation is the seed of the "
                "exact PatientID-blocked Freedman-Lane permutation in L (confirm.confirm_entry: seed + hit index, hits "
                "in (context, unit) order). H3a 6100, H3b 6200, H1 6300, H1S 6400, H2 6500, BRIDGE-L 6600, CTRL 6700 + "
                "registry index; checked against every registered or documented seed, no collision. "
                "Record: docs/registry/AMENDMENTS.md.")
REASON_BENCH = ("PI ruling A3(b) on the lockbox dry check, " + STATEMENT + ". The multiverse score within_line_rank needs "
                "genome-wide within-line ranks in L, but the registered BENCH read covers complex members only, so the "
                "entry is internally inconsistent; the PI chose not to widen the lockbox read and within_line_rank is "
                "declared NOT RUN (reported as a limitation). Record: docs/registry/AMENDMENTS.md.")
BENCH_NOT_RUN_REASON = ("needs genome-wide within-line ranks in L; the registered BENCH read covers complex members "
                        "only (amendment A3(b), PI ruling: the read is not widened; reported as a limitation)")


def git(*a) -> str:
    return subprocess.run(["git", "-C", str(ROOT), *a], capture_output=True, text=True, check=True).stdout


def ctrl_ids() -> list:
    return sorted(p.stem for p in REG.glob("CTRL_*.yaml"))


def planned_seeds(committed: dict) -> dict:
    """{entry: lockbox_permutation} with the +SHIFT rule applied per block on a collision."""
    out = dict(LOCKBOX_SEEDS)
    for k in ctrl_ids():
        out[k] = CTRL_BASE + int(committed[k]["seeds"]["permutation"]) - CTRL_REGISTERED_BASE
    return out


def registered_seeds(entries: dict) -> set:
    vals = set()
    for e in entries.values():
        for v in (e.get("seeds") or {}).values():
            if isinstance(v, int):
                vals.add(v)
        cal = e.get("calibration") or {}
        if isinstance(cal.get("seed"), int):
            vals.add(cal["seed"])
    return vals


def documented_seeds(commit: str = AMENDS) -> set:
    """Integers named as seeds in DESIGN.md, RESEARCH_TREE.md and docs/registry/*.md as of the amended commit (before
    the A3 record itself), e.g. 'seed 3000', 'seed is 7000', 'seeds 4000+'."""
    vals = set()
    names = [l for l in git("ls-tree", "--name-only", commit, "docs/registry/").split() if l.endswith(".md")]
    for rel in ["DESIGN.md", "RESEARCH_TREE.md", *sorted(names)]:
        for m in re.finditer(r"seeds?\b[^0-9\n]{0,24}(\d{3,5})", git("show", f"{commit}:{rel}"), flags=re.I):
            vals.add(int(m.group(1)))
    return vals


def collisions(seeds: dict, taken: set) -> list:
    return sorted(k for k, v in seeds.items() if v in taken)


def resolve(committed: dict) -> tuple[dict, list]:
    """Apply the PI's +50-per-block rule until no block collides; returns (seeds, shifted blocks)."""
    seeds = planned_seeds(committed)
    taken = registered_seeds(committed) | documented_seeds()
    blocks = {k: (k if not k.startswith("CTRL_") else "CTRL") for k in seeds}
    shifted = []
    for _ in range(20):
        bad = {blocks[k] for k in collisions(seeds, taken)}
        if not bad:
            break
        for b in sorted(bad):
            shifted.append(b)
            for k in seeds:
                if blocks[k] == b:
                    seeds[k] += SHIFT
    return seeds, shifted


def _insert_after_block(lines: list, header: str, new: list) -> list:
    i = lines.index(header)
    j = i + 1
    while j < len(lines) and (lines[j].startswith("  ") or not lines[j].strip()):
        j += 1
    return lines[:j] + new + lines[j:]


def _amendment_item(reason: str, paths: list) -> list:
    out = ["  - id: A3", f"    timestamp_utc: {TIMESTAMP}", f"    amends_commit: {AMENDS}", "    reason: >"]
    words, cur = reason.split(), "     "
    for w in words:
        if len(cur) + len(w) + 1 > 118:
            out.append(cur)
            cur = "     "
        cur += " " + w
    out.append(cur)
    out.append("    changed_paths: [" + ", ".join(paths) + "]")
    return out


def amend_text(text: str, seed: int | None, bench: bool) -> tuple[str, list]:
    lines = text.rstrip("\n").split("\n")
    paths = []
    if seed is not None:
        lines = _insert_after_block(lines, "seeds:", [f"  lockbox_permutation: {seed}"])
        paths.append("seeds.lockbox_permutation")
    if bench:
        lines = _insert_after_block(lines, "multiverse:", ["  not_run: [within_line_rank]",
                                                           f"  not_run_reason: \"{BENCH_NOT_RUN_REASON}\""])
        paths += ["multiverse.not_run", "multiverse.not_run_reason"]
    reason = REASON_BENCH if bench else REASON_SEEDS
    item = _amendment_item(reason, paths)
    st = lines.index("status: open") if "status: open" in lines else lines.index("status: closed")
    if "amendments:" in lines:
        lines = lines[:st] + item + lines[st:]
    else:
        lines = lines[:st] + ["amendments:"] + item + lines[st:]
    return "\n".join(lines) + "\n", paths


def plan() -> dict:
    """{entry: (new_text, declared paths, committed entry, new entry)}; raises on any validation failure."""
    committed = {p.stem: registry.load_yaml(git("show", f"{AMENDS}:docs/registry/{p.name}"))
                 for p in sorted(REG.glob("*.yaml"))}
    seeds, shifted = resolve(committed)
    if shifted:
        print(f"seed blocks shifted by +{SHIFT}: {shifted}")
    out = {}
    for k in sorted(committed):
        is_bench = k == "BENCH_complex_vs_gene"
        if k not in seeds and not is_bench:
            continue
        old_txt = git("show", f"{AMENDS}:docs/registry/{k}.yaml")
        new_txt, paths = amend_text(old_txt, seeds.get(k), is_bench)
        new = registry.load_yaml(new_txt)
        errs = registry.validate(new, k)
        if errs:
            raise SystemExit(f"{k}: {errs}")
        if registry.amendment_changes(committed[k], new):
            raise SystemExit(f"{k}: undeclared changes {registry.amendment_changes(committed[k], new)}")
        if registry.changed_top_paths(committed[k], new, paths) != set(paths):
            raise SystemExit(f"{k}: declared paths that did not change")
        out[k] = (new_txt, paths, committed[k], new)
    return out


def applied_commit() -> str:
    """The first commit whose H3a entry carries the A3 record."""
    log = git("log", "--reverse", "--format=%H", "-S", "  - id: A3", "--", "docs/registry/H3a_lineage_atlas.yaml")
    return log.split()[0] if log.split() else "HEAD"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="verify the working tree equals the planned amendment")
    a = ap.parse_args()
    p = plan()
    if a.check:                            # against the commit that applied A3 (A4+ may edit the same files later)
        c = applied_commit()
        bad = [k for k, (t, *_) in p.items() if git("show", f"{c}:docs/registry/{k}.yaml") != t]
        print("A3 applied" if not bad else f"differs: {bad}")
        raise SystemExit(1 if bad else 0)
    for k, (t, paths, _, new) in p.items():
        (REG / f"{k}.yaml").write_text(t)
        print(k, paths, new.get("seeds", {}).get("lockbox_permutation"))


if __name__ == "__main__":
    main()
