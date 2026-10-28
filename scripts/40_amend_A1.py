#!/usr/bin/env python
"""FZ-X amendment A1 (PI decision, before any H1 run): split the H1 event list.

* H1_within_complex_GI keeps the LoF-type and CN-type (hom-del) events of the frozen FZ-X list (headline panel);
  its contexts, complex_set_hash and family_m.value are recomputed.
* H1S_silencing_GI (new, standing exploratory, never headlined) holds the silencing-type events with the same
  estimand, score, classes, covariates, TOST / matched-set rules and seeds, its own BH family.
Both entries carry an `amendments` record; registry.amendment_changes(FZ-X version, new) must be empty and the
declared paths must all have changed. Inputs: the FZ-X commit's files only (event list, entry text, manifest,
lockbox-testable list) - event counts and types, no dependency value, statistic or lockbox outcome.

  .venv/bin/python scripts/40_amend_A1.py --timestamp 2026-10-28T22:54:50Z
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import fzx, registry  # noqa: E402

FZX = ROOT / "results" / "part2" / "fzx"
REG = ROOT / "docs" / "registry"
H1, H1S = "H1_within_complex_GI", "H1S_silencing_GI"
PRIMARY_TYPES = ("LoF", "CN")
REASON = ("PI decision A1. DESIGN v3.1 H1 anticipated a panel of about 10-20 events; the registered silencing clause "
          "(TPM < 1 where the lineage median is >= 3) yields 318 events, mostly lineage-restricted non-expression. "
          "Such events are confounded with lineage composition, a non-expressed gene is trivially non-dependent, and "
          "at m = 366 the intended LoF/CN panel would be diluted. Decided from event counts and types only (no "
          "dependency values, association statistics or lockbox outcomes read). Verbatim rationale: "
          "docs/registry/AMENDMENTS.md.")


def git(*a) -> str:
    return subprocess.run(["git", "-C", str(ROOT), *a], capture_output=True, text=True, check=True).stdout


def sha(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def fzx_commit() -> str:
    """The commit that created the FZ-X manifest (the freeze this amendment applies to)."""
    c = git("log", "--diff-filter=A", "--format=%H", "--", "results/part2/fzx/FZX_manifest.json").split()
    if len(c) != 1:
        raise SystemExit(f"cannot identify the FZ-X commit: {c}")
    return c[0]


def units_of(ev: pd.DataFrame) -> list:
    return sorted([f"{r.event}|N" for r in ev.itertuples() if r.n_N >= 2] +
                  [f"{r.event}|P" for r in ev.itertuples() if r.n_P >= 1])


def amendment_block(c: str, ts: str, paths: list, extra: dict | None = None) -> list:
    lines = ["amendments:", "  - id: A1", f"    timestamp_utc: {ts}", f"    amends_commit: {c}"]
    lines += fzx.folded("reason", REASON, indent=4)
    lines.append("    changed_paths: [" + ", ".join(fzx.yaml_scalar(p) for p in paths) + "]")
    for k, v in (extra or {}).items():
        lines.append(f"    {k}: {fzx.yaml_scalar(v)}")
    return lines


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--timestamp", required=True, help="UTC time of the PI decision record, YYYY-MM-DDTHH:MM:SSZ")
    a = ap.parse_args()
    c = fzx_commit()
    man0 = json.loads(git("show", f"{c}:results/part2/fzx/FZX_manifest.json"))
    for f in ("results/part2/fzx/h1_events.tsv", "results/part2/fzx/h1_contexts.tsv"):
        if sha(ROOT / f) != man0["files"][f]:
            raise SystemExit(f"{f} differs from the FZ-X manifest")
    ev = pd.read_csv(FZX / "h1_events.tsv", sep="\t")
    ctx = pd.read_csv(FZX / "h1_contexts.tsv", sep="\t")
    ev["entry"] = [H1 if t in PRIMARY_TYPES else H1S for t in ev.type]
    if set(ev.groupby("gene").entry.nunique()) != {1}:
        raise SystemExit("a gene has events in both entries")
    ev_p, ev_s = ev[ev.entry == H1], ev[ev.entry == H1S]
    u_p, u_s = units_of(ev_p), units_of(ev_s)

    # ---- frozen lists of the amendment
    ev[["index", "event", "unit", "gene", "type", "locus", "testable", "entry"]].to_csv(
        FZX / "h1_event_split_A1.tsv", sep="\t", index=False)
    (FZX / "units_h1_components_A1.txt").write_text("\n".join(u_p) + "\n")
    (FZX / "units_h1s_components.txt").write_text("\n".join(u_s) + "\n")
    lb = pd.read_csv(pd.io.common.StringIO(git("show", f"{c}:results/part2/fzx/lockbox_testable.tsv")), sep="\t",
                     keep_default_na=False)
    to_s = set(ev_s.event)
    lb.loc[(lb.entry == H1) & lb.event.isin(to_s), "entry"] = H1S
    lb.loc[lb.entry == H1S, "direction"] = "H1"
    lb.to_csv(FZX / "lockbox_testable.tsv", sep="\t", index=False)

    # ---- registry: H1 (amended) and H1S (new), both from the FZ-X version of the H1 entry
    rel = f"docs/registry/{H1}.yaml"
    t0 = git("show", f"{c}:{rel}")
    e0 = registry.load_yaml(t0)
    items = e0["context_family"]["contexts"]
    genes_p = set(ev_p.gene)
    it_p = [x for x in items if x["id"].split(":", 1)[1] in genes_p]
    it_s = [x for x in items if x["id"].split(":", 1)[1] not in genes_p]
    if sorted(x["type"] for x in it_p if x["type"] not in PRIMARY_TYPES) or len(it_p) + len(it_s) != len(items):
        raise SystemExit("context split inconsistent with event types")

    def common(t, its, units):
        t = fzx.replace_sub(t, "context_family", "contexts", ["  contexts:"] + fzx.yaml_block(its, 4))
        t = fzx.replace_top(t, "complex_set_hash", [f"complex_set_hash: {registry.set_hash(units)}"])
        lines = t.split("\n")
        k = next(i for i, l in enumerate(lines) if l.startswith("family_m:"))
        v = next(i for i in range(k, len(lines)) if lines[i].startswith("  value:"))
        lines[v] = f"  value: {len(units)}"
        return "\n".join(lines)

    paths_p = ["context_family.contexts", "complex_set_hash", "family_m.value"]
    t1 = common(t0, it_p, u_p)
    t1 = fzx.insert_before_top(t1, "status", amendment_block(c, a.timestamp, paths_p))
    t1 = t1.replace("# H1 within-complex genetic interactions (DESIGN v3.1 2.6). Drafted at FZ-R.",
                    "# H1 within-complex genetic interactions (DESIGN v3.1 2.6). Drafted at FZ-R; amended A1 (LoF / "
                    "hom-del events only;\n# silencing events in H1S_silencing_GI).", 1)

    paths_s = ["id", "hypothesis", "context_family.name", "context_family.contexts", "complex_set_hash",
               "family_m.value", "kill_rule.rule", "lockbox.target_list", "success_criterion", "standing", "headline"]
    ts = common(t0, it_s, u_s)
    ts = fzx.replace_top(ts, "id", [f"id: {H1S}"])
    ts = fzx.replace_top(ts, "hypothesis", fzx.folded("hypothesis", (
        "EXPLORATORY companion of H1 (amendment A1; never headlined; labelled exploratory wherever reported). "
        "Silencing-type loss of a gene g that is a member of, or attached to, a frozen coarse complex C rewires "
        "dependency on the rest of C: paralog compensation, residual sensitisation, collapse, or no change (neutral "
        "by TOST). Estimand, score, classes, covariates, TOST and matched-set rules are those of "
        "H1_within_complex_GI; the silencing-type events of the FZ-X event list form their own BH family.")))
    ts = fzx.replace_sub(ts, "context_family", "name", [
        '  name: "H1S silencing-type loss events (complex C, silenced gene g), FZ-X list split by amendment A1"'])
    ts = fzx.replace_sub(ts, "kill_rule", "rule", fzx.folded("rule", (
        "H1S results never enter any kill tally or main claim. The H1 kill-rule verdict (controls PC6, PC7, PC8, "
        "PC12, PC15) is carried as a label: if H1 is killed by controls, H1S is reported as killed too. The run and "
        "the lockbox read still happen and are reported as exploratory."), indent=2))
    ts = fzx.replace_sub(ts, "lockbox", "target_list", [f"  target_list: results/part2/lockbox_targets/{H1S}.json"])
    ts = fzx.replace_top(ts, "success_criterion", fzx.folded("success_criterion", (
        "None (exploratory). Per event: class by the H1 rules; lockbox confirmations (or the labelled rare-context "
        "route) are run and reported for every non-neutral call, but H1S confirmations do not enter the kill rule, "
        "the strongest-findings list or the main claims. The full panel is reported, labelled exploratory.")))
    ts = fzx.insert_before_top(ts, "status", ["standing: exploratory", "headline: false"] +
                               amendment_block(c, a.timestamp, paths_s, {"derived_from": f"{H1} at {c}"}))
    ts = ts.split("\n", 2)
    ts = ("# H1S silencing-type within-complex genetic interactions: EXPLORATORY entry created by FZ-X amendment A1 "
          "from\n# H1_within_complex_GI (docs/registry/AMENDMENTS.md). Never headlined.\n" + ts[2])

    for k, txt, ref_paths in ((H1, t1, paths_p), (H1S, ts, paths_s)):
        new = registry.load_yaml(txt)
        errs = registry.validate(new, k)
        bad = registry.amendment_changes(e0, new)
        changed = registry.changed_top_paths(e0, new, ref_paths)
        if errs or bad or changed != set(ref_paths) or registry.placeholders(new):
            raise SystemExit(f"{k}: validation {errs} / undeclared changes {bad} / declared but unchanged "
                             f"{set(ref_paths) - changed}")
        if not txt.isascii():
            raise SystemExit(f"{k}: non-ASCII")
        (REG / f"{k}.yaml").write_text(txt)

    # ---- manifest
    man = json.loads(json.dumps(man0))                    # rebuilt from the FZ-X version (idempotent)
    for f in ("h1_event_split_A1.tsv", "units_h1_components_A1.txt", "units_h1s_components.txt",
              "lockbox_testable.tsv"):
        man["files"][f"results/part2/fzx/{f}"] = sha(FZX / f)
    for k in (H1, H1S):
        man["registry"][k] = dict(filled=man0["registry"][H1]["filled"] if k == H1 else [],
                                  amendment="A1", sha256=sha(REG / f"{k}.yaml"))
    s = man["summary"]
    h1_fzx = dict(s["m"]["H1"])
    s["m"]["H1"] = dict(m=len(u_p), n_events=len(ev_p), n_N_axis=int((ev_p.n_N >= 2).sum()),
                        n_P_axis=int((ev_p.n_P >= 1).sum()), complex_set_hash=registry.set_hash(u_p))
    s["m"]["H1S"] = dict(m=len(u_s), n_events=len(ev_s), n_N_axis=int((ev_s.n_N >= 2).sum()),
                         n_P_axis=int((ev_s.n_P >= 1).sum()), complex_set_hash=registry.set_hash(u_s),
                         standing="exploratory")
    s["lockbox_testable"] = {d: dict(testable=int((g.lockbox_testable.astype(str) == "True").sum()), total=int(len(g)))
                             for d, g in lb.groupby("entry")}
    man["complex_set_hashes"]["h1_components"] = registry.set_hash(u_p)
    man["complex_set_hashes"]["h1s_components"] = registry.set_hash(u_s)
    man["code"]["scripts/40_amend_A1.py"] = sha(ROOT / "scripts/40_amend_A1.py")
    man["code"]["scripts/39_freeze_fzx.py (after A1: refuses to overwrite an amended freeze)"] = \
        sha(ROOT / "scripts/39_freeze_fzx.py")
    man["code"]["src/part2/registry.py"] = sha(ROOT / "src/part2/registry.py")
    man["code"]["src/part2/fzx.py"] = sha(ROOT / "src/part2/fzx.py")
    man["amendments"] = [dict(
        id="A1", timestamp_utc=a.timestamp, amends_commit=c, record="docs/registry/AMENDMENTS.md",
        decision="H1 keeps LoF / CN events; silencing events moved to the exploratory entry H1S_silencing_GI",
        H1_at_fzx=h1_fzx, H1_after=s["m"]["H1"], H1S=s["m"]["H1S"],
        inputs_read="FZ-X event list (counts and types) only; no dependency values, association statistics or "
                    "lockbox outcomes")]
    (FZX / "FZX_manifest.json").write_text(json.dumps(man, indent=1) + "\n")
    print(json.dumps(dict(H1=s["m"]["H1"], H1S=s["m"]["H1S"], fzx_commit=c), indent=1))


if __name__ == "__main__":
    main()
