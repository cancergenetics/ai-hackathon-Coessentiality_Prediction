#!/usr/bin/env python
"""FZ-X amendment A2 (PI ruling on the FZ-X review M1 and note b, before any H1 run).

Rule: an H1 event is GENETIC (H1_within_complex_GI, primary) iff >= 50% of its A loss models carry a LoF call
(passing the functional check) or a hom-del call (combined); otherwise it is SILENCING (H1S_silencing_GI,
exploratory). Applied mechanically to all FZ-X events from the frozen per-gene loss-kind counts
(h1_contexts.tsv: n_A_LoF, n_A_homdel, n_A; each A loss model carries exactly one kind, hom-del first, then LoF,
then silencing). Genetic events carry a descriptive sub-label (CN-majority if hom-del calls are >= 50% of the genetic
loss models, else LoF-majority) for reporting only. The analysis specification of every event (covariates, N / P
components; FZ-X rule-6 `type`) is unchanged.
Also (note b): H1S gets its own direction label `H1S` and an unambiguous fdr_rule; the stale `expected_events` and
`implementation_note` texts of both entries are replaced.

Inputs: the A1-amended outputs (scripts/41_freeze_fzx.py re-applies 42 then 43). Event counts and types only.
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
TIMESTAMP = "2026-10-28T23:28:31Z"                      # UTC time of the PI ruling record (docs/registry/AMENDMENTS.md)
GENETIC_MIN = 0.5
REASON = ("PI ruling A2 on FZ-X review M1. DESIGN 2.2 counts hom-del as LoF, and A1 intended genetic loss as the "
          "primary panel and silencing as exploratory. An event is genetic (H1) iff >= 50% of its A loss models carry "
          "a LoF or hom-del call (combined); otherwise silencing (H1S). Applied mechanically to all FZ-X events; "
          "genetic events keep a descriptive LoF-majority / CN-majority sub-label for reporting only; covariates and "
          "components stay as frozen at FZ-X. Also: H1S direction label and fdr_rule wording, stale text removed "
          "(review note b). Decided from event counts and types only (no dependency values, association statistics "
          "or lockbox outcomes read). Record: docs/registry/AMENDMENTS.md.")
PATHS = {H1: ["context_family.contexts", "context_family.expected_events", "complex_set_hash", "family_m.value",
              "implementation_note"],
         H1S: ["direction", "context_family.contexts", "context_family.expected_events", "complex_set_hash",
               "family_m.value", "fdr_rule", "implementation_note"]}


def git(*a) -> str:
    return subprocess.run(["git", "-C", str(ROOT), *a], capture_output=True, text=True, check=True).stdout


def sha(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def a1_commit() -> str:
    """The commit that created H1S_silencing_GI (amendment A1): the committed state A2 amends."""
    c = git("log", "--diff-filter=A", "--format=%H", "--", f"docs/registry/{H1S}.yaml").split()
    if len(c) != 1:
        raise SystemExit(f"cannot identify the A1 commit: {c}")
    return c[0]


def units_of(ev: pd.DataFrame) -> list:
    return sorted([f"{r.event}|N" for r in ev.itertuples() if r.n_N >= 2] +
                  [f"{r.event}|P" for r in ev.itertuples() if r.n_P >= 1])


def classify(ctx: pd.DataFrame) -> pd.DataFrame:
    """Per lost gene: genetic share, loss class and descriptive sub-label (A2 rule)."""
    g = ctx.n_A_LoF + ctx.n_A_homdel
    out = pd.DataFrame({"context": ctx.context, "gene": ctx.gene, "fzx_type": ctx.type, "n_A": ctx.n_A,
                        "n_A_LoF": ctx.n_A_LoF, "n_A_homdel": ctx.n_A_homdel, "n_A_silencing": ctx.n_A_silencing})
    out["genetic_share"] = (g / ctx.n_A).round(4)
    out["loss_class"] = ["genetic" if x >= GENETIC_MIN else "silencing" for x in g / ctx.n_A]
    out["genetic_sublabel"] = ["" if c != "genetic" else ("CN-majority" if h >= 0.5 * t else "LoF-majority")
                               for c, h, t in zip(out.loss_class, ctx.n_A_homdel, g)]
    return out


def a2_item(c: str, ts: str, paths: list) -> list:
    lines = ["  - id: A2", f"    timestamp_utc: {ts}", f"    amends_commit: {c}"]
    lines += fzx.folded("reason", REASON, indent=4)
    lines.append("    changed_paths: [" + ", ".join(fzx.yaml_scalar(p) for p in paths) + "]")
    return lines


def apply(timestamp: str = TIMESTAMP):
    c = a1_commit()
    man = json.loads((FZX / "FZX_manifest.json").read_text())
    if [x["id"] for x in man.get("amendments", [])] != ["A1"]:
        raise SystemExit("A2 applies on top of A1 only; run scripts/41_freeze_fzx.py (it re-applies A1 and A2)")
    for f in ("h1_events.tsv", "h1_contexts.tsv", "h1_event_split_A1.tsv"):
        if sha(FZX / f) != man["files"][f"results/part2/fzx/{f}"]:
            raise SystemExit(f"{f} differs from the manifest")
    ev = pd.read_csv(FZX / "h1_events.tsv", sep="\t")
    ctx = pd.read_csv(FZX / "h1_contexts.tsv", sep="\t")
    s1 = pd.read_csv(FZX / "h1_event_split_A1.tsv", sep="\t")
    cl = classify(ctx).set_index("gene")
    ev["loss_class"] = cl.loss_class.reindex(ev.gene).values
    ev["genetic_sublabel"] = cl.genetic_sublabel.reindex(ev.gene).values
    ev["genetic_share"] = cl.genetic_share.reindex(ev.gene).values
    ev["entry_A1"] = s1.set_index("event").entry.reindex(ev.event).values
    ev["entry"] = [H1 if x == "genetic" else H1S for x in ev.loss_class]
    moved = ev[ev.entry != ev.entry_A1]
    ev_p, ev_s = ev[ev.entry == H1], ev[ev.entry == H1S]
    u_p, u_s = units_of(ev_p), units_of(ev_s)
    ev[["index", "event", "unit", "gene", "type", "loss_class", "genetic_sublabel", "genetic_share", "n_A",
        "locus", "testable", "entry_A1", "entry"]].rename(columns={"type": "fzx_type"}).to_csv(
        FZX / "h1_event_split_A2.tsv", sep="\t", index=False)
    (FZX / "units_h1_components_A2.txt").write_text("\n".join(u_p) + "\n")
    (FZX / "units_h1s_components_A2.txt").write_text("\n".join(u_s) + "\n")
    lb = pd.read_csv(FZX / "lockbox_testable.tsv", sep="\t", keep_default_na=False)
    is_h1 = lb.entry.isin([H1, H1S])
    lb.loc[is_h1, "entry"] = ev.set_index("event").entry.reindex(lb.event[is_h1]).values
    lb.loc[lb.entry == H1S, "direction"] = "H1S"
    lb.to_csv(FZX / "lockbox_testable.tsv", sep="\t", index=False)

    # ---- registry
    texts = {k: (REG / f"{k}.yaml").read_text() for k in (H1, H1S)}
    olds = {k: registry.load_yaml(t) for k, t in texts.items()}
    for k, e in olds.items():
        if [x["id"] for x in e.get("amendments", [])] != ["A1"]:
            raise SystemExit(f"{k}: expected the A1 state")
    items = {x["id"]: x for k in (H1, H1S) for x in olds[k]["context_family"]["contexts"]}
    if sorted(items) != sorted(ctx.context):
        raise SystemExit("registry H1/H1S contexts differ from h1_contexts.tsv")

    def item(cid):
        x = dict(items[cid])
        r = cl.loc[cid.split(":", 1)[1]]
        x["loss_class"] = r.loss_class
        x["genetic_share"] = float(r.genetic_share)
        if r.loss_class == "genetic":
            x["genetic_sublabel"] = r.genetic_sublabel
        return x

    n_tot = len(ev)
    new_texts = {}
    for k, evk, units in ((H1, ev_p, u_p), (H1S, ev_s, u_s)):
        t = texts[k]
        cids = [f"H1loss:{g}" for g in ctx.gene if f"H1loss:{g}" in set("H1loss:" + evk.gene)]
        t = fzx.replace_sub(t, "context_family", "contexts", ["  contexts:"] + fzx.yaml_block([item(x) for x in cids], 4))
        t = fzx.replace_sub(t, "context_family", "expected_events", [
            f'  expected_events: "FZ-X enumeration rule: {n_tot} events; amendments A1/A2 (docs/registry/AMENDMENTS.md):'
            f' {len(ev_p)} genetic events in {H1}, {len(ev_s)} silencing events in {H1S}; this entry: {len(evk)}"'])
        t = fzx.replace_top(t, "complex_set_hash", [f"complex_set_hash: {registry.set_hash(units)}"])
        lines = t.split("\n")
        i = next(j for j, l in enumerate(lines) if l.startswith("family_m:"))
        v = next(j for j in range(i, len(lines)) if lines[j].startswith("  value:"))
        lines[v] = f"  value: {len(units)}"
        t = "\n".join(lines)
        t = fzx.replace_top(t, "implementation_note", fzx.folded("implementation_note", (
            "Implemented and tested: h1.assoc_inputs builds the per-context score frames and the (context, unit) "
            "pairs, and run_association(pairs=...) tests each loss column only against its own event components, "
            "with BH pooled over the counted tests of this entry (tests/test_part2_engine.py, "
            "tests/part2_synth_pipeline.py).")))
        if k == H1S:
            t = fzx.replace_top(t, "direction", ["direction: H1S"])
            t = fzx.replace_top(t, "fdr_rule", fzx.folded("fdr_rule", (
                "BH at q = 0.05 pooled over every counted test of THIS entry only (H1S_silencing_GI: all its events, "
                "both axes; m = family_m.value). It is a family of its own, separate from H1_within_complex_GI. The "
                "class conditions (contrast, breadth, matched null, TOST) are applied after BH and cannot rescue a "
                "non-rejected test.")))
        t = fzx.insert_before_top(t, "status", a2_item(c, timestamp, PATHS[k]))
        new = registry.load_yaml(t)
        errs = registry.validate(new, k)
        bad = registry.amendment_changes(olds[k], new)
        changed = registry.changed_top_paths(olds[k], new, PATHS[k])
        if errs or bad or changed != set(PATHS[k]) or not t.isascii():
            raise SystemExit(f"{k}: validation {errs} / undeclared {bad} / unchanged {set(PATHS[k]) - changed}")
        new_texts[k] = t
    for k, t in new_texts.items():
        (REG / f"{k}.yaml").write_text(t)

    # ---- manifest
    for f in ("h1_event_split_A2.tsv", "units_h1_components_A2.txt", "units_h1s_components_A2.txt",
              "lockbox_testable.tsv"):
        man["files"][f"results/part2/fzx/{f}"] = sha(FZX / f)
    for k in (H1, H1S):
        man["registry"][k] = dict(man["registry"][k], amendment="A1, A2", sha256=sha(REG / f"{k}.yaml"))
    s = man["summary"]

    def msum(e):
        return dict(m=len(units_of(e)), n_events=len(e), n_N_axis=int((e.n_N >= 2).sum()),
                    n_P_axis=int((e.n_P >= 1).sum()), complex_set_hash=registry.set_hash(units_of(e)))

    before = dict(H1=dict(s["m"]["H1"]), H1S=dict(s["m"]["H1S"]))
    s["m"]["H1"] = dict(msum(ev_p), sublabels={k: int(v) for k, v in ev_p.genetic_sublabel.value_counts().items()})
    s["m"]["H1S"] = dict(msum(ev_s), standing="exploratory", direction="H1S")
    s["lockbox_testable"] = {d: dict(testable=int((g.lockbox_testable.astype(str) == "True").sum()), total=int(len(g)))
                             for d, g in lb.groupby("entry")}
    man["complex_set_hashes"]["h1_components"] = registry.set_hash(u_p)
    man["complex_set_hashes"]["h1s_components"] = registry.set_hash(u_s)
    man["code"]["scripts/43_amend_A2.py"] = sha(ROOT / "scripts/43_amend_A2.py")
    for f in ("src/part2/registry.py", "src/part2/fzx.py", "src/part2/cscore.py"):
        man["code"][f] = sha(ROOT / f)
    man["amendments"].append(dict(
        id="A2", timestamp_utc=timestamp, amends_commit=c, record="docs/registry/AMENDMENTS.md",
        script="scripts/43_amend_A2.py",
        decision="event is genetic (H1) iff >= 50% of A loss models carry a LoF or hom-del call; else silencing (H1S)",
        moved=[dict(event=r.event, gene=r.gene, fzx_type=r.type, genetic_share=float(r.genetic_share),
                    sublabel=r.genetic_sublabel, from_entry=r.entry_A1, to_entry=r.entry)
               for r in moved.itertuples(index=False)],
        before=before, H1_after=s["m"]["H1"], H1S_after=s["m"]["H1S"],
        inputs_read="FZ-X per-gene loss-kind counts (h1_contexts.tsv) and event list only; no dependency values, "
                    "association statistics or lockbox outcomes"))
    man["files"] = dict(sorted(man["files"].items()))
    (FZX / "FZX_manifest.json").write_text(json.dumps(man, indent=1) + "\n")
    print(json.dumps(dict(moved=man["amendments"][-1]["moved"], H1=s["m"]["H1"], H1S=s["m"]["H1S"]), indent=1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--timestamp", default=TIMESTAMP)
    apply(ap.parse_args().timestamp)


if __name__ == "__main__":
    main()
