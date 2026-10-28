"""Pre-registered controls (DESIGN v3.1 §2.7): ONE run after FZ-C and FZ-X; writes results/part2/controls/
(controls_summary.tsv, kill_rule.json) and refuses to run a second time. Commit the outputs (and trials.jsonl)
before any direction's discovery run.

Usage (after every CTRL_* entry is filled and committed):
    .venv/bin/python scripts/34_run_controls.py --map results/part2/fzx/zavg_engine_membership.tsv \
        --context-table data_cache/part2/contexts/context_table.parquet

Both files must be listed with their SHA-256 in the committed freeze manifests named in the entries'
frozen_inputs (run_association checks this): the engine-format map in results/part2/fzx/FZX_manifest.json (FZ-X),
the context table in results/part2/contexts/checksums.json. Dependency data: A models only, via assoc.load_A_effects.
NOT RUN by the engine author (no dependency data may be touched before FZ-X).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from part2 import assoc, controls, covariates as covmod, cscore  # noqa: E402
from part2.cmap import FrozenMap  # noqa: E402


def homdel_shares(ctx: pd.DataFrame, A) -> dict:
    """LoF context -> fraction of its positive A lines that are hom-del calls of the gene (CN window rule)."""
    out = {}
    for c in ctx.columns:
        if c.startswith(("LoF:", "LoF_inclMSI:")):
            g = c.split(":", 1)[1]
            hd = f"CN_homdel:{g}"
            pos = ctx.loc[A, c] == 1
            if hd in ctx.columns and pos.any():
                out[c] = float((ctx.loc[A, hd][pos] == 1).mean())
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--map", required=True)
    ap.add_argument("--context-table", required=True)
    a = ap.parse_args()
    proj = assoc.Project()
    fmap = FrozenMap.from_membership(pd.read_csv(a.map, sep="\t"))
    ctx = pd.read_parquet(a.context_table)
    E_A = assoc.load_A_effects(proj=proj)
    A = E_A.index
    cov = covmod.source_table(ctx, A)
    annot = cscore.load_annotation(E_A.columns)
    tab, kill = controls.evaluate_controls(E_A, fmap, ctx.reindex(A), cov, annot, inputs=[a.context_table, a.map],
                                           proj=proj, homdel_share=homdel_shares(ctx, A))
    print(json.dumps(kill, indent=1))


if __name__ == "__main__":
    main()
