"""Before/after check for the h1 'not evaluable' fix (PI ruling R2): the synthetic e2e pipeline
(tests/part2_synth_pipeline.run_pipeline, seed 0, 10,000 permutations, K = 1000) and the direction runner on the same
synthetic project must give identical outputs, compared by SHA-256 of every output file (volatile fields such as
runtimes, commit heads and timestamps removed). No DepMap value is read.

Usage: .venv/bin/python scripts/verify/45_h1_fix_equivalence.py <out.json>
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import tempfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))

import part2_synth_pipeline as P  # noqa: E402
from part2 import assoc, synth  # noqa: E402

spec = importlib.util.spec_from_file_location("run_direction", ROOT / "scripts" / "44_run_direction.py")
R = importlib.util.module_from_spec(spec)
sys.modules["run_direction"] = R
spec.loader.exec_module(R)

VOLATILE = {"runtime_s", "peak_rss_bytes", "classify_runtime_s", "head", "time", "wall_s", "trial_line",
            "full_table", "full_table_sha256"}


def _clean(x):
    if isinstance(x, dict):
        return {k: _clean(v) for k, v in x.items() if k not in VOLATILE}
    if isinstance(x, list):
        return [_clean(v) for v in x]
    return x


def _sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def main():
    out = {}
    sim = synth.simulate(0)
    with tempfile.TemporaryDirectory() as td:
        e2e = P.run_pipeline(sim, Path(td) / "e2e")
        out["e2e_pipeline_json"] = _sha(json.dumps(_clean(e2e), sort_keys=True, default=str).encode())
        loss, kind, events = P.h1_inputs(sim)
        rep = P.build_repo(Path(td) / "runner", sim, events)
        proj = rep["proj"]
        tab = pd.read_parquet(proj.root / P.REL_CTX)
        E_A = assoc.load_A_effects(proj=proj)
        A = E_A.index
        env = R.Env(proj=proj, fmap=sim.fmap, ctx=tab.loc[A], E_A=E_A, cov=tab.loc[A], annot=sim.annot,
                    inputs=rep["inputs"], blocks=tab.PatientID.loc[A], strata=tab.OncotreeLineage.loc[A],
                    strict=False, write=True)
        for aid in ("H3a_lineage_atlas", "H3b_driver_family"):
            R.run_h3(env, aid)
        P.synth_commit(proj.root, "H3 committed")
        R.run_h2(env)
        ps = 10 * E_A.clip(-3, 3)[sorted(sim.fmap.all_genes(1))].corr()
        R.run_h1(env, "H1_within_complex_GI", events, loss.loc[A], ps)
        P.synth_commit(proj.root, "H1/H2 committed")
        R.run_sensitivity_cn(env, "H1_within_complex_GI", "ARID1A", events, loss.loc[A], ps, label="ARID1A_CN")
        d = proj.results_dir / "directions"
        files = {}
        for p in sorted(d.rglob("*")):
            if not p.is_file():
                continue
            rel = str(p.relative_to(d))
            if p.suffix == ".json":
                files[rel] = _sha(json.dumps(_clean(json.loads(p.read_text())), sort_keys=True).encode())
            else:
                files[rel] = _sha(p.read_bytes())
        out["runner_files"] = files
        cls = pd.read_csv(d / "H1_within_complex_GI" / "classes.tsv", sep="\t")
        out["h1_classes"] = cls[["event", "class_complex_wide", "class_paralog"]].values.tolist()
    txt = json.dumps(out, indent=1, sort_keys=True)
    print(txt)
    Path(sys.argv[1]).write_text(txt + "\n")


if __name__ == "__main__":
    main()
