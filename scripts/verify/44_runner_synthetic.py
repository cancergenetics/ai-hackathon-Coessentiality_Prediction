"""Synthetic check of scripts/44_run_direction.py (no DepMap value is read).

Builds the throw-away synthetic project of tests/part2_synth_pipeline.py (part2.synth.simulate, seed 0), runs the
runner's H3a -> H3b -> (commit) -> H2 -> H1 -> SMARCB1-style CN sensitivity (on the synthetic ARID1A event) with the
registered 10,000 permutations and K = 1000, and checks:
  * the runner's H3a / H3b discoveries equal those of the reference pipeline (tests/part2_synth_pipeline.discovery)
    when the calibration response does not fire (the runner adds no analysis logic);
  * FDR / power metrics and H1 class accuracy against the planted truth;
  * the single-look guard refuses a second run; outputs and trial lines are written.
Usage: .venv/bin/python scripts/verify/44_runner_synthetic.py [out.json]
"""
from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))

import part2_synth_pipeline as P  # noqa: E402
from part2 import assoc, registry, synth  # noqa: E402

spec = importlib.util.spec_from_file_location("run_direction", ROOT / "scripts" / "44_run_direction.py")
R = importlib.util.module_from_spec(spec)
sys.modules["run_direction"] = R
spec.loader.exec_module(R)


def env_for(sim, rep):
    proj = rep["proj"]
    tab = pd.read_parquet(proj.root / P.REL_CTX)
    E_A = assoc.load_A_effects(proj=proj)
    A = E_A.index
    return R.Env(proj=proj, fmap=sim.fmap, ctx=tab.loc[A], E_A=E_A, cov=tab.loc[A], annot=sim.annot,
                 inputs=rep["inputs"], blocks=tab.PatientID.loc[A], strata=tab.OncotreeLineage.loc[A],
                 strict=False, write=True)


def main():
    out = {}
    t0 = time.time()
    sim = synth.simulate(0)
    loss, kind, events = P.h1_inputs(sim)
    with tempfile.TemporaryDirectory() as td:
        rep = P.build_repo(Path(td) / "runner", sim, events)
        env = env_for(sim, rep)
        A = env.E_A.index
        t = time.time()
        r3 = {aid: R.run_h3(env, aid) for aid in ("H3a_lineage_atlas", "H3b_driver_family")}
        out["H3_runtime_s"] = round(time.time() - t, 1)
        for aid, key in (("H3a_lineage_atlas", "H3a"), ("H3b_driver_family", "H3b")):
            out[key] = dict(metrics=P._h3_metrics(r3[aid]["res"], sim, key),
                            calibration={k: r3[aid]["calibration"][k] for k in
                                         ("n_pseudo", "n_tests", "frac_p_lt_005", "response_fires", "bb_rejections")},
                            m=r3[aid]["summary"]["m_evaluable"], screened=r3[aid]["summary"]["n_screened_perm"])
        # reference pipeline on a fresh copy of the same synthetic project
        rep2 = P.build_repo(Path(td) / "reference", sim, events)
        _, _, ref = P.discovery(sim, rep2["proj"], rep2["inputs"], n_perm=R.N_PERM)
        for aid in r3:
            a = set(map(tuple, r3[aid]["res"].loc[r3[aid]["res"].discovery, ["context", "unit"]].values))
            b = set(map(tuple, ref[aid]["res"].loc[ref[aid]["res"].discovery, ["context", "unit"]].values))
            fires = r3[aid]["calibration"]["response_fires"]
            out[f"{aid}_equals_reference"] = (a == b) if not fires else f"response fired; |runner|={len(a)} |ref|={len(b)}"
        P.synth_commit(env.proj.root, "H3a/H3b discovery runs committed")
        t = time.time()
        r2 = R.run_h2(env)
        out["H2"] = dict(metrics=P._metrics(r2["res"], sim.truth["H2"], sim.truth.get("H2_expected_sibling", set()),
                                            planted=sim.truth.get("H2_planted", set())),
                         n_eligible_families=r2["summary"]["n_eligible_families"], m=r2["summary"]["m_evaluable"],
                         runtime_s=round(time.time() - t, 1))
        t = time.time()
        Ec = env.E_A.clip(-3, 3)
        ps = 10 * Ec[sorted(sim.fmap.all_genes(1))].corr()
        r1 = R.run_h1(env, "H1_within_complex_GI", events, loss.loc[A], ps)
        correct = n_axes = 0
        for r in r1["classes"].itertuples(index=False):
            exp = sim.truth["H1"].get(r.event, {})
            for ax, col in (("N", "class_complex_wide"), ("P", "class_paralog")):
                if ax in exp:
                    n_axes += 1
                    correct += getattr(r, col) == exp[ax]
        out["H1"] = dict(accuracy=correct / max(n_axes, 1), n_axes=n_axes, m=r1["summary"]["m_evaluable"],
                         class_counts=r1["summary"]["class_counts"], runtime_s=round(time.time() - t, 1))
        P.synth_commit(env.proj.root, "H2 / H1 discovery runs committed")
        sens = R.run_sensitivity_cn(env, "H1_within_complex_GI", "ARID1A", events, loss.loc[A], ps, label="ARID1A_CN")
        out["sensitivity"] = sens[["event", "covariate_set", "n_complete", "class_paralog", "beta_P"]].to_dict("records")
        prim = r1["classes"].set_index("event")
        s0 = sens[sens.covariate_set.str.startswith("primary")].set_index("event")
        out["sensitivity_primary_reproduces_classes"] = bool(
            all(abs(s0.loc[e, "beta_P"] - prim.loc[e, "beta_P"]) < 1e-12 and
                s0.loc[e, "class_paralog"] == prim.loc[e, "class_paralog"] for e in s0.index))
        env.strict = True
        refused = {}
        for aid in ("H3a_lineage_atlas", "H2_submodule_divergence", "H1_within_complex_GI"):
            try:
                R.guard_single_look(env, aid)
                refused[aid] = False
            except registry.RegistryError:
                refused[aid] = True
        out["second_run_refused"] = refused
        d = env.proj.results_dir / "directions"
        out["files"] = sorted(str(p.relative_to(d)) for p in d.rglob("*") if p.is_file())
        out["trial_lines"] = [json.loads(l)["analysis_id"] for l in env.proj.trials.read_text().splitlines()]
    out["wall_s"] = round(time.time() - t0, 1)
    txt = json.dumps(out, indent=1, default=str)
    print(txt)
    if len(sys.argv) > 1:
        Path(sys.argv[1]).write_text(txt + "\n")


if __name__ == "__main__":
    main()
