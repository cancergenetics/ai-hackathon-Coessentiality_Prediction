#!/usr/bin/env python
"""Registered per-hit robustness (LOSO, leave-top-2-out, breadth) and the final strongest-findings list (DESIGN v3.1
§2.4-2.5; registry `success_criterion` / `multiverse.robustness` of H3a_lineage_atlas, H3b_driver_family,
H1_within_complex_GI; rule and LOSO units in src/part2/strongest.py).

Computed AFTER the lockbox read (the rule was registered before it; disclosed in results/part2/robustness/README.md).
A side only: the outcome matrix comes from assoc.load_A_effects (A models); no lockbox read, persisted read or loader
is touched. L summary statistics shown next to each hit are copied from the committed confirmation tables.
Not a trial: no discovery, confirmation, registry entry or ledger is changed.

Modes:
  --candidates   committed files only (no outcome data): writes results/part2/robustness/candidates.tsv, the hit
                 list named by INTENT.md.
  --run          refuses unless INTENT.md and candidates.tsv are committed and unchanged; recomputes the candidate
                 list and checks it equals the committed one; re-fits every hit on A with the discovery run's score,
                 removals and covariates (scripts/44_run_direction.py real_env; full-fit beta / t must reproduce the
                 committed discovery table), then LOSO, leave-top-2-out and breadth; applies the rule mechanically.
Outputs (results/part2/robustness/): robustness_per_hit.tsv, strongest_findings.tsv, summary.json.

Usage (PART1_DIR set, PYTHONHASHSEED=0):
    .venv/bin/python scripts/56_strongest_findings.py --candidates
    .venv/bin/python scripts/56_strongest_findings.py --run
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from part2 import covariates as covmod, cscore, h1, registry, robust, strongest  # noqa: E402

OUT = ROOT / "results" / "part2" / "robustness"
DIR = ROOT / "results" / "part2" / "directions"
LB = ROOT / "results" / "part2" / "lockbox_confirm"
H3 = {"H3a": "H3a_lineage_atlas", "H3b": "H3b_driver_family"}
H1_ID = "H1_within_complex_GI"
TOL = 1e-6


def _tsv(p: Path) -> pd.DataFrame:
    return pd.read_csv(p, sep="\t")


def _load_script(name: str, mod: str):
    if mod in sys.modules:
        return sys.modules[mod]
    spec = importlib.util.spec_from_file_location(mod, ROOT / "scripts" / name)
    m = importlib.util.module_from_spec(spec)
    sys.modules[mod] = m
    spec.loader.exec_module(m)
    return m


def membership() -> pd.DataFrame:
    return _tsv(ROOT / "results" / "part2" / "fzx" / "zavg_engine_membership.tsv")


def candidates() -> pd.DataFrame:
    """Every H3a / H3b discovery and every H1 eligible call with the rule's non-robustness criteria (committed files
    only). `run` = robustness is computed for the row (all candidates, plus the H1 eligible calls)."""
    calls = _tsv(LB / "primary_relabelled_R7.tsv.gz")
    P = _load_script("part2_report.py", "part2_report")
    rows = []
    for short, aid in H3.items():
        d = _tsv(DIR / aid / "discoveries.tsv")
        x = strongest.h3_candidates(d, calls[calls.entry == aid])
        spec = P.discoveries(aid)[["context", "unit", "beats_all", "control"]]
        x = x.merge(spec, on=["context", "unit"], how="left", validate="one_to_one")
        x["writer_specific_set"] = (x.primary_call.eq("confirmed") & x.beats_all.astype(bool) & ~x.locus
                                    & ~x.defect_sensitive & x.control.ne("control readout"))
        x.insert(0, "direction", short)
        x["run"] = x.candidate
        rows.append(x)
    h1c = _tsv(LB / "h1_events_primary_vs_corrected_R7.tsv")
    conf = _tsv(LB / H1_ID / "confirmation.tsv")
    x = strongest.h1_candidates(_tsv(DIR / H1_ID / "discoveries.tsv"), h1c, _tsv(ROOT / "results/part2/fzx/h1_events.tsv"),
                                conf)
    x.insert(0, "direction", "H1")
    x["run"] = True                      # every H1 eligible call (task scope), candidate or not
    rows.append(x)
    cols = ["direction", "context", "unit", "beta", "t", "q_final", "tier", "locus", "primary_call", "route",
            "defect_sensitive", "discovery_q_lt_q", "lockbox_pass", "candidate", "run", "beats_all", "control",
            "writer_specific_set", "headline_eligible", "label_final"]
    out = pd.concat([r.reindex(columns=cols) for r in rows], ignore_index=True)
    return out.sort_values(["direction", "context", "unit"], kind="mergesort").reset_index(drop=True)


def write_candidates():
    c = candidates()
    OUT.mkdir(parents=True, exist_ok=True)
    c.to_csv(OUT / "candidates.tsv", sep="\t", index=False, float_format="%.6g")
    s = {d: dict(rows=int((c.direction == d).sum()), candidates=int(c[c.direction == d].candidate.sum()),
                 run=int(c[c.direction == d].run.astype(bool).sum()),
                 writer_specific_set=int(c[c.direction == d].writer_specific_set.astype(str).eq("True").sum()))
         for d in ("H3a", "H3b", "H1")}
    print(json.dumps(s, indent=1))


# ------------------------------------------------------------------------------------------------ run
def _check_committed():
    for p in (OUT / "INTENT.md", OUT / "candidates.tsv"):
        registry.committed_blob(ROOT, p)
    fresh = candidates()
    fresh.to_csv(buf := Path(OUT / ".candidates_check.tsv"), sep="\t", index=False, float_format="%.6g")
    same = buf.read_bytes() == (OUT / "candidates.tsv").read_bytes()
    buf.unlink()
    if not same:
        raise RuntimeError("recomputed candidate list differs from the committed candidates.tsv")
    return fresh


def _members_h3(fm, unit: str, removal, annot) -> list:
    mem = fm.members[unit]
    rem = set(removal.removed(mem, annot))
    return [g for g in mem if g not in rem]


def run():
    t0 = time.time()
    cand = _check_committed()
    todo = cand[cand.run.astype(bool)]
    D44 = _load_script("44_run_direction.py", "run_direction44")
    env = D44.real_env()                                         # A models only (hash-checked), frozen map / contexts
    A = env.E_A.index
    rows = []
    timing = {}
    for short, aid in list(H3.items()) + [("H1", H1_ID)]:
        t1 = time.time()
        hits = todo[todo.direction == short]
        entry, _ = registry.load_entry(env.proj.registry_dir, aid, root=env.proj.root)
        disc = D44.committed_discovery(env, aid).set_index(["context", "unit"])
        groups = env.cov[strongest.LOSO_UNIT[short]].reindex(A)
        if short == "H1":
            events = D44.load_events(aid, env.proj.root)
            loss_A = D44.load_loss_A(env.proj.root, A)
            Y1, ctx_all, _, comps = h1.assoc_inputs(env.E_A, events, loss_A)
            covs, _ = covmod.per_context(entry, list(ctx_all.columns), env.cov)
            ev_of = {e.id: e for e in events}
            score = lambda c, u: Y1[c][u]                                                  # noqa: E731
            members = lambda c, u: (ev_of[u.rsplit("|", 1)[0]].P if u.endswith("|P")      # noqa: E731
                                    else list(comps[u.rsplit("|", 1)[0]].Nprime))
        else:
            ctx_all = env.ctx
            ctxs = sorted(set(hits.context))
            covs, _ = covmod.per_context(entry, ctxs, env.cov)
            cache = {}

            def score(c, u):
                rem = cscore.removal_for_context(c, env.homdel.get(c))
                if rem not in cache:
                    cache[rem] = cscore.unit_scores(env.E_A, env.fmap, env.fmap.coarse(), rem, env.annot)[0]
                return cache[rem][u]

            members = lambda c, u: _members_h3(env.fmap, u, cscore.removal_for_context(c, env.homdel.get(c)),  # noqa
                                               env.annot)
        for r in hits.itertuples(index=False):
            d = disc.loc[(r.context, r.unit)]
            alt = d["alternative"]
            y, c, Z = score(r.context, r.unit), ctx_all[r.context].reindex(A), covs[r.context]
            out = strongest.hit_robustness(y, c, Z, groups, alt)
            if not (abs(out["beta_A_refit"] - float(d["beta"])) <= TOL * max(1, abs(float(d["beta"])))
                    and abs(out["t_A_refit"] - float(d["t"])) <= 1e-4 * max(1, abs(float(d["t"])))):
                raise RuntimeError(f"{short} {r.context} {r.unit}: A refit {out['beta_A_refit']}/{out['t_A_refit']} "
                                   f"does not reproduce the discovery {d['beta']}/{d['t']}")
            mem = members(r.context, r.unit)
            if short != "H1":            # H1 breadth is the engine's (classes.tsv; aligned N axis), not re-derived
                out.update(robust.breadth(env.E_A, mem, c, Z, np.sign(out["beta_A_refit"])))
            out.update(direction=short, context=r.context, unit=r.unit, alternative=alt,
                       loso_unit=strongest.LOSO_UNIT[short], members=";".join(sorted(mem)))
            rows.append(out)
        timing[short] = round(time.time() - t1, 1)
        print(f"{short}: {len(hits)} hits in {timing[short]} s", flush=True)
    rob = pd.DataFrame(rows)
    x = cand.merge(rob, on=["direction", "context", "unit"], how="left", validate="one_to_one")
    x = x[x.run.astype(bool)]
    x = strongest.apply_rule(x)
    x = x.merge(_l_stats(), on=["direction", "context", "unit"], how="left", validate="one_to_one")
    _write(x, timing, time.time() - t0)


def _l_stats() -> pd.DataFrame:
    """L summary statistics as committed by the single read (no lockbox data is opened)."""
    out = []
    for short, aid in list(H3.items()) + [("H1", H1_ID)]:
        c = _tsv(LB / aid / "confirmation.tsv")[["context", "unit", "beta_L", "ci95_L_lo", "ci95_L_hi", "p_L",
                                                 "n_context_L", "loso_same_sign", "loso_folds", "lomo_same_sign",
                                                 "demeter2_same_sign"]]
        c = c.rename(columns={"loso_same_sign": "rare_loso_same_sign", "loso_folds": "rare_loso_folds",
                              "lomo_same_sign": "rare_lomo_same_sign", "demeter2_same_sign": "rare_demeter2_same_sign"})
        c.insert(0, "direction", short)
        out.append(c)
    return pd.concat(out, ignore_index=True)


def _write(x: pd.DataFrame, timing: dict, wall: float):
    OUT.mkdir(parents=True, exist_ok=True)
    per = ["direction", "context", "unit", "members", "tier", "beta", "t", "q_final", "alternative", "primary_call",
           "route", "defect_sensitive", "locus", "candidate", "loso_unit", "loso_folds", "loso_pass", "loso_skipped",
           "loso_robust", "top2_beta", "top2_p", "top2_same_sign", "top2_robust", "breadth", "breadth_n",
           "n_A", "n_context_A", "strongest", "rank", "writer_specific_set", "beats_all", "control",
           "headline_eligible", "label_final"]
    x.sort_values(["direction", "context", "unit"], kind="mergesort")[per].to_csv(
        OUT / "robustness_per_hit.tsv", sep="\t", index=False, float_format="%.6g")
    s = x[x.strongest].sort_values("rank")
    cols = ["rank", "direction", "context", "unit", "members", "tier", "tier_class", "beta", "t", "q_final", "route",
            "beta_L", "ci95_L_lo", "ci95_L_hi", "p_L", "n_context_L", "loso_unit", "loso_folds", "loso_pass",
            "loso_skipped", "top2_beta", "top2_p", "top2_robust", "breadth", "breadth_n", "writer_specific_set",
            "control", "label_final"]
    s[cols].rename(columns={"beta": "beta_A", "t": "t_A"}).to_csv(OUT / "strongest_findings.tsv", sep="\t",
                                                                    index=False, float_format="%.6g")
    summ = {"rule": "discovery q < 0.05 AND primary lockbox-confirmed (or rare-context route passed) AND not "
                    "defect-sensitive (R7) AND not locus AND LOSO-robust (every evaluable fold same sign, p < 0.05); "
                    "H1 also non-neutral and headline-eligible (E2); ranked by tier class then |t_A|",
            "loso_unit": strongest.LOSO_UNIT, "per_direction": {}, "runtime_s": round(wall, 1), "timing_s": timing,
            "peak_rss_bytes": int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)}
    for d in ("H3a", "H3b", "H1"):
        g = x[x.direction == d]
        cg = g[g.candidate.astype(bool)]
        summ["per_direction"][d] = dict(
            computed=int(len(g)), candidates=int(len(cg)), loso_robust=int(cg.loso_robust.sum()),
            top2_robust=int(cg.top2_robust.astype(bool).sum()),
            loso_and_top2=int((cg.loso_robust & cg.top2_robust.astype(bool)).sum()),
            strongest=int(g.strongest.sum()), strongest_units=int(g[g.strongest].unit.nunique()),
            loso_fail_any_fold=int((cg.loso_pass < cg.loso_folds).sum()),
            loso_no_evaluable_fold=int((cg.loso_folds == 0).sum()),
            writer_specific_set=int(g.writer_specific_set.astype(str).eq("True").sum()),
            writer_specific_set_strongest=int((g.writer_specific_set.astype(str).eq("True") & g.strongest).sum()),
            strongest_per_tier={k: int(v) for k, v in g[g.strongest].tier.value_counts().sort_index().items()})
    (OUT / "summary.json").write_text(json.dumps(summ, indent=1, sort_keys=True, default=str) + "\n")
    print(json.dumps(summ, indent=1, default=str))


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--candidates", action="store_true")
    g.add_argument("--run", action="store_true")
    a = ap.parse_args()
    write_candidates() if a.candidates else run()


if __name__ == "__main__":
    main()
