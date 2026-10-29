"""Member-resolved profiles for H1 / H1S (registry score_rule: "the member-resolved profile (beta_j with 95% CI for
every member)", an "also reported" quantity that the discovery runs did not export). Construction set A only.

Descriptive export: NOT a trial, writes no trials.jsonl line, changes no engine logic, result or registry entry.

It reuses the engine's own per-member fit `h1.member_profile` (the one `h1.classify` uses for the breadth criterion)
with the same events, loss column, covariates (`covariates.per_context`), complete-case rule, [-3, 3] clipping and
baselines (`h1.baselines` / `Components.b`) as scripts/44_run_direction.py. Evaluable axes are taken from the committed
`tests.tsv.gz` (the run_association table `classify` consulted): an axis is exported only when it exists in the engine
(N: |N'| >= 2; P: |P| >= 1) and its counted test is `evaluable`. No lockbox model, matrix or manifest is read: the
effect matrix comes from `assoc.load_A_effects` (A only, hash-checked in `real_env`).

Columns of results/part2/directions/<ID>/member_profiles.tsv (one row per event x axis x member):
  event, gene, type, link, axis (N = complex-wide non-paralog members, P = in-complex paralogs), member (DepMap column
  key), in_aligned_set (N: member of N'; P: |b_j| >= 0.1), n (complete-case models with the member observed),
  n_loss (of those, loss = 1), beta, se, df, ci95_lo, ci95_hi (HC3 OLS, Satterthwaite df, raw Chronos scale; beta < 0
  = member more dependent in loss lines), baseline (b_j: median over reference A models), aligned_beta (ALIGN = +1,
  E1: sign(b_j) * beta_j, > 0 = dependency grows in its baseline direction; blank when |b_j| < 0.1) with its CI.

Consistency assertion (before anything is written): the per-event breadth fractions recomputed from the export
(share of N' members with aligned_beta > 0 / < 0) equal N_frac_aligned_pos / N_frac_aligned_neg of the committed
classes.tsv, and the set of events with an exported N axis equals the set with a breadth fraction in classes.tsv.

Usage (project root, PART1_DIR set, PYTHONHASHSEED=0):
    .venv/bin/python scripts/55_export_member_profiles.py [--id H1_within_complex_GI | --id H1S_silencing_GI]
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import covariates as covmod, cscore, h1  # noqa: E402

IDS = ("H1_within_complex_GI", "H1S_silencing_GI")
FRAC_ATOL = 1.5e-6          # classes.tsv stores 6 significant digits


def _runner():
    spec = importlib.util.spec_from_file_location("run_direction_44", ROOT / "scripts" / "44_run_direction.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod                          # dataclasses in the runner need its module registered
    spec.loader.exec_module(mod)
    return mod


def axis_profile(E_A: pd.DataFrame, event: h1.Event, comp: h1.Components, axis: str, ctx: pd.Series,
                 cov: pd.DataFrame | None) -> pd.DataFrame:
    """Per-member rows of one evaluable event axis, from the engine's own `h1.member_profile` fit."""
    if axis == "N":
        genes, b = list(event.N), comp.b
        aligned = set(comp.Nprime)
    else:
        genes, b = list(event.P), h1.baselines(E_A, event.P, ctx)
        aligned = {g for g in genes if abs(b.get(g, 0.0)) >= h1.MIN_BASELINE}
    prof = h1.member_profile(E_A, genes, ctx, cov, b)
    X, ok = h1._xy(ctx, cov, E_A.index)
    Y = E_A.reindex(columns=genes).astype(float).clip(-cscore.CLIP, cscore.CLIP).values[ok]
    obs = ~np.isnan(Y)
    in_set = np.array([g in aligned for g in genes])
    sgn = np.sign(prof["baseline"].values)
    a_lo, a_hi = sgn * prof.ci95_lo.values, sgn * prof.ci95_hi.values
    return pd.DataFrame({
        "event": event.id, "gene": event.gene, "type": event.type, "link": event.link, "axis": axis,
        "member": genes, "in_aligned_set": in_set, "n": obs.sum(axis=0),
        "n_loss": (obs * X[:, 1][:, None]).sum(axis=0).astype(int),
        "beta": prof.beta.values, "se": prof.se.values, "df": prof.df.values,
        "ci95_lo": prof.ci95_lo.values, "ci95_hi": prof.ci95_hi.values,
        "baseline": prof.baseline.values,
        "aligned_beta": np.where(in_set, prof.aligned_beta.values, np.nan),
        "aligned_ci95_lo": np.where(in_set, np.minimum(a_lo, a_hi), np.nan),
        "aligned_ci95_hi": np.where(in_set, np.maximum(a_lo, a_hi), np.nan)})


def build_profiles(E_A: pd.DataFrame, events: list, loss_A: pd.DataFrame, covs: dict,
                   tests: pd.DataFrame) -> pd.DataFrame:
    """Profiles of every engine axis (N: |N'| >= 2, P: |P| >= 1) of every event on A, with `counted_evaluable`
    from the committed run_association table `tests`. The caller exports only the evaluable ones."""
    evaluable = dict(zip(tests.unit, tests.evaluable.astype(bool)))
    parts = []
    for e in events:
        c = h1.component_scores(E_A, e, loss_A[e.gene])
        ctx = loss_A[e.gene].rename(e.context)
        for axis in ("N", "P"):
            if axis in c.y:
                p = axis_profile(E_A, e, c, axis, ctx, covs.get(e.context))
                p.insert(5, "counted_evaluable", evaluable.get(f"{e.id}|{axis}", False))
                parts.append(p)
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def breadth_from_profiles(prof: pd.DataFrame) -> pd.DataFrame:
    """Per-event N-axis breadth recomputed from the export, exactly as h1.classify does over N' members."""
    n = prof[(prof.axis == "N") & prof.in_aligned_set]
    g = n.groupby("event", sort=False).aligned_beta
    return pd.DataFrame({"N_frac_aligned_pos": g.apply(lambda s: float((s > 0).mean())),
                         "N_frac_aligned_neg": g.apply(lambda s: float((s < 0).mean()))})


def check_breadth(full: pd.DataFrame, classes: pd.DataFrame) -> dict:
    """Assert export-implied breadth == committed classes.tsv breadth (N_frac_aligned_pos / _neg, tolerance = the 6
    significant digits stored). `full` holds every engine axis with `counted_evaluable`.
      (1) every event with a committed breadth is reproduced from the export, with equal fractions;
      (2) every event with an EVALUABLE exported N axis has a committed breadth (none missing).
    The committed H1 run predates the R2 engine fix, so classes.tsv also holds breadth for N axes whose counted test
    is not evaluable (e.g. MAP2K4: 11 loss models); those are reproduced under (1) but not exported. H1S ran with the
    fix: its classes.tsv has breadth for evaluable axes only, and the two sets coincide."""
    ref = classes.set_index("event")[["N_frac_aligned_pos", "N_frac_aligned_neg"]].dropna(how="all")
    allax = breadth_from_profiles(full)
    evax = breadth_from_profiles(full[full.counted_evaluable])
    miss = sorted(set(ref.index) - set(allax.index))
    assert not miss, f"committed breadth not reproducible from the export for {miss}"
    extra = sorted(set(evax.index) - set(ref.index))
    assert not extra, f"evaluable exported N axes without committed breadth: {extra}"
    mine = allax.loc[ref.index]
    for col in mine.columns:
        d = (mine[col] - ref[col]).abs()
        bad = pd.DataFrame({"export": mine[col], "classes": ref[col]})[d > FRAC_ATOL]
        assert (d <= FRAC_ATOL).all(), f"{col} mismatch:\n{bad}"
    return dict(committed_events=len(ref), evaluable_exported=len(evax),
                committed_not_evaluable=sorted(set(ref.index) - set(evax.index)))


def export_one(env, run, aid: str, loss_A: pd.DataFrame) -> dict:
    t0 = time.time()
    od = env.proj.results_dir / "directions" / aid
    classes = pd.read_csv(od / "classes.tsv", sep="\t")
    tests = pd.read_csv(od / "tests.tsv.gz", sep="\t")
    entry, _ = run.entry_of(env, aid)
    events = run.load_events(aid, env.proj.root)
    ctxs = sorted({e.context for e in events})
    covs, _ = covmod.per_context(entry, ctxs, env.cov)
    assert list(classes.event) == [e.id for e in events], "classes.tsv events differ from the frozen event list"
    full = build_profiles(env.E_A, events, loss_A, covs, tests)
    chk = check_breadth(full, classes)                   # asserts before anything is written
    prof = full[full.counted_evaluable].drop(columns="counted_evaluable").reset_index(drop=True)
    path = od / "member_profiles.tsv"
    prof.to_csv(path, sep="\t", index=False, float_format="%.6g")
    dropped = full[~full.counted_evaluable].groupby("axis").event.nunique().to_dict()
    return dict(id=aid, rows=len(prof), events=prof.event.nunique(),
                axes=prof.groupby("axis").event.nunique().to_dict(),
                breadth_check=chk,
                axes_not_evaluable_not_exported=dropped, path=str(path), seconds=round(time.time() - t0, 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--id", choices=IDS)
    a = ap.parse_args()
    run = _runner()
    env = run.real_env()
    loss_A = run.load_loss_A(env.proj.root, env.E_A.index)
    for aid in ([a.id] if a.id else IDS):
        print(export_one(env, run, aid, loss_A))


if __name__ == "__main__":
    main()
