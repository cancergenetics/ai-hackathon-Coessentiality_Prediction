"""POST HOC, DESCRIPTIVE complex-specificity counts (not part of the registered error control; no trial line).

For every primary discovery of H3a, H3b and H2: does its effect lie outside the central 95% (numpy linear 2.5% and
97.5% quantiles) of the effects of its own matched random sets in the same context?
  H3a / H3b   the registered matched random gene sets (R = 10 per coarse unit; calibration D8, scripts/46)
  H2          the registered random member splits of the same parent (R = 20 per tested (context, sub-module))

Two versions:
  t  from the COMMITTED calibration outputs only (<ID>/calibration_D8_tests.tsv.gz stores t and p, not beta);
  beta  the committed calibration does not store beta, so the same sets are regenerated deterministically with the
        committed code and seeds (calib.matched_random_sets / calib.random_split_scores, sub-streams as in
        scripts/46_run_calibrations.py) and re-fitted; the regenerated t must equal the committed t (max |diff| is
        reported; the run stops if it exceeds 1e-3 relative), so the beta quantiles come from exactly the committed sets.
A-side data only (load_A_effects); no lockbox read; nothing registered is changed.

Writes results/part2/directions/specificity.json and specificity_discoveries.tsv.gz.
Usage: .venv/bin/python scripts/49_specificity_counts.py
"""
from __future__ import annotations

import gzip
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

from part2 import calib, covariates as covmod, cscore, h2  # noqa: E402

spec = importlib.util.spec_from_file_location("run_calibrations", ROOT / "scripts" / "46_run_calibrations.py")
C = importlib.util.module_from_spec(spec)
sys.modules["run_calibrations"] = C
spec.loader.exec_module(C)
R = C.R

LABEL = "post hoc descriptive specificity filter (not registered error control)"
DIRS = ROOT / "results" / "part2" / "directions"
CAL = {"H3a_lineage_atlas": "matched_random_set", "H3b_driver_family": "matched_random_set",
       "H2_submodule_divergence": "random_split"}


def committed_null(aid: str) -> pd.DataFrame:
    with gzip.open(DIRS / aid / "calibration_D8_tests.tsv.gz", "rt") as fh:
        t = pd.read_csv(fh, sep="\t")
    return t[t.calibration == CAL[aid]]


def outside(x: float, null: np.ndarray) -> tuple:
    v = np.asarray(null, float)
    v = v[np.isfinite(v)]
    if len(v) < 2 or not np.isfinite(x):
        return np.nan, np.nan, len(v), None
    lo, hi = np.quantile(v, [0.025, 0.975])
    return lo, hi, len(v), bool(x < lo or x > hi)


def regen_h3(env, aid: str) -> pd.DataFrame:
    """beta, t of every regenerated matched set (context, set, unit)."""
    entry, _ = R.entry_of(env, aid)
    seed = int(entry["calibration"]["seed"])
    ctxs = R.primary_contexts(entry)
    units = env.fmap.coarse()
    covs, _ = covmod.per_context(entry, ctxs, env.cov)
    sets, _ = calib.matched_random_sets(env.E_A, env.fmap, units, calib.clip_essentiality(env.E_A),
                                        R.load_pair_scores(env.proj.root), C.R_SETS, seed)
    of = {k: k.split("#")[0] for k in sets}
    rows, cache = [], {}
    for c in ctxs:
        rem = cscore.removal_for_context(c, env.homdel.get(c))
        if rem not in cache:
            cache[rem] = calib.score_sets(env.E_A, sets, rem, env.annot)
        f = calib.test_columns(cache[rem], env.ctx[c], covs[c], C.alt_of(entry, c))
        rows.append(pd.DataFrame(dict(context=c, set=f.index, unit=[of[k] for k in f.index], beta=f.beta.values,
                                      t=f.t.values)))
    return pd.concat(rows, ignore_index=True)


def regen_h2(env, aid: str, disc: pd.DataFrame) -> pd.DataFrame:
    entry, _ = R.entry_of(env, aid)
    seed = int(entry["calibration"]["seed"])
    ctxs = R.primary_contexts(entry)
    subs = h2.h2_units(env.fmap)
    cov2, _ = covmod.per_context(entry, ctxs, env.cov)
    rows = []
    for r in disc.itertuples(index=False):
        ci, si = ctxs.index(r.context), subs.index(r.unit)
        rem = cscore.removal_for_context(r.context, env.homdel.get(r.context))
        Ys = calib.random_split_scores(env.E_A, env.fmap, r.unit, rem, env.annot, C.R_SPLITS,
                                       np.random.default_rng([seed, si, ci]))
        f = calib.test_columns(Ys, env.ctx[r.context], cov2[r.context], r.alternative)
        rows.append(pd.DataFrame(dict(context=r.context, set=f.index, unit=r.unit, beta=f.beta.values, t=f.t.values)))
    return pd.concat(rows, ignore_index=True)


def main():
    t0 = time.time()
    env = R.real_env()
    out, per = dict(label=LABEL, quantiles="numpy linear 2.5% / 97.5% of the finite matched-set values; outside = "
                                           "strictly below the 2.5% or above the 97.5% quantile"), []
    for aid in CAL:
        real = R.committed_discovery(env, aid)
        disc = real[real.discovery].copy()
        cm = committed_null(aid)
        rg = regen_h2(env, aid, disc) if aid.startswith("H2") else regen_h3(env, aid)
        chk = cm.merge(rg, on=["context", "set", "unit"], how="inner", suffixes=("_committed", "_regen"))
        fin = np.isfinite(chk.t_committed) & np.isfinite(chk.t_regen)
        rel = (np.abs(chk.t_committed - chk.t_regen) / np.maximum(1.0, np.abs(chk.t_committed)))[fin]
        nan_mismatch = int((np.isfinite(chk.t_committed) != np.isfinite(chk.t_regen)).sum())
        if aid.startswith("H2"):
            need = cm.merge(disc[["context", "unit"]], on=["context", "unit"])
        else:
            need = cm
        if len(chk) != len(need) or nan_mismatch or (len(rel) and rel.max() > 1e-3):
            raise RuntimeError(f"{aid}: regenerated sets differ from the committed calibration "
                               f"({len(chk)} vs {len(need)} rows, {nan_mismatch} NaN mismatches, max rel {rel.max()})")
        gt = {k: g.t.values for k, g in cm.groupby(["context", "unit"])}
        gb = {k: g.beta.values for k, g in rg.groupby(["context", "unit"])}
        for r in disc.itertuples(index=False):
            k = (r.context, r.unit)
            tlo, thi, nt, ot = outside(r.t, gt.get(k, []))
            blo, bhi, nb, ob = outside(r.beta, gb.get(k, []))
            per.append(dict(analysis_id=aid, context=r.context, unit=r.unit, beta=r.beta, t=r.t, n_sets=nb,
                            beta_q025=blo, beta_q975=bhi, outside_beta=ob, t_q025=tlo, t_q975=thi, outside_t=ot,
                            same_side_beta=None if ob is None else bool((r.beta > bhi) if r.beta > 0 else (r.beta < blo))))
        d = pd.DataFrame([p for p in per if p["analysis_id"] == aid])
        has = d.outside_beta.notna()
        ctx_tab = d.groupby("context").agg(discoveries=("unit", "size"), with_sets=("outside_beta", lambda v: int(v.notna().sum())),
                                           outside_beta=("outside_beta", lambda v: int((v == True).sum())),  # noqa: E712
                                           outside_t=("outside_t", lambda v: int((v == True).sum())))  # noqa: E712
        out[aid] = dict(null=CAL[aid], n_discoveries=int(len(d)), n_with_sets=int(has.sum()),
                        n_without_sets=int((~has).sum()), outside_beta=int((d.outside_beta == True).sum()),  # noqa: E712
                        outside_beta_same_side=int((d.same_side_beta == True).sum()),  # noqa: E712
                        outside_t=int((d.outside_t == True).sum()),  # noqa: E712
                        regenerated_check=dict(rows_compared=int(len(chk)), max_rel_t_diff=float(rel.max()) if len(rel) else None),
                        per_context=ctx_tab.astype(int).to_dict("index"))
        print(aid, {k: v for k, v in out[aid].items() if k != "per_context"})
    out["runtime_s"] = round(time.time() - t0, 1)
    out["peak_rss_bytes"] = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    R.tsv_gz(pd.DataFrame(per), DIRS / "specificity_discoveries.tsv.gz")
    (DIRS / "specificity.json").write_text(json.dumps(out, indent=1, sort_keys=True, default=float) + "\n")
    print(f"wall {out['runtime_s']} s; peak RSS {out['peak_rss_bytes'] / 1e9:.2f} GB")


if __name__ == "__main__":
    main()
