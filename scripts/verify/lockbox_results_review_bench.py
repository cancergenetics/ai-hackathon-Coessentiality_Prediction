#!/usr/bin/env python
"""Independent BENCH recompute (results/part2/verify/lockbox_results_review.md, check 5).

The BENCH L matrix is obtained only through assoc.resume_lockbox_read(relook="R6") (the driver's R6/R7 resume path;
loader guarded; no ledger line). Delta(C, c) = signed HC3 t of S_C minus that of g* (clipped at +/-3), my own HC3 (FWL)
fits; PatientID-blocked bootstrap over L models with the registered seed (5100) and the same draw order as the
driver, so the replicates are the executed ones. Two designs per replicate:
  executed : covariates one-hot assembled over all L rows, then subset to the replicate rows (what bench_delta does);
  assembled: one-hot assembled over the replicate rows (the engine convention applied inside the resample).
Outputs: <out>/bench_recompute.json, bench_pairs_recompute.tsv (per-pair point Delta under the executed design).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lockbox_results_review as R  # noqa: E402
from part2 import assoc, covariates as covmod, cscore, data, paths  # noqa: E402

OUTC = R.OUTC
BENCH = "BENCH_complex_vs_gene"


def hc3_t_many(Y, x, Z):
    """Two-sided HC3 t of x for every column of Y (n x k, complete); NaN when the design cannot be fitted."""
    n = len(x)
    nc = x.sum()
    if nc < 2 or n - nc < 2:
        return np.full(Y.shape[1], np.nan), "few"
    Qz = R._basis(Z)
    xt = x - Qz @ (Qz.T @ x)
    sxx = xt @ xt
    if sxx < 1e-10 * max(1.0, x @ x):
        return np.full(Y.shape[1], np.nan), "span"
    h = (Qz ** 2).sum(1) + xt ** 2 / sxx
    if np.any(h > 1 - 1e-10):
        return np.full(Y.shape[1], np.nan), "lev1"
    c = xt / sxx
    b = c @ Y
    E = Y - Qz @ (Qz.T @ Y) - np.outer(xt, b)
    se = np.sqrt(((c / (1 - h))[:, None] ** 2 * E ** 2).sum(0))
    return b / se, ""


def design(names, src, asm_ids, rows_ids, ctx_vals):
    cov = R.my_covariates(names, src, asm_ids, rows_ids).reset_index(drop=True)
    x, Z, ok, _ = R.build_design(pd.Series(ctx_vals).reset_index(drop=True), cov)
    return x, Z, ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--B", type=int, default=2000)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    data.load_gene_effect = R._mmap_gene_effect
    drv = R.load_driver()
    Aset = set(assoc.a_models(assoc.Project()))
    oi = assoc.Project.__init__

    def init(self, root=paths.PROJECT, outcome_loader=None):
        oi(self, root=root, outcome_loader=outcome_loader or R.guarded_loader(Aset))
    assoc.Project.__init__ = init
    world = drv.real_world()
    world.proj.outcome_loader = R.guarded_loader(Aset)
    eB, _ = drv.entry(world, BENCH)
    pairs = pd.read_csv(OUTC / BENCH / "bench_pairs.tsv", sep="\t")
    n0 = len(R.LOADER_CALLS)
    E = assoc.resume_lockbox_read(BENCH, proj=world.proj, relook="R6")
    access = dict(call="assoc.resume_lockbox_read(BENCH, relook='R6')", n_models=int(E.shape[0]),
                  loader_calls_during=len(R.LOADER_CALLS) - n0, any_A_rows=bool(set(E.index) & Aset))
    L = E.index
    src = world.cov
    # per-pair L outcomes (full L rows)
    yc, yg = {}, {}
    for r in pairs.itertuples(index=False):
        rem = cscore.removal_for_context(r.context, world.homdel.get(r.context))
        yc[(r.context, r.unit)] = cscore.unit_scores(E, world.fmap, [r.unit], rem, world.annot)[0][r.unit].values
        yg[(r.context, r.unit)] = E[r.gstar].astype(float).clip(-3, 3).values
    ctxs = sorted(set(pairs.context))
    byc = {c: pairs.index[pairs.context == c].values for c in ctxs}
    names = {c: covmod.covariates_for(eB, c)[0] for c in ctxs}
    ctxv = {c: world.ctx[c].reindex(L).values for c in ctxs}
    sign = np.sign(pairs.beta_A.values)

    def deltas(rows, mode, stats):
        d = np.full(len(pairs), np.nan)
        ids = L[rows]
        for c in ctxs:
            asm = L if mode == "executed" else pd.Index(pd.unique(ids))
            x, Z, ok = design(names[c], src, asm, ids, ctxv[c][rows])
            idx = byc[c]
            YC = np.column_stack([yc[(pairs.context[i], pairs.unit[i])][rows][ok] for i in idx])
            YG = np.column_stack([yg[(pairs.context[i], pairs.unit[i])][rows][ok] for i in idx])
            Y = np.column_stack([YC, YG])
            good = ~np.isnan(Y).any(1)
            if not good.all():
                stats["y_nan_rows"] += 1
            tt, why = hc3_t_many(Y[good], x[good], Z[good])
            if why:
                stats[why] = stats.get(why, 0) + 1
            k = len(idx)
            d[idx] = sign[idx] * (tt[:k] - tt[k:])
        return d

    st0 = {"y_nan_rows": 0}
    d0 = deltas(np.arange(len(L)), "executed", st0)
    terc = pairs.tercile.values
    out = dict(access=access, point=dict(mean=float(np.nanmean(d0)), max_abs_diff_vs_committed=float(
        np.nanmax(np.abs(d0 - pairs.delta.values))), **{f"mean_{t}": float(np.nanmean(d0[terc == t]))
                                                        for t in ("low", "mid", "high")}))
    pairs.assign(delta_recomputed=d0).to_csv(a.out / "bench_pairs_recompute.tsv", sep="\t", index=False,
                                             float_format="%.6g")
    bl = src.PatientID.reindex(L)
    bl = bl.where(bl.notna(), pd.Series(L, index=L)).astype(str)
    groups = {b: np.flatnonzero(bl.values == b) for b in pd.unique(bl.values)}
    gnames = list(groups)
    rng = np.random.default_rng(int(eB["seeds"]["bootstrap"]))
    boots = {"executed": [], "assembled": []}
    stats = {"executed": {"y_nan_rows": 0}, "assembled": {"y_nan_rows": 0}}
    n_rep_with_nan = {"executed": 0, "assembled": 0}
    for b in range(a.B):
        pick = rng.integers(len(gnames), size=len(gnames))
        rows = np.concatenate([groups[gnames[i]] for i in pick])
        for mode in ("executed", "assembled"):
            d = deltas(rows, mode, stats[mode])
            n_rep_with_nan[mode] += int(np.isnan(d).any())
            boots[mode].append(d)
        if b % 200 == 0:
            print(f"[{time.time()-t0:.0f}s] rep {b}", flush=True)
    for mode in boots:
        Bm = np.array(boots[mode])
        s = dict(ci95=np.nanquantile(np.nanmean(Bm, 1), [0.025, 0.975]).tolist(), stats=stats[mode],
                 n_rep_with_any_nan_pair=n_rep_with_nan[mode], frac_pair_rep_nan=float(np.isnan(Bm).mean()))
        for t in ("low", "mid", "high"):
            m = terc == t
            s[f"ci95_{t}"] = np.nanquantile(np.nanmean(Bm[:, m], 1), [0.025, 0.975]).tolist()
        out[f"bootstrap_{mode}"] = s
    if a.B:
        diff = np.array(boots["executed"]) - np.array(boots["assembled"])
        out["executed_minus_assembled"] = dict(max_abs=float(np.nanmax(np.abs(diff))),
                                               frac_pair_rep_differs=float((np.abs(np.nan_to_num(diff, nan=9)) >
                                                                            1e-9).mean()))
    out["B"] = a.B
    out["runtime_s"] = round(time.time() - t0, 1)
    out["loader_calls"] = R.LOADER_CALLS
    (a.out / "bench_recompute.json").write_text(json.dumps(out, indent=1, default=str) + "\n")
    print(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
