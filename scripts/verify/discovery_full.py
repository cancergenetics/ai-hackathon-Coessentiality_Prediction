"""Full analytic recompute of every H3a and H3b test (8,441 + 7,340), vectorised per (context, complete-case mask);
own HC3 (numpy) and Bell-McCaffrey Satterthwaite df; BB on the recomputed analytic p. Does NOT import src/part2.

    .venv/bin/python scripts/verify/discovery_full.py --out <scratch>/discovery_review      # ~1-2 min
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as ss

sys.path.insert(0, str(Path(__file__).resolve().parent))
import discovery_lib as L  # noqa: E402
from discovery_review import Annot, H3A_COV, DIRS, h3b_removed, parse_registry_contexts, sym, ent, tests  # noqa: E402


def batch_fit(ctx, Ymat, Z):
    """Ymat: models x units with identical NaN mask. Returns beta, se, df, n."""
    y0 = pd.Series(Ymat.iloc[:, 0].values, index=Ymat.index)
    X, _, ok, _, keep = L.design(ctx, y0, Z)
    j = keep.index("context")
    Y = Ymat[ok.values].values
    G = np.linalg.inv(X.T @ X) @ X.T
    h = np.einsum("ij,ji->i", X, G)
    B = G @ Y
    Res = Y - X @ B
    w = (G[j] / (1 - h)) ** 2
    se = np.sqrt(w @ (Res ** 2))
    M = np.eye(len(h)) - X @ G
    df = np.sum(w * (1 - h)) ** 2 / (w @ (M ** 2) @ w)
    return B[j], se, df, int(ok.sum())


def run(kind, t, E, mem, ct, src, covs_of, removed_of):
    out = []
    for c, g in t.groupby("context", sort=False):
        Z = L.cov_frame(covs_of(c), src)
        ys = {}
        for u in g.unit:
            rem = removed_of(c, mem[u])
            kept = [m for m in mem[u] if m not in rem]
            ys[u] = L.score(E, kept) if len(kept) >= 2 else pd.Series(np.nan, index=E.index)
        Y = pd.DataFrame(ys)
        cc = ct[c].notna() & Z.notna().all(axis=1)
        masks = Y[cc].isna().apply(lambda col: hash(col.values.tobytes()))
        for _, us in masks.groupby(masks):
            b, se, df, n = batch_fit(ct[c], Y[list(us.index)], Z)
            for u, bb_, s_ in zip(us.index, b, se):
                out.append(dict(context=c, unit=u, n_mine=n, beta_mine=bb_, se_mine=s_, df_mine=df))
    r = pd.DataFrame(out)
    r["p_mine"] = 2 * ss.t.sf(np.abs(r.beta_mine / r.se_mine), r.df_mine)
    return t.merge(r, on=["context", "unit"], how="left")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    E = L.effects_A()
    units, mem, _ = L.membership()
    src = L.cov_source(list(E.index))
    ct = L.context_table().reindex(E.index)
    an = Annot()
    sym2ent = {sym(c): ent(c) for c in E.columns}
    h3b_cov = parse_registry_contexts(L.ROOT / "docs/registry/H3b_driver_family.yaml")
    rep = {}
    for kind in ("H3a", "H3b"):
        t = tests(DIRS[kind])
        if kind == "H3a":
            r = run(kind, t, E, mem, ct, src, lambda c: H3A_COV, lambda c, m: set())
        else:
            r = run(kind, t, E, mem, ct, src, lambda c: h3b_cov[c], lambda c, m: h3b_removed(c, m, an, sym2ent))
        rel = lambda a_, b_: ((r[a_] - r[b_]).abs() / r[b_].abs().clip(lower=1e-12))  # noqa: E731
        r["ok"] = (r.n == r.n_mine) & (rel("beta_mine", "beta") < 1e-4) & (rel("se_mine", "se") < 1e-4) & \
                  (rel("df_mine", "df") < 1e-4) & (rel("p_mine", "p") < 1e-3)
        rej, R, _ = L.bb(r.p_mine.values, r.context.values, 0.05)
        scr, _, _ = L.bb(r.p_mine.values, r.context.values, 0.1)
        bad = r[~r.ok]
        rep[kind] = dict(n_tests=len(r), n_match=int(r.ok.sum()), n_mismatch=len(bad),
                         mismatches=bad[["context", "unit", "n", "n_mine", "beta", "beta_mine", "p", "p_mine"]]
                         .head(30).to_dict("records"),
                         bb_reject_on_own_p=int(rej.sum()), equal_to_committed_reject=bool((rej == r.reject.values).all()),
                         screen_on_own_p=int(scr.sum()), equal_to_committed_screen=bool((scr == r.p_perm.notna().values).all()))
        r.to_csv(out / f"full_recompute_{kind}.tsv.gz", sep="\t", index=False, float_format="%.6g")
    (out / "full_recompute.json").write_text(json.dumps(rep, indent=1, default=str) + "\n")
    print(json.dumps(rep, indent=1, default=str))


if __name__ == "__main__":
    main()
