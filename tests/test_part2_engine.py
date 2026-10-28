"""Association engine on toy data with known answers, in a throw-away git repository: per-context score/covariate
frames, evaluability (excluded from m), (context, unit) pairs, BH / Benjamini-Bogomolov / stage-wise BB levels,
permutation re-entry, outputs (data_cache table, counts-only summary, trials without hits), robustness, pseudo-contexts,
lockbox one-sided test, zero-target close, target-list / close rules."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))

from part2 import assoc, confirm, registry, robust, stats, synth  # noqa: E402
from part2.registry import RegistryError  # noqa: E402
from test_part2_assoc import A_MODELS, GENES, L_MODELS, commit_all, proj, sh  # noqa: E402,F401

UNITS = [f"U{i}" for i in range(12)]


def add(proj, aid, direction="H3b", contexts=("c1", "c2"), units=UNITS, fdr=None, sign="two-sided", ctx_items=None):
    fdr = fdr or dict(method="benjamini_bogomolov", families="contexts", q=0.05)
    items = ctx_items or [dict(id=c) for c in contexts]
    txt = synth.entry_yaml(aid, direction, items, units, "S_C_mean_winsor3", ["library"], fdr,
                           ["results/part2/contexts/checksums.json"], sign=sign)
    (proj.registry_dir / f"{aid}.yaml").write_text(txt)
    commit_all(proj.root, f"register {aid}")


def data(seed=0, eff=0.6):
    rng = np.random.default_rng(seed)
    n = len(A_MODELS)
    C = pd.DataFrame({"c1": (rng.random(n) < 0.2).astype(float), "c2": (rng.random(n) < 0.03).astype(float)},
                     index=A_MODELS)
    Y = pd.DataFrame(rng.normal(size=(n, len(UNITS))), index=A_MODELS, columns=UNITS)
    Y1 = Y.copy()
    Y1["U0"] += -eff * C.c1.values                                      # planted for c1 only
    Y2 = Y.copy()
    Y2["U5"] = np.nan                                                   # U5 not evaluable for c2 (removals)
    Z = pd.DataFrame({"library:x": (rng.random(n) < 0.3).astype(float)}, index=A_MODELS)
    return {"c1": Y1, "c2": Y2}, C, {"c1": Z, "c2": Z.iloc[:, :0]}


INP = ["data_cache/ctx.txt"]


def test_per_context_frames_evaluability_and_outputs(proj):
    add(proj, "H3b_e")
    Y, C, Z = data()
    res = assoc.run_association("H3b_e", Y, C, Z, inputs=INP, proj=proj, min_context_n=12, min_ref_n=12)
    assert len(res) == 2 * len(UNITS)
    c2 = res[res.context == "c2"]
    assert not c2.evaluable.any() and (c2.reason.str.contains("below") | (c2.reason == "no score")).all()
    assert c2.set_index("unit").reason["U5"] == "no score"                          # all-NaN column for c2
    assert res.q_adj[~res.evaluable].isna().all()                                    # excluded from m
    hit = res[(res.context == "c1") & (res.unit == "U0")].iloc[0]
    m = C.c1.notna()
    X = np.column_stack([np.ones(m.sum()), C.c1[m], Z["c1"]["library:x"][m]])
    r = stats.ols_hc3(Y["c1"]["U0"][m].values, X, 1)
    assert np.isclose(r["t"], hit.t) and hit.discovery
    bb = stats.benjamini_bogomolov(res.p[res.evaluable].values, res.context[res.evaluable].values, 0.05)
    assert bb["M"] == 1 and res.discovery.sum() == bb["reject"].sum() == 1
    # outputs: full table in data_cache, counts-only summary, trials without unit / context names
    tab = pd.read_parquet(proj.root / "data_cache/part2/assoc/H3b_e/discovery.parquet")
    assert len(tab) == len(res)
    summ = json.loads((proj.root / "results/part2/assoc/H3b_e/discovery_summary.json").read_text())
    assert summ["n_discovery"] == 1 and summ["n_not_evaluable"] == len(UNITS) and "discoveries" not in summ
    tr = proj.trials.read_text()
    assert '"U0"' not in tr and '"c1"' not in tr and json.loads(tr.splitlines()[-1])["n_evaluable"] == len(UNITS)
    res2 = assoc.run_association("H3b_e", Y, C, Z, inputs=INP, proj=proj, min_context_n=12, min_ref_n=12,
                                 tag="released", release_hits=True)
    summ = json.loads((proj.root / "results/part2/assoc/H3b_e/released_summary.json").read_text())
    assert [d["unit"] for d in summ["discoveries"]] == ["U0"] and len(res2) == len(res)


def test_pairs_and_bh_with_permutation_reentry(proj):
    add(proj, "H1_e", direction="H1", fdr=dict(method="bh", q=0.05))
    Y, C, Z = data()
    pairs = [("c1", "U0"), ("c1", "U1"), ("c2", "U2")]
    res = assoc.run_association("H1_e", Y, C, Z, inputs=INP, proj=proj, pairs=pairs, n_perm=1000,
                                blocks=pd.Series(A_MODELS, index=A_MODELS),
                                strata=pd.Series(np.repeat(["a", "b", "c"], 100), index=A_MODELS))
    assert set(zip(res.context, res.unit)) == set(pairs)
    r0 = res.set_index("unit").loc["U0"]
    assert r0.p_perm < 0.01 and r0.discovery and "q_adj_perm" in res
    with pytest.raises(ValueError, match="pair unit"):
        assoc.run_association("H1_e", Y, C, Z, inputs=INP, proj=proj, pairs=[("c1", "nope")])


def test_stagewise_bb(proj):
    subs = ["P1.S1", "P1.S2", "P2.S1", "P2.S2", "P3.S1", "P3.S2"]
    parent_of = {s: s.split(".")[0] for s in subs}
    add(proj, "H3a_p", direction="H3a", units=["P1", "P2", "P3"])
    add(proj, "H2_e", direction="H2", units=subs,
        fdr=dict(method="stagewise_bb", parent_directions=["H3a"], q=0.05))
    rng = np.random.default_rng(1)
    n = len(A_MODELS)
    C = pd.DataFrame({"c1": (rng.random(n) < 0.25).astype(float), "c2": (rng.random(n) < 0.25).astype(float)},
                     index=A_MODELS)
    P = pd.DataFrame(rng.normal(size=(n, 3)), index=A_MODELS, columns=["P1", "P2", "P3"])
    P["P1"] -= 0.8 * C.c1
    S = pd.DataFrame(rng.normal(size=(n, 6)), index=A_MODELS, columns=subs)
    S["P1.S1"] -= 0.8 * C.c1
    S["P3.S1"] -= 0.8 * C.c1                                            # strong, but parent P3 not discovered
    with pytest.raises(RegistryError, match="parent_selection"):
        assoc.run_association("H2_e", S, C, None, inputs=INP, proj=proj)
    resP = assoc.run_association("H3a_p", P, C, None, inputs=INP, proj=proj)
    sel = assoc.parent_selection_from({"H3a": resP}, parent_of, q=0.05)
    assert sel["selected"] == {("c1", "P1")}
    R, m = int(resP.discovery.sum()), int(resP.evaluable.sum())
    assert sel["level"]["c1"] == pytest.approx(0.05 * R / m)
    with pytest.raises(RegistryError, match="committed discovery run"):          # trials not committed yet
        assoc.run_association("H2_e", S, C, None, inputs=INP, proj=proj, parent_selection=sel)
    commit_all(proj.root, "H3a discovery committed")
    res = assoc.run_association("H2_e", S, C, None, inputs=INP, proj=proj, parent_selection=sel)
    assert set(zip(res.context, res.unit)) == {("c1", "P1.S1"), ("c1", "P1.S2")}   # only eligible tests run
    rj, _ = stats.bh(res.p.values, sel["level"]["c1"])
    assert res.reject.tolist() == rj.tolist() and res.set_index("unit").loc["P1.S1"].discovery
    assert not ((res.unit == "P3.S1").any())                                         # no open sub-module screen


def test_robustness_and_pseudo_contexts():
    rng = np.random.default_rng(2)
    n = 600
    lin = pd.Series(np.repeat(list("abcdef"), 100))
    ctx = pd.Series((rng.random(n) < 0.3).astype(float))
    y_broad = pd.Series(rng.normal(size=n) - 0.5 * ctx)
    y_one = pd.Series(rng.normal(size=n) - 2.0 * ctx * (lin == "a"))             # driven by one lineage
    assert robust.loso(y_broad, ctx, None, lin)["loso_robust"]
    assert not robust.loso(y_one, ctx, None, lin)["loso_robust"]
    y_out = pd.Series(rng.normal(size=n))
    top = ctx[ctx == 1].index[:2]
    y_out[top] = -40.0                                                           # two extreme context models
    assert not robust.leave_top2_out(y_out, ctx, None)["top2_robust"]
    assert robust.leave_top2_out(y_broad, ctx, None)["top2_robust"]
    E = pd.DataFrame({f"g{i}": rng.normal(size=n) - (0.6 * ctx if i < 3 else 0) for i in range(6)})
    assert robust.breadth(E, list(E.columns), ctx, None, -1)["breadth"] == pytest.approx(0.5, abs=0.17)
    blocks = pd.Series([f"p{i // 2}" for i in range(n)])
    ps = robust.pseudo_contexts(ctx, 20, 7, strata=lin, blocks=blocks)
    for k in ps:
        for l in "abcdef":
            got, want = ps[k][lin == l].sum(), ctx[lin == l].sum()
            assert want <= got <= want + 1                                       # matched per stratum (blocks)
        assert all(ps[k][blocks == b].nunique() == 1 for b in blocks.unique()[:30])  # whole blocks


def test_lockbox_test_and_zero_target_close(proj):
    rng = np.random.default_rng(3)
    nL = len(L_MODELS)
    c = pd.Series((rng.random(nL) < 0.3).astype(float), index=L_MODELS)
    y = pd.Series(rng.normal(size=nL), index=L_MODELS) - 0.8 * c
    t = confirm.lockbox_test(y, c, None, beta_A=-0.5, n_perm=0)
    assert t["confirmed"] and t["beta_L"] < 0 and t["ci95_L"][1] < 0
    assert not confirm.lockbox_test(y, c, None, beta_A=+0.5, n_perm=0)["confirmed"]  # opposite A sign
    c_small = pd.Series(0.0, index=L_MODELS)
    c_small.iloc[:4] = 1
    t2 = confirm.lockbox_test(y, c_small, None, beta_A=-0.5, n_perm=1000,
                              blocks=pd.Series(L_MODELS, index=L_MODELS))
    assert not t2["testable"] and np.isfinite(t2["p_L_perm"])                    # < 5 L models; perm p used
    # zero-target close: ledger line, nothing read, then closed; a second target list is refused
    add(proj, "H3b_z")
    with pytest.raises(RegistryError, match="ledger"):
        confirm.close_entry(proj, "H3b_z")

    def boom(*a, **k):
        raise AssertionError("outcome loader must not be called for a zero-target close")
    proj.outcome_loader = boom
    confirm.write_target_list(proj, "H3b_z", [])
    with pytest.raises(RegistryError, match="already exists"):
        confirm.write_target_list(proj, "H3b_z", ["G1 (1)"])
    out = confirm.confirm_entry("H3b_z", pd.DataFrame(columns=["context", "unit", "beta", "alternative"]),
                                None, pd.DataFrame(), None, proj=proj)
    assert out.empty
    confirm.close_entry(proj, "H3b_z")
    assert confirm.is_closed(proj, "H3b_z")
    led = [json.loads(l) for l in proj.ledger.read_text().splitlines()]
    assert led[-1]["analysis_id"] == "H3b_z" and led[-1]["n_genes"] == 0
    with pytest.raises(RegistryError, match="closed"):
        assoc.run_association("H3b_z", *data()[:2], inputs=INP, proj=proj)
