"""Lockbox driver (scripts/48_lockbox_confirm.py), BENCH secondary / BRIDGE-L additions (src/part2/bench.py) and
amendment A3 (scripts/49_amend_A3.py), tested before any real lockbox read.

Synthetic world (tests/part2_synth_lockbox.py): a planted lockbox confirmation / failure pattern, two rare contexts
(LoF:KMT2D, fusion:EWSR1::FLI1|ERG with < 5 L models), synthetic DEMETER2, a GLS comparator map; A-side runs by the
real runner. Checks: planted confirmations reproduced and the A-only effect fails; the rare-context route fires for
KMT2D and the EWSR1 fusion only; BENCH / BRIDGE-L refused out of order (driver gate and assoc.read_lockbox); re-reads
refused; ledger order; E2 labels never promoted. No DepMap data is read."""
import dataclasses
import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))

import part2_synth_lockbox as S  # noqa: E402
from part2 import assoc, bench, confirm, registry  # noqa: E402

REG = ROOT / "docs" / "registry"


def _script(name):
    return S._script(name)


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    base = tmp_path_factory.mktemp("lockbox_driver")
    world, info = S.build_world(base / "repo")
    D = info["D"]
    pre = D.preflight(world)
    # an untouched copy (A side committed, nothing read) for the out-of-order checks
    shutil.copytree(base / "repo", base / "copy", symlinks=True)
    wc = dataclasses.replace(world, proj=assoc.Project(root=base / "copy", outcome_loader=world.proj.outcome_loader),
                             log=[])
    out = D.run_all(world)
    return dict(world=world, info=info, D=D, pre=pre, out=out, copy=wc)


def _conf(world, aid):
    t = pd.read_csv(world.out(aid) / "confirmation.tsv", sep="\t")
    for c in ("eligible", "lockbox_testable", "confirmed_final"):
        t[c] = t[c].astype(str).str.lower() == "true"
    return t


def test_preflight_has_no_gaps_and_routes_rare_contexts(runs):
    pre = runs["pre"]
    assert pre["gaps"] == []
    rare = {c for v in pre["entries"].values() for c in v.get("rare_contexts", [])}
    assert rare == set(S.RARE)


def test_planted_confirmations_reproduced(runs):
    w, nm = runs["world"], runs["info"]["named"]
    t = _conf(w, "H3b_driver_family").set_index(["context", "unit"])
    for c, role in S.CONFIRM:
        assert t.loc[(c, nm[role]), "confirmed_final"], (c, role)
    c, u, _ = S.FAIL
    assert (c, u) in t.index, "planted A-only effect not discovered in A"
    assert not t.loc[(c, u), "confirmed_final"]
    assert t.loc[(c, u), "status"] == "not confirmed"
    # H1: every engine non-neutral call planted in A and L is confirmed on its defining statistic
    h = _conf(w, "H1_within_complex_GI")
    el = h[h.eligible]
    assert len(el) >= 3 and el.confirmed_final.all(), el[["unit", "class_engine", "p_L"]]
    assert set(el.axis + ":" + el.class_engine) >= {"P:PARALOG COMPENSATION", "N:COLLAPSE", "N:RESIDUAL SENSITISATION"}
    # neutral calls are never confirmed but carry L effects with 90% CIs
    ne = h[~h.eligible]
    assert not ne.confirmed_final.any() and ne.ci90_L_lo.notna().all()
    # controls planted in both folds confirm
    for k in S.CTRL_INDEX:
        assert _conf(w, k).confirmed_final.all(), k


def test_rare_route_fires_for_kmt2d_and_ewsr1_only(runs):
    w = runs["world"]
    seen = set()
    for p in sorted(w.out().glob("*/confirmation.tsv")):
        t = _conf(w, p.parent.name)
        assert set(t.context[~t.lockbox_testable]) <= set(S.RARE), (p.parent.name, set(t.context))
        rr = t[~t.lockbox_testable & t.eligible]
        assert rr.rare_route_pass.notna().all()
        assert t[t.lockbox_testable].rare_route_pass.isna().all()
        assert (rr.status.str.startswith("not lockbox-testable")).all()
        assert not rr.confirmed_final.any()                  # labelled, never 'confirmed'
        seen |= set(rr.context)
    assert seen == set(S.RARE)
    h3b = _conf(w, "H3b_driver_family")
    for c, u in ((c, v[0]) for c, v in S.RARE.items()):
        r = h3b[(h3b.context == c) & (h3b.unit == u)]
        assert len(r) == 1 and bool(r.rare_route_pass.iloc[0]), (c, u)


def test_bench_bridge_refused_out_of_order(runs):
    D, wc = runs["D"], runs["copy"]
    with pytest.raises(assoc.FirewallError):
        D.bench_stage(wc)
    with pytest.raises(assoc.FirewallError):
        D.bridge_stage(wc)
    assert not wc.proj.ledger.exists()
    # machinery refusal at the read: a committed BENCH target list is still refused
    confirm.write_target_list(wc.proj, "BENCH_complex_vs_gene", ["X (1)"])
    with pytest.raises(assoc.FirewallError, match="before all confirmatory reads"):
        assoc.read_lockbox("BENCH_complex_vs_gene", confirm.target_path(wc.proj, "BENCH_complex_vs_gene"),
                           proj=wc.proj)
    # after one confirmation only: still refused, and the refusal names the missing entries
    order = D.confirm_order(wc.proj)
    D.run_entry(wc, D.build_spec(wc, order[0]))
    with pytest.raises(assoc.FirewallError, match="missing"):
        D.bridge_stage(wc)
    led = [json.loads(l) for l in wc.proj.ledger.read_text().splitlines()]
    assert [r["analysis_id"] for r in led] == [order[0]]


def test_rereads_refused_and_ledger_order(runs):
    w, D = runs["world"], runs["D"]
    led = [json.loads(l) for l in w.proj.ledger.read_text().splitlines()]
    ids = [r["analysis_id"] for r in led]
    assert ids == D.confirm_order(w.proj) + ["BENCH_complex_vs_gene", "BRIDGE-L_confirmation_rate"]
    assert ids[:4] == sorted(S.CTRL_INDEX, key=S.CTRL_INDEX.get)
    for aid in ids:
        assert confirm.is_closed(w.proj, aid), aid
    with pytest.raises(assoc.FirewallError, match="already read"):
        assoc.read_lockbox("H3a_lineage_atlas", confirm.target_path(w.proj, "H3a_lineage_atlas"), proj=w.proj)
    with pytest.raises(D.GapError, match="ledger line already exists"):
        D.run_all(w)
    with pytest.raises(assoc.FirewallError):
        D.bench_stage(w)


def test_e2_labels_never_promoted(runs):
    h = _conf(runs["world"], "H1_within_complex_GI")
    e2 = h[h.label_final == S.E2_LABEL]
    assert len(e2) >= 1 and e2.confirmed_final.all()          # confirmed on the registered statistic ...
    assert (e2.headline_eligible.astype(str) == "False").all()  # ... but the label and headline flag are unchanged
    ev = pd.read_csv(runs["world"].out() / "H1_within_complex_GI_events.tsv", sep="\t")
    assert (ev[ev.event.isin(e2.event)].label_final == S.E2_LABEL).all()


def test_bench_and_bridge_outputs(runs):
    out = runs["out"]
    b, g = out["bench"], out["bridge"]
    assert np.isfinite(b["primary"]["mean"]) and len(b["primary"]["ci95"]) == 2
    assert b["multiverse"]["within_line_rank"].startswith("NOT RUN")
    assert set(b["multiverse"]) >= {"eigengene", "gstar_largest_beta", "T1_T2_pairs"}
    assert np.isfinite(b["secondary"]["r_complex"]) and len(b["secondary"]["r_random"]) == 5
    assert b["distinct"]["n_distinct_gstar"] >= 1
    assert np.isfinite(g["primary"]["diff"]) and len(g["primary"]["ci95"]) == 2
    assert set(g["multiverse"]) >= {"unweighted", "matched_gene_coverage", "T1_T2_clusters"}
    rep = out["report"]
    assert rep["H3b_driver_family"]["overall"]["n_rare"] == 2
    # the GLS discovery pass is two trial lines (one BB family each), committed before the BRIDGE-L read
    w = runs["world"]
    tr = [json.loads(l) for l in w.proj.trials.read_text().splitlines()]
    gl = [r for r in tr if r["analysis_id"] == "BRIDGE-L_confirmation_rate"]
    assert [r["tag"] for r in gl] == ["discovery_H3a", "discovery_H3b"]


# ------------------------------------------------------------------------------------------------ unit pieces
def test_bridge_rate_uses_map_column():
    """Regression: `h.map` is DataFrame.map (a method) in pandas >= 2.1; the map column must be indexed."""
    h = pd.DataFrame(dict(map=["zavg"] * 4 + ["gls"] * 4, size=[4, 4, 10, 10] * 2,
                          confirmed=[True, True, True, False, True, False, False, False]))
    r = bench.bridge_rate(h, B=50, seed=1)
    assert r["rate_zavg"] == pytest.approx(0.75) and r["rate_gls"] == pytest.approx(0.25)
    u = bench.bridge_unweighted(h, B=20, seed=1)
    assert u["diff"] == pytest.approx(0.5) and len(u["ci95"]) == 2


def test_effect_rsd_and_al_correlation():
    rng = np.random.default_rng(0)
    n = 400
    idx = pd.Index([f"m{i}" for i in range(n)])
    c = pd.Series((rng.random(n) < 0.3).astype(float), idx)
    y = pd.Series(0.5 * c.values + rng.normal(0, 2.0, n), idx)
    e = bench.effect_rsd(y, c, None)
    assert e == pytest.approx(0.25, abs=0.12)
    pairs = pd.DataFrame(dict(context=["c"] * 3, unit=["u1", "u2", "u3"]))
    ctx = pd.DataFrame(dict(c=c))
    eff = {"u1": 0.2, "u2": 0.5, "u3": 1.0}
    sc = lambda cc, u: pd.Series(eff[u] * c.values + rng.normal(0, 1, n), idx)  # noqa: E731
    r = bench.al_correlation(pairs, sc, sc, ctx, ctx, {}, {}, method="pearson")
    assert r["r_complex"] > 0.8
    with pytest.raises(ValueError):
        bench.al_correlation(pairs, sc, sc, ctx, ctx, {}, {}, method=None)


def test_real_registry_gaps_are_refused():
    D = _script("48_lockbox_confirm.py")
    eb = registry.load_yaml((REG / "BENCH_complex_vs_gene.yaml").read_text())
    with pytest.raises(D.GapError):
        D.check_bench_rulings(eb.get("estimand_secondary"), (eb.get("multiverse") or {}).get("ci"))


# ------------------------------------------------------------------------------------------------ amendment A3
def _show(commit, path):
    r = subprocess.run(["git", "-C", str(ROOT), "show", f"{commit}:{path}"], capture_output=True, text=True)
    return r.stdout if r.returncode == 0 else None


def test_amendment_A3_seeds_and_bench_not_run():
    A3 = _script("49_amend_A3.py")
    entries = registry.load_all(REG)
    seeds = {k: e.get("seeds", {}).get("lockbox_permutation") for k, e in entries.items()}
    exp = dict(A3.LOCKBOX_SEEDS)
    for k, e in entries.items():
        if k.startswith("CTRL_"):
            exp[k] = 6700 + int(e["seeds"]["permutation"]) - 4400
    assert {k: v for k, v in seeds.items() if v is not None} == exp
    assert exp["CTRL_PC1"] == 6700 and exp["CTRL_PC22"] == 6722 and exp["CTRL_PC19b"] == 6719
    other = A3.registered_seeds({k: dict(e, seeds={a: b for a, b in e.get("seeds", {}).items()
                                                   if a != "lockbox_permutation"}) for k, e in entries.items()})
    assert not set(exp.values()) & (other | A3.documented_seeds())
    assert len(set(exp.values())) == len(exp)
    eb = entries["BENCH_complex_vs_gene"]
    assert eb["multiverse"]["not_run"] == ["within_line_rank"]
    for k in list(exp) + ["BENCH_complex_vs_gene"]:
        rec = [a for a in entries[k]["amendments"] if a["id"] == "A3"]
        assert len(rec) == 1 and rec[0]["amends_commit"] == A3.AMENDS and rec[0]["timestamp_utc"] == A3.TIMESTAMP, k
        assert "before any lockbox read; no L values read" in rec[0]["reason"], k
        assert registry.validate(entries[k], k) == [], k
        old_txt = _show(A3.AMENDS, f"docs/registry/{k}.yaml")
        if old_txt is None:
            pytest.skip("no history with the amended commit")
        old = registry.load_yaml(old_txt)
        assert registry.amendment_changes(old, entries[k]) == [], k
        assert registry.changed_top_paths(old, entries[k], rec[0]["changed_paths"]) == set(rec[0]["changed_paths"]), k
    txt = (REG / "AMENDMENTS.md").read_text()
    assert "## A3" in txt and A3.AMENDS in txt and A3.TIMESTAMP in txt
    assert "before any lockbox read; no L values read" in txt
    # the script re-derives exactly the committed text
    if _show(A3.AMENDS, "docs/registry/H3a_lineage_atlas.yaml") is not None:
        plan = A3.plan()
        assert all((REG / f"{k}.yaml").read_text() == t for k, (t, *_) in plan.items())
