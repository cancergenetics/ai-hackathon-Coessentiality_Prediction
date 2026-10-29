"""Lockbox driver (scripts/51_lockbox_confirm.py), BENCH secondary / BRIDGE-L additions (src/part2/bench.py) and
amendment A3 (scripts/52_amend_A3.py), tested before any real lockbox read.

Synthetic world (tests/part2_synth_lockbox.py): a planted lockbox confirmation / failure pattern, two rare contexts
(LoF:KMT2D, fusion:EWSR1::FLI1|ERG with < 5 L models), synthetic DEMETER2, a GLS comparator map; A-side runs by the
real runner. Checks: planted confirmations reproduced and the A-only effect fails; the rare-context route fires for
KMT2D and the EWSR1 fusion only; BENCH / BRIDGE-L refused out of order (driver gate and assoc.read_lockbox); re-reads
refused; ledger order; E2 labels never promoted. No DepMap data is read."""
import dataclasses
import types
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
    # A4/A5: the A-L random-set correlation NOT RUN, the distinct-complex count reported (A only); BENCH multiverse
    # point estimates only (primary keeps its CIs); BRIDGE-L multiverse with bootstrap CIs (A5)
    assert b["secondary"].startswith("NOT RUN") and b["distinct"]["n_distinct_gstar"] >= 1
    assert not (runs["world"].out("BENCH_complex_vs_gene") / "al_correlation_pairs.tsv").exists()
    for k in ("eigengene", "gstar_largest_beta", "T1_T2_pairs"):
        assert np.isfinite(b["multiverse"][k]["mean"]) and not any(x.startswith("ci95") for x in b["multiverse"][k])
    assert np.isfinite(g["primary"]["diff"]) and len(g["primary"]["ci95"]) == 2
    assert set(g["multiverse"]) >= {"unweighted", "matched_gene_coverage", "T1_T2_clusters"}
    assert all(len(v["ci95"]) == 2 for k, v in g["multiverse"].items() if "note" not in v)
    rep = out["report"]
    assert rep["H3b_driver_family"]["overall"]["n_rare"] == 2
    # the GLS discovery pass is two trial lines (one BB family each), committed before the BRIDGE-L read
    w = runs["world"]
    tr = [json.loads(l) for l in w.proj.trials.read_text().splitlines()]
    gl = [r for r in tr if r["analysis_id"] == "BRIDGE-L_confirmation_rate"]
    assert [r["tag"] for r in gl] == ["discovery_H3a", "discovery_H3b"]


def test_rare_route_with_A_alignment_and_missing_demeter2(runs):
    """A5 (G3): an eligible N-axis call on a non-testable event takes the rare route with the A CRISPR alignment
    (baselines and N' from A) on DEMETER2, with LOSO and leave-one-model-out. N5: when DEMETER2 is not available the
    route is NOT EVALUABLE, not 'failed'. A side only (synthetic H1 events are all lockbox-testable, so the routes are
    forced here)."""
    D, w = runs["D"], runs["world"]
    sp = D.build_spec(w, "H1_within_complex_GI")
    h = sp.hits
    n_i = h.index[h.eligible & (h.axis == "N")][0]
    p_i = h.index[h.eligible & (h.axis == "P")][0]
    h.loc[[n_i, p_i], "lockbox_testable"] = False
    sp.needs_fit = [h.unit[n_i]]
    rr = D.rare_routes(w, sp)
    for i in (n_i, p_i):
        assert rr[i]["demeter2_same_sign"] is not None and rr[i]["loso_folds"] > 0 and rr[i]["rare_route_note"] == ""
        assert isinstance(rr[i]["rare_route_pass"], bool)
    eid = h.unit[n_i].rsplit("|", 1)[0]
    import part2.h1 as h1m
    y = sp.score_D2(h.context[n_i], h.unit[n_i])            # = aligned score on DEMETER2 with the A parameters
    ev = {e.id: e for e in w.events["H1_within_complex_GI"]}[eid]
    lossA = w.loss.reindex(w.A).rename(columns=lambda c: c.split(":", 1)[1])
    cA = h1m.component_scores(w.E_A, ev, lossA[ev.gene])
    pd.testing.assert_series_equal(y, h1m.aligned_score(w.D2_A(), cA.Nprime, cA.b))
    sp.score_D2 = lambda c, u: None                         # gene not in DEMETER2
    rr = D.rare_routes(w, sp)
    assert rr[n_i]["rare_route_pass"] is None and rr[n_i]["rare_route_note"] == D.D2_NA
    row = h.loc[n_i].copy()
    row["rare_route_pass"], row["testable_L"], row["confirmed_final"] = None, False, False
    assert D._status(row) == f"{D.RARE_LABEL}: rare-context route {D.D2_NA}"


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


def test_real_registry_bench_rulings():
    """Without A4 the BENCH secondary gap is refused; with A4 (secondary NOT RUN, multiverse point estimates) the
    guard passes."""
    D = _script("51_lockbox_confirm.py")
    eb = registry.load_yaml((REG / "BENCH_complex_vs_gene.yaml").read_text())
    ci = (eb.get("multiverse") or {}).get("ci")
    with pytest.raises(D.GapError):
        D.check_bench_rulings(eb.get("estimand_secondary"), ci)
    with pytest.raises(D.GapError):
        D.check_bench_rulings(None, None, eb.get("estimand_secondary_not_run"))
    D.check_bench_rulings(eb.get("estimand_secondary"), ci, eb.get("estimand_secondary_not_run"))


# ------------------------------------------------------------------------------------------------ amendment A3
def _show(commit, path):
    r = subprocess.run(["git", "-C", str(ROOT), "show", f"{commit}:{path}"], capture_output=True, text=True)
    return r.stdout if r.returncode == 0 else None


def _frozen_entries():
    """Registry entries as frozen before the read (the post-read close changes only `status`)."""
    return {k: dict(e, status="open") for k, e in registry.load_all(REG).items()}


def _frozen_text(k):
    return (REG / f"{k}.yaml").read_text().replace("\nstatus: closed\n", "\nstatus: open\n")


def test_amendment_A3_seeds_and_bench_not_run():
    A3 = _script("52_amend_A3.py")
    entries = _frozen_entries()
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
        plan, c = A3.plan(), A3.applied_commit()
        assert all(_show(c, f"docs/registry/{k}.yaml") == t for k, (t, *_) in plan.items())


# ------------------------------------------------------------------------------------------------ amendment A4
def test_amendment_A4_rulings():
    """A4 as applied (its own commit); A5 later amends the same fields (test_amendment_A5_rulings)."""
    A4 = _script("53_amend_A4.py")
    c = A4.applied_commit()
    if _show(c, "docs/registry/BENCH_complex_vs_gene.yaml") is None:
        pytest.skip("no history with the A4 commit")
    entries = {p.stem: registry.load_yaml(_show(c, f"docs/registry/{p.name}")) for p in sorted(REG.glob("*.yaml"))}
    eb, eg, es = (entries[k] for k in ("BENCH_complex_vs_gene", "BRIDGE-L_confirmation_rate", "H1S_silencing_GI"))
    assert len(eb["estimand_secondary_not_run"]) == 2 and "estimand_secondary" not in eb
    assert eb["multiverse"]["ci"] == eg["multiverse"]["ci"] == "point estimate only"
    assert es["lockbox"]["rare_route_fitted_parameters"] == "not evaluable"
    assert sorted(es["lockbox"]["rare_route_not_evaluable"]) == sorted(A4.H1S_NOT_EVALUABLE)
    # every other entry is untouched by A4
    assert [k for k, e in entries.items() if any(a["id"] == "A4" for a in e.get("amendments", []))] == \
        sorted(A4.REASONS)
    for k in A4.REASONS:
        rec = [a for a in entries[k]["amendments"] if a["id"] == "A4"]
        assert len(rec) == 1 and rec[0]["amends_commit"] == A4.AMENDS and rec[0]["timestamp_utc"] == A4.TIMESTAMP
        assert "before any lockbox read; no L values read" in rec[0]["reason"], k
        assert registry.validate(entries[k], k) == [], k
        old_txt = _show(A4.AMENDS, f"docs/registry/{k}.yaml")
        if old_txt is None:
            pytest.skip("no history with the amended commit")
        old = registry.load_yaml(old_txt)
        assert registry.amendment_changes(old, entries[k]) == [], k
        assert registry.changed_top_paths(old, entries[k], rec[0]["changed_paths"]) == set(rec[0]["changed_paths"])
        assert _show(A4.applied_commit(), f"docs/registry/{k}.yaml") == A4.plan()[k][0], k
    txt = (REG / "AMENDMENTS.md").read_text()
    assert "## A4" in txt and A4.AMENDS in txt and A4.TIMESTAMP in txt
    a4 = txt[txt.index("## A4"):]
    assert "before any lockbox read; no L values read" in a4 and "accepted literal readings" in a4.lower()


# ------------------------------------------------------------------------------------------------ amendment A5
def test_amendment_A5_rulings():
    A5 = _script("54_amend_A5.py")
    entries = _frozen_entries()
    eb, eg, es = (entries[k] for k in ("BENCH_complex_vs_gene", "BRIDGE-L_confirmation_rate", "H1S_silencing_GI"))
    assert len(eb["estimand_secondary_not_run"]) == 1 and "A-L correlation" in eb["estimand_secondary_not_run"][0]
    assert "non-member" in eb["estimand_secondary_not_run_reason"] and "Erratum" in eb["estimand_secondary_not_run_reason"]
    assert eb["multiverse"]["ci"] == "point estimate only" and eg["multiverse"]["ci"] == "bootstrap as primary"
    assert es["lockbox"]["rare_route_fitted_parameters"] == "A CRISPR parameters"
    assert "rare_route_not_evaluable" not in es["lockbox"]
    assert sorted(es["lockbox"]["rare_route_fitted_hits"]) == sorted(A5.H1S_FITTED)
    assert [k for k, e in entries.items() if any(a["id"] == "A5" for a in e.get("amendments", []))] == \
        sorted(A5.REASONS)
    for k in A5.REASONS:
        rec = [a for a in entries[k]["amendments"] if a["id"] == "A5"]
        assert len(rec) == 1 and rec[0]["amends_commit"] == A5.AMENDS and rec[0]["timestamp_utc"] == A5.TIMESTAMP
        assert "before any lockbox read; no L values read" in rec[0]["reason"], k
        assert registry.validate(entries[k], k) == [], k
        old_txt = _show(A5.AMENDS, f"docs/registry/{k}.yaml")
        if old_txt is None:
            pytest.skip("no history with the amended commit")
        old = registry.load_yaml(old_txt)
        assert registry.amendment_changes(old, entries[k]) == [], k
        assert registry.changed_top_paths(old, entries[k], rec[0]["changed_paths"]) == set(rec[0]["changed_paths"])
        assert _frozen_text(k) == A5.plan()[k][0], k
    txt = (REG / "AMENDMENTS.md").read_text()
    a5 = txt[txt.index("## A5"):]
    assert A5.AMENDS in a5 and A5.TIMESTAMP in a5 and "before any lockbox read; no L values read" in a5
    for s in ("Accepted readings", "heme_vs_breast", "PC17", "within lineage", "all of C's members", "N11"):
        assert s in a5, s


# ------------------------------------------------------------------------------------------------ MUST-1 / MUST-2
def test_target_list_records_eligibility_and_route(runs):
    w = runs["world"]
    rec = json.loads(confirm.target_path(w.proj, "H1_within_complex_GI").read_text())
    h = rec["hits"]
    assert all({"eligible", "route", "label_final", "headline_eligible", "event", "axis"} <= set(x) for x in h)
    t = _conf(w, "H1_within_complex_GI")
    assert [(x["unit"], x["eligible"], x["route"]) for x in h] == list(zip(t.unit, t.eligible, t.route))
    rec = json.loads(confirm.target_path(w.proj, "H3b_driver_family").read_text())
    assert {x["route"] for x in rec["hits"] if x["context"] in S.RARE} == {"not lockbox-testable"}


def test_e2_status_not_promotable(runs):
    h = _conf(runs["world"], "H1_within_complex_GI")
    e2 = h[h.label_final == S.E2_LABEL]
    assert len(e2) and (e2.status == "confirmed (not promotable: E2 label or not headline-eligible)").all()
    assert not e2.promotable.astype(str).eq("True").any()
    ok = h[h.confirmed_final & (h.label_final == h.class_engine) & (h.headline_eligible.astype(str) == "True")]
    assert len(ok) and ok.promotable.astype(str).eq("True").all() and (ok.status == "confirmed").all()


def test_reads_are_persisted_with_committed_hash(runs):
    w = runs["world"]
    sc = [json.loads(l) for l in (w.proj.results_dir / "lockbox_reads.jsonl").read_text().splitlines()]
    led = [json.loads(l) for l in w.proj.ledger.read_text().splitlines()]
    assert [r["analysis_id"] for r in sc] == [r["analysis_id"] for r in led if r["n_genes"] > 0]
    for r in sc:
        assert assoc._sha_file(w.proj.root / r["file"]) == r["sha256"]
    assert assoc._committed_text(w.proj, w.proj.results_dir / "lockbox_reads.jsonl").count("\n") == len(sc)


def test_crash_after_read_resume_without_loader(runs):
    """MUST-2a: a crash after a read leaves the persisted matrix with its committed hash; close_entry by hand is refused;
    --resume finishes the entry from the persisted matrix without calling the loader; a tampered file is refused."""
    D, base = runs["D"], runs["world"]
    root = base.proj.root.parent / "resume"
    shutil.copytree(runs["copy"].proj.root, root, symlinks=True)        # copy: one entry read and closed already
    calls = []
    inner = base.proj.outcome_loader
    A = set(base.A)

    def loader(ids, genes, source):
        if not set(ids) <= A:
            calls.append(source)
        return inner(ids, genes, source)

    w = dataclasses.replace(base, proj=assoc.Project(root=root, outcome_loader=loader), log=[], _d2=None)
    aid = D.confirm_order(w.proj)[1]
    sp = D.build_spec(w, aid)
    good = sp.score_L

    def boom(E_L, c, u):
        raise RuntimeError("simulated crash after the read")
    sp.score_L = boom
    with pytest.raises(RuntimeError, match="simulated crash"):
        D.run_entry(w, sp)
    assert calls == ["crispr"]
    with pytest.raises(registry.RegistryError, match="not committed"):
        confirm.close_entry(w.proj, aid)                       # cannot be closed by hand
    with pytest.raises(assoc.FirewallError):
        D.run_entry(w, D.build_spec(w, aid))                   # no resume flag: refused
    wr = dataclasses.replace(w, resume=True)
    rec = [json.loads(l) for l in (root / "results/part2/lockbox_reads.jsonl").read_text().splitlines()][-1]
    f = root / rec["file"]
    raw = f.read_bytes()
    f.write_bytes(raw + b"x")
    with pytest.raises(assoc.FirewallError, match="SHA-256"):
        D.run_entry(wr, D.build_spec(wr, aid))
    f.write_bytes(raw)
    out = D.run_entry(wr, D.build_spec(wr, aid))
    assert calls == ["crispr"] and confirm.is_closed(w.proj, aid) and len(out)
    ref = _conf(runs["world"], aid)                          # the uncrashed run on the same synthetic L
    assert list(out.confirmed_final) == list(ref.confirmed_final)
    assert json.loads((w.out(aid) / "summary.json").read_text())["resumed"] is True
    assert good is not None


def test_close_entry_refused_without_output(runs):
    D, base = runs["D"], runs["world"]
    root = base.proj.root.parent / "noclose"
    shutil.copytree(runs["copy"].proj.root, root, symlinks=True)
    w = dataclasses.replace(base, proj=assoc.Project(root=root, outcome_loader=base.proj.outcome_loader), log=[])
    aid = D.confirm_order(w.proj)[1]
    sp = D.build_spec(w, aid)
    units = {u: sp.members_of(u) for u in sorted(set(sp.hits.unit))}
    confirm.write_target_list(w.proj, aid, sorted({g for v in units.values() for g in v}),
                              hits=sp.hits[["context", "unit", "beta", "alternative"]], units=units)
    assoc.read_lockbox(aid, confirm.target_path(w.proj, aid), proj=w.proj)
    with pytest.raises(registry.RegistryError, match="not committed"):
        confirm.close_entry(w.proj, aid)


def test_bench_delta_fixed_terciles():
    rng = np.random.default_rng(3)
    n = 120
    idx = pd.Index([f"m{i}" for i in range(n)])
    E = pd.DataFrame(rng.standard_normal((n, 8)), idx, [f"g{j}" for j in range(8)])
    ctx = pd.DataFrame(dict(c=(rng.random(n) < 0.3).astype(float)), idx)
    pairs = pd.DataFrame(dict(context=["c"] * 4, unit=["u0", "u1", "u2", "u3"], beta_A=[-1.0] * 4,
                              gstar=["g4", "g5", "g6", "g7"], breadth=[0.1, 0.2, 0.3, 0.4]))
    sL = lambda E_, c, u: E_.iloc[:, int(u[1]):int(u[1]) + 2].mean(axis=1)  # noqa: E731
    fixed = ["high", "high", "low", "low"]
    r = bench.bench_delta(pairs, sL, E, ctx, None, B=5, seed=1, tercile=fixed)
    assert list(r["table"].tercile) == fixed and "mean_mid" not in r["summary"]
    with pytest.raises(ValueError):
        bench.bench_delta(pairs, sL, E, ctx, None, B=2, seed=1, tercile=fixed[:3])


# ------------------------------------------------------------------------------------------------ R6
def test_r6_l_covariates_assembled_on_L_rows():
    """R6: an L subset without the all-model reference lineage and with a single-model lineage. Assembled over all
    models and then subset (the defect), the single model becomes its own group and HC3 refuses (leverage 1);
    assembled over the L rows (the fix, as the A side does), the engine merges it into the reference and the fit
    works."""
    from part2 import covariates as covmod, stats
    from part2.robust import fit_one
    D = _script("51_lockbox_confirm.py")
    rng = np.random.default_rng(0)
    A = [f"A{i}" for i in range(60)]
    L = [f"L{i}" for i in range(41)]
    lin = ["Aaa"] * 20 + ["Bbb"] * 20 + ["Ccc"] * 20 + ["Bbb"] * 20 + ["Ccc"] * 20 + ["Zzz"]
    cov = pd.DataFrame(dict(OncotreeLineage=lin), index=A + L)
    ctx = pd.Series(np.r_[rng.random(60) < 0.3, rng.random(41) < 0.3].astype(float), index=A + L)
    ctx[L[:8]] = np.nan                                           # a context with NaN labels for some L models
    e = {"id": "toy", "covariates": ["lineage_FE"], "context_family": {"contexts": [{"id": "c"}]}}
    old = covmod.per_context(e, ["c"], cov)[0]["c"].reindex(L)
    X, _, _ = assoc._design(ctx.reindex(L), old)
    with pytest.raises(ValueError, match="leverage"):
        stats.HC3Design(X, 1)
    new = D.l_covariates(types.SimpleNamespace(cov=cov), e, ["c"], L)["c"]
    X2, _, _ = assoc._design(ctx.reindex(L), new)
    stats.HC3Design(X2, 1)                                        # no error
    y = pd.Series(rng.standard_normal(41), index=L)
    assert not fit_one(y, ctx.reindex(L), old)["ok"] and fit_one(y, ctx.reindex(L), new)["ok"]


def test_r6_relook_path_needs_committed_ruling_and_reproduces(runs):
    """The relook reads persisted matrices of CLOSED entries only with a committed PI note and intent naming the
    ruling, never calls the loader, writes no ledger line, and (synthetic: the fix changes no design there) reproduces
    every committed fit."""
    D, base = runs["D"], runs["world"]
    root = base.proj.root.parent / "relook"
    shutil.copytree(base.proj.root, root, symlinks=True)
    calls = []
    inner = base.proj.outcome_loader
    A = set(base.A)

    def loader(ids, genes, source):
        if not set(ids) <= A:
            calls.append(source)
        return inner(ids, genes, source)

    w = dataclasses.replace(base, proj=assoc.Project(root=root, outcome_loader=loader), log=[], _d2=None)
    aid = D.confirm_order(w.proj)[0]
    with pytest.raises(assoc.FirewallError, match="closed"):
        assoc.resume_lockbox_read(aid, proj=w.proj)
    with pytest.raises(assoc.FirewallError, match="committed"):
        assoc.resume_lockbox_read(aid, proj=w.proj, relook="R6")
    for rel, txt in ((assoc.RELOOK_README, '"R6 (PI, synthetic). Relook authorised."\n'),
                     (assoc.RELOOK_INTENT, "R6 intent (synthetic)\n")):
        (root / rel).parent.mkdir(parents=True, exist_ok=True)
        (root / rel).write_text(txt)
    confirm.commit_paths(w.proj, [root / assoc.RELOOK_README, root / assoc.RELOOK_INTENT], "R6 note + intent")
    with pytest.raises(assoc.FirewallError):
        assoc.resume_lockbox_read(aid, proj=w.proj, relook="R7")
    n_led = len(w.proj.ledger.read_text().splitlines())
    chk = D.relook(w, "R6", "check")
    res = D.check_identical(chk)
    assert res["identical"] and res["n"] > 0 and res["n_confirmed_status_differs"] == 0, res
    assert calls == [] and len(w.proj.ledger.read_text().splitlines()) == n_led
