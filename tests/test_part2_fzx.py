"""FZ-X freeze (DESIGN v3.1 §0 F4; scripts/39_freeze_fzx.py, src/part2/fzx.py).

1. Toy checks of the FZ-X helpers: map conversion, evaluability rule, control-unit rules, registry text fill.
2. The committed freeze: manifest hashes, unit sets / complex-set hashes and m in every registry entry, H1 event list
   bookkeeping, lockbox-testable list, and that the FZ-X registry commit only filled placeholders (git history).
3. scripts/31_score_blocks.py refuses any model outside A (hard guard, not only the run-set convention).
No dependency value is read.
"""
import ast
import copy
import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import fzx, h1, registry  # noqa: E402
from part2.cmap import FrozenMap  # noqa: E402
from part2.cscore import Annotation, Removal  # noqa: E402

FZX = ROOT / "results" / "part2" / "fzx"
MAN = FZX / "FZX_manifest.json"
REG = ROOT / "docs" / "registry"


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# ------------------------------------------------------------------------------------------------ 1. toy checks
def _toy_map_files():
    cl = pd.DataFrame(dict(cluster_id=["Z1", "Z1.1", "Z1.2", "Z2"], level=[1, 2, 2, 1], parent=[np.nan, "Z1", "Z1", np.nan],
                           size=[5, 2, 2, 3], locus=[False, False, False, True]))
    mem = pd.DataFrame(dict(cluster_id=["Z1"] * 5 + ["Z1.1"] * 2 + ["Z1.2"] * 2 + ["Z2"] * 3,
                            level=[1] * 5 + [2] * 4 + [1] * 3,
                            symbol=["A", "B", "C", "D", "E", "A", "B", "C", "D", "X", "Y", "W"],
                            entrez=[1, 2, 3, 4, 5, 1, 2, 3, 4, 10, 11, 12]))
    att = pd.DataFrame(dict(cluster=["Z1", "Z1.1", "Z2", "Z1"], symbol=["G", "G", "H", "K"], entrez=[7, 7, 8, 9],
                            locus_attach=[False, False, False, True]))
    tiers = pd.DataFrame(dict(cluster_id=["Z1", "Z2"], tier=["T1", "T3-neg"]))
    return cl, mem, att, tiers


def test_engine_membership_conversion():
    cl, mem, att, tiers = _toy_map_files()
    em = fzx.engine_membership(cl, mem, att, tiers, column_of={i: f"S{i} ({i})" for i in range(1, 13)})
    fm = FrozenMap.from_membership(em)
    assert fm.coarse() == ["Z1", "Z2"] and fm.submodules("Z1") == ["Z1.1", "Z1.2"]
    assert fm.members["Z1"] == ["S1 (1)", "S2 (2)", "S3 (3)", "S4 (4)", "S5 (5)"]   # DepMap column keys by Entrez
    # attached: level-1 and not locus-flagged only (the sub-module and the locus attachment are dropped)
    assert fm.attached == {"Z1": ["S7 (7)"], "Z2": ["S8 (8)"]}
    assert fm.locus["Z2"] and fm.tier["Z1"] == "T1" and fm.parent["Z1.2"] == "Z1"
    with pytest.raises(KeyError):
        fzx.engine_membership(cl, mem, att, tiers, column_of={1: "S1 (1)"})
    bad = cl.copy()
    bad.loc[0, "size"] = 6
    with pytest.raises(ValueError):
        fzx.engine_membership(bad, mem)


def test_m_rule_counts_and_reasons():
    rng = np.random.default_rng(0)
    n = 60
    ids = [f"m{i}" for i in range(n)]
    ctx = pd.DataFrame({"c": np.r_[np.ones(14), np.zeros(46)], "c_small": np.r_[np.ones(11), np.zeros(49)]},
                       index=ids)
    ctx.iloc[0, 0] = np.nan                                         # NaN context -> excluded (13 left)
    cov = pd.DataFrame({"ScreenType3DO": rng.integers(0, 2, n).astype(float)}, index=ids)
    cov.iloc[1, 0] = np.nan                                         # incomplete covariate -> excluded (12 left)
    fm = FrozenMap(members={"U1": ["a (1)", "b (2)", "c (3)"], "U2": ["d (4)", "e (5)"]}, level={"U1": 1, "U2": 1})
    obs = pd.DataFrame(True, index=ids, columns=["a (1)", "b (2)", "c (3)", "d (4)", "e (5)"])
    obs.loc["m2", ["a (1)", "b (2)"]] = False                       # 1/3 observed in m2 -> U1 score NaN there
    entry = {"covariates": ["ScreenType3DO"], "context_family": {"contexts": [
        {"id": "c", "definition": "x"}, {"id": "c_small", "definition": "y"}]}}
    rem = {"c": Removal(), "c_small": Removal(genes=("d",), paralogs=False)}
    ev = fzx.evaluability(entry, ctx, cov, fm, ["U1", "U2"], Annotation(), obs, rem.get)
    r = ev.set_index(["context", "unit"])
    assert r.loc[("c", "U2"), "n_context"] == 12 and r.loc[("c", "U2"), "evaluable"]
    assert r.loc[("c", "U1"), "n_context"] == 11 and not r.loc[("c", "U1"), "evaluable"]   # m2 drops below 12
    assert "below 12/12" in r.loc[("c", "U1"), "reason"]
    assert r.loc[("c_small", "U2"), "reason"] == "< 2 members after removals"
    assert not r.loc[("c_small", "U1"), "evaluable"]
    s = fzx.context_summary(ev, 2).set_index("context")
    assert s.loc["c", "n_evaluable"] == 1 and not s.loc["c_small", "family_member"]


def test_parse_unit_rule_and_control_unit():
    r = fzx.parse_unit_rule("frozen coarse cluster containing SOX10 (else the one containing MITF); 'not evaluable "
                            "at complex level' if SOX10 is unclustered")
    assert r == dict(level=1, targets=["SOX10", "MITF"], partner=None)
    r2 = fzx.parse_unit_rule("frozen sub-module containing BRD9 (else BICRA); evaluable at complex level only if BRD9 "
                             "(or BICRA) and a cBAF-specific subunit (ARID1A, ARID1B or DPF2) fall in different "
                             "frozen sub-modules of the same frozen coarse complex; otherwise 'not evaluable'")
    assert r2 == dict(level=2, targets=["BRD9", "BICRA"], partner=["ARID1A", "ARID1B", "DPF2"])
    fm = FrozenMap(members={"C1": ["MITF (1)", "T2 (2)", "T3 (3)"], "C2": ["BRD9 (4)", "ARID1A (5)", "Q (6)", "R (7)"],
                            "C2.1": ["BRD9 (4)", "Q (6)"], "C2.2": ["ARID1A (5)", "R (7)"]},
                   level={"C1": 1, "C2": 1, "C2.1": 2, "C2.2": 2}, parent={"C1": "", "C2": "", "C2.1": "C2",
                                                                          "C2.2": "C2"})
    e1 = {"id": "CTRL_X", "score": "S_C_mean_winsor3", "control": {
        "design_id": "PCx", "complex_unit_rule": "frozen coarse cluster containing SOX10 (else the one containing "
                                                 "MITF); 'not evaluable at complex level' if ..."}}
    out = fzx.control_unit(e1, fm, Annotation(), Removal())
    assert out["unit"] == "C1" and out["target_used"] == "MITF" and out["structural_ok"]
    out = fzx.control_unit(e1, fm, Annotation(), Removal(genes=("T2", "T3"), paralogs=False))
    assert out["unit"] == fzx.NOT_EVALUABLE and "< 2 members" in out["reason"]
    e2 = dict(e1, score="submodule_vs_sibling", control=dict(e1["control"], complex_unit_rule=(
        "frozen sub-module containing BRD9 (else BICRA); evaluable at complex level only if BRD9 (or BICRA) and a "
        "cBAF-specific subunit (ARID1A, ARID1B or DPF2) fall in different frozen sub-modules of the same frozen "
        "coarse complex")))
    out = fzx.control_unit(e2, fm, Annotation(), Removal())
    assert out["unit"] == "C2.1" and out["partner_units"] == "C2.2" and out["structural_ok"]
    fm2 = copy.deepcopy(fm)
    fm2.members["C2.2"] = ["R (7)", "Q2 (8)"]
    out = fzx.control_unit(e2, fm2, Annotation(), Removal())
    assert out["unit"] == fzx.NOT_EVALUABLE and "different sub-module" in out["reason"]
    an = Annotation(paralogs={"T2": {"LOST"}, "LOST": {"T2"}})
    e3 = dict(e1, score="member_resolved", control=dict(e1["control"], complex_unit_rule=(
        "frozen coarse cluster containing T2; 'not evaluable at complex level' if T2 is unclustered")))
    out = fzx.control_unit(e3, fm, an, Removal(), dict(lost=("LOST",), target="T2"))
    assert out["structural_ok"] and out["readout_members"] == "T2 (2)"
    out = fzx.control_unit(e3, fm, Annotation(), Removal(), dict(lost=("LOST",), target="T2"))
    assert out["unit"] == fzx.NOT_EVALUABLE and "lacks T2" in out["reason"]


def test_fill_text_keeps_everything_else():
    txt = (REG / "CTRL_PC1.yaml").read_text()
    pend = _committed_pending_version("CTRL_PC1")
    if pend is None:
        pytest.skip("no git history for the FZ-R version")
    old = registry.load_yaml(pend)
    fills = {"complex_set_hash": "0123456789abcdef", "control.complex_unit": "Z0001", "family_m.value": 2,
             "frozen_inputs[2]": "a/b.json", "frozen_inputs[3]": "c/d.json"}
    new_txt = fzx.fill_text(pend, fills)
    new = registry.load_yaml(new_txt)
    assert registry.validate(new, "CTRL_PC1") == [] and registry.placeholders(new) == []
    assert registry.fill_only_changes(old, new) == [] and "pending" not in new
    assert new["control"]["complex_unit"] == "Z0001" and new["family_m"]["value"] == 2
    assert [l for l in new_txt.splitlines() if l.startswith("#")] == [l for l in pend.splitlines() if l.startswith("#")]
    with pytest.raises(ValueError):
        fzx.fill_text(pend, {"no_such_field": 1})
    assert txt                                                       # filled file exists


# ------------------------------------------------------------------------------------------------ 2. the freeze
@pytest.fixture(scope="module")
def man():
    return json.loads(MAN.read_text())


@pytest.fixture(scope="module")
def entries():
    return registry.load_all(REG)


def test_manifest_hashes(man):
    assert man["freeze"] == "FZ-X"
    for p, h in man["files"].items():
        f = ROOT / p
        if p.startswith("data_cache/"):
            if not f.exists():
                continue                                             # gitignored cache (clean clone)
        assert sha(f) == h, p
    listed = {p for p in man["files"] if p.startswith("results/part2/fzx/")}
    on_disk = {f"results/part2/fzx/{p.name}" for p in FZX.iterdir() if p.is_file()} - {
        "results/part2/fzx/FZX_manifest.json", "results/part2/fzx/README.md"}
    assert listed == on_disk
    fzc = json.loads((ROOT / "results/part2/construction/FZC_manifest.json").read_text())["sha256"]
    for p in ("results/part2/construction/zavg_map/membership.tsv", "results/part2/construction/zavg_map/clusters.tsv",
              "results/part2/construction/gls_map/membership.tsv"):
        assert man["files"][p] == fzc[p], p                          # the FZ-C map, unchanged


def test_unit_sets_hashes_and_m(man, entries):
    z = (FZX / "units_zavg_coarse.txt").read_text().split()
    g = (FZX / "units_gls_coarse.txt").read_text().split()
    sub = list(pd.read_csv(FZX / "h2_submodules.tsv", sep="\t").unit)
    h1u = (FZX / "units_h1_components_A1.txt").read_text().split()
    h1s = (FZX / "units_h1s_components.txt").read_text().split()
    cl = pd.read_csv(ROOT / "results/part2/construction/zavg_map/clusters.tsv", sep="\t")
    assert z == sorted(cl.cluster_id[cl.level == 1]) and len(z) == 367
    em = FrozenMap.from_membership(pd.read_csv(FZX / "zavg_engine_membership.tsv", sep="\t"))
    from part2 import h2
    assert sub == h2.h2_units(em)
    for k in ("H3a_lineage_atlas", "H3b_driver_family", "BENCH_complex_vs_gene"):
        assert entries[k]["complex_set_hash"] == registry.set_hash(z) == man["complex_set_hashes"]["zavg_coarse"], k
    assert entries["BRIDGE-L_confirmation_rate"]["complex_set_hash"] == registry.set_hash(g)
    assert entries["H2_submodule_divergence"]["complex_set_hash"] == registry.set_hash(sub)
    assert entries["H1_within_complex_GI"]["complex_set_hash"] == registry.set_hash(h1u) == \
        man["complex_set_hashes"]["h1_components"]
    assert entries["H1S_silencing_GI"]["complex_set_hash"] == registry.set_hash(h1s) == \
        man["complex_set_hashes"]["h1s_components"]
    assert entries["H1_within_complex_GI"]["family_m"]["value"] == len(h1u) == man["summary"]["m"]["H1"]["m"]
    assert entries["H1S_silencing_GI"]["family_m"]["value"] == len(h1s) == man["summary"]["m"]["H1S"]["m"]
    ev = pd.read_csv(FZX / "evaluability_tests.tsv.gz", sep="\t")
    for k in ("H3a_lineage_atlas", "H3b_driver_family"):
        e = ev[ev.entry == k]
        prim = [c["id"] for c in entries[k]["context_family"]["contexts"] if c.get("role", "primary") == "primary"]
        assert sorted(set(e.context)) == sorted(prim) and len(e) == len(prim) * len(z), k
        assert entries[k]["family_m"]["value"] == int(e.evaluable.sum()), k
        assert entries[k]["family_m"]["n_units"] == len(z), k
    b = ev[ev.entry == "BRIDGE-L_confirmation_rate"]
    assert entries["BRIDGE-L_confirmation_rate"]["family_m"]["value"]["total"] == int(b.evaluable.sum())
    for e in entries.values():
        assert registry.placeholders(e) == [] and "pending" not in e, e["id"]


def test_h1_event_list(man, entries):
    ev = pd.read_csv(FZX / "h1_events.tsv", sep="\t", keep_default_na=False)
    ctx = pd.read_csv(FZX / "h1_contexts.tsv", sep="\t")
    assert list(ev["index"]) == list(range(len(ev)))
    assert list(ev.event) == sorted(ev.event, key=lambda x: (x.split("|")[0], x.split("|")[1]))
    for r in ev.itertuples(index=False):
        assert r.event == f"{r.unit}|{r.gene}"
        assert r.n_A >= 12 and r.ref_A >= 12, r.event                       # rule 3
        assert r.n_N >= 2 or r.n_P >= 1, r.event                            # rule 5
        assert r.testable == (r.n_L >= 5 and r.ref_L >= 5), r.event          # rule 7
        assert r.type in ("LoF", "CN", "silencing")
        assert r.covariates.split(";") == h1.COVARIATES[r.type]
    units = sorted([f"{r.event}|N" for r in ev.itertuples() if r.n_N >= 2] +
                   [f"{r.event}|P" for r in ev.itertuples() if r.n_P >= 1])
    assert units == (FZX / "units_h1_components.txt").read_text().split()       # the FZ-X list (366)
    # amendment A1: LoF / CN events -> H1, silencing events -> H1S; together exactly the FZ-X list
    split = pd.read_csv(FZX / "h1_event_split_A1.tsv", sep="\t")
    assert list(split.event) == list(ev.event)
    assert (split.entry == np.where(ev.type.isin(["LoF", "CN"]), "H1_within_complex_GI", "H1S_silencing_GI")).all()
    u1 = (FZX / "units_h1_components_A1.txt").read_text().split()
    u2 = (FZX / "units_h1s_components.txt").read_text().split()
    assert sorted(u1 + u2) == units and not set(u1) & set(u2)
    assert entries["H1_within_complex_GI"]["family_m"]["value"] + entries["H1S_silencing_GI"]["family_m"]["value"] \
        == len(units)
    assert man["summary"]["h1"]["n_events"] == len(ev)
    assert sorted(ctx.context) == sorted(f"H1loss:{g}" for g in set(ev.gene))
    # events link only to coarse units; attachments only through non-locus level-1 attachments
    em = pd.read_csv(FZX / "zavg_engine_membership.tsv", sep="\t")
    lev = dict(zip(em.unit_id, em.level))
    assert all(lev[u] == 1 for u in ev.unit)
    att = set(zip(em.unit_id[em.role == "attached"], em.gene[em.role == "attached"].str.split(" \\(").str[0]))
    for r in ev[ev.link == "attached"].itertuples():
        assert (r.unit, r.gene) in att, r.event


def test_lockbox_testable_list(entries):
    lb = pd.read_csv(FZX / "lockbox_testable.tsv", sep="\t", keep_default_na=False)
    assert set(lb.entry) == set(entries)
    for r in lb.itertuples(index=False):
        assert str(r.lockbox_testable) == str(int(r.n_L) >= 5 and int(r.ref_L) >= 5), (r.entry, r.context)
    split = pd.read_csv(FZX / "h1_event_split_A1.tsv", sep="\t")
    for k in ("H1_within_complex_GI", "H1S_silencing_GI"):
        assert sorted(lb.event[lb.entry == k]) == sorted(split.event[split.entry == k]), k


def _committed_pending_version(eid: str):
    """Text of the registry entry in the last commit that still had `pending:` (the FZ-R draft), or None."""
    rel = f"docs/registry/{eid}.yaml"
    r = subprocess.run(["git", "-C", str(ROOT), "log", "--format=%H", "--", rel], capture_output=True, text=True)
    if r.returncode != 0:
        return None
    for c in r.stdout.split():
        t = subprocess.run(["git", "-C", str(ROOT), "show", f"{c}:{rel}"], capture_output=True, text=True)
        if t.returncode == 0 and "\npending:" in t.stdout:
            return t.stdout
    return None


def _show(commit: str, rel: str):
    t = subprocess.run(["git", "-C", str(ROOT), "show", f"{commit}:{rel}"], capture_output=True, text=True)
    return t.stdout if t.returncode == 0 else None


def _fzx_commit(man):
    return man["amendments"][0]["amends_commit"]


def test_fzx_fill_changed_only_placeholders(entries, man):
    """FZ-R draft -> FZ-X commit: fill-only for every entry (placeholders and `pending` only)."""
    n = 0
    c = _fzx_commit(man)
    for k in entries:
        old_txt = _committed_pending_version(k)
        fzx_txt = _show(c, f"docs/registry/{k}.yaml")
        if old_txt is None or fzx_txt is None:
            continue
        bad = registry.fill_only_changes(registry.load_yaml(old_txt), registry.load_yaml(fzx_txt))
        assert bad == [], (k, bad)
        n += 1
    if n == 0:
        pytest.skip("no history with the FZ-R drafts")


def test_amendment_A1_is_the_only_non_fill_change(entries, man):
    """FZ-X commit -> now: A1 is the only change; H1 and H1S differ from the FZ-X H1 entry only in A1's declared
    paths (all of which changed); every other entry is identical to its FZ-X version (status may close later)."""
    c = _fzx_commit(man)
    h1_fzx_txt = _show(c, "docs/registry/H1_within_complex_GI.yaml")
    if h1_fzx_txt is None:
        pytest.skip("no history with the FZ-X commit")
    h1_fzx = registry.load_yaml(h1_fzx_txt)
    assert "amendments" not in h1_fzx and len(man["amendments"]) == 1 and man["amendments"][0]["id"] == "A1"
    for k, e in entries.items():
        am = e.get("amendments", [])
        if k in ("H1_within_complex_GI", "H1S_silencing_GI"):
            assert [a["id"] for a in am] == ["A1"] and am[0]["amends_commit"] == c, k
            assert registry.amendment_changes(h1_fzx, e) == [], k
            declared = am[0]["changed_paths"]
            assert registry.changed_top_paths(h1_fzx, e, declared) == set(declared), k
        else:
            assert am == [], k
            old = registry.load_yaml(_show(c, f"docs/registry/{k}.yaml"))
            assert [p for p in registry.amendment_changes(old, e) if p != "status"] == [], k
    assert set(entries["H1_within_complex_GI"]["amendments"][0]["changed_paths"]) == {
        "context_family.contexts", "complex_set_hash", "family_m.value"}
    txt = (REG / "AMENDMENTS.md").read_text()
    assert "## A1" in txt and c in txt and man["amendments"][0]["timestamp_utc"] in txt


def test_amendment_machinery_toy():
    base = {"id": "X", "a": 1, "b": {"c": [1, 2]}, "status": "open"}
    new = copy.deepcopy(base)
    new["a"] = 2
    new["b"]["c"] = [1, 3]
    assert registry.amendment_changes(base, new) == ["a", "b.c[1]"]                 # undeclared -> reported
    new["amendments"] = [dict(id="A9", timestamp_utc="2026-01-01T00:00:00Z", amends_commit="0" * 40,
                              reason="r", changed_paths=["a"])]
    assert registry.amendment_changes(base, new) == ["b.c[1]"]
    new["amendments"][0]["changed_paths"] = ["a", "b.c"]
    assert registry.amendment_changes(base, new) == []
    assert registry.amendment_changes(new, new) == []
    assert registry.fill_only_changes(base, new) != []                               # fills stay strict
    bad = {"amendments": [dict(id="A9", timestamp_utc="yesterday", amends_commit="abc", reason="", changed_paths=[])]}
    errs = registry._validate_amendments(bad)
    assert any("timestamp_utc" in x for x in errs) and any("amends_commit" in x for x in errs)
    assert any("reason" in x for x in errs) and any("changed_paths" in x for x in errs)


# ------------------------------------------------------------------------------------------------ 3. scripts/31 guard
def _load_31():
    spec = importlib.util.spec_from_file_location("score_blocks_31", ROOT / "scripts" / "31_score_blocks.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_score_blocks_guard_is_called_before_scoring():
    tree = ast.parse((ROOT / "scripts" / "31_score_blocks.py").read_text())
    main = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "main")
    calls = [(n.lineno, n.func.attr if isinstance(n.func, ast.Attribute) else getattr(n.func, "id", ""))
             for n in ast.walk(main) if isinstance(n, ast.Call)]
    guard = min(l for l, f in calls if f == "assert_A_only")
    first_use = min(l for l, f in calls if f in ("score_block", "S_mask"))
    assert guard < first_use
    # the guard wraps the run-set model list itself
    src = ast.get_source_segment((ROOT / "scripts" / "31_score_blocks.py").read_text(), main)
    assert "ids = assert_A_only(runsets.models(a.set))" in src


def test_score_blocks_refuses_non_A_models(monkeypatch):
    if not (ROOT / "results/part2/lockbox/A_models.txt").exists():
        pytest.skip("no A list")
    mod = _load_31()
    from part2 import assoc
    A = assoc.a_models()
    assert mod.assert_A_only(A[:5]) == A[:5]
    with pytest.raises(SystemExit, match="outside the construction set A"):
        mod.assert_A_only(A[:5] + ["ACH-999999"])
    with pytest.raises(SystemExit):
        mod.assert_A_only([])
    # end to end: even if the run-set code returned a non-A model, main() stops before any scoring
    monkeypatch.setattr(mod.consensus, "ensure_hashseed", lambda: None)
    monkeypatch.setattr(mod.part1, "assert_locked", lambda: {"commit": "test"})
    monkeypatch.setattr(mod.runsets, "models", lambda name: list(A[:3]) + ["ACH-999999"])

    def boom(*a, **k):
        raise AssertionError("scoring reached with a non-A model")

    monkeypatch.setattr(mod.scores, "score_block", boom)
    monkeypatch.setattr(mod.runsets, "S_mask", boom)
    monkeypatch.setattr(sys, "argv", ["31_score_blocks.py", "--set", "A"])
    with pytest.raises(mod.NotAModelsError):
        mod.main()
