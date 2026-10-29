"""Independent pre-read review of the lockbox machinery: adversarial synthetic checks.

Synthetic data only (tests/part2_synth_lockbox.py): no DepMap value and no lockbox model list of the real project is
touched. Every test works on its own copy of a synthetic repository in which the A side is committed and nothing has
been read. The outcome loader is wrapped so that every call for non-A (lockbox) rows is counted.

Properties that hold today are plain tests. Gaps found in review were written as the desired property and marked
``xfail(strict=True)``; the five markers were removed once the fixes made them pass (MUST-1, S1, S2).
See results/part2/verify/lockbox_code_review.md.
"""
import dataclasses
import json
import shutil
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))

import part2_synth_lockbox as S  # noqa: E402
from part2 import assoc, confirm, registry  # noqa: E402
from part2.registry import RegistryError  # noqa: E402

REFUSALS = (RegistryError, RuntimeError, assoc.FreezeError, assoc.FirewallError)


@pytest.fixture(scope="module")
def base(tmp_path_factory):
    d = tmp_path_factory.mktemp("lockbox_review")
    world, info = S.build_world(d / "repo")
    return dict(tmp=d, root=d / "repo", world=world, info=info, D=info["D"])


def fresh(base, name):
    """A copy of the untouched synthetic repo; returns (world, list of lockbox loader calls)."""
    dst = base["tmp"] / name
    shutil.copytree(base["root"], dst, symlinks=True)
    calls = []
    A = set(base["world"].A)
    inner = base["world"].proj.outcome_loader

    def loader(model_ids, genes, source):
        ids = list(model_ids)
        if not set(ids) <= A:
            calls.append(dict(source=source, n_models=len(ids), n_genes=None if genes is None else len(genes)))
        return inner(ids, genes, source)

    w = dataclasses.replace(base["world"], proj=assoc.Project(root=dst, outcome_loader=loader), log=[], _d2=None)
    return w, calls


def _committed_ledger(proj):
    return [json.loads(l) for l in assoc._committed_text(proj, proj.ledger).splitlines() if l.strip()]


# ------------------------------------------------------------------------------------------------ crash safety
def test_crash_after_read_is_fail_closed_and_never_rereads(base):
    """A crash after the L values are loaded (here: the L score raises) leaves a committed ledger line, an open entry
    and no outputs. Every later path refuses: no second loader call, the driver cannot resume, BENCH stays blocked."""
    D = base["D"]
    w, calls = fresh(base, "crash_after_read")
    aid = D.confirm_order(w.proj)[0]
    sp = D.build_spec(w, aid)

    def boom(E_L, c, u):
        raise RuntimeError("simulated crash after the read")

    sp.score_L = boom
    with pytest.raises(RuntimeError, match="simulated crash"):
        D.run_entry(w, sp)
    assert len(calls) == 1 and calls[0]["source"] == "crispr"
    assert [r["analysis_id"] for r in _committed_ledger(w.proj)] == [aid]          # detectable from git alone
    assert not confirm.is_closed(w.proj, aid)
    assert not (w.out(aid) / "confirmation.tsv").exists()
    with pytest.raises(assoc.FirewallError):
        D.run_entry(w, D.build_spec(w, aid))
    with pytest.raises(assoc.FirewallError, match="already read"):
        assoc.read_lockbox(aid, confirm.target_path(w.proj, aid), proj=w.proj)
    gaps = D.preflight(w)["gaps"]
    assert any("ledger line already exists" in g for g in gaps) and any("target list already exists" in g for g in gaps)
    with pytest.raises(D.GapError):
        D.run_all(w)
    with pytest.raises(assoc.FirewallError):
        D.bench_stage(w)
    assert len(calls) == 1                                                         # never a silent re-read


def test_ledger_commit_failure_rolls_back_and_reads_nothing(base):
    """If the ledger commit fails (pre-commit hook rejects it), the ledger is rolled back and no L value is loaded;
    once the commit can succeed, exactly one read happens."""
    D = base["D"]
    w, calls = fresh(base, "ledger_commit_fails")
    aid = D.confirm_order(w.proj)[0]
    confirm.write_target_list(w.proj, aid, ["G0001 (1)"])
    hook = w.proj.root / ".git" / "hooks" / "pre-commit"
    hook.parent.mkdir(parents=True, exist_ok=True)
    hook.write_text("#!/bin/sh\ngit diff --cached --name-only | grep -q lockbox_ledger && exit 1\nexit 0\n")
    hook.chmod(0o755)
    with pytest.raises(assoc.FreezeError, match="nothing read"):
        assoc.read_lockbox(aid, confirm.target_path(w.proj, aid), proj=w.proj)
    assert calls == [] and not w.proj.ledger.exists()
    assert registry.git(w.proj.root, "ls-files", "--", w.proj.rel(w.proj.ledger)).stdout.strip() == ""
    hook.unlink()
    assoc.read_lockbox(aid, confirm.target_path(w.proj, aid), proj=w.proj)
    assert len(calls) == 1
    with pytest.raises(assoc.FirewallError, match="already read"):
        assoc.read_lockbox(aid, confirm.target_path(w.proj, aid), proj=w.proj)
    assert len(calls) == 1


def test_process_killed_between_ledger_append_and_commit_blocks_every_read(base):
    """An appended but uncommitted ledger line (process killed before the commit) blocks that entry in the driver and
    every read of every other entry, with no loader call."""
    D = base["D"]
    w, calls = fresh(base, "killed_before_commit")
    a0, a1 = D.confirm_order(w.proj)[:2]
    w.proj.ledger.parent.mkdir(parents=True, exist_ok=True)
    w.proj.ledger.write_text(json.dumps(dict(analysis_id=a0, source="crispr")) + "\n")
    with pytest.raises(assoc.FirewallError):
        D.run_entry(w, D.build_spec(w, a0))
    with pytest.raises(assoc.FreezeError, match="uncommitted"):
        D.run_entry(w, D.build_spec(w, a1))
    assert calls == []


# ------------------------------------------------------------------------------------------------ intent integrity
def test_target_list_tampered_after_intent(base):
    """Uncommitted edits to the committed target list are refused before the ledger is written. A *committed* edit
    is not prevented, but the ledger records the blob that was read, so it differs from the intent commit's blob."""
    D = base["D"]
    w, calls = fresh(base, "target_tamper")
    aid = D.confirm_order(w.proj)[0]
    sp = D.build_spec(w, aid)
    units = {u: sp.members_of(u) for u in sorted(set(sp.hits.unit))}
    genes = sorted({g for v in units.values() for g in v})
    tp = confirm.write_target_list(w.proj, aid, genes, hits=sp.hits[["context", "unit", "beta", "alternative"]],
                                   units=units)
    intent_blob = registry.committed_blob(w.proj.root, tp)
    rec = json.loads(tp.read_text())
    rec["genes"] = rec["genes"] + ["EXTRA (999999)"]
    tp.write_text(json.dumps(rec, indent=1, sort_keys=True) + "\n")
    with pytest.raises(assoc.FreezeError, match="committed"):
        assoc.read_lockbox(aid, tp, proj=w.proj)
    assert calls == [] and not w.proj.ledger.exists()
    confirm.commit_paths(w.proj, [tp], "tamper: target list widened after intent")
    assoc.read_lockbox(aid, tp, proj=w.proj)
    led = _committed_ledger(w.proj)
    assert led[0]["targets_blob"] != intent_blob and led[0]["n_genes"] == len(genes) + 1
    assert calls[0]["n_genes"] == len(genes) + 1


def test_hits_must_match_committed_target_list(base):
    D = base["D"]
    w, calls = fresh(base, "hits_tamper")
    aid = D.confirm_order(w.proj)[0]
    sp = D.build_spec(w, aid)
    h4 = sp.hits[["context", "unit", "beta", "alternative"]]
    units = {u: sp.members_of(u) for u in sorted(set(sp.hits.unit))}
    genes = sorted({g for v in units.values() for g in v})
    confirm.write_target_list(w.proj, aid, genes, hits=h4, units=units)
    flipped = h4.assign(beta=-h4.beta)                     # hypothesis sign changed after the intent commit
    with pytest.raises(REFUSALS):
        confirm.confirm_entry(aid, flipped, sp.score_L, sp.ctx_all, sp.cov_all, proj=w.proj, n_perm=0, write=False)
    assert calls == []


@pytest.mark.parametrize("which", ["classes", "testable"])
def test_uncommitted_eligibility_or_route_inputs_refused(base, which):
    D = base["D"]
    w, calls = fresh(base, f"input_tamper_{which}")
    if which == "classes":
        aid = "H1_within_complex_GI"
        p = w.proj.results_dir / "directions" / aid / "classes.tsv"
        t = pd.read_csv(p, sep="\t", keep_default_na=False)
        t["class_complex_wide"] = "COLLAPSE"                # every N-axis call becomes 'eligible'
    else:
        aid = "H3b_driver_family"
        p = w.proj.root / S._script("51_lockbox_confirm.py").TESTABLE_REL
        t = pd.read_csv(p, sep="\t", keep_default_na=False)
        t.loc[t.context.isin(list(S.RARE)), "lockbox_testable"] = True     # rare contexts re-routed to 'lockbox'
    t.to_csv(p, sep="\t", index=False)
    with pytest.raises(REFUSALS):
        D.run_entry(w, D.build_spec(w, aid))
    assert calls == []


# ------------------------------------------------------------------------------------------------ L into A
def test_loader_returning_lockbox_rows_for_A_is_refused(base):
    sim = base["info"]["sim"]
    w, _ = fresh(base, "leak_loader")
    L_id = sorted(sim.fold.index[sim.fold == "L"])[0]
    inner = w.proj.outcome_loader

    def leaky(model_ids, genes, source):
        return inner(list(model_ids) + [L_id], genes, source)

    with pytest.raises(assoc.FirewallError, match="non-A"):
        assoc.load_A_effects(proj=assoc.Project(root=w.proj.root, outcome_loader=leaky))


def test_driver_refuses_world_with_lockbox_rows_in_A(base):
    sim = base["info"]["sim"]
    D = base["D"]
    w, calls = fresh(base, "leak_world")
    L_id = sorted(sim.fold.index[sim.fold == "L"])[0]
    w = dataclasses.replace(w, E_A=pd.concat([w.E_A, sim.E.loc[[L_id], w.E_A.columns]]))
    pf = None
    try:
        pf = D.preflight(w)
    except REFUSALS:
        return
    assert pf["gaps"], "an L model in World.E_A passed the preflight"


def test_second_source_read_of_same_entry_refused(base):
    D = base["D"]
    w, calls = fresh(base, "second_source")
    aid = D.confirm_order(w.proj)[0]
    confirm.write_target_list(w.proj, aid, ["G0001 (1)"])
    assoc.read_lockbox(aid, confirm.target_path(w.proj, aid), proj=w.proj)
    with pytest.raises(assoc.FirewallError):
        assoc.read_lockbox(aid, confirm.target_path(w.proj, aid), proj=w.proj, source="demeter2")
    assert len(calls) == 1


# ------------------------------------------------------------------------------------------------ re-check (A5, resume)
def _crash_first_control(D, w):
    aid = D.confirm_order(w.proj)[0]
    sp = D.build_spec(w, aid)

    def boom(E_L, c, u):
        raise RuntimeError("simulated crash after the read")

    sp.score_L = boom
    with pytest.raises(RuntimeError, match="simulated crash"):
        D.run_entry(w, sp)
    return aid


def test_persist_hash_commit_failure_is_fail_closed(base):
    """The read happened but its hash could not be committed (pre-commit hook rejects lockbox_reads.jsonl): no resume,
    no hand close, no second loader call."""
    D = base["D"]
    w, calls = fresh(base, "persist_commit_fails")
    hook = w.proj.root / ".git" / "hooks" / "pre-commit"
    hook.parent.mkdir(parents=True, exist_ok=True)
    hook.write_text("#!/bin/sh\ngit diff --cached --name-only | grep -q lockbox_reads && exit 1\nexit 0\n")
    hook.chmod(0o755)
    aid = D.confirm_order(w.proj)[0]
    with pytest.raises(assoc.FreezeError, match="hash could not be committed"):
        D.run_entry(w, D.build_spec(w, aid))
    hook.unlink()
    assert len(calls) == 1
    wr = dataclasses.replace(w, resume=True)
    with pytest.raises(REFUSALS):
        D.run_entry(wr, D.build_spec(wr, aid))
    with pytest.raises(REFUSALS):
        assoc.resume_lockbox_read(aid, proj=w.proj)
    with pytest.raises(RegistryError):
        confirm.close_entry(w.proj, aid)
    assert len(calls) == 1


def test_resume_refused_when_committed_inputs_change_after_the_read(base):
    """After a crash, committing a changed hypothesis input (a control's A beta) makes --resume refuse: the hits no
    longer match the committed target list. The loader is never called again."""
    D = base["D"]
    w, calls = fresh(base, "resume_inputs_changed")
    aid = _crash_first_control(D, w)
    p = w.proj.results_dir / "controls" / "controls_summary.tsv"
    t = pd.read_csv(p, sep="\t")
    m = (t.entry == aid) & t.beta.notna()
    t.loc[m, "beta"] = t.loc[m, "beta"] * 1.5
    t.to_csv(p, sep="\t", index=False)
    confirm.commit_paths(w.proj, [p], "tamper: control A beta changed after the read")
    wr = dataclasses.replace(w, resume=True)
    with pytest.raises(RegistryError, match="differ from the committed target list"):
        D.run_entry(wr, D.build_spec(wr, aid))
    assert len(calls) == 1 and not confirm.is_closed(w.proj, aid)


def test_resume_refused_after_close_and_for_unread_entries(base):
    D = base["D"]
    w, calls = fresh(base, "resume_closed")
    a0, a1 = D.confirm_order(w.proj)[:2]
    D.run_entry(w, D.build_spec(w, a0))
    assert confirm.is_closed(w.proj, a0) and len(calls) == 1
    with pytest.raises(assoc.FirewallError, match="closed"):
        assoc.resume_lockbox_read(a0, proj=w.proj)
    with pytest.raises(assoc.FirewallError, match="no committed ledger line"):
        assoc.resume_lockbox_read(a1, proj=w.proj)
    assert len(calls) == 1
