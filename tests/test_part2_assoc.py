"""Registry + gated association runner + lockbox reader, exercised in a throw-away git repository with toy data:
refusals (unregistered / uncommitted / modified / closed entries, unfrozen inputs, non-A models, unregistered
contexts, wrong unit set), the trials ledger (counts only), a planted effect, and the lockbox ledger rules
(committed target list, one read per entry, confirmatory-before-benchmark order, ledger restored on commit failure)."""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import assoc, registry  # noqa: E402
from part2.assoc import FirewallError, FreezeError  # noqa: E402
from part2.registry import RegistryError  # noqa: E402

UNITS = [f"C{i}" for i in range(20)]
MODELS = [f"ACH-{i:06d}" for i in range(400)]
L_MODELS = MODELS[300:]
A_MODELS = MODELS[:300]
GENES = [f"G{i} ({i})" for i in range(30)]


def entry_yaml(aid, direction="H3b", status="open", contexts=("ctxA", "ctxB"), units=UNITS):
    ctx = "\n".join(f"    - id: {c}\n      definition: \"toy context {c}\"" for c in contexts)
    return f"""# toy registry entry
id: {aid}
direction: {direction}
hypothesis: |
  Toy hypothesis for {aid}.
  Second line.
sign: "+"
context_family:
  name: toy
  contexts:
{ctx}
complex_set_hash: {registry.set_hash(units)}
score: S_C_mean_winsor3
covariates: [lineage_FE, library]
fdr_family:
  method: benjamini_bogomolov
  families: contexts
  q: 0.05
success_criterion: "q < 0.05 in A and one-sided p < 0.05 in L"
stopping_rule: "single run; lockbox read once"
frozen_inputs:
  - results/part2/contexts/checksums.json
status: {status}
"""


def sh(root, *args):
    return subprocess.run(["git", "-C", str(root), *args], capture_output=True, text=True, check=True)


def commit_all(root, msg="x"):
    sh(root, "add", "-A")
    sh(root, "commit", "-q", "-m", msg)


def toy_loader(model_ids, genes, source):
    rng = np.random.default_rng(0)
    df = pd.DataFrame(rng.normal(size=(len(MODELS), len(GENES))), index=MODELS, columns=GENES)
    return df.reindex([m for m in model_ids]).reindex(columns=list(genes) if genes is not None else GENES)


@pytest.fixture
def proj(tmp_path):
    root = tmp_path / "repo"
    root.mkdir()
    sh(root, "init", "-q")
    sh(root, "config", "user.email", "t@example.org")
    sh(root, "config", "user.name", "toy")
    lb = root / "results" / "part2" / "lockbox"
    lb.mkdir(parents=True)
    tA, tL = "\n".join(A_MODELS) + "\n", "\n".join(L_MODELS) + "\n"
    (lb / "A_models.txt").write_text(tA)
    (lb / "lockbox_models.txt").write_text(tL)
    (lb / "manifest.json").write_text(json.dumps(dict(files={"A_models.txt": assoc._sha_text(tA),
                                                             "lockbox_models.txt": assoc._sha_text(tL)})))
    cdir = root / "results" / "part2" / "contexts"
    cdir.mkdir(parents=True)
    inp = root / "data_cache" / "ctx.txt"
    inp.parent.mkdir()
    inp.write_text("frozen context table\n")
    (cdir / "checksums.json").write_text(json.dumps(dict(files={"data_cache/ctx.txt": assoc._sha_file(inp)})))
    (root / "docs" / "registry").mkdir(parents=True)
    commit_all(root, "init")
    return assoc.Project(root=root, outcome_loader=toy_loader)


def add_entry(proj, aid, commit=True, **kw):
    p = proj.registry_dir / f"{aid}.yaml"
    p.write_text(entry_yaml(aid, **kw))
    if commit:
        commit_all(proj.root, f"register {aid}")
    return p


def toy_data(effect=-1.0, seed=1):
    rng = np.random.default_rng(seed)
    n = len(A_MODELS)
    C = pd.DataFrame({"ctxA": (rng.random(n) < 0.15).astype(float), "ctxB": (rng.random(n) < 0.2).astype(float)},
                     index=A_MODELS)
    C.iloc[:5, 1] = np.nan
    Y = pd.DataFrame(rng.normal(size=(n, len(UNITS))), index=A_MODELS, columns=UNITS)
    Y["C0"] += effect * C["ctxA"].values                      # planted: C0 more dependent in ctxA
    Y.iloc[10:13, 4] = np.nan
    Z = pd.DataFrame({"cov1": rng.normal(size=n)}, index=A_MODELS)
    return Y, C, Z


def test_yaml_subset_parser():
    e = registry.load_yaml(entry_yaml("T1"))
    assert registry.validate(e, "T1") == []
    assert e["hypothesis"] == "Toy hypothesis for T1.\nSecond line.\n"
    assert e["covariates"] == ["lineage_FE", "library"] and e["fdr_family"]["q"] == 0.05
    assert e["context_family"]["contexts"][1] == {"id": "ctxB", "definition": "toy context ctxB"}
    with pytest.raises(RegistryError):
        registry.load_yaml("a: {b: 1}\n")
    with pytest.raises(RegistryError):
        registry.load_yaml("a: 1\na: 2\n")
    bad = registry.load_yaml(entry_yaml("T1").replace("direction: H3b", "direction: H9"))
    assert any("direction" in x for x in registry.validate(bad, "T1"))


def test_registry_readme_example_validates():
    txt = (ROOT / "docs" / "registry" / "README.md").read_text()
    block = txt.split("```yaml\n", 1)[1].split("```", 1)[0]
    e = registry.load_yaml(block)
    assert registry.validate(e, e["id"]) == []


def test_run_association_refusals(proj):
    Y, C, Z = toy_data()
    inputs = ["data_cache/ctx.txt"]
    with pytest.raises(RegistryError, match="not registered"):
        assoc.run_association("H3b_toy", Y, C, Z, inputs=inputs, proj=proj)
    add_entry(proj, "H3b_toy", commit=False)
    with pytest.raises(RegistryError, match="not tracked"):
        assoc.run_association("H3b_toy", Y, C, Z, inputs=inputs, proj=proj)
    commit_all(proj.root)
    p = proj.registry_dir / "H3b_toy.yaml"
    p.write_text(p.read_text().replace("Second line.", "Edited line."))
    with pytest.raises(RegistryError, match="uncommitted"):
        assoc.run_association("H3b_toy", Y, C, Z, inputs=inputs, proj=proj)
    sh(proj.root, "checkout", "--", "docs/registry/H3b_toy.yaml")
    # unfrozen input: not in a manifest / changed after freeze
    (proj.root / "data_cache" / "other.txt").write_text("x")
    with pytest.raises(FreezeError, match="not listed"):
        assoc.run_association("H3b_toy", Y, C, Z, inputs=["data_cache/other.txt"], proj=proj)
    (proj.root / "data_cache" / "ctx.txt").write_text("edited\n")
    with pytest.raises(FreezeError, match="changed since freeze"):
        assoc.run_association("H3b_toy", Y, C, Z, inputs=inputs, proj=proj)
    (proj.root / "data_cache" / "ctx.txt").write_text("frozen context table\n")
    # lockbox model in the data
    Y2 = pd.concat([Y, Y.iloc[:1].rename(index={A_MODELS[0]: L_MODELS[0]})])
    with pytest.raises(FirewallError, match="outside A"):
        assoc.run_association("H3b_toy", Y2, C, Z, inputs=inputs, proj=proj)
    with pytest.raises(RegistryError, match="not registered"):
        assoc.run_association("H3b_toy", Y, C.assign(ctxZ=0.0), Z, inputs=inputs, proj=proj)
    with pytest.raises(RegistryError, match="complex_set_hash"):
        assoc.run_association("H3b_toy", Y.iloc[:, :-1], C, Z, inputs=inputs, proj=proj)
    assert not proj.trials.exists()                          # nothing was run
    add_entry(proj, "H3b_closed", status="closed")
    with pytest.raises(RegistryError, match="closed"):
        assoc.run_association("H3b_closed", Y, C, Z, inputs=inputs, proj=proj)


def test_run_association_known_answer_and_trials(proj):
    add_entry(proj, "H3b_toy")
    Y, C, Z = toy_data(effect=-1.0)
    res = assoc.run_association("H3b_toy", Y, C, Z, inputs=["data_cache/ctx.txt"], proj=proj, n_perm=999,
                                blocks=pd.Series(A_MODELS, index=A_MODELS),
                                strata=pd.Series(np.repeat(["x", "y"], 150), index=A_MODELS))
    assert len(res) == 2 * len(UNITS)
    hit = res[(res.context == "ctxA") & (res.unit == "C0")].iloc[0]
    assert hit.reject and hit.beta < 0 and hit.alternative == "less" and hit.p_perm < 0.01
    assert res.reject.sum() == 1                               # only the planted effect
    # complete-case bookkeeping
    assert res[(res.context == "ctxB") & (res.unit == "C4")].n.iloc[0] == 300 - 5 - 3
    # HC3 df/statistic equal a direct fit
    from part2 import stats
    m = C["ctxA"].notna()
    X = np.column_stack([np.ones(m.sum()), C["ctxA"][m], Z["cov1"][m]])
    r = stats.ols_hc3(Y["C0"][m].values, X, 1, "less")
    assert np.isclose(r["t"], hit.t) and np.isclose(r["p"], hit.p)
    recs = [json.loads(l) for l in proj.trials.read_text().splitlines()]
    assert len(recs) == 1 and recs[0]["n_tests"] == 40 and recs[0]["n_reject"] == 1
    assert recs[0]["registry_blob"] and recs[0]["inputs"]["data_cache/ctx.txt"]
    assert not any(u in proj.trials.read_text() for u in ('"C0"', "ctxA\""))  # no hits in the ledger


def write_targets(proj, aid, commit=True):
    p = proj.root / "results" / "part2" / "lockbox_targets" / f"{aid}.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(dict(analysis_id=aid, genes=GENES[:5])))
    if commit:
        commit_all(proj.root, f"targets {aid}")
    return p


def test_read_lockbox_rules(proj):
    add_entry(proj, "H3b_toy")
    t = write_targets(proj, "H3b_toy", commit=False)
    with pytest.raises(FreezeError, match="committed"):
        assoc.read_lockbox("H3b_toy", t, proj=proj)
    commit_all(proj.root)
    df = assoc.read_lockbox("H3b_toy", t, proj=proj)
    assert sorted(df.index) == L_MODELS and list(df.columns) == GENES[:5]
    led = [json.loads(l) for l in proj.ledger.read_text().splitlines()]
    assert len(led) == 1 and led[0]["analysis_id"] == "H3b_toy"
    assert "lockbox read: H3b_toy" in sh(proj.root, "log", "-1", "--format=%s").stdout
    assert not sh(proj.root, "status", "--porcelain").stdout.strip()   # ledger committed
    with pytest.raises(FirewallError, match="one read"):
        assoc.read_lockbox("H3b_toy", t, proj=proj)
    # benchmark only after every confirmatory entry has read
    add_entry(proj, "H1_toy", direction="H1")
    add_entry(proj, "bench_toy", direction="benchmark")
    tb = write_targets(proj, "bench_toy")
    with pytest.raises(FirewallError, match="before all confirmatory"):
        assoc.read_lockbox("bench_toy", tb, proj=proj)
    assoc.read_lockbox("H1_toy", write_targets(proj, "H1_toy"), proj=proj)
    with pytest.raises(FirewallError, match="closed"):           # ledger lines exist but entries are still open
        assoc.read_lockbox("bench_toy", tb, proj=proj)
    from part2 import confirm
    confirm.close_entry(proj, "H3b_toy")
    confirm.close_entry(proj, "H1_toy")
    assert confirm.is_closed(proj, "H1_toy")
    assert len(assoc.read_lockbox("bench_toy", tb, proj=proj)) == len(L_MODELS)


def test_read_lockbox_restores_ledger_on_commit_failure(proj):
    add_entry(proj, "H3b_toy")
    t = write_targets(proj, "H3b_toy")
    lock = Path(sh(proj.root, "rev-parse", "--git-dir").stdout.strip())
    lock = (lock if lock.is_absolute() else proj.root / lock) / "index.lock"
    lock.write_text("")
    try:
        with pytest.raises(FreezeError, match="nothing read"):
            assoc.read_lockbox("H3b_toy", t, proj=proj)
    finally:
        lock.unlink()
    assert not proj.ledger.exists()
    assert not sh(proj.root, "status", "--porcelain").stdout.strip()


def test_load_A_effects_is_A_only(proj):
    df = assoc.load_A_effects(genes=GENES[:3], proj=proj)
    assert sorted(df.index) == A_MODELS
