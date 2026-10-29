"""Lockbox firewall (DESIGN v3.1 §0 F2/F5; round2_ai_failure M-C1): static checks over Part 2 code.

1. No Part 2 file (src/part2, scripts/3x, Part 2 tests other than the firewall tests) names hu.MAP or test_eval.
2. Outcome data (CRISPR gene effect, ScreenGeneEffect, DEMETER2) is named only in the loader modules
   (data.py / paths.py), the construction code (scores.py, genefilter.py, pipeline.py, scripts/31 for the Part 1
   reproduction run) and assoc.py; inside assoc.py only in ``_default_outcome_loader``, and the outcome loader is
   called only from ``load_A_effects`` and ``read_lockbox``.
3. The lockbox model list is named only in paths.py, scripts/33 (its writer) and, in assoc.py, only inside
   ``read_lockbox``. ``load_A_effects`` reads A_models.txt only.
4. Context / lockbox / registry / statistics code never touches outcome values (incl. ``universe_models()``, which
   parses the gene-effect matrix) and never names the lockbox list.
New Part 2 files fall under these rules automatically; extending an allowlist needs a reviewed edit of this test.
"""
import ast
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = sorted((ROOT / "src" / "part2").glob("*.py"))
SCRIPTS = sorted((ROOT / "scripts").glob("[34]*_*.py"))
TESTS = sorted(p for p in (ROOT / "tests").glob("test_part2_*.py") if "firewall" not in p.name)
FILES = SRC + SCRIPTS


def rel(p):
    return str(p.relative_to(ROOT))


OUTCOME = [r"load_gene_effect", r"GENE_EFFECT_CSV", r"CRISPRGeneEffect\.csv", r"ScreenGeneEffect", r"DEMETER2_CSV",
           r"D2_combined", r"\.subset\(", r"universe_models\("]
OUTCOME_ALLOWED = {"src/part2/data.py", "src/part2/paths.py", "src/part2/scores.py", "src/part2/genefilter.py",
                   "src/part2/assoc.py", "src/part2/pipeline.py"}
# FZ-C: scripts/31 lost its all-model / U options and is no longer on the allowlist; construction model sets come
# only from src/part2/runsets.py, which reads A through assoc.a_models (manifest-verified).
LOCKBOX = [r"LOCKBOX_MODELS", r"lockbox_models"]
LOCKBOX_ALLOWED = {"src/part2/paths.py", "scripts/33_make_lockbox.py", "src/part2/assoc.py"}
NO_OUTCOME = ["src/part2/contexts.py", "src/part2/lockbox.py", "src/part2/registry.py", "src/part2/stats.py",
              "scripts/32_build_contexts.py", "scripts/33_make_lockbox.py"]


def _hits(files, pats):
    out = []
    for f in files:
        txt = f.read_text()
        for p in pats:
            for m in re.finditer(p, txt):
                out.append((rel(f), txt[: m.start()].count("\n") + 1, p))
    return out


def test_no_humap_or_test_eval_anywhere_in_part2():
    pats = [r"huMAP", r"HuMAP", r"test_eval", r"Test_positive", r"Test_negative", r"Train_positive", r"Train_negative"]
    hits = [h for h in _hits(FILES + TESTS, pats) if h[0] != "src/part2/__init__.py"]
    assert not hits, hits


def test_outcome_access_only_in_allowed_modules():
    hits = [h for h in _hits(FILES, OUTCOME) if h[0] not in OUTCOME_ALLOWED]
    assert not hits, hits
    bad = _hits([ROOT / f for f in NO_OUTCOME], OUTCOME)
    assert not bad, bad


def test_lockbox_list_named_only_in_allowed_places():
    hits = [h for h in _hits(FILES, LOCKBOX) if h[0] not in LOCKBOX_ALLOWED]
    assert not hits, hits
    bad = _hits([ROOT / f for f in NO_OUTCOME if "33_" not in f], LOCKBOX)
    assert not bad, bad


def _top_level_with(tree, pred):
    """Names of top-level statements (functions/classes; None for other statements) containing a node with pred."""
    return [getattr(node, "name", None) for node in tree.body if any(pred(n) for n in ast.walk(node))]


def test_assoc_structure():
    path = ROOT / "src" / "part2" / "assoc.py"
    text = path.read_text()
    tree = ast.parse(text)

    def is_lockbox_str(n):
        return isinstance(n, ast.Constant) and isinstance(n.value, str) and "lockbox_models" in n.value

    def is_outcome(n):
        return isinstance(n, ast.Attribute) and n.attr in ("load_gene_effect", "DEMETER2_CSV", "GENE_EFFECT_CSV")

    def calls_loader(n):
        return isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute) and n.func.attr == "outcome_loader"

    assert _top_level_with(tree, is_lockbox_str) == ["read_lockbox"]
    assert _top_level_with(tree, is_outcome) == ["_default_outcome_loader"]
    assert sorted(_top_level_with(tree, calls_loader)) == ["load_A_effects", "read_lockbox"]
    fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "load_A_effects")
    body = ast.get_source_segment(text, fn)
    assert "a_models(" in body and "_model_list" not in body and "lockbox_models" not in body


def test_construction_runs_use_A_only():
    s31 = (ROOT / "scripts" / "31_score_blocks.py").read_text()
    assert "runsets.models(" in s31 and "choices=runsets.run_names()" in s31
    for bad in ("universe_models", "crispr_model_ids", "load_gene_effect", '"all"', '"U"', "lockbox_models",
                "LOCKBOX"):
        assert bad not in s31, bad
    rs = (ROOT / "src" / "part2" / "runsets.py").read_text()
    assert "a_models()" in rs
    for bad in ("universe_models", "universe_model_ids", "crispr_model_ids", "lockbox_models", "LOCKBOX",
                "load_gene_effect"):
        assert bad not in rs, bad


# ------------------------------------------------------------------------------------------------ persisted reads (R1)
# The crash-recovery copies of lockbox reads (data_cache/part2/lockbox_reads/, hash ledger results/part2/
# lockbox_reads.jsonl) are a second route to L values. Only the reader / resume code in assoc.py, the confirmation step
# that resumes (confirm.py) and the lockbox driver may name them; any other module or script under src/ or scripts/
# (verify scripts included) fails. Tests are exempt (they exercise the mechanism on synthetic repositories).
PERSISTED = [r"lockbox_reads", r"resume_lockbox_read", r"_persist_read", r"_reads_dir", r"_reads_sidecar"]
PERSISTED_ALLOWED = {"src/part2/assoc.py", "src/part2/confirm.py", "scripts/51_lockbox_confirm.py"}
ALL_CODE = sorted((ROOT / "src").rglob("*.py")) + sorted((ROOT / "scripts").rglob("*.py"))


def test_persisted_lockbox_reads_only_in_reader_resume_and_driver():
    hits = [h for h in _hits(ALL_CODE, PERSISTED) if h[0] not in PERSISTED_ALLOWED]
    assert not hits, hits
    tree = ast.parse((ROOT / "src" / "part2" / "assoc.py").read_text())
    docs = {id(n.body[0].value) for n in ast.walk(tree)
            if isinstance(n, (ast.Module, ast.FunctionDef, ast.ClassDef)) and n.body
            and isinstance(n.body[0], ast.Expr) and isinstance(n.body[0].value, ast.Constant)}

    def names_dir(n):                                          # a string naming the directory (docstrings excluded)
        return isinstance(n, ast.Constant) and isinstance(n.value, str) and "lockbox_reads" in n.value \
            and id(n) not in docs

    def calls(name):
        return lambda n: isinstance(n, ast.Call) and getattr(n.func, "id", getattr(n.func, "attr", None)) == name

    assert sorted(x for x in _top_level_with(tree, names_dir) if x) == ["_reads_dir", "_reads_sidecar"]
    assert _top_level_with(tree, calls("_persist_read")) == ["read_lockbox"]
    assert _top_level_with(tree, calls("_reads_dir")) == ["_persist_read"]
    assert sorted(_top_level_with(tree, calls("_reads_sidecar"))) == ["_persist_read", "resume_lockbox_read"]
    fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "resume_lockbox_read")
    assert not any(isinstance(n, ast.Attribute) and n.attr == "outcome_loader" for n in ast.walk(fn))
