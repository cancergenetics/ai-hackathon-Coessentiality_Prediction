"""Firewall F2/F3: Part 2 construction code never reads hu.MAP files, test-evaluation results, curated complex
resources or validation-only data (BioPlex / proteomics / HGNC gene_group) for selection."""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FILES = sorted(list((ROOT / "src" / "part2").glob("*.py")) + list((ROOT / "scripts").glob("[3-6]*_*.py")))
FORBIDDEN = [r"huMAP", r"HuMAP", r"test_eval", r"Train_positive", r"Train_negative", r"Test_positive",
             r"Test_negative", r"gene_group", r"BioPlex", r"Proteomics", r"protein_quant", r"STRING", r"CORUM",
             r"ComplexPortal", r"testeval"]


def test_no_forbidden_inputs():
    assert FILES
    hits = []
    for f in FILES:
        txt = f.read_text()
        for pat in FORBIDDEN:
            for m in re.finditer(pat, txt):
                line = txt[: m.start()].count("\n") + 1
                hits.append(f"{f.relative_to(ROOT)}:{line}:{pat}")
    # mentions inside docstrings that describe the firewall itself are allowed only in __init__.py
    hits = [h for h in hits if not h.startswith("src/part2/__init__.py")]
    # CCLE proteomics is a context input for the pre-registered LoF functional check (DESIGN v3.1 §2.2) and a later
    # validation input; it may be named only where its path is defined and in the context builder.
    # (contexts.py reaches it only through paths.CCLE_PROTEOMICS_XLSX)
    allowed = {("src/part2/paths.py", "Proteomics"), ("src/part2/paths.py", "protein_quant")}
    # FZ-C validation (tiers + bridge) is the only Part 2 code that may name BioPlex; it runs only after the map
    # commit (validation.require_frozen) - see test_construction_never_touches_validation_inputs.
    allowed |= {(f, "BioPlex") for f in VALIDATION}
    # post hoc annotation column of the descriptive summary (after tiers; labelled POSTHOC; never used upstream)
    allowed |= {("scripts/38_describe_maps.py", "gene_group")}
    hits = [h for h in hits if (h.split(":")[0], h.rsplit(":", 1)[1]) not in allowed]
    assert not hits, hits


VALIDATION = {"src/part2/validation.py", "scripts/37_validate_maps.py"}
CONSTRUCTION = ["src/part2/consensus.py", "src/part2/graph.py", "src/part2/selection.py", "src/part2/pipeline.py",
                "src/part2/construct.py", "src/part2/attach.py", "src/part2/runsets.py", "src/part2/genefilter.py",
                "src/part2/scores.py", "src/part2/locus.py", "src/part2/toy.py", "scripts/31_score_blocks.py",
                "scripts/34_run_A_blocks.sh", "scripts/35_toy_decisions.py", "scripts/36_build_maps.py"]


def test_construction_never_touches_validation_inputs():
    """FZ-C: maps are built from score blocks only; BioPlex / proteomics are read by validation code only."""
    pats = [r"BioPlex", r"[Pp]roteomic", r"PROTEOMICS", r"protein_quant", r"read_proteomics", r"validation",
            r"BIOPLEX"]
    hits = []
    for f in CONSTRUCTION:
        txt = (ROOT / f).read_text()
        for pat in pats:
            for m in re.finditer(pat, txt):
                hits.append(f"{f}:{txt[: m.start()].count(chr(10)) + 1}:{pat}")
    assert not hits, hits


def test_validation_requires_frozen_maps():
    txt = (ROOT / "scripts" / "37_validate_maps.py").read_text()
    main = txt[txt.index("def main("):]
    first_read = min(i for i in (main.find("load_bioplex("), main.find("read_proteomics(")) if i >= 0)
    assert 0 <= main.find("require_frozen(") < first_read
