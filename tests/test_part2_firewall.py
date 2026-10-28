"""Firewall F2/F3: Part 2 construction code never reads hu.MAP files, test-evaluation results, curated complex
resources or validation-only data (BioPlex / proteomics / HGNC gene_group) for selection."""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FILES = sorted(list((ROOT / "src" / "part2").glob("*.py")) + list((ROOT / "scripts").glob("3*_*.py")))
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
    assert not hits, hits
