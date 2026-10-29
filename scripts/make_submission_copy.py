#!/usr/bin/env python3
"""Make the clean submission copy of the manuscript and report word counts.

Reads manuscript/main.md and manuscript/supplementary.md, strips every HTML comment (the source comments that
scripts/check_manuscript_numbers.py uses), normalises the blank lines left behind, rewrites figure paths so that
they resolve from manuscript/submission/, and writes manuscript/submission/{main,supplementary}.md. The source files
stay the reference; the copy is regenerated, never edited by hand.

Word counts (printed, and written to manuscript/submission/word_counts.json) follow Nature Communications
conventions: the abstract; the main text = Introduction + Results + Discussion (headings excluded; figure legends,
Methods, references and tables excluded); Methods separately. A word is a whitespace-separated token that contains
a letter or a digit; Markdown table rows are excluded from the main-text count.

Usage: python scripts/make_submission_copy.py
"""
from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "manuscript"
OUT = SRC / "submission"
COMMENT = re.compile(r"<!--.*?-->", re.S)


def strip(md: str) -> str:
    md = COMMENT.sub("", md)
    md = re.sub(r"[ \t]+\n", "\n", md)          # trailing spaces left by removed inline comments
    md = re.sub(r" {2,}", " ", md)
    md = re.sub(r"\n{3,}", "\n\n", md)
    return md.strip() + "\n"


def words(text: str) -> int:
    lines = [ln for ln in text.splitlines() if not ln.lstrip().startswith(("#", "|"))]
    return sum(1 for tok in " ".join(lines).split() if re.search(r"[A-Za-z0-9]", tok))


def section(md: str, start: str, ends: tuple[str, ...]) -> str:
    i = md.find(start)
    if i < 0:
        return ""
    j = min([k for k in (md.find(e, i + len(start)) for e in ends) if k > 0] or [len(md)])
    return md[i + len(start):j]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    counts = {}
    for name in ("main.md", "supplementary.md"):
        clean = strip((SRC / name).read_text())
        clean = clean.replace("](figures/", "](../figures/").replace("`manuscript/figures/", "`../figures/")
        (OUT / name).write_text(clean)
        assert "<!--" not in clean and "-->" not in clean, name
        if name == "main.md":
            figs = sorted((SRC / "figures").glob("fig*.pdf"))
            listing = "Figure files (PDF; PNG alongside): " + "; ".join(
                f"Fig. {f.name[3]} `../figures/{f.name}`" for f in figs) + ".\n\n"
            clean = clean.replace("## Figure legends\n\n", "## Figure legends\n\n" + listing, 1)
            (OUT / name).write_text(clean)
            counts["abstract"] = words(section(clean, "## Abstract", ("\n## ", "\n---")))
            body = "".join(section(clean, h, ("\n## Methods", "\n## Discussion", "\n## Results", "\n---\n\n## M"))
                           for h in ("## Introduction", "## Results", "## Discussion"))
            counts["main_text_intro_results_discussion"] = words(body)
            counts["methods"] = words(section(clean, "## Methods", ("\n## Data availability",)))
            counts["figure_legends"] = words(section(clean, "## Figure legends", ("\n## Tables",)))
        else:
            counts["supplementary"] = words(clean)
    (OUT / "word_counts.json").write_text(json.dumps(counts, indent=1) + "\n")
    for k, v in counts.items():
        print(f"{k:38s} {v:6d}")


if __name__ == "__main__":
    main()
