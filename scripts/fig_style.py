"""Shared style and loaders for the manuscript figures and tables (scripts/fig_*.py, scripts/tab_manuscript.py).

Every number shown in a figure or table is read from a committed result file; nothing is hardcoded.
The palette is Okabe-Ito (colourblind-safe). Output is deterministic: fixed seed, fixed hash salt,
no creation date in the PDF metadata.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "manuscript" / "figures"
TAB_DIR = ROOT / "manuscript" / "tables"
SEED = 0

# Committed inputs (paths relative to the repository root)
TEST_RESULTS = "results/test_eval/test_results.json"
SCORE_RECORDS = "results/test_eval/score_records/test"
CONSTR = "results/part2/construction"
BRIDGE = f"{CONSTR}/bridge/bridge.json"
TIERS = f"{CONSTR}/tiers/summary.json"
SUMMARY = f"{CONSTR}/summary/summary.json"
SELECTION = f"{CONSTR}/selection/selection.json"
VARIAB = f"{CONSTR}/variability/variability.json"
JAC_ZAVG = f"{CONSTR}/variability/zavg_BEF_per_complex_jaccard.tsv"
JAC_GLS = f"{CONSTR}/variability/gls_per_complex_jaccard.tsv"
CONTROLS = "results/part2/controls/controls_summary.tsv"
KILL = "results/part2/controls/kill_rule.json"
CONTROL_UNITS = "results/part2/fzx/control_units.tsv"

# Okabe-Ito
OI = {
    "black": "#000000", "orange": "#E69F00", "skyblue": "#56B4E9", "green": "#009E73",
    "yellow": "#F0E442", "blue": "#0072B2", "vermillion": "#D55E00", "purple": "#CC79A7", "grey": "#999999",
}
METHODS = ["zavg_BEF", "B", "gls_signed", "spearman_signed", "pearson_signed", "spec_k20_tinf_std"]
BASELINES = ["gls_signed", "spearman_signed", "pearson_signed"]
LABEL = {
    "zavg_BEF": "zavg_BEF (locked)", "B": "B (secondary candidate)", "gls_signed": "GLS",
    "spearman_signed": "Spearman", "pearson_signed": "Pearson",
    "spec_k20_tinf_std": "PC-removed Pearson (post hoc ref.)",
}
SHORT = {"zavg_BEF": "zavg_BEF", "B": "B", "gls_signed": "GLS", "spearman_signed": "Spearman",
         "pearson_signed": "Pearson", "spec_k20_tinf_std": "PC-rm Pearson"}
COLOR = {
    "zavg_BEF": OI["blue"], "B": OI["skyblue"], "gls_signed": OI["vermillion"],
    "spearman_signed": OI["green"], "pearson_signed": OI["orange"], "spec_k20_tinf_std": OI["purple"],
    "zavg": OI["blue"], "gls": OI["vermillion"], "pearson": OI["orange"],
}
MARKER = {"zavg_BEF": "o", "B": "s", "gls_signed": "D", "spearman_signed": "^", "pearson_signed": "v",
          "spec_k20_tinf_std": "P"}
TIER_COLOR = {"T1": OI["blue"], "T2": OI["skyblue"], "T3-neg": OI["grey"], "T3-nc": "#DDDDDD"}


def setup() -> None:
    np.random.seed(SEED)
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
        "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6, "axes.linewidth": 0.6,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "axes.spines.top": False,
        "axes.spines.right": False, "svg.hashsalt": "manuscript", "pdf.fonttype": 42,
        "figure.dpi": 150, "savefig.dpi": 300,
    })


def load_json(rel: str) -> dict:
    return json.loads((ROOT / rel).read_text())


def panel_label(ax, s: str) -> None:
    ax.text(-0.14, 1.06, s, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="left")


def save(fig, stem: str) -> list[Path]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    out = []
    for ext, kw in (("pdf", {"metadata": {"CreationDate": None, "ModDate": None, "Producer": None}}),
                    ("png", {"metadata": {"Software": None}})):
        p = FIG_DIR / f"{stem}.{ext}"
        fig.savefig(p, bbox_inches="tight", **kw)
        out.append(p)
    plt.close(fig)
    return out


def fmt_ci(point: float, lo: float, hi: float, nd: int = 3, signed: bool = False) -> str:
    f = f"{{:+.{nd}f}}" if signed else f"{{:.{nd}f}}"
    return f"{f.format(point)} [{f.format(lo)}, {f.format(hi)}]"
