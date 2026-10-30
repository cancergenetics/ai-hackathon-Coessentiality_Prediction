"""Static journey timeline (figures/journey.png) from data/journey.json, for the Markdown version of JOURNEY.md.
The HTML page draws the same data interactively (template.html). Usage: .venv/bin/python docs/explainer/make_journey_figure.py
"""
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
LANES = [("human", "Human input"), ("decision", "Decisions"), ("result", "Results"), ("failure", "Failures"),
         ("fix", "Fixes"), ("review", "Reviews")]
COL = {"human": "#1b9e77", "failure": "#d95f02"}
BLUE = "#2a78d6"


def ts(s):
    s = s if len(s) > 16 else s + ":00"
    return datetime.fromisoformat(s.replace("Z", "")).replace(tzinfo=timezone.utc)


def main():
    data = json.loads((HERE / "data/journey.json").read_text())
    ev = data["events"]
    fig, ax = plt.subplots(figsize=(12, 4.4), dpi=150)
    t0, t1 = ts("2026-10-28T14:00"), ts("2026-10-29T19:00")
    for i, p in enumerate(data["phases"]):
        a, b = ts(p["start"]), ts(p["end"])
        if i % 2 == 0:
            ax.axvspan(a, b, color="#eef2f8", zorder=0)
        ax.text(a + (b - a) / 2, len(LANES) - 0.45 + (0.32 if i % 2 else 0), p["name"], ha="center", va="bottom",
                fontsize=7.5, color="#5a5f69")
    for li, (k, label) in enumerate(LANES):
        y = len(LANES) - 1 - li
        ax.axhline(y, color="#d9dce2", lw=0.8, zorder=1)
        for e in [e for e in ev if e["kind"] == k]:
            x, tp, c = ts(e["datetime"]), e.get("tp"), COL.get(k, BLUE)
            if tp:
                ax.scatter(x, y, s=150, color=c, zorder=3, edgecolor="white", linewidth=0.8)
                ax.text(x, y, str(tp), ha="center", va="center", fontsize=7, color="white", fontweight="bold", zorder=4)
            else:
                ax.scatter(x, y, s=40, facecolor="white", edgecolor=c, linewidth=1.6, zorder=3)
    ax.set_yticks(range(len(LANES)))
    ax.set_yticklabels([l for _, l in LANES][::-1], fontsize=9)
    ax.set_ylim(-0.6, len(LANES) + 0.25)
    ax.set_xlim(t0, t1)
    ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, 3), tz=timezone.utc))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d Oct\n%H:%M", tz=timezone.utc))
    ax.tick_params(axis="x", labelsize=8)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    fig.suptitle("The project over 28–29 October 2026 (UTC). Numbered dots: the eight turning points.", fontsize=10,
                 x=0.01, ha="left", color="#14161a")
    fig.tight_layout()
    out = HERE / "figures/journey.png"
    fig.savefig(out, facecolor="white")
    print("wrote", out)


if __name__ == "__main__":
    main()
