"""Pre-registration gate and look ledger (docs/TEST_PROTOCOL.md §8).

Freeze (made by the PI, not by an agent): the protocol header reads "FROZEN" and the protocol contains no "draft"
string; ``scripts/22_freeze_protocol.py`` writes docs/TEST_PROTOCOL.FROZEN listing the SHA-256 of the protocol,
fitted_rules.json, every harness file, the venv package list and the two test files (which must equal DATA.md); the
PI commits it. ``require_frozen()`` re-verifies every hash, that the marker and all listed files are committed and
unmodified (src/, scripts/, tests/, docs/TEST_PROTOCOL*, fitted rules), and that the protocol is not a draft.

Look ledger: results/test_eval/looks.jsonl is git-tracked. Two steps: ``--intent`` writes an 'intent' record (the PI
commits it), ``--run`` requires that committed intent and refuses a second look found on disk or anywhere in git
history unless --relook-reason is given, in which case the run is automatically non-confirmatory.
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

from testeval import paths
from testeval.paths import PROJECT, PROTOCOL

MARKER = PROJECT / "docs" / "TEST_PROTOCOL.FROZEN"
CLEAN_PATHS = ("src", "scripts", "tests", "docs/TEST_PROTOCOL.md", "docs/TEST_PROTOCOL.FROZEN",
               "results/test_eval/fitted_rules.json")
DATA_MD = PROJECT / "DATA.md"


def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def _git(*a) -> tuple[int, str]:
    r = subprocess.run(["git", "-C", str(PROJECT), *a], capture_output=True, text=True)
    return r.returncode, r.stdout.strip()


def rel(p: Path) -> str:
    return str(Path(p).resolve().relative_to(PROJECT.resolve()))


def harness_files() -> list[Path]:
    """Every file whose content defines the look (sorted)."""
    fs = sorted((PROJECT / "src" / "testeval").glob("*.py"))
    fs += sorted(p for p in (PROJECT / "scripts").iterdir() if re.match(r"^(19|20|20b|21|22)_", p.name))
    fs += [PROJECT / "tests" / "test_testeval.py"]
    return fs


def pip_freeze_sha256() -> str:
    out = subprocess.run([sys.executable, "-m", "pip", "freeze", "--all"], capture_output=True, text=True).stdout
    lines = sorted(l for l in out.splitlines() if l.strip())
    return hashlib.sha256("\n".join(lines).encode()).hexdigest()


def data_md_test_hashes(text: str | None = None) -> dict:
    """SHA-256 of the two test files as recorded in DATA.md."""
    text = DATA_MD.read_text() if text is None else text
    out = {}
    for name in ("Test_positive_pairs.txt", "Test_negative_pairs.txt"):
        m = re.search(r"huMAP/" + re.escape(name) + r".*?\b([0-9a-f]{64})\b", text)
        if not m:
            raise RuntimeError(f"DATA.md has no SHA-256 for {name}")
        out[name] = m.group(1)
    return out


def check_protocol_text(text: str) -> None:
    """The frozen protocol must say FROZEN in its first heading block and contain no 'draft' anywhere."""
    if re.search(r"draft", text, re.IGNORECASE):
        raise SystemExit("docs/TEST_PROTOCOL.md still contains the word 'draft': not a frozen protocol")
    head = "\n".join(text.splitlines()[:6])
    if "FROZEN" not in head:
        raise SystemExit("docs/TEST_PROTOCOL.md header does not read FROZEN")


def build_marker(protocol_text: str | None = None) -> dict:
    """Content of docs/TEST_PROTOCOL.FROZEN (computed by scripts/22_freeze_protocol.py, run by the PI)."""
    text = PROTOCOL.read_text() if protocol_text is None else protocol_text
    check_protocol_text(text)
    pos, neg = paths.test_files()
    want = data_md_test_hashes()
    got = {pos.name: sha256_file(pos), neg.name: sha256_file(neg)}
    if got != want:
        raise SystemExit(f"test files differ from DATA.md: {got} != {want}")
    files = {rel(PROTOCOL): sha256_file(PROTOCOL), rel(paths.FITTED_RULES): sha256_file(paths.FITTED_RULES)}
    files.update({rel(p): sha256_file(p) for p in harness_files()})
    return dict(kind="TEST_PROTOCOL.FROZEN", version=1, files=files, test_files=got, pip_freeze_sha256=pip_freeze_sha256(),
                part1_locked_commit=paths.LOCKED_COMMIT)


def verify_marker(marker: dict, current: dict) -> list[str]:
    """Differences between a marker and the recomputed state (empty list = OK)."""
    errs = []
    for k in ("files", "test_files"):
        a, b = marker.get(k, {}), current.get(k, {})
        for f in sorted(set(a) | set(b)):
            if a.get(f) != b.get(f):
                errs.append(f"{k}:{f}: marker {a.get(f)} != current {b.get(f)}")
    for k in ("pip_freeze_sha256", "part1_locked_commit"):
        if marker.get(k) != current.get(k):
            errs.append(f"{k}: marker {marker.get(k)} != current {current.get(k)}")
    return errs


def git_state(paths_=CLEAN_PATHS) -> dict:
    _, head = _git("rev-parse", "HEAD")
    _, dirty = _git("status", "--porcelain", "--untracked-files=all", "--", *paths_)
    rc, _ = _git("ls-files", "--error-unmatch", rel(MARKER)) if MARKER.exists() else (1, "")
    _, marker_commit = _git("log", "-n", "1", "--format=%H", "--", rel(MARKER))
    _, proto_commit = _git("log", "-n", "1", "--format=%H", "--", rel(PROTOCOL))
    return dict(head=head, dirty=[l for l in dirty.splitlines() if l.strip()], marker_tracked=rc == 0,
                marker_commit=marker_commit or None, protocol_last_commit=proto_commit or None,
                protocol_sha256=sha256_file(PROTOCOL) if PROTOCOL.exists() else None,
                marker_sha256=sha256_file(MARKER) if MARKER.exists() else None)


def require_frozen() -> dict:
    """Raise SystemExit unless the protocol is frozen (marker committed, every hash equal, clean tree)."""
    if not MARKER.exists():
        raise SystemExit("docs/TEST_PROTOCOL.FROZEN missing: the PI freezes the protocol (scripts/22_freeze_protocol.py "
                         "+ commit) before any test scoring or look")
    st = git_state()
    if not st["marker_tracked"] or not st["marker_commit"]:
        raise SystemExit("docs/TEST_PROTOCOL.FROZEN is not committed")
    if st["dirty"]:
        raise SystemExit(f"uncommitted changes in frozen paths: {st['dirty']}")
    marker = json.loads(MARKER.read_text())
    errs = verify_marker(marker, build_marker())
    if errs:
        raise SystemExit("frozen state differs from docs/TEST_PROTOCOL.FROZEN:\n  " + "\n  ".join(errs))
    st["marker"] = marker
    return st


# ------------------------------------------------------------------------------------------------- ledger
def append_jsonl(path: Path, rec: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rec = dict(rec, time=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"))
    with open(path, "a") as fh:
        fh.write(json.dumps(rec, default=str) + "\n")


def read_jsonl(path: Path) -> list:
    if not Path(path).exists():
        return []
    return [json.loads(l) for l in Path(path).read_text().splitlines() if l.strip()]


def history_commits(path: Path) -> list[str]:
    """Commits (any ref) that touched ``path``."""
    _, out = _git("log", "--all", "--format=%H", "--", rel(path))
    return [l for l in out.splitlines() if l.strip()]


def prior_look_evidence(looks: Path, results: Path) -> dict:
    """Everything that shows a look (started / done / results) already happened, on disk or in git history."""
    recs = read_jsonl(looks)
    ran = [r for r in recs if r.get("status") in ("started", "done")]
    hist_results = history_commits(results)
    hist_ran = []
    for c in history_commits(looks):
        rc, txt = _git("show", f"{c}:{rel(looks)}")
        if rc == 0 and any(json.loads(l).get("status") in ("started", "done") for l in txt.splitlines() if l.strip()):
            hist_ran.append(c)
    return dict(ledger_run_records=len(ran), results_on_disk=Path(results).exists(),
                results_in_history=hist_results, ledger_runs_in_history=hist_ran,
                any=bool(ran or Path(results).exists() or hist_results or hist_ran))
