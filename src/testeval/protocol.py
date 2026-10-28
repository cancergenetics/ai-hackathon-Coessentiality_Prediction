"""Pre-registration gate: the test look may run only from a clean commit that contains docs/TEST_PROTOCOL.md."""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import subprocess
from pathlib import Path

from testeval.paths import PROJECT, PROTOCOL


def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def _git(*a) -> tuple[int, str]:
    r = subprocess.run(["git", "-C", str(PROJECT), *a], capture_output=True, text=True)
    return r.returncode, r.stdout.strip()


def git_state(paths=("src", "scripts", "docs/TEST_PROTOCOL.md")) -> dict:
    _, head = _git("rev-parse", "HEAD")
    _, dirty = _git("status", "--porcelain", "--", *paths)
    rc, _ = _git("ls-files", "--error-unmatch", str(PROTOCOL.relative_to(PROJECT)))
    _, proto_commit = _git("log", "-n", "1", "--format=%H", "--", str(PROTOCOL.relative_to(PROJECT)))
    return dict(head=head, dirty=[l for l in dirty.splitlines() if l.strip()], protocol_tracked=rc == 0,
                protocol_last_commit=proto_commit or None,
                protocol_sha256=sha256_file(PROTOCOL) if PROTOCOL.exists() else None)


def require_frozen() -> dict:
    """Raise unless TEST_PROTOCOL.md is committed and src/, scripts/ and the protocol are unmodified."""
    st = git_state()
    if not st["protocol_tracked"] or not st["protocol_last_commit"]:
        raise SystemExit("docs/TEST_PROTOCOL.md is not committed: freeze the protocol (git commit) before the look")
    if st["dirty"]:
        raise SystemExit(f"uncommitted changes in src/, scripts/ or the protocol: {st['dirty']}")
    return st


def append_jsonl(path: Path, rec: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rec = dict(rec, time=dt.datetime.now().isoformat(timespec="seconds"))
    with open(path, "a") as fh:
        fh.write(json.dumps(rec, default=str) + "\n")


def read_jsonl(path: Path) -> list:
    if not Path(path).exists():
        return []
    return [json.loads(l) for l in Path(path).read_text().splitlines() if l.strip()]
