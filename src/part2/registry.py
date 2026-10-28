"""Machine-readable analysis registry (FZ-R; DESIGN v3.1 §0 F4/F5): docs/registry/<id>.yaml.

The pinned environment has no YAML library, so this module parses a strict YAML subset (documented in
docs/registry/README.md): block mappings ``key: value``, block sequences ``- item`` (scalars or mappings),
flow sequences ``[a, b]``, block scalars ``|`` / ``>``, single/double-quoted strings, ints, floats,
true/false/null, ``#`` comments. Anything else (anchors, tags, multi-document streams, flow mappings) is rejected.

An entry is usable by run_association / read_lockbox only if it validates AND its file is committed in git with no
uncommitted modification (the committed blob SHA is recorded in every trial / ledger line).
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

DIRECTIONS = ("H3a", "H3b", "H1", "H2", "controls", "benchmark", "bridge")
CONFIRMATORY = ("controls", "H3a", "H3b", "H1", "H2")       # lockbox order: these before benchmark / bridge
SIGNS = ("+", "-", "two-sided", "per-context")
SCORES = ("S_C_mean_winsor3", "eigengene", "within_line_rank", "expression_masked", "pan_essential_referenced",
          "gene", "submodule_vs_sibling", "member_resolved")
STATUSES = ("open", "closed")
REQUIRED = {
    "id": str, "direction": str, "hypothesis": str, "sign": str, "context_family": dict, "complex_set_hash": str,
    "score": str, "covariates": list, "fdr_family": dict, "success_criterion": str, "stopping_rule": str,
    "frozen_inputs": list, "status": str,
}


class RegistryError(Exception):
    pass


# ------------------------------------------------------------------------------------------------ YAML subset
_NUM = re.compile(r"^[-+]?(\d+\.?\d*([eE][-+]?\d+)?|\.\d+([eE][-+]?\d+)?)$")


def _strip_comment(line: str) -> str:
    out, q = [], None
    for i, ch in enumerate(line):
        if q:
            if ch == q:
                q = None
        elif ch in "'\"":
            q = ch
        elif ch == "#" and (i == 0 or line[i - 1] in " \t"):
            break
        out.append(ch)
    return "".join(out).rstrip()


def _scalar(s: str):
    s = s.strip()
    if s == "" or s in ("null", "~", "Null", "NULL"):
        return None
    if s[0] in "'\"":
        if len(s) < 2 or s[-1] != s[0]:
            raise RegistryError(f"unterminated string: {s}")
        body = s[1:-1]
        return body.replace("''", "'") if s[0] == "'" else bytes(body, "utf-8").decode("unicode_escape")
    if s in ("true", "True", "TRUE"):
        return True
    if s in ("false", "False", "FALSE"):
        return False
    if s[0] == "[":
        if s[-1] != "]":
            raise RegistryError(f"unterminated flow sequence: {s}")
        inner = s[1:-1].strip()
        return [] if not inner else [_scalar(x) for x in _split_flow(inner)]
    if s[0] in "{&*!|>%@`":
        raise RegistryError(f"unsupported YAML construct: {s}")
    if _NUM.match(s):
        return int(s) if re.match(r"^[-+]?\d+$", s) else float(s)
    return s


def _split_flow(s: str):
    out, cur, q = [], [], None
    for ch in s:
        if q:
            cur.append(ch)
            if ch == q:
                q = None
        elif ch in "'\"":
            q = ch
            cur.append(ch)
        elif ch == ",":
            out.append("".join(cur))
            cur = []
        elif ch in "[]{}":
            raise RegistryError("nested flow collections are not supported")
        else:
            cur.append(ch)
    out.append("".join(cur))
    return out


def _key_value(s: str):
    m = re.match(r"^([A-Za-z_][A-Za-z0-9_.\-]*|\"[^\"]+\"|'[^']+')\s*:(\s+(.*))?$", s)
    if not m:
        return None
    k = m.group(1).strip("\"'")
    return k, (m.group(3) or "").strip()


def load_yaml(text: str):
    raw = text.splitlines()
    if any(l.strip() == "---" for l in raw[1:]) or any("\t" in l[: len(l) - len(l.lstrip())] for l in raw):
        raise RegistryError("multi-document streams and tab indentation are not supported")
    lines = []                                               # (indent, content, raw line index)
    for i, l in enumerate(raw):
        if l.strip() == "---" and i == 0:
            continue
        c = _strip_comment(l)
        if c.strip():
            lines.append((len(c) - len(c.lstrip(" ")), c.strip(), i))
    pos = [0]

    def block_scalar(style: str, parent_indent: int, start_raw: int):
        # consume raw lines after start_raw that are blank or indented deeper than parent
        buf, ind = [], None
        j = start_raw + 1
        while j < len(raw):
            l = raw[j]
            if l.strip() == "":
                buf.append("")
                j += 1
                continue
            li = len(l) - len(l.lstrip(" "))
            if li <= parent_indent:
                break
            ind = li if ind is None else min(ind, li)
            buf.append(l)
            j += 1
        while pos[0] < len(lines) and lines[pos[0]][2] < j:
            pos[0] += 1
        body = [b[ind:] if b else "" for b in buf] if ind is not None else []
        while body and body[-1] == "":
            body.pop()
        if style.startswith("|"):
            return "\n".join(body) + "\n"
        paras, cur = [], []
        for b in body:
            if b == "":
                paras.append(" ".join(cur))
                cur = []
            else:
                cur.append(b.strip())
        paras.append(" ".join(cur))
        return "\n".join(paras).strip() + "\n"

    def value_after(v: str, indent: int, raw_i: int):
        if v in ("|", ">", "|-", ">-", "|+", ">+"):
            return block_scalar(v, indent, raw_i)
        if v == "":
            if pos[0] < len(lines) and (lines[pos[0]][0] > indent or
                                        (lines[pos[0]][0] == indent and lines[pos[0]][1].startswith("- "))):
                return parse_block(lines[pos[0]][0])
            return None
        return _scalar(v)

    def parse_block(indent: int):
        first = lines[pos[0]][1]
        if first.startswith("- ") or first == "-":
            out = []
            while pos[0] < len(lines) and lines[pos[0]][0] == indent and (lines[pos[0]][1].startswith("- ")
                                                                           or lines[pos[0]][1] == "-"):
                _, c, ri = lines[pos[0]]
                item = c[1:].strip()
                pos[0] += 1
                kv = _key_value(item) if item else None
                if kv is not None and not item.startswith(("'", '"')):
                    # mapping inside a list item: first key inline, further keys at indent + 2
                    d = {}
                    k, v = kv
                    d[k] = value_after(v, indent + 2, ri)
                    while pos[0] < len(lines) and lines[pos[0]][0] == indent + 2 and not lines[pos[0]][1].startswith("- "):
                        _, c2, ri2 = lines[pos[0]]
                        kv2 = _key_value(c2)
                        if kv2 is None:
                            raise RegistryError(f"line {ri2 + 1}: expected 'key: value'")
                        pos[0] += 1
                        if kv2[0] in d:
                            raise RegistryError(f"line {ri2 + 1}: duplicate key {kv2[0]}")
                        d[kv2[0]] = value_after(kv2[1], indent + 2, ri2)
                    out.append(d)
                elif item == "":
                    out.append(parse_block(lines[pos[0]][0]) if pos[0] < len(lines) and lines[pos[0]][0] > indent
                               else None)
                else:
                    out.append(_scalar(item))
            return out
        d = {}
        while pos[0] < len(lines) and lines[pos[0]][0] == indent:
            _, c, ri = lines[pos[0]]
            if c.startswith("- "):
                break
            kv = _key_value(c)
            if kv is None:
                raise RegistryError(f"line {ri + 1}: expected 'key: value', got {c!r}")
            pos[0] += 1
            if kv[0] in d:
                raise RegistryError(f"line {ri + 1}: duplicate key {kv[0]}")
            d[kv[0]] = value_after(kv[1], indent, ri)
        return d

    if not lines:
        return None
    out = parse_block(lines[0][0])
    if pos[0] != len(lines):
        raise RegistryError(f"line {lines[pos[0]][2] + 1}: unexpected indentation")
    return out


# ------------------------------------------------------------------------------------------------ validation
def validate(entry: dict, stem: str | None = None) -> list[str]:
    errs = []
    if not isinstance(entry, dict):
        return ["entry is not a mapping"]
    for k, t in REQUIRED.items():
        if k not in entry or entry[k] is None:
            errs.append(f"missing field: {k}")
        elif not isinstance(entry[k], t):
            errs.append(f"{k}: expected {t.__name__}")
    if errs:
        return errs
    if stem is not None and entry["id"] != stem:
        errs.append(f"id {entry['id']!r} != file name {stem!r}")
    if not re.match(r"^[A-Za-z0-9_.\-]+$", entry["id"]):
        errs.append("id: only [A-Za-z0-9_.-]")
    if entry["direction"] not in DIRECTIONS:
        errs.append(f"direction must be one of {DIRECTIONS}")
    if entry["sign"] not in SIGNS:
        errs.append(f"sign must be one of {SIGNS}")
    if entry["score"] not in SCORES:
        errs.append(f"score must be one of {SCORES}")
    if entry["status"] not in STATUSES:
        errs.append(f"status must be one of {STATUSES}")
    if not re.match(r"^[0-9a-f]{16,64}$", str(entry["complex_set_hash"])):
        errs.append("complex_set_hash: 16-64 lowercase hex characters (part2.registry.set_hash of the tested units)")
    if not all(isinstance(c, str) and c for c in entry["covariates"]):
        errs.append("covariates: list of non-empty strings")
    if not all(isinstance(c, str) and c for c in entry["frozen_inputs"]):
        errs.append("frozen_inputs: list of manifest paths")
    cf = entry["context_family"]
    if not isinstance(cf.get("name"), str):
        errs.append("context_family.name missing")
    ctxs = cf.get("contexts")
    if not isinstance(ctxs, list) or not ctxs:
        errs.append("context_family.contexts: non-empty list")
    else:
        ids = []
        for i, c in enumerate(ctxs):
            if not isinstance(c, dict) or not isinstance(c.get("id"), str) or not isinstance(c.get("definition"), str):
                errs.append(f"context_family.contexts[{i}]: needs id and definition")
                continue
            ids.append(c["id"])
            if entry["sign"] == "per-context" and c.get("sign") not in ("+", "-", "two-sided"):
                errs.append(f"context {c['id']}: sign required when sign is per-context")
        if len(set(ids)) != len(ids):
            errs.append("context ids not unique")
    ff = entry["fdr_family"]
    if ff.get("method") not in ("benjamini_bogomolov", "bh", "none"):
        errs.append("fdr_family.method: benjamini_bogomolov | bh | none")
    if not isinstance(ff.get("q"), (int, float)) or not 0 < ff["q"] < 1:
        errs.append("fdr_family.q in (0, 1)")
    if ff.get("method") == "benjamini_bogomolov" and ff.get("families") not in ("contexts", "complexes"):
        errs.append("fdr_family.families: contexts | complexes")
    return errs


def set_hash(units) -> str:
    """Hash of a set of tested units (complex IDs / genes): SHA-256 (first 16 hex) of sorted, newline-joined strings."""
    import hashlib
    u = sorted(set(map(str, units)))
    return hashlib.sha256("\n".join(u).encode()).hexdigest()[:16]


def contexts_of(entry: dict) -> dict:
    """context id -> {'definition', 'sign'} (sign resolved from the entry when not per-context)."""
    return {c["id"]: dict(definition=c["definition"], sign=c.get("sign", entry["sign"]))
            for c in entry["context_family"]["contexts"]}


# ------------------------------------------------------------------------------------------------ git
def git(root, *args, check=True) -> subprocess.CompletedProcess:
    return subprocess.run(["git", "-C", str(root), *args], capture_output=True, text=True, check=check)


def committed_blob(root, path) -> str:
    """Blob SHA of `path` at HEAD; raises RegistryError if untracked or modified relative to HEAD."""
    root = Path(root)
    rel = str(Path(path).resolve().relative_to(root.resolve()))
    if git(root, "ls-files", "--error-unmatch", "--", rel, check=False).returncode != 0:
        raise RegistryError(f"{rel} is not tracked by git (commit it first)")
    if git(root, "status", "--porcelain", "--", rel).stdout.strip():
        raise RegistryError(f"{rel} has uncommitted changes")
    r = git(root, "rev-parse", f"HEAD:{rel}", check=False)
    if r.returncode != 0:
        raise RegistryError(f"{rel} is not in HEAD")
    return r.stdout.strip()


def load_entry(registry_dir, analysis_id: str, root=None, require_committed: bool = True) -> tuple[dict, str | None]:
    if not re.match(r"^[A-Za-z0-9_.\-]+$", str(analysis_id)):
        raise RegistryError(f"bad analysis id {analysis_id!r}")
    p = Path(registry_dir) / f"{analysis_id}.yaml"
    if not p.exists():
        raise RegistryError(f"analysis {analysis_id!r} is not registered ({p} missing)")
    entry = load_yaml(p.read_text())
    errs = validate(entry, analysis_id)
    if errs:
        raise RegistryError(f"{p.name}: " + "; ".join(errs))
    blob = committed_blob(root, p) if require_committed else None
    return entry, blob


def load_all(registry_dir) -> dict:
    out = {}
    for p in sorted(Path(registry_dir).glob("*.yaml")):
        e = load_yaml(p.read_text())
        errs = validate(e, p.stem)
        if errs:
            raise RegistryError(f"{p.name}: " + "; ".join(errs))
        out[p.stem] = e
    return out
