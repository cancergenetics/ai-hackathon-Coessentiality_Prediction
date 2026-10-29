#!/usr/bin/env python3
"""Check every number in a manuscript Markdown file against the source named in its HTML comment.

Convention (manuscript/main.md): a sentence or short run of sentences is followed by an HTML comment naming the
committed file(s) its numbers come from, optionally with JSON key paths, TSV rows/columns or a quoted section, e.g.

    ... AP 0.333 (0.256-0.410) ... <!-- TR: results.itt.bootstrap.cluster.stat.zavg_BEF|ap -->
    ... 601 of 1,145 ... <!-- CMP: H3a_lineage_atlas overall -->
    ... L -0.202 ... <!-- results/part2/robustness/strongest_findings.tsv rows Z0035: beta_L, ci95_L_lo -->

Aliases are declared in the first HTML comment of the file as lines "ALIAS = path".

For each comment the script takes the "claim text" (the text of the same paragraph since the previous comment),
extracts every number, collects every numeric value reachable from the cited sources, and calls a number
MATCHED if some source value equals it at the displayed precision (absolute value; x100 also tried for
percentages; mantissa compared for "a x 10^-b"). The check is deliberately lenient about *which* key supplies a
value, so a MATCHED number is supported by the cited source, not proven to come from the intended key.

A comment may add arithmetic for derived numbers: "calc: 3642/2862 - 1" adds that value to the candidates.

Statuses: matched (with the support kind: key, row, calc, file or text); UNMATCHED (cited source has no value that rounds to the number); NO-SOURCE (number in a
paragraph segment with no following comment); UNRESOLVED (the comment names no readable numeric source, e.g. only
"git log"). Numbers that are identifiers (gene symbols, unit IDs, figure/section/table numbers, commits, dates,
times, versions, citation years) are removed before extraction.

Usage:  python scripts/check_manuscript_numbers.py [manuscript/main.md] [--out results/verify_manuscript/number_check.tsv]
Reads committed files only (no dataset, no lockbox read). Exit status 1 if any number is UNMATCHED.
"""
from __future__ import annotations

import argparse
import csv
import fnmatch
import glob
import gzip
import io
import json
import math
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SUP = {"⁻": "-", "⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4", "⁵": "5", "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9"}
EXTS = (".json", ".jsonl", ".tsv", ".tsv.gz", ".md", ".yaml", ".yml", ".py", ".txt", ".csv")
SEARCH_DIRS = ("results", "manuscript", "docs", "scripts", "src")

# ----------------------------------------------------------------------------------------------- number parsing
NUM_RE = re.compile(
    r"(?<![\w.\-−/])(?P<sign>[+\-−±]?)(?P<int>\d{1,3}(?:,\d{3})+|\d+)(?:\.(?P<dec>\d+))?"
    r"(?:\s*[×x]\s*10(?P<exp>[⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+))?(?P<pct>\s*%)?(?![\w])")
POW_RE = re.compile(r"(?<![\d×x]\s)(?<![\d×])\b10(?P<exp>[⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+)")

STRIP_PATTERNS = [
    r"\bParts?\s+\d+(?:\s+and\s+\d+)?\b",
    r"`[^`]*`",                                                     # code spans
    r"\[[^\]]*\d{4}[^\]]*\]",                                       # citations [Name 2019; ...]
    r"\[M\d+\]",
    r"\(?Figs?\.\s*\d+[a-z]?(?:\s*[–,\-]\s*\d*[a-z])*\)?",          # figure refs
    r"Supplementary (?:Notes?|Tables?|Information)\s*S?\d+(?:\s*[–\-]\s*S?\d+)?",
    r"\bTables?\s+S?\d+[a-z]?(?:\s*[–\-]\s*S?\d+)?",
    r"§+\s*\d+(?:\.\d+)?(?:\s*[–\-]\s*\d+(?:\.\d+)?)?",
    r"\b(?:CTRL_)?PC\d+[ab]?(?:/PC\d+)?\b",
    r"\bZ\d{4}(?:\.\d+)?\b",
    r"\b[0-9a-f]{7,40}\b(?<=[a-f]\w{6})|\b(?=\w*[a-f])(?=\w*\d)[0-9a-f]{7,40}\b",   # commit hashes
    r"\b\d{1,2}:\d{2}(?::\d{2})?\b",                                # times
    r"\b\d{1,2}\s+(?:Oct|October|Nov|November)\b(?:\s+20\d\d)?",  # dates
    r"\b20\d\d-\d\d-\d\d\b",
    r"\b(?:19|20)\d\d\b",                                           # years
    r"\b\d{2}Q\d\b",                                                # DepMap release
    r"\b[A-Za-z]*[A-Z][A-Za-z]*\d+[A-Za-z0-9]*(?:[-:/][A-Za-z0-9]+)*\b",  # gene symbols, H3a, T1, C3, A1, chrY...
    r"\b\d+[pq]\d+(?:\.\d+)?\b",                                    # cytobands 9p21
    r"\b(?:hu\.MAP|BioPlex|Ensembl|Ensembl release|release|v|DEMETER|log)\s*\d+(?:\.\d+)?\b",
    r"\b(?:95|90)\s*%\s*(?=(?:CI|cluster|bootstrap|percentile|confidence|HC3|PatientID|protein|interval))",
    r"\b(?:95|97\.5|2\.5)(?:th)\b",
    r"\b\d+(?:st|nd|rd|th)\b",
    r"\bseeds?\s+\d+(?:\s*[–\-]\s*\d+)?",
    r"\b(?:one|two)-sided\b",
    r"\blevel[- ]\d\b", r"\bdegree-\d\b", r"\bR\d\b",
    r"\b[A-Z]\d?[–\-][A-Z]\d+\b",                                   # A3–A5 style ranges
]
STRIP_RE = re.compile("|".join(f"(?:{p})" for p in STRIP_PATTERNS))


def sup_to_int(s: str) -> int:
    return int("".join(SUP[c] for c in s))


def extract_numbers(text: str):
    """Return [(raw, value, decimals, is_pct, exp)] for numbers in text after removing identifiers."""
    clean = STRIP_RE.sub(lambda m: " " * len(m.group(0)), text)
    out = []
    taken = []
    for m in NUM_RE.finditer(clean):
        raw = m.group(0).strip()
        ip = m.group("int").replace(",", "")
        dec = m.group("dec") or ""
        val = float(ip + ("." + dec if dec else ""))
        exp = sup_to_int(m.group("exp")) if m.group("exp") else None
        if exp is None and ip == "10" and not dec:
            continue  # handled by POW_RE if it is a power
        out.append((raw, val, len(dec), bool(m.group("pct")), exp))
        taken.append(m.span())
    for m in POW_RE.finditer(clean):
        if any(a <= m.start() < b for a, b in taken):
            continue
        out.append((m.group(0), 1.0, 0, False, sup_to_int(m.group("exp"))))
    return out


# ----------------------------------------------------------------------------------------------- source values
def flatten(obj, path=()):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from flatten(v, path + (str(k),))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from flatten(v, path + (str(i),))
    else:
        yield path, obj


def num(v):
    if isinstance(v, bool) or v is None:
        return None
    if isinstance(v, (int, float)):
        return float(v) if math.isfinite(float(v)) else None
    if isinstance(v, str):
        s = v.strip().replace("−", "-").replace(",", "")
        try:
            f = float(s)
            return f if math.isfinite(f) else None
        except ValueError:
            return None
    return None


def cell_nums(c):
    """Numbers in a TSV cell: a plain number, or a bracketed list such as "[-0.59, 0.005]"."""
    v = num(c)
    if v is not None:
        return [v]
    if isinstance(c, str) and c.strip().startswith("["):
        return [float(x) for x in re.findall(r"-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?", c)]
    return []


def split_top(s: str, sep: str):
    """Split s on sep outside quotes, braces, brackets and parentheses."""
    out, depth, q, cur = [], 0, False, ""
    for c in s:
        if c == '"':
            q = not q
        elif not q and c in "{[(":
            depth += 1
        elif not q and c in "}])":
            depth = max(0, depth - 1)
        if c == sep and depth == 0 and not q:
            out.append(cur)
            cur = ""
        else:
            cur += c
    out.append(cur)
    return [x.strip() for x in out if x.strip()]


def expand_braces(seg: str):
    m = re.search(r"\{([^{}]*)\}", seg)
    if not m:
        return [seg]
    inner = m.group(1)
    alts = ["*"] if ".." in inner else inner.split(",")
    res = []
    for a in alts:
        res += expand_braces(seg[:m.start()] + a.strip() + seg[m.end():])
    return res


class Source:
    def __init__(self, path):
        self.path = path
        self.kind = None
        self.data = None
        self.header = None
        self.rows = None
        self.text = None
        p = path.lower()
        try:
            if p.endswith(".json"):
                self.kind, self.data = "json", json.load(open(path))
            elif p.endswith(".jsonl"):
                self.kind = "json"
                self.data = [json.loads(line) for line in open(path) if line.strip()]
            elif p.endswith(".tsv") or p.endswith(".tsv.gz") or p.endswith(".csv"):
                fh = gzip.open(path, "rt") if p.endswith(".gz") else open(path, newline="")
                rd = csv.reader(fh, delimiter="," if p.endswith(".csv") else "\t")
                rows = list(rd)
                self.kind, self.header, self.rows = "tsv", rows[0], rows[1:]
            else:
                self.kind, self.text = "text", open(path, errors="replace").read()
        except Exception as e:  # noqa: BLE001
            self.kind, self.text = "error", str(e)

    def num_leaves(self):
        """(path, value) for every numeric JSON leaf; numbers inside string leaves (rules, notes) count too."""
        out = []
        for p, x in flatten(self.data):
            v = num(x)
            if v is not None:
                out.append((p, v))
            elif isinstance(x, str):
                out += [(p, n[1] * (10.0 ** n[4] if n[4] is not None else 1)) for n in extract_numbers_loose(x)]
        return out

    # all numeric values
    def all_values(self):
        if self.kind == "json":
            return [v for _, v in self.num_leaves()]
        if self.kind == "tsv":
            return [v for r in self.rows for x in r for v in cell_nums(x)]
        if self.kind == "text":
            return [n[1] * (10 ** n[4] if n[4] is not None else 1) for n in extract_numbers_loose(self.text)]
        return []

    def json_values(self, keys):
        vals, unresolved, prev_parent = [], [], ()
        leaves = self.num_leaves()
        for key in keys:
            key = re.sub(r"\([^)]*\)", "", key)          # drop "(0.2814)" annotations
            key = re.split(r"\s+=|\s+\(|\s+—|\s+\"", key)[0].strip()
            key = key.split(" ")[0] if " " in key and "\"" not in key else key
            if not key or key in ("idem",):
                continue
            segs = [s.strip('"') for s in split_top(key, ".")]
            segs = [re.sub(r"\[[^\]]*\]$", "", s) for s in segs]
            found = self._match(leaves, segs)
            if not found and prev_parent:
                found = self._match(leaves, list(prev_parent) + segs)
            if not found:
                last = segs[-1] if segs else ""
                pats = expand_braces(last)
                found = [v for p, v in leaves
                         if any(fnmatch.fnmatchcase(".".join(p), f"*{pt}*") for pt in pats)]
            if found:
                vals += found
                if len(segs) > 1:
                    prev_parent = tuple(segs[:-1])
            else:
                unresolved.append(key)
        return vals, unresolved

    @staticmethod
    def _match(leaves, segs):
        pats = [expand_braces(s) for s in segs]
        out = []
        for p, v in leaves:
            if len(p) < len(pats):
                continue
            # allow list indices to be skipped: match pattern segments against path with list indices removed
            q = [s for s in p if not s.isdigit()] if not any(s.isdigit() for s in segs) else list(p)
            if len(q) < len(pats):
                continue
            if all(any(fnmatch.fnmatchcase(q[i], a) for a in pats[i]) for i in range(len(pats))):
                out.append(v)
        return out

    def tsv_values(self, selector: str):
        if not self.rows:
            return []
        sel, cols = selector, []
        m = re.match(r"(.*):\s+(.+)$", selector)   # "rows X Y: col1, col2" (colon followed by a space)
        if m and any(c.strip() in self.header for c in split_top(m.group(2), ",")):
            sel, cols = m.group(1), [c.strip() for c in split_top(m.group(2), ",")]
        tokens = [t.strip(" ,;()\"'") for t in re.split(r"[\s,]+", sel)]
        tokens = [t for t in tokens if t and t.lower() not in ("row", "rows", "and", "with", "column", "columns")]
        n = len(self.rows)
        specific = []
        for t in tokens:
            hits = sum(1 for r in self.rows if any(c == t or (len(t) >= 5 and t in c) for c in r))
            if 0 < hits <= max(3, 0.2 * n):
                specific.append(t)
        rows = self.rows
        if specific:
            rows = [r for r in self.rows if any(c == t or (len(t) >= 5 and t in c) for t in specific for c in r)]
        idx = [self.header.index(c) for c in cols if c in self.header] or range(len(self.header))
        return [v for r in rows for i in idx if i < len(r) for v in cell_nums(r[i])]


def extract_numbers_loose(text):
    """All numbers in a text source (Markdown, YAML, ...), including '7/14', '1e-4' and '10⁻⁴' forms."""
    out = []
    for m in re.finditer(r"(?<![\d.])[-−]?(\d{1,3}(?:,\d{3})+|\d+)(\.\d+)?(?:[eE]([+-]?\d+))?"
                         r"(?:\s*[×x]\s*10([⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+))?", text):
        v = float(m.group(1).replace(",", "") + (m.group(2) or ""))
        if m.group(3):
            e = int(m.group(3))
            if abs(e) > 300:   # e.g. hex or ID strings that look like exponents
                continue
            v *= 10.0 ** e
        exp = sup_to_int(m.group(4)) if m.group(4) else None
        out.append((m.group(0), v, len(m.group(2) or "") - 1 if m.group(2) else 0, False, exp))
    for m in re.finditer(r"\b10([⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+)", text):
        out.append((m.group(0), 1.0, 0, False, sup_to_int(m.group(1))))
    return out


# ----------------------------------------------------------------------------------------------- comments
def parse_aliases(md: str):
    first = re.search(r"<!--(.*?)-->", md, re.S)
    al = {}
    if first:
        for m in re.finditer(r"^\s*([A-Z][A-Z0-9]*)\s*=\s*(\S+)", first.group(1), re.M):
            al[m.group(1)] = m.group(2)
        # older prose form: '"TR" below abbreviates results/...'
        for m in re.finditer(r'"([A-Z][A-Z0-9]*)"(?:\s+below)?\s+abbreviates\s+(\S+?)[;,)]?(?=\s|$)', first.group(1)):
            al.setdefault(m.group(1), m.group(2))
        m = re.search(r'"RVD"\s+its output folder\s+(\S+?)\.?(?=\s|$)', first.group(1))
        if m:
            al.setdefault("RVD", m.group(1).rstrip("."))
    return al


class Resolver:
    def __init__(self, aliases):
        self.aliases = aliases
        self.cache = {}
        self.last_dir = None

    def expand(self, tok):
        tok = tok.strip().rstrip(":,")
        for _ in range(5):  # aliases may refer to other aliases (CMP = LB/...)
            head, _, rest = tok.partition("/")
            if head not in self.aliases:
                break
            tok = self.aliases[head].rstrip("/") + ("/" + rest if rest else "")
        return tok

    def is_file_token(self, tok):
        t = tok.rstrip(":,")
        head = t.split("/")[0]
        if head in self.aliases:
            return True
        return t.lower().endswith(EXTS) or ("/" in t and os.path.exists(os.path.join(ROOT, self.expand(t))))

    def files(self, tok):
        t = self.expand(tok)
        pats = expand_braces(t)
        out = []
        for p in pats:
            full = os.path.join(ROOT, p)
            g = sorted(glob.glob(full))
            if g:
                out += [x for x in g if os.path.isfile(x)]
                continue
            if "/" not in p:  # bare basename: last-used dir, then unique search
                cands = []
                if self.last_dir:
                    c = os.path.join(self.last_dir, p)
                    if os.path.isfile(c):
                        cands = [c]
                if not cands:
                    for d in SEARCH_DIRS:
                        cands += glob.glob(os.path.join(ROOT, d, "**", p), recursive=True)
                if len(cands) >= 1:
                    out.append(sorted(cands, key=len)[0])
        if out:
            self.last_dir = os.path.dirname(out[0])
        return out

    def load(self, path):
        if path not in self.cache:
            self.cache[path] = Source(path)
        return self.cache[path]

    def values(self, comment: str):
        """Return (values, sources, problems)."""
        vals, srcs, probs = [], [], []
        cur_files = []
        for spec in split_top(comment.replace("\n", " "), ";"):
            m = re.match(r"calc:\s*(.+)$", spec)
            if m:
                expr = m.group(1).replace("−", "-").replace("×", "*")
                for e in split_top(expr, ","):
                    if re.fullmatch(r"[\d\s.+\-*/()]+", e):
                        try:
                            vals.append((float(eval(e, {"__builtins__": {}})), "calc"))  # noqa: S307 (digits/operators only)
                        except Exception:  # noqa: BLE001
                            probs.append(f"bad calc {e}")
                continue
            parts = spec.split(None, 1)
            first = parts[0] if parts else ""
            rest = parts[1] if len(parts) > 1 else ""
            if first.endswith(":") and not self.is_file_token(first):
                pass
            if self.is_file_token(first):
                cur_files = self.files(first.rstrip(":"))
                if first.endswith(":"):
                    rest = ":" + rest
                if not cur_files:
                    probs.append(f"file not found: {first}")
                    continue
            elif not cur_files:
                continue  # e.g. "git log", "DESIGN reasoning" with no readable file
            else:
                rest = spec if spec.startswith(":") else (":" + spec if self.load(cur_files[0]).kind == "json"
                                                         else spec)
            for f in cur_files:
                s = self.load(f)
                srcs.append(os.path.relpath(f, ROOT))
                if s.kind == "json":
                    r = rest.strip()
                    if r.startswith(":"):
                        keys = split_top(r[1:], ",")
                        v, unres = s.json_values(keys)
                        vals += [(x, "key") for x in v]
                        if unres:
                            probs.append(f"{os.path.relpath(f, ROOT)}: unresolved keys {unres}")
                            vals += [(x, "file") for x in s.all_values()]  # lenient fallback, flagged
                    else:
                        vals += [(x, "file") for x in s.all_values()]
                elif s.kind == "tsv":
                    r = rest.strip().lstrip(":")
                    vals += [(x, "row") for x in s.tsv_values(r)] if r else [(x, "file") for x in s.all_values()]
                elif s.kind == "text":
                    vals += [(x, "text") for x in s.all_values()]
                else:
                    probs.append(f"unreadable {f}")
        return vals, srcs, probs


# ----------------------------------------------------------------------------------------------- matching
RANK = {"key": 0, "row": 0, "calc": 0, "file": 1, "text": 2}


def matches(n, vals):
    """Return the most specific support kind of a matching source value, or None."""
    raw, x, d, pct, exp = n
    tol = 0.5 * 10 ** (-d) + 1e-9
    best = None
    for v, kind in vals:
        a = abs(v)
        ok = False
        if exp is not None:
            ok = a != 0 and abs(a / 10 ** exp - x) <= tol * 1.0001
        elif abs(a - x) <= tol or ((pct or d == 0 or x > 1) and abs(a * 100 - x) <= tol):
            ok = True
        if ok and (best is None or RANK[kind] < RANK[best]):
            best = kind
            if RANK[kind] == 0:
                break
    return best


def segments(md: str):
    """Yield (line_no, claim_text, comment) for each non-header comment, plus trailing uncited text."""
    body_start = re.search(r"<!--.*?-->", md, re.S).end()
    pos = body_start
    for m in re.finditer(r"<!--(.*?)-->", md[body_start:], re.S):
        a, b = m.start() + body_start, m.end() + body_start
        text = md[pos:a]
        paras = re.split(r"\n\s*\n", text)
        # claim text is the last paragraph; if it is empty (comment on its own line), the one before
        claim = paras[-1] if paras[-1].strip() else (paras[-2] if len(paras) > 1 else "")
        uncited = paras[:-1] if paras[-1].strip() else paras[:-2]
        for u in uncited:
            yield md.count("\n", 0, pos) + 1, u, None
        yield md.count("\n", 0, a) + 1, claim, m.group(1)
        pos = b
    yield md.count("\n", 0, pos) + 1, md[pos:], None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("md", nargs="?", default="manuscript/main.md")
    ap.add_argument("--out", default="results/verify_manuscript/number_check.tsv")
    ap.add_argument("--stop-at", default="## Figure legends",
                    help="heading at which checking stops (legends, tables list and references are not checked)")
    args = ap.parse_args()
    md_full = open(os.path.join(ROOT, args.md)).read()
    stop = md_full.find(args.stop_at)
    md = md_full if stop < 0 else md_full[:stop]
    res = Resolver(parse_aliases(md_full))
    rows = []
    for line, claim, comment in segments(md):
        # drop headings and table separator noise
        claim_txt = "\n".join(l for l in claim.splitlines() if not l.lstrip().startswith("#"))
        nums = extract_numbers(claim_txt)
        if not nums:
            continue
        if comment is None:
            for n in nums:
                rows.append((line, n[0], "NO-SOURCE", "", "", ""))
            continue
        vals, srcs, probs = res.values(comment)
        for n in nums:
            kind = matches(n, vals) if vals else None
            st = "UNRESOLVED" if not vals else ("matched" if kind else "UNMATCHED")
            rows.append((line, n[0], st, kind or "", ",".join(sorted(set(srcs))), "; ".join(probs)))
    os.makedirs(os.path.dirname(os.path.join(ROOT, args.out)), exist_ok=True)
    with open(os.path.join(ROOT, args.out), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["line", "number", "status", "support", "sources", "problems"])
        w.writerows(rows)
    from collections import Counter
    c = Counter(r[2] for r in rows)
    print(f"{args.md}: {len(rows)} numbers checked (text before '{args.stop_at}')")
    for k in ("matched", "UNMATCHED", "UNRESOLVED", "NO-SOURCE"):
        print(f"  {k:10s} {c.get(k, 0)}")
    sup = Counter(r[3] for r in rows if r[2] == "matched")
    print("  matched by: " + ", ".join(f"{k} {sup[k]}" for k in ("key", "row", "calc", "file", "text") if sup[k])
          + "  (key/row/calc = the cited JSON key, TSV row or arithmetic; file = any value in a cited data file;"
            " text = any number in a cited Markdown/YAML document)")
    for r in rows:
        if r[2] != "matched":
            print(f"  line {r[0]:4d}  {r[2]:10s}  {r[1]!r:24s} {r[4][:90]}  {r[5][:120]}")
    print(f"report: {args.out}")
    return 1 if c.get("UNMATCHED", 0) else 0


if __name__ == "__main__":
    sys.exit(main())
