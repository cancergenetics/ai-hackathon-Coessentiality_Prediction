"""FZ-X freeze (DESIGN v3.1 §0 F4): family enumeration (m), unit sets, control units, H1 events, lockbox-testable list.

Everything here is computed from labels, the frozen FZ-C maps, Ensembl 116 / HGNC annotation and, for the
complete-case part of the registered evaluability rule only, the *missingness pattern* of A's gene-effect matrix
(``observed`` = a boolean A x gene frame; no effect value enters any function of this module). Lockbox (L) models
enter only as label counts (n_L / ref_L), with L = U minus A taken from the context table index.

* ``engine_membership``  FZ-C map files (clusters.tsv, membership.tsv, attachments.tsv, tiers) -> the FrozenMap
                         membership table used by the association engine (cmap.FrozenMap.from_membership). Attached
                         genes: level-1 attachments that are not locus-flagged (DECISIONS D3: only these may define
                         H1 events).
* ``evaluability``       the registered evaluability rule of H3a / H3b / BRIDGE-L, per (context, unit): >= 12
                         context and >= 12 reference A models after complete-case (context, covariates, score
                         observed) and >= 2 members after the context's removals, and a full-rank design (as in
                         assoc.run_association). m = number of evaluable tests.
* ``control_unit``       CTRL_* `control.complex_unit` from the entry's `complex_unit_rule` text.
* ``h1_events``          the H1 enumeration rule (h1.loss_calls + h1.enumerate_events) and the per-gene contexts.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

from part2 import covariates as covmod, h1, stats
from part2.assoc import _design
from part2.cmap import FrozenMap, symbol
from part2.cscore import MIN_FRAC, Annotation, Removal

NOT_EVALUABLE = "not evaluable at complex level"
MIN_N = 12
MIN_MEMBERS = 2
MIN_L = 5


# ------------------------------------------------------------------------------------------------ map conversion
def gene_key(sym, entrez) -> str:
    return f"{sym} ({int(entrez)})"


def engine_membership(clusters: pd.DataFrame, membership: pd.DataFrame, attachments: pd.DataFrame | None = None,
                      tiers: pd.DataFrame | None = None, column_of: dict | None = None) -> pd.DataFrame:
    """FZ-C files -> engine membership table (unit_id, level, parent, gene, role, locus, tier); IDs kept as frozen.
    gene = the DepMap gene-effect column key of the Entrez ID (column_of: Entrez -> 'SYMBOL (ENTREZ)', Part 1
    genes.tsv depmap_column); the FZ-C files carry Part 1 (HGNC) symbols, which differ for a few renamed genes."""
    cl = clusters.set_index("cluster_id")

    def key(sym, ent):
        if column_of is None:
            return gene_key(sym, ent)
        if int(ent) not in column_of:
            raise KeyError(f"Entrez {ent} ({sym}) has no DepMap column")
        return column_of[int(ent)]

    tier = tiers.set_index("cluster_id").tier if tiers is not None else pd.Series(dtype=str)

    def meta(u):
        par = cl.parent.get(u)
        return dict(level=int(cl.level[u]), parent="" if pd.isna(par) else str(par),
                    locus=bool(str(cl.locus[u]) == "True"), tier=str(tier.get(u, "")) if u in tier.index else "")

    rows = [dict(unit_id=r.cluster_id, gene=key(r.symbol, r.entrez), role="member", **meta(r.cluster_id))
            for r in membership.itertuples(index=False)]
    if attachments is not None and len(attachments):
        a = attachments[(~attachments.locus_attach.astype(str).eq("True"))]
        a = a[[int(cl.level[u]) == 1 for u in a.cluster]]
        mem = set(zip(membership.cluster_id, membership.symbol))
        for r in a.itertuples(index=False):
            if (r.cluster, r.symbol) in mem:
                raise ValueError(f"attached gene {r.symbol} is a member of {r.cluster}")
            rows.append(dict(unit_id=r.cluster, gene=key(r.symbol, r.entrez), role="attached", **meta(r.cluster)))
    out = pd.DataFrame(rows, columns=["unit_id", "level", "parent", "gene", "role", "locus", "tier"])
    if set(out.unit_id) != set(cl.index):
        raise ValueError("membership and clusters.tsv disagree on the unit set")
    sizes = out[out.role == "member"].groupby("unit_id").size()
    if not (sizes.reindex(cl.index) == cl["size"]).all():
        raise ValueError("member counts differ from clusters.tsv size")
    return out.sort_values(["unit_id", "role", "gene"], kind="mergesort").reset_index(drop=True)


# ------------------------------------------------------------------------------------------------ evaluability
def observed_score(observed: pd.DataFrame, genes, min_frac: float = MIN_FRAC) -> np.ndarray:
    """cscore.mean_score is non-NaN iff >= min_frac of its members are observed (and it has members)."""
    k = len(genes)
    if k == 0:
        return np.zeros(len(observed), bool)
    return observed.reindex(columns=list(genes), fill_value=False).values.sum(1) >= min_frac * k


class _DesignCache:
    """Full-rank / leverage check of the HC3 design on a row subset (labels only), cached by row mask."""

    def __init__(self, X):
        self.X, self.memo = X, {}

    def ok(self, rmask: np.ndarray) -> str:
        key = rmask.tobytes()
        if key not in self.memo:
            try:
                stats.HC3Design(self.X[rmask], 1)
                self.memo[key] = ""
            except (np.linalg.LinAlgError, ValueError) as e:
                self.memo[key] = f"design: {e}"
        return self.memo[key]


def evaluability(entry: dict, contexts_A: pd.DataFrame, cov_table: pd.DataFrame, fmap: FrozenMap, units,
                 annot: Annotation, observed: pd.DataFrame, removal_of, min_n: int = MIN_N,
                 min_members: int = MIN_MEMBERS, primary_only: bool = True) -> pd.DataFrame:
    """One row per (context, unit) of the entry's contexts x units with n, n_context, n_ref, evaluable, reason.
    Mirrors assoc.run_association (complete-case over context + covariates + score; >= min_n / min_n; HC3 design)
    and cscore.unit_scores (members after removals, 50% observed rule). removal_of(context) -> Removal."""
    items = entry["context_family"]["contexts"]
    ctxs = [c["id"] for c in items if not primary_only or c.get("role", "primary") == "primary"]
    rows = []
    for c in ctxs:
        names, dropped = covmod.covariates_for(entry, c)
        Z = covmod.assemble(names, cov_table.reindex(contexts_A.index))
        X, kept, ok = _design(contexts_A[c], Z)
        dc = _DesignCache(X)
        obs = observed.reindex(contexts_A.index).loc[ok]
        rem = removal_of(c)
        for u in units:
            mem = fmap.members[u]
            rm = set(rem.removed(mem, annot))
            use = [g for g in mem if g not in rm and g in observed.columns]
            base = dict(context=c, unit=u, n_members=len(mem), n_removed=len(rm), n_used=len(use),
                        covariates=";".join(names), covariates_dropped_by_rule=";".join(dropped))
            if len(use) < min_members:
                rows.append(dict(base, n=0, n_context=0, n_ref=0, evaluable=False,
                                 reason=f"< {min_members} members after removals"))
                continue
            rmask = observed_score(obs, use)
            n_ctx = int(X[rmask, 1].sum())
            n_ref = int(rmask.sum()) - n_ctx
            reason = ""
            if not rmask.any():
                reason = "no score"
            elif n_ctx < min_n or n_ref < min_n:
                reason = f"n_context {n_ctx} / n_ref {n_ref} below {min_n}/{min_n}"
            else:
                reason = dc.ok(rmask)
            rows.append(dict(base, n=int(rmask.sum()), n_context=n_ctx, n_ref=n_ref, evaluable=not reason,
                             reason=reason))
    return pd.DataFrame(rows)


def context_summary(ev: pd.DataFrame, n_units: int) -> pd.DataFrame:
    """Per context: complete-case n_context (max over units), evaluable / non-evaluable tests."""
    g = ev.groupby("context", sort=False)
    out = pd.DataFrame({"n_context_cc": g.n_context.max(), "n_ref_cc": g.n_ref.max(),
                        "n_tests": g.size(), "n_evaluable": g.evaluable.sum().astype(int)})
    out["n_not_evaluable"] = n_units - out.n_evaluable
    out["family_member"] = out.n_context_cc >= MIN_N
    return out.reset_index()


# ------------------------------------------------------------------------------------------------ controls
_CONTAINING = re.compile(r"frozen (coarse cluster|sub-module) containing ([A-Za-z0-9]+)"
                         r"(?: \(else (?:the one containing )?([A-Za-z0-9]+)\))?")


def _genes_in(text: str) -> list:
    return re.findall(r"\b[A-Z][A-Z0-9]+\b", text)


def parse_unit_rule(rule: str) -> dict:
    """complex_unit_rule text -> {'level': 1 | 2, 'targets': [primary, fallback...], 'partner': [...] or None}."""
    m = _CONTAINING.search(rule)
    if not m:
        raise ValueError(f"cannot parse complex_unit_rule: {rule!r}")
    level = 1 if m.group(1) == "coarse cluster" else 2
    targets = [m.group(2)] + ([m.group(3)] if m.group(3) else [])
    partner = None
    if level == 2:
        p = re.search(r"\) and (.*?) fall in different frozen sub-modules of the same frozen coarse complex", rule)
        if not p:
            raise ValueError(f"sub-module rule without a partner clause: {rule!r}")
        partner = _genes_in(p.group(1))
    return dict(level=level, targets=targets, partner=partner)


def control_unit(entry: dict, fmap: FrozenMap, annot: Annotation, removal: Removal | None,
                 member_resolved: dict | None = None) -> dict:
    """Frozen unit for a CTRL entry's complex-level readout (labels and map only).
    member_resolved: {'lost': tuple of symbols or 'chrY', 'target': symbol} for member_resolved readouts."""
    ctl = entry["control"]
    rule = parse_unit_rule(ctl["complex_unit_rule"])
    out = dict(entry=entry["id"], design_id=ctl["design_id"], score=entry["score"], rule_level=rule["level"],
               rule_targets=";".join(rule["targets"]), unit=NOT_EVALUABLE, structural_ok=False, reason="",
               target_used="", n_members=0, n_used=0, readout_members="")
    for t in rule["targets"]:
        us = fmap.units_of_symbol(t, rule["level"])
        if len(us) > 1:
            raise ValueError(f"{t} in several level-{rule['level']} units {us}")
        if us:
            out.update(unit=us[0], target_used=t)
            break
    if out["target_used"] == "":
        out["reason"] = f"{'/'.join(rule['targets'])} unclustered at level {rule['level']}"
        return out
    u = out["unit"]
    if rule["level"] == 2:
        parent = fmap.parent[u]
        partner_units = sorted({s for p in rule["partner"] for s in fmap.units_of_symbol(p, 2)
                                if fmap.parent.get(s) == parent and s != u})
        out["partner_units"] = ";".join(partner_units)
        if not partner_units:
            out.update(unit=NOT_EVALUABLE, reason=f"no {'/'.join(rule['partner'])} in a different sub-module of {parent}")
            return out
        rm = set(removal.removed(fmap.members[parent], annot)) if removal else set()
        s_mem = [g for g in fmap.members[u] if g not in rm]
        sib = [g for g in fmap.sibling_complement(u) if g not in rm]
        out.update(n_members=len(fmap.members[u]), n_used=len(s_mem), readout_members=";".join(s_mem),
                   sibling_members=";".join(sib))
        if len(s_mem) < MIN_MEMBERS or len(sib) < MIN_MEMBERS:
            out.update(unit=NOT_EVALUABLE, reason="< 2 members in sub-module or sibling complement after removals")
            return out
        out["structural_ok"] = True
        return out
    mem = fmap.members[u]
    rm = set(removal.removed(mem, annot)) if removal else set()
    use = [g for g in mem if g not in rm]
    out.update(n_members=len(mem), n_used=len(use), readout_members=";".join(use))
    if len(use) < MIN_MEMBERS:
        out.update(unit=NOT_EVALUABLE, reason="< 2 members after removals")
        return out
    if member_resolved is not None:
        lost, target = member_resolved["lost"], member_resolved["target"]
        if lost == "chrY":
            P = [g for g in use if annot.x_gametolog(symbol(g))]
        else:
            P = [g for g in use if any(annot.is_paralog(symbol(g), l) for l in lost)]
        out["readout_members"] = ";".join(P)
        if target not in {symbol(g) for g in P}:
            out.update(unit=NOT_EVALUABLE, reason=f"paralog component lacks {target}")
            return out
    out["structural_ok"] = True
    return out


# ------------------------------------------------------------------------------------------------ H1
def h1_contexts(events, loss: pd.DataFrame, kind: pd.DataFrame, fold: pd.Series) -> pd.DataFrame:
    """One registered context per lost gene with >= 1 event: counts, type, covariates (covariate-of-context rule)."""
    rows = []
    for g in sorted({e.gene for e in events}):
        evs = [e for e in events if e.gene == g]
        e0 = evs[0]
        cid = f"H1loss:{g}"
        bad = covmod.context_functions(cid)
        covs = [c for c in h1.COVARIATES[e0.type] if c not in bad]
        isA = fold.reindex(loss.index).eq("A").values
        k = kind[g][(loss[g].values == 1) & isA]
        rows.append(dict(context=cid, gene=g, type=e0.type, n_A=e0.n_A, ref_A=e0.ref_A, n_L=e0.n_L, ref_L=e0.ref_L,
                         lockbox_testable=bool(e0.testable), n_A_LoF=int((k == "LoF").sum()),
                         n_A_homdel=int((k == "homdel").sum()), n_A_silencing=int((k == "silencing").sum()),
                         covariates=";".join(covs), events=";".join(e.id for e in evs)))
    return pd.DataFrame(rows)


def h1_candidates(loss: pd.DataFrame, kind: pd.DataFrame, fold: pd.Series, events, min_A: int = MIN_N) -> pd.DataFrame:
    """Every linked gene with >= 1 loss model in A: counts and why it does / does not give an event."""
    isA = fold.reindex(loss.index).eq("A").values
    isL = fold.reindex(loss.index).eq("L").values
    by_gene = {}
    for e in events:
        by_gene.setdefault(e.gene, []).append(e.id)
    rows = []
    for g in loss.columns:
        v = loss[g].values
        nA, rA = int(((v == 1) & isA).sum()), int(((v == 0) & isA).sum())
        if nA == 0:
            continue
        k = kind[g][(v == 1) & isA]
        cand = nA >= min_A and rA >= min_A
        rows.append(dict(gene=g, n_A=nA, ref_A=rA, n_L=int(((v == 1) & isL).sum()), ref_L=int(((v == 0) & isL).sum()),
                         n_A_LoF=int((k == "LoF").sum()), n_A_homdel=int((k == "homdel").sum()),
                         n_A_silencing=int((k == "silencing").sum()), candidate=cand,
                         events=";".join(by_gene.get(g, [])),
                         outcome=("event" if g in by_gene else
                                  "candidate, no evaluable axis (|N| < 2 and |P| = 0) in any linked complex" if cand
                                  else f"< {min_A} loss or reference models in A")))
    return pd.DataFrame(rows).sort_values(["candidate", "n_A"], ascending=[False, False]).reset_index(drop=True)


# ------------------------------------------------------------------------------------------------ registry fill
def yaml_scalar(v) -> str:
    """Render a scalar in the registry YAML subset (strings double-quoted, ASCII only)."""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        return repr(float(v))
    s = str(v)
    if not s.isascii():
        raise ValueError(f"non-ASCII registry value {s!r}")
    if re.fullmatch(r"[0-9a-f]{16,64}", s) or re.fullmatch(r"[A-Za-z0-9_./-]+", s) and not re.fullmatch(
            r"[-+]?[0-9.eE]+|true|false|null|FZ-[CX]", s):
        return s
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def yaml_block(v, indent: int) -> list:
    """Lines for a block value (mapping / list of mappings / flow list of scalars) at `indent`."""
    pad = " " * indent
    if isinstance(v, dict):
        out = []
        for k, x in v.items():
            if isinstance(x, (dict,)) or (isinstance(x, list) and x and isinstance(x[0], dict)):
                out.append(f"{pad}{k}:")
                out += yaml_block(x, indent + 2)
            elif isinstance(x, list):
                out.append(f"{pad}{k}: [" + ", ".join(yaml_scalar(y) for y in x) + "]")
            else:
                out.append(f"{pad}{k}: {yaml_scalar(x)}")
        return out
    if isinstance(v, list):
        out = []
        for item in v:
            sub = yaml_block(item, indent + 2)
            out.append(pad + "- " + sub[0].lstrip())
            out += sub[1:]
        return out
    raise TypeError(type(v))


def fill_text(text: str, fills: dict) -> str:
    """Replace the freeze placeholders of a registry entry's text by `fills` ({pending path: value}) and drop the
    `pending` block; everything else (comments, order, wording) is kept byte for byte. Supported paths: top-level
    scalars, '<block>.<key>' one level deep, 'frozen_inputs[i]', 'context_family.contexts' (list of mappings)."""
    lines = text.split("\n")
    out, top, item_i, done = [], None, -1, set()
    skip_pending = False
    for line in lines:
        if skip_pending:
            if line.startswith("  - ") or line.strip() == "":
                if line.strip() == "":
                    skip_pending = False
                    out.append(line)
                continue
            skip_pending = False
        m_top = re.match(r"^([A-Za-z_][A-Za-z0-9_]*):(.*)$", line)
        if m_top:
            top, item_i = m_top.group(1), -1
            if top == "pending":
                skip_pending = True
                continue
            if top in fills and re.match(r"^\s*FZ-[CX]\s*(#.*)?$", m_top.group(2)):
                out.append(f"{top}: {yaml_scalar(fills[top])}")
                done.add(top)
                continue
        m_sub = re.match(r"^  ([A-Za-z_][A-Za-z0-9_]*):\s*FZ-[CX]\s*(#.*)?$", line)
        if m_sub and top is not None and f"{top}.{m_sub.group(1)}" in fills:
            path = f"{top}.{m_sub.group(1)}"
            v = fills[path]
            if isinstance(v, (dict, list)):
                out.append(f"  {m_sub.group(1)}:")
                out += yaml_block(v, 4)
            else:
                out.append(f"  {m_sub.group(1)}: {yaml_scalar(v)}")
            done.add(path)
            continue
        if top == "frozen_inputs" and line.startswith("  - "):
            item_i += 1
            path = f"frozen_inputs[{item_i}]"
            m_it = re.match(r"^  - FZ-[CX](\s+#.*)?$", line)
            if m_it and path in fills:
                out.append(f"  - {fills[path]}" + (m_it.group(1) or ""))
                done.add(path)
                continue
        out.append(line)
    missing = set(fills) - done
    if missing:
        raise ValueError(f"placeholders not found for {sorted(missing)}")
    return "\n".join(out)


# ------------------------------------------------------------------------------------------------ amendments (text)
def _block_end(lines, i, indent):
    j = i + 1
    while j < len(lines):
        l = lines[j]
        if l.strip() and (len(l) - len(l.lstrip(" "))) <= indent:
            break
        j += 1
    return j


def replace_top(text: str, key: str, new_lines: list) -> str:
    """Replace the top-level field `key` (its line and indented continuation) by new_lines."""
    lines = text.split("\n")
    i = next((k for k, l in enumerate(lines) if l.startswith(f"{key}:")), None)
    if i is None:
        raise ValueError(f"no top-level field {key}")
    return "\n".join(lines[:i] + list(new_lines) + lines[_block_end(lines, i, 0):])


def replace_sub(text: str, top: str, key: str, new_lines: list) -> str:
    """Replace `top.key` (two-space indented field and its continuation) by new_lines."""
    lines = text.split("\n")
    i = next((k for k, l in enumerate(lines) if l.startswith(f"{top}:")), None)
    if i is None:
        raise ValueError(f"no top-level field {top}")
    end = _block_end(lines, i, 0)
    s = next((k for k in range(i + 1, end) if lines[k].startswith(f"  {key}:")), None)
    if s is None:
        raise ValueError(f"no field {top}.{key}")
    return "\n".join(lines[:s] + list(new_lines) + lines[_block_end(lines, s, 2):])


def insert_before_top(text: str, key: str, new_lines: list) -> str:
    lines = text.split("\n")
    i = next(k for k, l in enumerate(lines) if l.startswith(f"{key}:"))
    return "\n".join(lines[:i] + list(new_lines) + lines[i:])


def folded(key: str, text: str, indent: int = 0, width: int = 116) -> list:
    """`key: >` block scalar lines (ASCII), wrapped at `width`."""
    import textwrap
    pad = " " * (indent + 2)
    if not text.isascii():
        raise ValueError("non-ASCII registry text")
    return [" " * indent + f"{key}: >"] + textwrap.wrap(text, width=width, initial_indent=pad, subsequent_indent=pad)
