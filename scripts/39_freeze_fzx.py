#!/usr/bin/env python
"""FZ-X freeze (DESIGN v3.1 §0 F4): unit sets + complex-set hashes, family enumeration (m per family), control units,
H1 event list, lockbox-testable list; fills the registry placeholders (fill-only edits) and writes the manifest.

Inputs (all verified against their committed manifests before use):
  FZ-C map files (results/part2/construction/FZC_manifest.json), tiers (results/part2/construction/tiers/),
  the context table (results/part2/contexts/checksums.json), A (results/part2/lockbox/manifest.json), the omics
  label matrices cached by scripts/32 (LoF, hom-del, damaging, expression), Ensembl 116 paralogs / coordinates, HGNC.
Dependency data: NONE of its values. The registered evaluability rule needs complete-case n per (context, unit), so
the *missingness pattern* of A's gene-effect matrix is read (assoc.load_A_effects -> notna()) and the values are
discarded at once. Lockbox models enter only as label counts (L = U minus A over the context table index).

Outputs: results/part2/fzx/ (lists + FZX_manifest.json + README.md), data_cache/part2/fzx/h1_loss_table.parquet
(H1 loss columns over U; gitignored, SHA-256 in the manifest), filled docs/registry/*.yaml.

  PYTHONHASHSEED=0 PART1_DIR=... .venv/bin/python scripts/39_freeze_fzx.py
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import assoc, contexts as cx, controls, data, covariates as covmod, cscore, fzx, h1, h2, paths, registry  # noqa: E402
from part2.cmap import FrozenMap, symbol  # noqa: E402

CON = paths.RESULTS / "construction"
FZC_MANIFEST = CON / "FZC_manifest.json"
OUT = paths.RESULTS / "fzx"
CACHE_OUT = paths.CACHE / "fzx"
MANIFEST = OUT / "FZX_manifest.json"


def sha(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def rel(p) -> str:
    return str(Path(p).resolve().relative_to(ROOT.resolve())) if Path(p).resolve().is_relative_to(ROOT.resolve()) \
        else str(Path(p).relative_to(paths.PROJECT))


def verify_fzc() -> dict:
    man = json.loads(FZC_MANIFEST.read_text())["sha256"]
    for f in ("clusters.tsv", "membership.tsv", "attachments.tsv"):
        p = f"results/part2/construction/zavg_map/{f}"
        if sha(ROOT / p) != man[p]:
            raise SystemExit(f"{p} differs from the FZ-C manifest")
    for f in ("clusters.tsv", "membership.tsv"):
        p = f"results/part2/construction/gls_map/{f}"
        if sha(ROOT / p) != man[p]:
            raise SystemExit(f"{p} differs from the FZ-C manifest")
    r = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "--", "results/part2/construction"],
                       capture_output=True, text=True)
    if r.stdout.strip():
        raise SystemExit("results/part2/construction has uncommitted changes")
    return man


def verify_context_table() -> pd.DataFrame:
    ck = json.loads((paths.CONTEXT_RESULTS / "checksums.json").read_text())["files"]
    p = paths.CONTEXT_CACHE / "context_table.parquet"
    if sha(p) != ck["data_cache/part2/contexts/context_table.parquet"]:
        raise SystemExit("context table differs from results/part2/contexts/checksums.json")
    for k in ("lof_matrix.parquet", "cn_homdel.parquet"):
        if sha(paths.CONTEXT_CACHE / k) != ck[f"data_cache/part2/contexts/{k}"]:
            raise SystemExit(f"{k} differs from checksums.json")
    return cx.load_table(p)


def load_map(name: str, attachments: bool) -> tuple[pd.DataFrame, FrozenMap]:
    d = CON / f"{name}_map"
    cl = pd.read_csv(d / "clusters.tsv", sep="\t")
    mem = pd.read_csv(d / "membership.tsv", sep="\t")
    att = pd.read_csv(d / "attachments.tsv", sep="\t") if attachments else None
    tiers = pd.read_csv(CON / "tiers" / f"{name}_tiers.tsv", sep="\t")[["cluster_id", "tier"]]
    g = data.protein_genes()
    em = fzx.engine_membership(cl, mem, att, tiers, column_of=dict(zip(g.entrez.astype(int), g.depmap_column)))
    return em, FrozenMap.from_membership(em)


def homdel_shares(ctx: pd.DataFrame, A, ctx_ids) -> pd.DataFrame:
    """LoF context -> fraction of its positive A lines that are hom-del calls of the gene (cn_homdel cache)."""
    hd = pd.read_parquet(paths.CONTEXT_CACHE / "cn_homdel.parquet")
    sym2e = omics_sym2e()
    rows = []
    for c in ctx_ids:
        if not c.startswith(("LoF:", "LoF_inclMSI:", "PTEN_LoF")):
            continue
        g = "PTEN" if c.startswith("PTEN_LoF") else c.split(":", 1)[1]
        ent = sym2e.get(g)
        pos = ctx.loc[A, c] == 1
        col = str(ent) if ent is not None else None
        share = float((hd.loc[A, col][pos.values] == 1).mean()) if col in hd.columns and pos.any() else np.nan
        rows.append(dict(context=c, gene=g, n_pos_A=int(pos.sum()), homdel_share=share,
                         window_removal=bool(share >= 0.5) if np.isfinite(share) else False))
    return pd.DataFrame(rows)


def fzr_draft(p: Path) -> str:
    """The committed FZ-R draft of a registry entry: the newest committed version that still has `pending:`
    (HEAD before the FZ-X commit; its parent afterwards, so the script can be re-run to verify the freeze)."""
    r = rel(p)
    log = subprocess.run(["git", "-C", str(ROOT), "log", "--format=%H", "--", r], capture_output=True, text=True,
                         check=True).stdout.split()
    for c in log:
        t = subprocess.run(["git", "-C", str(ROOT), "show", f"{c}:{r}"], capture_output=True, text=True)
        if t.returncode == 0 and "\npending:" in t.stdout:
            return t.stdout
    raise SystemExit(f"no committed FZ-R draft (with pending) of {r}")


def omics_sym2e() -> dict:
    """Symbol -> Entrez from the column keys of the damaging-mutation matrix (as in contexts.build)."""
    df, _ = cx.read_omics(paths.OMICS_DAMAGING, [])
    out = {}
    for c in df.columns:
        s, e = cx._entrez_of(c)
        if e is not None:
            out.setdefault(s, e)
    return out


def omics_by_symbol(path, U, ent_of: dict) -> pd.DataFrame:
    df, _ = cx.read_omics(path, U)
    e = cx._by_entrez(df).reindex(U)
    cols = {s: e[k] if k in e.columns else pd.Series(np.nan, index=e.index) for s, k in ent_of.items()}
    return pd.DataFrame(cols, index=U).astype(float)


def h1_build(ctx, fold, fmap: FrozenMap, annot):
    """H1 enumeration rule 1-8 (labels, frozen map, Ensembl 116, coordinates)."""
    U = list(ctx.index)
    genes = sorted({g for u in fmap.coarse() for g in fmap.members[u]} |
                   {g for u, att in fmap.attached.items() for g in att})
    ent_of = {}
    for k in genes:
        s, e = symbol(k), int(k.split("(")[-1].rstrip(")"))
        if s in ent_of and ent_of[s] != e:
            raise SystemExit(f"symbol {s} with two Entrez IDs")
        ent_of[s] = e
    lof_all = pd.read_parquet(paths.CONTEXT_CACHE / "lof_matrix.parquet")
    hd_all = pd.read_parquet(paths.CONTEXT_CACHE / "cn_homdel.parquet")

    def by_sym(df):
        return pd.DataFrame({s: (df[str(e)] if str(e) in df.columns else pd.Series(np.nan, index=df.index))
                             for s, e in ent_of.items()}).reindex(U).astype(float)

    lof_incl = by_sym(lof_all)
    msi = ctx["MSI_H"].reindex(U).values
    lof = lof_incl.copy()
    for s in lof.columns:
        if s not in cx.MMR_GENES:
            lof[s] = np.where(msi == 0, lof_incl[s].values, np.nan)
    homdel = by_sym(hd_all)
    dam = omics_by_symbol(paths.OMICS_DAMAGING, U, ent_of)
    expr = omics_by_symbol(paths.OMICS_EXPR, U, ent_of)
    # LoF functional check: the FZ-LB-era table, plus the same check for linked genes with MSS LoF lines not in it
    func = pd.read_csv(paths.CONTEXT_RESULTS / "lof_functional_check.tsv", sep="\t")
    in_tab = set(func.entrez.astype(int))
    mss = (msi == 0)
    need = [s for s in lof.columns if ent_of[s] not in in_tab and np.nansum((lof_incl[s].values == 1) & mss) > 0]
    extra = pd.DataFrame()
    if need:
        ex_df, _ = cx.read_omics(paths.OMICS_EXPR, U)
        ex_e = cx._by_entrez(ex_df).reindex(U)
        lof_sub = pd.DataFrame({ent_of[s]: lof_incl[s].values for s in need}, index=pd.Index(U, name="ModelID"))
        ann = pd.DataFrame({"symbol": need, "entrez": [ent_of[s] for s in need]},
                           index=pd.Index([ent_of[s] for s in need], dtype=np.int64))
        sym2e = pd.Series([ent_of[s] for s in need], index=need)
        extra = cx.lof_functional_check(lof_sub, ctx.reindex(U), ex_e, sym2e, ann, min_lof=0)
        extra["source"] = "FZ-X (same check, genes not in lof_functional_check.tsv)"
    passed = set(func.loc[func["pass"].astype(bool), "entrez"].astype(int))
    if len(extra):
        passed |= set(extra.loc[extra["pass"].astype(bool), "entrez"].astype(int))
    lof_pass = {s for s in lof.columns if ent_of[s] in passed}
    lineage = ctx["OncotreeLineage"].reindex(U)
    loss, kind = h1.loss_calls(lof, homdel, dam, expr, lineage, lof_pass=lof_pass)
    events = h1.enumerate_events(fmap, loss, kind, annot, fold)
    return dict(loss=loss, kind=kind, events=events, extra_check=extra, ent_of=ent_of, n_genes=len(ent_of),
                n_lof_pass=len(lof_pass))


def control_rows(entries, fmap, annot, ctx, A, cov_table, observed, shares) -> pd.DataFrame:
    rows = []
    for k, e in sorted(entries.items()):
        if e["direction"] != "controls":
            continue
        ctl = e["control"]
        ctxs = [c["id"] for c in e["context_family"]["contexts"] if c.get("role", "primary") == "primary"]
        cid = ctxs[0]
        base = dict(entry=k, design_id=ctl["design_id"], context=cid, kill_direction=e["kill_rule"]["direction"],
                    level=ctl["level"], gene_units=";".join(ctl["gene_units"]))
        if ctl["level"] == "G":
            rows.append(dict(base, unit="", structural_ok=False, reason="gene level only", complex_evaluable=False,
                             family_m=len(ctl["gene_units"]),
                             complex_set_hash=registry.set_hash(ctl["gene_units"])))
            continue
        spec = controls.CONTROL_SPECS[ctl["design_id"]]
        sh = shares.get(cid)
        rem = controls._removal(spec, sh)
        mr = dict(lost=spec["lost"], target=spec["target"]) if e["score"] == "member_resolved" else None
        cu = fzx.control_unit(e, fmap, annot, rem, mr)
        n = n_ctx = n_ref = 0
        ev_reason = cu["reason"]
        if cu["structural_ok"]:
            names, _ = covmod.covariates_for(e, cid)
            Z = covmod.assemble(names, cov_table.reindex(A))
            X, _, ok = assoc._design(ctx.loc[A, cid], Z)
            obs = observed.reindex(A).loc[ok]
            rmask = fzx.observed_score(obs, cu["readout_members"].split(";"))
            if e["score"] == "submodule_vs_sibling":
                rmask &= fzx.observed_score(obs, cu["sibling_members"].split(";"))
            n, n_ctx = int(rmask.sum()), int(X[rmask, 1].sum())
            n_ref = n - n_ctx
            ev_reason = "" if (n_ctx >= fzx.MIN_N and n_ref >= fzx.MIN_N) else \
                f"n_context {n_ctx} / n_ref {n_ref} below 12/12 after complete-case"
            if not ev_reason:
                ev_reason = fzx._DesignCache(X).ok(rmask)
        complex_ev = bool(cu["structural_ok"] and not ev_reason)
        units = list(ctl["gene_units"]) + ([cu["unit"]] if cu["unit"] in fmap.members else [])
        rows.append(dict(base, unit=cu["unit"], target_used=cu["target_used"], rule_targets=cu["rule_targets"],
                         n_members=cu["n_members"], n_used=cu["n_used"], readout_members=cu["readout_members"],
                         partner_units=cu.get("partner_units", ""), structural_ok=cu["structural_ok"],
                         reason=ev_reason, n_cc=n, n_context_cc=n_ctx, n_ref_cc=n_ref, complex_evaluable=complex_ev,
                         family_m=len(ctl["gene_units"]) + int(complex_ev), complex_set_hash=registry.set_hash(units)))
    return pd.DataFrame(rows)


def lockbox_rows(entries, h1ctx: pd.DataFrame, h3b_drop: set) -> pd.DataFrame:
    rows = []
    for k, e in sorted(entries.items()):
        cs = e["context_family"]["contexts"]
        if e["direction"] == "H1":                       # one row per frozen event (the contexts are per gene)
            for r in h1ctx.itertuples(index=False):
                for ev in r.events.split(";"):
                    rows.append(dict(entry=k, direction=e["direction"], context=r.context, event=ev, role="primary",
                                     n_L=r.n_L, ref_L=r.ref_L, lockbox_testable=r.lockbox_testable,
                                     route="lockbox" if r.lockbox_testable else "rare-context route (LOSO + DEMETER2 "
                                     "on A + leave-one-model-out), labelled 'not lockbox-testable'"))
            continue
        for c in cs:
            role = c.get("role", "primary")
            if role != "primary":
                continue
            t = bool(c.get("n_L", 0) >= fzx.MIN_L and c.get("ref_L", 0) >= fzx.MIN_L)
            in_family = not (e["direction"] in ("H3b",) and c["id"] in h3b_drop)
            rows.append(dict(entry=k, direction=e["direction"], context=c["id"], event="", role=role,
                             n_L=c.get("n_L"), ref_L=c.get("ref_L"), lockbox_testable=t,
                             route=("lockbox" if t else "rare-context route (LOSO + DEMETER2 on A + "
                                    "leave-one-model-out), labelled 'not lockbox-testable'")
                             + ("" if in_family else "; context left the family after complete-case")))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-registry", action="store_true", help="compute and write lists only (no registry fill)")
    a = ap.parse_args()
    t0 = time.time()
    timings = {}
    if MANIFEST.exists() and json.loads(MANIFEST.read_text()).get("amendments"):
        raise SystemExit("FZ-X has been amended (docs/registry/AMENDMENTS.md); re-running would overwrite the amended "
                         "lists and entries. To reproduce FZ-X, run this script on a checkout of the FZ-X commit; "
                         "amendments are re-applied by their own scripts (e.g. scripts/40_amend_A1.py).")
    fzc = verify_fzc()
    ctx = verify_context_table()
    A = assoc.a_models()
    U = list(ctx.index)
    if not set(A) <= set(U):
        raise SystemExit("A is not a subset of the context table index (U)")
    fold = pd.Series(np.where(pd.Index(U).isin(A), "A", "L"), index=U)
    em_z, fz = load_map("zavg", attachments=True)
    em_g, fg = load_map("gls", attachments=False)
    timings["maps"] = time.time() - t0

    keys = sorted(set(em_z.gene) | set(em_g.gene))
    # ---- missingness pattern of A's gene effects (values discarded immediately)
    E = assoc.load_A_effects()
    observed = E.notna()
    all_cols = list(E.columns)
    del E
    missing_cols = [k for k in keys if k not in all_cols]
    if missing_cols:
        raise SystemExit(f"map genes without a gene-effect column: {missing_cols}")
    observed = observed.reindex(A)
    timings["observed_mask"] = time.time() - t0
    annot = cscore.load_annotation(all_cols)          # as scripts/34_run_controls.py: every gene-effect column
    timings["annotation"] = time.time() - t0

    fzr_text = {p.stem: fzr_draft(p) for p in sorted(paths.REGISTRY_DIR.glob("*.yaml"))}
    entries = {k: registry.load_yaml(t) for k, t in fzr_text.items()}   # the committed FZ-R drafts (pending)
    for k, e in entries.items():
        errs = registry.validate(e, k)
        if errs:
            raise SystemExit(f"{k}: FZ-R draft does not validate: {errs}")
    cov_table = covmod.source_table(ctx, A)
    ctxA = ctx.reindex(A)

    # ---- unit sets
    coarse_z, coarse_g = fz.coarse(), fg.coarse()
    h2u = h2.h2_units(fz)
    OUT.mkdir(parents=True, exist_ok=True)
    CACHE_OUT.mkdir(parents=True, exist_ok=True)
    em_z.to_csv(OUT / "zavg_engine_membership.tsv", sep="\t", index=False)
    em_g.to_csv(OUT / "gls_engine_membership.tsv", sep="\t", index=False)
    (OUT / "units_zavg_coarse.txt").write_text("\n".join(coarse_z) + "\n")
    (OUT / "units_gls_coarse.txt").write_text("\n".join(coarse_g) + "\n")
    h2tab = pd.DataFrame([dict(unit=s, parent=fz.parent[s], n_members=len(fz.members[s]),
                               n_siblings_in_parent=len(fz.submodules(fz.parent[s])),
                               parent_size=len(fz.members[fz.parent[s]]), locus=fz.locus.get(s, False),
                               tier=fz.tier.get(s, ""), members=";".join(symbol(g) for g in fz.members[s]))
                          for s in h2u])
    h2tab.to_csv(OUT / "h2_submodules.tsv", sep="\t", index=False)

    # ---- removals that depend on hom-del share (LoF contexts)
    all_ctx_ids = sorted({c["id"] for e in entries.values() if isinstance(e["context_family"]["contexts"], list)
                          for c in e["context_family"]["contexts"]})
    sh = homdel_shares(ctx, A, [c for c in all_ctx_ids if c in ctx.columns])
    sh.to_csv(OUT / "lof_homdel_shares.tsv", sep="\t", index=False)
    shares = dict(zip(sh.context, sh.homdel_share))

    def removal_of(c):
        return cscore.removal_for_context(c, shares.get(c))

    # ---- evaluability and m (H3a, H3b on the zavg map; BRIDGE-L on the GLS map)
    evs, summ = {}, {}
    for eid, fmap, units in (("H3a_lineage_atlas", fz, coarse_z), ("H3b_driver_family", fz, coarse_z),
                             ("BRIDGE-L_confirmation_rate", fg, coarse_g)):
        ev = fzx.evaluability(entries[eid], ctxA, cov_table, fmap, units, annot, observed, removal_of)
        evs[eid] = ev
        summ[eid] = fzx.context_summary(ev, len(units))
    timings["evaluability"] = time.time() - t0
    s3b = summ["H3b_driver_family"]
    h3b_drop = set(s3b.context[~s3b.family_member])
    ev_all = pd.concat([v.assign(entry=k) for k, v in evs.items()], ignore_index=True)
    ev_all[["entry", "context", "unit", "n_members", "n_removed", "n_used", "n", "n_context", "n_ref", "evaluable",
            "reason", "covariates", "covariates_dropped_by_rule"]].to_csv(OUT / "evaluability_tests.tsv.gz",
                                                                          sep="\t", index=False,
                                                                          compression={"method": "gzip", "mtime": 0})

    def m_of(eid, drop=()):
        ev = evs[eid]
        return int(ev.evaluable[~ev.context.isin(drop)].sum())

    m = {"H3a_lineage_atlas": m_of("H3a_lineage_atlas"),
         "H3b_driver_family": m_of("H3b_driver_family", h3b_drop)}
    br = evs["BRIDGE-L_confirmation_rate"]
    h3a_ids = {c["id"] for c in entries["H3a_lineage_atlas"]["context_family"]["contexts"]}
    m_bridge = {"H3a": int(br.evaluable[br.context.isin(h3a_ids)].sum()),
                "H3b": int(br.evaluable[~br.context.isin(h3a_ids) & ~br.context.isin(h3b_drop)].sum())}
    m_bridge["total"] = m_bridge["H3a"] + m_bridge["H3b"]

    # ---- context tables (final lists with n_A / n_L and the evaluability outcome)
    def ctx_table(eid):
        items = [c for c in entries[eid]["context_family"]["contexts"] if c.get("role", "primary") == "primary"]
        s = summ[eid].set_index("context")
        rows = []
        for c in items:
            r = s.loc[c["id"]]
            rows.append(dict(context=c["id"], definition=c["definition"], n_A=c.get("n_A"), ref_A=c.get("ref_A"),
                             n_L=c.get("n_L"), ref_L=c.get("ref_L"), lockbox_testable=c.get("lockbox_testable"),
                             n_context_cc=int(r.n_context_cc), n_ref_cc=int(r.n_ref_cc), n_tests=int(r.n_tests),
                             n_evaluable=int(r.n_evaluable), n_not_evaluable=int(r.n_not_evaluable),
                             in_family=bool(r.family_member),
                             outcome=("in family" if r.family_member else
                                      "left the family: < 12 A context models after complete-case (exploratory)")))
        return pd.DataFrame(rows)

    t3a, t3b = ctx_table("H3a_lineage_atlas"), ctx_table("H3b_driver_family")
    t3a.to_csv(OUT / "h3a_contexts.tsv", sep="\t", index=False)
    t3b_ex = pd.DataFrame([dict(context=x["id"], n_A=x["n_A"], outcome="excluded at FZ-R: " + x["reason"])
                           for x in entries["H3b_driver_family"]["context_family"].get("excluded", [])])
    pd.concat([t3b, t3b_ex], ignore_index=True).to_csv(OUT / "h3b_contexts.tsv", sep="\t", index=False)
    tbr = ctx_table("BRIDGE-L_confirmation_rate")
    tbr.to_csv(OUT / "bridgeL_gls_contexts.tsv", sep="\t", index=False)

    # ---- controls
    ctab = control_rows(entries, fz, annot, ctx, A, cov_table, observed, shares)
    ctab.to_csv(OUT / "control_units.tsv", sep="\t", index=False)

    # ---- H1 events (labels only)
    hb = h1_build(ctx, fold, fz, annot)
    events = hb["events"]
    etab = h1.events_table(events)
    etab["axes"] = [";".join(e.axes) for e in events]
    etab.to_csv(OUT / "h1_events.tsv", sep="\t", index=False)
    h1ctx = fzx.h1_contexts(events, hb["loss"], hb["kind"], fold)
    h1ctx.to_csv(OUT / "h1_contexts.tsv", sep="\t", index=False)
    cand_tab = fzx.h1_candidates(hb["loss"], hb["kind"], fold, events)
    cand_tab.to_csv(OUT / "h1_candidates.tsv", sep="\t", index=False)
    if len(hb["extra_check"]):
        hb["extra_check"].to_csv(OUT / "h1_lof_check_extra.tsv", sep="\t", index=False)
    h1_units = h1.units_of(events)
    (OUT / "units_h1_components.txt").write_text("\n".join(h1_units) + ("\n" if h1_units else ""))
    genes_ev = sorted({e.gene for e in events})
    loss_tab = pd.DataFrame({f"H1loss:{g}": hb["loss"][g] for g in genes_ev}, index=pd.Index(U, name="ModelID"))
    loss_tab.to_parquet(CACHE_OUT / "h1_loss_table.parquet")
    timings["h1"] = time.time() - t0

    # ---- lockbox-testable list
    lb = lockbox_rows(entries, h1ctx, h3b_drop)
    lb.to_csv(OUT / "lockbox_testable.tsv", sep="\t", index=False)

    # ---- registry fills
    hashes = dict(zavg_coarse=registry.set_hash(coarse_z), gls_coarse=registry.set_hash(coarse_g),
                  h2_submodules=registry.set_hash(h2u), h1_components=registry.set_hash(h1_units))
    fzc_rel, fzx_rel = rel(FZC_MANIFEST), rel(MANIFEST)
    fills = {}
    for k, e in entries.items():
        pend = list(e.get("pending", []))
        if not pend:
            continue
        f = {}
        for p in pend:
            if p.startswith("frozen_inputs["):
                i = int(p[len("frozen_inputs["):-1])
                cur = e["frozen_inputs"][i]
                f[p] = fzc_rel if cur == "FZ-C" else fzx_rel
        d = e["direction"]
        if d == "controls":
            r = ctab.set_index("entry").loc[k]
            if "control.complex_unit" in pend:
                f["control.complex_unit"] = r.unit
            if "complex_set_hash" in pend:
                f["complex_set_hash"] = r.complex_set_hash
            if "family_m.value" in pend:
                f["family_m.value"] = int(r.family_m)
        elif k in ("H3a_lineage_atlas", "H3b_driver_family"):
            f.update({"complex_set_hash": hashes["zavg_coarse"], "family_m.n_units": len(coarse_z),
                      "family_m.value": m[k]})
        elif k == "H2_submodule_divergence":
            f["complex_set_hash"] = hashes["h2_submodules"]
        elif k == "BENCH_complex_vs_gene":
            f.update({"complex_set_hash": hashes["zavg_coarse"],
                      "family_m.value": "set by the committed H3a and H3b discovery runs (estimation, no test family)"})
        elif k == "BRIDGE-L_confirmation_rate":
            f.update({"complex_set_hash": hashes["gls_coarse"], "family_m.value": dict(m_bridge)})
        elif k == "H1_within_complex_GI":
            items = []
            for r in h1ctx.itertuples(index=False):
                items.append(dict(
                    id=r.context,
                    definition=(f"loss({r.gene}) = 1 if LoF (functional check passed) / hom-del / silencing by the "
                                f"enumeration rule; 0 = reference (no damaging mutation, no hom-del, expressed); NaN "
                                f"otherwise; {r.type}-type event"),
                    sign="two-sided", type=r.type, covariates=r.covariates.split(";"),
                    n_A=int(r.n_A), ref_A=int(r.ref_A), n_L=int(r.n_L), ref_L=int(r.ref_L),
                    lockbox_testable=bool(r.lockbox_testable), events=r.events.split(";")))
            f.update({"context_family.contexts": items, "complex_set_hash": hashes["h1_components"],
                      "family_m.value": len(h1_units)})
        if set(f) != set(pend):
            raise SystemExit(f"{k}: fills {sorted(f)} != pending {sorted(pend)}")
        fills[k] = f

    reg_report = {}
    if not a.no_registry:
        for k, f in sorted(fills.items()):
            p = paths.REGISTRY_DIR / f"{k}.yaml"
            old_txt = fzr_text[k]
            old = registry.load_yaml(old_txt)
            new_txt = fzx.fill_text(old_txt, f)
            new = registry.load_yaml(new_txt)
            errs = registry.validate(new, k)
            bad = registry.fill_only_changes(old, new)
            if errs or bad or registry.placeholders(new):
                raise SystemExit(f"{k}: fill failed validation {errs} / non-fill changes {bad}")
            p.write_text(new_txt)
            reg_report[k] = dict(filled=sorted(f), fill_only_changes=[], sha256=sha(p))
        if set(reg_report) != set(entries):
            raise SystemExit(f"entries without a fill: {sorted(set(entries) - set(reg_report))}")

    # ---- summary + manifest
    summary = dict(
        n_U=len(U), n_A=len(A), n_L=int((fold == "L").sum()),
        zavg=dict(n_coarse=len(coarse_z), n_submodules=len(fz.submodules()), n_h2_units=len(h2u),
                  n_h2_parents=len(h2.eligible_parents(fz)), n_attached_level1_nonlocus=int((em_z.role == "attached").sum()),
                  complex_set_hash=hashes["zavg_coarse"], h2_set_hash=hashes["h2_submodules"]),
        gls=dict(n_coarse=len(coarse_g), complex_set_hash=hashes["gls_coarse"]),
        m=dict(H3a=dict(M=int(t3a.in_family.sum()), m=m["H3a_lineage_atlas"],
                        tests=int(len(evs["H3a_lineage_atlas"]))),
               H3b=dict(M=int(t3b.in_family.sum()), m=m["H3b_driver_family"], left_family=sorted(h3b_drop),
                        tests=int(len(evs["H3b_driver_family"]))),
               H1=dict(m=len(h1_units), n_events=len(events),
                       n_N_axis=sum(len(e.N) >= 2 for e in events), n_P_axis=sum(len(e.P) >= 1 for e in events)),
               H2="set by the committed H3a and H3b discovery runs (stage-wise)",
               BRIDGE_L_gls=m_bridge,
               controls={r.entry: int(r.family_m) for r in ctab.itertuples(index=False)}),
        h1=dict(n_linked_genes=hb["n_genes"], n_candidates=int(cand_tab.candidate.sum()),
                n_events=len(events), n_event_genes=len(genes_ev),
                by_type={t: int(sum(e.type == t for e in events)) for t in ("LoF", "CN", "silencing")},
                by_link={l: int(sum(e.link == l for e in events)) for l in ("member", "attached")},
                n_locus_events=int(sum(e.locus for e in events)),
                n_lockbox_testable_events=int(sum(e.testable for e in events)),
                n_lof_functional_extra=int(len(hb["extra_check"]))),
        controls=dict(complex_evaluable=sorted(ctab.loc[ctab.complex_evaluable, "design_id"]),
                      not_evaluable={r.design_id: r.reason for r in ctab.itertuples(index=False)
                                     if r.level != "G" and not r.complex_evaluable}),
        lockbox_testable={d: dict(testable=int(g.lockbox_testable.sum()), total=int(len(g)))
                          for d, g in lb.groupby("entry")},
    )
    files = {}
    for p in sorted(OUT.iterdir()):
        if p.is_file() and p.name not in ("FZX_manifest.json", "README.md"):
            files[rel(p)] = sha(p)
    files[rel(CACHE_OUT / "h1_loss_table.parquet")] = sha(CACHE_OUT / "h1_loss_table.parquet")
    for f in ("zavg_map/clusters.tsv", "zavg_map/membership.tsv", "zavg_map/attachments.tsv", "gls_map/clusters.tsv",
              "gls_map/membership.tsv", "tiers/zavg_tiers.tsv", "tiers/gls_tiers.tsv"):
        files[f"results/part2/construction/{f}"] = sha(CON / f)
    head = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    man = dict(
        freeze="FZ-X", design="DESIGN.md v3.1 §0 F4, §2.2, §2.5-2.7", created=time.strftime("%Y-%m-%dT%H:%M:%S"),
        parent_commit=head, fzc_map_commit="27fa9b275a83fd76a9a2f35107fc0b4f3aaa9ac4",
        inputs={rel(FZC_MANIFEST): sha(FZC_MANIFEST),
                "results/part2/contexts/checksums.json": sha(paths.CONTEXT_RESULTS / "checksums.json"),
                "results/part2/lockbox/manifest.json": sha(paths.LOCKBOX_MANIFEST),
                "results/part2/contexts/lof_functional_check.tsv": sha(paths.CONTEXT_RESULTS / "lof_functional_check.tsv"),
                "paralogs_ensembl116": sha(paths.PARALOGS), "coordinates_ensembl116": sha(paths.GENE_COORDS)},
        code={"scripts/39_freeze_fzx.py": sha(ROOT / "scripts/39_freeze_fzx.py"),
              "src/part2/fzx.py": sha(ROOT / "src/part2/fzx.py"), "src/part2/h1.py": sha(ROOT / "src/part2/h1.py"),
              "src/part2/cscore.py": sha(ROOT / "src/part2/cscore.py"),
              "scripts/31_score_blocks.py": sha(ROOT / "scripts/31_score_blocks.py")},
        complex_set_hashes=hashes,
        membership_sha256={"zavg (FZ-C membership.tsv)": fzc["results/part2/construction/zavg_map/membership.tsv"],
                           "gls (FZ-C membership.tsv)": fzc["results/part2/construction/gls_map/membership.tsv"],
                           "zavg engine membership": files[rel(OUT / "zavg_engine_membership.tsv")],
                           "gls engine membership": files[rel(OUT / "gls_engine_membership.tsv")]},
        files=files, registry=reg_report, summary=summary,
        runtime_s=round(time.time() - t0, 1), timings_s={k: round(v, 1) for k, v in timings.items()})
    MANIFEST.write_text(json.dumps(man, indent=1, default=lambda x: x.item() if hasattr(x, "item") else str(x)) + "\n")
    print(json.dumps(dict(summary=summary, runtime_s=man["runtime_s"]), indent=1, default=str))


if __name__ == "__main__":
    main()
