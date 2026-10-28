#!/usr/bin/env python
"""Independent FZ-X review (results/part2/verify/fzx_review.md).

Re-derives from raw inputs, WITHOUT the producer's fzx / h1 / cscore / covariates / contexts code:
  * unit sets and set hashes (zavg / GLS coarse, H2 sub-modules, H1 / H1S components);
  * H3a / H3b / BRIDGE-L context counts (n_A, ref_A, n_L, ref_L) from the context table, the A list and L = U \\ A;
  * evaluability and m per family: complete-case over context + registered covariates + score observed (>= 50% of
    the used members observed in A's gene-effect MISSINGNESS mask; values discarded at once), removals re-implemented
    from the registry removed_members text with Ensembl-116 paralogs keyed by Ensembl gene ID (Entrez -> ENSG);
  * control units (target clustered at level 1 / in a sub-module);
  * the H1 enumeration rule 1-8 from the omics label matrices (LoF, hom-del, damaging, log1p TPM), the LoF functional
    check tables, the FZ-C attachment table and Ensembl 116 (paralogs by ENSG, midpoints for +/-10 Mb), plus the
    effect of the producer's two resolutions (LoF check on genes not in the FZ-LB table; level-1 non-locus
    attachments only).
Only the registry YAML parser (part2.registry.load_yaml) is imported from the producer's package.
No dependency value, association statistic, trials/ledger entry or lockbox outcome is read or written.

  PART1_DIR=... .venv/bin/python scripts/verify/fzx_independent.py   (writes results/part2/verify/fzx/)
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from part2.registry import load_yaml  # noqa: E402  (parser only)

DS = ROOT / "datasets"
CC = ROOT / "data_cache" / "part2" / "contexts"
CON = ROOT / "results" / "part2" / "construction"
FZX = ROOT / "results" / "part2" / "fzx"
REG = ROOT / "docs" / "registry"
OUT = ROOT / "results" / "part2" / "verify" / "fzx"
WIN = 10_000_000


def shash(units) -> str:
    return hashlib.sha256("\n".join(sorted(set(map(str, units)))).encode()).hexdigest()[:16]


def ent_of_col(c):
    m = re.match(r"^(.*) \((\d+)\)$", str(c))
    return int(m.group(2)) if m else None


def by_entrez(df):
    e = [ent_of_col(c) for c in df.columns]
    keep = [i for i, x in enumerate(e) if x is not None]
    out = df.iloc[:, keep].copy()
    out.columns = [e[i] for i in keep]
    return out.loc[:, ~pd.Index(out.columns).duplicated()]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rep = {}
    # ------------------------------------------------------------------ labels, A, L
    ctx = pd.read_parquet(CC / "context_table.parquet")
    U = list(ctx.index)
    A = (ROOT / "results/part2/lockbox/A_models.txt").read_text().split()
    A = [a for a in A if a.strip()]
    assert set(A) <= set(U)
    isA = pd.Index(U).isin(A)
    isL = ~isA
    rep["n_U_A_L"] = [len(U), int(isA.sum()), int(isL.sum())]

    # ------------------------------------------------------------------ maps
    def mp(name):
        cl = pd.read_csv(CON / f"{name}_map/clusters.tsv", sep="\t")
        mem = pd.read_csv(CON / f"{name}_map/membership.tsv", sep="\t")
        return cl, mem
    clz, memz = mp("zavg")
    clg, memg = mp("gls")
    att = pd.read_csv(CON / "zavg_map/attachments.tsv", sep="\t")
    coarse_z = sorted(clz.cluster_id[clz.level == 1])
    coarse_g = sorted(clg.cluster_id[clg.level == 1])
    sub = clz[clz.level == 2]
    nsub = sub.groupby("parent").size()
    h2 = sorted(sub.cluster_id[sub.parent.isin(nsub.index[nsub >= 2])])
    rep["units"] = dict(zavg=[len(coarse_z), shash(coarse_z)], gls=[len(coarse_g), shash(coarse_g)],
                        h2=[len(h2), shash(h2), int((nsub >= 2).sum())], n_submodules=int(len(sub)))
    members_z = memz.groupby("cluster_id").entrez.apply(lambda s: sorted(map(int, s))).to_dict()
    members_g = memg.groupby("cluster_id").entrez.apply(lambda s: sorted(map(int, s))).to_dict()
    sym = dict(zip(memz.entrez.astype(int), memz.symbol)) | dict(zip(memg.entrez.astype(int), memg.symbol)) | \
        dict(zip(att.entrez.astype(int), att.symbol))
    locus_z = dict(zip(clz.cluster_id, clz.locus.astype(str).eq("True")))

    # ------------------------------------------------------------------ Ensembl 116 by ENSG
    co = pd.read_csv(DS / "Paralogs/ensembl116_gene_coordinates.tsv", sep="\t", dtype=str)
    co.columns = ["ensg", "name", "entrez", "chrom", "start", "end"]
    co = co.dropna(subset=["entrez"])
    e2g = co.groupby(co.entrez.astype(int)).ensg.apply(set).to_dict()
    prim = [str(i) for i in range(1, 23)] + ["X", "Y"]
    cp = co[co.chrom.isin(prim)].copy()
    cp["s"], cp["e"] = cp.start.astype(np.int64), cp.end.astype(np.int64)
    gp = cp.groupby(cp.entrez.astype(int)).agg(n=("chrom", "nunique"), chrom=("chrom", "first"), lo=("s", "min"),
                                                 hi=("e", "max"))
    gp = gp[(gp.n == 1) & (gp.hi - gp.lo <= WIN)]
    mid = {int(k): (r.chrom, (r.lo + r.hi) // 2) for k, r in gp.iterrows()}
    hg = pd.read_csv(DS / "HGNC/hgnc_complete_set.txt", sep="\t", dtype=str, usecols=["symbol", "entrez_id", "location"])
    arm_of = {}
    for r in hg.dropna(subset=["entrez_id", "location"]).itertuples():
        m = re.match(r"^(\d+|X|Y)([pq])", r.location)
        if m:
            arm_of[int(r.entrez_id)] = m.group(1) + m.group(2)
    par = pd.read_csv(DS / "Paralogs/ensembl116_human_paralogues.tsv", sep="\t", dtype=str)
    pc = list(par.columns)
    par = par[par[pc[6]].isin(["within_species_paralog", "other_paralog"])]
    par = par[np.minimum(par[pc[4]].astype(float), par[pc[5]].astype(float)) >= 20]
    pset = set(zip(par[pc[0]], par[pc[2]])) | set(zip(par[pc[2]], par[pc[0]]))

    def is_par(a: int, b: int) -> bool:
        return any((x, y) in pset for x in e2g.get(a, ()) for y in e2g.get(b, ()) if x != y)

    def near(a: int, b: int) -> bool:
        if a in mid and b in mid:
            return mid[a][0] == mid[b][0] and abs(mid[a][1] - mid[b][1]) <= WIN
        return bool(arm_of.get(a)) and arm_of.get(a) == arm_of.get(b)

    hgnc_sym2e = {r.symbol: int(r.entrez_id) for r in hg.dropna(subset=["entrez_id"]).itertuples()}
    chrom_of = {k: v[0] for k, v in mid.items()}

    # ------------------------------------------------------------------ registry entries (parser only)
    entries = {p.stem: load_yaml(p.read_text()) for p in sorted(REG.glob("*.yaml"))}

    # ------------------------------------------------------------------ context counts
    rows = []
    for k in ("H3a_lineage_atlas", "H3b_driver_family", "BRIDGE-L_confirmation_rate", "H2_submodule_divergence"):
        for c in entries[k]["context_family"]["contexts"]:
            v = ctx[c["id"]].values
            got = dict(n_A=int(((v == 1) & isA).sum()), ref_A=int(((v == 0) & isA).sum()),
                       n_L=int(((v == 1) & isL).sum()), ref_L=int(((v == 0) & isL).sum()))
            reg = {x: c.get(x) for x in got}
            t = got["n_L"] >= 5 and got["ref_L"] >= 5
            rows.append(dict(entry=k, context=c["id"], role=c.get("role", "primary"), **got,
                             registry_match=got == reg, testable=t, registry_testable=c.get("lockbox_testable"),
                             testable_match=t == c.get("lockbox_testable")))
    cc = pd.DataFrame(rows)
    cc.to_csv(OUT / "context_counts_recomputed.tsv", sep="\t", index=False)
    rep["context_counts_all_match"] = bool(cc.registry_match.all() and cc.testable_match.all())
    lin = ctx["OncotreeLineage"] if "OncotreeLineage" in ctx.columns else None
    models = pd.read_csv(DS / "DepMAP/Model.csv", index_col=0)
    linA = models.OncotreeLineage.reindex(A).value_counts()
    rep["H3a_lineages_ge12_A"] = sorted(linA.index[linA >= 12])
    h3a_ids = [c["id"] for c in entries["H3a_lineage_atlas"]["context_family"]["contexts"]]
    rep["H3a_registry_ids"] = h3a_ids
    rep["H3a_lineage_list_matches"] = sorted(f"lineage:{x}" for x in linA.index[linA >= 12]) + ["melanoma"] == \
        sorted(h3a_ids)
    h3b = entries["H3b_driver_family"]["context_family"]
    rep["H3b_primary"] = [(c["id"], c["n_A"], c["n_L"], c["lockbox_testable"]) for c in h3b["contexts"]
                          if c.get("role", "primary") == "primary"]
    rep["H3b_excluded_fzr"] = [(x["id"], x["n_A"]) for x in h3b.get("excluded", [])]
    exc_check = []
    for x in h3b.get("excluded", []):
        if x["id"] in ctx.columns:
            v = ctx[x["id"]].values
            exc_check.append((x["id"], int(((v == 1) & isA).sum()), x["n_A"]))
        else:
            exc_check.append((x["id"], None, x["n_A"]))
    rep["H3b_excluded_recount"] = exc_check

    # ------------------------------------------------------------------ missingness mask of A's gene effects
    ge = pd.read_csv(DS / "DepMAP/CRISPRGeneEffect.csv", index_col=0)
    ge = ge.reindex([a for a in A if a in ge.index])
    obs = by_entrez(ge.notna())
    del ge
    obs = obs.reindex(A).fillna(False).astype(bool)
    rep["A_with_gene_effect"] = int(obs.any(axis=1).sum())

    # ------------------------------------------------------------------ covariates (source columns only)
    conf = pd.read_csv(DS / "DepMAP/CRISPRConfounders.csv", index_col=0).reindex(A)
    mA = models.reindex(A)
    ctxA = ctx.reindex(A)
    libcols = [c for c in conf.columns if c.startswith("Library")]

    def cov_ok(names, cid):
        drop = set()
        lineage_ctx = cid.startswith("lineage:") or cid == "melanoma"
        if lineage_ctx:
            drop.add("lineage_FE")
        if cid == "MSI_H":
            drop |= {"MSIScore", "log_burden"}
        if cid.startswith("WGD"):
            drop |= {"Ploidy", "WGD"}
        if "TP53" in cid:
            drop.add("TP53_altered")
        ok = np.ones(len(A), bool)
        for n in names:
            if n in drop:
                continue
            if n == "lineage_FE":
                s = mA.OncotreeLineage.notna().values
            elif n == "library":
                s = conf[libcols].notna().all(axis=1).values
            elif n in ("ScreenType3DO", "CasActivity"):
                s = conf[n].notna().values
            elif n == "NNMD":
                s = conf["ScreenNNMD"].notna().values
            elif n == "ModelType":
                s = mA.ModelType.notna().values
            else:
                s = ctxA[n].notna().values
            ok &= s
        return ok

    # ------------------------------------------------------------------ removals from registry text
    def homdel_share(g):
        e = hgnc_sym2e.get(g)
        hd = pd.read_parquet(CC / "cn_homdel.parquet", columns=[str(e)]) if e is not None else None
        pos = ctxA[f"LoF:{g}"].values == 1
        if hd is None or str(e) not in hd.columns or not pos.any():
            return np.nan
        return float((hd.reindex(A)[str(e)].values[pos] == 1).mean())

    shares = {}

    def removal(cid):
        if cid.startswith("lineage:") or cid in ("melanoma", "MSI_H", "WGD"):
            return dict(genes=[], cn=[], chrY=False)
        if cid in ("LoY", "sex_female"):
            return dict(genes=[], cn=[], chrY=True)
        if cid == "9p21_loss":
            return dict(genes=["MTAP", "CDKN2A", "CDKN2B"], cn=["MTAP"], chrY=False)
        if cid in ("MYC_amp", "MYCN_amp"):
            g = cid.split("_")[0]
            return dict(genes=[g], cn=[g], chrY=False)
        if cid == "NRF2_active":
            return dict(genes=["KEAP1", "CUL3", "NFE2L2"], cn=[], chrY=False)
        if cid == "TP53_altered":
            return dict(genes=["TP53"], cn=[], chrY=False)
        if cid.startswith("hotspot:"):
            return dict(genes=[cid.split(":")[1]], cn=[], chrY=False)
        if cid.startswith("LoF:"):
            g = cid.split(":")[1]
            if g not in shares:
                shares[g] = homdel_share(g)
            return dict(genes=[g], cn=[g] if (shares[g] or 0) >= 0.5 else [], chrY=False)
        if cid.startswith("fusion:"):
            return dict(genes=[x for x in re.split(r"::|\|", cid[7:]) if x], cn=[], chrY=False)
        raise KeyError(cid)

    def removed(mem, rm):
        ge_ = [hgnc_sym2e[g] for g in rm["genes"] if g in hgnc_sym2e]
        cn_ = [hgnc_sym2e[g] for g in rm["cn"] if g in hgnc_sym2e]
        out = set()
        for x in mem:
            if x in ge_ or any(is_par(x, g) for g in ge_) or any(near(x, g) for g in cn_):
                out.add(x)
            if rm["chrY"] and (chrom_of.get(x) == "Y" or x in xgam):
                out.add(x)
        return out

    ysY = [y for y in chrom_of if chrom_of[y] == "Y"]
    ensY = set().union(*[e2g.get(y, set()) for y in ysY])
    xgam = {x for x in chrom_of if chrom_of[x] == "X" and any((a, b) in pset for a in e2g.get(x, ()) for b in ensY)}

    obs_cols = set(obs.columns)
    O = obs.values
    col_ix = {c: i for i, c in enumerate(obs.columns)}

    def evaluate(entry_id, units, members):
        e = entries[entry_id]
        out = []
        for c in e["context_family"]["contexts"]:
            if c.get("role", "primary") != "primary":
                continue
            cid = c["id"]
            cv = ctxA[cid].values
            base = (~np.isnan(cv.astype(float))) & cov_ok(c.get("covariates", e["covariates"]), cid)
            rm = removal(cid)
            for u in units:
                mem = members[u]
                rmv = removed(mem, rm)
                use = [x for x in mem if x not in rmv and x in obs_cols]
                if len(use) < 2:
                    out.append(dict(entry=entry_id, context=cid, unit=u, n_used=len(use), n=0, n_context=0,
                                    n_ref=0, evaluable=False))
                    continue
                sc = O[:, [col_ix[x] for x in use]].sum(1) >= 0.5 * len(use)
                r = base & sc
                nc = int((r & (cv == 1)).sum())
                nr = int(r.sum()) - nc
                out.append(dict(entry=entry_id, context=cid, unit=u, n_used=len(use), n=int(r.sum()), n_context=nc,
                                n_ref=nr, evaluable=nc >= 12 and nr >= 12))
        return pd.DataFrame(out)

    ev = pd.concat([evaluate("H3a_lineage_atlas", coarse_z, members_z),
                    evaluate("H3b_driver_family", coarse_z, members_z),
                    evaluate("BRIDGE-L_confirmation_rate", coarse_g, members_g)], ignore_index=True)
    ev.to_csv(OUT / "evaluability_recomputed.tsv.gz", sep="\t", index=False,
              compression={"method": "gzip", "mtime": 0})
    pe = pd.read_csv(FZX / "evaluability_tests.tsv.gz", sep="\t")
    j = ev.merge(pe, on=["entry", "context", "unit"], suffixes=("", "_prod"), how="outer", indicator=True)
    rep["evaluability_rows"] = dict(mine=len(ev), producer=len(pe), both=int((j._merge == "both").sum()))
    diffs = {}
    for col in ("n_used", "n", "n_context", "n_ref", "evaluable"):
        d = j[j[col] != j[f"{col}_prod"]]
        diffs[col] = int(len(d))
        if len(d):
            d.head(50).to_csv(OUT / f"evaluability_diff_{col}.tsv", sep="\t", index=False)
    rep["evaluability_diffs"] = diffs
    h3a_set = set(h3a_ids)
    m = dict(H3a=int(ev.evaluable[ev.entry == "H3a_lineage_atlas"].sum()),
             H3b=int(ev.evaluable[ev.entry == "H3b_driver_family"].sum()),
             BRIDGE_H3a=int(ev.evaluable[(ev.entry == "BRIDGE-L_confirmation_rate") & ev.context.isin(h3a_set)].sum()),
             BRIDGE_H3b=int(ev.evaluable[(ev.entry == "BRIDGE-L_confirmation_rate") & ~ev.context.isin(h3a_set)].sum()))
    rep["m"] = m
    s3b = ev[ev.entry == "H3b_driver_family"].groupby("context").n_context.max().sort_values()
    rep["H3b_min_ncontext_cc"] = s3b.head(6).to_dict()
    rep["H3b_contexts_below12_cc"] = sorted(s3b.index[s3b < 12])
    rep["bridge_nonevaluable"] = ev[(ev.entry == "BRIDGE-L_confirmation_rate") & ~ev.evaluable][
        ["context", "unit", "n_used"]].values.tolist()
    rep["lof_homdel_shares"] = {k: (None if not np.isfinite(v) else round(v, 3)) for k, v in shares.items()}

    # ------------------------------------------------------------------ controls: target clustered?
    def units_with(g, level):
        e = hgnc_sym2e.get(g)
        cl = clz.set_index("cluster_id")
        return sorted(u for u, mem in members_z.items() if e in mem and int(cl.level[u]) == level)
    rep["control_targets"] = {g: dict(level1=units_with(g, 1), level2=units_with(g, 2))
                              for g in ("SMARCA2", "STAG1", "KIF18A", "CDK6", "CCND3", "CDK4", "CCND1", "BRD9",
                                        "BICRA", "ARID1A", "ARID1B", "DPF2", "VPS4A", "WRN", "MDM2", "SHOC2", "PRMT5",
                                        "PELO", "NFE2L2", "DDX3X", "PIK3CB", "PIK3CA", "SOX10", "MYB")}

    # ------------------------------------------------------------------ H1 enumeration (rules 1-8), independent
    lof_all = pd.read_parquet(CC / "lof_matrix.parquet")
    hd_all = pd.read_parquet(CC / "cn_homdel.parquet")
    raw = CC / "raw"
    dam = by_entrez(pd.read_parquet(next(raw.glob("OmicsSomaticMutationsMatrixDamaging.*.parquet"))).drop(
        columns="SequencingID")).reindex(U)
    expr = by_entrez(pd.read_parquet(next(raw.glob("OmicsExpressionTPMLogp1HumanProteinCodingGenes.*.parquet"))).drop(
        columns="SequencingID")).reindex(U)
    func = pd.read_csv(ROOT / "results/part2/contexts/lof_functional_check.tsv", sep="\t")
    extra = pd.read_csv(FZX / "h1_lof_check_extra.tsv", sep="\t")
    pass_orig = set(func.entrez[func["pass"].astype(bool)].astype(int))
    in_orig = set(func.entrez.astype(int))
    pass_extra = set(extra.entrez[extra["pass"].astype(bool)].astype(int))
    msi = ctx["MSI_H"].reindex(U).values
    mmr = {hgnc_sym2e[g] for g in ("MLH1", "MSH2", "MSH6", "PMS2", "MLH3", "MSH3", "PMS1", "EPCAM")}
    lineage = models.OncotreeLineage.reindex(U)

    # linked genes and links
    lvl = dict(zip(clz.cluster_id, clz.level))
    links_member = {}
    for u in coarse_z:
        for x in members_z[u]:
            links_member.setdefault(x, []).append(u)
    attv = att.copy()
    attv["lvl"] = [lvl[c] for c in attv.cluster]
    attv["locus_attach"] = attv.locus_attach.astype(str).eq("True")
    rep["attachments"] = dict(total=len(attv), level1=int((attv.lvl == 1).sum()),
                              level1_nonlocus=int(((attv.lvl == 1) & ~attv.locus_attach).sum()),
                              level1_locus=int(((attv.lvl == 1) & attv.locus_attach).sum()),
                              level2=int((attv.lvl == 2).sum()))

    def enumerate_(att_rows, pass_set, label):
        links = {x: [(u, "member") for u in us] for x, us in links_member.items()}
        for r in att_rows.itertuples():
            links.setdefault(int(r.entrez), [])
            if (r.cluster, "member") not in links[int(r.entrez)]:
                links[int(r.entrez)].append((r.cluster, "attached"))
        genes = sorted(links)
        def col(df, g):
            return df[str(g)].reindex(U).values.astype(float) if str(g) in df.columns else np.full(len(U), np.nan)
        evs, loss_tab = [], {}
        E = expr.reindex(columns=[g for g in genes if g in expr.columns])
        med = pd.DataFrame(np.tile(E.median().values, (len(U), 1)), index=U, columns=E.columns)
        for l, idx in lineage.groupby(lineage).groups.items():
            sub_ = E.loc[idx]
            lm = sub_.median().where(sub_.notna().sum() >= 5)
            ok = lm.notna()
            med.loc[idx, ok[ok].index] = np.tile(lm[ok].values, (len(idx), 1))
        for g in genes:
            lof = col(lof_all, g)
            if g not in mmr:
                lof = np.where(msi == 0, lof, np.nan)
            hd = col(hd_all, g)
            dm = dam[g].values.astype(float) if g in dam.columns else np.full(len(U), np.nan)
            ex = E[g].values.astype(float) if g in E.columns else np.full(len(U), np.nan)
            if g in E.columns:
                mg = med[g].values
            else:
                mg = np.full(len(U), np.nan)
            a = (lof == 1) & (g in pass_set)
            b = hd == 1
            c = (ex < 1) & (mg >= 3) & ~b
            one = a | b | c
            zero = (dm == 0) & (hd == 0) & ~(ex < 1)
            loss = np.where(one, 1.0, np.where(zero, 0.0, np.nan))
            nA, rA = int(((loss == 1) & isA).sum()), int(((loss == 0) & isA).sum())
            if nA < 12 or rA < 12:
                continue
            la = (loss == 1) & isA
            sh_lof, sh_hd = (a & la).sum() / la.sum(), (b & la).sum() / la.sum()
            sh_lof_excl = (a & ~b & la).sum() / la.sum()
            typ = "CN" if sh_hd >= 0.5 else ("LoF" if sh_lof_excl >= 0.5 else "silencing")
            typ_incl = "CN" if (sh_hd >= 0.5 and not sh_lof > 0.5) else ("LoF" if sh_lof >= 0.5 else
                                                                        ("CN" if sh_hd >= 0.5 else "silencing"))
            nL, rL = int(((loss == 1) & isL).sum()), int(((loss == 0) & isL).sum())
            for u, how in links[g]:
                rest = [x for x in members_z[u] if x != g]
                P = [x for x in rest if is_par(x, g)]
                N = [x for x in rest if x not in P and not (typ == "CN" and near(x, g))]
                if len(N) < 2 and len(P) < 1:
                    continue
                evs.append(dict(event=f"{u}|{sym.get(g, g)}", entrez=g, unit=u, link=how, type=typ,
                                type_inclusive_LoF=typ_incl, n_P=len(P), n_N=len(N), n_A=nA, ref_A=rA, n_L=nL,
                                ref_L=rL, testable=nL >= 5 and rL >= 5, locus=locus_z[u],
                                n_A_lof=int((a & la).sum()), n_A_homdel=int((b & la).sum()),
                                n_A_sil=int((c & la).sum()),
                                genetic_share=round(float(((a | b) & la).sum() / la.sum()), 6), lof_check=("orig" if g in in_orig else
                                                                        "extra" if g in set(extra.entrez) else "none"),
                                lof_pass=g in pass_set))
        return pd.DataFrame(evs)

    att_l1 = attv[(attv.lvl == 1) & ~attv.locus_attach]
    mine = enumerate_(att_l1, pass_orig | pass_extra, "registered")
    mine.to_csv(OUT / "h1_events_recomputed.tsv", sep="\t", index=False)
    prod = pd.read_csv(FZX / "h1_events.tsv", sep="\t", keep_default_na=False)
    prod["key"] = prod.unit + "|" + prod.gene
    mine["key"] = mine.unit + "|" + [p.split("|")[1] for p in mine.event]
    # compare on (unit, entrez) to avoid symbol-version issues
    gk = pd.read_csv(FZX / "zavg_engine_membership.tsv", sep="\t")
    dep2e = {k.split(" (")[0]: int(k.split("(")[-1].rstrip(")")) for k in gk.gene}
    prod["entrez"] = [dep2e.get(g, -1) for g in prod.gene]
    jj = mine.merge(prod, on=["unit", "entrez"], how="outer", suffixes=("", "_prod"), indicator=True)
    rep["h1"] = dict(mine=len(mine), producer=len(prod), both=int((jj._merge == "both").sum()),
                     only_mine=jj[jj._merge == "left_only"].event.tolist(),
                     only_producer=jj[jj._merge == "right_only"].event_prod.tolist(),
                     by_type_mine=mine.type.value_counts().to_dict(),
                     type_inclusive_differs=mine[mine.type != mine.type_inclusive_LoF].event.tolist())
    b = jj[jj._merge == "both"]
    rep["h1_field_diffs"] = {c: b[b[c] != b[f"{c}_prod"]].event.tolist()[:20]
                             for c in ("type", "n_P", "n_N", "n_A", "ref_A", "n_L", "ref_L", "testable", "locus", "link")}
    # hashes / m for H1 and H1S from my events
    def units_of(df):
        return sorted([f"{r.unit}|{r.g}|N" for r in df.assign(g=[e.split('|')[1] for e in df.event]).itertuples()
                       if r.n_N >= 2] +
                      [f"{r.unit}|{r.g}|P" for r in df.assign(g=[e.split('|')[1] for e in df.event]).itertuples()
                       if r.n_P >= 1])
    p_ = mine[mine.type.isin(["LoF", "CN"])]
    s_ = mine[mine.type == "silencing"]
    rep["h1_m_hash"] = dict(all=[len(units_of(mine)), shash(units_of(mine))],
                            H1=[len(p_), len(units_of(p_)), shash(units_of(p_))],
                            H1S=[len(s_), len(units_of(s_)), shash(units_of(s_))],
                            H1_testable=int(p_.testable.sum()), H1S_testable=int(s_.testable.sum()),
                            locus=int(mine.locus.sum()), locus_silencing=int((mine.locus & (mine.type == "silencing")).sum()),
                            P_axis=int((mine.n_P >= 1).sum()), H1_P_axis=int((p_.n_P >= 1).sum()))
    # sensitivity of the producer's resolutions
    alt1 = enumerate_(att_l1, pass_orig, "no extra LoF check")
    alt2 = enumerate_(attv[attv.lvl == 1], pass_orig | pass_extra, "incl locus attachments")
    def summ(df):
        return dict(n=len(df), by_type=df.type.value_counts().to_dict())
    rep["sensitivity"] = dict(
        extra_lof_check_ignored=dict(**summ(alt1),
                                     events_changed=sorted(set(map(tuple, alt1[["event", "type", "n_A"]].values.tolist()))
                                                           ^ set(map(tuple, mine[["event", "type", "n_A"]].values.tolist())))[:40]),
        locus_attachments_included=dict(**summ(alt2),
                                        extra_events=sorted(set(alt2.event) - set(mine.event))[:40]))
    # 10 spot checks
    pick = (list(mine[mine.type.isin(["LoF", "CN"])].sort_values("event").event[::5][:5]) +
            list(mine[mine.type == "silencing"].sample(5, random_state=11).event))
    rep["spot_checks"] = mine[mine.event.isin(pick)][
        ["event", "link", "type", "n_A", "ref_A", "n_L", "ref_L", "n_A_lof", "n_A_homdel", "n_A_sil", "lof_check",
         "lof_pass", "n_P", "n_N", "locus"]].to_dict("records")
    # the 8 DepMap/HGNC symbol mismatches
    g8 = pd.read_csv(Path(__import__("os").environ.get("PART1_DIR", "")) / "outputs/all_pairs/genes.tsv", sep="\t")
    g8 = g8[g8.depmap_column.notna()]
    g8["dsym"] = g8.depmap_column.str.split(" \\(").str[0]
    mm = g8[(g8.dsym != g8.symbol) & g8.entrez.astype(int).isin(set(memz.entrez) | set(memg.entrez) | set(att.entrez))]
    ens_names = set(par[pc[1]]) | set(par[pc[3]])
    rep["symbol_mismatches"] = [dict(hgnc=r.symbol, depmap=r.dsym, entrez=int(r.entrez),
                                     depmap_in_ensembl_paralogs=r.dsym in ens_names,
                                     hgnc_in_ensembl_paralogs=r.symbol in ens_names) for r in mm.itertuples()]
    # ---- delta review: amendment A2 (event genetic iff >= 50% of A loss models carry a LoF or hom-del call)
    g_ = mine[mine.genetic_share >= 0.5]
    s2_ = mine[mine.genetic_share < 0.5]
    a2p = pd.read_csv(FZX / "h1_event_split_A2.tsv", sep="\t", keep_default_na=False) \
        if (FZX / "h1_event_split_A2.tsv").exists() else None
    a2 = dict(H1=[len(g_), len(units_of(g_)), shash(units_of(g_)), int(g_.testable.sum())],
              H1S=[len(s2_), len(units_of(s2_)), shash(units_of(s2_)), int(s2_.testable.sum())],
              moved_vs_A1=sorted(set(g_.event) ^ set(p_.event)),
              near_threshold=mine[(mine.genetic_share >= 0.3) & (mine.genetic_share < 0.7)][
                  ["event", "type", "genetic_share"]].values.tolist())
    if a2p is not None:
        pm = dict(zip(a2p.event, a2p.entry))
        mm_ = {e: ("H1_within_complex_GI" if gs >= 0.5 else "H1S_silencing_GI")
               for e, gs in zip(mine.event, mine.genetic_share)}
        a2["entry_mismatches_vs_producer"] = sorted(e for e in mm_ if pm.get(e) != mm_[e])
        a2["share_mismatches_vs_producer"] = sorted(
            e for e, gs in zip(mine.event, mine.genetic_share)
            if abs(float(dict(zip(a2p.event, a2p.genetic_share)).get(e, -1)) - gs) > 1e-3)
        for k, f in (("H1_units_file", "units_h1_components_A2.txt"), ("H1S_units_file", "units_h1s_components_A2.txt")):
            a2[k] = shash((FZX / f).read_text().split())
    rep["A2"] = a2
    # gametologs (my derivation: X-linked gene with an Ensembl-116 paralogue, by ENSG, located on chrY)
    rep["xgam_in_maps"] = sorted(sym.get(x, x) for x in xgam
                                 if any(x in members_z[u] for u in members_z) or any(x in members_g[u] for u in members_g))
    pd.Series(sorted(xgam)).to_csv(OUT / "x_gametologs_entrez.tsv", sep="\t", index=False, header=["entrez"])
    (OUT / "fzx_independent.json").write_text(json.dumps(rep, indent=1, default=str) + "\n")
    print(json.dumps(rep, indent=1, default=str)[:20000])


if __name__ == "__main__":
    main()
