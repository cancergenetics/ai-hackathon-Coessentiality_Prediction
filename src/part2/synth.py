"""Synthetic data with known answers for the association engine (DESIGN v3.1 §3 'planted context effects ... H1
classifier, FDR control and HC3/permutation calibration'). No DepMap value is read.

simulate(seed, null=...) -> Synth: models x genes effects with
  * 90 coarse complexes (sizes 6-16, distinct essentiality levels, member loadings on a complex factor), noise genes;
    one parent with two planted sub-modules (CDK6/CCND3 vs CDK4/CCND1) and three more parents with two sub-modules;
  * lineages (8, two primary diseases each), PatientID blocks (~20% of patients with 2 lines), library batch,
    ScreenType3DO / ModelType, MSIScore, burden, Ploidy, WGD, TP53 status, NNMD, CasActivity;
  * genotype contexts correlated with lineage (KRAS enriched in Pancreas/Bowel) and a Pancreas lineage effect on a
    complex that is not a KRAS effect (confounding that lineage_FE must absorb);
  * planted lineage effects (H3a), genotype effects (H3b), a sub-module divergence (H2), H1 loss events:
    paralog compensation (ARID1A -> ARID1B), collapse, residual sensitisation, neutral (incl. a neutral paralog);
  * a PatientID-blocked 25% lockbox split.
tests/part2_synth_pipeline.py builds a throw-away project (freeze manifests, context table, frozen map, registry
entries written with entry_yaml) and runs the whole gated pipeline end to end, exactly as on real data.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from part2 import registry
from part2.cmap import FrozenMap, symbol
from part2.cscore import Annotation

LINEAGES = ["Myeloid", "Lymphoid", "Breast", "Lung", "Skin", "Bowel", "CNS/Brain", "Pancreas"]
H3B_CONTEXTS = ["hotspot:KRAS", "hotspot:BRAF", "LoF:PTEN", "MSI_H", "TP53_altered", "WGD"]
MUT_COV = ["lineage_FE", "library", "ScreenType3DO", "ModelType", "log_burden", "MSIScore"]
BASE_COV = ["lineage_FE", "library", "ScreenType3DO", "ModelType"]
H3B_COV = {"hotspot:KRAS": MUT_COV, "hotspot:BRAF": MUT_COV, "LoF:PTEN": MUT_COV, "TP53_altered": MUT_COV,
           "MSI_H": BASE_COV, "WGD": BASE_COV + ["TP53_altered"]}
H3A_COV = ["library", "ScreenType3DO", "ModelType"]


@dataclass
class Synth:
    E: pd.DataFrame
    meta: pd.DataFrame
    contexts: pd.DataFrame
    fmap: FrozenMap
    annot: Annotation
    fold: pd.Series
    lof: pd.DataFrame
    truth: dict
    named: dict = field(default_factory=dict)       # role -> unit id


def _key(sym: str, i: int) -> str:
    return f"{sym} ({900000 + i})"


def simulate(seed: int = 0, n_models: int = 900, null: bool = False, lineage_null: bool | None = None,
             n_complexes: int = 90, n_noise: int = 400, effect: float = 0.3) -> Synth:
    """null=True: no planted genotype, H1 or H2 effects; lineage_null (default = null): no lineage effects either
    (else the planted lineage effects incl. the confounded Pancreas effect stay)."""
    lineage_null = null if lineage_null is None else lineage_null
    rng = np.random.default_rng(seed)
    n = n_models
    # ---------------------------------------------------------------- models
    n_pat = int(n / 1.2)
    pat = np.concatenate([np.arange(n_pat), rng.choice(n_pat, n - n_pat, replace=False)])
    rng.shuffle(pat)
    lin_of_pat = rng.choice(len(LINEAGES), n_pat, p=np.full(len(LINEAGES), 1 / len(LINEAGES)))
    lin = np.array(LINEAGES)[lin_of_pat[pat]]
    dis = np.array([f"{l} D{int(d)}" for l, d in zip(lin, rng.integers(0, 2, n))])
    ids = [f"ACH-{i:06d}" for i in range(n)]
    lib_ky = (rng.random(n) < 0.3).astype(float)
    st3d = (rng.random(n) < 0.06).astype(float)
    msi = (rng.random(n) < np.where(np.isin(lin, ["Bowel"]), 0.25, 0.04)).astype(float)
    msis = np.where(msi == 1, rng.uniform(25, 60, n), rng.uniform(0, 15, n))
    burden = np.log1p(rng.poisson(np.where(msi == 1, 800, 150)))
    tp53 = (rng.random(n) < 0.65).astype(float)
    tp53[rng.random(n) < 0.02] = np.nan
    wgd = (rng.random(n) < 0.45).astype(float)                # independent of TP53 (clean truth)
    ploidy = np.where(wgd == 1, rng.normal(3.6, 0.3, n), rng.normal(2.1, 0.2, n))
    nowgs = rng.random(n) < 0.1
    wgd[nowgs], ploidy[nowgs] = np.nan, np.nan
    kras = (rng.random(n) < np.select([lin == "Pancreas", lin == "Bowel"], [0.6, 0.35], 0.08)).astype(float)
    braf = (rng.random(n) < np.where(lin == "Skin", 0.45, 0.04)).astype(float)   # independent of KRAS given lineage
    pten = (rng.random(n) < 0.1).astype(float)
    pten[msi == 1] = np.nan                                   # primary LoF: MSI-H lines NaN
    meta = pd.DataFrame(dict(PatientID=[f"P{p}" for p in pat], OncotreeLineage=lin, OncotreePrimaryDisease=dis,
                             LibraryAvana=1 - lib_ky, LibraryKY=lib_ky, ScreenType3DO=st3d,
                             ModelType=np.where(st3d == 1, "Organoid", "Cell line"), MSIScore=msis,
                             log_burden=burden, Ploidy=ploidy, WGD=wgd, TP53_altered=tp53,
                             ScreenNNMD=rng.normal(-3, 0.5, n), CasActivity=np.where(rng.random(n) < 0.1, np.nan,
                                                                                     rng.uniform(40, 90, n))),
                        index=ids)
    ctx = pd.DataFrame({f"lineage:{l}": (lin == l).astype(float) for l in LINEAGES}, index=ids)
    ctx["heme"] = np.isin(lin, ["Myeloid", "Lymphoid"]).astype(float)
    ctx["heme_vs_breast"] = np.where(np.isin(lin, ["Myeloid", "Lymphoid"]), 1.0, np.where(lin == "Breast", 0.0, np.nan))
    ctx["hotspot:KRAS"], ctx["hotspot:BRAF"], ctx["LoF:PTEN"] = kras, braf, pten
    ctx["MSI_H"], ctx["TP53_altered"], ctx["WGD"] = msi, tp53, wgd
    # ---------------------------------------------------------------- genes and complexes
    sizes = rng.integers(6, 17, n_complexes)
    named_members = {  # role -> (symbols placed first in the complex)
        "WRN": ["WRN"], "MYB": ["MYB"], "KRASC": ["SHOC2", "KRAS"], "PTENC": ["PIK3CB"], "TP53C": ["MDM2"],
        "WGDC": ["KIF18A"], "SKIN": [], "LUNG": [], "PANC": [], "CDK": ["CDK6", "CCND3", "CDK4", "CCND1"],
        "EVA": ["ARID1A", "ARID1B"], "EVB": ["GLOSB"], "EVC": ["GLOSC"], "EVD": ["STAG2", "STAG1"],
        "EVE": ["GLOSE"], "SUB2": [], "SUB3": [], "SUB4": [], "LOCUS": [],
    }
    roles = list(named_members)
    mu_c = np.linspace(-1.2, -0.1, n_complexes)[rng.permutation(n_complexes)]
    role_of = {r: i for i, r in enumerate(roles)}
    mu_c[role_of["EVB"]] = -1.0
    mu_c[role_of["EVC"]] = -0.5
    mu_c[role_of["EVA"]], mu_c[role_of["EVD"]], mu_c[role_of["EVE"]] = -0.6, -0.7, -0.8
    sizes[[role_of[r] for r in ("EVA", "EVB", "EVC", "EVD", "EVE")]] = [9, 8, 8, 9, 8]
    sizes[role_of["CDK"]] = 12
    keys, members, gi = [], {}, 0
    named = {}
    for c in range(n_complexes):
        u = f"C{c:04d}"
        syms = list(named_members[roles[c]]) if c < len(roles) else []
        while len(syms) < sizes[c]:
            syms.append(f"G{gi:04d}")
            gi += 1
        mk = []
        for s in syms:
            mk.append(_key(s, len(keys)))
            keys.append(mk[-1])
        members[u] = mk
        if c < len(roles):
            named[roles[c]] = u
    noise = []
    for _ in range(n_noise):
        noise.append(_key(f"N{gi:04d}", len(keys)))
        keys.append(noise[-1])
        gi += 1
    for s in ("MTAP", "MYC"):
        noise.append(_key(s, len(keys)))
        keys.append(noise[-1])
    K = len(keys)
    kidx = {k: i for i, k in enumerate(keys)}
    # ---------------------------------------------------------------- effects
    E = rng.normal(0, 0.2, (n, K))
    mu_g = rng.normal(-0.1, 0.1, K)
    for ci, (u, mk) in enumerate(members.items()):
        f = rng.normal(0, 1, n)
        for k in mk:
            j = kidx[k]
            mu_g[j] = mu_c[ci] + rng.normal(0, 0.08)
            E[:, j] += rng.uniform(0.15, 0.35) * f
    # sub-modules: CDK (CDK6/CCND3 half vs CDK4/CCND1 half) + three neutral parents with two halves
    sub = {}
    for r in ("CDK", "SUB2", "SUB3", "SUB4"):
        u = named[r]
        mk = members[u]
        if r == "CDK":
            s1 = [mk[0], mk[1]] + mk[4:8]
            s2 = [mk[2], mk[3]] + mk[8:12]
        else:
            h = len(mk) // 2
            s1, s2 = mk[:h], mk[h:]
        for s_i, s in enumerate((s1, s2), 1):
            fs = rng.normal(0, 1, n)
            for k in s:
                E[:, kidx[k]] += 0.2 * fs
            sub[f"{u}.S{s_i}"] = (u, s)
    E += mu_g
    E += np.outer(lib_ky, rng.normal(0, 0.08, K))              # library batch
    truth = dict(H3a=set(), H3b=set(), H2=set(), H1={})
    cidx = {k: ci for ci, mk in enumerate(members.values()) for k in mk}
    sig = np.zeros((n, n_complexes))                           # planted lineage/genotype signal on S_C

    def shift(rows, genes, d, track=True):
        for k in genes:
            E[rows, kidx[k]] += d
            if track:
                ci = cidx[k]
                sig[rows, ci] += d / len(members[f"C{ci:04d}"])

    if not lineage_null:
        heme = np.isin(lin, ["Myeloid", "Lymphoid"])
        shift(heme, members[named["MYB"]], -0.45)
        truth["H3a"] |= {("lineage:Myeloid", named["MYB"]), ("lineage:Lymphoid", named["MYB"])}
        shift(lin == "Skin", members[named["SKIN"]], -effect - 0.05)
        truth["H3a"].add(("lineage:Skin", named["SKIN"]))
        shift(lin == "Lung", members[named["LUNG"]], effect)
        truth["H3a"].add(("lineage:Lung", named["LUNG"]))
        shift(lin == "Pancreas", members[named["PANC"]], -0.35)      # confounded with KRAS, not a KRAS effect
        truth["H3a"].add(("lineage:Pancreas", named["PANC"]))
        cdk6_sub = sub[f"{named['CDK']}.S1"][1]
        if not null:
            shift(heme, cdk6_sub, -0.6)                              # H2: CDK6/CCND3 sub-module in heme
            truth["H3a"] |= {("lineage:Myeloid", named["CDK"]), ("lineage:Lymphoid", named["CDK"])}
            truth["H2_planted"] = {("lineage:Myeloid", f"{named['CDK']}.S1"),
                                   ("lineage:Lymphoid", f"{named['CDK']}.S1")}
            # one-vs-rest truth: every lineage whose share of the heme shift differs from the rest by >= 0.03; the
            # sibling sub-module (S2, residualised on a complement containing S1) diverges by construction
            for l in LINEAGES:
                m = lin == l
                if abs(-0.6 * (heme[m].mean() - heme[~m].mean())) >= 0.03:
                    truth["H2"].add((f"lineage:{l}", f"{named['CDK']}.S1"))
                    truth.setdefault("H2_expected_sibling", set()).add((f"lineage:{l}", f"{named['CDK']}.S2"))
    if not null:
        for c, r, d in (("MSI_H", "WRN", -0.5), ("hotspot:KRAS", "KRASC", -effect), ("LoF:PTEN", "PTENC", -effect),
                        ("TP53_altered", "TP53C", effect), ("WGD", "WGDC", -effect)):
            rows = (ctx[c].values == 1)
            shift(rows, [k for k in members[named[r]] if symbol(k) != "KRAS"], d)
            truth["H3b"].add((c, named[r]))
    # ---------------------------------------------------------------- H1 loss events (LoF, ~25% of lines)
    lof_syms = ["ARID1A", "GLOSB", "GLOSC", "STAG2", "GLOSE"]
    lof = pd.DataFrame({s: (rng.random(n) < (0.3 if s == "STAG2" else 0.25)).astype(float) for s in lof_syms},
                       index=ids)
    if not null:
        L = lof.values.astype(bool)
        a1b = [k for k in members[named["EVA"]] if symbol(k) == "ARID1B"]
        shift(L[:, 0], a1b, -0.8, False)                                  # paralog compensation
        for k in members[named["EVB"]]:                             # collapse toward 0
            if symbol(k) != "GLOSB":
                E[L[:, 1], kidx[k]] += 0.6 * abs(mu_g[kidx[k]])
        shift(L[:, 2], [k for k in members[named["EVC"]] if symbol(k) != "GLOSC"], -0.4, False)   # sensitisation
    for s, r in zip(lof_syms, ("EVA", "EVB", "EVC", "EVD", "EVE")):  # lost gene's own dependency vanishes
        k = next(k for k in members[named[r]] if symbol(k) == s)
        E[lof[s].values == 1, kidx[k]] = rng.normal(0, 0.15, int((lof[s] == 1).sum()))
    eva, evb, evc, evd, eve = (f"{named[r]}|{s}" for r, s in zip(("EVA", "EVB", "EVC", "EVD", "EVE"), lof_syms))
    if null:
        truth["H1"] = {eva: dict(N="NEUTRAL", P="NEUTRAL"), evb: dict(N="NEUTRAL"), evc: dict(N="NEUTRAL"),
                       evd: dict(N="NEUTRAL", P="NEUTRAL"), eve: dict(N="NEUTRAL")}
    else:
        truth["H1"] = {eva: dict(N="NEUTRAL", P="PARALOG COMPENSATION"), evb: dict(N="COLLAPSE"),
                       evc: dict(N="RESIDUAL SENSITISATION"), evd: dict(N="NEUTRAL", P="NEUTRAL"),
                       eve: dict(N="NEUTRAL")}
    E = E.astype(np.float32).astype(float)
    E[rng.random(E.shape) < 0.003] = np.nan
    Edf = pd.DataFrame(E, index=ids, columns=keys)
    # one-vs-rest truth for H3a from the planted signal: |mean(lineage) - mean(rest)| >= 0.03 on S_C; complexes with
    # no planted lineage/genotype signal at all are the null units (a discovery there is a false positive)
    units = list(members)
    truth["H3a_planted"] = set(truth["H3a"])
    for l in LINEAGES:
        m = lin == l
        diff = sig[m].mean(0) - sig[~m].mean(0)
        truth["H3a"] |= {(f"lineage:{l}", units[ci]) for ci in np.flatnonzero(np.abs(diff) >= 0.03)}
    truth["null_units"] = {units[ci] for ci in np.flatnonzero(np.abs(sig).sum(0) == 0)}
    # ---------------------------------------------------------------- map, annotation, fold
    level = {u: 1 for u in members}
    parent = {u: "" for u in members}
    allm = dict(members)
    for s, (p, mk) in sub.items():
        allm[s], level[s], parent[s] = list(mk), 2, p
    fmap = FrozenMap(members=allm, level=level, parent=parent, locus={named["LOCUS"]: True})
    par = {"ARID1A": {"ARID1B"}, "ARID1B": {"ARID1A"}, "STAG2": {"STAG1"}, "STAG1": {"STAG2"},
           "KRAS": {"NRAS"}, "NRAS": {"KRAS"}}
    chrom, pos = {}, {}
    for k in keys:
        s = symbol(k)
        chrom[s] = str(rng.integers(1, 23))
        pos[s] = int(rng.integers(1, 200_000_000))
    annot = Annotation(paralogs=par, chrom=chrom, pos=pos)
    pats = pd.Series(meta.PatientID.unique())
    Lp = set(pats.sample(frac=0.25, random_state=int(seed) + 3000))
    fold = pd.Series(np.where(meta.PatientID.isin(Lp), "L", "A"), index=ids)
    return Synth(E=Edf, meta=meta, contexts=ctx, fmap=fmap, annot=annot, fold=fold, lof=lof, truth=truth, named=named)


# ------------------------------------------------------------------------------------------------ registry text
def _q(s) -> str:
    return '"' + str(s).replace("\\", "\\\\").replace('"', '\\"') + '"'


def entry_yaml(aid: str, direction: str, contexts: list, units, score: str, covariates: list, fdr: dict,
               frozen_inputs: list, sign: str = "two-sided", extra: str = "", unit_hash: str | None = None) -> str:
    lines = [f"id: {aid}", f"direction: {direction}", "hypothesis: |", f"  Synthetic entry {aid}.",
             f"sign: {_q(sign)}", "context_family:", "  name: synthetic", "  contexts:"]
    for c in contexts:
        lines.append(f"    - id: {_q(c['id'])}")
        lines.append(f"      definition: {_q(c.get('definition', 'synthetic context ' + c['id']))}")
        if "covariates" in c:
            lines.append("      covariates: [" + ", ".join(c["covariates"]) + "]")
        if "type" in c:
            lines.append(f"      type: {_q(c['type'])}")
    lines.append(f"complex_set_hash: {unit_hash or registry.set_hash(units)}")
    lines.append(f"score: {score}")
    lines.append("covariates: [" + ", ".join(covariates) + "]")
    lines.append("fdr_family:")
    for k, v in fdr.items():
        lines.append(f"  {k}: [" + ", ".join(v) + "]" if isinstance(v, list) else f"  {k}: {v}")
    lines += ["success_criterion: \"synthetic\"", "stopping_rule: \"synthetic\"", "frozen_inputs:"]
    lines += [f"  - {p}" for p in frozen_inputs]
    lines.append("status: open")
    return "\n".join(lines) + "\n" + extra


def ctrl_yaml(aid, design_id, ctx_id, sign, gene_units, complex_unit, score, direction, covs, finputs, ctx_type=""):
    units = list(gene_units) + ([complex_unit] if complex_unit else [])
    extra = ["control:", f"  design_id: {design_id}", "  level: \"G+C\"", "  gene_units:"]
    extra += [f"    - {_q(g)}" for g in gene_units]
    extra += [f"  complex_unit: {_q(complex_unit or 'not evaluable at complex level')}",
              "kill_rule:", f"  direction: {direction}", "  counted: true", "  level: C"]
    c = dict(id=ctx_id)
    if ctx_type:
        c["type"] = ctx_type
    return entry_yaml(aid, "controls", [c], units, score, covs, dict(method="none", q=0.05), finputs, sign=sign,
                      extra="\n".join(extra) + "\n")

