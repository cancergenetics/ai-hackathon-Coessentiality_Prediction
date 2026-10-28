#!/usr/bin/env python
"""Descriptive summary of the frozen FZ-C zavg map (after tiers): counts, size distribution, gene coverage by
A-essentiality class, locus-flag rate, tiers, attachments, and the largest T1 non-locus complexes by member symbols.

POST HOC column only: HGNC `gene_group` (most frequent group among members and its share), clearly labelled as a
post hoc annotation; it is never used for construction, selection, tiers or the bridge. No curated complex resource
is read. Writes results/part2/construction/summary/.

  PYTHONHASHSEED=0 .venv/bin/python scripts/38_describe_maps.py
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from part2 import data, genefilter, paths, runsets  # noqa: E402

CON = paths.RESULTS / "construction"


def posthoc_group(members: str, sym2grp: pd.Series) -> tuple[str, float]:
    ms = members.split(";")
    g = [x for m in ms for x in str(sym2grp.get(m, "")).split("|") if x and x != "nan"]
    if not g:
        return "", 0.0
    s = pd.Series(g).value_counts()
    return s.index[0], round(int(s.iloc[0]) / len(ms), 3)


def describe(name: str, d: Path, tiers: pd.DataFrame, sym2grp, S_ess: pd.Series) -> dict:
    tab = pd.read_csv(d / "clusters.tsv", sep="\t", keep_default_na=False)
    memb = pd.read_csv(d / "membership.tsv", sep="\t")
    tab["locus"] = tab.locus.astype(str) == "True"
    l1, l2 = tab[tab.level == 1], tab[tab.level == 2]
    m1 = memb[memb.level == 1]
    cov = {c: dict(in_S=int((S_ess == c).sum()), clustered=int((m1.ess_class == c).sum())) for c in ("pan", "selective", "other")}
    for c in cov:
        cov[c]["frac"] = round(cov[c]["clustered"] / max(cov[c]["in_S"], 1), 3)
    out = dict(n_level1=int(len(l1)), n_level2=int(len(l2)), n_parents_with_submodules=int(l2.parent.nunique()),
               genes_level1=int(l1["size"].sum()), genes_S=int(len(S_ess)),
               coverage_by_class=cov,
               size_level1=dict(zip(["min", "q25", "median", "q75", "max"], np.percentile(l1["size"], [0, 25, 50, 75, 100]).tolist())),
               size_level2=dict(zip(["min", "q25", "median", "q75", "max"], np.percentile(l2["size"], [0, 25, 50, 75, 100]).tolist())) if len(l2) else {},
               size_bins_level1={b: int(((l1["size"] >= lo) & (l1["size"] <= hi)).sum()) for b, (lo, hi) in
                                 {"3-4": (3, 4), "5-8": (5, 8), "9-16": (9, 16), "17-60": (17, 60)}.items()},
               locus_rate_level1=round(float(l1.locus.mean()), 3), locus_rate_level2=round(float(l2.locus.mean()), 3) if len(l2) else None)
    if tiers is not None:
        t = tiers.merge(tab[["cluster_id", "parent"]], on="cluster_id")
        out["tiers"] = {f"level{lev}": t[t.level == lev].tier.value_counts().to_dict() for lev in (1, 2)}
        out["tiers_nonlocus_level1"] = t[(t.level == 1) & ~t.locus.astype(bool)].tier.value_counts().to_dict()
        top = t[(t.level == 1) & (t.tier == "T1") & ~t.locus.astype(bool)].sort_values(["size", "q_bioplex"], ascending=[False, True])
        rows = []
        for _, r in top.head(25).iterrows():
            g, share = posthoc_group(r.members, sym2grp)
            subs = tab[tab.parent == r.cluster_id]
            rows.append(dict(cluster_id=r.cluster_id, size=int(r["size"]), bioplex_edges=int(r.edges),
                             bioplex_obs_members=int(r.n_obs), null_mean_edges=round(float(r.null_mean_edges), 2),
                             q_bioplex=float(r.q_bioplex), n_submodules=int(len(subs)), members=r.members,
                             POSTHOC_hgnc_gene_group=g, POSTHOC_group_share=share))
        pd.DataFrame(rows).to_csv(CON / "summary" / f"{name}_top_T1_complexes.tsv", sep="\t", index=False)
        out["top_T1_nonlocus"] = [dict(cluster_id=r["cluster_id"], size=r["size"], members=r["members"],
                                       POSTHOC_hgnc_gene_group=r["POSTHOC_hgnc_gene_group"]) for r in rows[:10]]
    att = d / "attachments.tsv"
    if att.exists():
        a = pd.read_csv(att, sep="\t")
        out["attachments"] = dict(pairs=int(len(a)), genes=int(a.gene_index.nunique()),
                                  pairs_not_locus=int((~a.locus_attach).sum()),
                                  genes_not_locus=int(a.loc[~a.locus_attach, "gene_index"].nunique()))
    return out


def main():
    (CON / "summary").mkdir(parents=True, exist_ok=True)
    S = runsets.S_mask()
    _, st = genefilter.dependency_bearing_for_models(list(runsets.A()), return_stats=True)
    fb = st["frac_below"][S]
    S_ess = pd.Series(np.where(fb >= 0.9, "pan", np.where(fb >= 0.05, "selective", "other")))
    h = pd.read_csv(paths.HGNC_FILE, sep="\t", dtype=str, usecols=["symbol", "gene_group"])   # POST HOC only
    sym2grp = h.set_index("symbol").gene_group
    res = {}
    for name in ("zavg", "gls"):
        tp = CON / "tiers" / f"{name}_tiers.tsv"
        tiers = pd.read_csv(tp, sep="\t", keep_default_na=False) if tp.exists() else None
        if tiers is not None:
            for c in ("edges", "n_obs", "q_bioplex", "null_mean_edges"):
                tiers[c] = pd.to_numeric(tiers[c], errors="coerce")
            tiers["locus"] = tiers.locus.astype(str) == "True"
        res[name] = describe(name, CON / f"{name}_map", tiers, sym2grp, S_ess)
    res["note"] = ("POSTHOC_hgnc_gene_group: HGNC gene_group of the most frequent group among members (post hoc "
                   "annotation; not used for construction, selection, tiers or the bridge)")
    (CON / "summary" / "summary.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps({k: {kk: v for kk, v in res[k].items() if kk != "top_T1_nonlocus"} for k in ("zavg", "gls")}, indent=1, default=float))
    for r in res["zavg"].get("top_T1_nonlocus", []):
        print(r["cluster_id"], r["size"], r["POSTHOC_hgnc_gene_group"], r["members"][:200])


if __name__ == "__main__":
    main()
