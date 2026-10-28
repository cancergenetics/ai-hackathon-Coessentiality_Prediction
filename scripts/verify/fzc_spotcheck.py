#!/usr/bin/env python
"""Reviewer spot-check of FZ-C clusters: (a) 5 random level-1 T1 zavg complexes (gene lists, BioPlex edges, locus
flag); (b) 5 random locus-flagged zavg clusters with the flag recomputed independently from Ensembl 116 coordinates
(by Entrez ID; gene midpoint; primary chromosomes) and HGNC `location` arm; (c) the locus flag recomputed for ALL
clusters of both maps and compared with the committed flag.

The rule recomputed (DESIGN §2.3): locus iff >= 50% of members within one 10-Mb window on one chromosome, or >= 50%
on one chromosome arm. Genes without coordinates / arm count only in the denominator. Multi-chromosome Entrez IDs
or genes spanning > 10 Mb have no position (as documented in src/part2/locus.py).

  .venv/bin/python scripts/verify/fzc_spotcheck.py --out <dir> [--seed 6100]
"""
import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CON = ROOT / "results" / "part2" / "construction"
W = 10_000_000
PRIMARY = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}


def coords():
    c = pd.read_csv(ROOT / "datasets" / "Paralogs" / "ensembl116_gene_coordinates.tsv", sep="\t", dtype=str)
    c.columns = ["ensg", "name", "entrez", "chrom", "start", "end"]
    c = c.dropna(subset=["entrez"])
    c = c[c.chrom.isin(PRIMARY)]
    c[["start", "end"]] = c[["start", "end"]].astype(np.int64)
    g = c.groupby(c.entrez.astype(np.int64)).agg(n=("chrom", "nunique"), chrom=("chrom", "first"),
                                                 lo=("start", "min"), hi=("end", "max"))
    ok = (g.n == 1) & (g.hi - g.lo <= W)
    g["mid"] = np.where(ok, (g.lo + g.hi) // 2, -1)
    g.loc[~ok, "chrom"] = ""
    return g


def arms():
    h = pd.read_csv(ROOT / "datasets" / "HGNC" / "hgnc_complete_set.txt", sep="\t", dtype=str,
                    usecols=["symbol", "entrez_id", "location"])
    h = h.dropna(subset=["entrez_id"])
    h["arm"] = h.location.fillna("").map(lambda s: (lambda m: m.group(1) + m.group(2) if m else "")(
        re.match(r"^(\d+|X|Y)([pq])", s)))
    return h.set_index(h.entrez_id.astype(np.int64)).arm


def flag(entrez, C, AR):
    n = len(entrez)
    cc = C.reindex(entrez)
    best_w, best_chr = 0, ""
    for ch, grp in cc[(cc.chrom.fillna("") != "") & (cc.mid >= 0)].groupby("chrom"):
        x = np.sort(grp.mid.values)
        j, k = 0, 0
        for i in range(len(x)):
            while x[i] - x[j] > W:
                j += 1
            if i - j + 1 > best_w:
                best_w, best_chr = i - j + 1, ch
    a = AR.reindex(entrez).fillna("")
    a = a[a != ""]
    vc = a.value_counts()
    arm_n, arm = (int(vc.iloc[0]), vc.index[0]) if len(vc) else (0, "")
    return dict(locus=bool(best_w / n >= 0.5 or arm_n / n >= 0.5), window_share=round(best_w / n, 3),
                window_chr=best_chr, arm_share=round(arm_n / n, 3), top_arm=arm)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=6100)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    C, AR = coords(), arms()
    rep = {}
    for m, d in (("zavg", "zavg_map"), ("gls", "gls_map")):
        tab = pd.read_csv(CON / d / "clusters.tsv", sep="\t", keep_default_na=False)
        rec = []
        for r in tab.itertuples():
            f = flag([int(x) for x in str(r.entrez).split(";")], C, AR)
            rec.append(dict(cluster_id=r.cluster_id, committed=str(r.locus) == "True", **f,
                            committed_window=r.window_share, committed_arm=r.arm_share))
        rec = pd.DataFrame(rec)
        rec.to_csv(out / f"{m}_locus_recomputed.tsv", sep="\t", index=False)
        dis = rec[rec.committed != rec.locus]
        rep[m] = dict(n=len(rec), agree=int((rec.committed == rec.locus).sum()), discordant=dis.cluster_id.tolist(),
                      window_share_maxabsdiff=float((rec.window_share - rec.committed_window).abs().max()),
                      arm_share_maxabsdiff=float((rec.arm_share - rec.committed_arm).abs().max()))
    rng = np.random.default_rng(a.seed)
    tz = pd.read_csv(CON / "tiers" / "zavg_tiers.tsv", sep="\t")
    ct = pd.read_csv(CON / "zavg_map" / "clusters.tsv", sep="\t", keep_default_na=False).set_index("cluster_id")
    rz = pd.read_csv(out / "zavg_locus_recomputed.tsv", sep="\t").set_index("cluster_id")
    t1 = tz[(tz.level == 1) & (tz.tier == "T1")].cluster_id.values
    lf = ct[(ct.locus.astype(str) == "True")].index.values
    pick = dict(T1=sorted(rng.choice(t1, 5, replace=False).tolist()), locus=sorted(rng.choice(lf, 5, replace=False).tolist()))
    rows = []
    for kind, ids in pick.items():
        for c in ids:
            t = tz.set_index("cluster_id").loc[c]
            rows.append(dict(kind=kind, cluster_id=c, size=int(ct.loc[c, "size"]), tier=t.tier, edges=t.edges,
                             n_obs=t.n_obs, q_bioplex=t.q_bioplex, locus_committed=ct.loc[c, "locus"],
                             locus_recomputed=bool(rz.loc[c, "locus"]), top_arm=rz.loc[c, "top_arm"],
                             arm_share=rz.loc[c, "arm_share"], window_chr=rz.loc[c, "window_chr"],
                             window_share=rz.loc[c, "window_share"], members=ct.loc[c, "members"]))
    pd.DataFrame(rows).to_csv(out / "spotcheck.tsv", sep="\t", index=False)
    rep["picked"] = pick
    (out / "spotcheck.json").write_text(json.dumps(rep, indent=1, default=str))
    print(json.dumps(rep, indent=1, default=str))
    print(pd.DataFrame(rows).to_string())


if __name__ == "__main__":
    main()
