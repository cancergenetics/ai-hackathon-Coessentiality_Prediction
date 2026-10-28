"""Locus-artefact flag (DESIGN v3.1 §2.3): a cluster is a "locus cluster" if >= 50% of its members lie within
10 Mb (Ensembl 116 gene coordinates) or on one chromosome arm (HGNC cytogenetic location).

* Coordinates: datasets/Paralogs/ensembl116_gene_coordinates.tsv joined on NCBI (Entrez) gene ID; only primary
  assembly chromosomes (1-22, X, Y, MT) are used; position = gene midpoint. An Entrez ID with primary loci on more
  than one chromosome, or spanning > 10 Mb, has no position (counts in the denominator only).
* "Within 10 Mb": the largest number of members on one chromosome whose midpoints fit in a window of `window_bp`.
* Arm: '19q13.43' -> '19q' (HGNC `location`); band ('19q13') is reported for interpretation.
Genes without a parsable location / position count in the denominator but never in a shared arm or window.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

_ARM = re.compile(r"^(\d+|X|Y)([pq])")
_BAND = re.compile(r"^((?:\d+|X|Y)[pq]\d+)")
PRIMARY = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
WINDOW_BP = 10_000_000


def parse_location(loc) -> tuple[str, str]:
    """'19q13.43' -> ('19q', '19q13'); unparsable / missing -> ('', '')."""
    if not isinstance(loc, str):
        return "", ""
    a, b = _ARM.match(loc), _BAND.match(loc)
    return (a.group(1) + a.group(2) if a else ""), (b.group(1) if b else "")


def entrez_positions(path=None) -> pd.DataFrame:
    """Entrez -> (chrom, pos = midpoint bp) from the Ensembl coordinates table (primary chromosomes only)."""
    from part2 import paths
    c = pd.read_csv(path or paths.GENE_COORDS, sep="\t", dtype=str)
    c.columns = ["ensg", "name", "entrez", "chrom", "start", "end"]
    c = c.dropna(subset=["entrez"])
    c = c[c.chrom.isin(PRIMARY)].copy()
    c["start"], c["end"] = c.start.astype(np.int64), c.end.astype(np.int64)
    g = c.groupby("entrez").agg(nchrom=("chrom", "nunique"), chrom=("chrom", "first"), lo=("start", "min"),
                                hi=("end", "max"))
    ok = (g.nchrom == 1) & (g.hi - g.lo <= WINDOW_BP)
    g["pos"] = np.where(ok, (g.lo + g.hi) // 2, -1)
    g.loc[~ok, "chrom"] = ""
    g.index = g.index.astype(np.int64)
    return g[["chrom", "pos"]]


def gene_locations(genes: pd.DataFrame | None = None, hgnc_path=None, coords_path=None) -> pd.DataFrame:
    """Per genes.tsv row: location, arm, band (HGNC via hgnc_id), chrom, pos (Ensembl via Entrez; '' / -1 if none)."""
    from part2 import data, paths
    genes = data.protein_genes() if genes is None else genes
    h = pd.read_csv(hgnc_path or paths.HGNC_FILE, sep="\t", dtype=str, usecols=["hgnc_id", "location"])
    loc = genes[["hgnc_id"]].merge(h, on="hgnc_id", how="left").location.values
    ab = [parse_location(l) for l in loc]
    pos = entrez_positions(coords_path).reindex(genes.entrez.astype(np.int64).values)
    return pd.DataFrame(dict(location=loc, arm=[a for a, _ in ab], band=[b for _, b in ab],
                             chrom=pos.chrom.fillna("").values, pos=pos.pos.fillna(-1).astype(np.int64).values),
                        index=genes.index)


def _max_share(labels):
    labels = [l for l in labels if l]
    if not labels:
        return "", 0
    s = pd.Series(labels).value_counts()
    return s.index[0], int(s.iloc[0])


def max_in_window(chrom: np.ndarray, pos: np.ndarray, window_bp: int = WINDOW_BP) -> int:
    """Largest number of genes on one chromosome whose positions fit in a window of window_bp."""
    best = 0
    ok = (np.asarray(pos) >= 0) & (np.asarray(chrom) != "")
    chrom, pos = np.asarray(chrom)[ok], np.asarray(pos)[ok]
    for ch in np.unique(chrom):
        p = np.sort(pos[chrom == ch])
        j = np.searchsorted(p, p + window_bp, side="right")
        best = max(best, int((j - np.arange(len(p))).max()))
    return best


def locus_flag(members, arm: np.ndarray, band: np.ndarray, frac: float = 0.5, chrom=None, pos=None,
               window_bp: int = WINDOW_BP) -> dict:
    """Arm rule (HGNC) OR window rule (Ensembl coordinates, when chrom/pos are given)."""
    members = np.asarray(members)
    n = len(members)
    a_lab, a_n = _max_share(arm[members])
    b_lab, b_n = _max_share(band[members])
    w_n = max_in_window(np.asarray(chrom)[members], np.asarray(pos)[members], window_bp) if chrom is not None else 0
    flag = a_n / n >= frac or w_n / n >= frac
    return dict(top_arm=a_lab, arm_share=a_n / n, top_band=b_lab, band_share=b_n / n, window_share=w_n / n,
                locus=bool(flag))
