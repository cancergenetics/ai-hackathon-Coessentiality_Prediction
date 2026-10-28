"""Complex scores (DESIGN v3.1 §2.4) and the per-context member removals.

Primary score  S_C(m) = mean over the fixed members of C (after the context's removals) of the gene effect, each value
clipped to [-3, 3] before averaging; NaN members are skipped; S_C(m) is NaN when fewer than 50% of the remaining
members are observed in m. A unit with fewer than 2 remaining members is not evaluable for that context.

Removals for a context (``Removal``): the context gene(s), their paralogs inside C (Ensembl 116, H1 paralog
definition), and every member within +/-10 Mb of a copy-number context gene (Ensembl 116 midpoints; HGNC arm as the
fallback when a gene has no coordinate). Activating target genes of two controls are retained (PC11 NFE2L2, PC19b
PIK3CA). Sex-chromosome contexts remove chrY genes and X-linked members that are paralogs of a chrY gene.

Multiverse scores (reported per hit, never selected):
* ``eigengene``: members clipped, standardised with the fitting rows' mean / SD (A), NaN -> 0, projected on the first
  right singular vector (sign: loadings sum > 0, so a lower score = more dependent, like the mean). Units: SD.
  The fitted (mu, sd, v) are re-used unchanged on other rows (lockbox).
* ``within_line_rank``: mean over members of the within-line percentile rank of the effect among all genes
  (``ranks`` = E_all.rank(axis=1, pct=True); low = more dependent).
* ``expression_masked``: member values set to NaN where log1p TPM < 1 in that line (not expressed), then the
  primary mean (50% rule over all remaining members).
* ``pan_essential_referenced``: primary score divided by -median(clipped effect of a pan-essential reference set) in
  the same line (screen-strength control; NaN when that median is > -0.05).
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from part2.cmap import FrozenMap, symbol

CLIP = 3.0
MIN_FRAC = 0.5
WINDOW_BP = 10_000_000
METHODS = ("S_C_mean_winsor3", "eigengene", "within_line_rank", "expression_masked", "pan_essential_referenced")
LINEAGE_LIKE = ("melanoma", "heme", "heme_vs_breast")


# ------------------------------------------------------------------------------------------------ annotation
@dataclass
class Annotation:
    """Gene annotation keyed by symbol: paralogs (symmetric), chromosome, midpoint, arm (fallback)."""
    paralogs: dict = field(default_factory=dict)       # symbol -> set of symbols
    chrom: dict = field(default_factory=dict)          # symbol -> '1'..'22','X','Y' ('' unknown)
    pos: dict = field(default_factory=dict)            # symbol -> midpoint bp (-1 unknown)
    arm: dict = field(default_factory=dict)            # symbol -> '9p' ('' unknown)

    def is_paralog(self, a: str, b: str) -> bool:
        return b in self.paralogs.get(a, ()) or a in self.paralogs.get(b, ())

    def near(self, a: str, b: str, window: int = WINDOW_BP) -> bool:
        ca, cb = self.chrom.get(a, ""), self.chrom.get(b, "")
        pa, pb = self.pos.get(a, -1), self.pos.get(b, -1)
        if ca and cb and pa >= 0 and pb >= 0:
            return ca == cb and abs(int(pa) - int(pb)) <= window
        aa, ab = self.arm.get(a, ""), self.arm.get(b, "")
        return bool(aa) and aa == ab

    def x_gametolog(self, a: str) -> bool:
        """X-linked gene with an Ensembl paralogue on chrY. The partner's chromosome comes from the annotation, which
        load_annotation fills for every Ensembl gene (not only screened ones; FZ-X review M2: RBMX -> RBMY1*)."""
        return self.chrom.get(a) == "X" and any(self.chrom.get(p) == "Y" for p in self.paralogs.get(a, ()))

    def gametolog(self, a: str) -> bool:
        """X-Y gametolog in either direction: chrX gene with a chrY paralogue, or chrY gene with a chrX paralogue."""
        ca = self.chrom.get(a)
        other = {"X": "Y", "Y": "X"}.get(ca)
        return other is not None and any(self.chrom.get(p) == other for p in self.paralogs.get(a, ()))

    @staticmethod
    def paralogs_from_table(df: pd.DataFrame, min_identity: float = 20.0,
                            types=("within_species_paralog", "other_paralog")) -> dict:
        """Ensembl 116 paralog table (BioMart columns as in datasets/Paralogs) -> symmetric symbol -> set map.
        Rule (H1 paralog_definition): homology type in `types` (gene_split excluded) and min(%id t, %id q) >= 20."""
        c = list(df.columns)
        a, b, idt, idq, typ = c[1], c[3], c[4], c[5], c[6]
        d = df[df[typ].isin(types)]
        d = d[np.minimum(d[idt].astype(float), d[idq].astype(float)) >= min_identity]
        out: dict = {}
        for x, y in zip(d[a].astype(str), d[b].astype(str)):
            if x and y and x != "nan" and y != "nan" and x != y:
                out.setdefault(x, set()).add(y)
                out.setdefault(y, set()).add(x)
        return out


def load_annotation(gene_keys, min_identity: float = 20.0) -> Annotation:
    """Annotation for DepMap gene keys from Ensembl 116 (paralogues, coordinates via Entrez) and HGNC arms."""
    from part2 import paths
    from part2.locus import entrez_positions, parse_location
    keys = list(gene_keys)
    sym = [symbol(k) for k in keys]
    ent = pd.to_numeric(pd.Series([str(k).split("(")[-1].rstrip(")") for k in keys]), errors="coerce")
    pos = entrez_positions().reindex(ent.fillna(-1).astype(np.int64).values)
    coords = pd.read_csv(paths.GENE_COORDS, sep="\t", dtype=str)
    par_tab = pd.read_csv(paths.PARALOGS, sep="\t", dtype=str)
    par, ens_chrom = paralogs_by_ensembl_id(par_tab, coords, dict(zip(sym, ent)), min_identity)
    chrom = dict(ens_chrom)                       # every Ensembl gene, so unscreened paralogs (e.g. RBMY1A1) count
    chrom.update(zip(sym, pos.chrom.fillna("").values))
    h = pd.read_csv(paths.HGNC_FILE, sep="\t", dtype=str, usecols=["entrez_id", "location"]).dropna()
    arm_e = {str(e): parse_location(l)[0] for e, l in zip(h.entrez_id, h.location)}
    return Annotation(paralogs=par, chrom=chrom,
                      pos=dict(zip(sym, pos.pos.fillna(-1).astype(np.int64).values)),
                      arm={s: arm_e.get(str(int(e)) if pd.notna(e) else "", "") for s, e in zip(sym, ent)})


def paralogs_by_ensembl_id(par_tab: pd.DataFrame, coords: pd.DataFrame, entrez_of: dict,
                           min_identity: float = 20.0,
                           types=("within_species_paralog", "other_paralog")) -> tuple[dict, dict]:
    """H1 paralog rule applied on Ensembl gene IDs (FZ-X review note a / M2): Ensembl 116 paralogue pairs (homology
    type in `types`, min(%id t, %id q) >= min_identity) are named by the DepMap symbol of the gene's Entrez ID
    (Ensembl coordinates file: ENSG -> Entrez) when the gene is a DepMap column, else by its Ensembl gene name.
    So renamed genes (DepMap TSLIG3A = HGNC/Ensembl FAM98A) keep their paralogs, and unscreened partners keep their
    chromosome. Returns (symmetric name -> set of names, name -> chromosome for every named Ensembl gene)."""
    c = coords.copy()
    c.columns = ["ensg", "name", "entrez", "chrom", "start", "end"]
    by_entrez = {}
    for s, e in entrez_of.items():
        if pd.notna(e):
            by_entrez.setdefault(str(int(e)), []).append(s)
    names = {}
    for g, n, e in zip(c.ensg, c.name, c.entrez):
        dm = by_entrez.get(str(e)) if isinstance(e, str) and e else None
        names.setdefault(g, set()).update(dm if dm else ([n] if isinstance(n, str) and n else []))
    ens_chrom = {}
    for g, ch in zip(c.ensg, c.chrom):
        for n in names.get(g, ()):
            ens_chrom.setdefault(n, ch)
    col = list(par_tab.columns)
    ga, na, gb, nb, idt, idq, typ = col[0], col[1], col[2], col[3], col[4], col[5], col[6]
    d = par_tab[par_tab[typ].isin(types)]
    d = d[np.minimum(d[idt].astype(float), d[idq].astype(float)) >= min_identity]
    out: dict = {}
    for a, an, b, bn in zip(d[ga], d[na], d[gb], d[nb]):
        A = names.get(a) or ({an} if isinstance(an, str) and an and an != "nan" else set())
        B = names.get(b) or ({bn} if isinstance(bn, str) and bn and bn != "nan" else set())
        for x in A:
            for y in B:
                if x != y:
                    out.setdefault(x, set()).add(y)
                    out.setdefault(y, set()).add(x)
    return out, ens_chrom


# ------------------------------------------------------------------------------------------------ removals
@dataclass(frozen=True)
class Removal:
    genes: tuple = ()          # context gene symbols: removed with their in-complex paralogs
    cn_genes: tuple = ()       # members within +/- window of these genes are removed
    retain: tuple = ()         # never removed (activating target genes of PC11 / PC19b)
    paralogs: bool = True
    chrY: bool = False         # chrY genes and X-linked gametologs removed (LoY / sex contexts)
    ygenes: bool = False       # chrY genes only (PC13: the X gametologs are the readout)
    extra: tuple = ()          # further symbols removed without their paralogs

    def removed(self, members, annot: Annotation, window: int = WINDOW_BP) -> list:
        out = []
        for g in members:
            s = symbol(g)
            if s in self.retain:
                continue
            hit = s in self.genes or s in self.extra
            hit = hit or (self.paralogs and any(annot.is_paralog(s, c) for c in self.genes))
            hit = hit or any(annot.near(s, c, window) for c in self.cn_genes)
            hit = hit or (self.chrY and (annot.chrom.get(s) == "Y" or annot.x_gametolog(s)))
            hit = hit or (self.ygenes and annot.chrom.get(s) == "Y")
            if hit:
                out.append(g)
        return out


def is_lineage_context(cid: str, item: dict | None = None) -> bool:
    if item is not None and "lineage" in str(item.get("type", "")).lower():
        return True
    return cid.startswith("lineage:") or cid in LINEAGE_LIKE


def removal_for_context(cid: str, homdel_share: float | None = None) -> Removal:
    """Machine-readable form of the registry `removed_members` rules for H3a/H3b/H2 contexts (registry text is the
    specification; tests/test_part2_cscore.py checks each rule against it). homdel_share: fraction of the context's
    positive A lines that are hom-del calls of the gene (LoF contexts: +/-10 Mb window iff >= 0.5)."""
    if is_lineage_context(cid) or cid in ("MSI_H", "WGD", "WGD_anyprofile"):
        return Removal()
    if cid in ("LoY", "sex_female"):
        return Removal(chrY=True)
    if cid in ("9p21_loss", "CN_homdel:CDKN2A"):
        return Removal(genes=("MTAP", "CDKN2A", "CDKN2B"), cn_genes=("MTAP",))
    if cid in ("MYC_amp", "MYCN_amp"):
        g = cid.split("_")[0]
        return Removal(genes=(g,), cn_genes=(g,))
    if cid in ("NRF2_active", "LoF:KEAP1"):
        return Removal(genes=("KEAP1", "CUL3", "NFE2L2"))
    if cid in ("TP53_altered", "LoF:TP53", "hotspot:TP53"):        # registry: no +/-10 Mb clause for TP53
        return Removal(genes=("TP53",))
    if cid == "PIK3CA_hotspot_noPTENLoF":
        return Removal(genes=("PIK3CA",))
    if cid == "PTEN_LoF_noPIK3CAhs":
        cid = "LoF:PTEN"
    m = re.match(r"^(hotspot|LoF|LoF_inclMSI|CN_homdel|CN_amp|CN_highamp|CN_loss):([A-Za-z0-9.\-]+)$", cid)
    if m:
        kind, g = m.groups()
        if kind == "hotspot":
            return Removal(genes=(g,))
        if kind.startswith("LoF"):
            return Removal(genes=(g,), cn_genes=(g,) if (homdel_share or 0) >= 0.5 else ())
        return Removal(genes=(g,), cn_genes=(g,))
    if cid.startswith("fusion:"):
        genes = tuple(dict.fromkeys(x for x in re.split(r"::|\|", cid[len("fusion:"):]) if x))
        return Removal(genes=genes)
    raise KeyError(f"no removal rule for context {cid!r}; add one (registry removed_members is the specification)")


# ------------------------------------------------------------------------------------------------ scores
def _clipped(E: pd.DataFrame, members) -> pd.DataFrame:
    return E.reindex(columns=list(members)).astype(float).clip(-CLIP, CLIP)


def _enough(V: pd.DataFrame, min_frac: float = MIN_FRAC) -> np.ndarray:
    k = V.shape[1]
    return (V.notna().sum(axis=1).values >= min_frac * k) & (k > 0)


def mean_score(E: pd.DataFrame, members, min_frac: float = MIN_FRAC) -> pd.Series:
    V = _clipped(E, members)
    s = V.mean(axis=1, skipna=True)
    s[~_enough(V, min_frac)] = np.nan
    return s


def score(E: pd.DataFrame, members, method: str = "S_C_mean_winsor3", *, ranks: pd.DataFrame | None = None,
          expr: pd.DataFrame | None = None, pan_genes=None, fit: dict | None = None,
          min_frac: float = MIN_FRAC) -> tuple[pd.Series, dict]:
    """One unit's score over the rows of E (models x gene keys). Returns (score, fitted parameters)."""
    members = list(members)
    if method in ("S_C_mean_winsor3", "gene_mean"):
        return mean_score(E, members, min_frac), {}
    if method == "eigengene":
        V = _clipped(E, members)
        ok = _enough(V, min_frac)
        if fit is None:
            mu, sd = V.mean(), V.std(ddof=1).replace(0, np.nan)
            Z = ((V - mu) / sd).fillna(0.0)
            _, _, vt = np.linalg.svd(Z.values[ok], full_matrices=False)
            v = vt[0] * (1 if vt[0].sum() >= 0 else -1)
            fit = dict(mu=mu, sd=sd, v=pd.Series(v, index=V.columns))
        Z = ((V - fit["mu"]) / fit["sd"]).fillna(0.0)
        s = pd.Series(Z.values @ fit["v"].reindex(V.columns).fillna(0).values, index=E.index)
        s[~ok] = np.nan
        return s, fit
    if method == "within_line_rank":
        if ranks is None:
            raise ValueError("within_line_rank needs ranks = E_all.rank(axis=1, pct=True)")
        R = ranks.reindex(index=E.index, columns=members)
        s = R.mean(axis=1)
        s[~_enough(R, min_frac)] = np.nan
        return s, {}
    if method == "expression_masked":
        if expr is None:
            raise ValueError("expression_masked needs expr (log1p TPM, models x gene keys)")
        V = _clipped(E, members)
        X = expr.reindex(index=E.index, columns=members)
        V = V.mask(X < 1.0)
        s = V.mean(axis=1)
        s[~_enough(V, min_frac)] = np.nan
        return s, {}
    if method == "pan_essential_referenced":
        if pan_genes is None:
            raise ValueError("pan_essential_referenced needs pan_genes")
        ref = _clipped(E, list(pan_genes)).median(axis=1)
        s = mean_score(E, members, min_frac) / (-ref)
        s[~(ref < -0.05)] = np.nan
        return s, {}
    raise ValueError(f"unknown score method {method!r}")


def gene_unit(E: pd.DataFrame, unit: str) -> pd.Series:
    """Gene-level readout: 'SYMBOL (ENTREZ)' -> raw effect; 'A (x) minus B (y)' -> effect A - effect B."""
    if " minus " in unit:
        a, b = unit.split(" minus ")
        return (E[a.strip()] - E[b.strip()]).rename(unit)
    return E[unit].astype(float).rename(unit)


def unit_scores(E: pd.DataFrame, fmap: FrozenMap, units, removal: Removal, annot: Annotation,
                method: str = "S_C_mean_winsor3", min_members: int = 2, fits: dict | None = None,
                **kw) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Scores of `units` for one context. Returns (models x units scores; NaN column = not evaluable),
    info table (n_members, n_removed, n_used, evaluable, removed), fitted parameters per unit."""
    cols, info, fitted = {}, [], {}
    for u in units:
        mem = fmap.members[u]
        rem = set(removal.removed(mem, annot))
        use = [g for g in mem if g not in rem and g in E.columns]
        ok = len(use) >= min_members
        if ok:
            s, f = score(E, use, method, fit=(fits or {}).get(u), **kw)
            fitted[u] = f
        else:
            s = pd.Series(np.nan, index=E.index)
        cols[u] = s
        info.append(dict(unit=u, n_members=len(mem), n_removed=len(rem), n_used=len(use), evaluable=ok,
                         removed=";".join(sorted(rem))))
    return pd.DataFrame(cols, index=E.index), pd.DataFrame(info), fitted
