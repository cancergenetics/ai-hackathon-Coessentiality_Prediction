"""Context table (DESIGN v3.1 §2.2; FZ-X inputs). Built from DepMap Omics + Model.csv + Ensembl 116 coordinates + HGNC
locations only. No gene-effect (dependency) value is read anywhere in this module: the model universe U comes from
the ModelID column of the CRISPR gene-effect file only (``data.universe_model_ids``).

Conventions
-----------
* Index: ModelID over U (1,345 models). Each Omics matrix is reduced to the row flagged ``IsDefaultEntryForModel ==
  'Yes'`` (OmicsProfiles default profile of that data type); a model without such a row is NaN for every context
  that needs that data type. Uniqueness of the default row per model and data type is asserted.
* Binary contexts are float (1 / 0 / NaN). Logic is three-valued (Kleene): an OR is 1 if any term is 1, 0 if all
  terms are 0, NaN otherwise; an AND is 0 if any term is 0. No imputation: complete-case per test downstream.
* Genes are matched across matrices by Entrez ID (column key 'SYMBOL (ENTREZ)').

Definitions (one place; mirrored in the summary TSV)
---------------------------------------------------
* MSI-H: MSIScore > 20 (OmicsGlobalSignatures, default DNA profile). log_burden = log1p(row sum of the damaging matrix).
* CN: r = OmicsCNGeneWGS (ploidy-relative linear ratio), clipped at 10; a = r x Ploidy (signature of the same WGS
  profile). hom-del: a < 0.5 (r < 0.25 when Ploidy is NaN) AND expression log1p(TPM) < 1 where measured;
  partial loss: 0.5 <= a <= Ploidy - 0.75; loss (any): a <= Ploidy - 0.75 or hom-del; amp: r >= 2; high-amp: r >= 3.
* LoF(g): damaging == 2 OR (damaging >= 1 AND a < Ploidy - 0.75) OR hom-del OR (male line AND g on chrX outside the
  PARs AND damaging >= 1). Primary LoF contexts set MSI-H models to NaN, except for MMR genes (MLH1, MSH2, MSH6,
  PMS2, MLH3, MSH3, PMS1, EPCAM); ``LoF_inclMSI:`` columns keep them.
* Hotspot(g): hotspot matrix >= 1 (separate from LoF).
* TP53 status: 1 = altered (hotspot or LoF), 0 = WT (no damaging, no hotspot, no hom-del), NaN = ambiguous (single
  damaging hit without loss / hom-del / hotspot) or unknown.
* Fusions: OmicsFusionFiltered (default RNA profile); Gene1/Gene2 'SYMBOL (ENSG...)' parsed to symbols; unordered
  pairs from the pre-registered list; self-fusions excluded; models without a default RNA profile are NaN.
* Sex from Model.csv (Unknown -> NaN). LoY: male line with summed chrY (male-specific region) log1p TPM below a
  threshold set at the antimode of the bimodal male distribution (Gaussian KDE; see ``loy_threshold``).
* NRF2-active: LoF(KEAP1) OR LoF(CUL3) OR hotspot(NFE2L2) (MSI-H lines: 1 only through the hotspot, else NaN).
* Lineage: OncotreeLineage one-hot; melanoma = OncotreePrimaryDisease contains 'Melanoma' (any site).
"""
from __future__ import annotations

import json
import logging
import re
import zipfile
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from part2 import data, paths
from part2.locus import entrez_positions, parse_location

log = logging.getLogger("part2.contexts")

META = ["SequencingID", "ModelID", "ModelConditionID", "IsDefaultEntryForModel", "IsDefaultEntryForMC"]
MMR_GENES = ("MLH1", "MSH2", "MSH6", "PMS2", "MLH3", "MSH3", "PMS1", "EPCAM")
# Named LoF contexts: DESIGN §2.1 lockbox list + §2.6 H3b + §2.7 controls (+ MMR, KEAP1/CUL3 for NRF2, TP53, NF1).
LOF_GENES = ("ARID1A", "SMARCA4", "STAG2", "SMARCB1", "PTEN", "RB1", "KMT2D", "CREBBP", "EP300", "KEAP1", "CUL3",
             "TP53", "NF1", "STK11", "CDKN2A", "MTAP", "FOCAD", "VPS4B") + MMR_GENES
HOTSPOT_GENES = ("KRAS", "NRAS", "HRAS", "BRAF", "PIK3CA", "NFE2L2", "TP53")
HOMDEL_GENES = ("MTAP", "CDKN2A", "FOCAD", "VPS4B", "PTEN", "RB1", "SMARCB1")
LOSS_GENES = ("VPS4B",)
AMP_GENES = ("MYC", "MYCN")
# Pre-registered fusion list (DESIGN §2.2): context name -> list of unordered symbol pairs.
FUSIONS = {
    "EWSR1::FLI1": [("EWSR1", "FLI1")],
    "EWSR1::ERG": [("EWSR1", "ERG")],
    "EWSR1::FLI1|ERG": [("EWSR1", "FLI1"), ("EWSR1", "ERG")],
    "SS18::SSX": [("SS18", "SSX1"), ("SS18", "SSX2"), ("SS18", "SSX4")],
    "KMT2A-r": [("KMT2A", p) for p in ("MLLT3", "AFF1", "MLLT1", "MLLT10", "AFDN", "ELL")],
    "PAX3|7::FOXO1": [("PAX3", "FOXO1"), ("PAX7", "FOXO1")],
    "BCR::ABL1": [("BCR", "ABL1")],
}
# Expected minimum counts in U (docs/reviews/round2_biology.md item 6); the build fails if a count is below these.
FUSION_MIN_U = {"EWSR1::FLI1|ERG": 19, "SS18::SSX": 5, "KMT2A-r": 13, "PAX3|7::FOXO1": 7, "BCR::ABL1": 12}
PAR = ((0, 2_781_479), (155_701_383, 156_030_895))     # GRCh38 X PAR1 / PAR2


# ----------------------------------------------------------------------------------------------- three-valued logic
def k_or(*xs) -> np.ndarray:
    xs = [np.asarray(x, dtype=np.float64) for x in xs]
    any1 = np.zeros(np.broadcast(*xs).shape, bool)
    all0 = np.ones_like(any1)
    for x in xs:
        any1 |= (x == 1)
        all0 &= (x == 0)
    return np.where(any1, 1.0, np.where(all0, 0.0, np.nan))


def k_and(*xs) -> np.ndarray:
    xs = [np.asarray(x, dtype=np.float64) for x in xs]
    any0 = np.zeros(np.broadcast(*xs).shape, bool)
    all1 = np.ones_like(any0)
    for x in xs:
        any0 |= (x == 0)
        all1 &= (x == 1)
    return np.where(any0, 0.0, np.where(all1, 1.0, np.nan))


def cmp(x, op, thr) -> np.ndarray:
    """Comparison returning 1/0 with NaN propagated."""
    x = np.asarray(x, dtype=np.float64)
    with np.errstate(invalid="ignore"):
        r = {"<": x < thr, "<=": x <= thr, ">": x > thr, ">=": x >= thr, "==": x == thr}[op]
    return np.where(np.isnan(x) | np.isnan(np.asarray(thr, dtype=np.float64)), np.nan, r.astype(np.float64))


# ------------------------------------------------------------------------------------------------------- loaders
def _entrez_of(col: str):
    m = re.match(r"^(.*) \((\d+)\)$", col)
    return (m.group(1), int(m.group(2))) if m else (col, None)


def read_omics(path, universe: list[str], cache: bool = True) -> tuple[pd.DataFrame, pd.Series]:
    """Default-profile rows of an Omics matrix for the universe: (values models x 'SYMBOL (ENTREZ)', SequencingID).
    Cached as parquet keyed by the CSV's SHA-256."""
    h = data.sha256_file(path)
    pq = paths.CONTEXT_CACHE / "raw" / f"{path.stem}.{h[:16]}.parquet"
    if cache and pq.exists():
        df = pd.read_parquet(pq)
    else:
        df = pd.read_csv(path, index_col=0, low_memory=False)
        df = df[df.IsDefaultEntryForModel == "Yes"]
        if df.ModelID.duplicated().any():
            raise ValueError(f"{path.name}: >1 default row for {df.ModelID[df.ModelID.duplicated()].tolist()[:5]}")
        df = df.set_index("ModelID")
        vals = [c for c in df.columns if c not in META]
        df = pd.concat([df[["SequencingID"]], df[vals].astype(np.float32)], axis=1)
        if cache:
            pq.parent.mkdir(parents=True, exist_ok=True)
            df.to_parquet(pq)
    df = df.reindex([m for m in universe if m in df.index])
    return df.drop(columns="SequencingID"), df["SequencingID"]


def read_signatures(universe: list[str]) -> pd.DataFrame:
    """OmicsGlobalSignatures: the default DNA (wes/wgs) profile row per model (RNA rows carry only PAM50)."""
    prof = pd.read_csv(paths.OMICS_PROFILES).set_index("SequencingID")
    s = pd.read_csv(paths.OMICS_SIGNATURES, index_col=0)
    s["DataType"] = prof.DataType.reindex(s.SequencingID).values
    s = s[(s.IsDefaultEntryForModel == "Yes") & s.DataType.isin(["wes", "wgs"])]
    if s.ModelID.duplicated().any():
        raise ValueError("OmicsGlobalSignatures: >1 default DNA row per model")
    return s.set_index("ModelID").reindex(universe)


def read_signatures_any(universe: list[str]) -> pd.DataFrame:
    """Per model, the first DNA-profile signature row with non-NaN Ploidy: default first, then WGS before WES, then
    latest SequencingDate (sensitivity columns only)."""
    prof = pd.read_csv(paths.OMICS_PROFILES).set_index("SequencingID")
    s = pd.read_csv(paths.OMICS_SIGNATURES, index_col=0)
    s["DataType"] = prof.DataType.reindex(s.SequencingID).values
    s["date"] = prof.SequencingDate.reindex(s.SequencingID).values
    s = s[s.DataType.isin(["wes", "wgs"]) & s.Ploidy.notna()].copy()
    s["_o1"] = (s.IsDefaultEntryForModel != "Yes").astype(int)
    s["_o2"] = (s.DataType != "wgs").astype(int)
    s = s.sort_values(["ModelID", "_o1", "_o2", "date"], ascending=[True, True, True, False])
    return s.drop_duplicates("ModelID").set_index("ModelID").reindex(universe)


def read_fusions(universe: list[str]) -> pd.DataFrame:
    f = pd.read_csv(paths.OMICS_FUSIONS, index_col=0)
    f = f[f.IsDefaultEntryForModel == "Yes"].copy()
    for c in ("Gene1", "Gene2"):
        # 'SYMBOL (ENSG...)'; also 'SYMBOL (.)' (no Ensembl gene) and viral 'NC_...@... (NC_...)' entries
        parsed = f[c].astype(str).str.extract(r"^\s*(\S+?)\s*\(([^()]*)\)\s*$")
        if parsed[0].isna().any():
            raise ValueError(f"fusion {c}: {int(parsed[0].isna().sum())} entries do not parse as 'SYMBOL (ID)'")
        ens = parsed[1].str.startswith("ENSG")
        if ens.mean() < 0.5:                                  # Gene2: ~79% ENSG, rest "(.)"
            raise ValueError(f"fusion {c}: only {ens.mean():.1%} of entries carry an Ensembl gene ID")
        f[c + "_symbol"] = parsed[0]
    return f[f.ModelID.isin(universe)]


def read_proteomics(cache: bool = True) -> pd.DataFrame:
    """CCLE proteomics (Nusinow 2020) normalized protein quant, models (ModelID via Model.csv CCLEName) x gene
    symbol (Gene_Symbol as published; the row with most measured lines per symbol). Replicate TenPx runs of the same
    line are averaged. Parsed from the xlsx with the standard library (no openpyxl in the pinned env)."""
    path = paths.CCLE_PROTEOMICS_XLSX
    h = data.sha256_file(path)
    pq = paths.CONTEXT_CACHE / "raw" / f"ccle_protein.{h[:16]}.parquet"
    if cache and pq.exists():
        return pd.read_parquet(pq)
    rows = _xlsx_rows(path, "Normalized Protein Expression")
    header = rows[0]
    body = rows[1:]
    df = pd.DataFrame(body, columns=header)
    sym = df["Gene_Symbol"].astype(str)
    valcols = [c for c in header if isinstance(c, str) and "_TenPx" in c and not c.endswith("_Peptides")]
    X = df[valcols].apply(pd.to_numeric, errors="coerce").astype(np.float64)
    X.index = sym.values
    X = X.loc[X.index.notna() & (X.index != "nan")]
    nobs = X.notna().sum(1)
    X = X.assign(_n=nobs.values).sort_values("_n", ascending=False)
    X = X[~X.index.duplicated(keep="first")].drop(columns="_n")
    ccle = [re.sub(r"_TenPx\d+$", "", c) for c in valcols]
    models = data.load_models()
    name2id = pd.Series(models.index, index=models.CCLEName).dropna()
    name2id = name2id[~name2id.index.duplicated()]
    T = X.T
    T.index = [name2id.get(c, None) for c in ccle]
    T = T[pd.Index(T.index).notna()]
    T = T.groupby(level=0).mean()
    T.columns = T.columns.astype(str)
    if cache:
        pq.parent.mkdir(parents=True, exist_ok=True)
        T.to_parquet(pq)
    return T


def _xlsx_rows(path, sheet_name: str) -> list[list]:
    """Minimal streaming xlsx reader (shared strings + numeric cells) for one sheet."""
    from lxml import etree
    ns = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    rns = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
    z = zipfile.ZipFile(path)
    wb = etree.fromstring(z.read("xl/workbook.xml"))
    rid = next(s.get(rns + "id") for s in wb.iter(ns + "sheet") if s.get("name") == sheet_name)
    rels = etree.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    target = next(r.get("Target") for r in rels if r.get("Id") == rid)
    target = "xl/" + target.lstrip("/").removeprefix("xl/")
    sst = []
    with z.open("xl/sharedStrings.xml") as fh:
        for _, si in etree.iterparse(fh, tag=ns + "si"):
            sst.append("".join(t.text or "" for t in si.iter(ns + "t")))
            si.clear()

    def colidx(ref):
        letters = re.match(r"[A-Z]+", ref).group(0)
        n = 0
        for ch in letters:
            n = n * 26 + ord(ch) - 64
        return n - 1

    out, width = [], 0
    with z.open(target) as fh:
        for _, row in etree.iterparse(fh, tag=ns + "row"):
            vals = {}
            for c in row.iter(ns + "c"):
                v = c.find(ns + "v")
                if v is None or v.text is None:
                    continue
                t = c.get("t")
                vals[colidx(c.get("r"))] = sst[int(v.text)] if t == "s" else (v.text if t == "str" else float(v.text))
            if vals:
                width = max(width, max(vals) + 1)
                out.append(vals)
            row.clear()
    return [[r.get(i) for i in range(width)] for r in out]


# ------------------------------------------------------------------------------------------------ gene annotation
def gene_annotation(cols) -> pd.DataFrame:
    """'SYMBOL (ENTREZ)' -> symbol, entrez, chrom, pos (Ensembl 116 midpoint), location / band / arm (HGNC)."""
    se = [_entrez_of(c) for c in cols]
    g = pd.DataFrame({"symbol": [s for s, _ in se], "entrez": [e for _, e in se]}, index=pd.Index(cols, name="col"))
    pos = entrez_positions()
    g["chrom"] = pos.chrom.reindex(g.entrez.fillna(-1).astype(np.int64)).fillna("").values
    g["pos"] = pos.pos.reindex(g.entrez.fillna(-1).astype(np.int64)).fillna(-1).astype(np.int64).values
    h = pd.read_csv(paths.HGNC_FILE, sep="\t", dtype=str, usecols=["entrez_id", "location"]).dropna()
    h = h.drop_duplicates("entrez_id").set_index("entrez_id").location
    g["location"] = h.reindex(g.entrez.astype("Int64").astype(str)).values
    ab = [parse_location(l) for l in g.location]
    g["arm"] = [a for a, _ in ab]
    g["band"] = [b for _, b in ab]
    x = g.chrom == "X"
    in_par = np.zeros(len(g), bool)
    for lo, hi in PAR:
        in_par |= (g.pos.values >= lo) & (g.pos.values <= hi)
    g["x_nonpar"] = x & ~in_par
    return g


def _by_entrez(df: pd.DataFrame) -> pd.DataFrame:
    """Re-key a 'SYMBOL (ENTREZ)' matrix by Entrez (drop unparsable / duplicate keys)."""
    e = [_entrez_of(c)[1] for c in df.columns]
    keep = [i for i, x in enumerate(e) if x is not None]
    out = df.iloc[:, keep]
    out.columns = pd.Index([e[i] for i in keep], name="entrez")
    return out.loc[:, ~out.columns.duplicated()]


# ----------------------------------------------------------------------------------------------- LoY threshold
def loy_threshold(male_values: np.ndarray, grid_n: int = 2048) -> dict:
    """Antimode of the bimodal male chrY-expression distribution: Gaussian KDE (Scott bandwidth); the threshold is the
    density minimum between the two highest local maxima. Raises if the density is not bimodal."""
    from scipy.stats import gaussian_kde
    v = np.asarray(male_values, float)
    v = v[~np.isnan(v)]
    kde = gaussian_kde(v)
    grid = np.linspace(v.min(), v.max(), grid_n)
    d = kde(grid)
    peaks = [i for i in range(1, grid_n - 1) if d[i] >= d[i - 1] and d[i] > d[i + 1]]
    if len(peaks) < 2:
        raise ValueError("male chrY expression is not bimodal")
    top = sorted(sorted(peaks, key=lambda i: -d[i])[:2])
    lo, hi = top
    j = lo + int(np.argmin(d[lo:hi + 1]))
    return dict(threshold=float(grid[j]), low_mode=float(grid[lo]), high_mode=float(grid[hi]),
                density_at_threshold=float(d[j]), bandwidth_factor=float(kde.factor), n_male=int(len(v)))


# ------------------------------------------------------------------------------------------------ segments
def merge_segments(E: pd.DataFrame, ann: pd.DataFrame, label: str, jaccard: float = 0.9, min_n: int = 1,
                   gene_frac: float = 0.5) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Merge adjacent genes (Ensembl order within chromosome) whose binary event vectors over models have Jaccard >=
    `jaccard` (computed on models observed for both genes; two event-free genes are not merged into an event run, i.e.
    runs are built only over genes with >= 1 event). A model carries the segment event if it has the event in >=
    gene_frac of the segment's observed genes. Returns (segment table, models x segments matrix) for segments whose
    event count >= min_n. Segments are named '<band or band range> <label>' (HGNC cytoband)."""
    ann = ann.reindex(E.columns)
    ok = (ann.pos.values >= 0) & (ann.chrom.values != "") & (np.nansum(E.values == 1, 0) >= 1)
    cols = E.columns[ok]
    a = ann.loc[cols].copy()
    a["_c"] = pd.Categorical(a.chrom, [str(i) for i in range(1, 23)] + ["X", "Y", "MT"], ordered=True)
    a = a.sort_values(["_c", "pos"])
    V = E[a.index].values
    obs = ~np.isnan(V)
    one = (V == 1)
    rows, mats = [], {}
    start = 0
    n = len(a)
    breaks = np.ones(max(n - 1, 0), bool)
    if n > 1:
        both = obs[:, :-1] & obs[:, 1:]
        inter = (one[:, :-1] & one[:, 1:] & both).sum(0)
        union = ((one[:, :-1] | one[:, 1:]) & both).sum(0)
        same_chrom = a.chrom.values[:-1] == a.chrom.values[1:]
        with np.errstate(invalid="ignore", divide="ignore"):
            J = np.where(union > 0, inter / np.maximum(union, 1), 0.0)
        breaks = ~(same_chrom & (J >= jaccard))
    bounds = [0] + [i + 1 for i in np.flatnonzero(breaks)] + [n]
    names = {}
    for s, e in zip(bounds[:-1], bounds[1:]):
        if e <= s:
            continue
        idx = a.index[s:e]
        sub_one, sub_obs = one[:, s:e], obs[:, s:e]
        nobs = sub_obs.sum(1)
        frac = np.where(nobs > 0, sub_one.sum(1) / np.maximum(nobs, 1), np.nan)
        vec = np.where(nobs > 0, (frac >= gene_frac).astype(float), np.nan)
        npos = int(np.nansum(vec))
        if npos < min_n:
            continue
        b0, b1 = a.band.iloc[s] or a.location.iloc[s], a.band.iloc[e - 1] or a.location.iloc[e - 1]
        loc0 = a.location.iloc[s] if isinstance(a.location.iloc[s], str) else b0
        loc1 = a.location.iloc[e - 1] if isinstance(a.location.iloc[e - 1], str) else b1
        locus = loc0 if loc0 == loc1 else f"{loc0}-{loc1}"
        genes_lab = a.symbol.iloc[s] if e - s == 1 else f"{a.symbol.iloc[s]}-{a.symbol.iloc[e - 1]}"
        nm = f"{locus} {label} {genes_lab}"
        names[nm] = names.get(nm, 0) + 1
        if names[nm] > 1:
            nm = f"{nm} ({names[nm]})"
        mats[nm] = vec
        rows.append(dict(segment=nm, event=label, chrom=a.chrom.iloc[s], start_pos=int(a.pos.iloc[s]),
                         end_pos=int(a.pos.iloc[e - 1]), n_genes=int(e - s), first_gene=a.symbol.iloc[s],
                         last_gene=a.symbol.iloc[e - 1], genes=";".join(a.symbol.iloc[s:e]), n_events=npos,
                         n_observed=int(np.sum(~np.isnan(vec)))))
    tab = pd.DataFrame(rows)
    mat = pd.DataFrame(mats, index=E.index)
    return tab, mat


# ------------------------------------------------------------------------------------------------ builder
@dataclass
class ContextBuild:
    table: pd.DataFrame                    # ModelID x named contexts (+ meta columns)
    definitions: dict                      # context -> (kind, definition)
    lof: pd.DataFrame                      # genome-wide LoF incl. MSI-H (models x entrez)
    cn_events: dict                        # 'homdel' / 'loss' / 'partial' / 'amp' / 'highamp' -> models x entrez
    segments: pd.DataFrame                 # segment table (all event types)
    segment_matrix: pd.DataFrame           # models x segment
    functional: pd.DataFrame               # LoF functional check per gene
    log: dict = field(default_factory=dict)


def _hk_mask(ann: pd.DataFrame) -> np.ndarray:
    """Ribosomal-protein genes (RPL*/RPS*, excluding pseudo/kinase/like genes) as the 'housekeeping' CN-scale set."""
    return ann.symbol.str.match(r"^RP[LS]\d+[AXY]?\d*$").values


def build(universe: list[str] | None = None, lof_check_min: int = 8) -> ContextBuild:
    U = universe if universe is not None else data.universe_model_ids()
    Ui = pd.Index(U, name="ModelID")
    models = data.load_models().reindex(Ui)
    L: dict = {}
    defs: dict = {}
    T: dict = {}                                         # name -> Series (assembled into a DataFrame at the end)

    def add(name, values, kind, definition):
        T[name] = pd.Series(np.asarray(values, dtype=np.float64) if kind != "categorical" else np.asarray(values),
                            index=Ui)
        defs[name] = (kind, definition)

    # ---- inputs
    sig = read_signatures(U)
    dam, dam_seq = read_omics(paths.OMICS_DAMAGING, U)
    hot, hot_seq = read_omics(paths.OMICS_HOTSPOT, U)
    cn, cn_seq = read_omics(paths.OMICS_CN, U)
    ex, ex_seq = read_omics(paths.OMICS_EXPR, U)
    L["n_models"] = dict(U=len(U), damaging=len(dam), hotspot=len(hot), cn_wgs=len(cn), expression=len(ex),
                         signatures_dna=int(sig.SequencingID.notna().sum()))
    # the CN profile and the signature (Ploidy) profile must be the same WGS sequencing run
    same = (sig.SequencingID.reindex(cn_seq.index) == cn_seq)
    L["cn_ploidy_same_profile"] = dict(n=int(len(cn_seq)), same=int(same.sum()))
    if same.mean() < 0.99:
        raise ValueError("Ploidy profile differs from the CN profile for >1% of models")
    dam_e, hot_e = _by_entrez(dam).reindex(Ui), _by_entrez(hot).reindex(Ui)
    cn_e, ex_e = _by_entrez(cn).reindex(Ui), _by_entrez(ex).reindex(Ui)
    allcols = list(dict.fromkeys(list(dam.columns) + list(cn.columns) + list(ex.columns)))
    ann = gene_annotation(allcols)
    ann = ann.dropna(subset=["entrez"]).drop_duplicates("entrez").set_index("entrez", drop=False)
    ann.index = ann.index.astype(np.int64)
    sym2e = pd.Series(ann.entrez.values.astype(np.int64), index=ann.symbol.values)
    sym2e = sym2e[~sym2e.index.duplicated()]

    import time as _t
    _t0 = _t.time()

    def stage(name):
        log.info("stage %s at %.1f s", name, _t.time() - _t0)

    stage("inputs loaded")
    # ---- signatures
    msi = sig.MSIScore.values.astype(float)
    ploidy = sig.Ploidy.values.astype(float)
    add("MSIScore", msi, "continuous", "OmicsGlobalSignatures MSIScore (default DNA profile)")
    add("MSI_H", cmp(msi, ">", 20), "binary", "MSIScore > 20")
    add("log_burden", np.log1p(dam_e.sum(1, min_count=1).values), "continuous",
        "log1p(row sum of OmicsSomaticMutationsMatrixDamaging, default DNA profile)")
    add("Ploidy", ploidy, "continuous", "OmicsGlobalSignatures Ploidy (default DNA profile)")
    add("WGD", sig.WGD.values.astype(float), "binary", "OmicsGlobalSignatures WGD (default DNA profile)")
    # sensitivity only (not primary): Ploidy / WGD borrowed from the model's other DNA profiles when the default
    # profile has none (default first, then non-default WGS, then WES; latest SequencingDate first)
    alt = read_signatures_any(U)
    add("Ploidy_anyprofile", np.where(np.isnan(ploidy), alt.Ploidy.values, ploidy), "continuous",
        "SENSITIVITY: Ploidy of the default DNA profile, else of another DNA profile of the model")
    add("WGD_anyprofile", np.where(np.isnan(sig.WGD.values.astype(float)), alt.WGD.values, sig.WGD.values), "binary",
        "SENSITIVITY: WGD of the default DNA profile, else of another DNA profile of the model")
    msi_h = T["MSI_H"].values

    stage('CN')
    # ---- CN
    r = np.minimum(cn_e.values.astype(np.float64), 10.0)
    P = ploidy[:, None]
    a = r * P
    genes_cn = cn_e.columns
    exr = ex_e.reindex(columns=genes_cn).values.astype(np.float64)
    deep = np.where(np.isnan(P), cmp(r, "<", 0.25), cmp(a, "<", 0.5))
    deep = np.where(np.isnan(r), np.nan, deep)
    lowexp = np.where(np.isnan(exr), 1.0, cmp(exr, "<", 1.0))           # expression check where measured
    homdel = k_and(deep, lowexp)
    loss_a = cmp(a - (P - 0.75), "<=", 0.0)
    partial = k_and(cmp(a, ">=", 0.5), loss_a)
    loss_any = k_or(loss_a, homdel)
    amp, high = cmp(r, ">=", 2.0), cmp(r, ">=", 3.0)
    # sex chromosomes: r is relative to the autosomal ploidy, so a male X (or Y) copy sits at r ~ 0.5 and every male
    # X gene would read as 'loss'. Male X non-PAR / Y genes: partial and any-loss events -> NaN (hom-del kept: a < 0.5
    # is a real deletion of the single copy). chrY genes in non-male lines: all events NaN.
    is_male = (models.Sex == "Male").values[:, None]
    an_cn = ann.reindex(genes_cn)
    xg = an_cn.x_nonpar.fillna(False).values.astype(bool)[None, :]
    yg_cn = (an_cn.chrom == "Y").values[None, :]
    hemi = is_male & (xg | yg_cn)
    noY = ~is_male & yg_cn
    partial = np.where(hemi | noY, np.nan, partial)
    loss_any = np.where(hemi | noY, np.nan, loss_any)
    homdel = np.where(noY, np.nan, homdel)
    amp = np.where(noY, np.nan, amp)
    high = np.where(noY, np.nan, high)
    L["sex_chromosome_rule"] = dict(n_x_nonpar_genes=int(xg.sum()), n_y_genes=int(yg_cn.sum()))
    ev = {k: pd.DataFrame(v.astype(np.float32), index=Ui, columns=genes_cn)
          for k, v in dict(homdel=homdel, partial=partial, loss=loss_any, amp=amp, highamp=high).items()}
    # CN-scale check (logged): r is ploidy-relative if the median r of housekeeping genes ~1 in WGD and non-WGD lines
    hk = _hk_mask(ann.reindex(genes_cn).fillna({"symbol": ""}))
    wgd = sig.WGD.values
    med_all = np.nanmedian(np.where(np.isnan(r), np.nan, r), axis=1)
    med_hk = np.nanmedian(r[:, hk], axis=1) if hk.any() else np.full(len(Ui), np.nan)
    has = ~np.isnan(med_all)
    L["cn_scale_check"] = {
        "n_housekeeping_genes(RPL/RPS)": int(hk.sum()),
        **{f"median_r_{lab}_{grp}": float(np.nanmedian(v[has & (wgd == w)]))
           for lab, v in (("allgenes", med_all), ("housekeeping", med_hk)) for grp, w in (("WGD", 1), ("nonWGD", 0))},
        "median_abs_CN_housekeeping_over_ploidy": float(np.nanmedian((med_hk * ploidy / ploidy)[has])),
        "n_models_with_CN": int(has.sum()), "n_models_with_CN_and_Ploidy_NaN": int((has & np.isnan(ploidy)).sum()),
        "fraction_of_cells_clipped_at_10": float(np.nanmean(cn_e.values > 10)),
    }
    log.info("CN-scale check: %s", L["cn_scale_check"])

    stage('LoF')
    # ---- LoF (genome-wide, genes in the damaging matrix)
    genes_dam = dam_e.columns
    D = dam_e.values.astype(np.float64)
    a_d = pd.DataFrame(a, index=Ui, columns=genes_cn).reindex(columns=genes_dam).values
    hd_d = ev["homdel"].reindex(columns=genes_dam).values.astype(np.float64)
    male = (models.Sex == "Male").values.astype(float)[:, None]
    xnp = ann.x_nonpar.reindex(genes_dam).fillna(False).values.astype(float)[None, :]
    lof = k_or(cmp(D, "==", 2), k_and(cmp(D, ">=", 1), cmp(a_d - (P - 0.75), "<", 0.0)), hd_d,
               k_and(male * xnp, cmp(D, ">=", 1)))
    ydam = (ann.chrom.reindex(genes_dam) == "Y").values[None, :]
    lof = np.where(ydam & (male == 0), np.nan, lof)                 # chrY genes: not applicable in non-male lines
    lof_df = pd.DataFrame(lof.astype(np.float32), index=Ui, columns=genes_dam)
    L["lof_male_x_rule_genes"] = int(xnp.sum())

    def e(sym):
        if sym not in sym2e.index:
            raise KeyError(f"gene symbol {sym} not in the Omics matrices")
        return int(sym2e[sym])

    def col(df, sym):
        k = e(sym)
        if k not in df.columns:
            raise KeyError(f"{sym} ({k}) missing")
        return df[k].values.astype(np.float64)

    for g in LOF_GENES:
        v = col(lof_df, g)
        add(f"LoF_inclMSI:{g}", v, "binary", f"LoF({g}) incl. MSI-H lines")
        prim = v if g in MMR_GENES else np.where(msi_h == 0, v, np.nan)
        add(f"LoF:{g}", prim, "binary", f"LoF({g}); MSI-H lines NaN" + (" (MMR gene: kept)" if g in MMR_GENES else ""))
    for g in HOTSPOT_GENES:
        add(f"hotspot:{g}", cmp(col(hot_e, g), ">=", 1), "binary", f"hotspot matrix {g} >= 1")
    for g in HOMDEL_GENES:
        add(f"CN_homdel:{g}", col(ev["homdel"], g), "binary", f"hom-del {g} (a<0.5 & expr<1 where measured)")
    for g in LOSS_GENES:
        add(f"CN_loss:{g}", col(ev["loss"], g), "binary", f"partial or hom loss {g} (a <= Ploidy-0.75 or hom-del)")
    for g in AMP_GENES:
        add(f"CN_amp:{g}", col(ev["amp"], g), "binary", f"amp {g} (r >= 2)")
        add(f"CN_highamp:{g}", col(ev["highamp"], g), "binary", f"high-level amp {g} (r >= 3)")

    stage('TP53 status')
    # ---- TP53 status
    d53, h53 = col(dam_e, "TP53"), T["hotspot:TP53"].values
    alt = k_or(h53, T["LoF_inclMSI:TP53"].values)
    wt = k_and(cmp(d53, "==", 0), cmp(h53, "==", 0), 1 - col(ev["homdel"], "TP53"))
    tp53 = np.where(alt == 1, 1.0, np.where(wt == 1, 0.0, np.nan))
    add("TP53_altered", tp53, "binary", "1 = TP53 hotspot or LoF; 0 = no damaging/hotspot/hom-del; NaN = ambiguous "
        "single hit or unknown")
    L["tp53"] = dict(altered=int(np.sum(tp53 == 1)), wt=int(np.sum(tp53 == 0)),
                     ambiguous_damaging_hit=int(np.sum((d53 >= 1) & np.isnan(tp53))),
                     unknown_other=int(np.sum((~(d53 >= 1)) & np.isnan(tp53))))

    # ---- composite mutation contexts
    ras = k_or(T["hotspot:KRAS"].values, T["hotspot:NRAS"].values)
    add("RAS_hotspot", ras, "binary", "KRAS or NRAS hotspot (PC4)")
    braf = T["hotspot:BRAF"].values
    pc5 = np.where((braf == 1) & (ras == 0), 1.0, np.where((ras == 1) & (braf == 0), 0.0, np.nan))
    add("BRAF_vs_RAS_hotspot", pc5, "binary", "1 = BRAF hotspot (no RAS hotspot), 0 = KRAS/NRAS hotspot (no BRAF); "
        "others NaN (PC5)")
    nfe = T["hotspot:NFE2L2"].values
    nrf2_all = k_or(T["LoF_inclMSI:KEAP1"].values, T["LoF_inclMSI:CUL3"].values, nfe)
    nrf2 = np.where(msi_h == 0, nrf2_all, np.where(nfe == 1, 1.0, np.nan))
    add("NRF2_active", nrf2, "binary", "LoF KEAP1 or LoF CUL3 or NFE2L2 hotspot; MSI-H lines 1 only via hotspot, else NaN")
    pten, pik = T["LoF:PTEN"].values, T["hotspot:PIK3CA"].values
    add("PTEN_LoF_noPIK3CAhs", np.where(pik == 1, np.nan, pten), "binary", "LoF PTEN; lines with PIK3CA hotspot NaN (PC19a)")
    add("PIK3CA_hotspot_noPTENLoF", np.where(pten == 1, np.nan, pik), "binary",
        "PIK3CA hotspot; lines with PTEN LoF NaN (PC19b)")

    # ---- CN composite contexts
    add("9p21_loss", T["CN_homdel:MTAP"].values, "binary", "9p21 loss = MTAP hom-del (PC9; segment reported)")
    add("FOCAD_del_or_MSI_H", k_or(T["CN_homdel:FOCAD"].values, msi_h), "binary", "FOCAD hom-del or MSI-H (PC10)")
    add("VPS4B_loss", T["CN_loss:VPS4B"].values, "binary", "VPS4B partial or hom loss (PC12)")
    add("MYC_amp", T["CN_amp:MYC"].values, "binary", "MYC r >= 2")
    add("MYCN_amp", T["CN_amp:MYCN"].values, "binary", "MYCN r >= 2")

    stage('fusions')
    # ---- fusions
    fus = read_fusions(U)
    rna = pd.Series(np.isin(Ui, ex.index), index=Ui)
    fus = fus[fus.Gene1_symbol != fus.Gene2_symbol]
    pairs = fus.apply(lambda x: frozenset((x.Gene1_symbol, x.Gene2_symbol)), axis=1)
    fcounts = {}
    for name, plist in FUSIONS.items():
        want = {frozenset(p) for p in plist}
        hit = set(fus.ModelID[pairs.isin(want)])
        v = np.where(rna.values, np.isin(Ui, list(hit)).astype(float), np.nan)
        add(f"fusion:{name}", v, "binary", "unordered " + ", ".join("::".join(p) for p in plist)
            + " (OmicsFusionFiltered default RNA profile; no RNA -> NaN)")
        fcounts[name] = int(np.nansum(v))
        hit_no_rna = len(hit - set(Ui[rna.values]))
        if hit_no_rna:
            raise ValueError(f"fusion {name}: {hit_no_rna} hits in models without a default RNA expression profile")
    L["fusion_counts_U"] = fcounts
    for name, mn in FUSION_MIN_U.items():
        if fcounts[name] == 0 or fcounts[name] < mn:
            raise AssertionError(f"fusion {name}: {fcounts[name]} models in U < expected minimum {mn}")
    smb = T["LoF:SMARCB1"].values
    add("SMARCB1_LoF_or_SS18SSX", k_or(smb, T["fusion:SS18::SSX"].values), "binary", "LoF SMARCB1 or SS18::SSX (PC15/PC22)")
    npm1c = cmp(col(dam_e, "NPM1"), ">=", 1)
    kmt = T["fusion:KMT2A-r"].values
    add("KMT2A_r_vs_nonNPM1c", np.where((kmt == 0) & (npm1c == 1), np.nan, kmt), "binary",
        "KMT2A-r; reference excludes NPM1-mutant (damaging >= 1) lines (PC16)")

    stage('sex, LoY')
    # ---- sex, LoY
    sex = models.Sex
    female = np.where(sex == "Female", 1.0, np.where(sex == "Male", 0.0, np.nan))
    add("sex_female", female, "binary", "Model.csv Sex: Female 1, Male 0, Unknown NaN")
    yg = ann.index[(ann.chrom == "Y").values]
    ygs = set(yg)
    ycols = [c for c in ex_e.columns if c in ygs]
    ysum = ex_e[ycols].sum(1, min_count=1).values.astype(float)
    add("chrY_expr_sum", ysum, "continuous", f"sum of log1p TPM over {len(ycols)} chrY (male-specific) protein-coding genes")
    thr = loy_threshold(ysum[(female == 0) & ~np.isnan(ysum)])
    fem_v = ysum[(female == 1) & ~np.isnan(ysum)]
    thr.update(n_chrY_genes=len(ycols), chrY_genes=sorted(ann.symbol.reindex(ycols).tolist()),
               female_p99=float(np.percentile(fem_v, 99)), female_max=float(fem_v.max()),
               n_male_below=int(np.sum((female == 0) & (ysum < thr["threshold"]))),
               n_male_above=int(np.sum((female == 0) & (ysum >= thr["threshold"]))),
               n_female_above=int(np.sum((female == 1) & (ysum >= thr["threshold"]))))
    if "DDX3Y" in set(ann.symbol):
        ddx3y = ex_e[e("DDX3Y")].values
        thr["agreement_with_DDX3Y_lt_1_in_males"] = float(np.mean(
            ((ysum < thr["threshold"]) == (ddx3y < 1))[(female == 0) & ~np.isnan(ysum)]))
    L["loy"] = thr
    loy = np.where(female == 0, cmp(ysum, "<", thr["threshold"]), np.nan)
    add("LoY", loy, "binary", f"male lines with chrY_expr_sum < {thr['threshold']:.3f} (KDE antimode); female/unknown NaN")
    add("female_or_LoY", np.where(female == 1, 1.0, np.where(female == 0, loy, np.nan)), "binary",
        "female line or LoY male line (PC13); intact-Y males 0")

    # ---- lineage
    lin = models.OncotreeLineage
    add("OncotreeLineage", lin.values, "categorical", "Model.csv OncotreeLineage")
    add("PatientID", models.PatientID.values, "categorical", "Model.csv PatientID (blocking unit)")
    for l in sorted(lin.dropna().unique()):
        add(f"lineage:{l}", np.where(lin.isna(), np.nan, (lin == l).astype(float)), "binary", f"OncotreeLineage == {l}")
    heme = np.where(lin.isna(), np.nan, lin.isin(["Myeloid", "Lymphoid"]).astype(float))
    add("heme", heme, "binary", "OncotreeLineage Myeloid or Lymphoid (PC20)")
    add("heme_vs_breast", np.where(heme == 1, 1.0, np.where(lin == "Breast", 0.0, np.nan)), "binary",
        "1 = Myeloid/Lymphoid, 0 = Breast, others NaN (PC21)")
    pdis = models.OncotreePrimaryDisease
    add("melanoma", np.where(pdis.isna(), np.nan, pdis.str.contains("Melanoma", na=False).astype(float)), "binary",
        "OncotreePrimaryDisease contains 'Melanoma' (any site; PC17)")

    stage('segments')
    # ---- segments (genome-wide CN events)
    segs, smats = [], []
    for k in ("homdel", "loss", "amp", "highamp"):
        tab, mat = merge_segments(ev[k], ann, k, jaccard=0.9, min_n=1)
        segs.append(tab)
        smats.append(mat)
    seg_tab = pd.concat(segs, ignore_index=True)
    seg_mat = pd.concat(smats, axis=1)
    for g, k in (("MTAP", "homdel"), ("CDKN2A", "homdel"), ("FOCAD", "homdel"), ("VPS4B", "loss"), ("MYC", "amp"),
                 ("MYCN", "amp")):
        hit = seg_tab[(seg_tab.event == k) & seg_tab.genes.str.split(";").apply(lambda s: g in s)]
        L.setdefault("segment_of", {})[f"{g} {k}"] = (hit[["segment", "n_genes", "n_events", "first_gene", "last_gene"]]
                                                      .to_dict("records"))
    L["segments"] = {k: dict(n_segments=int((seg_tab.event == k).sum()),
                             n_ge8=int(((seg_tab.event == k) & (seg_tab.n_events >= 8)).sum()),
                             n_ge8_multigene=int(((seg_tab.event == k) & (seg_tab.n_events >= 8)
                                                  & (seg_tab.n_genes > 1)).sum()),
                             median_genes_per_segment=float(seg_tab.n_genes[seg_tab.event == k].median()))
                     for k in ("homdel", "loss", "amp", "highamp")}

    stage('LoF functional check')
    # ---- LoF functional check
    table = pd.DataFrame(T, index=Ui)
    func = lof_functional_check(lof_df, table, ex_e, sym2e, ann, min_lof=lof_check_min)
    L["lof_functional_check"] = dict(n_genes=int(len(func)), n_pass=int(func["pass"].sum()),
                                     named_fail=func.loc[func.named & ~func["pass"], "gene"].tolist())
    return ContextBuild(table, defs, lof_df, ev, seg_tab, seg_mat, func, L)


def lof_functional_check(lof_df, T, ex_e, sym2e, ann, min_lof: int = 8, alpha: float = 0.05) -> pd.DataFrame:
    """Pre-registered LoF functional check (DESIGN §2.2): for each gene with >= min_lof MSS LoF lines in U (plus every
    named LoF gene), compare the gene's own expression (log1p TPM) and protein (CCLE proteomics, where measured) in
    MSS LoF lines vs MSS non-LoF lines. One-sided Mann-Whitney U (LoF lower). Pass = (expression p < alpha and median
    difference < 0) OR (protein p < alpha, median difference < 0, >= 3 LoF lines with protein). Never uses dependency."""
    from scipy.stats import mannwhitneyu
    mss = T["MSI_H"].values == 0
    L = lof_df.values
    nlof = np.nansum((L == 1) & mss[:, None], 0)
    named = {int(sym2e[g]) for g in LOF_GENES if g in sym2e.index}
    cand = [c for c, n in zip(lof_df.columns, nlof) if n >= min_lof or c in named]
    prot = read_proteomics()
    prot = prot.reindex(T.index)
    sym = ann.symbol
    rows = []
    for g in cand:
        s = sym.get(g, str(g))
        v = lof_df[g].values
        lo, ref = (v == 1) & mss, (v == 0) & mss
        rec = dict(gene=s, entrez=int(g), named=g in named, n_lof_mss=int(lo.sum()), n_ref_mss=int(ref.sum()))
        for lab, M in (("expr", ex_e.get(g)), ("prot", prot.get(s))):
            if M is None:
                rec.update({f"{lab}_n_lof": 0, f"{lab}_n_ref": 0, f"{lab}_meddiff": np.nan, f"{lab}_p": np.nan})
                continue
            x = M.values.astype(float)
            a_, b_ = x[lo & ~np.isnan(x)], x[ref & ~np.isnan(x)]
            p = mannwhitneyu(a_, b_, alternative="less").pvalue if len(a_) >= 2 and len(b_) >= 2 else np.nan
            rec.update({f"{lab}_n_lof": len(a_), f"{lab}_n_ref": len(b_),
                        f"{lab}_meddiff": float(np.median(a_) - np.median(b_)) if len(a_) and len(b_) else np.nan,
                        f"{lab}_p": float(p)})
        ok_e = rec["expr_p"] < alpha and rec["expr_meddiff"] < 0
        ok_p = rec["prot_n_lof"] >= 3 and rec["prot_p"] < alpha and rec["prot_meddiff"] < 0
        rec.update(expr_pass=bool(ok_e), prot_pass=bool(ok_p), **{"pass": bool(ok_e or ok_p)})
        rows.append(rec)
    return pd.DataFrame(rows).sort_values(["named", "n_lof_mss"], ascending=[False, False]).reset_index(drop=True)


def summary(build: ContextBuild) -> pd.DataFrame:
    rows = []
    for c, (kind, d) in build.definitions.items():
        v = build.table[c]
        if kind == "categorical":
            rows.append(dict(context=c, kind=kind, n_U=len(v), n_nonNaN=int(v.notna().sum()), n_pos="",
                             definition=d))
            continue
        x = v.values.astype(float)
        rows.append(dict(context=c, kind=kind, n_U=len(x), n_nonNaN=int((~np.isnan(x)).sum()),
                         n_pos=int(np.nansum(x == 1)) if kind == "binary" else "", definition=d))
    return pd.DataFrame(rows)


def content_hash(df: pd.DataFrame) -> str:
    import hashlib
    h = hashlib.sha256()
    h.update(json.dumps([str(c) for c in df.columns]).encode())
    h.update(json.dumps([str(i) for i in df.index]).encode())
    h.update(pd.util.hash_pandas_object(df, index=True).values.tobytes())
    return h.hexdigest()


def load_table(path=None) -> pd.DataFrame:
    return pd.read_parquet(path or paths.CONTEXT_CACHE / "context_table.parquet")
