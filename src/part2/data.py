"""DepMap data layer: gene effect (models x genes, Entrez-keyed), Model.csv, CRISPRConfounders, model subsets,
PatientID-blocked subsamples / splits.

Gene keys: the DepMap column 'SYMBOL (ENTREZ)'; the Part 2 gene universe is Part 1's protein-gene table
(genes.tsv: index, symbol, entrez, ..., low_coverage, depmap_column, depmap_index), so every Part 2 array indexed
by gene uses genes.tsv ``index`` order (= DepMap column order restricted to protein genes).

The gene-effect matrix is parsed exactly as Part 1 (``coessential.data.load_depmap``: pandas, float32) and cached
as .npy in data_cache/part2 keyed by the CSV's SHA-256. Subsets keep the CSV row order whatever order model IDs
are given in (row order changes floating-point results), so a subset is a function of the *set* of model IDs.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
import pandas as pd

from part2 import paths


def sha256_file(path, chunk: int = 1 << 24) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def model_set_hash(model_ids) -> str:
    """SHA-256 (first 16 hex) of the sorted, newline-joined model IDs."""
    ids = sorted(set(map(str, model_ids)))
    if len(ids) != len(list(model_ids)):
        raise ValueError("duplicate model IDs")
    return hashlib.sha256("\n".join(ids).encode()).hexdigest()[:16]


@dataclass
class GeneEffect:
    X: np.ndarray            # float32, models x all DepMap genes (CSV order)
    models: pd.Index         # ModelID, CSV row order
    columns: pd.Index        # 'SYMBOL (ENTREZ)'
    csv_sha256: str

    @property
    def entrez(self) -> np.ndarray:
        return np.array([c.rsplit("(", 1)[1].rstrip(")") for c in self.columns], dtype=np.int64)

    def rows(self, model_ids) -> np.ndarray:
        """Row indices (ascending CSV order) of the given model IDs; error on unknown / duplicate IDs."""
        ids = pd.Index(list(map(str, model_ids)))
        if ids.has_duplicates:
            raise ValueError("duplicate model IDs")
        miss = ids.difference(self.models)
        if len(miss):
            raise KeyError(f"{len(miss)} model IDs not in CRISPRGeneEffect: {list(miss[:5])}")
        return np.sort(self.models.get_indexer(ids))

    def subset(self, model_ids) -> np.ndarray:
        """Models x all-genes float32 matrix for the given model set (CSV row order)."""
        r = self.rows(model_ids)
        return self.X if len(r) == len(self.models) else np.ascontiguousarray(self.X[r])


@lru_cache(maxsize=1)
def load_gene_effect(use_cache: bool = True) -> GeneEffect:
    csv = paths.GENE_EFFECT_CSV
    h = sha256_file(csv)
    npy = paths.CACHE / f"CRISPRGeneEffect.{h[:16]}.f32.npy"
    meta = npy.with_suffix(".json")
    if use_cache and npy.exists() and meta.exists():
        m = json.loads(meta.read_text())
        return GeneEffect(np.load(npy), pd.Index(m["models"]), pd.Index(m["columns"]), h)
    df = pd.read_csv(csv, index_col=0)                     # == coessential.data.load_depmap
    X = df.values.astype(np.float32)
    ge = GeneEffect(X, pd.Index(df.index.astype(str)), df.columns, h)
    if use_cache:
        paths.CACHE.mkdir(parents=True, exist_ok=True)
        np.save(npy, X)
        meta.write_text(json.dumps(dict(csv=str(csv), sha256=h, models=list(ge.models), columns=list(ge.columns))))
    return ge


@lru_cache(maxsize=1)
def protein_genes() -> pd.DataFrame:
    """Part 1 genes.tsv (protein-gene universe; low_coverage = Part 1's NaN >= 10% flag on all models)."""
    return pd.read_csv(paths.GENES_TSV, sep="\t")


@lru_cache(maxsize=1)
def load_models() -> pd.DataFrame:
    return pd.read_csv(paths.MODEL_CSV, low_memory=False).set_index("ModelID")


@lru_cache(maxsize=1)
def load_confounders() -> pd.DataFrame:
    return pd.read_csv(paths.CONFOUNDERS_CSV, index_col=0)


@lru_cache(maxsize=1)
def crispr_model_ids() -> tuple[str, ...]:
    """ModelIDs of CRISPRGeneEffect.csv in CSV row order, reading ONLY the first (ModelID) column: no gene-effect
    value is parsed. Used by context / lockbox code, which must not touch dependency values."""
    return tuple(pd.read_csv(paths.GENE_EFFECT_CSV, usecols=[0]).iloc[:, 0].astype(str))


def universe_model_ids() -> list[str]:
    """U (same rule as :func:`universe_models`) from the ModelID column only (no effect values read)."""
    ids = pd.Index(crispr_model_ids())
    m = load_models().reindex(ids)
    if m.index.isna().any() or len(m) != len(ids):
        raise KeyError("CRISPR models missing from Model.csv")
    eng = m.EngineeredModel.notna() | m.CulturedResistanceDrug.notna()
    return list(ids[~eng.values])


def universe_models() -> list[str]:
    """U = CRISPRGeneEffect models minus engineered / drug-adapted models (Model.csv EngineeredModel non-empty;
    covers all CulturedResistanceDrug models). DESIGN §2.1: 1,361 - 16 = 1,345."""
    ge = load_gene_effect()
    m = load_models().reindex(ge.models)
    if m.index.isna().any():
        raise KeyError("CRISPR models missing from Model.csv")
    eng = m.EngineeredModel.notna() | m.CulturedResistanceDrug.notna()
    return list(ge.models[~eng.values])


# ------------------------------------------------------------------------------------------ blocked sampling
def _groups(model_ids, block_col: str = "PatientID") -> pd.Series:
    m = load_models().reindex(list(model_ids))
    g = m[block_col].astype(object)
    g = g.where(g.notna(), pd.Series(["__model__" + i for i in m.index], index=m.index))   # no patient: own block
    return g


def blocked_subsample(model_ids, frac: float, seed: int, block_col: str = "PatientID",
                      blocks: pd.Series | None = None) -> list[str]:
    """Subsample ~frac of the models by whole blocks (all models of a patient together). Blocks are visited in a
    seeded random order and added while the running model count stays <= round(frac * n) (blocks that would
    overshoot are skipped), so the result has at most round(frac * n) models and is a deterministic function of
    (model set, frac, seed). Returned in sorted order."""
    ids = sorted(map(str, model_ids))
    g = (blocks.reindex(ids) if blocks is not None else _groups(ids, block_col))
    sizes = g.value_counts().sort_index()
    order = np.random.default_rng(seed).permutation(len(sizes))
    target = int(round(frac * len(ids)))
    take, n = [], 0
    for i in order:
        if n + sizes.iloc[i] <= target:
            take.append(sizes.index[i])
            n += int(sizes.iloc[i])
        if n == target:
            break
    return sorted(g.index[g.isin(take)])


def blocked_split(model_ids, frac: float, seed: int, strata: pd.Series | None = None,
                  block_col: str = "PatientID", blocks: pd.Series | None = None) -> tuple[list[str], list[str]]:
    """Two-way split (first part ~frac) by whole blocks, optionally stratified: within each stratum (e.g.
    OncotreeLineage; a block is assigned to the stratum of its first model in sorted order) blocks are shuffled
    with rng(seed) and assigned to part 1 until the stratum's part-1 count reaches round(frac * stratum size)."""
    ids = sorted(map(str, model_ids))
    g = (blocks.reindex(ids) if blocks is not None else _groups(ids, block_col))
    st = (strata.reindex(ids).fillna("__NA__") if strata is not None else pd.Series("all", index=ids))
    blk_stratum = st.groupby(g).first()
    blk_size = g.value_counts()
    rng = np.random.default_rng(seed)
    part1 = set()
    for s in sorted(blk_stratum.unique()):
        bl = np.array(sorted(blk_stratum.index[blk_stratum == s]))
        bl = bl[rng.permutation(len(bl))]
        tot = int(blk_size[bl].sum())
        target = round(frac * tot)
        n = 0
        for b in bl:
            if n >= target:
                break
            part1.add(b)
            n += int(blk_size[b])
    a = sorted(g.index[g.isin(part1)])
    b = sorted(set(ids) - set(a))
    return a, b
