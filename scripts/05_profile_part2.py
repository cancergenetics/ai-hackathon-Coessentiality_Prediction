#!/usr/bin/env python
"""Part 2 data profile (descriptive only; no labels, no huMAP files, no curated complex DBs).

Reads (read-only):
  datasets/DepMAP/{CRISPRGeneEffect,ScreenGeneEffect(first column only),CRISPRConfounders,Model,ModelCondition}.csv
  datasets/HGNC/hgnc_complete_set.txt            (not needed; gene universe comes from Part 1 genes.tsv)
  <Part1>/outputs/all_pairs/{manifest.json,genes.tsv,zavg_BEF.npy (mmap),top_pairs.parquet}
Writes:
  results/profile/part2_profile.json, results/profile/part2_profile.md
Run:  .venv/bin/python scripts/05_profile_part2.py
"""
import json, os, re, resource, sys, time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

SEED = 0
np.random.seed(SEED)

ROOT = Path(__file__).resolve().parents[1]
DM = ROOT / "datasets" / "DepMAP"
P1 = Path(os.environ.get("PART1_DIR", ROOT.parent / "coessentiality_protein complexes")) / "outputs" / "all_pairs"
OUT = ROOT / "results" / "profile"
OUT.mkdir(parents=True, exist_ok=True)

CUTS = [3, 4, 5, 6, 8]
PE_FRAC = 0.9      # pan-essential: fraction of models with effect < -0.5 >= 0.9
EFF_THR = -0.5
MIN_LIN = (15, 30)
BLOCK = 1024

J = {"seed": SEED, "cutoffs": CUTS, "definitions": {
    "pan_essential": f"fraction of models with gene effect < {EFF_THR} >= {PE_FRAC}",
    "pan_essential_alt": "mean gene effect < -0.5",
}}
TIM, MEM = {}, {}


def rss_gb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / (1024 ** 3 if sys.platform == "darwin" else 1024 ** 2)


class step:
    def __init__(self, name): self.name = name
    def __enter__(self): self.t = time.perf_counter(); print(f"[{self.name}] ...", flush=True)
    def __exit__(self, *a):
        TIM[self.name] = round(time.perf_counter() - self.t, 2)
        MEM[self.name] = round(rss_gb(), 2)
        print(f"[{self.name}] {TIM[self.name]}s  peakRSS={MEM[self.name]}GB", flush=True)


def vc(s, n=None):
    v = s.fillna("NA").value_counts()
    return {str(k): int(x) for k, x in (v.head(n) if n else v).items()}


def eta2(values, groups):
    """one-way eta^2 (between SS / total SS), NaN-aware over rows."""
    ok = ~np.isnan(values)
    v, g = values[ok], np.asarray(groups)[ok]
    if len(v) < 3: return np.nan
    gm = v.mean(); tot = ((v - gm) ** 2).sum()
    if tot == 0: return np.nan
    bet = 0.0
    for u in np.unique(g):
        x = v[g == u]; bet += len(x) * (x.mean() - gm) ** 2
    return float(bet / tot)


# ------------------------------------------------------------------ 0. load
with step("load_depmap_matrix"):
    ge = pd.read_csv(DM / "CRISPRGeneEffect.csv", index_col=0, dtype={0: str})
    ge = ge.astype(np.float32)
    X = ge.values
    model_ids = ge.index.to_numpy()
    cols = ge.columns.to_numpy()
    J["crispr_matrix"] = {"n_models": int(X.shape[0]), "n_genes": int(X.shape[1]),
                          "nan_fraction": float(np.isnan(X).mean()),
                          "genes_with_any_nan": int(np.isnan(X).any(0).sum())}
    screens = pd.read_csv(DM / "ScreenGeneEffect.csv", usecols=[0]).iloc[:, 0]
    scr_model = screens.str.extract(r"SC-(\d+)")[0].map(lambda x: f"ACH-{x}")
    J["screen_gene_effect"] = {"n_screens": int(len(screens)),
                               "n_unique_models_from_SC_id": int(scr_model.nunique()),
                               "frac_SC_id_mapping_to_CRISPR_model": float(scr_model.isin(model_ids).mean()),
                               "screens_per_model_counts": {str(k): int(v) for k, v in scr_model.value_counts().value_counts().items()},
                               "note": "rows are screens (SC-xxxxxx.CDnn); not linkable to ModelCondition IDs, only to models via numeric id"}

# ------------------------------------------------------------------ 1. models
with step("models"):
    M = pd.read_csv(DM / "Model.csv")
    J["model_csv"] = {"n_rows": int(len(M)), "n_columns": int(M.shape[1])}
    Mc = M[M.ModelID.isin(model_ids)].set_index("ModelID").loc[model_ids].copy()
    J["models"] = {"n_crispr": int(len(Mc)), "n_crispr_in_model_csv": int(len(Mc)),
                   "n_model_csv_without_crispr": int(len(M) - len(Mc)),
                   "n_duplicate_patient": int(Mc.PatientID.duplicated().sum())}
    lin = Mc.OncotreeLineage.fillna("NA")
    lc = lin.value_counts()
    J["models"]["by_lineage"] = {k: int(v) for k, v in lc.items()}
    J["models"]["n_lineages"] = int(len(lc))
    for t in MIN_LIN:
        J["models"][f"lineages_ge{t}"] = {k: int(v) for k, v in lc[lc >= t].items()}
        J["models"][f"n_models_in_lineages_ge{t}"] = int(lc[lc >= t].sum())
    pd_ = Mc.OncotreePrimaryDisease.fillna("NA").value_counts()
    st = Mc.OncotreeSubtype.fillna("NA").value_counts()
    J["models"]["n_primary_disease"] = int(len(pd_)); J["models"]["n_subtype"] = int(len(st))
    J["models"]["primary_disease_ge15"] = {k: int(v) for k, v in pd_[pd_ >= 15].items()}
    J["models"]["primary_disease_ge30"] = {k: int(v) for k, v in pd_[pd_ >= 30].items()}
    J["models"]["n_primary_disease_ge10"] = int((pd_ >= 10).sum())
    J["models"]["subtype_ge10"] = {k: int(v) for k, v in st[st >= 10].items()}
    J["models"]["n_subtype_ge10"] = int((st >= 10).sum())
    J["models"]["n_subtype_ge15"] = int((st >= 15).sum())
    # per-lineage: primary-disease/subtype counts within lineages >=15
    J["models"]["subtypes_within_lineage_ge15"] = {
        l: {k: int(v) for k, v in Mc[lin == l].OncotreeSubtype.value_counts().head(6).items()} for l in lc[lc >= 15].index}

    # annotation fields
    ann = {}
    for c in ["Sex", "PrimaryOrMetastasis", "GrowthPattern", "AgeCategory", "PediatricModelType", "ModelType",
              "EngineeredModel", "PatientRace", "Stage", "PatientTreatmentStatus", "SampleCollectionSite"]:
        ann[c] = vc(Mc[c], 12)
    ann["Age"] = {"n_nonmissing": int(Mc.Age.notna().sum()), **{k: float(v) for k, v in Mc.Age.describe().items()}}
    ann["SourceType_top"] = vc(Mc.SourceType, 10)
    ann["missing_fraction_crispr_models"] = {c: round(float(Mc[c].isna().mean()), 3) for c in Mc.columns
                                            if Mc[c].isna().mean() > 0.0}
    ann["EngineeredModelDetails"] = vc(Mc.EngineeredModelDetails.dropna())
    J["model_annotations"] = ann
    # Sex x lineage balance, metastasis x lineage (lineages >=30)
    big = lc[lc >= 30].index
    J["models"]["sex_by_lineage_ge30"] = {l: vc(Mc[lin == l].Sex) for l in big}
    J["models"]["primary_vs_met_by_lineage_ge30"] = {l: vc(Mc[lin == l].PrimaryOrMetastasis) for l in big}
    J["models"]["pediatric_by_lineage_ge30"] = {l: {"pediatric": int((Mc[lin == l].PediatricModelType.astype(str) == "True").sum()),
                                                    "n": int((lin == l).sum())} for l in big}

    # ---- molecular annotation inside Model.csv
    msf = Mc.ModelSubtypeFeatures.fillna(""); psf = Mc.PatientSubtypeFeatures.fillna("")
    feat = (msf + " | " + psf).str.strip(" |")
    has_any = feat != ""
    mol = {"n_ModelSubtypeFeatures": int((msf != "").sum()), "n_PatientSubtypeFeatures": int((psf != "").sum()),
           "n_any_feature": int(has_any.sum()), "frac_any_feature": round(float(has_any.mean()), 3),
           "n_distinct_ModelSubtypeFeatures": int(msf[msf != ""].nunique()),
           "n_distinct_PatientSubtypeFeatures": int(psf[psf != ""].nunique()),
           "ModelSubtypeFeatures_top": vc(msf[msf != ""], 40),
           "PatientSubtypeFeatures_top": vc(psf[psf != ""], 15)}
    nb = lambda pat: r"(?<![A-Za-z0-9])(?:" + pat + r")(?![A-Za-z0-9])"
    kw = {
        "KRAS": nb("KRAS"), "NRAS": nb("NRAS"), "HRAS": nb("HRAS"), "BRAF": nb("BRAF"), "EGFR": nb("EGFR"),
        "PIK3CA": nb("PIK3CA"), "PTEN": nb("PTEN"), "TP53": nb("TP53"), "CDKN2A": nb("CDKN2A"), "RB1": nb("RB1"),
        "IDH1/2": nb("IDH1|IDH2"), "MYC (not MYCN)": nb("MYC"), "MYCN": nb("MYCN"), "ERBB2/HER2": nb("ERBB2|HER2|HER2\\+"),
        "ER/PR (ER+, PR+)": nb("ER\\+|PR\\+"), "TNBC": nb("TNBC"), "MSI": nb("MSI"), "MGMT": nb("MGMT"),
        "NF1/NF2": nb("NF1|NF2"), "APC/CTNNB1": nb("APC|CTNNB1"), "KMT2A (MLL)": nb("KMT2A|MLL"),
        "BCR-ABL1": nb("BCR(?:::|-)ABL1?"), "EWSR1/EWS fusion": nb("EWS(?:R1)?(?:::|-)(?:FLI1?|ERG)"),
        "'::' fusion (any)": r"::", "amplified/amplification": r"[Aa]mplif", "deletion/del/loss": r"[Dd]elet|\(del\)|[Ll]oss",
        "mutation/mutant (text)": r"[Mm]utat|[Mm]utant",
    }
    kwc = {}
    for k, p in kw.items():
        hit = feat.str.contains(p, regex=True)
        kwc[k] = {"n_models": int(hit.sum()), "lineages": vc(lin[hit], 4)}
    mol["keyword_counts_in_feature_fields"] = kwc
    # keywords in subtype / disease names (derived label proxies)
    nm = (Mc.OncotreeSubtype.fillna("") + " | " + Mc.OncotreePrimaryDisease.fillna(""))
    nmk = {}
    for k, p in {"IDH-mutant": r"IDH-[Mm]utant", "IDH-wildtype": r"IDH-[Ww]ildtype", "BCR-ABL1+": r"BCR-ABL1\+",
                 "ALK positive/negative": r"ALK (?:Positive|Negative)", "NUT": r"NUT ", "HER2/ER/TNBC in name": r"HER2|ER\+|Triple|TNBC",
                 "MSI/MSS in name": r"MSI|Microsatellite", "p53/MYCN in name": r"MYCN|p53", "Ewing": r"Ewing",
                 "Melanoma (BRAF/NRAS-driven proxy only)": r"Melanoma"}.items():
        nmk[k] = int(nm.str.contains(p, regex=True).sum())
    mol["keyword_counts_in_subtype_names"] = nmk
    # token frequency: gene-like tokens in feature fields
    toks = pd.Series([t for s in feat[has_any] for t in set(re.findall(r"(?<![A-Za-z0-9])[A-Z][A-Z0-9]{2,}(?:::[A-Z][A-Z0-9]+)?(?![a-z])", s))])
    mol["gene_like_token_top"] = {k: int(v) for k, v in toks.value_counts().head(30).items()}
    mol["n_tokens_ge10_models"] = int((toks.value_counts() >= 10).sum())
    mol["verdict"] = ("Model.csv contains no per-gene mutation/copy-number/expression matrix; only free-text subtype features "
                      "for a minority of models (see counts). Usable binary contexts with >=10 models are limited to a handful "
                      "(MSI, HER2+, TNBC, MYC-amplified, BCR::ABL1, EWS fusion, MGMT status, IDH-mutant/GBM subtype names), "
                      "mostly nested within single lineages.")
    J["molecular_annotation"] = mol

# ------------------------------------------------------------------ 2. confounders / conditions
with step("confounders"):
    C = pd.read_csv(DM / "CRISPRConfounders.csv", index_col=0)
    Cc = C.loc[model_ids]
    conf = {"n_rows_total": int(len(C)), "n_rows_crispr_models": int(len(Cc)),
            "all_crispr_models_present": bool(set(model_ids) <= set(C.index)),
            "n_columns": int(C.shape[1]), "columns": list(C.columns)}
    conf["missing_fraction_all_rows"] = {c: round(float(C[c].isna().mean()), 3) for c in C.columns if C[c].isna().any()}
    conf["missing_fraction_crispr_models"] = {c: round(float(Cc[c].isna().mean()), 3) for c in C.columns}
    libc = ["LibraryAvana", "LibraryHumagneCD", "LibraryKY", "LibraryTKOv3", "LibraryBrunello"]
    conf["library_counts"] = {c: int(Cc[c].sum()) for c in libc}
    conf["screentype_counts"] = {c: int(Cc[c].sum()) for c in ["ScreenType2DS", "ScreenType3DO"]}
    conf["library_onehot_rowsum_values"] = {str(k): int(v) for k, v in Cc[libc].sum(axis=1).value_counts().items()}
    # pan-essential-strength per model computed later; placeholders
    lib_lab = Cc[libc].idxmax(1).str.replace("Library", "")
    lib_lab[Cc[libc].sum(axis=1) == 0] = "none"
    conf["lineage_library_cramers_v_all"] = None
    ct = pd.crosstab(lin, lib_lab)
    ctb = ct.loc[lc[lc >= 15].index]
    chi2, p, dof, _ = stats.chi2_contingency(ctb.values)
    conf["lineage_ge15_x_library"] = {"chi2": float(chi2), "p": float(p), "dof": int(dof),
                                      "cramers_v": float(np.sqrt(chi2 / (ctb.values.sum() * (min(ctb.shape) - 1))))}
    conf["library_share_by_lineage_ge15"] = {l: {k: round(float(v), 2) for k, v in (ctb.loc[l] / ctb.loc[l].sum()).items()} for l in ctb.index}
    ct2 = pd.crosstab(lin, Cc.ScreenType3DO.astype(int))
    conf["frac_3DO_by_lineage_ge15"] = {l: round(float(ct2.loc[l].get(1, 0) / ct2.loc[l].sum()), 2) for l in lc[lc >= 15].index}
    conts = ["ScreenNNMD", "ScreenROCAUC", "ScreenFPR", "ScreenMedianEssentialDepletion", "ScreenMedianNonessentialDepletion",
             "ScreenMADEssentials", "ScreenMADNonessentials", "CasActivity", "ScreenDoublingTime"]
    big15 = lin.isin(lc[lc >= 15].index).values
    conf["eta2_lineage_ge15"] = {c: round(eta2(Cc[c].values[big15].astype(float), lin.values[big15]), 3) for c in conts}
    conf["eta2_lineage_ge15"].update({c: round(eta2(Cc[c].values[big15].astype(float), lin.values[big15]), 3) for c in libc[:3] + ["ScreenType3DO"]})
    conf["spearman_NNMD_vs"] = {c: round(float(stats.spearmanr(Cc.ScreenNNMD, Cc[c], nan_policy="omit")[0]), 2)
                                for c in ["ScreenROCAUC", "CasActivity", "ScreenDoublingTime", "ScreenMADNonessentials", "LibraryAvana", "LibraryKY", "LibraryHumagneCD"]}
    # growth pattern agreement with Model.csv
    gp = Cc[[c for c in C.columns if c.startswith("GrowthPattern")]].idxmax(1).str.replace("GrowthPattern", "")
    conf["growthpattern_agrees_with_Model_csv"] = float((gp.values == Mc.GrowthPattern.values).mean())
    # ModelCondition
    MC = pd.read_csv(DM / "ModelCondition.csv")
    MCc = MC[MC.ModelID.isin(model_ids)]
    src = MCc.groupby("ModelID").DataSource.apply(lambda s: "+".join(sorted(set(s))))
    mcd = {"n_rows": int(len(MC)), "n_rows_crispr_models": int(len(MCc)), "datasource_counts_crispr": vc(MCc.DataSource),
           "models_by_datasource_set": {k: int(v) for k, v in src.value_counts().items()},
           "conditions_per_model": {str(k): int(v) for k, v in MCc.groupby("ModelID").size().value_counts().items()},
           "missing_fraction": {c: round(float(MC[c].isna().mean()), 3) for c in MC.columns},
           "n_with_parent_condition": int(MC.ParentModelConditionID.notna().sum()),
           "CellFormat": vc(MCc.CellFormat), "n_with_prescreen_or_anchor_drug": int(MC.AnchorDrug.notna().sum() + MC.PrescreenTreatmentDrug.notna().sum()),
           "growth_media_distinct": int(MCc.GrowthMedia.nunique())}
    # Sanger-only vs library KY
    sset = src.reindex(model_ids)
    mcd["datasource_vs_library_KY"] = {k: {"n": int((sset == k).sum()), "frac_KY": round(float(Cc.LibraryKY[(sset == k).values].mean()), 2)}
                                       for k in sset.dropna().unique()}
    mcd["lineage_ge15_x_datasource_set_p"] = float(stats.chi2_contingency(pd.crosstab(lin[big15], sset[big15].fillna("none")).values)[1])
    mcd["usefulness"] = ("ModelCondition adds Broad/Sanger source (screen-site batch proxy, coincides with library) and culture medium; "
                         "no drug-treatment conditions among CRISPR models; one row per screen condition (several per model).")
    conf["model_condition"] = mcd
    J["confounders"] = conf

# ------------------------------------------------------------------ 4. essentiality landscape (needed also for graph step)
with step("essentiality"):
    n_models = X.shape[0]
    frac = np.nanmean(X < EFF_THR, 0)  # NaN compare False; adjust by valid count
    valid = (~np.isnan(X)).sum(0)
    frac = (X < EFF_THR).sum(0) / np.maximum(valid, 1)
    mean = np.nanmean(X, 0); sd = np.nanstd(X, 0)
    minv = np.nanmin(X, 0)
    q05 = np.nanpercentile(X, 5, 0); q50 = np.nanpercentile(X, 50, 0)
    skew = np.asarray(stats.skew(X, axis=0, nan_policy="omit"))
    kurt = np.asarray(stats.kurtosis(X, axis=0, nan_policy="omit"))  # excess
    nn = valid.astype(float)
    bc = (skew ** 2 + 1) / (kurt + 3 * (nn - 1) ** 2 / np.maximum((nn - 2) * (nn - 3), 1))
    pan = frac >= PE_FRAC
    pan_alt = mean < -0.5
    bins = [(0, 0.0, "0 (never < -0.5)"), (0.0, 0.01, "(0,1%)"), (0.01, 0.05, "[1,5%)"), (0.05, 0.25, "[5,25%)"),
            (0.25, 0.5, "[25,50%)"), (0.5, 0.9, "[50,90%)"), (0.9, 1.0001, ">=90% (pan-essential)")]
    dist = {}
    for lo, hi, lab in bins:
        m_ = (frac == 0) if lab.startswith("0") else ((frac > lo) & (frac < hi)) if lo == 0.0 else ((frac >= lo) & (frac < hi))
        dist[lab] = int(m_.sum())
    lineage_eta = np.full(X.shape[1], np.nan)
    gl = lin.values; keep = np.isin(gl, lc[lc >= 15].index.to_numpy())
    Xk = X[keep]; gk = gl[keep]
    tot_mean = np.nanmean(Xk, 0); tot_ss = np.nansum((Xk - tot_mean) ** 2, 0)
    bet = np.zeros(X.shape[1])
    for u in np.unique(gk):
        s = Xk[gk == u]; cnt = (~np.isnan(s)).sum(0)
        bet += cnt * (np.nanmean(s, 0) - tot_mean) ** 2
    lineage_eta = bet / np.where(tot_ss > 0, tot_ss, np.nan)
    sel_skew = (~pan) & (frac >= 0.02) & (skew < -1)
    bimod = (~pan) & (frac >= 0.05) & (frac <= 0.95) & (bc > 0.555)
    ess = {"n_genes": int(X.shape[1]),
           "pan_essential_frac_ge90": int(pan.sum()), "pan_essential_mean_lt_-0.5": int(pan_alt.sum()),
           "pan_both_definitions": int((pan & pan_alt).sum()),
           "frac_lt_-0.5_distribution": dist,
           "genes_with_effect_lt_-1_in_ge5_models": int(((X < -1).sum(0) >= 5).sum()),
           "genes_with_effect_lt_-1_in_ge15_models": int(((X < -1).sum(0) >= 15).sum()),
           "selective_skew_lt_-1_and_frac_ge2pct_not_pan": int(sel_skew.sum()),
           "bimodal_BC_gt_0.555_frac5to95_not_pan": int(bimod.sum()),
           "frac_gt_0.9_in_pan_median_mean": float(np.median(mean[pan])),
           "lineage_eta2_ge0.2": int((lineage_eta >= 0.2).sum()), "lineage_eta2_ge0.3": int((lineage_eta >= 0.3).sum()),
           "lineage_eta2_median": float(np.nanmedian(lineage_eta)),
           "lineage_eta2_ge0.2_and_selective": int(((lineage_eta >= 0.2) & ~pan & (frac >= 0.02)).sum()),
           "n_genes_effect_lt_-0.5_in_ge15_models_not_pan": int((((X < EFF_THR).sum(0) >= 15) & ~pan).sum())}
    # per-model strength: mean effect over pan-essential genes
    pe_str = np.nanmean(X[:, pan], 1)
    J["essentiality"] = ess
    conf = J["confounders"]
    conf["pan_essential_mean_effect_per_model"] = {"mean": float(pe_str.mean()), "sd": float(pe_str.std()),
                                                   "eta2_lineage_ge15": round(eta2(pe_str[big15], lin.values[big15]), 3),
                                                   "spearman_with_NNMD": round(float(stats.spearmanr(pe_str, Cc.ScreenNNMD)[0]), 2)}
    # which lineages have extreme strengths
    by_l = pd.Series(pe_str, index=lin.values).groupby(level=0).agg(["mean", "count"])
    conf["pan_essential_mean_effect_by_lineage_ge15"] = {l: round(float(by_l.loc[l, "mean"]), 3) for l in lc[lc >= 15].index}
    # model-wise: correlations of confounders with the dominant PC of the effect matrix (cheap, 1361 models)
    Z = np.where(np.isnan(X), np.nanmean(X, 0), X); Z = Z - Z.mean(0)
    u, s_, vt = np.linalg.svd(Z, full_matrices=False)
    pcs = u[:, :5] * s_[:5]
    ev = (s_ ** 2) / (s_ ** 2).sum()
    conf["effect_matrix_PCs"] = {"var_explained_PC1_5": [round(float(x), 4) for x in ev[:5]]}
    conf["effect_matrix_PCs"]["PC_spearman_with"] = {
        c: [round(float(stats.spearmanr(pcs[:, i], Cc[c], nan_policy="omit")[0]), 2) for i in range(3)]
        for c in ["ScreenNNMD", "ScreenMedianEssentialDepletion", "LibraryKY", "LibraryAvana", "CasActivity"]}
    conf["effect_matrix_PCs"]["PC_eta2_lineage_ge15"] = [round(eta2(pcs[big15, i], lin.values[big15]), 2) for i in range(5)]
    del Z, u, vt

# ------------------------------------------------------------------ 3. genome-wide scores
with step("score_matrix_scan"):
    man = json.load(open(P1 / "manifest.json"))
    G = pd.read_csv(P1 / "genes.tsv", sep="\t")
    n = len(G)
    assert (G.depmap_column.values == cols[G.depmap_index.values]).all(), "genes.tsv depmap_index mismatch"
    lowcov = G.low_coverage.values.astype(bool)
    pe_g = pan[G.depmap_index.values]
    selected_g = (~pe_g) & (frac[G.depmap_index.values] >= 0.02)
    fr_g = frac[G.depmap_index.values]
    cls_g = np.where(fr_g >= PE_FRAC, 0, np.where(fr_g >= 0.05, 1, np.where(fr_g >= 0.01, 2, 3)))  # 0 pan,1 sel(5-90%),2 weak(1-5%),3 ~never(<1%)
    CLS_NAMES = ["pan(>=90%)", "selective(5-90%)", "weak(1-5%)", "never(<1%)"]
    ecls = {c: np.zeros((4, 4), np.int64) for c in CUTS}
    S = np.load(P1 / "zavg_BEF.npy", mmap_mode="r")
    assert S.shape == (n, n)
    deg = {c: np.zeros(n, np.int64) for c in CUTS}
    cnt = {c: 0 for c in CUTS}; cnt_lc = {c: 0 for c in CUTS}; cnt_pe2 = {c: 0 for c in CUTS}
    cnt_pe1 = {c: 0 for c in CUTS}; cnt_nonpe2 = {c: 0 for c in CUTS}; cnt_sel2 = {c: 0 for c in CUTS}
    edges_b = (-8.0, 40.0, 9600)  # 0.005 bins
    hist = np.zeros(edges_b[2], np.int64); n_below = 0; n_above = 0
    for i0 in range(0, n, BLOCK):
        blk = np.array(S[i0:i0 + BLOCK], dtype=np.float32)  # rows i0:i1 x all
        np.nan_to_num(blk, copy=False, nan=-100.0)
        h, _ = np.histogram(blk, bins=edges_b[2], range=edges_b[:2]); hist += h
        n_below += int((blk < edges_b[0]).sum()) - int((blk == -100.0).sum()); n_above += int((blk > edges_b[1]).sum())
        r_lc = lowcov[i0:i0 + BLOCK][:, None]; r_pe = pe_g[i0:i0 + BLOCK][:, None]; r_sel = selected_g[i0:i0 + BLOCK][:, None]
        for c in CUTS:
            m_ = blk >= c
            deg[c][i0:i0 + BLOCK] = m_.sum(1)
            cnt[c] += int(m_.sum())
            cnt_lc[c] += int((m_ & (r_lc | lowcov[None, :])).sum())
            cnt_pe2[c] += int((m_ & r_pe & pe_g[None, :]).sum())
            cnt_pe1[c] += int((m_ & (r_pe | pe_g[None, :])).sum())
            cnt_sel2[c] += int((m_ & r_sel & selected_g[None, :]).sum())
            rr, cc_ = np.nonzero(m_)
            ecls[c] += np.bincount(cls_g[i0:i0 + BLOCK][rr] * 4 + cls_g[cc_], minlength=16).reshape(4, 4)
        print(f"   block {i0}/{n}", end="\r", flush=True)
    # quantiles from histogram (full matrix, includes mirrored duplicates; symmetric so same distribution)
    tot = hist.sum(); cs = np.cumsum(hist)
    qs = [0.0001, 0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999, 0.9999, 0.99999, 0.999999]
    centers = edges_b[0] + (np.arange(edges_b[2]) + 0.5) * (edges_b[1] - edges_b[0]) / edges_b[2]
    quant = {str(q): float(centers[np.searchsorted(cs, q * tot)]) for q in qs}
    mq = man["scores"]["zavg_BEF_over_all_pairs"]["quantiles"]
    sm = {"n_offdiag_values_scanned": int(tot), "quantiles_hist(0.005 bins)": quant,
          "manifest_quantiles_for_comparison": {k: v for k, v in mq.items() if k in {str(q) for q in qs}},
          "mean_sd_manifest": [man["scores"]["zavg_BEF_over_all_pairs"]["mean"], man["scores"]["zavg_BEF_over_all_pairs"]["sd"]],
          "n_values_outside_[-8,40]": int(n_below + n_above)}
    sm["pairs_ge"] = {str(c): cnt[c] // 2 for c in CUTS}
    sm["pairs_ge_manifest_check"] = {k: v for k, v in man["scores"]["zavg_BEF_over_all_pairs"]["count_ge"].items()}
    sm["frac_edges_low_coverage_any"] = {str(c): round(cnt_lc[c] / max(cnt[c], 1), 4) for c in CUTS}
    sm["baseline_frac_pairs_low_coverage_any"] = round(1 - (1 - lowcov.mean()) ** 2, 4)
    sm["frac_edges_both_pan_essential"] = {str(c): round(cnt_pe2[c] / max(cnt[c], 1), 4) for c in CUTS}
    sm["frac_edges_any_pan_essential"] = {str(c): round(cnt_pe1[c] / max(cnt[c], 1), 4) for c in CUTS}
    sm["frac_edges_both_nonpan_with_frac_ge2pct"] = {str(c): round(cnt_sel2[c] / max(cnt[c], 1), 4) for c in CUTS}
    sm["n_edges_both_nonpan_with_frac_ge2pct"] = {str(c): cnt_sel2[c] // 2 for c in CUTS}
    sm["n_genes_pan_essential_in_universe"] = int(pe_g.sum())
    sm["n_genes_low_coverage"] = int(lowcov.sum())
    dg = {}
    for c in CUTS:
        d = deg[c]
        dg[str(c)] = {"genes_ge1": int((d >= 1).sum()), "genes_ge3": int((d >= 3).sum()), "genes_ge10": int((d >= 10).sum()),
                      "genes_ge50": int((d >= 50).sum()),
                      "median_deg_among_ge1": float(np.median(d[d >= 1])) if (d >= 1).any() else 0,
                      "p90": float(np.percentile(d[d >= 1], 90)) if (d >= 1).any() else 0,
                      "p99": float(np.percentile(d[d >= 1], 99)) if (d >= 1).any() else 0, "max": int(d.max()),
                      "frac_lowcov_genes_ge1": round(float((d[lowcov] >= 1).mean()), 3),
                      "frac_pan_genes_ge1": round(float((d[pe_g] >= 1).mean()), 3),
                      "frac_nonpan_genes_ge1": round(float((d[~pe_g] >= 1).mean()), 3),
                      "top_hubs": [(G.symbol.values[j], int(d[j])) for j in np.argsort(-d)[:8]]}
    sm["degree"] = dg
    sm["gene_classes_by_frac_lt_-0.5"] = {CLS_NAMES[i]: int((cls_g == i).sum()) for i in range(4)}
    sm["edge_class_matrix_row_col_order"] = CLS_NAMES
    sm["edges_by_class_pair"] = {}
    sm["frac_genes_deg_ge1_by_class"] = {}
    sm["median_deg_by_class"] = {}
    for c in CUTS:
        E_ = ecls[c]; E_ = (E_ + E_.T)  # symmetric double count -> keep as ordered counts / 2 for unordered
        und = {f"{CLS_NAMES[i]}-{CLS_NAMES[j]}": int((E_[i, j] if i != j else E_[i, i]) // 2) for i in range(4) for j in range(i, 4)}
        # ecls counts both orientations already (full matrix); unordered pair (i,j) i!=j = ecls[i,j]; i==j = ecls[i,i]/2
        und = {f"{CLS_NAMES[i]}-{CLS_NAMES[j]}": int(ecls[c][i, j] if i != j else ecls[c][i, i] // 2) for i in range(4) for j in range(i, 4)}
        sm["edges_by_class_pair"][str(c)] = und
        sm["frac_genes_deg_ge1_by_class"][str(c)] = {CLS_NAMES[i]: round(float((deg[c][cls_g == i] >= 1).mean()), 3) for i in range(4)}
        sm["median_deg_by_class"][str(c)] = {CLS_NAMES[i]: float(np.median(deg[c][cls_g == i])) for i in range(4)}
    d5 = deg[5]
    sm["hubs_deg_ge50_at_5_class_counts"] = {CLS_NAMES[i]: int(((d5 >= 50) & (cls_g == i)).sum()) for i in range(4)}
    sm["spearman_deg5_vs_gene_sd"] = round(float(stats.spearmanr(d5, sd[G.depmap_index.values])[0]), 3)
    sm["spearman_deg5_vs_frac"] = round(float(stats.spearmanr(d5, fr_g)[0]), 3)
    sm["degree_at_5_hist"] = {lab: int(((deg[5] >= lo) & (deg[5] < hi)).sum())
                              for lab, lo, hi in [("0", 0, 1), ("1", 1, 2), ("2-3", 2, 4), ("4-9", 4, 10), ("10-24", 10, 25),
                                                  ("25-99", 25, 100), (">=100", 100, 10 ** 9)]}
    J["score_matrix"] = sm

# ------------------------------------------------------------------ graph / components / leiden
with step("graph_components"):
    T = pd.read_parquet(P1 / "top_pairs.parquet", columns=["entrez_a", "entrez_b", "zavg_BEF", "low_coverage_any"])
    e2i = dict(zip(G.entrez.values, np.arange(n)))
    ia = T.entrez_a.map(e2i).values; ib = T.entrez_b.map(e2i).values
    assert not (pd.isna(ia).any() or pd.isna(ib).any())
    ia = ia.astype(np.int64); ib = ib.astype(np.int64); zz = T.zavg_BEF.values.astype(np.float64)
    J["score_matrix"]["top_pairs_rows"] = int(len(T)); J["score_matrix"]["top_pairs_matches_scan_ge3"] = bool(len(T) == cnt[3] // 2)
    comp = {}
    for c in [3, 4, 5, 6]:
        k = zz >= c
        a, b = ia[k], ib[k]
        A = coo_matrix((np.ones(len(a)), (a, b)), shape=(n, n))
        nc, lab = connected_components(A, directed=False)
        sizes = np.bincount(lab); deg_ = np.bincount(np.r_[a, b], minlength=n)
        sizes_nz = np.sort(sizes[np.unique(lab[deg_ > 0])])[::-1]
        comp[str(c)] = {"n_edges": int(k.sum()), "n_nodes_nonisolated": int((deg_ > 0).sum()), "n_components_nonisolated": int(len(sizes_nz)),
                        "largest": [int(x) for x in sizes_nz[:5]], "n_components_size_ge3": int((sizes_nz >= 3).sum()),
                        "n_components_size_ge10": int((sizes_nz >= 10).sum()), "n_components_size_2": int((sizes_nz == 2).sum()),
                        "frac_nodes_in_largest": round(float(sizes_nz[0] / (deg_ > 0).sum()), 3)}
    J["components"] = comp

with step("leiden_z5"):
    import igraph as ig, leidenalg as la
    k = zz >= 5
    a, b, w = ia[k], ib[k], zz[k]
    nodes = np.unique(np.r_[a, b]); remap = {int(x): i for i, x in enumerate(nodes)}
    g = ig.Graph(n=len(nodes), edges=[(remap[int(x)], remap[int(y)]) for x, y in zip(a, b)])
    g.es["w"] = w.tolist()
    pe_nodes = pe_g[nodes]
    res = {"graph": {"n_nodes": int(g.vcount()), "n_edges": int(g.ecount())}, "runs": []}
    runs = [("RBConfiguration(weighted z)", la.RBConfigurationVertexPartition, r, "w") for r in (0.25, 0.5, 1.0, 2.0, 5.0)] + \
           [("CPM(unweighted)", la.CPMVertexPartition, r, None) for r in (0.05, 0.1, 0.25, 0.5)]
    for name, cls, r, wt in runs:
        t0 = time.perf_counter()
        part = la.find_partition(g, cls, weights=wt, resolution_parameter=r, seed=SEED, n_iterations=2)
        dt = time.perf_counter() - t0
        mem = np.array(part.membership); sz = np.bincount(mem); sz_s = np.sort(sz)[::-1]
        pe_frac_c = np.bincount(mem, weights=pe_nodes.astype(float)) / sz
        big_ = sz >= 3
        res["runs"].append({"method": name, "resolution": r, "runtime_s": round(dt, 2), "n_clusters": int(len(sz)),
                            "singletons": int((sz == 1).sum()), "pairs": int((sz == 2).sum()), "size_ge3": int(big_.sum()),
                            "size_ge5": int((sz >= 5).sum()), "size_ge10": int((sz >= 10).sum()), "size_ge50": int((sz >= 50).sum()),
                            "size_ge100": int((sz >= 100).sum()), "median_size_ge3": float(np.median(sz[big_])) if big_.any() else 0,
                            "largest5": [int(x) for x in sz_s[:5]], "frac_nodes_in_clusters_ge3": round(float(sz[big_].sum() / g.vcount()), 3),
                            "frac_clusters_ge3_mostly_pan_ess(>=80%)": round(float((pe_frac_c[big_] >= 0.8).mean()), 3) if big_.any() else 0,
                            "modularity_RB1": round(float(g.modularity(mem, weights="w")), 3)})
        print("   ", res["runs"][-1]["method"], r, res["runs"][-1]["n_clusters"], res["runs"][-1]["runtime_s"], flush=True)
    J["leiden_z5"] = res

J["timings_s"] = TIM; J["peak_rss_gb_cumulative_after_step"] = MEM
J["peak_rss_gb_final"] = round(rss_gb(), 2)
J["total_runtime_s"] = round(sum(TIM.values()), 1)
J["env"] = {"python": sys.version.split()[0], "numpy": np.__version__, "pandas": pd.__version__}
J["inputs_not_used"] = "datasets/huMAP/* never opened; no curated complex DB"
J["part1_manifest_summary"] = {"created": man["created"], "n_protein_genes": man["counts"]["n_protein_genes"],
                               "n_low_coverage": man["counts"]["n_protein_genes_low_coverage"], "code_commit": man["code"]["commit"]}


def clean(o):
    if isinstance(o, dict): return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [clean(v) for v in o]
    if isinstance(o, (np.floating, float)): return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.integer): return int(o)
    if isinstance(o, np.bool_): return bool(o)
    return o


J = clean(J)
json.dump(J, open(OUT / "part2_profile.json", "w"), indent=1)

# ------------------------------------------------------------------ markdown
def tbl(header, rows):
    s = "| " + " | ".join(header) + " |\n|" + "|".join(["---"] * len(header)) + "|\n"
    for r in rows: s += "| " + " | ".join(str(x) for x in r) + " |\n"
    return s


m = J["models"]; mo = J["molecular_annotation"]; cf = J["confounders"]; sm = J["score_matrix"]; es = J["essentiality"]
L = []
L.append("# Part 2 data profile (descriptive; no labels)\n")
L.append(f"Script `scripts/05_profile_part2.py` (seed {SEED}); total {J['total_runtime_s']} s, peak RSS {J['peak_rss_gb_final']} GB. "
         "No huMAP files or curated complex DBs read.\n")
L.append("## 1. Models\n")
L.append(f"CRISPRGeneEffect: **{J['crispr_matrix']['n_models']} models x {J['crispr_matrix']['n_genes']} genes** "
         f"(NaN {100*J['crispr_matrix']['nan_fraction']:.1f}%); all join to Model.csv ({m['n_model_csv_without_crispr']} Model.csv rows have no screen). "
         f"{m['n_lineages']} lineages, {m['n_primary_disease']} primary diseases ({m['n_primary_disease_ge10']} with >=10), "
         f"{m['n_subtype']} subtypes ({m['n_subtype_ge10']} with >=10). ScreenGeneEffect has {J['screen_gene_effect']['n_screens']} screens "
         f"(rows `SC-id`, {J['screen_gene_effect']['n_unique_models_from_SC_id']} models; not linkable to ModelCondition).\n")
L.append(tbl(["Lineage (n>=15)", "n", "Lineage (n>=15)", "n"],
             [[a[0], a[1], b[0] if b else "", b[1] if b else ""] for a, b in
              zip(list(m["lineages_ge15"].items())[::2], list(m["lineages_ge15"].items())[1::2] + [None])]))
L.append(f"\nLineages >=15 models: {len(m['lineages_ge15'])} ({m['n_models_in_lineages_ge15']} models); >=30: {len(m['lineages_ge30'])} "
         f"({m['n_models_in_lineages_ge30']} models). Primary diseases >=30: {len(m['primary_disease_ge30'])}; >=15: {len(m['primary_disease_ge15'])}; "
         f"subtypes >=15: {m['n_subtype_ge15']}.\n")
an = J["model_annotations"]
L.append(tbl(["Field", "Distribution (CRISPR models)"],
             [["Sex", an["Sex"]], ["PrimaryOrMetastasis", an["PrimaryOrMetastasis"]], ["AgeCategory", an["AgeCategory"]],
              ["PediatricModelType", an["PediatricModelType"]], ["GrowthPattern", an["GrowthPattern"]], ["ModelType", an["ModelType"]],
              ["EngineeredModel (flag to exclude)", an["EngineeredModel"]],
              ["Age", f"n={an['Age']['n_nonmissing']}, median {an['Age'].get('50%')}, range {an['Age'].get('min')}-{an['Age'].get('max')}"]]))
L.append("\n**Molecular annotation in Model.csv.** No mutation/CN/expression columns. Only free text: "
         f"ModelSubtypeFeatures n={mo['n_ModelSubtypeFeatures']} ({mo['n_distinct_ModelSubtypeFeatures']} distinct), PatientSubtypeFeatures "
         f"n={mo['n_PatientSubtypeFeatures']}; any feature: {mo['n_any_feature']} models ({100*mo['frac_any_feature']:.0f}%). Keyword hits (feature fields):\n")
kc = mo["keyword_counts_in_feature_fields"]
keys = [k for k in kc if kc[k]["n_models"] >= 3]
L.append(tbl(["Keyword", "n models", "top lineages"], [[k, kc[k]["n_models"], ", ".join(f"{a}:{b}" for a, b in kc[k]["lineages"].items())] for k in keys]))
L.append(f"\nSubtype-name proxies: {mo['keyword_counts_in_subtype_names']}. "
         "Conclusion: only a few binary contexts reach >=10 models (MSI, HER2+/TNBC/ER, MYC-amp, BCR::ABL1, EWS fusions, IDH-mutant/GBM names, MGMT), "
         "mostly lineage-nested; KRAS/BRAF/TP53 status etc. cannot be derived without an external mutation file (candidate DATA_REQUEST).\n")
L.append("## 2. Confounders and ModelCondition\n")
L.append(f"CRISPRConfounders: {cf['n_rows_total']} rows x {cf['n_columns']} cols; covers all CRISPR models ({cf['missing_fraction_crispr_models']['ScreenNNMD']*100:.0f}% missing for "
         f"quality cols among CRISPR models; CasActivity {100*cf['missing_fraction_crispr_models']['CasActivity']:.0f}%, DoublingTime {100*cf['missing_fraction_crispr_models']['ScreenDoublingTime']:.0f}% missing). "
         f"Libraries: {cf['library_counts']}; screen type {cf['screentype_counts']}. Lineage(>=15) x library Cramer's V = {cf['lineage_ge15_x_library']['cramers_v']:.2f} (p={cf['lineage_ge15_x_library']['p']:.1e}).\n")
e = cf["eta2_lineage_ge15"]
L.append(tbl(["Covariate", "eta^2 vs lineage"], [[k, v] for k, v in e.items()] +
             [["pan-essential mean effect (derived)", cf["pan_essential_mean_effect_per_model"]["eta2_lineage_ge15"]]]))
pcs = cf["effect_matrix_PCs"]
L.append(f"\nEffect-matrix PCA: var explained PC1-5 = {pcs['var_explained_PC1_5']}; lineage eta^2 of PC1-5 = {pcs['PC_eta2_lineage_ge15']}; "
         f"PC1-3 Spearman with NNMD {pcs['PC_spearman_with']['ScreenNNMD']}, with KY library {pcs['PC_spearman_with']['LibraryKY']}. "
         f"Pan-essential strength per model vs NNMD rho={cf['pan_essential_mean_effect_per_model']['spearman_with_NNMD']}. "
         f"ModelCondition: {cf['model_condition']['n_rows']} rows, Broad/Sanger/other = {cf['model_condition']['datasource_counts_crispr']}; "
         f"models by source set {cf['model_condition']['models_by_datasource_set']}; Sanger~library KY {cf['model_condition']['datasource_vs_library_KY']}. "
         "Useful only as site/medium batch labels.\n")
L.append("## 3. Genome-wide score matrix (zavg_BEF, n=%d genes, %d low-coverage)\n" % (J["part1_manifest_summary"]["n_protein_genes"], sm["n_genes_low_coverage"]))
q = sm["quantiles_hist(0.005 bins)"]
L.append(tbl(["Quantile"] + list(q.keys()), [["z"] + [f"{v:.2f}" for v in q.values()]]))
rows = []
for c in CUTS:
    d = sm["degree"][str(c)]
    rows.append([f">={c}", f"{sm['pairs_ge'][str(c)]:,}", d["genes_ge1"], d["genes_ge3"], d["genes_ge10"], d["median_deg_among_ge1"], d["max"],
                 sm["frac_edges_low_coverage_any"][str(c)], sm["frac_edges_both_pan_essential"][str(c)], sm["frac_edges_both_nonpan_with_frac_ge2pct"][str(c)]])
L.append("\n" + tbl(["z", "pairs", "genes deg>=1", ">=3", ">=10", "median deg", "max deg", "edge frac low-cov", "edge frac both pan-ess", "edge frac both non-pan"], rows))
L.append(f"\nBaseline low-coverage-any pair fraction {sm['baseline_frac_pairs_low_coverage_any']}. Top hubs z>=5: {sm['degree']['5']['top_hubs']}. "
         f"Degree histogram at z>=5: {sm['degree_at_5_hist']}.\n")
ec = sm["edges_by_class_pair"]
L.append("\nGene classes by frac(effect<-0.5): " + str(sm["gene_classes_by_frac_lt_-0.5"]) + ". Edges by class pair:\n")
pairs_ = list(ec["5"].keys())
L.append(tbl(["z"] + pairs_, [[c] + [ec[str(c)][k] for k in pairs_] for c in CUTS]))
L.append(f"\nMedian degree by class at z>=5: {sm['median_deg_by_class']['5']}; hubs (deg>=50 at z>=5) by class: {sm['hubs_deg_ge50_at_5_class_counts']}; "
         f"Spearman(deg, gene SD)={sm['spearman_deg5_vs_gene_sd']}, (deg, frac<-0.5)={sm['spearman_deg5_vs_frac']}.\n")
cp = J["components"]
L.append(tbl(["z", "edges", "nodes", "components", "largest 3 sizes", "components >=3", "frac nodes in giant"],
             [[c, v["n_edges"], v["n_nodes_nonisolated"], v["n_components_nonisolated"], v["largest"][:3], v["n_components_size_ge3"], v["frac_nodes_in_largest"]] for c, v in cp.items()]))
lr = J["leiden_z5"]
L.append(f"\nLeiden on z>=5 graph ({lr['graph']['n_nodes']} nodes, {lr['graph']['n_edges']} edges; seed 0, 2 iterations):\n")
L.append(tbl(["Method", "res", "s", "clusters", "single", "size>=3", ">=10", ">=50", "median(>=3)", "largest", "pan-ess-dominated (>=80%)"],
             [[r["method"].split("(")[0], r["resolution"], r["runtime_s"], r["n_clusters"], r["singletons"], r["size_ge3"], r["size_ge10"], r["size_ge50"],
               r["median_size_ge3"], r["largest5"][:2], r["frac_clusters_ge3_mostly_pan_ess(>=80%)"]] for r in lr["runs"]]))
L.append("\n## 4. Essentiality landscape\n")
L.append(tbl(["Quantity", "Genes"],
             [["pan-essential (frac<-0.5 >= 90%)", es["pan_essential_frac_ge90"]], ["pan-essential (mean<-0.5)", es["pan_essential_mean_lt_-0.5"]],
              ["effect<-1 in >=15 models", es["genes_with_effect_lt_-1_in_ge15_models"]],
              ["non-pan, effect<-0.5 in >=15 models", es["n_genes_effect_lt_-0.5_in_ge15_models_not_pan"]],
              ["selective: skew<-1, frac>=2%, non-pan", es["selective_skew_lt_-1_and_frac_ge2pct_not_pan"]],
              ["bimodal (BC>0.555), frac 5-95%, non-pan", es["bimodal_BC_gt_0.555_frac5to95_not_pan"]],
              ["lineage eta^2 >=0.2 / >=0.3", f"{es['lineage_eta2_ge0.2']} / {es['lineage_eta2_ge0.3']}"]] +
             [[f"frac<-0.5 in {k}", v] for k, v in es["frac_lt_-0.5_distribution"].items()]))
L.append("\n## 5. Timings (s) and peak RSS (GB, cumulative max)\n")
L.append(tbl(["Step", "s", "peak RSS GB"], [[k, v, J["peak_rss_gb_cumulative_after_step"][k]] for k, v in TIM.items()]))
open(OUT / "part2_profile.md", "w").write("\n".join(L))
print("done", J["total_runtime_s"], "s")
