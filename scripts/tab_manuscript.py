"""Manuscript tables, built only from committed result files (no hardcoded numbers).

 Table 1  test metrics (results/test_eval/test_results.json)
 Table 2  map summary (results/part2/construction/{summary,tiers,bridge,variability})
 Table S1 pre-registered controls (results/part2/controls/controls_summary.tsv, kill_rule.json)

Writes Markdown and TSV to manuscript/tables/.
Usage: python scripts/tab_manuscript.py
"""
from __future__ import annotations

import pandas as pd

import fig_style as S


def ci(st: dict, nd: int = 3) -> str:
    return S.fmt_ci(st["point"], st["ci_low"], st["ci_high"], nd)


def md_table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for _, r in df.iterrows():
        out.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join(out)


def write(stem: str, title: str, notes: list[str], df: pd.DataFrame, extra: str = "") -> None:
    S.TAB_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(S.TAB_DIR / f"{stem}.tsv", sep="\t", index=False)
    txt = f"# {title}\n\n{md_table(df)}\n\n{extra}" + "".join(f"- {n}\n" for n in notes)
    (S.TAB_DIR / f"{stem}.md").write_text(txt)


def pfmt(p: float, bound: float) -> str:
    """Bootstrap p; 0 means no replicate on the other side of 0, shown as its bound (revision round 1, minor 17)."""
    return f"< {bound:.3g}" if p == 0 else f"{p:.3g}"


def table1() -> None:
    d = S.load_json(S.TEST_RESULTS)
    rows = []
    for sname, skey in (("ITT", "itt"), ("S1", "mapped_only")):
        r = d["results"][skey]
        st = r["bootstrap"]["cluster"]["stat"]
        for m in S.METHODS:
            mm = r["metrics"][m]
            iso = mm["calibration"]["isotonic"]
            rows.append({
                "set": sname, "method": S.LABEL[m],
                "AP [95% CI]": ci(st[f"{m}|ap"]),
                "AUROC [95% CI]": ci(st[f"{m}|auroc"]),
                "macro per-cluster AP [95% CI]": ci(st[f"{m}|macro_cluster_ap"]),
                "AP w/o largest comp. [95% CI]": ci(st[f"{m}|ap_no_largest_component"]),
                "AP w/o 10 largest": f"{mm['ap_no_top10_clusters']:.3f}",
                "P, rule (b) [95% CI]": ci(st[f"{m}|precision_rule_b"]),
                "R, rule (b) [95% CI]": ci(st[f"{m}|recall_rule_b"]),
                "F1, rule (b) [95% CI]": ci(st[f"{m}|f1_rule_b"]),
                "Brier (isotonic) [95% CI]": ci(st[f"{m}|brier_isotonic"], 5),
                "Brier skill": f"{iso['bss_vs_fit_prevalence']:.3f}",
                "cal.-in-the-large": f"{iso['calibration_in_the_large']:.3f}",
                "cal. slope": f"{iso['calibration_slope']:.3f}",
            })
    df = pd.DataFrame(rows)
    # Delta block
    drows = []
    for sname, skey in (("ITT", "itt"), ("S1", "mapped_only")):
        r = d["results"][skey]
        cs = r["ci_statistic"]
        nb = r["bootstrap"]["cluster"]["n_boot"]
        for cand in ["zavg_BEF", "B"]:
            for b in S.BASELINES:
                v = r["bootstrap"]["cluster"]["diff"][f"{cand}|ap - {b}|ap"]
                lv = r["locked_verdict"]["pairs"][f"{cand}-{b}"]
                pc = cs["pairs"][f"{cand} vs {b}"]
                holm = cs["holm_adjusted_p_primary"].get(f"{cand} vs {b}") if cand == "zavg_BEF" else None
                drows.append({
                    "set": sname, "comparison": f"{S.SHORT[cand]} − {S.SHORT[b]}",
                    "ΔAP [95% CI]": S.fmt_ci(v["diff"], v["ci_low"], v["ci_high"], 3, True),
                    "relative gain": f"{v['rel_gain']:+.3f}",
                    "macro-AP ratio [95% CI]": S.fmt_ci(pc["non_degradation_ratio_ci"]["macro_cluster_ap"]["ratio"],
                                                        pc["non_degradation_ratio_ci"]["macro_cluster_ap"]["ci_low"],
                                                        pc["non_degradation_ratio_ci"]["macro_cluster_ap"]["ci_high"]),
                    "p (two-sided)": pfmt(v["p_two_sided"], 2 / nb),
                    "Holm p": "" if holm is None else pfmt(holm, len(S.BASELINES) * 2 / nb),
                    "locked rule passes": "yes" if lv["passes"] else "no",
                    "CI excludes 0": "yes" if pc["ci_excludes_0"] else "no",
                    "borderline": "yes" if pc["borderline"] else "no",
                })
    dd = pd.DataFrame(drows)
    dd.to_csv(S.TAB_DIR / "table1b_test_deltas.tsv", sep="\t", index=False)
    itt = d["results"]["itt"]
    notes = [
        f"Source: {S.TEST_RESULTS} (single confirmatory look; started commit {d['started_commit'][:7]}, "
        f"protocol {d['protocol']['commit'][:7]}).",
        f"ITT: {itt['n_pos']:,} positive / {itt['n_neg']:,} negative pairs; unmapped pairs kept at the floor score. "
        f"S1: pairs mapped under every method ({d['results']['mapped_only']['n_pos']:,} / "
        f"{d['results']['mapped_only']['n_neg']:,}).",
        f"CIs: paired cluster bootstrap over {itt['bootstrap']['cluster']['n_units']} test clusters, "
        f"B = {itt['bootstrap']['cluster']['n_boot']:,}, seed {itt['bootstrap']['cluster']['seed']}, 95% percentile.",
        "Rule (b): threshold maximising F1 on all hu.MAP train pairs (fitted before the freeze). "
        "Calibration: isotonic regression fitted on train pairs, applied unchanged; Brier skill against the train prevalence.",
        "The locked verdict (Part 1 'clearly beats' rule) is computed on S1; ITT is the primary metric. "
        "p = 0 means no bootstrap replicate on the other side of 0 (p < 1/2,000).",
        "PC-removed Pearson is a post hoc reference, not a candidate; B is the pre-declared secondary candidate.",
    ]
    extra = "## Table 1b: paired differences\n\n" + md_table(dd) + "\n\n## Notes\n\n"
    write("table1_test_metrics", "Table 1: test-set performance of the locked method and baselines", notes, df, extra)


def table2() -> None:
    summ = S.load_json(S.SUMMARY)
    tiers = S.load_json(S.TIERS)
    var = S.load_json(S.VARIAB)
    br = S.load_json(S.BRIDGE)
    rows = []

    def add(item, z, g):
        rows.append({"quantity": item, "zavg_BEF map": z, "GLS map": g})

    Z, G = summ["zavg"], summ["gls"]
    add("level-1 complexes", Z["n_level1"], G["n_level1"])
    add("genes in level-1 complexes (of S)", f"{Z['genes_level1']:,} / {Z['genes_S']:,}", f"{G['genes_level1']:,} / {G['genes_S']:,}")
    add("sub-modules (parents with any)", f"{Z['n_level2']} ({Z['n_parents_with_submodules']})",
        f"{G['n_level2']} ({G['n_parents_with_submodules']})")
    for key, lab in (("size_level1", "level-1 size, median (IQR; max)"),):
        add(lab, f"{Z[key]['median']:g} ({Z[key]['q25']:g}–{Z[key]['q75']:g}; {Z[key]['max']:g})",
            f"{G[key]['median']:g} ({G[key]['q25']:g}–{G[key]['q75']:g}; {G[key]['max']:g})")
    add("level-1 size bins 3–4 / 5–8 / 9–16 / 17–60", " / ".join(str(Z["size_bins_level1"][b]) for b in ["3-4", "5-8", "9-16", "17-60"]),
        " / ".join(str(G["size_bins_level1"][b]) for b in ["3-4", "5-8", "9-16", "17-60"]))
    for c, lab in (("pan", "pan-essential"), ("selective", "selective"), ("other", "other")):
        add(f"coverage of {lab} S genes",
            f"{Z['coverage_by_class'][c]['clustered']:,}/{Z['coverage_by_class'][c]['in_S']:,} ({Z['coverage_by_class'][c]['frac']:.1%})",
            f"{G['coverage_by_class'][c]['clustered']:,}/{G['coverage_by_class'][c]['in_S']:,} ({G['coverage_by_class'][c]['frac']:.1%})")
    add("locus-flagged, level 1 / level 2", f"{Z['locus_rate_level1']:.1%} / {Z['locus_rate_level2']:.1%}",
        f"{G['locus_rate_level1']:.1%} / {G['locus_rate_level2']:.1%}")
    for lvl in ("level1", "level2"):
        tz, tg = tiers["tiers"]["zavg"][lvl], tiers["tiers"]["gls"][lvl]
        add(f"tiers {lvl.replace('level', 'level ')}: T1 / T2 / T3-neg / T3-nc",
            f"{tz['T1']} / {tz['T2']} / {tz['T3-neg']} / {tz['T3-nc']}",
            f"{tg['T1']} / {tg['T2']} / {tg['T3-neg']} / {tg['T3-nc']}")
        add(f"expected false T1 ({lvl.replace('level', 'level ')})", f"{tz['expected_false_T1']:.1f}", f"{tg['expected_false_T1']:.1f}")
    nz, ng = Z["tiers_nonlocus_level1"], G["tiers_nonlocus_level1"]
    add("non-locus level 1: T1 / T2 / T3-neg / T3-nc", f"{nz['T1']} / {nz['T2']} / {nz['T3-neg']} / {nz['T3-nc']}",
        f"{ng['T1']} / {ng['T2']} / {ng['T3-neg']} / {ng['T3-nc']}")
    vz, vg = var["zavg_BEF"], var["gls"]
    az, ag = vz["size_bins"]["all"], vg["size_bins"]["all"]
    add("best-match Jaccard in bootstrap maps, pooled median (fraction ≥ 0.5)",
        f"{az['observed_pooled']['median']:.2f} ({az['observed_pooled']['frac_ge_0_5']:.2f})",
        f"{ag['observed_pooled']['median']:.2f} ({ag['observed_pooled']['frac_ge_0_5']:.2f})")
    add("random size-matched sets: median (P95)", f"{az['random_S_pooled']['median']:.2f} ({az['random_S_p95']:.2f})",
        f"{ag['random_S_pooled']['median']:.2f} ({ag['random_S_p95']:.2f})")
    add("complexes with median Jaccard > random P95", f"{az['frac_complexes_median_above_random_S_p95']:.1%}",
        f"{ag['frac_complexes_median_above_random_S_p95']:.1%}")
    add("pairwise ARI of bootstrap maps, median (IQR)",
        f"{vz['pairwise']['ari_both']['median']:.3f} ({vz['pairwise']['ari_both']['q25']:.3f}–{vz['pairwise']['ari_both']['q75']:.3f})",
        f"{vg['pairwise']['ari_both']['median']:.3f} ({vg['pairwise']['ari_both']['q25']:.3f}–{vg['pairwise']['ari_both']['q75']:.3f})")
    p = br["primary"]
    add("bridge: size-standardised BioPlex precision (matched coverage)", f"{p['zavg']:.3f}", f"{p['gls']:.3f}")
    df = pd.DataFrame(rows)
    pr = br["primary"]
    notes = [
        f"Sources: {S.SUMMARY}, {S.TIERS}, {S.VARIAB}, {S.BRIDGE}. Maps frozen at commit "
        f"{br['map_commits']['zavg_map'][:7]} (FZ-C).",
        f"Bridge primary difference zavg − GLS: {S.fmt_ci(pr['diff'], pr['ci95'][0], pr['ci95'][1], 3, True)} "
        f"(protein-node bootstrap, B = {pr['B']}, seed {pr['seed']}); decision: {br['decision']}.",
        "Coverage classes by A-essentiality: pan-essential = effect < −0.5 in ≥ 90% of A models, selective 5–90%.",
        "Tiers use BioPlex 3.0 and CCLE proteomics (A models only), read only after the maps were committed. "
        "T2 is not locus-aware (copy-number dosage makes neighbouring proteins co-abundant).",
        "Variability: 20 single-subsample maps (PatientID-blocked 80% subsamples of A). The frozen map is built from the "
        "same subsamples, so agreement with it is not independent; the random baseline is a floor, not a null.",
    ]
    write("table2_map_summary", "Table 2: complex maps on the construction set A", notes, df, "## Notes\n\n")


def table_s1() -> None:
    c = pd.read_csv(S.ROOT / S.CONTROLS, sep="\t")
    kill = S.load_json(S.KILL)
    rows = []
    for _, r in c.iterrows():
        rows.append({
            "control": r["design_id"], "direction": r["direction"], "counted": r["counted"],
            "context": r["context"], "unit": r["unit"], "level": r["level"],
            "n (context)": "" if pd.isna(r["n"]) else f"{int(r['n'])} ({int(r['n_context'])})",
            "β (Chronos)": "" if pd.isna(r["beta"]) else f"{r['beta']:+.3f}",
            "t": "" if pd.isna(r["t"]) else f"{r['t']:+.2f}",
            "p (one-sided)": "" if pd.isna(r["p"]) else f"{r['p']:.3g}",
            "alternative": "" if pd.isna(r["alternative"]) else r["alternative"],
            "evaluable for kill rule (C rows)": r["evaluable_kill"] if r["level"] == "C" else "",
            "recovered (C rows)": r["recovered"] if r["level"] == "C" else "",
            "p < 0.05 (G rows)": r["gene_level_p_lt_005"] if r["level"] == "G" else "",
            "note": "" if pd.isna(r["note"]) else r["note"],
        })
    df = pd.DataFrame(rows)
    order = {"PC%d" % i: i for i in range(1, 23)}
    df["_o"] = df["control"].map(lambda x: (order.get(x.rstrip("ab"), 99), x))
    df = df.sort_values(["_o", "level"]).drop(columns="_o")
    krows = [{"direction": k, "controls": ", ".join(v["controls"]), "evaluable": v["evaluable"],
              "recovered": v["recovered"], "verdict": v["verdict"]} for k, v in kill.items()]
    kd = pd.DataFrame(krows)
    kd.to_csv(S.TAB_DIR / "tableS1b_kill_rule.tsv", sep="\t", index=False)
    notes = [
        f"Sources: {S.CONTROLS}, {S.KILL} (single run, results commit 4b5f800 after intent e82b657).",
        "Model: OLS of the unit score on the context plus registered covariates, HC3 standard error, "
        "Satterthwaite df; one-sided test in the registered direction at α = 0.05; A models only (n ≤ 1,007).",
        "β < 0 with alternative 'less' means more dependent (more negative Chronos effect) in the context. "
        "PC2 (context TP53-altered, expected dependency in TP53-WT) and PC5 (BRAF vs RAS hotspot, expected SHOC2 "
        "dependency in RAS) are coded so that the registered alternative is 'greater'.",
        "Level C = frozen complex unit (counted for the kill rule when evaluable); G = gene-level readout (reported, never counted).",
        "Kill rule (DESIGN §2.6): a direction is dropped if fewer than 50% of its evaluable complex-level controls are recovered.",
        "PC6 and PC12 reduce to single paralog genes under the registered paralog-component readout; only PC15 is a "
        "multi-gene H1 readout (independent review, results/part2/verify/controls_review.md).",
    ]
    extra = "## Kill-rule tally\n\n" + md_table(kd) + "\n\n## Notes\n\n"
    write("tableS1_controls", "Supplementary Table S1: pre-registered positive controls", notes, df, extra)


if __name__ == "__main__":
    S.TAB_DIR.mkdir(parents=True, exist_ok=True)
    table1()
    table2()
    table_s1()
