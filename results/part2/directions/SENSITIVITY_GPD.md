# GPD-tail sensitivity: pre-declaration (PI ruling R3)

Committed **before** any of the quantities below is computed. This is a **sensitivity analysis**: it writes no
trial line, changes no registered output, and the registered `p_final` (permutation p with the GPD tail below 1e-4)
stays primary. Discoveries and q-values of every direction are unchanged.

## Why

With 10,000 permutations the empirical p cannot go below 1/10,001. Below 1e-4 the registered rule fits a generalised
Pareto tail, which can return p far below the analytic HC3 p (e.g. H3a: 435 of 1,145 discoveries have p_perm < p;
one 1e-24 against an analytic 1e-10). The extrapolated values must not drive ranking or reporting.

## Sensitivity (computed once per direction from the committed discovery tables)

For each of **H3a_lineage_atlas, H3b_driver_family, H1_within_complex_GI**, and, labelled the same way, for the
exploratory H1S_silencing_GI and H2_submodule_divergence once their single runs are committed:

1. `p_cons = max(p_final, p)` per evaluable test, where `p` is the analytic HC3/Satterthwaite p of the registered
   alternative. For tests not re-tested by permutation, `p_final = p`, so `p_cons = p`.
2. Re-run the entry's registered FDR procedure on `p_cons` with the engine function `assoc._fdr` (same families,
   q = 0.05): Benjamini–Bogomolov (families = contexts) for H3a/H3b; pooled BH for H1/H1S; stage-wise BB for H2 with
   the same parent selection (`assoc.parent_selection_from` the committed primary H3a/H3b discovery tables).
3. Report: number of primary discoveries, of conservative discoveries, their overlap (both / primary only /
   conservative only), overall and per context; for H1/H1S, the event axes whose discovery flag would change (class
   calls are not recomputed).
4. Input check: each discovery table is read from `data_cache/part2/assoc/<ID>/discovery.parquet` only if its
   SHA-256 equals `full_table_sha256` in the committed `results/part2/assoc/<ID>/discovery_summary.json`.

## Reporting tables (all directions)

`results/part2/directions/<ID>/discoveries_report.tsv`: the primary discoveries ranked by |t| (never by GPD p), with
β, 95% CI, analytic p, `p_perm_report` (`"<1e-4"` whenever the GPD branch was used, i.e. no permutation statistic
reached the observed one so the empirical p is 1/10,001; otherwise the empirical permutation p), the registered q, and
`discovery_cons` (the conservative flag from step 2). Script: `scripts/45_gpd_sensitivity.py`. Outputs:
`<ID>/gpd_sensitivity.json` and `gpd_sensitivity_summary.json` in this folder.

## Disclosure

Before this note was committed, the coordinator relayed (ruling R5) the independent reviewer's counts for the same
quantity: p = max(p_perm, p_analytic) loses 32 H3a, 31 H3b and 0 H1 discoveries. The reviewer's file was not opened.
The computation here is fixed above and is run mechanically; agreement or disagreement is reported as found.
