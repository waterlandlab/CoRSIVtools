# CoRSIVtools QC suite — design spec (v0.2.0)

- **Date:** 2026-10-12
- **Status:** Draft, awaiting review
- **Builds on:** `2026-10-12-corsivtools-design.md` (v0.1). All v0.1 conventions apply:
  SummarizedExperiment in, tidy DataFrames and ggplots out, `region_id` keys, user-facing
  errors, and no participant data in the package.

## Purpose
Standard array QC that a user runs before (or alongside) a CoRSIV analysis. It answers four questions:
1. Are any samples poor?
2. Are there batch effects, and which variables drive the main axes of variation?
3. Are blood cell proportions confounded with the comparison of interest?
4. Do known CoRSIV behaviours appear in the data, showing the data behave as expected?

## Inputs
- **Data:** an annotated SE from `annotate_probes()`, with sample metadata in `colData`.
- **Metadata variables:** chosen by the user, e.g. `group`, `chip`, `position`, `sex`, `age`.
- **Cell fractions:** optional `colData` columns. If none are given, `estimate_cell_fractions()`
  computes them.

## Functions

| Function | Computation | Returns |
|---|---|---|
| `qc_samples(se)` | Per sample: fraction of NA CpGs, median β, the β density curve. A sample is flagged if its missing fraction is > 3 MAD above the median, or > 0.05 absolute. | `list(table, plot_missing, plot_density)` |
| `qc_pca(se, vars, n_top = 10000, features = c("all", "corsiv"), n_pcs = 10)` | PCA (`prcomp`, centred, unscaled) on M-values of the `n_top` most variable CpGs. Only complete CpGs are used (CpGs with NA in any sample are dropped and their number is reported). `features = "corsiv"` restricts to CpGs in CoRSIV 2.0. | `list(scores, var_explained, plots)`; one PC1–PC2 and one PC3–PC4 plot per variable in `vars` |
| `qc_batch(pca, se, vars)` | For each of PCs 1–`n_pcs` and each variable: linear model `PC ~ var` (factor → F-test; numeric → slope t-test). | `list(table, plot)`; heatmap of −log10 P, cells annotated with the % variance of that PC |
| `estimate_cell_fractions(se, ref = c("cent12CT", "centDHSbloodDMC"))` | EpiDISH RPC on β. Each reference `cg` ID is matched to the SE's `probe_ids` by the base name before the first `_`. A reference CpG is used only when its base name maps to **exactly one** CpG coordinate in the SE. The numbers matched, ambiguous and missing are reported. EpiDISH is in Suggests. | `DataFrame` of fractions (samples × cell types); match counts in `metadata()` |
| `qc_cells(se, cells, group, pair = NULL)` | For each cell type, a group test: paired t-test if `pair` is given, else Wilcoxon (2 groups) or Kruskal–Wallis (> 2). Also each cell fraction vs each of PCs 1–10, via `qc_batch` logic. | `list(table, plot_box, plot_pc)` |
| `qc_spatc1l(se)` | The two CoRSIVs at the ends of SPATC1L (bundled hg38 coordinates: `chr21:46160701-46162600` and `chr21:46184001-46185600`). Probe-level correlation matrix (all CpGs in both regions, ordered by position) and the correlation of the two regional means. **Pass** if the regional-mean r < 0 and the median between-CoRSIV probe correlation is < 0. | `list(r_regional, median_between, pass, plot_corr, plot_scatter)` |
| `qc_vtrna2(se)` | CoRSIV `chr5:136079901-136081000` overlapping VTRNA2-1: per-sample mean β and the probe-level distribution. **Pass** if bimodal: a 2-component Gaussian mixture (`stats`-only EM, no new dependency) beats one component by BIC, and the component means differ by > 0.2. Reports each sample's mode (low / high). | `list(table, pass, plot_hist, plot_probes)` |
| `qc_report(se, vars, cells = NULL, group = NULL, pair = NULL)` | Runs all of the above, skipping steps whose inputs are missing. | `list(summary, results)`; `summary` has one row per check: check, status (pass / flag / skipped), key metric |

- **Control loci:** the coordinates are stored in `inst/extdata/qc_loci.csv` (locus, region_id, role,
  expected pattern). They are matched to the CoRSIV 2.0 coordinates, and a check is skipped with a
  message if the array has fewer than 3 CpGs in a control region.
- **Plots:** lab figure style. Correlation matrices use a ggplot tile heatmap (no corrplot
  dependency), red–white–blue with limits −1 to 1.

## Testing
- **Simulated data:** tests run on `simulate_corsiv_data()`. Without seeding control-locus behaviour
  the simulation can't reproduce SPATC1L/VTRNA2-1 patterns, so `qc_spatc1l`/`qc_vtrna2` are tested
  with small hand-built SEs at the bundled coordinates:
  - an anti-correlated pair → pass;
  - a positively correlated pair → flag;
  - a bimodal fixture → pass;
  - a unimodal fixture → flag.
- **`qc_batch`:** detects an injected chip effect.
- **`estimate_cell_fractions`:** ambiguous base names are excluded, and the counts are reported (fixture).
- **Real-data validation (Addanki), not in package tests:** SPATC1L r ≈ −0.77 (pass); VTRNA2-1
  bimodal (10 low / 64 high, pass).

## Out of scope
- An HTML report template (`qc_report` returns objects; the user saves figures).
- Sex prediction from chrX/chrY (EPIC-Xtra has little sex-chromosome coverage).
- Automatic batch correction. QC reports batch effects; adjustment stays in the user's `design`.

## Dependencies
- **No new Imports.**
- **EpiDISH** in Suggests, used only by `estimate_cell_fractions()`.
