# CoRSIVtools — design spec (v0.1)

- **Date:** 2026-10-12
- **Status:** Approved 2026-10-12
- **Repository:** `waterlandlab/CoRSIVtools` (public GitHub). Local clone at `~/code/CoRSIVtools`, outside OneDrive.
- **Licence:** MIT. Copyright holder: Baylor College of Medicine.
- **First user:** Baylor_Addanki_10022026 (EPIC-Xtra, 74 subjects, HCC case/control). The analysis lives in its own repo and pins the CoRSIVtools version with `renv`.

## 1. Purpose

CoRSIVtools is an R package for DNA methylation analysis **restricted to CoRSIVs**
(correlated regions of systemic interindividual variation). A lab provides:

1. a normalised beta matrix,
2. a sample metadata table, and
3. a probe manifest (built in for 450K, EPICv1 and EPICv2; supplied by the user for EPIC-Xtra or other custom arrays).

The package then:

- maps probes to CoRSIVs,
- summarises methylation per CoRSIV,
- tests for differential methylation with any covariate design, and
- compares CoRSIV results against a background.

The package grows by adding analysis modules (DMRcate, PCA, mixOmics, and other statistical methods). Every module works only on CoRSIV probes.

### Goals (v0.1)

- One working path from input to results: input → annotation → CoRSIV/background sets → region summaries → limma test → background comparison → plots.
- Works with any array, given a manifest.
- Every result is reproducible and records its provenance.
- Structured so it can later pass Bioconductor submission (SummarizedExperiment-based, BiocCheck-clean).

### Non-goals (v0.1), planned for later versions

- DMR calling (DMRcate), PCA/association screens, multi-omics integration (mixOmics).
- Estimating cell composition. v0.1 accepts cell fractions as covariates in `colData`.
- Sex-stratified chrX models, mixed models for repeated measures, beta regression, per-CpG test combination (e.g. ACAT).
- Whole-array analysis. **The rest of the array is only ever used as a background.**

## 2. Facts from the first dataset that shaped the design

All numbers come from Baylor_Addanki_10022026, measured 2026-10-12.

| Item | Value |
|---|---|
| Beta matrix | `data.frame`, 1,020,136 probes × 74 samples; IDs like `cg06185909_TC11` |
| Xtra manifest | 4-col BED, 1-bp, 0-based, hg38; 1,009,415 probes; every ID is in the betas |
| Beta IDs with no coordinates | 10,721 (7,967 `cg`, 2,294 `Ne*` controls, others) |
| EPICv2 Bioconductor annotation | covers only 4,249 Xtra IDs (0.4%), so a user-supplied manifest is required |
| CoRSIV 2.0 BED | 6,669 merged, non-overlapping, non-abutting regions; mean width 409 bp; 320 on chrX |
| CoRSIVs with ≥1 / ≥2 / ≥3 / ≥5 probes | 4,517 / 3,753 / 3,083 / 1,230 |
| Probes in CoRSIVs / outside CoRSIVs | 18,649 / 990,766 |
| chrX CoRSIVs covered | 18 / 320 |
| Probes falling in more than one CoRSIV | 0 |
| Original CoRSIV names | **not unique**: 6,590 distinct names for 6,669 regions |
| Replicate probes sharing a `cg` base name | 20,983. Of 4,246 IDs shared with EPICv2, 60 (`_TC21`/`_BC21`) map to a different position (mostly 1 bp off; a few to distant loci) |

## 3. Architecture

```
 from_matrix(betas, meta, sample_col)          → probe-level SummarizedExperiment
        │
 annotate_probes(se, manifest)                 → rowRanges = probe coordinates
        │                                         replicates collapsed by coordinate
 define_sets(se, corsivs, background, …)       → CoRSIVSets
        │   target:     CoRSIV regions + their probes
        │   background: "matched" | "array" | user GRanges
        │
 ┌──────┴───────────────────────────┐
 CpG level                          CoRSIV level
 test_dm(level = "cpg")             summarise_regions() → region-level SE
                                    test_dm(level = "corsiv")
 └──────┬───────────────────────────┘
        │   one tidy result table: target rows + background rows (tagged by draw)
 compare_to_background(res)
 plot_manhattan / plot_volcano / plot_corsiv / plot_background
 coverage_report(sets); write_bed(res)
```

### Rules for the architecture

1. **SummarizedExperiment is what passes between steps.** Covariates (cell fractions, ancestry PCs, batch) are `colData` columns and enter only through `design`.
2. **Only CoRSIV probes are analysed.** Non-CoRSIV probes are used only as a background, and background results never feed into target statistics.
3. **FDR is calculated within each set.** The target set, and each background draw, gets its own Benjamini–Hochberg adjustment.
4. **Heavy future back-ends (DMRcate, mixOmics) go in `Suggests`** and are loaded through `requireNamespace()`. The core install stays light.
5. **One package for now.** A module that grows large may move to a companion package later without breaking the core API.

## 4. Identifiers

- **`region_id`** (the CoRSIV ID for target regions) is the primary key: a UCSC-style 1-based string, `chr:start-end`. Example: BED `chr1 629001 635000` → `chr1:629002-635000`. The same format is used for background regions. Checked on CoRSIV 2.0: 6,669 regions give 6,669 unique IDs, and ID → coordinates round-trips exactly.
- **`source_corsivs`** keeps the original BED name (e.g. `black_chr1_38,…,white_chr1_15`) in full, as a column. It is never used as a key, because it isn't unique.
- **Probe IDs** are kept exactly as supplied. They are never truncated to the `cg` base name.
- **`write_bed()`** converts back to 0-based BED for bedtools. A unit test checks that BED → ID → BED returns identical coordinates.

## 5. Public API (v0.1)

| Function | Input → output | Notes |
|---|---|---|
| `from_matrix(betas, meta, sample_col)` | matrix/data.frame + data.frame → SE | Coerces the betas to a numeric matrix. Errors on: sample IDs not matching one-to-one, duplicated probe or sample IDs, values outside [0, 1], non-numeric columns. |
| `read_manifest(x)` | path (BED, Illumina CSV) or `"450K"`, `"EPICv1"`, `"EPICv2"` → `GRanges` (1 bp, hg38, `probe_id` column) | The built-ins need the matching Illumina annotation package (Suggests). The genome build must be hg38; other builds error. |
| `annotate_probes(se, manifest)` | SE + `GRanges` → SE with `rowRanges` | Drops probes without coordinates and reports how many. Collapses replicate probes **by identical `chr:pos`** (mean beta), keeping the contributing IDs in `rowData$probe_ids`. |
| `corsiv_regions(version = "2.0", chrX = FALSE)` | → `GRanges` with `region_id`, `source_corsivs` | Reads `inst/extdata/CoRSIV2.0_hg38.bed`. |
| `define_sets(se, corsivs = corsiv_regions(), background = "matched", n_draws = 100, seed = 1, flank = 1000, min_cpgs = 3)` | → `CoRSIVSets` (S4) | `background`: `"matched"`, `"array"`, or a `GRanges`. Holds the target regions, the background regions (with `draw`), the probe→region maps and the matching diagnostics. |
| `coverage_report(sets)` | → `DataFrame` | Number of regions with ≥1/2/3/5 probes, probes-per-region distribution, number dropped by `min_cpgs`, for both target and background. |
| `summarise_regions(se, sets, fun = "mean", max_na_frac = 0.2)` | → region-level SE (`rowRanges` = regions; `rowData$set`, `$draw`, `$n_cpg`) | Uses `min_cpgs` from `sets`. |
| `test_dm(x, design, coef, level = c("corsiv", "cpg"))` | SE (+ sets for `"cpg"`) → `DataFrame` | Schema in §6.4. |
| `compare_to_background(res, threshold = c(P = 0.05))` | → list of `DataFrame`s | Enrichment, distributions, variability (§7.3). |
| `plot_manhattan(res)`, `plot_volcano(res)`, `plot_corsiv(se, id, group)`, `plot_background(cmp)` | → `ggplot` | Lab figure style: palette and font sizes from the lab plotting guide. |
| `write_bed(res, file)` | → 0-based BED with `region_id` in the name column | For bedtools. |

## 6. Statistical defaults

### 6.1 Region definition
- `min_cpgs = 3`, counting distinct CpG coordinates after replicates are collapsed. On Xtra this keeps 3,083 of the 4,517 covered CoRSIVs. Can be changed; the Addanki analysis plan includes a sensitivity analysis at `min_cpgs = 2` (3,753 regions).
- chrX is excluded by default (`chrX = FALSE`). Including it would need sex-stratified models, which are out of scope for v0.1.

### 6.2 Region summary
- For each sample, the mean beta over the region's CpGs (`fun = "mean"`; `"median"` is also available).
- If a sample has fewer than `min_cpgs` non-NA CpGs in a region, that sample's region value is NA.
- A region is dropped if more than `max_na_frac` (0.2) of its samples are NA. The number dropped is reported.

### 6.3 Test
- limma on M-values: `M = log2(β / (1 − β))`, with β clamped to [1e-6, 1 − 1e-6].
- `lmFit` on `model.matrix(design, colData)`, then `eBayes(trend = TRUE, robust = TRUE)`.
- Effect size: `delta_beta`, the coefficient from the same design fitted on the beta scale with `lmFit` (for a two-group factor, this is the difference in group means after adjustment). P values come from the M-value model.
- `coef` must name exactly one column of the design matrix; otherwise it errors and lists the available names.
- Samples with NA in any design variable are dropped, and the count is reported.
- A design matrix that is not of full rank errors, naming the aliased columns.
- limma is fitted **separately for each (`set`, `draw`) block**, so target moderation and FDR never depend on the background.

### 6.4 Result schema (one row per region or CpG, per set and draw)
`id, region_id, set ("target"/"background"), draw (NA for target), chr, start, end, n_cpg, source_corsivs, delta_beta, logFC_M, AveM, t, P, FDR, sd_beta`.
At CoRSIV level, `id == region_id`. At CpG level, `id` is the CpG (`chr:pos`) and `region_id` is the region that contains it. `sd_beta` is the SD of beta across the samples used (needed for the variability comparison in §7.3).
FDR is BH within each (`set`, `draw`) group.

### 6.5 Provenance
`metadata(result)` records: CoRSIVtools version, CoRSIV set version, manifest md5 (or built-in name), `background` mode, `n_draws`, `seed`, `flank`, `min_cpgs`, `fun`, `max_na_frac`, the design formula, `coef`, the number of samples used, R and Bioconductor versions, and a timestamp.

## 7. Background

### 7.1 Modes

| Mode | Contents | Levels |
|---|---|---|
| `"matched"` (default) | `n_draws` sets of random control regions matched to the covered CoRSIVs | CoRSIV and CpG |
| `"array"` | every mapped probe outside the CoRSIVs ± `flank` | CpG only. Calling `test_dm(level = "corsiv")` with this mode errors with an explanation. |
| a `GRanges` | the user's own control regions; must not overlap the CoRSIVs | CoRSIV and CpG |

### 7.2 Matched sampling algorithm

For each draw `r` in `1..n_draws`, using seed `seed + r`, and for each target CoRSIV `i` (width `w_i`, probe count `k_i`, CpG density `d_i`, chromosome type `c_i`):

1. **Exclusion zone:** every CoRSIV ± `flank` bp, plus the control regions already accepted in this draw.
2. **Propose a candidate:** pick a random mapped probe on a chromosome of type `c_i` (autosome or chrX) outside the exclusion zone. Place a window of width `w_i` with a uniformly random offset so that it contains that probe.
3. **Accept the candidate only if all four criteria hold:**
   - **Width:** exactly `w_i` (true by construction).
   - **Probe count:** distinct probe coordinates in the window equal `k_i` if `k_i < 10`, or are within ±1 of `k_i` if `k_i ≥ 10`.
   - **Chromosome type:** same as `c_i`.
   - **CpG density:** CpGs per 100 bp, counted with `Biostrings` on `BSgenome.Hsapiens.UCSC.hg38`, within ±20% of `d_i` (absolute floor ±0.25 CpG/100 bp).
4. **Retry limit:** up to `max_tries = 1000` proposals per target. If it still fails, record the target as unmatched for that draw.
5. **Diagnostics** go in `CoRSIVSets`: for each draw, the matched/unmatched counts and the distributions of all four criteria for target versus control. `define_sets()` warns if any draw has more than 5% unmatched.

The same seed gives identical draws. The control regions themselves are stored, so the comparison can be checked and exported with `write_bed()`.

### 7.3 `compare_to_background()`

- **Enrichment:** the fraction of target regions with P < threshold, compared with each draw.
  - Empirical p = (1 + #draws with a fraction ≥ the target's) / (1 + `n_draws`).
  - Fisher's exact test of target versus the pooled background is also reported. It is the main test in `"array"` mode.
  - FDR < 0.05 counts are reported alongside.
- **Distributions:** quantiles of |delta_beta| and of P, and data for a QQ plot, for the target versus the background draws.
- **Variability between individuals:** per-region SD of beta across samples, target median versus the distribution over draws, with an empirical p. This also checks the data: CoRSIVs are expected to vary more between individuals than matched background regions.

## 8. Package layout

```
CoRSIVtools/
├── DESCRIPTION          License: MIT + file LICENSE
├── LICENSE              YEAR: 2026 / COPYRIGHT HOLDER: Baylor College of Medicine
├── NAMESPACE            roxygen2
├── NEWS.md
├── R/  input.R manifest.R annotate.R corsiv.R sets.R background.R
│       summarise.R test_dm.R compare.R plots.R provenance.R
├── inst/extdata/CoRSIV2.0_hg38.bed  (+ README: source, md5, licence)
├── inst/scripts/make_corsiv_bed.R
├── inst/CITATION
├── vignettes/CoRSIVtools.Rmd   (end-to-end on simulated data)
├── vignettes/methods.Rmd       (formal methods: model, assumptions, limitations, references)
├── tests/testthat/helper-simulate.R, test-*.R (one per R/ file)
└── .github/workflows/  R-CMD-check.yaml, BiocCheck.yaml, test-coverage.yaml
```

**Imports:** SummarizedExperiment, S4Vectors, GenomicRanges, IRanges, BiocGenerics, GenomeInfoDb, limma, BSgenome.Hsapiens.UCSC.hg38, Biostrings, ggplot2, withr, methods, stats, utils, tools.
**Also exported:** `simulate_corsiv_data()`, used by the tests and the vignette, plus the accessors `targetRegions()`, `backgroundRegions()` and `matchDiagnostics()` for `CoRSIVSets`.
**Suggests:** IlluminaHumanMethylation450kanno.ilmn12.hg19, IlluminaHumanMethylationEPICanno.ilm10b5.hg38, IlluminaHumanMethylationEPICv2anno.20a1.hg38, testthat (≥ 3e), knitr, rmarkdown, BiocStyle.
**450K:** the Bioconductor 450K annotation package is hg19-only, so `read_manifest("450K")` errors and points users to an hg38 450K manifest file (e.g. InfiniumAnnotation), which they then pass as a path.

## 9. Testing

- **Approach:** test-driven development, with testthat 3rd edition. All unit tests use **simulated data only**; no participant data is ever in the package.
- **`helper-simulate.R`** builds a small dataset: about 40 CoRSIV-like regions plus candidate background regions on a chr21 slice, 20 samples in two groups, with a known shift in Δβ added to a known subset of regions. Real chr21 coordinates are used so the BSgenome CpG-density step runs for real.

| Area | Tests |
|---|---|
| Correctness | The regions with an added shift are recovered (FDR < 0.05). Under the null with no shift, the fraction with P < 0.05 falls in a binomial 99% interval around 0.05 (fixed seed). |
| Matching | Every accepted control meets all four criteria. No control overlaps a CoRSIV ± flank or another control in the same draw. The same seed gives identical draws, and a different seed gives different ones. |
| IDs | BED → `region_id` → BED is identical. `source_corsivs` survives unchanged. Replicates are collapsed by coordinate (a fixture with two replicate IDs at different positions stays as 2 CpGs). |
| Inputs | Each `from_matrix()` error condition. Probes missing from the manifest are reported. A `coef` not in the design errors. |
| FDR | Target FDR is identical whichever background mode is used. |
| Provenance | Every field in §6.5 is present. |
| Plots | Each returns a `ggplot` and builds without error. |

**CI (GitHub Actions):** `R CMD check` on macOS and Ubuntu against Bioconductor release; BiocCheck; test coverage reported, target ≥ 80%.
**Real-data validation** happens in the Addanki repo (`code/01_corsiv_dm.R`), never in the package.

## 10. Versioning and release

- `0.1.0` during GitHub development, following semantic versioning. Each release gets a `NEWS.md` entry and a git tag.
- For Bioconductor submission the version becomes `0.99.0`, and all BiocCheck notes are resolved.

## 11. Example (the first user)

```r
library(CoRSIVtools)
load("data/raw/betas_Xtra_filtered.RData")
meta <- read.csv("data/raw/meta_data_files.csv", check.names = FALSE)
se   <- from_matrix(betas_Xtra_filtered, meta, sample_col = "<column matching colnames(betas)>")
se   <- annotate_probes(se, read_manifest("data/raw/XTRACoRSIV1_SS_20156551_A1.bed"))
sets <- define_sets(se, background = "matched", n_draws = 100, seed = 1)
reg  <- summarise_regions(se, sets)
res  <- test_dm(reg, design = ~ HCC + age + sex, coef = "HCCCase")
cmp  <- compare_to_background(res)
plot_manhattan(res); plot_volcano(res); plot_background(cmp)
```

The Addanki metadata has 1,008 rows but 74 array samples, and it's keyed by `Patient ID`, not by array sample ID. The join from metadata rows to array columns has to be built and checked in the Addanki repo before `from_matrix()` is called. This is a task for the user's analysis, not for the package.

## 12. Open items

1. **CoRSIV 2.0 citation (resolved 2026-10-12):** not yet published. `inst/CITATION` cites Gunasekara CJ et al., "A genomic atlas of systemic interindividual epigenetic variation in humans", *Genome Biology* 20:105 (2019), doi:10.1186/s13059-019-1708-1, and marks the CoRSIV 2.0 citation as forthcoming.
2. **Provenance of the bundled BED:** `inst/scripts/make_corsiv_bed.R` needs to record how `CoRSIV2.0_ESS_chrX_male_female_all_cooridnates_100bp_hg38.bedtools_merged_fixed.bed` was produced from the CoRSIV 2.0 release. The lab to supply the source steps.
3. **Sharing the vendor manifest:** the Xtra manifest is not bundled. The vignette tells users to get it from the vendor.
4. **Before going public:** check with Baylor's technology transfer office or Dr. Waterland that the code can be released openly.
