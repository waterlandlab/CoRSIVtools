# CoRSIVtools

<!-- badges: start -->
[![check](https://github.com/waterlandlab/CoRSIVtools/actions/workflows/check.yaml/badge.svg)](https://github.com/waterlandlab/CoRSIVtools/actions/workflows/check.yaml)
<!-- badges: end -->

DNA methylation analysis restricted to **CoRSIVs** (correlated regions of
systemic interindividual variation), for any Illumina methylation array. Given
a normalised beta matrix, sample metadata and a probe manifest, CoRSIVtools:

- maps probes to the bundled CoRSIV 2.0 regions (hg38). Replicate probes are
  merged by coordinate.
- summarises methylation per CoRSIV, or analyses each CoRSIV CpG.
- tests differential methylation with limma under any covariate or paired
  design.
- compares CoRSIV results with matched random control regions, the rest of the
  array, or your own control regions.
- records full provenance and draws Manhattan, volcano and per-CoRSIV plots.

**Documentation:** <https://waterlandlab.github.io/CoRSIVtools/>. It includes
the [workflow vignette](https://waterlandlab.github.io/CoRSIVtools/articles/CoRSIVtools.html)
and the [statistical methods](https://waterlandlab.github.io/CoRSIVtools/articles/methods.html).

## Installation

```r
install.packages("BiocManager")
BiocManager::install("waterlandlab/CoRSIVtools")
```

Requires R ≥ 4.5 and Bioconductor ≥ 3.22. `BSgenome.Hsapiens.UCSC.hg38`, about
800 MB, is installed as a dependency; it is used for CpG-density matching.

## Probe manifests

| Array | How to supply it |
|---|---|
| EPICv2 | `read_manifest("EPICv2")`. Needs `BiocManager::install("IlluminaHumanMethylationEPICv2anno.20a1.hg38")`. |
| EPICv1 | `read_manifest("EPICv1")`. Needs `remotes::install_github("achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38")`. |
| 450K, EPIC-Xtra, custom arrays | `read_manifest("path/to/manifest.bed")`, a 4-column hg38 BED (`chr, start, end, probe_id`), or an Illumina manifest CSV. |

## Quick start

```r
library(CoRSIVtools)
se   <- from_matrix(betas, meta, sample_col = "Sample_ID")    # betas: probes x samples
se   <- annotate_probes(se, read_manifest("EPICv2"))
sets <- define_sets(se)                                       # CoRSIVs only (default)
coverage_report(sets)                                         # how many CoRSIVs the array covers
reg  <- summarise_regions(se, sets)
res  <- test_dm(reg, design = ~ group + age + sex, coef = "groupCase")
cmp  <- compare_to_background(res)
plot_manhattan(res); plot_volcano(res); plot_background(cmp)
```

The [workflow vignette](https://waterlandlab.github.io/CoRSIVtools/articles/CoRSIVtools.html)
(`vignette("CoRSIVtools", package = "CoRSIVtools")`) walks through a full
analysis of a public EPIC dataset, GSE131433: ART vs natural conception, at
birth and in adulthood. It covers GEO cleaning, QC, matched backgrounds and
interpretation. The scripts that produce it are in `inst/scripts/`
(`gse131433_prepare.R`, `gse131433_analysis.R`). For a small example that runs
in seconds, use `simulate_corsiv_data()`.

**Seeds:** draw *r* of the matched background uses seed `seed + r`. For two
analyses to have independent background draws, their seeds must differ by at
least `n_draws`.

## Contributing / adding features

`main` is protected. Every change goes through a pull request, which needs the
`check` CI job to pass (about 5 minutes).

### Where things go

| What you're adding | File(s) |
|---|---|
| New function | `R/<topic>.R` with a roxygen `@export` block. Reuse an existing file if the topic fits (`plots.R`, `annotate.R`, `test_dm.R`, `qc_*.R`). |
| Its tests | `tests/testthat/test-<topic>.R`. Use small simulated data from `simulate_corsiv_data()` or `helper-qc.R`, never real data. |
| Website reference page | Add the function name to a section of `_pkgdown.yml`. **The website build fails if an exported function is missing.** |
| New dependency | `DESCRIPTION`. Put heavy or optional packages in **Suggests** and guard them with `.pkg_available()` (see EpiDISH in `qc_cells.R`). |
| Bundled small data | `inst/extdata/`. Record its source in `inst/extdata/README.md`. |
| User guide / methods text | `vignettes/CoRSIVtools.Rmd` (usage) or `vignettes/methods.Rmd` (formal methods). |
| Changelog | `NEWS.md`. Bump `Version:` in `DESCRIPTION`: 0.x.y for fixes, 0.(x+1).0 for features. |

### Workflow

```bash
git checkout main && git pull
git checkout -b feat/<short-name>          # e.g. feat/annotate-results

# edit R/, tests/, _pkgdown.yml, NEWS.md, DESCRIPTION
Rscript -e 'devtools::document()'          # regenerates man/ and NAMESPACE
Rscript -e 'devtools::test()'              # fast loop
Rscript -e 'devtools::check()'             # aim for 0 errors / 0 warnings / 0 notes

git add -A && git commit -m "Add <feature>: <why>"
git push -u origin feat/<short-name>
gh pr create --fill
```

`devtools::check()` needs pandoc. If it is not found, set `RSTUDIO_PANDOC` to a
pandoc directory, e.g. `RSTUDIO_PANDOC=/Applications/quarto/bin/tools/aarch64`
with Quarto on Apple Silicon.

Once CI passes, merge the PR on GitHub or with `gh pr merge --merge`. Merging to
`main` redeploys the documentation website automatically.

To use the new version in an analysis project:

```r
remotes::install_github("waterlandlab/CoRSIVtools")   # or devtools::install("path/to/CoRSIVtools")
```

### House rules

- **Keep defaults fast.** Expensive steps must be opt-in: background draws,
  large downloads, many permutations.
- Functions take and return a `SummarizedExperiment` or a plain `data.frame`,
  like the existing ones.
- Plots return ggplot objects in the house palette (see `plots.R`).
- Never commit participant-level data, including in tests or `inst/extdata`.
- For a larger feature, write a short design note in `docs/` before coding.

## Citation

`citation("CoRSIVtools")`. Please also cite the CoRSIV atlas: Gunasekara CJ,
et al. (2019) *Genome Biology* 20:105.

## Licence

MIT © Baylor College of Medicine.
