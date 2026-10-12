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

For a runnable example on simulated data, see
`vignette("CoRSIVtools", package = "CoRSIVtools")`.

**Seeds:** draw *r* of the matched background uses seed `seed + r`. For two
analyses to have independent background draws, their seeds must differ by at
least `n_draws`.

## Citation

`citation("CoRSIVtools")`. Please also cite the CoRSIV atlas: Gunasekara CJ,
et al. (2019) *Genome Biology* 20:105.

## Licence

MIT © Baylor College of Medicine.
