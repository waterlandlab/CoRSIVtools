# CoRSIVtools

CoRSIV-focused DNA methylation analysis for any Illumina array: map probes to
CoRSIVs, summarise per CoRSIV, test differential methylation with limma, and
compare against matched control regions or the rest of the array.

```r
# install.packages("BiocManager")
BiocManager::install("waterlandlab/CoRSIVtools")
```

```r
library(CoRSIVtools)
se   <- from_matrix(betas, meta, sample_col = "Sample_ID")
se   <- annotate_probes(se, read_manifest("EPICv2"))      # or a manifest file path
sets <- define_sets(se, background = "matched", n_draws = 100)
res  <- test_dm(summarise_regions(se, sets), ~ group + age + sex, coef = "groupCase")
cmp  <- compare_to_background(res)
```

See `vignette("CoRSIVtools")` and `vignette("methods", package = "CoRSIVtools")`.
Licence: MIT © Baylor College of Medicine.
