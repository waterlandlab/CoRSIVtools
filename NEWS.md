# CoRSIVtools 0.3.0

* New default plot style for every CoRSIVtools plot, in the manner of Nature Genetics figures:
  * large, plain sans-serif text that stays legible when a panel is shrunk into a multi-panel
    figure;
  * point sizes and line widths doubled;
  * legends on top;
  * short panel titles, with details in subtitles.
* Contributors: the style is set in `plots.R` via `.theme_corsiv()`, `.pt()` and `.lw()`.
* `qc_samples()` hides per-sample axis labels above 60 samples.
* Requires ggplot2 >= 3.5.0.

# CoRSIVtools 0.2.2

* **Bug fix:** `read_manifest("EPICv1")` returned hg19 coordinates.
  `IlluminaHumanMethylationEPICanno.ilm10b5.hg38` stores Illumina's hg19 `MAPINFO` in `Locations`
  despite its name. hg38 positions are now read from its `Other` table
  (`CHR_hg38`, `Start_hg38 + 1`); 99.9% of sampled cg probes land on a CpG in hg38. EPICv1 results
  from earlier versions overlapped hg19 probe positions with hg38 CoRSIVs and should be re-run.
  `"EPICv2"` was not affected.
* The workflow vignette now analyses a real public dataset, GSE131433 (ART vs natural conception,
  at birth and in adulthood), from precomputed results. Reproduce it with
  `inst/scripts/gse131433_prepare.R` and `inst/scripts/gse131433_analysis.R`.

# CoRSIVtools 0.2.1

* `qc_pca()` uses autosomal CpGs only by default (`autosomes = TRUE`). In mixed-sex cohorts,
  sex-chromosome CpGs otherwise dominate PC1 and hide batch effects.
* `qc_report()` gains `technical`: only these variables (e.g. chip, row, processing date) can flag
  the batch check. Biological variables such as sex or age no longer raise a batch flag.

# CoRSIVtools 0.2.0

* New QC suite, run via `qc_report()` or individually:
  * `qc_samples()`: per-sample missingness and beta distributions, with outlier samples flagged.
  * `qc_pca()` and `qc_batch()`: PCA coloured by metadata, and the association of each PC with batch
    and biological variables.
  * `qc_cells()` and `estimate_cell_fractions()`: blood cell composition against groups and PCs.
    Fractions are estimated with EpiDISH (in Suggests) when not supplied; reference CpGs are matched
    by unambiguous `cg` name only.
  * `qc_spatc1l()`: probe-level correlation of the CoRSIVs at the two ends of SPATC1L; expected
    anti-correlated.
  * `qc_vtrna2()`: bimodality of the VTRNA2-1 CoRSIV; expected bimodal.

# CoRSIVtools 0.1.2

* `define_sets()` now defaults to `background = "none"`: a CoRSIV-only analysis with no
  background sampling (and no BSgenome use). Pass `background = "matched"`, `"array"` or a
  `GRanges` to compare against a background. Calls that relied on the old default
  (`"matched"`) must now request it explicitly.

# CoRSIVtools 0.1.1

* Bug fix: the bundled CoRSIV 2.0 file stores 1-based, end-inclusive coordinates
  (whole 100-bp bins), but `corsiv_regions()` read it as 0-based BED and added 1 to
  every start. Each CoRSIV lost its first base (e.g. chr13:110812302-110814100
  instead of chr13:110812301-110814100), and array CpGs on that base were left out
  (36 CpGs on EPIC-Xtra). `region_id`s change accordingly; re-run analyses made
  with 0.1.0. `write_bed()` still writes standard 0-based BED.

# CoRSIVtools 0.1.0

* First release: from_matrix(), read_manifest(), annotate_probes(),
  corsiv_regions() (CoRSIV 2.0, hg38), define_sets() with matched / array /
  user backgrounds, coverage_report(), summarise_regions(), test_dm() at
  CoRSIV and CpG level, compare_to_background(), plots, write_bed().
* Background draws are seeded as `seed + r` for draw r (spec section 7.2).
  Two analyses share draws unless their base seeds differ by at least
  `n_draws`; use e.g. `seed = 1` and `seed = 1001` for a 100-draw sensitivity
  analysis.

## Known limitations / planned for 0.2

* Replace the `seed + r` scheme with independent per-draw random streams
  (e.g. L'Ecuyer-CMRG substreams derived from the base seed), so that
  different base seeds never share draws.
