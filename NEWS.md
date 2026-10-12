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
