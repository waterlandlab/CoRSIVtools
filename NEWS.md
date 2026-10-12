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
