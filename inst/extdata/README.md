# CoRSIV2.0_hg38.bed

- Contents: CoRSIV 2.0 regions (autosomes + chrX; ESS; male and female), hg38,
  4 tab-separated columns (`chr, start, end, source_corsivs`). **Coordinates are
  1-based and end-inclusive, not standard 0-based BED**: every region is a whole
  number of 100-bp bins (starts end in `01`, ends in `00`; e.g. chr13
  110812301-110814100 is 18 bins). `corsiv_regions()` reads them as such, and
  `write_bed()` exports standard 0-based BED for bedtools. Overlapping
  population-specific CoRSIVs were merged, apparently with `bedtools merge`
  (inferred from the original file name and the comma-separated, merged name
  column); the name column lists every merged source CoRSIV. The exact upstream
  commands are unconfirmed; the lab will supply them.
- Regions: 6,669 (320 on chrX). Non-overlapping and non-abutting.
- Original file name: `CoRSIV2.0_ESS_chrX_male_female_all_cooridnates_100bp_hg38.bedtools_merged_fixed.bed`
- md5: `ca1b9a6c50e78e9a43a92f2e56400fb2`
- Source: Waterland Lab, Baylor College of Medicine (CoRSIV 2.0; publication forthcoming).
- Licence: MIT, same as the package (approved 2026-10-12).
- How it was made: see `inst/scripts/make_corsiv_bed.R`.
