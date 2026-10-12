# CoRSIV2.0_hg38.bed

- Contents: CoRSIV 2.0 regions (autosomes + chrX; ESS; male and female), hg38,
  0-based BED, 4 columns (`chr, start, end, source_corsivs`). Overlapping
  population-specific CoRSIVs were merged with `bedtools merge`; the name column
  lists every merged source CoRSIV, comma-separated.
- Regions: 6,669 (320 on chrX). Non-overlapping and non-abutting.
- Original file name: `CoRSIV2.0_ESS_chrX_male_female_all_cooridnates_100bp_hg38.bedtools_merged_fixed.bed`
- md5: `ca1b9a6c50e78e9a43a92f2e56400fb2`
- Source: Waterland Lab, Baylor College of Medicine (CoRSIV 2.0; publication forthcoming).
- Licence: MIT, same as the package (approved 2026-10-12).
- How it was made: see `inst/scripts/make_corsiv_bed.R`.
