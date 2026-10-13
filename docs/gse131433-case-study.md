# GSE131433 case study: data preparation and analysis

This note records how the public dataset used in the package vignette was cleaned
and analysed. It is not part of the built package (`docs/` is Rbuildignored).

## Source

- GEO series **GSE131433** (CHART cohort, Murdoch Children's Research Institute).
  Novakovic B, et al. (2019) *Nat Commun* 10:3922 (PMID 31477727); follow-up
  PMID 37085889.
- Illumina Infinium MethylationEPIC (GPL21145). 440 arrays: Guthrie (newborn
  heel-prick) blood spots and adult whole blood from individuals conceived by
  assisted reproductive technology (ART) and naturally conceived controls.
- Files used, as downloaded from GEO and never modified:
  - `GSE131433_series_matrix.txt` (sample characteristics).
  - `GSE131433_BorisNovakovicCHART_matrix_processed.csv` (supplementary file;
    submitter-normalised SWAN betas and detection P values; 722,292 probes, so
    some probe filtering was already applied by the submitters).

## Preparation (`inst/scripts/gse131433_prepare.R`)

1. **Metadata.** Sample characteristics are parsed from the series matrix.
   Derived variables: `timepoint` (birth/adult), `conception` (natural/ART,
   reference = natural), `sex`, `art_subtype`, `donor`, and the Sentrix
   `slide` and `position` parsed from the IDAT file names (`array_row` = row of
   the position). Each sample title encodes timepoint and conception; these are
   cross-checked against the characteristics fields, and the script stops on
   any disagreement.
2. **Join.** Beta columns are matched to metadata on the full sample title
   (exact string match). Matching must be one-to-one or the script stops.
3. **Masking.** Beta values with detection P > 0.01 are set to NA.
4. **Probe filter.** Probes masked in > 5% of samples, and non-CpG probes
   (`ch.`/`rs` IDs), are dropped.
5. **Annotation.** Probes are mapped to hg38 with the EPIC v1 manifest
   (`IlluminaHumanMethylationEPICanno.ilm10b5.hg38`) through `read_manifest("EPICv1")`.
   Replicate probes at one coordinate are averaged by `annotate_probes()`.
   hg38 positions come from the package's `Other$CHR_hg38` / `Start_hg38 + 1`.
   Its `Locations` table is hg19; using it was a bug fixed in v0.2.2, and it
   left SPATC1L and VTRNA2-1 with no CpGs.

The submitters' SWAN normalisation is kept; no re-normalisation from IDATs is done.

## Analysis (`inst/scripts/gse131433_analysis.R`)

- **Cell composition.** EpiDISH RPC with the 7-type `centDHSbloodDMC`
  reference (Teschendorff et al. 2017). The adult-blood reference is also
  applied to newborn blood spots. It contains no nucleated red blood cells,
  so newborn estimates are approximate.
- **QC.** `qc_report()` on all arrays; technical variables are `slide` and
  `array_row`. Samples are excluded only above 5% masked CpGs; none were.
  The relative rule (> 3 MADs above the median) flags 89 samples, but median
  missingness is ~0.005%, so the rule is over-sensitive here. All flagged
  samples have < 2.3% missing.
- **Sets.** `define_sets(background = "matched", n_draws = 100, seed = 1)`.
  Both timepoints share the same sets and region summaries, which depend only on
  array coverage and on each sample's own values.
- **Model.** Fitted separately at birth and in adulthood. Each fit is a limma
  model on M-values of CoRSIV means:

  `M ~ conception + sex + array_row + B + NK + CD4T + CD8T + Mono`,

  testing `conceptionART`. Neutrophils and eosinophils are the omitted
  reference, since the fractions sum to ~1.
- **Why `array_row`.** Within each timepoint, chip row is confounded with
  conception: naturally conceived samples sit mostly on rows R05–R07 at birth
  and R06/R08 in adulthood (simulated Fisher P < 1e-4). Without it, 24% (birth)
  and 37% (adulthood) of matched background regions had P < 0.05 for ART. With
  it, about 3% do. Slide (56–57 levels) is not adjusted for; it is not
  associated with conception. Background draws are fitted separately; see
  `vignettes/methods.Rmd`.

## Limitations

- Adult age is not provided in the GEO record and cannot be adjusted for. Adults
  were sampled as young adults (see Novakovic et al. 2019 for the age range).
- The ART group pools fresh-embryo, frozen-embryo and GIFT conceptions.
- Analyses are cross-sectional per timepoint. The paired birth→adult structure
  (most donors have both arrays) is not used here.
- Newborn and adult samples differ in tissue (dried blood spot vs whole blood),
  so the two timepoints should not be compared directly.
