# Clean GEO series GSE131433 (CHART cohort; Novakovic et al. 2019, Nat Commun 10:3922)
# into a CoRSIVtools-ready, annotated SummarizedExperiment.
#
# Inputs, downloaded unchanged from GEO and never modified:
#   GSE131433_series_matrix.txt                 sample metadata (440 EPIC arrays)
#   BorisNovakovicCHART_matrix_processed.csv    SWAN-normalised betas + detection P (supplementary file)
#
# Cleaning steps (see vignettes/CoRSIVtools.Rmd, "Data"):
#   1. Parse sample characteristics; derive donor, timepoint, conception, sex, Sentrix slide/position.
#   2. Join betas to metadata on the exact sample title (checked one-to-one).
#   3. Mask beta values with detection P > DET_P (set to NA).
#   4. Drop probes masked in > MAX_PROBE_FAIL of samples and non-CpG (ch./rs) probes.
#   5. Map probes to hg38 (EPICv1 ilm10b5 manifest) and collapse to CpG coordinates.
#
# Usage: Rscript inst/scripts/gse131433_prepare.R <geo_dir> <out_rds>
# Needs: data.table, IlluminaHumanMethylationEPICanno.ilm10b5.hg38. Peak memory ~20 GB.

DET_P          <- 0.01   # detection P above which a beta value is masked
MAX_PROBE_FAIL <- 0.05   # drop probes masked in more than this fraction of samples

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: Rscript gse131433_prepare.R <geo_dir> <out_rds>")
geo_dir <- args[1]; out_rds <- args[2]
suppressPackageStartupMessages({
  library(CoRSIVtools)
  library(SummarizedExperiment)
})

# ---- 1. sample metadata from the series matrix -------------------------------------------------
sm <- readLines(file.path(geo_dir, "GSE131433_series_matrix.txt"))
field <- function(tag) {
  l <- grep(paste0("^!", tag, "\t"), sm, value = TRUE)
  lapply(l, function(x) scan(text = x, what = "", sep = "\t", quiet = TRUE)[-1])
}
chars <- field("Sample_characteristics_ch1")
char_tab <- do.call(cbind, lapply(chars, function(v) sub("^[^:]+: ", "", v)))
colnames(char_tab) <- vapply(chars, function(v) unique(sub(":.*$", "", v)), character(1))
stopifnot(identical(colnames(char_tab),
                    c("gender", "sample type", "art status", "art subtype", "time of collection")))
idat <- field("Sample_supplementary_file")[[1]]                      # .../GSM..._<slide>_<RxxCxx>_Grn.idat.gz
sentrix <- regmatches(idat, regexec("_([0-9]+)_(R[0-9]+C[0-9]+)_Grn\\.idat", idat))
meta <- data.frame(
  sample_id   = field("Sample_title")[[1]],
  gsm         = field("Sample_geo_accession")[[1]],
  donor       = field("Sample_source_name_ch1")[[1]],
  timepoint   = factor(ifelse(char_tab[, "time of collection"] == "birth", "birth", "adult"),
                       levels = c("birth", "adult")),
  tissue      = char_tab[, "sample type"],
  conception  = factor(ifelse(char_tab[, "art status"] == "ART", "ART", "natural"),
                       levels = c("natural", "ART")),
  art_subtype = char_tab[, "art subtype"],
  sex         = factor(char_tab[, "gender"], levels = c("F", "M")),
  slide       = vapply(sentrix, `[`, "", 2),
  position    = vapply(sentrix, `[`, "", 3),
  stringsAsFactors = FALSE
)
meta$array_row <- sub("C[0-9]+$", "", meta$position)
stopifnot(nrow(meta) == 440L, !anyDuplicated(meta$sample_id), !anyDuplicated(meta$gsm),
          all(nzchar(meta$slide)), !anyNA(meta$timepoint), !anyNA(meta$conception),
          # every title encodes its own timepoint and conception: cross-check against the characteristics
          identical(grepl("^neonate_", meta$sample_id), meta$timepoint == "birth"),
          identical(grepl("_non-ART_", meta$sample_id), meta$conception == "natural"),
          # one array per donor per timepoint
          !anyDuplicated(paste(meta$donor, meta$timepoint)))

# ---- 2. betas and detection P --------------------------------------------------------------------
csv <- file.path(geo_dir, "BorisNovakovicCHART_matrix_processed.csv")
hdr <- names(data.table::fread(csv, nrows = 0L))
stopifnot(identical(hdr[1:2], c("V1", "ID_REF")))
beta_cols <- hdr[-(1:2)][!grepl("_Detection_Pval$", hdr[-(1:2)])]
stopifnot(setequal(beta_cols, meta$sample_id), length(beta_cols) == nrow(meta))
dt <- data.table::fread(csv, select = c("ID_REF", beta_cols, paste0(beta_cols, "_Detection_Pval")),
                        showProgress = FALSE)
probe_ids <- dt$ID_REF
stopifnot(!anyDuplicated(probe_ids))
betas <- as.matrix(dt[, beta_cols, with = FALSE])
detp  <- as.matrix(dt[, paste0(beta_cols, "_Detection_Pval"), with = FALSE])
rm(dt); invisible(gc())
rownames(betas) <- probe_ids
n_raw_probes <- nrow(betas)

# ---- 3-4. mask and filter -------------------------------------------------------------------------
masked <- detp > DET_P
rm(detp); invisible(gc())
betas[masked] <- NA_real_
probe_fail <- rowMeans(masked)
keep <- probe_fail <= MAX_PROBE_FAIL & startsWith(probe_ids, "cg")
log <- list(n_raw_probes = n_raw_probes, n_masked_values = sum(masked),
            n_probes_failed = sum(probe_fail > MAX_PROBE_FAIL),
            n_probes_non_cg = sum(!startsWith(probe_ids, "cg")),
            sample_frac_masked = colMeans(masked))
rm(masked); invisible(gc())
betas <- betas[keep, , drop = FALSE]
log$n_probes_kept <- nrow(betas)

# ---- 5. CoRSIVtools objects -----------------------------------------------------------------------
se <- from_matrix(betas, meta, sample_col = "sample_id")
rm(betas); invisible(gc())
se <- annotate_probes(se, read_manifest("EPICv1"))
metadata(se)$prepare <- c(list(source = "GEO GSE131433", det_p = DET_P, max_probe_fail = MAX_PROBE_FAIL,
                               csv_md5 = unname(tools::md5sum(csv))), log)
saveRDS(se, out_rds)
message(sprintf("wrote %s: %d CpGs x %d samples", out_rds, nrow(se), ncol(se)))
