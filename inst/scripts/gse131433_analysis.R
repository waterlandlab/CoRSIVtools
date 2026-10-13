# CoRSIV analysis of GSE131433 (ART vs natural conception, at birth and in adulthood) used by
# the package vignette. Writes the small, aggregated outputs the vignette displays:
#   <out_dir>/gse131433_results.rds   tables and printed summaries (no beta matrices)
#   <out_dir>/figures/*.png           figures
#
# Usage: Rscript inst/scripts/gse131433_analysis.R <prepared_se.rds> <out_dir>
#   <prepared_se.rds> is written by inst/scripts/gse131433_prepare.R.
#   The vignette expects <out_dir> = vignettes/gse131433.
# Model: limma on CoRSIV means, ~ conception + sex + array_row + cell fractions, per timepoint.
# Needs: EpiDISH, BSgenome.Hsapiens.UCSC.hg38. Runtime ~ tens of minutes (100 matched draws).

N_DRAWS   <- 100       # matched background draws (>= 20 needed for empirical P)
SEED      <- 1
CELL_REF  <- "centDHSbloodDMC"
CELL_COVS <- c("B", "NK", "CD4T", "CD8T", "Mono")   # neutrophils + eosinophils are the implicit reference
FIG_DPI   <- 110
MAX_SAMPLE_MISSING <- 0.05   # exclude samples with more masked CpGs than this (qc_samples' absolute limit)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: Rscript gse131433_analysis.R <prepared_se.rds> <out_dir>")
in_rds <- args[1]; out_dir <- args[2]
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages({
  library(CoRSIVtools)
  library(SummarizedExperiment)
})
save_fig <- function(p, name, width = 5.5, height = 5.5) {   # small panels: house style uses large text
  ggplot2::ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = width, height = height,
                  dpi = FIG_DPI, bg = "white")
}
printed <- function(x) paste(utils::capture.output(print(x)), collapse = "\n")
out <- list()

se <- readRDS(in_rds)
out$prepare <- metadata(se)$prepare[setdiff(names(metadata(se)$prepare), "sample_frac_masked")]
out$se <- printed(se)
out$design_table <- as.data.frame(table(timepoint = se$timepoint, conception = se$conception))

# ---- QC --------------------------------------------------------------------------------------------
cells <- estimate_cell_fractions(se, ref = CELL_REF)
out$cells_meta <- metadata(cells)
for (cn in colnames(cells)) colData(se)[[cn]] <- cells[[cn]]

qc <- qc_report(se, vars = c("timepoint", "conception", "sex", "slide", "array_row"),
                technical = c("slide", "array_row"), cells = colnames(cells), group = "conception")
out$qc_summary <- qc$summary
out$qc_batch_top <- head(qc$results$batch$table[order(qc$results$batch$table$p), ], 10)
out$qc_cells <- qc$results$cells$table
out$qc_pca <- list(n_cpgs_used = qc$results$pca$n_cpgs_used, var_explained = qc$results$pca$var_explained)
save_fig(qc$results$samples$plot_density, "qc_density")
save_fig(qc$results$pca$plots$timepoint$pc12, "qc_pca_timepoint")
save_fig(qc$results$batch$plot, "qc_batch", width = 10, height = 6)
save_fig(qc$results$spatc1l$plot_scatter, "qc_spatc1l")
save_fig(qc$results$vtrna2$plot_hist, "qc_vtrna2")

# cell composition is checked within each timepoint: birth vs adult blood differ by design
out$cells_by_timepoint <- lapply(split(seq_len(ncol(se)), se$timepoint), function(i)
  qc_cells(se[, i], cells = CELL_COVS, group = "conception")$table)
save_fig(qc_cells(se[, se$timepoint == "birth"], cells = colnames(cells), group = "conception")$plot_box,
         "qc_cells_birth", width = 10, height = 8)

# qc_samples() also flags samples > 3 MADs above the median missingness. With very low median
# missingness that relative rule flags samples missing well under 1% of CpGs, so exclusion here
# uses the absolute limit only; the relative flags are reported.
st <- as.data.frame(qc$results$samples$table)
out$sample_missing <- list(n_flag_relative = sum(st$flag), max_frac_missing = max(st$frac_missing),
                           median_frac_missing = stats::median(st$frac_missing),
                           flagged_by_timepoint = table(se$timepoint[st$flag]))
excluded <- st$sample[st$frac_missing > MAX_SAMPLE_MISSING]
out$excluded <- excluded
se <- se[, !colnames(se) %in% excluded]
out$row_by_timepoint <- table(array_row = se$array_row, timepoint = se$timepoint)
out$slide_by_timepoint <- c(n_slides = length(unique(se$slide)),
                            n_mixed = sum(tapply(se$timepoint, se$slide, function(x) length(unique(x)) > 1)))
out$n_analysed <- as.data.frame(table(timepoint = se$timepoint, conception = se$conception))

# ---- sets and region summaries (shared by both timepoints) ----------------------------------------
t0 <- Sys.time()
sets <- define_sets(se, background = "matched", n_draws = N_DRAWS, seed = SEED)
out$sets <- printed(sets)
out$coverage <- coverage_report(sets)
out$match <- as.data.frame(matchDiagnostics(sets))
reg <- summarise_regions(se, sets)
out$time_sets <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

# ---- ART vs natural, separately at birth and in adulthood -----------------------------------------
# array row (chip position) is confounded with conception within each timepoint (naturally conceived
# samples were placed mostly on a few rows), so it is adjusted for; slide is not (56-57 levels).
out$row_by_conception <- lapply(split(seq_len(ncol(se)), se$timepoint), function(i)
  table(array_row = se$array_row[i], conception = se$conception[i]))
out$row_conception_p <- vapply(out$row_by_conception, function(t)
  withr::with_seed(SEED, stats::fisher.test(t, simulate.p.value = TRUE, B = 1e5)$p.value), numeric(1))
design <- stats::as.formula(paste("~ conception + sex + array_row +", paste(CELL_COVS, collapse = " + ")))
out$design <- deparse(design)
for (tp in levels(se$timepoint)) {
  r  <- reg[, reg$timepoint == tp]
  res <- test_dm(r, design = design, coef = "conceptionART")
  cmp <- compare_to_background(res)
  tg <- as.data.frame(res[res$set == "target", ])
  tg <- tg[order(tg$P), ]
  out[[tp]] <- list(n = ncol(r), n_target = sum(!is.na(tg$P)), n_fdr = sum(tg$FDR < 0.05, na.rm = TRUE),
                    top = head(tg[, c("region_id", "chr", "start", "end", "n_cpg", "delta_beta", "sd_beta",
                                      "P", "FDR")], 10),
                    enrichment = as.data.frame(cmp$enrichment), variability = as.data.frame(cmp$variability),
                    distributions = as.data.frame(cmp$distributions),
                    target = tg[, c("region_id", "chr", "start", "end", "n_cpg", "delta_beta", "sd_beta",
                                    "P", "FDR")],
                    provenance = metadata(res))
  save_fig(plot_manhattan(res), paste0(tp, "_manhattan"), width = 11, height = 5)
  save_fig(plot_volcano(res), paste0(tp, "_volcano"))
  save_fig(plot_background(cmp, "enrichment"), paste0(tp, "_bg_enrichment"))
  save_fig(plot_background(cmp, "variability"), paste0(tp, "_bg_variability"))
  save_fig(plot_corsiv(r, tg$region_id[1], "conception"), paste0(tp, "_top_corsiv"))
}

out$session <- printed(utils::sessionInfo())
out$run_date <- format(Sys.Date())
saveRDS(out, file.path(out_dir, "gse131433_results.rds"))
message("done: ", out_dir)
