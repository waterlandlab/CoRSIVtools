# Build a small annotated SE at given coordinates (QC fixtures).
make_qc_se <- function(chr, pos, betas, coldata = NULL, ids = NULL) {
  if (is.null(ids)) ids <- sprintf("cg%08d_TC11", seq_along(pos))
  rownames(betas) <- ids
  if (is.null(colnames(betas))) colnames(betas) <- sprintf("S%02d", seq_len(ncol(betas)))
  if (is.null(coldata)) coldata <- data.frame(sample_id = colnames(betas))
  coldata$sample_id <- colnames(betas)
  se <- from_matrix(betas, coldata, "sample_id")
  man <- GenomicRanges::GRanges(chr, IRanges::IRanges(pos, width = 1), probe_id = ids)
  S4Vectors::metadata(man)$manifest <- list(source = "fixture")
  suppressMessages(annotate_probes(se, man))
}

# Simulation with an injected chip effect on half of the samples.
sim_with_chip <- function(shift = 0.08) {
  s <- simulate_corsiv_data()
  b <- SummarizedExperiment::assay(s$se, "beta")
  chip <- rep(c("chipA", "chipB"), length.out = ncol(b))
  b[, chip == "chipB"] <- pmin(b[, chip == "chipB"] + shift, 0.999)
  s$se <- from_matrix(b, cbind(as.data.frame(SummarizedExperiment::colData(s$se)), chip = chip), "sample_id")
  s$a <- suppressMessages(annotate_probes(s$se, s$manifest))
  s
}
