.region_id <- function(gr) {
  sprintf("%s:%d-%d", as.character(GenomicRanges::seqnames(gr)),
          BiocGenerics::start(gr), BiocGenerics::end(gr))
}

.cpg_id <- function(gr) {
  sprintf("%s:%d", as.character(GenomicRanges::seqnames(gr)), BiocGenerics::start(gr))
}

.chr_type <- function(chr) {
  ifelse(chr %in% paste0("chr", 1:22), "autosome",
         ifelse(chr == "chrX", "chrX", "other"))
}

.stop_dups <- function(x, what) {
  d <- unique(x[duplicated(x)])
  if (length(d)) {
    stop(sprintf("duplicated %s (%d), e.g. %s", what, length(d),
                 paste(utils::head(d, 3), collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

.check_annotated <- function(se) {
  if (!is(se, "RangedSummarizedExperiment") ||
      !"probe_ids" %in% names(mcols(rowRanges(se)))) {
    stop("se must be the output of annotate_probes()", call. = FALSE)
  }
  invisible(TRUE)
}
