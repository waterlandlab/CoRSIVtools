#' Summarise CpG methylation per region (CoRSIV and background)
#'
#' @param se Output of [annotate_probes()].
#' @param sets Output of [define_sets()] (not `background = "array"`).
#' @param fun `"mean"` (default) or `"median"` across the region's CpGs.
#' @param max_na_frac Drop regions with more than this fraction of NA samples.
#' @return Region-level `RangedSummarizedExperiment`.
#' @export
#' @examples
#' sim <- simulate_corsiv_data()
#' se <- annotate_probes(sim$se, sim$manifest)
#' sets <- define_sets(se, corsivs = sim$corsivs, background = sim$controls)
#' reg <- summarise_regions(se, sets)
#' reg
summarise_regions <- function(se, sets, fun = c("mean", "median"), max_na_frac = 0.2) {
  fun <- match.arg(fun)
  .check_annotated(se)
  stopifnot(is(sets, "CoRSIVSets"))
  if (sets@mode == "array") {
    stop("background = 'array' is CpG level only; use test_dm(level = 'cpg', sets = sets)", call. = FALSE)
  }
  min_cpgs <- sets@params$min_cpgs
  bg <- sets@background; mcols(bg) <- mcols(bg)[, .std_cols]
  regions <- suppressWarnings(c(sets@target, bg))
  grp <- paste(regions$set, regions$draw)
  blocks <- split(seq_along(regions), factor(grp, levels = unique(grp)))
  b <- assay(se, "beta"); cpg <- rowRanges(se)
  vals <- do.call(rbind, lapply(blocks, function(i) .summarise_block(regions[i], cpg, b, fun, min_cpgs)))
  regions <- regions[unlist(blocks, use.names = FALSE)]
  na_frac <- rowMeans(is.na(vals))
  keep <- na_frac <= max_na_frac
  dropped <- vapply(c("target", "background"), function(s) sum(!keep & regions$set == s), integer(1))
  if (any(dropped)) {
    message(sprintf("dropped %d target and %d background regions with > %g NA samples",
                    dropped[["target"]], dropped[["background"]], max_na_frac))
  }
  regions <- regions[keep]; vals <- vals[keep, , drop = FALSE]
  rn <- ifelse(regions$set == "target", regions$region_id, paste0(regions$region_id, "|draw", regions$draw))
  rownames(vals) <- rn
  colnames(vals) <- colnames(se)
  out <- SummarizedExperiment(
    assays = list(beta = vals), rowRanges = regions, colData = colData(se),
    metadata = c(metadata(se), list(sets = c(sets@params, mode = sets@mode),
                                    summarise = list(fun = fun, max_na_frac = max_na_frac),
                                    level = "region", n_dropped_na = dropped))
  )
  rownames(out) <- rn
  out
}

.summarise_block <- function(reg, cpg, b, fun, min_cpgs) {
  out <- matrix(NA_real_, length(reg), ncol(b))
  h <- findOverlaps(reg, cpg)
  if (!length(h)) return(out)
  q <- queryHits(h)
  sub <- b[subjectHits(h), , drop = FALSE]
  ok <- !is.na(sub)
  cnt <- rowsum(ok * 1, q, reorder = TRUE)
  if (fun == "mean") {
    sub0 <- sub; sub0[!ok] <- 0
    v <- rowsum(sub0, q, reorder = TRUE) / cnt
  } else {
    v <- do.call(rbind, lapply(split(seq_along(q), q), function(i)
      apply(sub[i, , drop = FALSE], 2, stats::median, na.rm = TRUE)))
  }
  v[cnt < min_cpgs] <- NA_real_
  out[as.integer(rownames(cnt)), ] <- v
  out
}
