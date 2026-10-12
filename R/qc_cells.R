#' Estimate blood cell fractions with EpiDISH
#'
#' Reference CpG IDs (e.g. `cg00000029`) are matched to the SE's probe IDs by
#' the name before the first `_` (EPICv2/EPIC-Xtra IDs carry suffixes such as
#' `_TC11`). A reference CpG is used only when its name maps to exactly one CpG
#' coordinate in `se`; names mapping to several coordinates are excluded as
#' ambiguous. CpGs with any missing value are dropped. Counts are recorded in
#' `metadata()`.
#'
#' @param se Annotated SummarizedExperiment (from [annotate_probes()]).
#' @param ref EpiDISH reference: `"cent12CT"` (12 blood cell types, EPIC) or
#'   `"centDHSbloodDMC"` (7 blood cell types).
#' @return `DataFrame` of cell fractions (rows = samples).
#' @export
#' @examples
#' \dontrun{
#' f <- estimate_cell_fractions(se)
#' }
estimate_cell_fractions <- function(se, ref = c("cent12CT", "centDHSbloodDMC")) {
  ref <- match.arg(ref)
  .check_annotated(se)
  if (!.pkg_available("EpiDISH")) {
    stop("install EpiDISH to estimate cell fractions: BiocManager::install('EpiDISH')", call. = FALSE)
  }
  e <- new.env()
  utils::data(list = paste0(ref, ".m"), package = "EpiDISH", envir = e)
  R <- e[[paste0(ref, ".m")]]
  ids <- as.list(rowRanges(se)$probe_ids)
  long <- data.frame(row = rep(seq_along(ids), lengths(ids)), name = sub("_.*$", "", unlist(ids)))
  long <- unique(long)
  n_rows <- table(long$name)
  usable <- names(n_rows)[n_rows == 1L]
  matched <- intersect(rownames(R), usable)
  ambiguous <- intersect(rownames(R), names(n_rows)[n_rows > 1L])
  B <- assay(se, "beta")[long$row[match(matched, long$name)], , drop = FALSE]
  rownames(B) <- matched
  complete <- rowSums(is.na(B)) == 0
  B <- B[complete, , drop = FALSE]
  if (nrow(B) < 50L) {
    stop(sprintf("only %d reference CpGs available after matching; cannot estimate cell fractions", nrow(B)),
         call. = FALSE)
  }
  est <- EpiDISH::epidish(beta.m = B, ref.m = R[rownames(B), , drop = FALSE], method = "RPC")$estF
  out <- DataFrame(est, check.names = FALSE)
  rownames(out) <- colnames(se)
  metadata(out) <- list(ref = ref, n_ref = nrow(R), n_matched = length(matched), n_ambiguous = length(ambiguous),
                        n_missing = nrow(R) - length(matched) - length(ambiguous),
                        n_dropped_na = sum(!complete), n_used = nrow(B))
  out
}

#' Blood cell composition QC: confounding with the comparison groups
#'
#' Tests each cell fraction between groups: paired t-test when `pair` is given
#' (two groups), otherwise Wilcoxon rank-sum (two groups) or Kruskal-Wallis
#' (more). If `pca` is given, also tests each fraction against the PCs.
#'
#' @param se SummarizedExperiment with sample metadata.
#' @param cells Names of `colData` columns holding cell fractions, or a
#'   `DataFrame`/`data.frame` of fractions with rownames = `colnames(se)`
#'   (e.g. from [estimate_cell_fractions()]).
#' @param group `colData` column with the comparison groups.
#' @param pair Optional `colData` column identifying matched pairs.
#' @param pca Optional output of [qc_pca()].
#' @return List: `table` (cell, test, p, per-group means), `plot_box`, and
#'   `plot_pc` (if `pca` given).
#' @export
#' @examples
#' sim <- simulate_corsiv_data()
#' se <- annotate_probes(sim$se, sim$manifest)
#' se$Neu <- runif(ncol(se), 0.5, 0.7)
#' qc_cells(se, cells = "Neu", group = "group")$table
qc_cells <- function(se, cells, group, pair = NULL, pca = NULL) {
  cd <- as.data.frame(colData(se))
  if (is.character(cells)) {
    .check_vars(se, cells)
    frac <- cd[, cells, drop = FALSE]
  } else {
    frac <- as.data.frame(cells)
    if (!identical(rownames(frac), colnames(se))) {
      stop("rownames of cells must equal colnames(se)", call. = FALSE)
    }
  }
  .check_vars(se, c(group, pair))
  g <- factor(cd[[group]])
  rows <- lapply(names(frac), function(cn) {
    x <- frac[[cn]]
    if (!is.null(pair)) {
      if (nlevels(g) != 2L) stop("a paired test needs exactly two groups", call. = FALSE)
      pr <- cd[[pair]]
      x1 <- x[g == levels(g)[1]][order(pr[g == levels(g)[1]])]
      x2 <- x[g == levels(g)[2]][order(pr[g == levels(g)[2]])]
      if (!identical(sort(pr[g == levels(g)[1]]), sort(pr[g == levels(g)[2]]))) {
        stop("every pair must have exactly one sample in each group", call. = FALSE)
      }
      p <- stats::t.test(x2, x1, paired = TRUE)$p.value; test <- "paired t-test"
    } else if (nlevels(g) == 2L) {
      p <- stats::wilcox.test(x ~ g, exact = FALSE)$p.value; test <- "Wilcoxon rank-sum"
    } else {
      p <- stats::kruskal.test(x ~ g)$p.value; test <- "Kruskal-Wallis"
    }
    means <- tapply(x, g, mean, na.rm = TRUE)
    cbind(data.frame(cell = cn, test = test, p = p),
          as.data.frame(as.list(stats::setNames(as.numeric(means), paste0("mean_", names(means))))))
  })
  tab <- do.call(rbind, rows)
  long <- do.call(rbind, lapply(names(frac), function(cn) data.frame(cell = cn, group = g, fraction = frac[[cn]])))
  long$cell <- factor(long$cell, levels = names(frac))
  pb <- ggplot2::ggplot(long, ggplot2::aes(.data$group, .data$fraction)) +
    ggplot2::geom_boxplot(outlier.shape = NA, fill = "white") +
    ggplot2::geom_jitter(width = 0.15, height = 0, size = 1.5, alpha = 0.7, colour = "#3498db") +
    ggplot2::facet_wrap(~ cell, scales = "free_y") +
    ggplot2::labs(x = group, y = "Estimated fraction", title = "Blood cell composition by group") +
    .theme_corsiv()
  out <- list(table = tab, plot_box = pb)
  if (!is.null(pca)) {
    fr <- frac[as.character(pca$scores$sample), , drop = FALSE]
    out$plot_pc <- .assoc_heatmap(.pc_assoc(as.data.frame(pca$scores), fr, pca$var_explained),
                                  "Association of principal components with cell fractions")
  }
  out
}
