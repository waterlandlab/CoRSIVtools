#' PCA of methylation for batch-effect QC
#'
#' PCA (centred, unscaled) on M-values of the `n_top` most variable CpGs. CpGs
#' with a missing value in any sample are dropped (and counted).
#'
#' @param se Annotated SummarizedExperiment (from [annotate_probes()]).
#' @param vars `colData` columns to colour the PC plots by (e.g. group, chip, position, sex, age).
#' @param n_top Number of most variable CpGs used.
#' @param features `"all"` (all array CpGs, default) or `"corsiv"` (CpGs in `corsivs` only).
#' @param n_pcs Number of PCs returned.
#' @param corsivs CoRSIV `GRanges` used when `features = "corsiv"`.
#' @return List: `scores` (`DataFrame`: sample, PCs, `vars`), `var_explained`,
#'   `plots` (per variable: `pc12`, `pc34`), `n_cpgs_used`, `n_dropped_na`, `features`.
#' @export
#' @examples
#' sim <- simulate_corsiv_data()
#' se <- annotate_probes(sim$se, sim$manifest)
#' p <- qc_pca(se, vars = c("group", "age"), n_top = 500)
#' p$plots$group$pc12
qc_pca <- function(se, vars = character(), n_top = 10000, features = c("all", "corsiv"), n_pcs = 10,
                   corsivs = corsiv_regions()) {
  features <- match.arg(features)
  b <- .qc_beta(se)
  .check_vars(se, vars)
  if (features == "corsiv") b <- b[overlapsAny(rowRanges(se), corsivs), , drop = FALSE]
  complete <- rowSums(is.na(b)) == 0
  n_dropped <- sum(!complete)
  M <- .to_m(b[complete, , drop = FALSE])
  if (nrow(M) < 2L) stop("fewer than 2 complete CpGs available for PCA", call. = FALSE)
  rv <- rowSums((M - rowMeans(M))^2)
  top <- order(rv, decreasing = TRUE)[seq_len(min(n_top, nrow(M)))]
  p <- stats::prcomp(t(M[top, , drop = FALSE]), center = TRUE, scale. = FALSE)
  k <- min(n_pcs, ncol(p$x))
  ve <- (p$sdev^2 / sum(p$sdev^2))[seq_len(k)]
  scores <- data.frame(sample = colnames(b), p$x[, seq_len(k), drop = FALSE], row.names = NULL)
  cd <- as.data.frame(colData(se))
  for (v in vars) scores[[v]] <- cd[[v]]

  plot_pair <- function(v, i, j) {
    if (j > k) return(NULL)
    xi <- paste0("PC", i); yj <- paste0("PC", j)
    g <- ggplot2::ggplot(scores, ggplot2::aes(.data[[xi]], .data[[yj]], colour = .data[[v]])) +
      ggplot2::geom_point(size = 2.5, alpha = 0.85) +
      ggplot2::labs(x = sprintf("%s (%.1f%%)", xi, 100 * ve[i]), y = sprintf("%s (%.1f%%)", yj, 100 * ve[j]),
                    title = sprintf("PCA coloured by %s", v)) +
      .theme_corsiv()
    if (!is.numeric(scores[[v]]) && length(unique(scores[[v]])) > 12) g <- g + ggplot2::theme(legend.position = "none")
    g
  }
  plots <- stats::setNames(lapply(vars, function(v) list(pc12 = plot_pair(v, 1, 2), pc34 = plot_pair(v, 3, 4))), vars)
  list(scores = DataFrame(scores, check.names = FALSE), var_explained = ve, plots = plots,
       n_cpgs_used = length(top), n_dropped_na = as.integer(n_dropped), features = features)
}

#' Association of principal components with sample variables (batch QC)
#'
#' For each PC and variable: factor variables are tested with an F-test of
#' `PC ~ variable`, numeric variables with a correlation test.
#'
#' @param pca Output of [qc_pca()].
#' @param se The SummarizedExperiment passed to [qc_pca()].
#' @param vars `colData` columns to test.
#' @return List: `table` (`data.frame`: pc, variable, p, var_explained) and `plot`
#'   (heatmap of -log10 P).
#' @export
#' @examples
#' sim <- simulate_corsiv_data()
#' se <- annotate_probes(sim$se, sim$manifest)
#' p <- qc_pca(se, vars = c("group", "age"), n_top = 500)
#' qc_batch(p, se, c("group", "age"))$table
qc_batch <- function(pca, se, vars) {
  .check_vars(se, vars)
  cd <- as.data.frame(colData(se))[as.character(pca$scores$sample), vars, drop = FALSE]
  tab <- .pc_assoc(as.data.frame(pca$scores), cd, pca$var_explained)
  list(table = tab, plot = .assoc_heatmap(tab, "Association of principal components with sample variables"))
}
