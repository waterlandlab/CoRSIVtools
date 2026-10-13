#' Per-sample QC: missingness and beta distributions
#'
#' A sample is flagged if its fraction of missing (masked) CpGs is more than
#' 3 MADs above the median across samples, or above 0.05.
#'
#' @param se SummarizedExperiment with a `beta` assay (e.g. from [annotate_probes()]).
#' @param max_cpgs Maximum number of CpGs (randomly sampled, seed 1) used for the density curves.
#' @return List with `table` (`DataFrame`: sample, frac_missing, median_beta, flag),
#'   `plot_missing` and `plot_density` (ggplots).
#' @export
#' @examples
#' sim <- simulate_corsiv_data()
#' q <- qc_samples(annotate_probes(sim$se, sim$manifest))
#' q$table
qc_samples <- function(se, max_cpgs = 20000) {
  b <- .qc_beta(se)
  miss <- colMeans(is.na(b))
  med <- apply(b, 2, stats::median, na.rm = TRUE)
  spread <- stats::mad(miss)
  flag <- unname(miss > 0.05 | (spread > 0 & miss > stats::median(miss) + 3 * spread))
  tab <- DataFrame(sample = colnames(b), frac_missing = unname(miss), median_beta = unname(med), flag = unname(flag))

  pm <- ggplot2::ggplot(as.data.frame(tab), ggplot2::aes(stats::reorder(.data$sample, .data$frac_missing),
                                                        .data$frac_missing, colour = .data$flag)) +
    ggplot2::geom_point(size = .pt(2.5)) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "#bdc3c7", `TRUE` = "#e74c3c"), name = "Flagged") +
    ggplot2::geom_hline(yintercept = 0.05, colour = "grey", linetype = "dashed", linewidth = .lw(0.5)) +
    ggplot2::labs(x = "Sample", y = "Fraction missing", title = "Missingness") +
    .theme_corsiv() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, size = .small_text(7), vjust = 0.5))
  # sample names are unreadable beyond a few dozen samples at the house text size
  if (ncol(b) > 60L) pm <- pm + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                               axis.ticks.x = ggplot2::element_blank())

  rows <- if (nrow(b) > max_cpgs) withr::with_seed(1, sort(sample.int(nrow(b), max_cpgs))) else seq_len(nrow(b))
  dens <- do.call(rbind, lapply(seq_len(ncol(b)), function(j) {
    x <- b[rows, j]; x <- x[!is.na(x)]
    d <- stats::density(x, from = 0, to = 1, n = 256)
    data.frame(sample = colnames(b)[j], beta = d$x, density = d$y, flag = flag[j])
  }))
  pd <- ggplot2::ggplot(dens, ggplot2::aes(.data$beta, .data$density, group = .data$sample, colour = .data$flag)) +
    ggplot2::geom_line(alpha = 0.6, linewidth = .lw(0.5)) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "#bdc3c7", `TRUE` = "#e74c3c"), name = "Flagged") +
    ggplot2::labs(x = expression(beta), y = "Density", title = "Beta distributions") +
    .theme_corsiv()
  list(table = tab, plot_missing = pm, plot_density = pd)
}
