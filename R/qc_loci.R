.qc_loci <- function() {
  utils::read.csv(system.file("extdata", "qc_loci.csv", package = "CoRSIVtools", mustWork = TRUE),
                  stringsAsFactors = FALSE, check.names = FALSE)
}

.skipped <- function(msg) list(status = "skipped", message = msg)

#' QC: anti-correlated CoRSIVs at the two ends of SPATC1L
#'
#' Methylation of the CoRSIVs at the 5' (TSS) and 3' ends of SPATC1L
#' (`chr21:46184001-46185600` and `chr21:46160701-46162600`, hg38) is
#' anti-correlated across individuals. The check passes when the regional
#' means are negatively correlated and the median probe-level correlation
#' between the two CoRSIVs is negative. Skipped if either CoRSIV has fewer than
#' 3 CpGs on the array.
#'
#' @param se Annotated SummarizedExperiment (from [annotate_probes()]).
#' @return List: `status` (`"pass"`, `"flag"` or `"skipped"`), `r_regional`,
#'   `median_between`, `n_cpg`, `corr` (probe correlation matrix), `plot_corr`,
#'   `plot_scatter`.
#' @export
#' @examples
#' \dontrun{
#' q <- qc_spatc1l(se); q$status; q$plot_corr
#' }
qc_spatc1l <- function(se) {
  b <- .qc_beta(se)
  loci <- .qc_loci(); loci <- loci[loci$locus == "SPATC1L", ]
  three <- loci$region_id[grepl("^3'", loci$role)]; five <- loci$region_id[grepl("^5'", loci$role)]
  r3 <- .region_rows(se, three); r5 <- .region_rows(se, five)
  if (length(r3) < 3L || length(r5) < 3L) {
    return(.skipped(sprintf("fewer than 3 array CpGs in a SPATC1L CoRSIV (3' end: %d, 5' end: %d)",
                            length(r3), length(r5))))
  }
  rows <- c(r3, r5)
  lab <- sprintf("%s:%d", as.character(seqnames(rowRanges(se)))[rows], start(rowRanges(se))[rows])
  C <- stats::cor(t(b[rows, , drop = FALSE]), use = "pairwise.complete.obs")
  dimnames(C) <- list(lab, lab)
  between <- C[seq_along(r3), length(r3) + seq_along(r5)]
  m3 <- colMeans(b[r3, , drop = FALSE], na.rm = TRUE); m5 <- colMeans(b[r5, , drop = FALSE], na.rm = TRUE)
  r_reg <- stats::cor(m3, m5, use = "complete.obs")
  med_between <- stats::median(between, na.rm = TRUE)
  pass <- isTRUE(r_reg < 0 && med_between < 0)

  long <- data.frame(x = factor(rep(lab, times = length(lab)), levels = lab),
                     y = factor(rep(lab, each = length(lab)), levels = rev(lab)), r = as.vector(C))
  pc <- ggplot2::ggplot(long, ggplot2::aes(.data$x, .data$y, fill = .data$r)) +
    ggplot2::geom_tile(colour = "white", linewidth = .lw(0.5)) +
    ggplot2::scale_fill_gradient2(low = "#3498db", mid = "white", high = "#e74c3c", limits = c(-1, 1), guide = .colourbar(),
                                  name = "Pearson r") +
    ggplot2::geom_vline(xintercept = length(r3) + 0.5, linewidth = .lw(0.8)) +
    ggplot2::geom_hline(yintercept = length(r5) + 0.5, linewidth = .lw(0.8)) +
    ggplot2::labs(x = NULL, y = NULL, title = "SPATC1L end CoRSIVs",
                  subtitle = sprintf("3' %s (%d CpGs)\n5' %s (%d CpGs)\nmedian between-CoRSIV r = %.2f",
                                     three, length(r3), five, length(r5), med_between)) +
    .theme_corsiv() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, size = .small_text(7), vjust = 0.5),
                   axis.text.y = ggplot2::element_text(size = .small_text(7)))
  ps <- ggplot2::ggplot(data.frame(m3 = m3, m5 = m5), ggplot2::aes(.data$m3, .data$m5)) +
    ggplot2::geom_point(size = .pt(2.5), alpha = 0.85, colour = "#3498db") +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#333", linewidth = .lw(0.6)) +
    ggplot2::labs(x = expression("Mean" ~ beta * ", 3' CoRSIV"), y = expression("Mean" ~ beta * ", 5' CoRSIV"),
                  title = sprintf("SPATC1L: r = %.2f", r_reg),
                  subtitle = sprintf("Expected < 0\n3' %s\n5' %s", three, five)) +
    .theme_corsiv()
  list(status = if (pass) "pass" else "flag", r_regional = r_reg, median_between = med_between,
       n_cpg = c(three_prime = length(r3), five_prime = length(r5)), corr = C, plot_corr = pc, plot_scatter = ps)
}

# Two-component Gaussian mixture by EM (1-D); deterministic start at min/max.
.gmm2 <- function(x, iter = 500L, tol = 1e-8) {
  km <- stats::kmeans(x, centers = matrix(range(x), ncol = 1))
  mu <- as.numeric(km$centers); sd <- pmax(tapply(x, km$cluster, stats::sd), 0.01); pi <- as.numeric(table(km$cluster)) / length(x)
  sd[is.na(sd)] <- 0.01
  ll_old <- -Inf
  for (i in seq_len(iter)) {
    d <- cbind(pi[1] * stats::dnorm(x, mu[1], sd[1]), pi[2] * stats::dnorm(x, mu[2], sd[2]))
    tot <- rowSums(d); ll <- sum(log(tot))
    z <- d / tot
    nk <- colSums(z); pi <- nk / length(x)
    mu <- colSums(z * x) / nk
    sd <- pmax(sqrt(colSums(z * (x - rep(mu, each = length(x)))^2) / nk), 0.01)
    if (abs(ll - ll_old) < tol) break
    ll_old <- ll
  }
  list(mu = mu, sd = sd, pi = pi, loglik = ll, z = z)
}

#' QC: bimodal methylation at the VTRNA2-1 CoRSIV
#'
#' VTRNA2-1 (`chr5:136079901-136081000`, hg38) is a polymorphically imprinted
#' locus: most people are about half-methylated, while a minority are nearly
#' unmethylated. The check fits one- and two-component Gaussian models to the
#' per-sample mean beta. It passes when two components fit better by BIC and
#' their means differ by more than 0.2. Skipped if fewer than 3 CpGs are on the array.
#'
#' @param se Annotated SummarizedExperiment (from [annotate_probes()]).
#' @return List: `status`, `n_low`, `n_high`, `mu` (component means),
#'   `delta_bic` (BIC1 - BIC2; > 0 favours two modes), `table` (sample,
#'   mean_beta, mode), `plot_hist`, `plot_probes`.
#' @export
#' @examples
#' \dontrun{
#' q <- qc_vtrna2(se); q$status; q$plot_hist
#' }
qc_vtrna2 <- function(se) {
  b <- .qc_beta(se)
  reg <- .qc_loci()$region_id[.qc_loci()$locus == "VTRNA2-1"]
  rows <- .region_rows(se, reg)
  if (length(rows) < 3L) return(.skipped(sprintf("fewer than 3 array CpGs in the VTRNA2-1 CoRSIV (%d)", length(rows))))
  m <- colMeans(b[rows, , drop = FALSE], na.rm = TRUE)
  ok <- !is.na(m); x <- m[ok]; n <- length(x)
  ll1 <- sum(stats::dnorm(x, mean(x), max(stats::sd(x), 0.01), log = TRUE))
  g <- .gmm2(x)
  delta_bic <- (-2 * ll1 + 2 * log(n)) - (-2 * g$loglik + 5 * log(n))
  lo <- which.min(g$mu)
  mode <- rep(NA_character_, length(m))
  mode[ok] <- ifelse(max.col(g$z) == lo, "low", "high")
  pass <- isTRUE(delta_bic > 0 && abs(diff(g$mu)) > 0.2)
  tab <- DataFrame(sample = colnames(b), mean_beta = unname(m), mode = mode)

  ph <- ggplot2::ggplot(as.data.frame(tab[ok, ]), ggplot2::aes(.data$mean_beta, fill = .data$mode)) +
    ggplot2::geom_histogram(binwidth = 0.02, boundary = 0, colour = "white", linewidth = .lw(0.5)) +
    ggplot2::scale_fill_manual(values = c(low = "#3498db", high = "#e74c3c")) +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::labs(x = expression("Mean" ~ beta * ", VTRNA2-1 CoRSIV"), y = "Samples", title = "VTRNA2-1",
                  subtitle = sprintf("%d low, %d high (means %.2f / %.2f)\n%s, %d CpGs",
                                     sum(mode == "low", na.rm = TRUE), sum(mode == "high", na.rm = TRUE),
                                     min(g$mu), max(g$mu), reg, length(rows))) +
    .theme_corsiv()
  long <- data.frame(cpg = factor(rep(start(rowRanges(se))[rows], times = ncol(b))),
                     beta = as.vector(b[rows, , drop = FALSE]), mode = rep(mode, each = length(rows)))
  long <- long[!is.na(long$beta) & !is.na(long$mode), ]
  pp <- ggplot2::ggplot(long, ggplot2::aes(.data$cpg, .data$beta, colour = .data$mode)) +
    ggplot2::geom_jitter(width = 0.2, height = 0, size = .pt(1.2), alpha = 0.7) +
    ggplot2::scale_colour_manual(values = c(low = "#3498db", high = "#e74c3c")) +
    ggplot2::labs(x = "CpG position (chr5)", y = expression(beta), title = "VTRNA2-1 CpGs by mode") +
    .theme_corsiv() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, size = .small_text(7), vjust = 0.5))
  list(status = if (pass) "pass" else "flag", n_low = sum(mode == "low", na.rm = TRUE),
       n_high = sum(mode == "high", na.rm = TRUE), mu = sort(g$mu), delta_bic = delta_bic,
       table = tab, plot_hist = ph, plot_probes = pp)
}
