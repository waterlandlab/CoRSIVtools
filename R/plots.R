.col <- c(hyper = "#e74c3c", hypo = "#3498db", ns = "#bdc3c7", alt = "#7f8c8d")

.theme_corsiv <- function() {
  ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16),
                   axis.title = ggplot2::element_text(face = "bold", size = 15),
                   axis.text = ggplot2::element_text(size = 12),
                   legend.text = ggplot2::element_text(size = 11))
}

.subset_set <- function(res, set, draw = NULL) {
  d <- as.data.frame(res); d <- d[d$set == set & !is.na(d$P), , drop = FALSE]
  if (!nrow(d)) stop("no rows with set == '", set, "'", call. = FALSE)
  if (set == "background") {
    nd <- length(unique(d$draw))
    if (is.null(draw) && nd > 1L) {
      stop("`res` has ", nd, " background draws; pass `draw = <n>` to choose one.",
           call. = FALSE)
    }
    if (!is.null(draw)) {
      d <- d[!is.na(d$draw) & d$draw %in% draw, , drop = FALSE]
      if (!nrow(d)) stop("no background rows with draw == ", paste(draw, collapse = ", "),
                         "; check `draw`.", call. = FALSE)
    }
  }
  d
}

#' Manhattan plot of CoRSIV results
#' @param res Output of [test_dm()].
#' @param set `"target"` or `"background"`.
#' @param fdr FDR threshold; significant points in red, dashed line at the
#'   largest P with FDR < `fdr`.
#' @param draw Background draw to plot (integer); required when `set =
#'   "background"` and several draws are present, ignored for `"target"`.
#' @return `ggplot`.
#' @export
#' @examples
#' sim <- simulate_corsiv_data()
#' se <- annotate_probes(sim$se, sim$manifest)
#' sets <- define_sets(se, corsivs = sim$corsivs, background = sim$controls)
#' res <- test_dm(summarise_regions(se, sets), ~ group, coef = "groupB")
#' plot_manhattan(res)
plot_manhattan <- function(res, set = "target", fdr = 0.05, draw = NULL) {
  d <- .subset_set(res, set, draw)
  lv <- paste0("chr", c(1:22, "X", "Y"))
  d <- d[d$chr %in% lv, , drop = FALSE]
  d$chr <- factor(d$chr, levels = intersect(lv, unique(d$chr)))
  d <- d[order(d$chr, d$start), ]
  len <- tapply(as.numeric(d$end), d$chr, max)
  off <- stats::setNames(c(0, cumsum(len))[seq_along(len)], names(len))
  d$pos <- d$start + off[as.character(d$chr)]
  d$band <- factor(as.integer(d$chr) %% 2)
  d$sig <- d$FDR < fdr
  ax <- tapply(d$pos, d$chr, function(v) mean(range(v)))
  g <- ggplot2::ggplot(d, ggplot2::aes(.data$pos, -log10(.data$P))) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$band), size = 1.2, alpha = 0.6, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c(`0` = .col[["ns"]], `1` = .col[["alt"]])) +
    ggplot2::geom_point(data = d[d$sig, , drop = FALSE], colour = .col[["hyper"]], size = 2.5) +
    ggplot2::scale_x_continuous(breaks = ax, labels = sub("^chr", "", names(ax))) +
    ggplot2::labs(x = "Chromosome", y = expression(-log[10](P)),
                  title = sprintf("CoRSIV differential methylation (%s)", set)) +
    .theme_corsiv()
  if (any(d$sig)) {
    g <- g + ggplot2::geom_hline(yintercept = -log10(max(d$P[d$sig])), colour = "grey",
                                 linetype = "dashed", linewidth = 1, alpha = 0.7)
  }
  g
}

#' Volcano plot of CoRSIV results
#' @inheritParams plot_manhattan
#' @param delta Optional |delta_beta| reference lines.
#' @return `ggplot`.
#' @export
#' @examples
#' sim <- simulate_corsiv_data()
#' se <- annotate_probes(sim$se, sim$manifest)
#' sets <- define_sets(se, corsivs = sim$corsivs, background = sim$controls)
#' res <- test_dm(summarise_regions(se, sets), ~ group, coef = "groupB")
#' plot_volcano(res)
plot_volcano <- function(res, set = "target", fdr = 0.05, delta = NULL, draw = NULL) {
  d <- .subset_set(res, set, draw)
  d$cls <- ifelse(d$FDR < fdr & d$delta_beta > 0, "hyper",
                  ifelse(d$FDR < fdr & d$delta_beta < 0, "hypo", "ns"))
  d$cls <- factor(d$cls, levels = c("hyper", "hypo", "ns"))
  g <- ggplot2::ggplot(d, ggplot2::aes(.data$delta_beta, -log10(.data$P), colour = .data$cls)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$cls, alpha = .data$cls)) +
    ggplot2::scale_colour_manual(values = .col[c("hyper", "hypo", "ns")], drop = FALSE,
                                 labels = c("Hyper (FDR)", "Hypo (FDR)", "NS"), name = NULL) +
    ggplot2::scale_size_manual(values = c(hyper = 2.5, hypo = 2.5, ns = 1), guide = "none") +
    ggplot2::scale_alpha_manual(values = c(hyper = 0.85, hypo = 0.85, ns = 0.45), guide = "none") +
    ggplot2::geom_vline(xintercept = 0, colour = "black", linewidth = 0.5, alpha = 0.3) +
    ggplot2::labs(x = expression(Delta * beta), y = expression(-log[10](P)),
                  title = sprintf("CoRSIV volcano (%s)", set)) +
    .theme_corsiv()
  if (!is.null(delta)) {
    g <- g + ggplot2::geom_vline(xintercept = c(-delta, delta), colour = .col[["hyper"]],
                                 linetype = "dotted", alpha = 0.5)
  }
  g
}

#' Beta values of one CoRSIV by group
#' @param x Region-level SE from [summarise_regions()].
#' @param id Target `region_id`, e.g. `"chr1:629002-635000"`.
#' @param group `colData` column to group by.
#' @return `ggplot`.
#' @export
#' @examples
#' sim <- simulate_corsiv_data()
#' se <- annotate_probes(sim$se, sim$manifest)
#' sets <- define_sets(se, corsivs = sim$corsivs, background = sim$controls)
#' plot_corsiv(summarise_regions(se, sets), sim$dm_ids[1], "group")
plot_corsiv <- function(x, id, group) {
  i <- which(rowRanges(x)$region_id == id & rowRanges(x)$set == "target")
  if (length(i) != 1L) stop("target region '", id, "' not found in x", call. = FALSE)
  d <- data.frame(beta = as.numeric(assay(x, "beta")[i, ]), group = factor(colData(x)[[group]]))
  ggplot2::ggplot(d, ggplot2::aes(.data$group, .data$beta)) +
    ggplot2::geom_boxplot(outlier.shape = NA, fill = "white") +
    ggplot2::geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.85, colour = .col[["hypo"]]) +
    ggplot2::labs(x = group, y = "Mean beta across CoRSIV CpGs", title = id) +
    .theme_corsiv()
}

#' Target vs background distribution across draws
#' @param cmp Output of [compare_to_background()].
#' @param which `"enrichment"` (fraction P < threshold) or `"variability"`
#'   (median SD of beta).
#' @return `ggplot`.
#' @export
#' @examples
#' sim <- simulate_corsiv_data()
#' se <- annotate_probes(sim$se, sim$manifest)
#' sets <- define_sets(se, corsivs = sim$corsivs, background = sim$controls)
#' res <- test_dm(summarise_regions(se, sets), ~ group, coef = "groupB")
#' plot_background(compare_to_background(res), "variability")
plot_background <- function(cmp, which = c("enrichment", "variability")) {
  which <- match.arg(which)
  pd <- as.data.frame(cmp$per_draw)
  if (which == "enrichment") {
    v <- pd$frac_sig; obs <- cmp$enrichment$target_frac_sig
    xl <- sprintf("Fraction with P < %g", cmp$enrichment$threshold_P)
  } else {
    v <- pd$median_sd; obs <- cmp$variability$target_median_sd
    xl <- "Median SD of beta across samples"
  }
  ggplot2::ggplot(data.frame(v = v), ggplot2::aes(.data$v)) +
    ggplot2::geom_histogram(bins = max(5L, min(30L, length(v))), fill = .col[["ns"]], colour = "white") +
    ggplot2::geom_vline(xintercept = obs, colour = .col[["hyper"]], linewidth = 1.2) +
    ggplot2::labs(x = xl, y = "Background draws",
                  title = sprintf("CoRSIVs (red) vs %d background draw(s)", length(v))) +
    .theme_corsiv()
}
