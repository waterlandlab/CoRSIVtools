.col <- c(hyper = "#e74c3c", hypo = "#3498db", ns = "#bdc3c7", alt = "#7f8c8d")

# House plot style (Nature Genetics-like): plain sans-serif text, large enough to stay legible when
# a panel is shrunk into a multi-panel figure; points and lines at twice the usual size. Every plot
# takes text size from .theme_corsiv() and point sizes / line widths from .pt() / .lw().
.style <- list(base_size = 22, point = 2, line = 2)
.pt <- function(x) x * .style$point
.lw <- function(x) x * .style$line
.theme_corsiv <- function() {
  s <- .style$base_size
  ggplot2::theme_classic(base_size = s, base_family = "sans", base_line_size = .lw(0.6)) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = s, colour = "black"),
                   plot.title.position = "plot",
                   plot.subtitle = ggplot2::element_text(size = s * 0.8, colour = "black"),
                   axis.title = ggplot2::element_text(size = s, colour = "black"),
                   axis.text = ggplot2::element_text(size = s * 0.9, colour = "black"),
                   axis.ticks = ggplot2::element_line(colour = "black"),
                   axis.line = ggplot2::element_line(colour = "black"),
                   legend.position = "top", legend.justification = "left",
                   legend.title = ggplot2::element_text(size = s * 0.85),
                   legend.text = ggplot2::element_text(size = s * 0.85),
                   legend.key.spacing = ggplot2::unit(s * 0.3, "pt"),
                   legend.margin = ggplot2::margin(0, 0, 0, 0),
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(size = s * 0.9, face = "bold"))
}
# continuous colour/fill legend sized for the top legend position
.colourbar <- function() {
  ggplot2::guide_colourbar(theme = ggplot2::theme(legend.key.width = ggplot2::unit(.style$base_size * 5, "pt"),
                                                  legend.title.position = "top"))
}
# small text (dense tick labels, cell labels) scaled with the theme: n is the size at base_size 11
.small_text <- function(n) n * .style$base_size / 11

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

# chromosome tick labels; the small chromosomes after chr13 crowd at large text, so label every other one
.chr_labels <- function(chr) {
  lab <- sub("^chr", "", chr)
  crowded <- lab %in% as.character(c(14, 16, 18, 20, 22))
  if (length(chr) > 12L) lab[crowded] <- ""
  lab
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
    ggplot2::geom_point(ggplot2::aes(colour = .data$band), size = .pt(1.2), alpha = 0.6, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c(`0` = .col[["ns"]], `1` = .col[["alt"]])) +
    ggplot2::geom_point(data = d[d$sig, , drop = FALSE], colour = .col[["hyper"]], size = .pt(2.5)) +
    ggplot2::scale_x_continuous(breaks = ax, labels = .chr_labels(names(ax))) +
    ggplot2::labs(x = "Chromosome", y = expression(-log[10](P)),
                  title = sprintf("Differential methylation (%s)", set)) +
    .theme_corsiv()
  if (any(d$sig)) {
    g <- g + ggplot2::geom_hline(yintercept = -log10(max(d$P[d$sig])), colour = "grey",
                                 linetype = "dashed", linewidth = .lw(1), alpha = 0.7)
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
    ggplot2::geom_point(ggplot2::aes(size = .data$cls, alpha = .data$cls), show.legend = TRUE) +
    ggplot2::scale_colour_manual(values = .col[c("hyper", "hypo", "ns")], drop = FALSE,
                                 labels = c("Hyper", "Hypo", "n.s."), name = NULL) +
    ggplot2::scale_size_manual(values = c(hyper = .pt(2.5), hypo = .pt(2.5), ns = .pt(1)), guide = "none") +
    ggplot2::scale_alpha_manual(values = c(hyper = 0.85, hypo = 0.85, ns = 0.45), guide = "none") +
    # constant key glyphs, so classes with no points (e.g. no FDR hits) still show a key
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = .pt(2.5), alpha = 1))) +
    ggplot2::geom_vline(xintercept = 0, colour = "black", linewidth = .lw(0.5), alpha = 0.3) +
    ggplot2::labs(x = expression(Delta * beta), y = expression(-log[10](P)),
                  title = sprintf("CoRSIV volcano (%s)", set),
                  subtitle = sprintf("Coloured: FDR < %g", fdr)) +
    .theme_corsiv()
  if (!is.null(delta)) {
    g <- g + ggplot2::geom_vline(xintercept = c(-delta, delta), colour = .col[["hyper"]],
                                 linetype = "dotted", linewidth = .lw(0.5), alpha = 0.5)
  }
  g
}

#' Beta values of one CoRSIV by group
#' @param x Region-level SE from [summarise_regions()].
#' @param id Target `region_id`, e.g. `"chr1:629001-635000"`.
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
    ggplot2::geom_boxplot(outlier.shape = NA, fill = "white", linewidth = .lw(0.5)) +
    ggplot2::geom_jitter(width = 0.15, height = 0, size = .pt(2), alpha = 0.85, colour = .col[["hypo"]]) +
    ggplot2::labs(x = group, y = expression("Mean" ~ beta), title = id) +
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
    xl <- expression("Median SD of" ~ beta)
  }
  ggplot2::ggplot(data.frame(v = v), ggplot2::aes(.data$v)) +
    ggplot2::geom_histogram(bins = max(5L, min(30L, length(v))), fill = .col[["ns"]], colour = "white", linewidth = .lw(0.5)) +
    ggplot2::geom_vline(xintercept = obs, colour = .col[["hyper"]], linewidth = .lw(1.2)) +
    ggplot2::labs(x = xl, y = "Background draws", title = "CoRSIVs vs background",
                  subtitle = sprintf("Red: CoRSIVs; grey: %d matched draw(s)", length(v))) +
    .theme_corsiv()
}
