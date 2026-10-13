# Shared helpers for the QC suite.

.qc_beta <- function(se) {
  if (!is(se, "SummarizedExperiment") || !"beta" %in% assayNames(se)) {
    stop("se must be a SummarizedExperiment with a 'beta' assay (e.g. from annotate_probes())", call. = FALSE)
  }
  assay(se, "beta")
}

.to_m <- function(b) {
  bc <- pmin(pmax(b, 1e-6), 1 - 1e-6)
  log2(bc / (1 - bc))
}

.check_vars <- function(se, vars) {
  miss <- setdiff(vars, names(colData(se)))
  if (length(miss)) stop("variables not in colData(se): ", paste(miss, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

# Association of each PC with each variable: factor -> F-test of PC ~ var; numeric -> correlation test.
.pc_assoc <- function(scores, vars_df, var_explained) {
  pcs <- grep("^PC[0-9]+$", names(scores), value = TRUE)
  out <- list()
  for (v in names(vars_df)) {
    x <- vars_df[[v]]
    for (k in pcs) {
      y <- scores[[k]]
      ok <- !is.na(x) & !is.na(y)
      p <- NA_real_
      if (is.numeric(x)) {
        if (sum(ok) > 3 && stats::sd(x[ok]) > 0) p <- stats::cor.test(x[ok], y[ok])$p.value
      } else {
        f <- droplevels(factor(x[ok]))
        if (nlevels(f) >= 2 && nlevels(f) < sum(ok)) p <- stats::anova(stats::lm(y[ok] ~ f))[["Pr(>F)"]][1]
      }
      out[[length(out) + 1]] <- data.frame(pc = k, variable = v, p = p,
                                           var_explained = var_explained[match(k, pcs)])
    }
  }
  res <- do.call(rbind, out)
  res$pc <- factor(res$pc, levels = pcs)
  res
}

.assoc_heatmap <- function(tab, title) {
  tab$neglog10p <- -log10(pmax(tab$p, 1e-300))
  tab$label <- ifelse(is.na(tab$p), "", ifelse(tab$p < 0.001, "***", ifelse(tab$p < 0.01, "**",
                      ifelse(tab$p < 0.05, "*", ""))))
  pcs <- levels(tab$pc)
  ve <- tab$var_explained[match(pcs, as.character(tab$pc))]
  ggplot2::ggplot(tab, ggplot2::aes(.data$pc, .data$variable, fill = .data$neglog10p)) +
    ggplot2::geom_tile(colour = "white", linewidth = .lw(0.5)) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label), size = .small_text(4), colour = "#333") +
    ggplot2::scale_fill_gradient(low = "white", high = "#e74c3c", na.value = "grey90", guide = .colourbar(),
                                 name = expression(-log[10](P))) +
    ggplot2::scale_x_discrete(labels = sprintf("%s\n%.1f%%", pcs, 100 * ve)) +
    ggplot2::labs(x = NULL, y = NULL, title = title,
                  subtitle = "* P < 0.05, ** P < 0.01, *** P < 0.001; % = variance explained") +
    .theme_corsiv() + ggplot2::theme(axis.text.x = ggplot2::element_text(size = .small_text(7)))
}

# Rows of an annotated SE inside a "chr:start-end" region (1-based, inclusive).
.region_rows <- function(se, region_id) {
  m <- regmatches(region_id, regexec("^(.+):([0-9]+)-([0-9]+)$", region_id))[[1]]
  if (!m[2] %in% GenomeInfoDb::seqlevels(rowRanges(se))) return(integer())
  gr <- GRanges(m[2], IRanges(as.integer(m[3]), as.integer(m[4])))
  which(overlapsAny(rowRanges(se), gr))
}
