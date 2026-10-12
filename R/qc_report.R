#' Run the full QC suite
#'
#' Runs [qc_samples()], [qc_pca()] + [qc_batch()], [qc_cells()],
#' [qc_spatc1l()] and [qc_vtrna2()], and summarises each check as `pass`,
#' `flag` or `skipped`:
#' - **samples**: flag if any sample is flagged.
#' - **batch**: flag if any variable in `vars` other than `group` is
#'   associated with PC1-PC3 at P < 0.001.
#' - **cells**: flag if any cell fraction differs between groups at P < 0.05;
#'   skipped without `group` or cell fractions.
#' - **SPATC1L** / **VTRNA2-1**: as in the individual checks.
#'
#' @param se Annotated SummarizedExperiment (from [annotate_probes()]).
#' @param vars `colData` columns for PCA colouring and batch tests.
#' @param cells Cell-fraction `colData` columns or a `DataFrame` of fractions;
#'   if `NULL` and `estimate_cells = TRUE`, estimated with [estimate_cell_fractions()].
#' @param group `colData` column with the comparison groups (optional).
#' @param pair `colData` column identifying matched pairs (optional).
#' @param estimate_cells Estimate cell fractions when `cells` is `NULL` (needs EpiDISH).
#' @param n_top,features Passed to [qc_pca()].
#' @return List: `summary` (`data.frame`: check, status, metric) and `results`
#'   (the output of each check).
#' @export
#' @examples
#' sim <- simulate_corsiv_data()
#' se <- annotate_probes(sim$se, sim$manifest)
#' qc_report(se, vars = "group", n_top = 500, estimate_cells = FALSE)$summary
qc_report <- function(se, vars = character(), cells = NULL, group = NULL, pair = NULL, estimate_cells = TRUE,
                      n_top = 10000, features = c("all", "corsiv")) {
  res <- list()
  row <- function(check, status, metric) data.frame(check = check, status = status, metric = metric)

  res$samples <- qc_samples(se)
  nf <- sum(res$samples$table$flag)
  s1 <- row("samples", if (nf) "flag" else "pass", sprintf("%d of %d samples flagged", nf, ncol(se)))

  pca_vars <- unique(c(vars, group))
  res$pca <- qc_pca(se, vars = pca_vars, n_top = n_top, features = match.arg(features))
  res$batch <- qc_batch(res$pca, se, pca_vars)
  tech <- res$batch$table[res$batch$table$variable %in% setdiff(vars, group) &
                            res$batch$table$pc %in% c("PC1", "PC2", "PC3") & !is.na(res$batch$table$p), ]
  top <- res$batch$table[which.min(res$batch$table$p), ]
  s2 <- row("batch", if (any(tech$p < 0.001)) "flag" else "pass",
            sprintf("strongest: %s ~ %s (P = %.2g)", top$pc, top$variable, top$p))

  if (is.null(cells) && estimate_cells && .pkg_available("EpiDISH")) {
    cells <- tryCatch(estimate_cell_fractions(se), error = function(e) NULL)
  }
  if (!is.null(group) && !is.null(cells)) {
    res$cells <- qc_cells(se, cells, group, pair = pair, pca = res$pca)
    k <- sum(res$cells$table$p < 0.05)
    s3 <- row("cells", if (k) "flag" else "pass",
              sprintf("%d of %d cell types differ by %s (P < 0.05)", k, nrow(res$cells$table), group))
  } else {
    res$cells <- .skipped("needs group and cell fractions")
    s3 <- row("cells", "skipped", res$cells$message)
  }

  res$spatc1l <- qc_spatc1l(se)
  s4 <- row("SPATC1L", res$spatc1l$status,
            if (res$spatc1l$status == "skipped") res$spatc1l$message
            else sprintf("regional r = %.2f; median between-CoRSIV probe r = %.2f", res$spatc1l$r_regional,
                         res$spatc1l$median_between))
  res$vtrna2 <- qc_vtrna2(se)
  s5 <- row("VTRNA2-1", res$vtrna2$status,
            if (res$vtrna2$status == "skipped") res$vtrna2$message
            else sprintf("%d low / %d high; means %.2f / %.2f", res$vtrna2$n_low, res$vtrna2$n_high,
                         res$vtrna2$mu[1], res$vtrna2$mu[2]))
  list(summary = rbind(s1, s2, s3, s4, s5), results = res)
}
