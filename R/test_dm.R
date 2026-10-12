#' Differential methylation of CoRSIVs (and background) with limma
#'
#' Fits limma on M-values separately for the target set and each background
#' draw, so target moderation and FDR never depend on the background.
#' `delta_beta` is the same design fitted on the beta scale. Rows without
#' residual degrees of freedom, or whose own missing values leave a design
#' column non-estimable, are returned with NA statistics and are excluded
#' from the variance prior and the BH family.
#'
#' @param x Region-level SE from [summarise_regions()] (`level = "corsiv"`)
#'   or annotated SE from [annotate_probes()] (`level = "cpg"`).
#' @param design One-sided formula using `colData(x)` columns.
#' @param coef Name of the design-matrix column to test, e.g. `"groupB"`.
#' @param level `"corsiv"` or `"cpg"`.
#' @param sets [define_sets()] output; required for `level = "cpg"`.
#' @return `DataFrame`; see the spec for the column schema. Provenance in
#'   `metadata()`.
#' @export
#' @examples
#' sim <- simulate_corsiv_data()
#' se <- annotate_probes(sim$se, sim$manifest)
#' sets <- define_sets(se, corsivs = sim$corsivs, background = sim$controls)
#' res <- test_dm(summarise_regions(se, sets), ~ group + age, coef = "groupB")
#' head(res[res$set == "target", ])
test_dm <- function(x, design, coef, level = c("corsiv", "cpg"), sets = NULL) {
  level <- match.arg(level)
  if (level == "corsiv") {
    if (!identical(metadata(x)$level, "region")) {
      stop("level = 'corsiv' needs the output of summarise_regions()", call. = FALSE)
    }
    rr <- rowRanges(x)
    info <- data.frame(id = rr$region_id, region_id = rr$region_id, set = rr$set, draw = rr$draw,
                       chr = as.character(seqnames(rr)), start = start(rr), end = end(rr),
                       n_cpg = rr$n_cpg, source_corsivs = rr$source_corsivs)
    b <- assay(x, "beta")
    set_params <- metadata(x)$sets
  } else {
    if (!is(sets, "CoRSIVSets")) stop("level = 'cpg' needs sets = define_sets(...)", call. = FALSE)
    .check_annotated(x)
    built <- .cpg_rows(x, sets)
    info <- built$info; b <- built$b
    set_params <- c(sets@params, mode = sets@mode)
  }
  cd <- as.data.frame(colData(x))
  vars <- all.vars(design)
  miss <- setdiff(vars, names(cd))
  if (length(miss)) stop("design variables not in colData: ", paste(miss, collapse = ", "), call. = FALSE)
  keep <- stats::complete.cases(cd[, vars, drop = FALSE])
  if (any(!keep)) message(sprintf("%d samples dropped for missing design variables", sum(!keep)))
  cd <- cd[keep, , drop = FALSE]; b <- b[, keep, drop = FALSE]
  mm <- stats::model.matrix(design, data = cd)
  if (!coef %in% colnames(mm)) {
    stop(sprintf("coef '%s' is not a design-matrix column; available: %s", coef,
                 paste(colnames(mm), collapse = ", ")), call. = FALSE)
  }
  q <- qr(mm)
  if (q$rank < ncol(mm)) {
    stop("design matrix is not of full rank; aliased column(s): ",
         paste(colnames(mm)[q$pivot[(q$rank + 1):ncol(mm)]], collapse = ", "), call. = FALSE)
  }
  if (nrow(mm) <= ncol(mm)) {
    stop(sprintf(paste("only %d samples remain with complete design variables, but the design has %d columns,",
                       "so no residual degrees of freedom are left; use fewer covariates or fill missing values"),
                 nrow(mm), ncol(mm)), call. = FALSE)
  }
  grp <- paste(info$set, info$draw)
  blocks <- split(seq_len(nrow(info)), factor(grp, levels = unique(grp)))
  fits <- lapply(blocks, function(i) .fit_block(b[i, , drop = FALSE], mm, coef, info[i, , drop = FALSE]))
  n_nonest <- sum(vapply(fits, `[[`, integer(1), "n_nonest"))
  if (n_nonest > 0L) {
    message(sprintf("%d row(s) returned as NA because the full design was not estimable given their missing values",
                    n_nonest))
  }
  res <- do.call(rbind, lapply(fits, `[[`, "info"))
  rownames(res) <- NULL
  res <- DataFrame(res, check.names = FALSE)
  metadata(res) <- .provenance(metadata(x), set_params, design, coef, level, sum(keep))
  res
}

.fit_block <- function(b, mm, coef, info) {
  stat_cols <- c("delta_beta", "logFC_M", "AveM", "t", "P", "FDR")
  for (cc in stat_cols) info[[cc]] <- NA_real_
  info$sd_beta <- apply(b, 1, stats::sd, na.rm = TRUE)
  ok <- rowSums(!is.na(b)) > ncol(mm)          # rows with residual df > 0; others stay NA
  if (!any(ok)) return(list(info = info, n_nonest = 0L))
  bo <- b[ok, , drop = FALSE]
  bc <- pmin(pmax(bo, 1e-6), 1 - 1e-6)
  M <- log2(bc / (1 - bc))
  # rows whose own NA pattern leaves a design column non-estimable stay NA and are
  # excluded from the variance prior and the BH family
  est <- rowSums(is.na(.lmfit_quiet(M, mm)$coefficients)) == 0L
  if (!any(est)) return(list(info = info, n_nonest = sum(!est)))
  fit <- limma::eBayes(.lmfit_quiet(M[est, , drop = FALSE], mm), trend = TRUE, robust = TRUE)
  tt <- limma::topTable(fit, coef = coef, number = Inf, sort.by = "none")
  fb <- .lmfit_quiet(bo[est, , drop = FALSE], mm)
  idx <- which(ok)[est]
  info$delta_beta[idx] <- unname(fb$coefficients[, coef])
  info$logFC_M[idx] <- tt$logFC; info$AveM[idx] <- tt$AveExpr; info$t[idx] <- tt$t
  info$P[idx] <- tt$P.Value; info$FDR[idx] <- stats::p.adjust(tt$P.Value, "BH")
  list(info = info, n_nonest = sum(!est))
}

.lmfit_quiet <- function(y, mm) {
  withCallingHandlers(limma::lmFit(y, mm), warning = function(w) {
    if (grepl("Partial NA coefficients", conditionMessage(w))) invokeRestart("muffleWarning")
  })
}

.cpg_rows <- function(x, sets) {
  cpg <- rowRanges(x); b <- assay(x, "beta")
  one <- function(regions) {
    h <- findOverlaps(cpg, regions)
    i <- queryHits(h); r <- regions[subjectHits(h)]
    list(b = b[i, , drop = FALSE],
         info = data.frame(id = rownames(x)[i], region_id = r$region_id, set = r$set, draw = r$draw,
                           chr = as.character(seqnames(cpg))[i], start = start(cpg)[i],
                           end = end(cpg)[i], n_cpg = rep(1L, length(i)), source_corsivs = r$source_corsivs))
  }
  t <- one(sets@target); g <- one(sets@background)
  list(b = rbind(t$b, g$b), info = rbind(t$info, g$info))
}
