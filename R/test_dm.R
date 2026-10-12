#' Differential methylation of CoRSIVs (and background) with limma
#'
#' Fits limma on M-values separately for the target set and each background
#' draw, so target moderation and FDR never depend on the background.
#' `delta_beta` is the same design fitted on the beta scale.
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
  grp <- paste(info$set, info$draw)
  blocks <- split(seq_len(nrow(info)), factor(grp, levels = unique(grp)))
  res <- do.call(rbind, lapply(blocks, function(i) .fit_block(b[i, , drop = FALSE], mm, coef, info[i, , drop = FALSE])))
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
  if (!any(ok)) return(info)
  bo <- b[ok, , drop = FALSE]
  bc <- pmin(pmax(bo, 1e-6), 1 - 1e-6)
  M <- log2(bc / (1 - bc))
  fit <- limma::eBayes(limma::lmFit(M, mm), trend = TRUE, robust = TRUE)
  tt <- limma::topTable(fit, coef = coef, number = Inf, sort.by = "none")
  fb <- limma::lmFit(bo, mm)
  info$delta_beta[ok] <- unname(fb$coefficients[, coef])
  info$logFC_M[ok] <- tt$logFC; info$AveM[ok] <- tt$AveExpr; info$t[ok] <- tt$t
  info$P[ok] <- tt$P.Value; info$FDR[ok] <- tt$adj.P.Val
  info
}

.cpg_rows <- function(x, sets) {
  cpg <- rowRanges(x); b <- assay(x, "beta")
  one <- function(regions) {
    h <- findOverlaps(cpg, regions)
    i <- queryHits(h); r <- regions[subjectHits(h)]
    list(b = b[i, , drop = FALSE],
         info = data.frame(id = rownames(x)[i], region_id = r$region_id, set = r$set, draw = r$draw,
                           chr = as.character(seqnames(cpg))[i], start = start(cpg)[i],
                           end = end(cpg)[i], n_cpg = 1L, source_corsivs = r$source_corsivs))
  }
  t <- one(sets@target); g <- one(sets@background)
  list(b = rbind(t$b, g$b), info = rbind(t$info, g$info))
}
