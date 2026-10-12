#' How many CoRSIVs (and background regions) the array covers
#'
#' @param sets Output of [define_sets()].
#' @return `DataFrame`, one row for the supplied CoRSIVs and one for the
#'   background (pooled across draws).
#' @export
coverage_report <- function(sets) {
  stopifnot(is(sets, "CoRSIVSets"))
  row <- function(label, n, kept) {
    DataFrame(set = label, n_regions = length(n), n_ge1 = sum(n >= 1), n_ge2 = sum(n >= 2),
              n_ge3 = sum(n >= 3), n_ge5 = sum(n >= 5),
              median_cpg = if (any(n > 0)) stats::median(n[n > 0]) else NA_real_,
              max_cpg = max(c(0L, n)), min_cpgs = sets@params$min_cpgs, n_kept = kept)
  }
  rbind(row("corsiv", sets@corsivs$n_cpg, length(sets@target)),
        row("background", sets@background$n_cpg, length(sets@background)))
}
