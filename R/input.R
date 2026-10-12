#' Build a probe-level SummarizedExperiment from a beta matrix
#'
#' @param betas Numeric matrix or data.frame of beta values in \[0, 1\];
#'   rownames are probe IDs, colnames are sample IDs.
#' @param meta data.frame of sample metadata with exactly one row per column
#'   of `betas`.
#' @param sample_col Name of the column in `meta` holding the sample IDs
#'   that match `colnames(betas)`.
#' @return A `SummarizedExperiment` with assay `"beta"`.
#' @export
#' @examples
#' b <- matrix(runif(6), 3, 2, dimnames = list(paste0("p", 1:3), c("s1", "s2")))
#' from_matrix(b, data.frame(id = c("s1", "s2"), group = c("A", "B")), "id")
from_matrix <- function(betas, meta, sample_col) {
  if (is.data.frame(betas)) {
    num <- vapply(betas, is.numeric, logical(1))
    if (!all(num)) {
      stop("betas has non-numeric columns: ",
           paste(utils::head(names(betas)[!num], 5), collapse = ", "), call. = FALSE)
    }
    betas <- as.matrix(betas)
  }
  if (!is.matrix(betas) || !is.numeric(betas)) {
    stop("betas must be a numeric matrix or data.frame", call. = FALSE)
  }
  if (is.null(rownames(betas)) || is.null(colnames(betas))) {
    stop("betas must have probe IDs as rownames and sample IDs as colnames", call. = FALSE)
  }
  .stop_dups(rownames(betas), "probe IDs in betas")
  .stop_dups(colnames(betas), "sample IDs in betas")
  if (all(is.na(betas))) stop("betas is all NA", call. = FALSE)
  rng <- range(betas, na.rm = TRUE)
  if (rng[1] < 0 || rng[2] > 1) {
    stop(sprintf("betas must lie in [0, 1]; observed range [%g, %g]. Were M-values supplied?",
                 rng[1], rng[2]), call. = FALSE)
  }
  meta <- as.data.frame(meta, check.names = FALSE)
  if (!is.character(sample_col) || length(sample_col) != 1L || !sample_col %in% names(meta)) {
    stop("sample_col must name one column of meta", call. = FALSE)
  }
  ids <- as.character(meta[[sample_col]])
  .stop_dups(ids, sprintf("sample IDs in meta$%s", sample_col))
  in_b <- setdiff(colnames(betas), ids)
  in_m <- setdiff(ids, colnames(betas))
  if (length(in_b) || length(in_m)) {
    stop(sprintf(paste0("sample IDs must match one-to-one: %d in betas but not meta (e.g. %s); ",
                        "%d in meta but not betas (e.g. %s). Subset meta to the arrayed samples first."),
                 length(in_b), paste(utils::head(in_b, 3), collapse = ", "),
                 length(in_m), paste(utils::head(in_m, 3), collapse = ", ")), call. = FALSE)
  }
  meta <- meta[match(colnames(betas), ids), , drop = FALSE]
  rownames(meta) <- colnames(betas)
  SummarizedExperiment(
    assays = list(beta = betas),
    colData = DataFrame(meta, check.names = FALSE),
    rowData = DataFrame(probe_id = rownames(betas)),
    metadata = list(n_input_probes = nrow(betas))
  )
}
