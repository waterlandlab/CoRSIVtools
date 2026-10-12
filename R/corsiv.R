.corsiv_files <- c("2.0" = "CoRSIV2.0_hg38.bed")

#' CoRSIV regions bundled with the package
#'
#' @param version CoRSIV release. Currently only `"2.0"`.
#' @param chrX Include chrX CoRSIVs? Default `FALSE` (chrX needs
#'   sex-stratified models, not supported in v0.1).
#' @return `GRanges` (1-based) with `region_id` (`chr:start-end`) and
#'   `source_corsivs` (original merged CoRSIV names, not unique).
#' @export
#' @examples
#' cr <- corsiv_regions()
#' head(cr)
corsiv_regions <- function(version = "2.0", chrX = FALSE) {
  if (!version %in% names(.corsiv_files)) {
    stop("unknown CoRSIV version '", version, "'; available: ",
         paste(names(.corsiv_files), collapse = ", "), call. = FALSE)
  }
  path <- system.file("extdata", .corsiv_files[[version]], package = "CoRSIVtools", mustWork = TRUE)
  gr <- .read_region_bed(path)
  if (!chrX) gr <- gr[as.character(seqnames(gr)) != "chrX"]
  metadata(gr)$corsiv_version <- version
  gr
}

.read_region_bed <- function(path) {
  b <- utils::read.delim(path, header = FALSE, quote = "", comment.char = "",
                         colClasses = c("character", "integer", "integer", "character"))
  # The CoRSIV 2.0 file stores 1-based, end-inclusive coordinates (100-bp bins xx01-xx00),
  # not standard 0-based BED, so the start column is used as is.
  gr <- GRanges(b[[1]], IRanges(b[[2]], b[[3]]))
  gr$region_id <- .region_id(gr)
  gr$source_corsivs <- b[[4]]
  gr <- GenomeInfoDb::sortSeqlevels(gr)
  sort(gr)
}

#' Write regions or results as a 0-based BED file for bedtools
#'
#' @param x `GRanges` (name column = `region_id`) or a results table with
#'   columns `chr`, `start`, `end`, `id` (1-based, as returned by [test_dm()]).
#' @param file Output path.
#' @return `file`, invisibly.
#' @export
#' @examples
#' f <- tempfile(fileext = ".bed")
#' write_bed(head(corsiv_regions()), f)
#' readLines(f, 2)
write_bed <- function(x, file) {
  if (is(x, "GRanges")) {
    nm <- if (!is.null(x$region_id)) x$region_id else .region_id(x)
    df <- data.frame(as.character(seqnames(x)), start(x) - 1L, end(x), nm)
  } else {
    x <- as.data.frame(x)
    need <- c("chr", "start", "end", "id")
    if (!all(need %in% names(x))) {
      stop("x needs columns: ", paste(need, collapse = ", "), call. = FALSE)
    }
    df <- data.frame(x$chr, as.integer(x$start) - 1L, as.integer(x$end), x$id)
  }
  utils::write.table(df, file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  invisible(file)
}
