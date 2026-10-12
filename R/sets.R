.std_cols <- c("region_id", "source_corsivs", "n_cpg", "set", "draw")

#' Target CoRSIVs and background regions for a CoRSIV analysis
#'
#' Created by [define_sets()]. Use [targetRegions()], [backgroundRegions()]
#' and [matchDiagnostics()] to access the contents.
#' @slot target,background,corsivs `GRanges`.
#' @slot mode `"matched"`, `"array"` or `"user"`.
#' @slot params list of settings.
#' @slot diagnostics `DataFrame` of per-draw matching results.
#' @export
setClass("CoRSIVSets", representation(target = "GRanges", background = "GRanges",
                                      corsivs = "GRanges", mode = "character",
                                      params = "list", diagnostics = "DataFrame"))

#' @rdname CoRSIVSets-class
#' @param x A `CoRSIVSets` object.
#' @export
targetRegions <- function(x) x@target
#' @rdname CoRSIVSets-class
#' @export
backgroundRegions <- function(x) x@background
#' @rdname CoRSIVSets-class
#' @export
matchDiagnostics <- function(x) x@diagnostics

#' @rdname CoRSIVSets-class
#' @param object A `CoRSIVSets` object.
#' @export
setMethod("show", "CoRSIVSets", function(object) {
  cat("CoRSIVSets\n")
  cat(sprintf("  CoRSIV set: %s; %d supplied, %d targets with >= %d CpGs\n",
              object@params$corsiv_version, length(object@corsivs), length(object@target),
              object@params$min_cpgs))
  cat(sprintf("  background: %s; %d regions across %d draw(s)\n", object@mode,
              length(object@background), length(unique(object@background$draw))))
})

#' Define target CoRSIVs and the background to compare them with
#'
#' @param se Output of [annotate_probes()].
#' @param corsivs `GRanges` of CoRSIVs; default [corsiv_regions()].
#' @param background `"matched"` (random control regions matched on width,
#'   probe count, chromosome type and CpG density), `"array"` (every CpG
#'   outside CoRSIVs +/- `flank`; CpG level only), or a `GRanges` of
#'   user control regions (must not overlap CoRSIVs).
#' @param n_draws Number of matched draws.
#' @param seed Base seed; draw r uses `seed + r`.
#' @param flank bp excluded around every CoRSIV when choosing background.
#' @param min_cpgs Minimum distinct CpGs for a region to be analysed.
#' @param max_tries Maximum candidate proposals per target per draw.
#' @return A [CoRSIVSets-class] object.
#' @export
define_sets <- function(se, corsivs = corsiv_regions(), background = "matched", n_draws = 100,
                        seed = 1, flank = 1000, min_cpgs = 3, max_tries = 1000) {
  .check_annotated(se)
  cpg <- rowRanges(se)
  if (is.null(corsivs$region_id)) corsivs$region_id <- .region_id(corsivs)
  if (is.null(corsivs$source_corsivs)) corsivs$source_corsivs <- NA_character_
  # a no-common-seqlevels warning is superseded by the error below
  corsivs$n_cpg <- suppressWarnings(countOverlaps(corsivs, cpg))
  target <- corsivs[corsivs$n_cpg >= min_cpgs]
  if (!length(target)) {
    stop(sprintf(paste0("no CoRSIV has >= %d CpGs on this array. Check that chromosome names match ",
                        "(e.g. 'chr1' in both the manifest and the CoRSIVs)."), min_cpgs), call. = FALSE)
  }
  target$set <- "target"; target$draw <- NA_integer_
  mcols(target) <- mcols(target)[, .std_cols]
  excl <- reduce(granges(corsivs) + flank)
  diag <- DataFrame()
  if (is(background, "GRanges")) {
    if (any(overlapsAny(background, corsivs))) {
      stop("user background regions overlap CoRSIVs; remove the overlapping regions first", call. = FALSE)
    }
    bg <- granges(background)
    bg$region_id <- .region_id(bg); bg$source_corsivs <- NA_character_
    bg$n_cpg <- countOverlaps(bg, cpg); bg <- bg[bg$n_cpg >= min_cpgs]
    bg$set <- rep("background", length(bg)); bg$draw <- rep(1L, length(bg))
    mode <- "user"
  } else if (identical(background, "array")) {
    bg <- granges(cpg[!overlapsAny(cpg, excl)])
    mcols(bg) <- DataFrame(region_id = .cpg_id(bg), source_corsivs = NA_character_, n_cpg = 1L,
                           set = "background", draw = 1L)
    mode <- "array"
  } else if (identical(background, "matched")) {
    m <- .match_background(target, cpg, excl, n_draws, seed, max_tries)
    bg <- m$background; diag <- m$diagnostics
    mode <- "matched"
  } else {
    stop("background must be 'matched', 'array', or a GRanges of control regions", call. = FALSE)
  }
  new("CoRSIVSets", target = target, background = bg, corsivs = corsivs, mode = mode,
      params = list(corsiv_version = if (is.null(metadata(corsivs)$corsiv_version)) NA_character_
                                     else metadata(corsivs)$corsiv_version,
                    min_cpgs = min_cpgs, n_draws = if (mode == "matched") n_draws else 1L,
                    seed = seed, flank = flank, max_tries = max_tries),
      diagnostics = diag)
}
