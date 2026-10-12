#' Attach probe coordinates and collapse replicate probes by coordinate
#'
#' Probes without coordinates are dropped. Replicate probes at the same
#' chr:pos are averaged (NA-aware). Probes that share a cg base name but map
#' to different positions are kept separate.
#'
#' @param se Output of [from_matrix()].
#' @param manifest Output of [read_manifest()].
#' @return A `RangedSummarizedExperiment`, one row per CpG coordinate.
#' @export
#' @examples
#' sim <- simulate_corsiv_data()
#' se <- annotate_probes(sim$se, sim$manifest)
#' se
annotate_probes <- function(se, manifest) {
  if (!is(se, "SummarizedExperiment")) stop("se must be a SummarizedExperiment", call. = FALSE)
  if (!is(manifest, "GRanges") || is.null(manifest$probe_id)) {
    stop("manifest must be the output of read_manifest()", call. = FALSE)
  }
  m <- match(rownames(se), manifest$probe_id)
  n_miss <- sum(is.na(m))
  if (n_miss == nrow(se)) {
    stop("no probe IDs in se were found in the manifest; IDs must match exactly (e.g. 'cg06185909_TC11')",
         call. = FALSE)
  }
  if (n_miss) {
    message(sprintf("%d of %d probes have no coordinates in the manifest and were dropped",
                    n_miss, nrow(se)))
  }
  keep <- !is.na(m)
  b <- assay(se, "beta")[keep, , drop = FALSE]
  gr <- granges(manifest)[m[keep]]
  ids <- rownames(se)[keep]
  key <- .cpg_id(gr)
  f <- factor(key, levels = unique(key))
  ok <- !is.na(b)
  b0 <- b; b0[!ok] <- 0
  sums <- rowsum(b0, f, reorder = FALSE)
  cnt <- rowsum(ok * 1, f, reorder = FALSE)
  vals <- sums / cnt
  vals[cnt == 0] <- NA_real_
  gr <- gr[!duplicated(key)]
  mcols(gr) <- DataFrame(probe_ids = IRanges::CharacterList(split(ids, f)))
  rownames(vals) <- levels(f)
  gr <- GenomeInfoDb::sortSeqlevels(gr)
  out <- SummarizedExperiment(
    assays = list(beta = vals), rowRanges = gr, colData = colData(se),
    metadata = c(metadata(se), list(manifest = metadata(manifest)$manifest,
                                    n_unmapped = n_miss,
                                    n_replicates_collapsed = sum(duplicated(key))))
  )
  rownames(out) <- levels(f)
  out[order(rowRanges(out)), ]
}
