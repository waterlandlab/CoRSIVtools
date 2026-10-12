#' Simulate a small CoRSIV methylation dataset on real hg38 chr21 sequence
#'
#' Probes sit on real CpG positions in chr21:30,000,001-34,000,000, so
#' CpG-density matching works on genuine sequence. 400 probe clusters
#' (300-600 bp, 3-8 probes) are laid out one per 10-kb bin; `n_corsivs` of
#' them become CoRSIVs and the rest are returned as `controls` (usable as a
#' user background). CoRSIVs get high
#' between-individual variation (SD 0.10) vs other clusters (SD 0.03).
#' `n_dm` CoRSIVs are shifted by `shift` in group B. Singleton probes are added
#' between clusters. One CpG carries a same-position replicate probe and one
#' probe ID shares its cg base name with a probe elsewhere.
#'
#' @param n_samples Even number of samples, split equally into groups A and B.
#' @param n_corsivs Number of the 400 clusters that are CoRSIVs (<= 200).
#' @param n_dm Number of CoRSIVs with a group effect.
#' @param shift Beta-scale group effect.
#' @param seed Random seed.
#' @return `list(se, manifest, corsivs, controls, dm_ids)`.
#' @export
#' @examples
#' sim <- simulate_corsiv_data(n_samples = 10, n_corsivs = 20, n_dm = 5)
#' sim$se
#' sim$dm_ids
simulate_corsiv_data <- function(n_samples = 20, n_corsivs = 40, n_dm = 10, shift = 0.15, seed = 1) {
  stopifnot(n_samples %% 2 == 0, n_dm <= n_corsivs, n_corsivs <= 200)
  withr::with_seed(seed, {
    genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    s0 <- 30000001L; s1 <- 34000000L
    cg <- BiocGenerics::start(Biostrings::matchPattern("CG", genome[["chr21"]][s0:s1])) + s0 - 1L
    n_clust <- 400L
    bin <- floor((s1 - s0 + 1) / n_clust)
    starts <- integer(n_clust); widths <- integer(n_clust); probes <- vector("list", n_clust)
    for (j in seq_len(n_clust)) {
      for (attempt in 1:50) {
        w <- sample(c(300L, 400L, 500L, 600L), 1)
        st <- s0 + (j - 1L) * bin + sample.int(bin - w - 2000L, 1)
        inside <- cg[cg >= st & cg <= st + w - 1L]
        if (length(inside) >= 3L) break
      }
      starts[j] <- st; widths[j] <- w; probes[[j]] <- sort(sample(inside, min(sample(3:8, 1), length(inside))))
    }
    clusters <- GRanges("chr21", IRanges(starts, width = widths))
    corsiv_idx <- sort(sample.int(n_clust, n_corsivs))
    # singletons between clusters (array background); two spare positions for replicate fixtures
    out_cl <- cg[!IRanges::overlapsAny(GRanges("chr21", IRanges(cg, width = 1L)), clusters + 1000L)]
    single <- sort(sample(out_cl, 1001L))
    spare <- single[1001L]; single <- single[-1001L]
    pos <- c(unlist(probes), single)
    clust_of <- c(rep(seq_len(n_clust), lengths(probes)), rep(NA_integer_, length(single)))
    ids <- sprintf("cg%08d_TC11", seq_along(pos))
    # fixture 1: same-position replicate of probe 1 (collapsed by annotate_probes)
    # fixture 2: probe 2's base name at a different position (must stay a separate CpG)
    ids <- c(ids, sub("_TC11$", "_TC12", ids[1]), sub("_TC11$", "_TC21", ids[2]))
    pos <- c(pos, pos[1], spare)
    clust_of <- c(clust_of, clust_of[1], NA_integer_)
    n <- n_samples
    group <- factor(rep(c("A", "B"), each = n / 2), levels = c("A", "B"))
    in_cv <- !is.na(clust_of) & clust_of %in% corsiv_idx
    sd_r <- ifelse(in_cv, 0.10, 0.03)
    mu <- stats::runif(length(pos), 0.2, 0.8)
    reff <- matrix(stats::rnorm(n_clust * n), n_clust, n)
    b <- t(vapply(seq_along(pos), function(i) {
      r <- if (is.na(clust_of[i])) stats::rnorm(n) else reff[clust_of[i], ]
      mu[i] + sd_r[i] * r + stats::rnorm(n, 0, 0.02)
    }, numeric(n)))
    dm_idx <- corsiv_idx[sample.int(length(corsiv_idx), n_dm)]
    dm_rows <- which(clust_of %in% dm_idx)
    b[dm_rows, group == "B"] <- b[dm_rows, group == "B"] + shift
    b <- pmin(pmax(b, 0.001), 0.999)
    sids <- sprintf("S%02d", seq_len(n))
    dimnames(b) <- list(ids, sids)
    meta <- data.frame(sample_id = sids, group = group, age = round(stats::runif(n, 40, 70), 1))
    se <- from_matrix(b, meta, "sample_id")
    manifest <- GRanges("chr21", IRanges(pos, width = 1L), probe_id = ids)
    metadata(manifest)$manifest <- list(source = "simulated", file = NA_character_, build = "hg38",
                                        n = length(ids))
    corsivs <- clusters[corsiv_idx]
    corsivs$region_id <- .region_id(corsivs)
    corsivs$source_corsivs <- sprintf("sim_corsiv_%03d", seq_along(corsivs))
    metadata(corsivs)$corsiv_version <- "sim"
    controls <- clusters[-corsiv_idx]
    list(se = se, manifest = manifest, corsivs = corsivs, controls = controls,
         dm_ids = .region_id(clusters[dm_idx]))
  })
}
