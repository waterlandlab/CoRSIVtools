.hg38 <- function() BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

.cpg_density <- function(gr) {
  seqs <- Biostrings::getSeq(.hg38(), granges(gr))
  Biostrings::vcountPattern("CG", seqs) / width(gr) * 100
}

.in_bounds <- function(chr, s, e) {
  len <- unname(GenomeInfoDb::seqlengths(.hg38())[chr])
  !is.na(len) & s >= 1L & e <= len
}

.count_in_windows <- function(chr, s, e, pos_by_chr) {
  out <- integer(length(chr))
  for (ch in unique(chr)) {
    p <- pos_by_chr[[ch]]
    if (is.null(p)) next
    i <- which(chr == ch)
    out[i] <- findInterval(e[i], p) - findInterval(s[i] - 1L, p)
  }
  out
}

.match_background <- function(target, cpg, excl, n_draws, seed, max_tries, batch = 20L) {
  target$cpg_density <- .cpg_density(target)
  pos_by_chr <- split(start(cpg), as.character(seqnames(cpg)))
  pos_by_chr <- lapply(pos_by_chr, sort)
  pool <- granges(cpg[!overlapsAny(cpg, excl)])
  pool_type <- .chr_type(as.character(seqnames(pool)))
  no_pool <- !.chr_type(as.character(seqnames(target))) %in% unique(pool_type)
  if (any(no_pool)) {
    warning(sprintf("%d target(s) are on chromosome types with no background probes and will be unmatched",
                    sum(no_pool)), call. = FALSE)
  }
  out <- vector("list", n_draws); diag <- vector("list", n_draws)
  for (r in seq_len(n_draws)) {
    ctrl <- withr::with_seed(seed + r, .match_one_draw(target, pool, pool_type, pos_by_chr, excl,
                                                       max_tries, batch))
    ctrl$draw <- rep(as.integer(r), length(ctrl))
    out[[r]] <- ctrl
    diag[[r]] <- DataFrame(draw = as.integer(r), n_target = length(target), n_matched = length(ctrl),
                           n_unmatched = length(target) - length(ctrl))
  }
  diag <- do.call(rbind, diag)
  bad <- diag$n_unmatched / diag$n_target > 0.05
  if (any(bad)) {
    warning(sprintf("%d of %d draws left >5%% of targets unmatched (max %d unmatched); see matchDiagnostics()",
                    sum(bad), n_draws, max(diag$n_unmatched)), call. = FALSE)
  }
  bg <- unlist(GRangesList(out), use.names = FALSE)
  list(background = bg, diagnostics = diag)
}

.match_one_draw <- function(target, pool, pool_type, pos_by_chr, excl, max_tries, batch) {
  n <- length(target); w <- width(target); k <- target$n_cpg; d <- target$cpg_density
  tchr <- .chr_type(as.character(seqnames(target)))
  pool_chr <- as.character(seqnames(pool)); pool_pos <- start(pool)
  tries <- integer(n); done <- logical(n)
  for (ty in unique(tchr)) if (!any(pool_type == ty)) tries[tchr == ty] <- max_tries
  acc <- GRanges(); acc_t <- integer(); acc_k <- integer(); acc_d <- numeric()
  repeat {
    todo <- which(!done & tries < max_tries)
    if (!length(todo)) break
    nb <- pmin(batch, max_tries - tries[todo])
    ti <- rep(todo, nb)
    tries[todo] <- tries[todo] + nb
    anc <- integer(length(ti))
    for (ty in unique(tchr[ti])) {
      sel <- tchr[ti] == ty; cand <- which(pool_type == ty)
      anc[sel] <- cand[sample.int(length(cand), sum(sel), replace = TRUE)]
    }
    a_chr <- pool_chr[anc]
    s <- pool_pos[anc] - as.integer(floor(stats::runif(length(ti)) * w[ti]))
    e <- s + w[ti] - 1L
    keep <- .in_bounds(a_chr, s, e)
    cnt <- .count_in_windows(a_chr, s, e, pos_by_chr)
    keep <- keep & ifelse(k[ti] < 10L, cnt == k[ti], abs(cnt - k[ti]) <= 1L)
    if (!any(keep)) next
    ti <- ti[keep]; cnt <- cnt[keep]
    cand <- GRanges(a_chr[keep], IRanges(s[keep], e[keep]))
    keep <- !overlapsAny(cand, excl)
    if (length(acc)) keep <- keep & !overlapsAny(cand, acc)
    if (!any(keep)) next
    ti <- ti[keep]; cnt <- cnt[keep]; cand <- cand[keep]
    dens <- .cpg_density(cand)
    keep <- abs(dens - d[ti]) <= pmax(0.2 * d[ti], 0.25)
    keep <- keep & !duplicated(ti) & !done[ti]
    if (!any(keep)) next
    ti <- ti[keep]; cnt <- cnt[keep]; cand <- cand[keep]; dens <- dens[keep]
    h <- findOverlaps(cand, drop.self = TRUE, drop.redundant = TRUE)
    if (length(h)) {
      drop <- unique(pmax(queryHits(h), subjectHits(h)))
      ti <- ti[-drop]; cnt <- cnt[-drop]; cand <- cand[-drop]; dens <- dens[-drop]
    }
    acc <- c(acc, cand); acc_t <- c(acc_t, ti); acc_k <- c(acc_k, cnt); acc_d <- c(acc_d, dens)
    done[ti] <- TRUE
  }
  o <- order(acc_t)
  acc <- acc[o]
  mcols(acc) <- DataFrame(region_id = .region_id(acc), source_corsivs = rep(NA_character_, length(acc)),
                          n_cpg = acc_k[o], set = rep("background", length(acc)),
                          draw = rep(NA_integer_, length(acc)), matched_to = target$region_id[acc_t[o]],
                          cpg_density = acc_d[o], target_density = d[acc_t[o]])
  acc
}
