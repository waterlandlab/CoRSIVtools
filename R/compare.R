#' Compare CoRSIV results with the background
#'
#' @param res Output of [test_dm()].
#' @param p Nominal P threshold for enrichment.
#' @param fdr FDR threshold reported alongside.
#' @return Named list of `DataFrame`s: `enrichment`, `per_draw`,
#'   `distributions`, `qq`, `variability`.
#' @export
compare_to_background <- function(res, p = 0.05, fdr = 0.05) {
  res <- as.data.frame(res)
  tg <- res[res$set == "target", , drop = FALSE]
  bg <- res[res$set == "background", , drop = FALSE]
  if (!nrow(bg)) stop("res has no background rows; run define_sets() with a background", call. = FALSE)
  per_draw <- do.call(rbind, lapply(split(bg, bg$draw), function(d) data.frame(
    draw = d$draw[1], n = sum(!is.na(d$P)), n_sig = sum(d$P < p, na.rm = TRUE),
    n_fdr = sum(d$FDR < fdr, na.rm = TRUE), median_sd = stats::median(d$sd_beta, na.rm = TRUE))))
  per_draw$frac_sig <- per_draw$n_sig / per_draw$n
  nd <- nrow(per_draw)
  emp <- function(stat, obs) if (nd >= 20) (1 + sum(stat >= obs)) / (1 + nd) else NA_real_
  t_n <- sum(!is.na(tg$P)); t_sig <- sum(tg$P < p, na.rm = TRUE); t_frac <- t_sig / t_n
  b_n <- sum(!is.na(bg$P)); b_sig <- sum(bg$P < p, na.rm = TRUE)
  ft <- stats::fisher.test(matrix(c(t_sig, t_n - t_sig, b_sig, b_n - b_sig), 2))
  enrichment <- DataFrame(threshold_P = p, target_n = t_n, target_n_sig = t_sig, target_frac_sig = t_frac,
                          bg_mean_frac_sig = mean(per_draw$frac_sig),
                          bg_sd_frac_sig = if (nd > 1) stats::sd(per_draw$frac_sig) else NA_real_,
                          empirical_p = emp(per_draw$frac_sig, t_frac),
                          fisher_or = unname(ft$estimate), fisher_p = ft$p.value,
                          target_n_fdr = sum(tg$FDR < fdr, na.rm = TRUE),
                          bg_mean_n_fdr = mean(per_draw$n_fdr), n_draws = nd)
  probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  qf <- function(v) unname(stats::quantile(v, probs, na.rm = TRUE))
  distributions <- DataFrame(quantile = probs, target_abs_delta_beta = qf(abs(tg$delta_beta)),
                             bg_abs_delta_beta = qf(abs(bg$delta_beta)), target_P = qf(tg$P), bg_P = qf(bg$P))
  t_sd <- stats::median(tg$sd_beta, na.rm = TRUE)
  variability <- DataFrame(target_median_sd = t_sd, bg_mean_median_sd = mean(per_draw$median_sd),
                           empirical_p = emp(per_draw$median_sd, t_sd))
  list(enrichment = enrichment, per_draw = DataFrame(per_draw),
       distributions = distributions,
       qq = DataFrame(rbind(.qq(tg$P, "target"), .qq(bg$P, "background"))),
       variability = variability)
}

.qq <- function(p, set) {
  p <- sort(p[!is.na(p)])
  data.frame(set = rep(set, length(p)), expected = -log10(stats::ppoints(length(p))), observed = -log10(p))
}
