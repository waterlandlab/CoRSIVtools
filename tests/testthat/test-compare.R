fake_res <- function(n_draws = 20) {
  tg <- data.frame(set = "target", draw = NA_integer_, P = c(rep(0.01, 30), rep(0.5, 70)),
                   FDR = c(rep(0.04, 10), rep(0.5, 90)), delta_beta = 0.1, sd_beta = 0.1)
  bg <- do.call(rbind, lapply(seq_len(n_draws), function(d)
    data.frame(set = "background", draw = d, P = c(rep(0.01, 5), rep(0.5, 95)),
               FDR = 0.9, delta_beta = 0.01, sd_beta = 0.03)))
  S4Vectors::DataFrame(rbind(tg, bg))
}

test_that("enrichment, empirical p and Fisher test are computed correctly", {
  cmp <- compare_to_background(fake_res())
  e <- cmp$enrichment
  expect_identical(e$target_n_sig, 30L)
  expect_equal(e$target_frac_sig, 0.30)
  expect_equal(e$bg_mean_frac_sig, 0.05)
  expect_equal(e$empirical_p, 1 / 21)
  ft <- stats::fisher.test(matrix(c(30, 70, 100, 1900), 2))
  expect_equal(e$fisher_p, ft$p.value)
  expect_identical(e$target_n_fdr, 10L)
  expect_identical(nrow(cmp$per_draw), 20L)
  expect_equal(cmp$variability$target_median_sd, 0.1)
  expect_equal(cmp$variability$empirical_p, 1 / 21)
  expect_identical(sort(unique(cmp$qq$set)), c("background", "target"))
  expect_identical(cmp$distributions$quantile, c(0.05, 0.25, 0.5, 0.75, 0.95))
})

test_that("empirical p is NA with fewer than 20 draws; missing background errors", {
  expect_true(is.na(compare_to_background(fake_res(1))$enrichment$empirical_p))
  r <- fake_res(); r <- r[r$set == "target", ]
  expect_error(compare_to_background(r), "no background rows")
})
