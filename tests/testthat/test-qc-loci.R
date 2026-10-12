spatc1l_se <- function(sign = -1, n = 40, k = 5) {
  z <- withr::with_seed(4, stats::rnorm(n))
  eps <- function(s) withr::with_seed(s, matrix(stats::rnorm(2 * k * n, 0, 0.02), 2 * k, n))
  a <- matrix(rep(0.5 + 0.1 * z, each = k), k)
  b <- matrix(rep(0.5 + sign * 0.1 * z, each = k), k)
  betas <- pmin(pmax(rbind(a, b) + eps(5), 0.01), 0.99)
  pos <- c(seq(46160800, by = 300, length.out = k), seq(46184100, by = 300, length.out = k))
  make_qc_se("chr21", pos, betas)
}

vtrna2_se <- function(bimodal = TRUE, n = 40, k = 6) {
  m <- if (bimodal) c(rep(0.1, 10), rep(0.5, n - 10)) else rep(0.45, n)
  noise <- withr::with_seed(6, matrix(stats::rnorm(k * n, 0, 0.03), k, n))
  betas <- pmin(pmax(matrix(rep(m, each = k), k) + noise, 0.01), 0.99)
  make_qc_se("chr5", seq(136080000, by = 150, length.out = k), betas)
}

test_that("bundled QC loci use CoRSIV 2.0 region IDs", {
  l <- CoRSIVtools:::.qc_loci()
  expect_true(all(c("SPATC1L", "VTRNA2-1") %in% l$locus))
  expect_true(all(l$region_id %in% corsiv_regions(chrX = TRUE)$region_id))
})

test_that("qc_spatc1l passes when the two ends are anti-correlated and flags otherwise", {
  q <- qc_spatc1l(spatc1l_se(-1))
  expect_identical(q$status, "pass")
  expect_lt(q$r_regional, -0.8)
  expect_lt(q$median_between, 0)
  expect_identical(dim(q$corr), c(10L, 10L))
  expect_s3_class(q$plot_corr, "ggplot")
  expect_silent(ggplot2::ggplot_build(q$plot_corr))
  expect_s3_class(q$plot_scatter, "ggplot")
  expect_identical(qc_spatc1l(spatc1l_se(1))$status, "flag")
})

test_that("qc_spatc1l is skipped when the array lacks probes there", {
  se <- make_qc_se("chr21", c(46160800, 46184100), matrix(0.5, 2, 5))
  q <- qc_spatc1l(se)
  expect_identical(q$status, "skipped")
  expect_match(q$message, "fewer than 3")
})

test_that("qc_vtrna2 detects bimodality and assigns modes", {
  q <- qc_vtrna2(vtrna2_se(TRUE))
  expect_identical(q$status, "pass")
  expect_identical(q$n_low, 10L)
  expect_identical(sum(q$table$mode == "high"), 30L)
  expect_gt(q$delta_bic, 0)
  expect_s3_class(q$plot_hist, "ggplot")
  expect_silent(ggplot2::ggplot_build(q$plot_probes))
  expect_identical(qc_vtrna2(vtrna2_se(FALSE))$status, "flag")
})
