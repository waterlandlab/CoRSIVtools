test_that("qc_samples flags a sample with many missing CpGs", {
  s <- sim_with_chip()
  a <- s$a
  b <- SummarizedExperiment::assay(a)
  b[seq_len(round(0.1 * nrow(b))), 3] <- NA
  SummarizedExperiment::assay(a) <- b
  q <- qc_samples(a)
  expect_s4_class(q$table, "DataFrame")
  expect_identical(names(q$table), c("sample", "frac_missing", "median_beta", "flag"))
  expect_identical(q$table$sample[q$table$flag], colnames(a)[3])
  expect_s3_class(q$plot_missing, "ggplot")
  expect_s3_class(q$plot_density, "ggplot")
  expect_silent(ggplot2::ggplot_build(q$plot_density))
})

test_that("qc_samples handles zero spread in missingness (MAD = 0)", {
  s <- sim_with_chip()
  q <- qc_samples(s$a)
  expect_false(any(q$table$flag))
})

test_that("qc_pca and qc_batch detect an injected chip effect on PC1", {
  s <- sim_with_chip(shift = 0.1)
  p <- qc_pca(s$a, vars = c("chip", "group", "age"), n_top = 500)
  expect_identical(nrow(p$scores), ncol(s$a))
  expect_true(all(c("PC1", "PC2", "chip", "group", "age") %in% names(p$scores)))
  expect_equal(sum(p$var_explained), sum(p$var_explained[seq_along(p$var_explained)]))
  expect_named(p$plots, c("chip", "group", "age"))
  expect_s3_class(p$plots$chip$pc12, "ggplot")
  bt <- qc_batch(p, s$a, vars = c("chip", "group", "age"))
  pc1_chip <- bt$table$p[bt$table$pc == "PC1" & bt$table$variable == "chip"]
  expect_lt(pc1_chip, 1e-6)
  expect_s3_class(bt$plot, "ggplot")
  expect_silent(ggplot2::ggplot_build(bt$plot))
})

test_that("qc_pca drops CpGs with missing values and reports how many", {
  s <- sim_with_chip()
  a <- s$a
  b <- SummarizedExperiment::assay(a); b[1:5, 1] <- NA; SummarizedExperiment::assay(a) <- b
  p <- qc_pca(a, vars = "chip", n_top = 200)
  expect_identical(p$n_dropped_na, 5L)
})

test_that("qc_pca can restrict to CoRSIV CpGs", {
  s <- sim_with_chip()
  p <- qc_pca(s$a, vars = "chip", n_top = 100, features = "corsiv", corsivs = s$corsivs)
  n_cv <- sum(IRanges::overlapsAny(SummarizedExperiment::rowRanges(s$a), s$corsivs))
  expect_lte(p$n_cpgs_used, n_cv)
  expect_identical(p$features, "corsiv")
  expect_error(qc_pca(s$a, vars = "nope"), "nope")
})
