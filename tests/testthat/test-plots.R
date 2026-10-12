pl <- local({
  s <- simulate_corsiv_data()
  a <- suppressMessages(annotate_probes(s$se, s$manifest))
  x <- define_sets(a, s$corsivs, background = s$controls)
  reg <- summarise_regions(a, x, max_na_frac = 1)
  res <- test_dm(reg, ~ group, "groupB")
  list(reg = reg, res = res, cmp = compare_to_background(res), dm = s$dm_ids)
})

test_that("all plots return buildable ggplot objects", {
  for (g in list(plot_manhattan(pl$res), plot_volcano(pl$res, delta = 0.05),
                 plot_corsiv(pl$reg, pl$dm[1], "group"),
                 plot_background(pl$cmp, "enrichment"), plot_background(pl$cmp, "variability"))) {
    expect_s3_class(g, "ggplot")
    expect_silent(ggplot2::ggplot_build(g))
  }
})

test_that("plot_corsiv errors on an unknown region", {
  expect_error(plot_corsiv(pl$reg, "chr1:1-2", "group"), "not found")
})

.mk <- function(chr, start, P, fdr, db, set = "target", draw = NA_integer_) {
  data.frame(chr = chr, start = start, end = start + 100L, P = P, FDR = fdr,
             delta_beta = db, set = set, draw = draw)
}

test_that("manhattan base points are ordered along the genome", {
  x <- ggplot2::layer_data(plot_manhattan(pl$res), 1)$x
  expect_false(is.unsorted(x))
})

test_that("manhattan orders chromosomes numerically and labels without chr", {
  d <- .mk(c("chr10", "chr2", "chr10", "chr2"), c(500L, 900L, 100L, 100L),
           c(0.01, 0.2, 0.5, 0.03), c(0.02, 0.4, 0.6, 0.05), c(.1, -.1, .2, -.2))
  g <- plot_manhattan(d)
  ld <- ggplot2::layer_data(g, 1)
  grp <- d[order(factor(d$chr, levels = c("chr2", "chr10")), d$start), ]
  expect_true(all(ld$x[grp$chr == "chr10"] > max(ld$x[grp$chr == "chr2"])))
  b <- ggplot2::ggplot_build(g)$layout$panel_params[[1]]$x
  expect_equal(b$get_labels(), c("2", "10"))
  expect_true(b$breaks[1] < b$breaks[2])
})

test_that("volcano colour counts match FDR and sign", {
  t <- as.data.frame(pl$res); t <- t[t$set == "target" & !is.na(t$P), ]
  exp <- c(sum(t$FDR < 0.05 & t$delta_beta > 0), sum(t$FDR < 0.05 & t$delta_beta < 0),
           sum(!(t$FDR < 0.05 & t$delta_beta != 0)))
  ld <- ggplot2::layer_data(plot_volcano(pl$res), 1)
  got <- c(sum(ld$colour == "#e74c3c"), sum(ld$colour == "#3498db"), sum(ld$colour == "#bdc3c7"))
  expect_equal(got, exp)
  expect_equal(sum(got), nrow(t))
})

test_that("plot_corsiv has one point per sample and the group levels on x", {
  g <- plot_corsiv(pl$reg, pl$dm[1], "group")
  b <- ggplot2::ggplot_build(g)
  expect_equal(nrow(b$data[[2]]), ncol(pl$reg))
  expect_equal(b$layout$panel_params[[1]]$x$get_labels(),
               levels(factor(SummarizedExperiment::colData(pl$reg)$group)))
})

test_that("background plots require and honour a draw choice", {
  bg <- rbind(.mk("chr1", c(1L, 500L), c(.1, .2), c(.3, .4), c(.1, -.1), "background", 1L),
              .mk("chr1", c(1L, 500L, 900L), c(.01, .2, .3), c(.03, .4, .5), c(.2, -.2, .1),
                  "background", 2L))
  expect_error(plot_volcano(bg, set = "background"), "draw")
  expect_error(plot_manhattan(bg, set = "background"), "draw")
  expect_error(plot_volcano(bg, set = "background", draw = 9L), "draw")
  g <- plot_volcano(bg, set = "background", draw = 1L)
  expect_silent(ggplot2::ggplot_build(g))
  expect_equal(nrow(g$data), 2L)
  expect_equal(unique(g$data$draw), 1L)
})
