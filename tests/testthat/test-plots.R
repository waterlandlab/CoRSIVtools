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
