test_that("qc_report runs every check and summarises statuses", {
  s <- sim_with_chip(shift = 0.1)
  a <- s$a
  a$Neu <- withr::with_seed(3, stats::runif(ncol(a), 0.5, 0.7))
  r <- qc_report(a, vars = c("chip", "group"), cells = "Neu", group = "group", n_top = 300)
  expect_identical(r$summary$check, c("samples", "batch", "cells", "SPATC1L", "VTRNA2-1"))
  expect_true(all(r$summary$status %in% c("pass", "flag", "skipped")))
  expect_identical(r$summary$status[r$summary$check == "batch"], "flag")      # injected chip effect
  expect_identical(r$summary$status[r$summary$check == "SPATC1L"], "skipped") # simulation is chr21:30-34 Mb
  expect_identical(r$summary$status[r$summary$check == "VTRNA2-1"], "skipped")
  expect_named(r$results, c("samples", "pca", "batch", "cells", "spatc1l", "vtrna2"))
})

test_that("qc_report skips cell QC without a group", {
  s <- sim_with_chip()
  r <- qc_report(s$a, vars = "chip", n_top = 300, estimate_cells = FALSE)
  expect_identical(r$summary$status[r$summary$check == "cells"], "skipped")
})
