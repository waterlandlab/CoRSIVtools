setup_sim <- function(n_dm = 10, n_corsivs = 40, background = "array") {
  s <- simulate_corsiv_data(n_dm = n_dm, n_corsivs = n_corsivs)
  a <- suppressMessages(annotate_probes(s$se, s$manifest))
  list(s = s, a = a, x = define_sets(a, s$corsivs, background = background))
}
user_bg <- function(r) r$s$controls

test_that("CoRSIV-level test recovers simulated DM regions with correct effect size", {
  r <- setup_sim()
  x <- define_sets(r$a, r$s$corsivs, background = user_bg(r))
  reg <- summarise_regions(r$a, x, max_na_frac = 1)
  res <- test_dm(reg, design = ~ group + age, coef = "groupB")
  expect_s4_class(res, "DataFrame")
  expect_identical(names(res), c("id", "region_id", "set", "draw", "chr", "start", "end", "n_cpg",
                                 "source_corsivs", "delta_beta", "logFC_M", "AveM", "t", "P", "FDR", "sd_beta"))
  tg <- res[res$set == "target", ]
  hit <- tg$id[tg$FDR < 0.05]
  expect_gte(sum(r$s$dm_ids %in% hit), 8L)
  expect_lt(abs(mean(tg$delta_beta[tg$id %in% r$s$dm_ids]) - 0.15), 0.03)
  expect_identical(tg$id, tg$region_id)
})

test_that("under the null, P < 0.05 occurs at about the nominal rate", {
  r <- setup_sim(n_dm = 0, n_corsivs = 200)
  x <- define_sets(r$a, r$s$corsivs, background = user_bg(r))
  reg <- summarise_regions(r$a, x, max_na_frac = 1)
  tg <- test_dm(reg, ~ group, "groupB"); tg <- tg[tg$set == "target", ]
  ci <- stats::qbinom(c(0.005, 0.995), nrow(tg), 0.05) / nrow(tg)
  expect_gte(mean(tg$P < 0.05), ci[1]); expect_lte(mean(tg$P < 0.05), ci[2])
})

test_that("target statistics are identical whatever the background", {
  r <- setup_sim()
  xa <- define_sets(r$a, r$s$corsivs, background = "array")
  xu <- define_sets(r$a, r$s$corsivs, background = user_bg(r))
  ra <- test_dm(r$a, ~ group, "groupB", level = "cpg", sets = xa)
  ru <- test_dm(r$a, ~ group, "groupB", level = "cpg", sets = xu)
  ta <- ra[ra$set == "target", ]; tu <- ru[ru$set == "target", ]
  expect_identical(ta$id, tu$id)
  expect_equal(ta$P, tu$P); expect_equal(ta$FDR, tu$FDR)
  expect_true(all(ta$region_id %in% targetRegions(xa)$region_id))
  expect_true(all(ta$n_cpg == 1L))
})

test_that("bad designs give clear errors; NA covariates drop samples", {
  r <- setup_sim()
  expect_error(test_dm(r$a, ~ group, "groupB"), "summarise_regions")
  expect_error(test_dm(r$a, ~ group, "groupB", level = "cpg"), "sets")
  expect_error(test_dm(r$a, ~ group, "HCCCase", level = "cpg", sets = r$x), "groupB")
  expect_error(test_dm(r$a, ~ nope, "x", level = "cpg", sets = r$x), "nope")
  a2 <- r$a; a2$dup <- a2$group                                            # Review Focus 3
  expect_error(test_dm(a2, ~ group + dup, "groupB", level = "cpg", sets = r$x), "full rank")
  a3 <- r$a; a3$age[1:2] <- NA
  expect_message(res <- test_dm(a3, ~ group + age, "groupB", level = "cpg", sets = r$x), "2 samples dropped")
  expect_identical(metadata(res)$n_samples, 18L)
})

test_that("an all-NA row returns NA statistics instead of erroring", {  # Review Focus 2
  r <- setup_sim()
  a <- r$a
  i <- which(IRanges::overlapsAny(rowRanges(a), targetRegions(r$x)))[1]
  assay(a)[i, ] <- NA
  res <- test_dm(a, ~ group, "groupB", level = "cpg", sets = r$x)
  row <- res[res$id == rownames(a)[i], ]
  expect_identical(nrow(row), 1L)
  expect_true(is.na(row$P)); expect_true(is.na(row$delta_beta))
  expect_false(all(is.na(res$P)))
})

test_that("a row whose own NA pattern makes the design non-estimable is NA and outside the BH family", {
  r <- setup_sim()
  a <- r$a
  tg_rows <- which(IRanges::overlapsAny(rowRanges(a), targetRegions(r$x)))
  i <- tg_rows[1]
  assay(a)[i, a$group == "B"] <- NA                    # 10 non-NA samples, but groupB not estimable
  expect_no_warning(
    expect_message(res <- test_dm(a, ~ group, "groupB", level = "cpg", sets = r$x),
                   "1 row\\(s\\) returned as NA")
  )
  row <- res[res$id == rownames(a)[i], ]
  expect_true(all(is.na(unlist(as.data.frame(row)[c("delta_beta", "logFC_M", "AveM", "t", "P", "FDR")]))))
  tg <- res[res$set == "target", ]
  ok <- !is.na(tg$P)
  expect_identical(sum(!ok), 1L)
  expect_equal(tg$FDR[ok], stats::p.adjust(tg$P[ok], "BH"))
  expect_false(anyNA(res$P[res$set == "background"]))
})

test_that("too few complete samples for the design is an error, not all-NA output", {
  r <- setup_sim()
  a <- r$a; a$age[-c(1, 2, 11)] <- NA                 # 3 samples left (A, A, B), 3 design columns
  expect_error(suppressMessages(test_dm(a, ~ group + age, "groupB", level = "cpg", sets = r$x)),
               "3 samples remain")
})

test_that("provenance metadata records every setting", {
  r <- setup_sim()
  res <- test_dm(r$a, ~ group, "groupB", level = "cpg", sets = r$x)
  md <- metadata(res)
  expect_true(all(c("CoRSIVtools_version", "corsiv_version", "manifest", "background_mode", "n_draws",
                    "seed", "flank", "min_cpgs", "fun", "max_na_frac", "design", "coef", "level",
                    "n_samples", "R_version", "Bioc_version", "timestamp") %in% names(md)))
  expect_identical(md$background_mode, "array")
  expect_identical(md$design, "~group")
  expect_identical(md$manifest, "simulated")
  expect_identical(md$corsiv_version, "sim")
})
