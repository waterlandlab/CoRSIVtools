sim_sets <- function(n_draws = 3, seed = 1, ...) {
  s <- simulate_corsiv_data()
  a <- suppressMessages(annotate_probes(s$se, s$manifest))
  # The default simulation has CoRSIVs with CpG-island-like density and no feasible
  # control window on the simulated 4 Mb, so the ">5% unmatched" warning is expected here.
  x <- withCallingHandlers(
    define_sets(a, s$corsivs, background = "matched", n_draws = n_draws, seed = seed, ...),
    warning = function(w) if (grepl(">5% of targets unmatched", conditionMessage(w), fixed = TRUE))
      invokeRestart("muffleWarning"))
  list(s = s, a = a, x = x)
}

test_that(".count_in_windows counts positions in closed intervals", {
  p <- list(chr1 = c(10L, 20L, 30L), chr2 = 5L)
  expect_identical(CoRSIVtools:::.count_in_windows(c("chr1", "chr1", "chr2", "chr3"),
                                                   c(10L, 11L, 1L, 1L), c(20L, 29L, 5L, 9L), p),
                   c(2L, 1L, 1L, 0L))
})

test_that("every matched control satisfies all four criteria and exclusions", {
  r <- sim_sets()
  bg <- backgroundRegions(r$x); tg <- targetRegions(r$x)
  expect_gt(length(bg), 0L)
  t <- tg[match(bg$matched_to, tg$region_id)]
  expect_identical(width(bg), width(t))                                     # 1 width
  k <- t$n_cpg
  cnt <- GenomicRanges::countOverlaps(bg, rowRanges(r$a))
  expect_identical(cnt, bg$n_cpg)
  expect_true(all(ifelse(k < 10, cnt == k, abs(cnt - k) <= 1)))            # 2 probe count
  expect_identical(CoRSIVtools:::.chr_type(as.character(GenomicRanges::seqnames(bg))),
                   CoRSIVtools:::.chr_type(as.character(GenomicRanges::seqnames(t))))  # 3 chr type
  d_bg <- CoRSIVtools:::.cpg_density(bg); d_t <- CoRSIVtools:::.cpg_density(t)
  expect_equal(bg$cpg_density, d_bg)
  expect_true(all(abs(d_bg - d_t) <= pmax(0.2 * d_t, 0.25)))               # 4 CpG density
  expect_false(any(IRanges::overlapsAny(bg, r$s$corsivs + 1000)))
  for (d in unique(bg$draw)) {
    g <- bg[bg$draw == d]
    expect_identical(sum(GenomicRanges::countOverlaps(g, g)), length(g))   # no within-draw overlap
  }
  expect_true(all(bg$set == "background"))
  expect_true(all(is.na(bg$source_corsivs)))
})

test_that("same seed reproduces draws; different seed changes them", {
  a <- sim_sets(n_draws = 2, seed = 7); b <- sim_sets(n_draws = 2, seed = 7); c <- sim_sets(n_draws = 2, seed = 8)
  expect_identical(backgroundRegions(a$x)$region_id, backgroundRegions(b$x)$region_id)
  expect_false(identical(backgroundRegions(a$x)$region_id, backgroundRegions(c$x)$region_id))
})

test_that("diagnostics report per-draw matching", {
  r <- sim_sets(n_draws = 3)
  dg <- matchDiagnostics(r$x)
  expect_identical(dg$draw, 1:3)
  expect_identical(dg$n_matched + dg$n_unmatched, rep(length(targetRegions(r$x)), 3L))
  expect_identical(as.integer(table(backgroundRegions(r$x)$draw)), dg$n_matched)
})

test_that("a chromosome type with no background probes is unmatched with a warning, not an error", {  # Review Focus 4
  r <- sim_sets(n_draws = 1)
  tg <- targetRegions(r$x)
  fake <- GenomicRanges::GRanges("chrX", IRanges::IRanges(1e6, width = 400))
  fake$region_id <- "chrX:1000000-1000399"; fake$source_corsivs <- NA_character_
  fake$n_cpg <- 3L; fake$set <- "target"; fake$draw <- NA_integer_
  tg2 <- suppressWarnings(c(tg, fake))
  excl <- GenomicRanges::reduce(GenomicRanges::granges(r$s$corsivs) + 1000)
  expect_warning(expect_warning(m <- CoRSIVtools:::.match_background(tg2, rowRanges(r$a), excl, 1, 1, 50),
                                "no background probes and will be unmatched"),
                 ">5% of targets unmatched", fixed = TRUE)
  expect_identical(m$diagnostics$n_unmatched >= 1L, TRUE)
  expect_false("chrX:1000000-1000399" %in% m$background$matched_to)
})

test_that("windows past the chromosome end are rejected before getSeq", {  # Review Focus 5
  ok <- CoRSIVtools:::.in_bounds(c("chr21", "chr21"), c(46709900L, 100L), c(46710100L, 500L))
  expect_identical(ok, c(FALSE, TRUE))   # chr21 length is 46,709,983
})
