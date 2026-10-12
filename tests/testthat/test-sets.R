sim_annot <- function(...) {
  s <- simulate_corsiv_data(...)
  s$a <- suppressMessages(annotate_probes(s$se, s$manifest))
  s
}

test_that("targets are CoRSIVs with >= min_cpgs CpGs; user background must not overlap", {
  s <- sim_annot()
  user <- GenomicRanges::GRanges("chr21", IRanges::IRanges(c(1e6, 2e6), width = 500))
  x <- define_sets(s$a, s$corsivs, background = user, min_cpgs = 3)
  expect_s4_class(x, "CoRSIVSets")
  expect_identical(x@mode, "user")
  n <- GenomicRanges::countOverlaps(s$corsivs, rowRanges(s$a))
  expect_identical(length(targetRegions(x)), sum(n >= 3))
  expect_identical(names(mcols(targetRegions(x))), CoRSIVtools:::.std_cols)
  expect_true(all(is.na(targetRegions(x)$draw)))
  expect_true(all(targetRegions(x)$set == "target"))
  expect_error(define_sets(s$a, s$corsivs, background = s$corsivs[1:2]), "overlap")
  expect_error(define_sets(s$a, s$corsivs, background = "nope"), "background must be")
})

test_that("array background is every CpG outside CoRSIVs +/- flank", {
  s <- sim_annot()
  x <- define_sets(s$a, s$corsivs, background = "array", flank = 1000)
  bg <- backgroundRegions(x)
  expect_true(all(width(bg) == 1L))
  expect_false(any(IRanges::overlapsAny(bg, s$corsivs + 1000)))
  expect_identical(length(bg), sum(!IRanges::overlapsAny(rowRanges(s$a), s$corsivs + 1000)))
  expect_true(all(bg$draw == 1L))
  expect_identical(bg$region_id, CoRSIVtools:::.cpg_id(bg))
})

test_that("seqlevel-style mismatch gives a helpful error", {
  s <- sim_annot()
  cv <- s$corsivs; GenomeInfoDb::seqlevels(cv) <- sub("chr", "", GenomeInfoDb::seqlevels(cv))
  expect_error(define_sets(s$a, cv, background = "array"), "chr")
})

test_that("show() summarises the object", {
  s <- sim_annot()
  x <- define_sets(s$a, s$corsivs, background = "array")
  expect_output(show(x), "CoRSIVSets")
  expect_output(show(x), "array")
})
