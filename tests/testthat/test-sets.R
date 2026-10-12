sim_annot <- function(...) {
  s <- simulate_corsiv_data(...)
  s$a <- suppressMessages(annotate_probes(s$se, s$manifest))
  s
}

test_that("targets are CoRSIVs with >= min_cpgs CpGs; user background must not overlap", {
  s <- sim_annot()
  user <- GenomicRanges::GRanges("chr21", IRanges::IRanges(c(1e6, 2e6), width = 500))
  expect_warning(x <- define_sets(s$a, s$corsivs, background = user, min_cpgs = 3),
                 "none of the 2 user background regions")
  expect_length(backgroundRegions(x), 0L)
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

test_that("a non-empty user background is filtered and standardised", {
  s <- sim_annot()
  x <- define_sets(s$a, s$corsivs, background = s$controls, min_cpgs = 3)
  bg <- backgroundRegions(x)
  expect_identical(names(mcols(bg)), CoRSIVtools:::.std_cols)
  expect_true(length(bg) > 0L)
  expect_true(all(bg$n_cpg >= 3L))
  expect_true(all(bg$draw == 1L))
  expect_true(all(bg$set == "background"))
  expect_identical(length(bg), sum(GenomicRanges::countOverlaps(s$controls, rowRanges(s$a)) >= 3))
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

test_that("array background is restricted to the targets' chromosome types", {
  s <- simulate_corsiv_data()
  b <- SummarizedExperiment::assay(s$se, "beta")
  ex_ids <- c("cgX0001_TC11", "cgX0002_TC11", "cgY0001_TC11", "cgY0002_TC11")
  extra <- matrix(0.5, length(ex_ids), ncol(b), dimnames = list(ex_ids, colnames(b)))
  se <- from_matrix(rbind(b, extra), as.data.frame(SummarizedExperiment::colData(s$se)), "sample_id")
  m0 <- s$manifest
  man <- GenomicRanges::GRanges(c(as.character(GenomicRanges::seqnames(m0)), "chrX", "chrX", "chrY", "chrY"),
                                IRanges::IRanges(c(GenomicRanges::start(m0), 5e6, 6e6, 3e6, 4e6), width = 1L),
                                probe_id = c(m0$probe_id, ex_ids))
  S4Vectors::metadata(man) <- S4Vectors::metadata(m0)
  a <- suppressMessages(annotate_probes(se, man))
  x <- define_sets(a, s$corsivs, background = "array")
  bg <- backgroundRegions(x)
  expect_true(all(as.character(GenomicRanges::seqnames(bg)) == "chr21"))
  expect_identical(length(bg), sum(!IRanges::overlapsAny(rowRanges(a), s$corsivs + 1000) &
                                     as.character(GenomicRanges::seqnames(rowRanges(a))) == "chr21"))
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

test_that("background = 'none' (the default) gives a target-only analysis", {
  s <- simulate_corsiv_data()
  a <- suppressMessages(annotate_probes(s$se, s$manifest))
  x <- define_sets(a, s$corsivs)
  expect_identical(x@mode, "none")
  expect_length(backgroundRegions(x), 0L)
  expect_identical(names(mcols(backgroundRegions(x))), CoRSIVtools:::.std_cols)
  expect_identical(x@params$n_draws, 0L)
  expect_output(show(x), "none")
  cr <- coverage_report(x)
  expect_identical(cr$n_kept, c(length(targetRegions(x)), 0L))
  reg <- summarise_regions(a, x, max_na_frac = 1)
  expect_true(all(rowRanges(reg)$set == "target"))
  res <- test_dm(reg, ~ group, "groupB")
  expect_true(all(res$set == "target"))
  expect_identical(metadata(res)$background_mode, "none")
  cpg <- test_dm(a, ~ group, "groupB", level = "cpg", sets = x)
  expect_true(all(cpg$set == "target"))
  expect_error(compare_to_background(res), "no background rows")
})
