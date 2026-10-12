tiny <- function() {
  b <- rbind(c(.1, .2, .3, .4, .5), c(.3, .2, NA, .4, .5), c(.5, .2, NA, .4, .5),
             c(.2, .2, .2, .2, .2), c(.4, .4, .4, .4, .4), c(.6, .6, .6, NA, NA))
  dimnames(b) <- list(paste0("p", 1:6), paste0("s", 1:5))
  se <- from_matrix(b, data.frame(id = paste0("s", 1:5)), "id")
  man <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(10, 20, 30, 1010, 1020, 1030), width = 1),
                                probe_id = paste0("p", 1:6))
  a <- annotate_probes(se, man)
  corsiv <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 100))
  user <- GenomicRanges::GRanges("chr1", IRanges::IRanges(5001, 5100))
  list(a = a, corsiv = corsiv, user = user)
}

test_that("mean beta per region; samples with < min_cpgs non-NA CpGs are NA", {
  t <- tiny()
  x <- define_sets(t$a, t$corsiv, background = GenomicRanges::GRanges("chr1", IRanges::IRanges(1001, 1100)),
                   flank = 0, min_cpgs = 3)
  r <- summarise_regions(t$a, x, max_na_frac = 0.5)
  expect_identical(rownames(r), c("chr1:1-100", "chr1:1001-1100|draw1"))
  v <- assay(r)["chr1:1-100", ]
  expect_equal(unname(v[c(1, 2, 4, 5)]), c(0.3, 0.2, 0.4, 0.5))
  expect_true(is.na(v[3]))                                  # only 1 non-NA CpG in s3
  expect_identical(metadata(r)$level, "region")
  expect_identical(rowRanges(r)$set, c("target", "background"))
  expect_identical(metadata(r)$sets$mode, "user")
})

test_that("regions with too many NA samples are dropped and counted", {
  t <- tiny()
  x <- define_sets(t$a, t$corsiv, background = GenomicRanges::GRanges("chr1", IRanges::IRanges(1001, 1100)),
                   flank = 0, min_cpgs = 3)
  expect_message(r <- summarise_regions(t$a, x, max_na_frac = 0.2), "dropped")
  # target has 1/5 = 0.2 NA (kept); background has 2/5 = 0.4 NA (dropped)
  expect_identical(rownames(r), "chr1:1-100")
  expect_identical(metadata(r)$n_dropped_na[["background"]], 1L)
})

test_that("median option and array-mode error", {
  t <- tiny()
  x <- define_sets(t$a, t$corsiv, background = GenomicRanges::GRanges("chr1", IRanges::IRanges(1001, 1100)),
                   flank = 0, min_cpgs = 3)
  r <- summarise_regions(t$a, x, fun = "median", max_na_frac = 1)
  expect_equal(unname(assay(r)["chr1:1-100", 1]), 0.3)
  xa <- define_sets(t$a, t$corsiv, background = "array", flank = 0)
  expect_error(summarise_regions(t$a, xa), "CpG level only")
})
