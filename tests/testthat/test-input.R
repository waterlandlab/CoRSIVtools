b <- matrix(c(.1, .2, .3, .4, .5, .6), 3, 2,
            dimnames = list(c("p1", "p2", "p3"), c("s1", "s2")))

test_that("from_matrix builds an SE with metadata reordered to beta columns", {
  meta <- data.frame(id = c("s2", "s1"), grp = c("B", "A"))
  se <- from_matrix(b, meta, "id")
  expect_s4_class(se, "SummarizedExperiment")
  expect_identical(colnames(se), c("s1", "s2"))
  expect_identical(se$grp, c("A", "B"))
  expect_identical(assay(se, "beta"), b)
  expect_identical(rowData(se)$probe_id, c("p1", "p2", "p3"))
  expect_identical(metadata(se)$n_input_probes, 3L)
})

test_that("from_matrix accepts a numeric data.frame", {
  se <- from_matrix(as.data.frame(b), data.frame(id = c("s1", "s2")), "id")
  expect_identical(unname(assay(se, "beta")), unname(b))
})

test_that("from_matrix rejects invalid input with actionable messages", {
  meta <- data.frame(id = c("s1", "s2"))
  b_dup <- b; rownames(b_dup) <- c("p1", "p1", "p3")
  expect_error(from_matrix(b_dup, meta, "id"), "duplicated probe IDs")
  expect_error(from_matrix(b * 10, meta, "id"), "\\[0, 1\\]")
  expect_error(from_matrix(b, data.frame(id = c("s1", "s3")), "id"), "one-to-one")
  expect_error(from_matrix(b, data.frame(id = c("s1", "s2", "s9")), "id"), "one-to-one")
  expect_error(from_matrix(b, data.frame(id = c("s1", "s1")), "id"), "duplicated")
  expect_error(from_matrix(b, meta, "nope"), "sample_col")
  df <- as.data.frame(b); df$s2 <- as.character(df$s2)
  expect_error(from_matrix(df, meta, "id"), "non-numeric")
  expect_error(from_matrix(b * NA, meta, "id"), "all NA")
  b_nonames <- b; dimnames(b_nonames) <- NULL
  expect_error(from_matrix(b_nonames, meta, "id"), "rownames")
})

test_that(".region_id and .cpg_id give 1-based UCSC strings", {
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(629002, 635000))
  expect_identical(CoRSIVtools:::.region_id(gr), "chr1:629002-635000")
  expect_identical(CoRSIVtools:::.cpg_id(GenomicRanges::GRanges("chr2", IRanges::IRanges(5, 5))), "chr2:5")
  expect_identical(CoRSIVtools:::.chr_type(c("chr1", "chr22", "chrX", "chrY", "chrM")),
                   c("autosome", "autosome", "chrX", "other", "other"))
})
