mk <- function() {
  b <- rbind(cg1_TC11 = c(.2, .4), cg1_TC12 = c(.4, NA), cg2_TC11 = c(.5, .5),
             cg2_TC21 = c(.9, .9), cgZ_BC11 = c(.1, .1), ctl_Neg01 = c(.3, .3))
  colnames(b) <- c("s1", "s2")
  se <- from_matrix(b, data.frame(id = c("s1", "s2")), "id")
  man <- GenomicRanges::GRanges(c("chr2", "chr2", "chr1", "chr1", "chr10"),
                                IRanges::IRanges(c(100, 100, 50, 51, 5), width = 1),
                                probe_id = c("cg1_TC11", "cg1_TC12", "cg2_TC11", "cg2_TC21", "cgZ_BC11"))
  S4Vectors::metadata(man)$manifest <- list(source = "test")
  list(se = se, man = man)
}

test_that("replicates at the same coordinate are averaged; same name at different coordinates is not", {
  x <- mk()
  expect_message(a <- annotate_probes(x$se, x$man), "1 of 6 probes have no coordinates")
  expect_s4_class(a, "RangedSummarizedExperiment")
  expect_identical(rownames(a), c("chr1:50", "chr1:51", "chr2:100", "chr10:5"))
  expect_equal(unname(assay(a)["chr2:100", ]), c(.3, .4))
  expect_identical(as.list(rowRanges(a)$probe_ids)[[3]], c("cg1_TC11", "cg1_TC12"))
  expect_equal(unname(assay(a)["chr1:50", ]), c(.5, .5))
  expect_equal(unname(assay(a)["chr1:51", ]), c(.9, .9))
  expect_identical(metadata(a)$n_unmapped, 1L)
  expect_identical(metadata(a)$n_replicates_collapsed, 1L)
  expect_identical(metadata(a)$manifest$source, "test")
})

test_that("a CpG whose replicates are all NA in a sample stays NA", {
  x <- mk(); assay(x$se)["cg1_TC11", 2] <- NA
  a <- suppressMessages(annotate_probes(x$se, x$man))
  expect_true(is.na(assay(a)["chr2:100", 2]))
})

test_that("no overlap between se and manifest IDs is an error", {
  x <- mk(); x$man$probe_id <- paste0("other", seq_along(x$man))
  expect_error(annotate_probes(x$se, x$man), "no probe IDs")
})
