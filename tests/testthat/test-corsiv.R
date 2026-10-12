test_that("corsiv_regions returns the bundled CoRSIV 2.0 set with coordinate IDs", {
  all <- corsiv_regions(chrX = TRUE)
  expect_length(all, 6669L)
  expect_false(anyDuplicated(all$region_id) > 0)
  expect_identical(metadata(all)$corsiv_version, "2.0")
  i <- match("chr1:629001-635000", all$region_id)
  expect_identical(all$source_corsivs[i], paste0(
    "black_chr1_38,black_chr1_40,black_chr1_41,multi_racial_chr1_11,multi_racial_chr1_12,",
    "multi_racial_chr1_14,multi_racial_chr1_16,multi_racial_chr1_17,white_chr1_15"))
  auto <- corsiv_regions()
  expect_length(auto, 6669L - 320L)
  expect_false(any(as.character(GenomicRanges::seqnames(auto)) == "chrX"))
  expect_false(is.unsorted(auto))
})

test_that("the bundled file's start column is read as 1-based (whole 100-bp bins)", {
  # The CoRSIV 2.0 file stores 1-based, end-inclusive coordinates (bins xx01-xx00),
  # not standard 0-based BED: every region is a whole number of 100-bp bins.
  all <- corsiv_regions(chrX = TRUE)
  expect_true(all(width(all) %% 100L == 0L))
  expect_true(all(start(all) %% 100L == 1L))
  expect_true("chr13:110812301-110814100" %in% all$region_id)
  orig <- utils::read.delim(system.file("extdata", "CoRSIV2.0_hg38.bed", package = "CoRSIVtools"),
                            header = FALSE, quote = "", comment.char = "")
  expect_setequal(all$region_id, sprintf("%s:%d-%d", orig$V1, orig$V2, orig$V3))
})

test_that("write_bed writes standard 0-based BED for bedtools", {
  all <- corsiv_regions(chrX = TRUE)
  f <- tempfile(fileext = ".bed")
  write_bed(all, f)
  got <- utils::read.delim(f, header = FALSE)
  i <- match("chr13:110812301-110814100", got$V4)
  expect_identical(c(got$V2[i], got$V3[i]), c(110812300L, 110814100L))
  expect_true(all((got$V3 - got$V2) %% 100L == 0L))
})

test_that("write_bed accepts a results table", {
  res <- data.frame(chr = "chr2", start = 101L, end = 300L, id = "chr2:101-300")
  f <- tempfile(); write_bed(res, f)
  expect_identical(readLines(f), "chr2\t100\t300\tchr2:101-300")
})

test_that("unknown versions error", {
  expect_error(corsiv_regions("3.0"), "version")
})
