test_that("corsiv_regions returns the bundled CoRSIV 2.0 set with coordinate IDs", {
  all <- corsiv_regions(chrX = TRUE)
  expect_length(all, 6669L)
  expect_false(anyDuplicated(all$region_id) > 0)
  expect_identical(metadata(all)$corsiv_version, "2.0")
  i <- match("chr1:629002-635000", all$region_id)
  expect_identical(all$source_corsivs[i], paste0(
    "black_chr1_38,black_chr1_40,black_chr1_41,multi_racial_chr1_11,multi_racial_chr1_12,",
    "multi_racial_chr1_14,multi_racial_chr1_16,multi_racial_chr1_17,white_chr1_15"))
  auto <- corsiv_regions()
  expect_length(auto, 6669L - 320L)
  expect_false(any(as.character(GenomicRanges::seqnames(auto)) == "chrX"))
  expect_false(is.unsorted(auto))
})

test_that("BED -> region_id -> BED round-trips exactly", {
  all <- corsiv_regions(chrX = TRUE)
  f <- tempfile(fileext = ".bed")
  write_bed(all, f)
  orig <- readLines(system.file("extdata", "CoRSIV2.0_hg38.bed", package = "CoRSIVtools"))
  got <- readLines(f)
  key <- function(x) sub("\t[^\t]*$", "", x)   # chr, start, end only
  expect_setequal(key(got), key(orig))
  expect_true(all(sub("^[^\t]*\t[^\t]*\t[^\t]*\t", "", got) %in% all$region_id))
})

test_that("write_bed accepts a results table", {
  res <- data.frame(chr = "chr2", start = 101L, end = 300L, id = "chr2:101-300")
  f <- tempfile(); write_bed(res, f)
  expect_identical(readLines(f), "chr2\t100\t300\tchr2:101-300")
})

test_that("unknown versions error", {
  expect_error(corsiv_regions("3.0"), "version")
})
