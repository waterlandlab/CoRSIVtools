write_lines <- function(lines) { f <- tempfile(fileext = ".bed"); writeLines(lines, f); f }

test_that("BED manifests are read as 1-bp GRanges at the BED end coordinate", {
  f <- write_lines(c("chr4\t5891975\t5891976\tcg06185909_TC11",
                     "chr8\t118111029\t118111030\tcg18975462_BC11"))
  gr <- read_manifest(f)
  expect_s4_class(gr, "GRanges")
  expect_identical(start(gr), c(5891976L, 118111030L))
  expect_true(all(width(gr) == 1L))
  expect_identical(gr$probe_id, c("cg06185909_TC11", "cg18975462_BC11"))
  expect_identical(metadata(gr)$manifest$source, unname(tools::md5sum(f)))
  expect_identical(metadata(gr)$manifest$build, "hg38")
})

test_that("chromosome names without a chr prefix get one", {  # Review Focus 1
  gr <- read_manifest(write_lines("4\t10\t11\tcgA"))
  expect_identical(as.character(GenomicRanges::seqnames(gr)), "chr4")
})

test_that("extra BED columns are ignored", {
  gr <- read_manifest(write_lines("chr1\t10\t11\tcgA\t0\t+"))
  expect_identical(gr$probe_id, "cgA")
})

test_that("invalid manifests are rejected", {
  expect_error(read_manifest(write_lines("chr1\t10\t12\tcgA")), "1-bp")
  expect_error(read_manifest(write_lines(c("chr1\t10\t11\tcgA", "chr1\t20\t21\tcgA"))), "duplicated")
  expect_error(read_manifest("does_not_exist.bed"), "not found")
  expect_error(read_manifest(write_lines("chr1\t10\t11\tcgA"), build = "hg19"), "hg38")
  expect_error(read_manifest("450K"), "hg38 450K manifest")
  f <- tempfile(fileext = ".txt"); writeLines("x", f)
  expect_error(read_manifest(f), "unrecognised")
})

test_that("Illumina CSV manifests are parsed from the [Assay] block", {
  f <- tempfile(fileext = ".csv")
  writeLines(c("Illumina, Inc.,,", "[Heading],,,", "Descriptor File Name,x,,",
               "[Assay],,,",
               "IlmnID,CHR,MAPINFO,Genome_Build",
               "cgA,1,100,38", "cgB,chrX,200,38", "cgC,,0,38",
               "[Controls],,,", "123,Staining,Red,"), f)
  gr <- read_manifest(f)
  expect_identical(gr$probe_id, c("cgA", "cgB"))
  expect_identical(as.character(GenomicRanges::seqnames(gr)), c("chr1", "chrX"))
  expect_identical(start(gr), c(100L, 200L))
  g <- readLines(f); g <- sub(",38$", ",37", g); writeLines(g, f)
  expect_error(read_manifest(f), "Genome_Build")
})

test_that("built-in EPICv2 manifest loads hg38 coordinates", {
  skip_if_not_installed("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
  gr <- read_manifest("EPICv2")
  expect_gt(length(gr), 900000)
  i <- match("cg25324105_BC11", gr$probe_id)
  expect_identical(as.character(GenomicRanges::seqnames(gr))[i], "chr19")
  expect_identical(start(gr)[i], 37692358L)
  expect_identical(metadata(gr)$manifest$source, "EPICv2")
})
