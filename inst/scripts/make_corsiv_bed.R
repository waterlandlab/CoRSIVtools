# Provenance and validation of inst/extdata/CoRSIV2.0_hg38.bed
#
# Upstream generation (Waterland Lab): population-specific CoRSIV 2.0 calls
# (ESS, autosomes + chrX, male and female) at 100-bp resolution, merged with
#   bedtools merge -c 4 -o distinct
# producing CoRSIV2.0_ESS_chrX_male_female_all_cooridnates_100bp_hg38.bedtools_merged_fixed.bed.
# The upstream calling scripts live in https://github.com/waterlandlab/CoRSIV2.0 (scripts/).
#
# This script copies the merged file into the package and validates it.
# Usage: Rscript inst/scripts/make_corsiv_bed.R <path/to/merged.bed>
args <- commandArgs(trailingOnly = TRUE)
src <- args[1]
stopifnot(unname(tools::md5sum(src)) == "ca1b9a6c50e78e9a43a92f2e56400fb2")
b <- utils::read.delim(src, header = FALSE, quote = "", comment.char = "",
                       colClasses = c("character", "integer", "integer", "character"))
gr <- GenomicRanges::GRanges(b[[1]], IRanges::IRanges(b[[2]] + 1L, b[[3]]))
stopifnot(nrow(b) == 6669L,
          !any(GenomicRanges::countOverlaps(gr, gr, maxgap = 0L) > 1L),
          !anyDuplicated(sprintf("%s:%d-%d", b[[1]], b[[2]] + 1L, b[[3]])))
file.copy(src, "inst/extdata/CoRSIV2.0_hg38.bed", overwrite = TRUE)
