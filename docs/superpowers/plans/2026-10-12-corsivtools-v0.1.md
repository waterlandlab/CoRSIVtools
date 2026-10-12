# CoRSIVtools v0.1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build CoRSIVtools v0.1, an R package that takes a normalised beta matrix, sample metadata and a probe manifest, then:
- restricts the analysis to CoRSIV probes,
- summarises methylation per CoRSIV,
- tests for differential methylation with limma, and
- compares the results with a matched, whole-array or user-supplied background.

**Architecture:** Every step takes and returns a `SummarizedExperiment` (SE), except the final results, which are tidy `S4Vectors::DataFrame`s.
- A `CoRSIVSets` S4 object holds the target CoRSIVs and the background regions.
- Matched background sampling is vectorised rejection sampling: probe counts come from `findInterval` on sorted CpG positions, and CpG density comes from `BSgenome.Hsapiens.UCSC.hg38`.
- limma is fitted separately for each (set, draw) block, so target statistics never depend on the background.

**Tech Stack:** R ≥ 4.5, Bioconductor 3.22.
- Core Bioconductor classes: SummarizedExperiment, GenomicRanges, IRanges, S4Vectors, BiocGenerics, GenomeInfoDb.
- Analysis and genome: limma, Biostrings, BSgenome.Hsapiens.UCSC.hg38.
- Other: ggplot2, withr; roxygen2 and testthat (3rd edition) for development.

**Spec:** `docs/superpowers/specs/2026-10-12-corsivtools-design.md` (approved 2026-10-12). Read it before starting any task.

## Global Constraints

- Package name `CoRSIVtools`; version `0.1.0`; `License: MIT + file LICENSE`; copyright holder `Baylor College of Medicine`.
- Supported genome build: hg38 only.
- The primary key is `region_id` = `sprintf("%s:%d-%d", chr, start_1based, end)`. CpG IDs are `sprintf("%s:%d", chr, pos)`.
- Never truncate or rewrite probe IDs. The full original CoRSIV name is kept in `source_corsivs` and is never used as a key, because it isn't unique.
- Replicate probes are collapsed **by identical coordinate**, never by `cg` base name.
- Defaults: `min_cpgs = 3`, `n_draws = 100`, `seed = 1`, `flank = 1000`, `max_tries = 1000`, `fun = "mean"`, `max_na_frac = 0.2`, `chrX = FALSE`.
- M = log2(β/(1−β)) with β clamped to [1e-6, 1 − 1e-6]; `eBayes(trend = TRUE, robust = TRUE)`; BH FDR within each (set, draw) block.
- Matching rule: exact width; probe count exact if k < 10, else ±1; same chromosome type (autosome / chrX); CpG density (CpGs per 100 bp) within `max(0.2 * d, 0.25)`; no overlap with any CoRSIV ± `flank` or with another control in the same draw.
- No participant data anywhere in the package. Tests use `simulate_corsiv_data()` only.
- Figure style: white background, no top/right spines, bold axis titles of 14–16 pt, palette red `#e74c3c` / blue `#3498db` / grey `#bdc3c7`.
- Error messages are user-facing: name the argument, state the problem, and say how to fix it. Use `call. = FALSE`.
- Package root: `~/code/CoRSIVtools`. Run every command from there.

## Review Focus

These are failure modes the spec implies but no task's main tests exercise. Each has a test in the task named in brackets.

1. **The manifest uses `1`, `2`, … instead of `chr1`, `chr2`, …** Expected: `read_manifest()` adds the `chr` prefix, so overlaps with CoRSIVs aren't silently zero. [Task 2]
2. **A CpG or region row is all NA after filtering.** Expected: `test_dm()` returns that row with `P = NA` and doesn't error. [Task 9]
3. **The design matrix isn't full rank** (e.g. a covariate identical to the group). Expected: a clear error naming the aliased column, not a silent limma `NA` coefficient. [Task 9]
4. **A target's chromosome type has no background probes** (e.g. chrX included but no chrX probes outside the CoRSIVs). Expected: that target is counted as unmatched and `define_sets()` warns but doesn't error. [Task 6]
5. **A candidate window would run past the end of a chromosome.** Expected: the candidate is rejected before `getSeq()`, and nothing errors. [Task 6]

---

## File map

| File | Responsibility | Task |
|---|---|---|
| `DESCRIPTION`, `LICENSE`, `LICENSE.md`, `.Rbuildignore`, `R/CoRSIVtools-package.R` | Package metadata, imports | 1 |
| `R/utils.R` | ID helpers, chromosome type, argument checks | 1 |
| `R/input.R` | `from_matrix()` | 1 |
| `R/manifest.R` | `read_manifest()` | 2 |
| `R/annotate.R` | `annotate_probes()` | 3 |
| `R/corsiv.R`, `inst/extdata/*`, `inst/scripts/make_corsiv_bed.R` | `corsiv_regions()`, `write_bed()`, bundled BED | 4 |
| `R/simulate.R` | `simulate_corsiv_data()` | 5 |
| `R/sets.R` | `CoRSIVSets` class, accessors, `define_sets()` (user and array modes) | 5 |
| `R/background.R` | Matched sampling | 6 |
| `R/coverage.R` | `coverage_report()` | 7 |
| `R/summarise.R` | `summarise_regions()` | 8 |
| `R/test_dm.R`, `R/provenance.R` | `test_dm()` (both levels), provenance | 9 |
| `R/compare.R` | `compare_to_background()` | 10 |
| `R/plots.R` | Plots and theme | 11 |
| `vignettes/*`, `inst/CITATION`, `NEWS.md`, `README.md`, `.github/workflows/*` | Docs, CI, check-clean | 12 |
| — | Real-data smoke run on Addanki (no inference) | 13 |

---

### Task 1: Package scaffold and `from_matrix()`

**Files:**
- Create: `DESCRIPTION`, `LICENSE`, `LICENSE.md`, `.Rbuildignore`, `R/CoRSIVtools-package.R`, `R/utils.R`, `R/input.R`, `tests/testthat.R`, `tests/testthat/test-input.R`
- Modify: `.gitignore` (already exists)

**Interfaces:**
- Produces: `from_matrix(betas, meta, sample_col)` → `SummarizedExperiment`. It has assay `"beta"` (numeric matrix; rownames = probe IDs, colnames = sample IDs), `colData` = the metadata rows reordered to match the columns, `rowData$probe_id`, and `metadata(se)$n_input_probes`.
- Produces (internal, `R/utils.R`): `.region_id(gr)`, `.cpg_id(gr)`, `.chr_type(chr)` → `"autosome"` / `"chrX"` / `"other"`, `.stop_dups(x, what)`, `.check_annotated(se)`.

- [ ] **Step 1: Install the missing development dependencies**

```bash
Rscript -e 'BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg38","BiocStyle","BiocCheck"), update = FALSE, ask = FALSE)'
Rscript -e 'for (p in c("BSgenome.Hsapiens.UCSC.hg38","BiocStyle","BiocCheck")) cat(p, requireNamespace(p, quietly = TRUE), "\n")'
```
Expected: all three print `TRUE`.

- [ ] **Step 2: Write the package metadata files**

`DESCRIPTION`:
```
Package: CoRSIVtools
Type: Package
Title: CoRSIV-Focused DNA Methylation Analysis
Version: 0.1.0
Authors@R: c(
    person("Chathura J.", "Gunasekara", email = "cjgunase@mtu.edu",
           role = c("aut", "cre")),
    person("Baylor College of Medicine", role = "cph"))
Description: Restricts DNA methylation array analysis to correlated regions
    of systemic interindividual variation (CoRSIVs). It maps array probes to
    CoRSIVs for any array given a probe manifest, summarises methylation per
    CoRSIV, tests differential methylation with limma under arbitrary
    covariate designs, and compares CoRSIV results against matched random
    control regions, the rest of the array, or user-supplied control regions.
License: MIT + file LICENSE
Encoding: UTF-8
Roxygen: list(markdown = TRUE)
RoxygenNote: 7.3.2
Depends: R (>= 4.5.0)
Imports:
    BiocGenerics,
    Biostrings,
    BSgenome.Hsapiens.UCSC.hg38,
    GenomeInfoDb,
    GenomicRanges,
    ggplot2,
    IRanges,
    limma,
    methods,
    S4Vectors,
    stats,
    SummarizedExperiment,
    tools,
    utils,
    withr
Suggests:
    BiocStyle,
    IlluminaHumanMethylationEPICanno.ilm10b5.hg38,
    IlluminaHumanMethylationEPICv2anno.20a1.hg38,
    knitr,
    rmarkdown,
    testthat (>= 3.0.0)
Config/testthat/edition: 3
VignetteBuilder: knitr
biocViews: DNAMethylation, MethylationArray, DifferentialMethylation,
    Epigenetics, Software
URL: https://github.com/waterlandlab/CoRSIVtools
BugReports: https://github.com/waterlandlab/CoRSIVtools/issues
```

`LICENSE`:
```
YEAR: 2026
COPYRIGHT HOLDER: Baylor College of Medicine
```

`LICENSE.md`: run `Rscript -e 'usethis::use_mit_license("Baylor College of Medicine")'`. Then check that `DESCRIPTION` still says `License: MIT + file LICENSE` and that `LICENSE` has the two lines above. Restore them if usethis rewrote them differently.

`.Rbuildignore`:
```
^docs$
^\.github$
^LICENSE\.md$
^README\.Rmd$
^.*\.Rproj$
^\.Rproj\.user$
```

`tests/testthat.R`:
```r
library(testthat)
library(CoRSIVtools)
test_check("CoRSIVtools")
```

`R/CoRSIVtools-package.R`:
```r
#' CoRSIVtools: CoRSIV-focused DNA methylation analysis
#'
#' @keywords internal
#' @import methods
#' @import S4Vectors
#' @import IRanges
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @importFrom BiocGenerics start end width
#' @importFrom ggplot2 .data
"_PACKAGE"
```

- [ ] **Step 3: Write the failing tests**

`tests/testthat/test-input.R`:
```r
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
```

- [ ] **Step 4: Run the tests and confirm they fail**

Run: `Rscript -e 'devtools::test(filter = "input")'`
Expected: FAIL with `could not find function "from_matrix"` (or a load error because `R/input.R` doesn't exist yet).

- [ ] **Step 5: Implement**

`R/utils.R`:
```r
.region_id <- function(gr) {
  sprintf("%s:%d-%d", as.character(GenomicRanges::seqnames(gr)),
          BiocGenerics::start(gr), BiocGenerics::end(gr))
}

.cpg_id <- function(gr) {
  sprintf("%s:%d", as.character(GenomicRanges::seqnames(gr)), BiocGenerics::start(gr))
}

.chr_type <- function(chr) {
  ifelse(chr %in% paste0("chr", 1:22), "autosome",
         ifelse(chr == "chrX", "chrX", "other"))
}

.stop_dups <- function(x, what) {
  d <- unique(x[duplicated(x)])
  if (length(d)) {
    stop(sprintf("duplicated %s (%d), e.g. %s", what, length(d),
                 paste(utils::head(d, 3), collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

.check_annotated <- function(se) {
  if (!is(se, "RangedSummarizedExperiment") ||
      !"probe_ids" %in% names(mcols(rowRanges(se)))) {
    stop("se must be the output of annotate_probes()", call. = FALSE)
  }
  invisible(TRUE)
}
```

`R/input.R`:
```r
#' Build a probe-level SummarizedExperiment from a beta matrix
#'
#' @param betas Numeric matrix or data.frame of beta values in \[0, 1\];
#'   rownames are probe IDs, colnames are sample IDs.
#' @param meta data.frame of sample metadata with exactly one row per column
#'   of `betas`.
#' @param sample_col Name of the column in `meta` holding the sample IDs
#'   that match `colnames(betas)`.
#' @return A `SummarizedExperiment` with assay `"beta"`.
#' @export
#' @examples
#' b <- matrix(runif(6), 3, 2, dimnames = list(paste0("p", 1:3), c("s1", "s2")))
#' from_matrix(b, data.frame(id = c("s1", "s2"), group = c("A", "B")), "id")
from_matrix <- function(betas, meta, sample_col) {
  if (is.data.frame(betas)) {
    num <- vapply(betas, is.numeric, logical(1))
    if (!all(num)) {
      stop("betas has non-numeric columns: ",
           paste(utils::head(names(betas)[!num], 5), collapse = ", "), call. = FALSE)
    }
    betas <- as.matrix(betas)
  }
  if (!is.matrix(betas) || !is.numeric(betas)) {
    stop("betas must be a numeric matrix or data.frame", call. = FALSE)
  }
  if (is.null(rownames(betas)) || is.null(colnames(betas))) {
    stop("betas must have probe IDs as rownames and sample IDs as colnames", call. = FALSE)
  }
  .stop_dups(rownames(betas), "probe IDs in betas")
  .stop_dups(colnames(betas), "sample IDs in betas")
  if (all(is.na(betas))) stop("betas is all NA", call. = FALSE)
  rng <- range(betas, na.rm = TRUE)
  if (rng[1] < 0 || rng[2] > 1) {
    stop(sprintf("betas must lie in [0, 1]; observed range [%g, %g]. Were M-values supplied?",
                 rng[1], rng[2]), call. = FALSE)
  }
  meta <- as.data.frame(meta, check.names = FALSE)
  if (!is.character(sample_col) || length(sample_col) != 1L || !sample_col %in% names(meta)) {
    stop("sample_col must name one column of meta", call. = FALSE)
  }
  ids <- as.character(meta[[sample_col]])
  .stop_dups(ids, sprintf("sample IDs in meta$%s", sample_col))
  in_b <- setdiff(colnames(betas), ids)
  in_m <- setdiff(ids, colnames(betas))
  if (length(in_b) || length(in_m)) {
    stop(sprintf(paste0("sample IDs must match one-to-one: %d in betas but not meta (e.g. %s); ",
                        "%d in meta but not betas (e.g. %s). Subset meta to the arrayed samples first."),
                 length(in_b), paste(utils::head(in_b, 3), collapse = ", "),
                 length(in_m), paste(utils::head(in_m, 3), collapse = ", ")), call. = FALSE)
  }
  meta <- meta[match(colnames(betas), ids), , drop = FALSE]
  rownames(meta) <- colnames(betas)
  SummarizedExperiment(
    assays = list(beta = betas),
    colData = DataFrame(meta, check.names = FALSE),
    rowData = DataFrame(probe_id = rownames(betas)),
    metadata = list(n_input_probes = nrow(betas))
  )
}
```

- [ ] **Step 6: Generate NAMESPACE and run the tests**

Run: `Rscript -e 'devtools::document(); devtools::test(filter = "input")'`
Expected: all tests PASS. `NAMESPACE` exports `from_matrix`.

- [ ] **Step 7: Commit**

```bash
git add DESCRIPTION LICENSE LICENSE.md .Rbuildignore NAMESPACE R/ man/ tests/
git commit -m "Scaffold CoRSIVtools and add from_matrix() input validation"
```

---

### Task 2: `read_manifest()`

**Files:**
- Create: `R/manifest.R`, `tests/testthat/test-manifest.R`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `read_manifest(x, build = "hg38")` → `GRanges`. Width-1 ranges at the 1-based CpG position, `mcols` column `probe_id` (character), and `metadata(gr)$manifest = list(source = <md5 or built-in name>, file = <basename or NA>, build = "hg38", n = <length>)`.

- [ ] **Step 1: Write the failing tests**

`tests/testthat/test-manifest.R`:
```r
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
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `Rscript -e 'devtools::test(filter = "manifest")'`
Expected: FAIL with `could not find function "read_manifest"`.

- [ ] **Step 3: Implement `R/manifest.R`**

```r
.builtin_manifests <- c(EPICv1 = "IlluminaHumanMethylationEPICanno.ilm10b5.hg38",
                        EPICv2 = "IlluminaHumanMethylationEPICv2anno.20a1.hg38")

#' Read a probe manifest as hg38 CpG coordinates
#'
#' @param x Path to a BED (`chr, start, end, probe_id`; 1-bp, 0-based) or
#'   Illumina CSV manifest, or one of `"EPICv1"`, `"EPICv2"` for the
#'   Bioconductor hg38 annotation packages.
#' @param build Genome build of the manifest. Only `"hg38"` is supported.
#' @return `GRanges` of width 1 with column `probe_id`.
#' @export
read_manifest <- function(x, build = "hg38") {
  if (!is.character(x) || length(x) != 1L) stop("x must be a single path or built-in name", call. = FALSE)
  if (!identical(build, "hg38")) {
    stop("only hg38 manifests are supported; got build = '", build, "'", call. = FALSE)
  }
  if (identical(x, "450K")) {
    stop(paste0("No hg38 450K annotation package exists on Bioconductor. Supply an hg38 450K manifest ",
                "file (e.g. from https://zwdzwd.github.io/InfiniumAnnotation) as a path."), call. = FALSE)
  }
  if (x %in% names(.builtin_manifests)) {
    gr <- .manifest_builtin(.builtin_manifests[[x]])
    src <- x; file <- NA_character_
  } else {
    if (!file.exists(x)) stop("manifest file not found: ", x, call. = FALSE)
    gr <- if (grepl("\\.bed(\\.gz)?$", x, ignore.case = TRUE)) {
      .manifest_bed(x)
    } else if (grepl("\\.csv(\\.gz)?$", x, ignore.case = TRUE)) {
      .manifest_csv(x)
    } else {
      stop("unrecognised manifest format (expected .bed or .csv): ", x, call. = FALSE)
    }
    src <- unname(tools::md5sum(x)); file <- basename(x)
  }
  .stop_dups(gr$probe_id, "probe IDs in manifest")
  metadata(gr)$manifest <- list(source = src, file = file, build = build, n = length(gr))
  gr
}

.add_chr <- function(chr) ifelse(grepl("^chr", chr), chr, paste0("chr", chr))

.make_probe_gr <- function(chr, pos, id) {
  GRanges(.add_chr(as.character(chr)), IRanges(as.integer(pos), width = 1L), probe_id = as.character(id))
}

.manifest_bed <- function(path) {
  nf <- utils::count.fields(path, sep = "\t", quote = "", comment.char = "")[1]
  if (is.na(nf) || nf < 4) stop("BED manifest needs at least 4 columns: chr, start, end, probe_id", call. = FALSE)
  cc <- c("character", "integer", "integer", "character", rep("NULL", nf - 4))
  b <- utils::read.delim(path, header = FALSE, colClasses = cc, quote = "", comment.char = "")
  if (any(b[[3]] - b[[2]] != 1L)) {
    stop("BED manifest intervals must be 1-bp (end - start == 1)", call. = FALSE)
  }
  .make_probe_gr(b[[1]], b[[3]], b[[4]])
}

.manifest_csv <- function(path) {
  lines <- readLines(path)
  a <- grep("^\\[Assay\\]", lines)
  if (length(a) != 1L) stop("Illumina CSV manifest has no [Assay] section", call. = FALSE)
  ctrl <- grep("^\\[Controls\\]", lines)
  n <- if (length(ctrl)) ctrl[1] - a - 2L else -1L
  m <- utils::read.csv(path, skip = a, nrows = n, check.names = FALSE,
                       colClasses = "character", na.strings = c("", "NA"))
  need <- c("IlmnID", "CHR", "MAPINFO")
  if (!all(need %in% names(m))) {
    stop("Illumina CSV manifest is missing columns: ",
         paste(setdiff(need, names(m)), collapse = ", "), call. = FALSE)
  }
  if ("Genome_Build" %in% names(m)) {
    gb <- unique(stats::na.omit(m$Genome_Build))
    if (!all(gb %in% c("38", "GRCh38", "hg38"))) {
      stop("Illumina CSV Genome_Build is ", paste(gb, collapse = "/"), ", not hg38", call. = FALSE)
    }
  }
  pos <- suppressWarnings(as.integer(m$MAPINFO))
  ok <- !is.na(m$CHR) & !is.na(pos) & pos > 0L
  .make_probe_gr(m$CHR[ok], pos[ok], m$IlmnID[ok])
}

.manifest_builtin <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("install the annotation package first: BiocManager::install('", pkg, "')", call. = FALSE)
  }
  e <- new.env()
  utils::data("Locations", package = pkg, envir = e)
  L <- e$Locations
  ok <- !is.na(L$chr) & !is.na(L$pos)
  .make_probe_gr(L$chr[ok], L$pos[ok], rownames(L)[ok])
}
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `Rscript -e 'devtools::document(); devtools::test(filter = "manifest")'`
Expected: PASS. The EPICv2 test runs because the package is installed.

- [ ] **Step 5: Commit**

```bash
git add R/manifest.R tests/testthat/test-manifest.R NAMESPACE man/
git commit -m "Add read_manifest() for BED, Illumina CSV and built-in EPIC hg38 manifests"
```

---

### Task 3: `annotate_probes()` with replicate collapse by coordinate

**Files:**
- Create: `R/annotate.R`, `tests/testthat/test-annotate.R`

**Interfaces:**
- Consumes: `from_matrix()` output; `read_manifest()` output (GRanges with `probe_id`); `.cpg_id()`.
- Produces: `annotate_probes(se, manifest)` → `RangedSummarizedExperiment`:
  - one row per distinct CpG coordinate; rownames = `.cpg_id()`;
  - `rowRanges` of width 1, with mcols `probe_ids` (`CharacterList`);
  - assay `"beta"`, where replicates are averaged with `na.rm`;
  - rows sorted by `GenomeInfoDb::sortSeqlevels` order and position;
  - `metadata` adds `manifest`, `n_unmapped` and `n_replicates_collapsed`.

- [ ] **Step 1: Write the failing tests**

`tests/testthat/test-annotate.R`:
```r
mk <- function() {
  b <- rbind(cg1_TC11 = c(.2, .4), cg1_TC12 = c(.4, NA), cg2_TC11 = c(.5, .5),
             cg2_TC21 = c(.9, .9), cgZ_BC11 = c(.1, .1), ctl_Neg01 = c(.3, .3))
  colnames(b) <- c("s1", "s2")
  se <- from_matrix(b, data.frame(id = c("s1", "s2")), "id")
  man <- GenomicRanges::GRanges(c("chr2", "chr2", "chr1", "chr1", "chr10"),
                                IRanges::IRanges(c(100, 100, 50, 51, 5), width = 1),
                                probe_id = c("cg1_TC11", "cg1_TC12", "cg2_TC11", "cg2_TC21", "cgZ_BC11"))
  metadata(man)$manifest <- list(source = "test")
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
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `Rscript -e 'devtools::test(filter = "annotate")'`
Expected: FAIL with `could not find function "annotate_probes"`.

- [ ] **Step 3: Implement `R/annotate.R`**

```r
#' Attach probe coordinates and collapse replicate probes by coordinate
#'
#' Probes without coordinates are dropped. Replicate probes at the same
#' chr:pos are averaged (NA-aware). Probes that share a cg base name but map
#' to different positions are kept separate.
#'
#' @param se Output of [from_matrix()].
#' @param manifest Output of [read_manifest()].
#' @return A `RangedSummarizedExperiment`, one row per CpG coordinate.
#' @export
annotate_probes <- function(se, manifest) {
  if (!is(se, "SummarizedExperiment")) stop("se must be a SummarizedExperiment", call. = FALSE)
  if (!is(manifest, "GRanges") || is.null(manifest$probe_id)) {
    stop("manifest must be the output of read_manifest()", call. = FALSE)
  }
  m <- match(rownames(se), manifest$probe_id)
  n_miss <- sum(is.na(m))
  if (n_miss == nrow(se)) {
    stop("no probe IDs in se were found in the manifest; IDs must match exactly (e.g. 'cg06185909_TC11')",
         call. = FALSE)
  }
  if (n_miss) {
    message(sprintf("%d of %d probes have no coordinates in the manifest and were dropped",
                    n_miss, nrow(se)))
  }
  keep <- !is.na(m)
  b <- assay(se, "beta")[keep, , drop = FALSE]
  gr <- granges(manifest)[m[keep]]
  ids <- rownames(se)[keep]
  key <- .cpg_id(gr)
  f <- factor(key, levels = unique(key))
  ok <- !is.na(b)
  b0 <- b; b0[!ok] <- 0
  sums <- rowsum(b0, f, reorder = FALSE)
  cnt <- rowsum(ok * 1, f, reorder = FALSE)
  vals <- sums / cnt
  vals[cnt == 0] <- NA_real_
  gr <- gr[!duplicated(key)]
  mcols(gr) <- DataFrame(probe_ids = IRanges::CharacterList(split(ids, f)))
  rownames(vals) <- levels(f)
  gr <- GenomeInfoDb::sortSeqlevels(gr)
  out <- SummarizedExperiment(
    assays = list(beta = vals), rowRanges = gr, colData = colData(se),
    metadata = c(metadata(se), list(manifest = metadata(manifest)$manifest,
                                    n_unmapped = n_miss,
                                    n_replicates_collapsed = sum(duplicated(key))))
  )
  rownames(out) <- levels(f)
  out[order(rowRanges(out)), ]
}
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `Rscript -e 'devtools::document(); devtools::test(filter = "annotate")'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/annotate.R tests/testthat/test-annotate.R NAMESPACE man/
git commit -m "Add annotate_probes(): map probes to hg38 and collapse replicates by coordinate"
```

---

### Task 4: Bundled CoRSIV 2.0 BED, `corsiv_regions()`, `write_bed()`

**Files:**
- Create: `inst/extdata/CoRSIV2.0_hg38.bed` (copied), `inst/extdata/README.md`, `inst/scripts/make_corsiv_bed.R`, `R/corsiv.R`, `tests/testthat/test-corsiv.R`

**Interfaces:**
- Consumes: `.region_id()`.
- Produces:
  - `corsiv_regions(version = "2.0", chrX = FALSE)` → `GRanges` (1-based), sorted. mcols: `region_id` (character), `source_corsivs` (character). `metadata(gr)$corsiv_version = "2.0"`.
  - `write_bed(x, file)`: `x` is a GRanges (uses `region_id` if present) or a results table with `chr, start, end, id`. Writes a 0-based BED; returns `file` invisibly.
  - Internal: `.read_region_bed(path)` → GRanges as above.

- [ ] **Step 1: Copy the BED and record its provenance**

```bash
mkdir -p inst/extdata inst/scripts
SRC="$HOME/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/WORK/WaterlandLab/XTraCoRSIV_Addanki/Baylor_Addanki_10022026/Normalized_Data/data/raw/CoRSIV2.0_ESS_chrX_male_female_all_cooridnates_100bp_hg38.bedtools_merged_fixed.bed"
cp "$SRC" inst/extdata/CoRSIV2.0_hg38.bed
md5 -q inst/extdata/CoRSIV2.0_hg38.bed   # expect ca1b9a6c50e78e9a43a92f2e56400fb2
```

`inst/extdata/README.md`:
```markdown
# CoRSIV2.0_hg38.bed

- Contents: CoRSIV 2.0 regions (autosomes + chrX; ESS; male and female), hg38,
  0-based BED, 4 columns (`chr, start, end, source_corsivs`). Overlapping
  population-specific CoRSIVs were merged with `bedtools merge`; the name column
  lists every merged source CoRSIV, comma-separated.
- Regions: 6,669 (320 on chrX). Non-overlapping and non-abutting.
- Original file name: `CoRSIV2.0_ESS_chrX_male_female_all_cooridnates_100bp_hg38.bedtools_merged_fixed.bed`
- md5: `ca1b9a6c50e78e9a43a92f2e56400fb2`
- Source: Waterland Lab, Baylor College of Medicine (CoRSIV 2.0; publication forthcoming).
- Licence: MIT, same as the package (approved 2026-10-12).
- How it was made: see `inst/scripts/make_corsiv_bed.R`.
```

`inst/scripts/make_corsiv_bed.R`:
```r
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
```

- [ ] **Step 2: Write the failing tests**

`tests/testthat/test-corsiv.R`:
```r
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
```

- [ ] **Step 3: Run the tests and confirm they fail**

Run: `Rscript -e 'devtools::test(filter = "corsiv")'`
Expected: FAIL with `could not find function "corsiv_regions"`.

- [ ] **Step 4: Implement `R/corsiv.R`**

```r
.corsiv_files <- c("2.0" = "CoRSIV2.0_hg38.bed")

#' CoRSIV regions bundled with the package
#'
#' @param version CoRSIV release. Currently only `"2.0"`.
#' @param chrX Include chrX CoRSIVs? Default `FALSE` (chrX needs
#'   sex-stratified models, not supported in v0.1).
#' @return `GRanges` (1-based) with `region_id` (`chr:start-end`) and
#'   `source_corsivs` (original merged CoRSIV names, not unique).
#' @export
corsiv_regions <- function(version = "2.0", chrX = FALSE) {
  if (!version %in% names(.corsiv_files)) {
    stop("unknown CoRSIV version '", version, "'; available: ",
         paste(names(.corsiv_files), collapse = ", "), call. = FALSE)
  }
  path <- system.file("extdata", .corsiv_files[[version]], package = "CoRSIVtools", mustWork = TRUE)
  gr <- .read_region_bed(path)
  if (!chrX) gr <- gr[as.character(seqnames(gr)) != "chrX"]
  metadata(gr)$corsiv_version <- version
  gr
}

.read_region_bed <- function(path) {
  b <- utils::read.delim(path, header = FALSE, quote = "", comment.char = "",
                         colClasses = c("character", "integer", "integer", "character"))
  gr <- GRanges(b[[1]], IRanges(b[[2]] + 1L, b[[3]]))
  gr$region_id <- .region_id(gr)
  gr$source_corsivs <- b[[4]]
  gr <- GenomeInfoDb::sortSeqlevels(gr)
  sort(gr)
}

#' Write regions or results as a 0-based BED file for bedtools
#'
#' @param x `GRanges` (name column = `region_id`) or a results table with
#'   columns `chr`, `start`, `end`, `id` (1-based, as returned by [test_dm()]).
#' @param file Output path.
#' @return `file`, invisibly.
#' @export
write_bed <- function(x, file) {
  if (is(x, "GRanges")) {
    nm <- if (!is.null(x$region_id)) x$region_id else .region_id(x)
    df <- data.frame(as.character(seqnames(x)), start(x) - 1L, end(x), nm)
  } else {
    x <- as.data.frame(x)
    need <- c("chr", "start", "end", "id")
    if (!all(need %in% names(x))) {
      stop("x needs columns: ", paste(need, collapse = ", "), call. = FALSE)
    }
    df <- data.frame(x$chr, as.integer(x$start) - 1L, as.integer(x$end), x$id)
  }
  utils::write.table(df, file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  invisible(file)
}
```

- [ ] **Step 5: Run the tests and confirm they pass**

Run: `Rscript -e 'devtools::document(); devtools::test(filter = "corsiv")'`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add inst/ R/corsiv.R tests/testthat/test-corsiv.R NAMESPACE man/
git commit -m "Bundle CoRSIV 2.0 BED; add corsiv_regions() with coordinate IDs and write_bed()"
```

---

### Task 5: `simulate_corsiv_data()`, the `CoRSIVSets` class, and `define_sets()` (user and array modes)

**Files:**
- Create: `R/simulate.R`, `R/sets.R`, `tests/testthat/test-simulate.R`, `tests/testthat/test-sets.R`

**Interfaces:**
- Consumes: `from_matrix()`, `annotate_probes()`, `.region_id()`, `.cpg_id()`, `.check_annotated()`.
- Produces:
  - `simulate_corsiv_data(n_samples = 20, n_corsivs = 40, n_dm = 10, shift = 0.15, seed = 1)` → `list(se, manifest, corsivs, controls, dm_ids)`:
    - `se` is a probe-level SE from `from_matrix()` with `colData` columns `sample_id`, `group` (factor with levels `A`, `B`) and `age`;
    - `manifest` is a GRanges with `probe_id`;
    - `corsivs` is a GRanges with `region_id` and `source_corsivs`, and `metadata(corsivs)$corsiv_version = "sim"`;
    - `controls` is a GRanges of the non-CoRSIV clusters (each has 3–8 probes; use it as a user background);
    - `dm_ids` lists the `region_id`s that were shifted.
  - S4 class `CoRSIVSets`, with slots:
    - `target`: GRanges; mcols `region_id, source_corsivs, n_cpg, set, draw`;
    - `background`: GRanges; same five mcols, plus `matched_to`, `cpg_density`, `target_density` in matched mode;
    - `corsivs`: every supplied CoRSIV, with `n_cpg`;
    - `mode`: `"matched"`, `"array"` or `"user"`;
    - `params`: list of `corsiv_version, min_cpgs, n_draws, seed, flank, max_tries`;
    - `diagnostics`: DataFrame.
  - Accessors `targetRegions(x)`, `backgroundRegions(x)`, `matchDiagnostics(x)`.
  - `define_sets(se, corsivs = corsiv_regions(), background = "matched", n_draws = 100, seed = 1, flank = 1000, min_cpgs = 3, max_tries = 1000)` → `CoRSIVSets`. In this task the `"matched"` branch calls `.match_background()`, which Task 6 writes; until then that branch isn't tested.
  - Internal `.std_cols <- c("region_id", "source_corsivs", "n_cpg", "set", "draw")`.

- [ ] **Step 1: Write the failing tests**

`tests/testthat/test-simulate.R`:
```r
test_that("simulate_corsiv_data returns a consistent, reproducible dataset", {
  s <- simulate_corsiv_data()
  expect_named(s, c("se", "manifest", "corsivs", "controls", "dm_ids"))
  expect_length(s$controls, 360L)
  expect_false(any(IRanges::overlapsAny(s$controls, s$corsivs)))
  expect_identical(ncol(s$se), 20L)
  expect_identical(levels(s$se$group), c("A", "B"))
  expect_length(s$corsivs, 40L)
  expect_length(s$dm_ids, 10L)
  expect_true(all(s$dm_ids %in% s$corsivs$region_id))
  expect_true(all(rownames(s$se) %in% s$manifest$probe_id))
  a <- suppressMessages(annotate_probes(s$se, s$manifest))
  expect_true(all(GenomicRanges::countOverlaps(s$corsivs, rowRanges(a)) >= 3L))
  expect_identical(simulate_corsiv_data()$dm_ids, s$dm_ids)
  # same-position replicate collapsed; same base name at a different position kept separate
  expect_true(any(lengths(rowRanges(a)$probe_ids) == 2L))
  base2 <- sub("_TC11$", "", rownames(s$se)[2])
  hits <- vapply(as.list(rowRanges(a)$probe_ids), function(p) any(startsWith(p, base2)), logical(1))
  expect_identical(sum(hits), 2L)
})
```

`tests/testthat/test-sets.R`:
```r
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
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `Rscript -e 'devtools::test(filter = "simulate|sets")'`
Expected: FAIL with `could not find function "simulate_corsiv_data"`.

- [ ] **Step 3: Implement `R/simulate.R`**

```r
#' Simulate a small CoRSIV methylation dataset on real hg38 chr21 sequence
#'
#' Probes sit on real CpG positions in chr21:30,000,001-34,000,000, so
#' CpG-density matching works on genuine sequence. 400 probe clusters
#' (300-600 bp, 3-8 probes) are laid out one per 10-kb bin; `n_corsivs` of
#' them become CoRSIVs and the rest are returned as `controls` (usable as a
#' user background). CoRSIVs get high
#' between-individual variation (SD 0.10) vs other clusters (SD 0.03).
#' `n_dm` CoRSIVs are shifted by `shift` in group B. Singleton probes are added
#' between clusters. One CpG carries a same-position replicate probe and one
#' probe ID shares its cg base name with a probe elsewhere.
#'
#' @param n_samples Even number of samples, split equally into groups A and B.
#' @param n_corsivs Number of the 400 clusters that are CoRSIVs (<= 200).
#' @param n_dm Number of CoRSIVs with a group effect.
#' @param shift Beta-scale group effect.
#' @param seed Random seed.
#' @return `list(se, manifest, corsivs, controls, dm_ids)`.
#' @export
simulate_corsiv_data <- function(n_samples = 20, n_corsivs = 40, n_dm = 10, shift = 0.15, seed = 1) {
  stopifnot(n_samples %% 2 == 0, n_dm <= n_corsivs, n_corsivs <= 200)
  withr::with_seed(seed, {
    genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    s0 <- 30000001L; s1 <- 34000000L
    cg <- BiocGenerics::start(Biostrings::matchPattern("CG", genome[["chr21"]][s0:s1])) + s0 - 1L
    n_clust <- 400L
    bin <- floor((s1 - s0 + 1) / n_clust)
    starts <- integer(n_clust); widths <- integer(n_clust); probes <- vector("list", n_clust)
    for (j in seq_len(n_clust)) {
      for (attempt in 1:50) {
        w <- sample(c(300L, 400L, 500L, 600L), 1)
        st <- s0 + (j - 1L) * bin + sample.int(bin - w - 2000L, 1)
        inside <- cg[cg >= st & cg <= st + w - 1L]
        if (length(inside) >= 8L) break
      }
      starts[j] <- st; widths[j] <- w; probes[[j]] <- sort(sample(inside, sample(3:8, 1)))
    }
    clusters <- GRanges("chr21", IRanges(starts, width = widths))
    corsiv_idx <- sort(sample.int(n_clust, n_corsivs))
    # singletons between clusters (array background); two spare positions for replicate fixtures
    out_cl <- cg[!IRanges::overlapsAny(GRanges("chr21", IRanges(cg, width = 1L)), clusters + 1000L)]
    single <- sort(sample(out_cl, 1001L))
    spare <- single[1001L]; single <- single[-1001L]
    pos <- c(unlist(probes), single)
    clust_of <- c(rep(seq_len(n_clust), lengths(probes)), rep(NA_integer_, length(single)))
    ids <- sprintf("cg%08d_TC11", seq_along(pos))
    # fixture 1: same-position replicate of probe 1 (collapsed by annotate_probes)
    # fixture 2: probe 2's base name at a different position (must stay a separate CpG)
    ids <- c(ids, sub("_TC11$", "_TC12", ids[1]), sub("_TC11$", "_TC21", ids[2]))
    pos <- c(pos, pos[1], spare)
    clust_of <- c(clust_of, clust_of[1], NA_integer_)
    n <- n_samples
    group <- factor(rep(c("A", "B"), each = n / 2), levels = c("A", "B"))
    in_cv <- !is.na(clust_of) & clust_of %in% corsiv_idx
    sd_r <- ifelse(in_cv, 0.10, 0.03)
    mu <- stats::runif(length(pos), 0.2, 0.8)
    reff <- matrix(stats::rnorm(n_clust * n), n_clust, n)
    b <- t(vapply(seq_along(pos), function(i) {
      r <- if (is.na(clust_of[i])) stats::rnorm(n) else reff[clust_of[i], ]
      mu[i] + sd_r[i] * r + stats::rnorm(n, 0, 0.02)
    }, numeric(n)))
    dm_idx <- corsiv_idx[sample.int(length(corsiv_idx), n_dm)]
    dm_rows <- which(clust_of %in% dm_idx)
    b[dm_rows, group == "B"] <- b[dm_rows, group == "B"] + shift
    b <- pmin(pmax(b, 0.001), 0.999)
    sids <- sprintf("S%02d", seq_len(n))
    dimnames(b) <- list(ids, sids)
    meta <- data.frame(sample_id = sids, group = group, age = round(stats::runif(n, 40, 70), 1))
    se <- from_matrix(b, meta, "sample_id")
    manifest <- GRanges("chr21", IRanges(pos, width = 1L), probe_id = ids)
    metadata(manifest)$manifest <- list(source = "simulated", file = NA_character_, build = "hg38",
                                        n = length(ids))
    corsivs <- clusters[corsiv_idx]
    corsivs$region_id <- .region_id(corsivs)
    corsivs$source_corsivs <- sprintf("sim_corsiv_%03d", seq_along(corsivs))
    metadata(corsivs)$corsiv_version <- "sim"
    controls <- clusters[-corsiv_idx]
    list(se = se, manifest = manifest, corsivs = corsivs, controls = controls,
         dm_ids = .region_id(clusters[dm_idx]))
  })
}
```

- [ ] **Step 4: Implement `R/sets.R`**

```r
.std_cols <- c("region_id", "source_corsivs", "n_cpg", "set", "draw")

#' Target CoRSIVs and background regions for a CoRSIV analysis
#'
#' Created by [define_sets()]. Use [targetRegions()], [backgroundRegions()]
#' and [matchDiagnostics()] to access the contents.
#' @slot target,background,corsivs `GRanges`.
#' @slot mode `"matched"`, `"array"` or `"user"`.
#' @slot params list of settings.
#' @slot diagnostics `DataFrame` of per-draw matching results.
#' @export
setClass("CoRSIVSets", representation(target = "GRanges", background = "GRanges",
                                      corsivs = "GRanges", mode = "character",
                                      params = "list", diagnostics = "DataFrame"))

#' @rdname CoRSIVSets-class
#' @param x A `CoRSIVSets` object.
#' @export
targetRegions <- function(x) x@target
#' @rdname CoRSIVSets-class
#' @export
backgroundRegions <- function(x) x@background
#' @rdname CoRSIVSets-class
#' @export
matchDiagnostics <- function(x) x@diagnostics

#' @rdname CoRSIVSets-class
#' @param object A `CoRSIVSets` object.
#' @export
setMethod("show", "CoRSIVSets", function(object) {
  cat("CoRSIVSets\n")
  cat(sprintf("  CoRSIV set: %s; %d supplied, %d targets with >= %d CpGs\n",
              object@params$corsiv_version, length(object@corsivs), length(object@target),
              object@params$min_cpgs))
  cat(sprintf("  background: %s; %d regions across %d draw(s)\n", object@mode,
              length(object@background), length(unique(object@background$draw))))
})

#' Define target CoRSIVs and the background to compare them with
#'
#' @param se Output of [annotate_probes()].
#' @param corsivs `GRanges` of CoRSIVs; default [corsiv_regions()].
#' @param background `"matched"` (random control regions matched on width,
#'   probe count, chromosome type and CpG density), `"array"` (every CpG
#'   outside CoRSIVs +/- `flank`; CpG level only), or a `GRanges` of
#'   user control regions (must not overlap CoRSIVs).
#' @param n_draws Number of matched draws.
#' @param seed Base seed; draw r uses `seed + r`.
#' @param flank bp excluded around every CoRSIV when choosing background.
#' @param min_cpgs Minimum distinct CpGs for a region to be analysed.
#' @param max_tries Maximum candidate proposals per target per draw.
#' @return A [CoRSIVSets-class] object.
#' @export
define_sets <- function(se, corsivs = corsiv_regions(), background = "matched", n_draws = 100,
                        seed = 1, flank = 1000, min_cpgs = 3, max_tries = 1000) {
  .check_annotated(se)
  cpg <- rowRanges(se)
  if (is.null(corsivs$region_id)) corsivs$region_id <- .region_id(corsivs)
  if (is.null(corsivs$source_corsivs)) corsivs$source_corsivs <- NA_character_
  corsivs$n_cpg <- countOverlaps(corsivs, cpg)
  target <- corsivs[corsivs$n_cpg >= min_cpgs]
  if (!length(target)) {
    stop(sprintf(paste0("no CoRSIV has >= %d CpGs on this array. Check that chromosome names match ",
                        "(e.g. 'chr1' in both the manifest and the CoRSIVs)."), min_cpgs), call. = FALSE)
  }
  target$set <- "target"; target$draw <- NA_integer_
  mcols(target) <- mcols(target)[, .std_cols]
  excl <- reduce(granges(corsivs) + flank)
  diag <- DataFrame()
  if (is(background, "GRanges")) {
    if (any(overlapsAny(background, corsivs))) {
      stop("user background regions overlap CoRSIVs; remove the overlapping regions first", call. = FALSE)
    }
    bg <- granges(background)
    bg$region_id <- .region_id(bg); bg$source_corsivs <- NA_character_
    bg$n_cpg <- countOverlaps(bg, cpg); bg <- bg[bg$n_cpg >= min_cpgs]
    bg$set <- rep("background", length(bg)); bg$draw <- rep(1L, length(bg))
    mode <- "user"
  } else if (identical(background, "array")) {
    bg <- granges(cpg[!overlapsAny(cpg, excl)])
    mcols(bg) <- DataFrame(region_id = .cpg_id(bg), source_corsivs = NA_character_, n_cpg = 1L,
                           set = "background", draw = 1L)
    mode <- "array"
  } else if (identical(background, "matched")) {
    m <- .match_background(target, cpg, excl, n_draws, seed, max_tries)
    bg <- m$background; diag <- m$diagnostics
    mode <- "matched"
  } else {
    stop("background must be 'matched', 'array', or a GRanges of control regions", call. = FALSE)
  }
  new("CoRSIVSets", target = target, background = bg, corsivs = corsivs, mode = mode,
      params = list(corsiv_version = if (is.null(metadata(corsivs)$corsiv_version)) NA_character_
                                     else metadata(corsivs)$corsiv_version,
                    min_cpgs = min_cpgs, n_draws = if (mode == "matched") n_draws else 1L,
                    seed = seed, flank = flank, max_tries = max_tries),
      diagnostics = diag)
}
```

- [ ] **Step 5: Run the tests and confirm they pass**

Run: `Rscript -e 'devtools::document(); devtools::test(filter = "simulate|sets")'`
Expected: PASS. The tests don't exercise `.match_background()`. R CMD check will flag it as undefined until Task 6; that's expected.

- [ ] **Step 6: Commit**

```bash
git add R/simulate.R R/sets.R tests/testthat/test-simulate.R tests/testthat/test-sets.R NAMESPACE man/
git commit -m "Add simulate_corsiv_data(), CoRSIVSets class and define_sets() user/array backgrounds"
```

---

### Task 6: Matched background sampling

**Files:**
- Create: `R/background.R`, `tests/testthat/test-background.R`

**Interfaces:**
- Consumes: `target` GRanges (mcols `.std_cols`), `cpg` (sorted rowRanges of the annotated SE), `excl` (GRanges of CoRSIVs ± flank), `.chr_type()`, `.region_id()`.
- Produces:
  - `.match_background(target, cpg, excl, n_draws, seed, max_tries, batch = 20L)` → `list(background = GRanges, diagnostics = DataFrame)`.
    - `background` mcols: `.std_cols`, plus `matched_to` (target `region_id`), `cpg_density`, `target_density`.
    - `diagnostics` columns: `draw, n_target, n_matched, n_unmatched`.
  - Internal helpers: `.count_in_windows(chr, s, e, pos_by_chr)` → integer; `.cpg_density(gr)` → numeric (CpGs per 100 bp, hg38).

- [ ] **Step 1: Write the failing tests**

`tests/testthat/test-background.R`:
```r
sim_sets <- function(n_draws = 3, seed = 1, ...) {
  s <- simulate_corsiv_data()
  a <- suppressMessages(annotate_probes(s$se, s$manifest))
  list(s = s, a = a, x = define_sets(a, s$corsivs, background = "matched",
                                     n_draws = n_draws, seed = seed, ...))
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
  expect_warning(m <- CoRSIVtools:::.match_background(tg2, rowRanges(r$a), excl, 1, 1, 50),
                 "unmatched")
  expect_identical(m$diagnostics$n_unmatched >= 1L, TRUE)
  expect_false("chrX:1000000-1000399" %in% m$background$matched_to)
})

test_that("windows past the chromosome end are rejected before getSeq", {  # Review Focus 5
  ok <- CoRSIVtools:::.in_bounds(c("chr21", "chr21"), c(46709900L, 100L), c(46710100L, 500L))
  expect_identical(ok, c(FALSE, TRUE))   # chr21 length is 46,709,983
})
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `Rscript -e 'devtools::test(filter = "background")'`
Expected: FAIL with `object '.count_in_windows' not found`.

- [ ] **Step 3: Implement `R/background.R`**

```r
.hg38 <- function() BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

.cpg_density <- function(gr) {
  seqs <- Biostrings::getSeq(.hg38(), granges(gr))
  Biostrings::vcountPattern("CG", seqs) / width(gr) * 100
}

.in_bounds <- function(chr, s, e) {
  len <- GenomeInfoDb::seqlengths(.hg38())[chr]
  !is.na(len) & s >= 1L & e <= len
}

.count_in_windows <- function(chr, s, e, pos_by_chr) {
  out <- integer(length(chr))
  for (ch in unique(chr)) {
    p <- pos_by_chr[[ch]]
    if (is.null(p)) next
    i <- which(chr == ch)
    out[i] <- findInterval(e[i], p) - findInterval(s[i] - 1L, p)
  }
  out
}

.match_background <- function(target, cpg, excl, n_draws, seed, max_tries, batch = 20L) {
  target$cpg_density <- .cpg_density(target)
  pos_by_chr <- split(start(cpg), as.character(seqnames(cpg)))
  pos_by_chr <- lapply(pos_by_chr, sort)
  pool <- granges(cpg[!overlapsAny(cpg, excl)])
  pool_type <- .chr_type(as.character(seqnames(pool)))
  no_pool <- !.chr_type(as.character(seqnames(target))) %in% unique(pool_type)
  if (any(no_pool)) {
    warning(sprintf("%d target(s) are on chromosome types with no background probes and will be unmatched",
                    sum(no_pool)), call. = FALSE)
  }
  out <- vector("list", n_draws); diag <- vector("list", n_draws)
  for (r in seq_len(n_draws)) {
    ctrl <- withr::with_seed(seed + r, .match_one_draw(target, pool, pool_type, pos_by_chr, excl,
                                                       max_tries, batch))
    ctrl$draw <- rep(as.integer(r), length(ctrl))
    out[[r]] <- ctrl
    diag[[r]] <- DataFrame(draw = as.integer(r), n_target = length(target), n_matched = length(ctrl),
                           n_unmatched = length(target) - length(ctrl))
  }
  diag <- do.call(rbind, diag)
  bad <- diag$n_unmatched / diag$n_target > 0.05
  if (any(bad)) {
    warning(sprintf("%d of %d draws left >5%% of targets unmatched (max %d unmatched); see matchDiagnostics()",
                    sum(bad), n_draws, max(diag$n_unmatched)), call. = FALSE)
  }
  bg <- unlist(GRangesList(out), use.names = FALSE)
  list(background = bg, diagnostics = diag)
}

.match_one_draw <- function(target, pool, pool_type, pos_by_chr, excl, max_tries, batch) {
  n <- length(target); w <- width(target); k <- target$n_cpg; d <- target$cpg_density
  tchr <- .chr_type(as.character(seqnames(target)))
  pool_chr <- as.character(seqnames(pool)); pool_pos <- start(pool)
  tries <- integer(n); done <- logical(n)
  for (ty in unique(tchr)) if (!any(pool_type == ty)) tries[tchr == ty] <- max_tries
  acc <- GRanges(); acc_t <- integer(); acc_k <- integer(); acc_d <- numeric()
  repeat {
    todo <- which(!done & tries < max_tries)
    if (!length(todo)) break
    nb <- pmin(batch, max_tries - tries[todo])
    ti <- rep(todo, nb)
    tries[todo] <- tries[todo] + nb
    anc <- integer(length(ti))
    for (ty in unique(tchr[ti])) {
      sel <- tchr[ti] == ty; cand <- which(pool_type == ty)
      anc[sel] <- cand[sample.int(length(cand), sum(sel), replace = TRUE)]
    }
    a_chr <- pool_chr[anc]
    s <- pool_pos[anc] - as.integer(floor(stats::runif(length(ti)) * w[ti]))
    e <- s + w[ti] - 1L
    keep <- .in_bounds(a_chr, s, e)
    cnt <- .count_in_windows(a_chr, s, e, pos_by_chr)
    keep <- keep & ifelse(k[ti] < 10L, cnt == k[ti], abs(cnt - k[ti]) <= 1L)
    if (!any(keep)) next
    ti <- ti[keep]; cnt <- cnt[keep]
    cand <- GRanges(a_chr[keep], IRanges(s[keep], e[keep]))
    keep <- !overlapsAny(cand, excl)
    if (length(acc)) keep <- keep & !overlapsAny(cand, acc)
    if (!any(keep)) next
    ti <- ti[keep]; cnt <- cnt[keep]; cand <- cand[keep]
    dens <- .cpg_density(cand)
    keep <- abs(dens - d[ti]) <= pmax(0.2 * d[ti], 0.25)
    keep <- keep & !duplicated(ti) & !done[ti]
    if (!any(keep)) next
    ti <- ti[keep]; cnt <- cnt[keep]; cand <- cand[keep]; dens <- dens[keep]
    h <- findOverlaps(cand, drop.self = TRUE, drop.redundant = TRUE)
    if (length(h)) {
      drop <- unique(pmax(queryHits(h), subjectHits(h)))
      ti <- ti[-drop]; cnt <- cnt[-drop]; cand <- cand[-drop]; dens <- dens[-drop]
    }
    acc <- c(acc, cand); acc_t <- c(acc_t, ti); acc_k <- c(acc_k, cnt); acc_d <- c(acc_d, dens)
    done[ti] <- TRUE
  }
  o <- order(acc_t)
  acc <- acc[o]
  mcols(acc) <- DataFrame(region_id = .region_id(acc), source_corsivs = rep(NA_character_, length(acc)),
                          n_cpg = acc_k[o], set = rep("background", length(acc)),
                          draw = rep(NA_integer_, length(acc)), matched_to = target$region_id[acc_t[o]],
                          cpg_density = acc_d[o], target_density = d[acc_t[o]])
  acc
}
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `Rscript -e 'devtools::document(); devtools::test(filter = "background|sets")'`
Expected: PASS. If the criteria test fails on CpG density only, check that `.cpg_density()` is applied to the **target** GRanges in `.match_background()` before `.match_one_draw()`, and that `d` is read from `target$cpg_density`.

- [ ] **Step 5: Commit**

```bash
git add R/background.R tests/testthat/test-background.R NAMESPACE man/
git commit -m "Add matched background sampling on width, probe count, chromosome type and CpG density"
```

---

### Task 7: `coverage_report()`

**Files:**
- Create: `R/coverage.R`, `tests/testthat/test-coverage.R`

**Interfaces:**
- Consumes: `CoRSIVSets` (slots `corsivs`, `target`, `background`, `params`).
- Produces: `coverage_report(sets)` → `DataFrame` with one row each for `set = "corsiv"` and `"background"`, and columns `set, n_regions, n_ge1, n_ge2, n_ge3, n_ge5, median_cpg, max_cpg, min_cpgs, n_kept`.

- [ ] **Step 1: Write the failing test**

`tests/testthat/test-coverage.R`:
```r
test_that("coverage_report counts CoRSIV coverage exactly", {
  s <- simulate_corsiv_data()
  a <- suppressMessages(annotate_probes(s$se, s$manifest))
  x <- define_sets(a, s$corsivs, background = "array")
  cr <- coverage_report(x)
  expect_s4_class(cr, "DataFrame")
  expect_identical(cr$set, c("corsiv", "background"))
  n <- GenomicRanges::countOverlaps(s$corsivs, rowRanges(a))
  expect_identical(cr$n_regions[1], length(n))
  expect_identical(cr$n_ge3[1], sum(n >= 3))
  expect_identical(cr$n_ge5[1], sum(n >= 5))
  expect_identical(cr$n_kept[1], length(targetRegions(x)))
  expect_identical(cr$n_kept[2], length(backgroundRegions(x)))
  expect_identical(cr$min_cpgs, c(3, 3))
})
```

- [ ] **Step 2: Run the test and confirm it fails**

Run: `Rscript -e 'devtools::test(filter = "coverage")'`
Expected: FAIL with `could not find function "coverage_report"`.

- [ ] **Step 3: Implement `R/coverage.R`**

```r
#' How many CoRSIVs (and background regions) the array covers
#'
#' @param sets Output of [define_sets()].
#' @return `DataFrame`, one row for the supplied CoRSIVs and one for the
#'   background (pooled across draws).
#' @export
coverage_report <- function(sets) {
  stopifnot(is(sets, "CoRSIVSets"))
  row <- function(label, n, kept) {
    DataFrame(set = label, n_regions = length(n), n_ge1 = sum(n >= 1), n_ge2 = sum(n >= 2),
              n_ge3 = sum(n >= 3), n_ge5 = sum(n >= 5),
              median_cpg = if (any(n > 0)) stats::median(n[n > 0]) else NA_real_,
              max_cpg = max(c(0L, n)), min_cpgs = sets@params$min_cpgs, n_kept = kept)
  }
  rbind(row("corsiv", sets@corsivs$n_cpg, length(sets@target)),
        row("background", sets@background$n_cpg, length(sets@background)))
}
```

- [ ] **Step 4: Run the test and confirm it passes**

Run: `Rscript -e 'devtools::document(); devtools::test(filter = "coverage")'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/coverage.R tests/testthat/test-coverage.R NAMESPACE man/
git commit -m "Add coverage_report() for CoRSIV and background coverage on the array"
```

---

### Task 8: `summarise_regions()`

**Files:**
- Create: `R/summarise.R`, `tests/testthat/test-summarise.R`

**Interfaces:**
- Consumes: annotated SE; `CoRSIVSets` (`target`, `background`, `mode`, `params`); `.std_cols`.
- Produces: `summarise_regions(se, sets, fun = c("mean", "median"), max_na_frac = 0.2)` → `RangedSummarizedExperiment`, one row per region:
  - assay `"beta"`;
  - `rowRanges` mcols = `.std_cols`;
  - rownames: `region_id` for targets, `paste0(region_id, "|draw", draw)` for background;
  - `metadata` adds `sets = c(params, mode = mode)`, `summarise = list(fun, max_na_frac)`, `level = "region"` and `n_dropped_na` (named integer by set).

- [ ] **Step 1: Write the failing tests**

`tests/testthat/test-summarise.R`:
```r
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
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `Rscript -e 'devtools::test(filter = "summarise")'`
Expected: FAIL with `could not find function "summarise_regions"`.

- [ ] **Step 3: Implement `R/summarise.R`**

```r
#' Summarise CpG methylation per region (CoRSIV and background)
#'
#' @param se Output of [annotate_probes()].
#' @param sets Output of [define_sets()] (not `background = "array"`).
#' @param fun `"mean"` (default) or `"median"` across the region's CpGs.
#' @param max_na_frac Drop regions with more than this fraction of NA samples.
#' @return Region-level `RangedSummarizedExperiment`.
#' @export
summarise_regions <- function(se, sets, fun = c("mean", "median"), max_na_frac = 0.2) {
  fun <- match.arg(fun)
  .check_annotated(se)
  stopifnot(is(sets, "CoRSIVSets"))
  if (sets@mode == "array") {
    stop("background = 'array' is CpG level only; use test_dm(level = 'cpg', sets = sets)", call. = FALSE)
  }
  min_cpgs <- sets@params$min_cpgs
  bg <- sets@background; mcols(bg) <- mcols(bg)[, .std_cols]
  regions <- suppressWarnings(c(sets@target, bg))
  grp <- paste(regions$set, regions$draw)
  blocks <- split(seq_along(regions), factor(grp, levels = unique(grp)))
  b <- assay(se, "beta"); cpg <- rowRanges(se)
  vals <- do.call(rbind, lapply(blocks, function(i) .summarise_block(regions[i], cpg, b, fun, min_cpgs)))
  regions <- regions[unlist(blocks, use.names = FALSE)]
  na_frac <- rowMeans(is.na(vals))
  keep <- na_frac <= max_na_frac
  dropped <- vapply(c("target", "background"), function(s) sum(!keep & regions$set == s), integer(1))
  if (any(dropped)) {
    message(sprintf("dropped %d target and %d background regions with > %g NA samples",
                    dropped[["target"]], dropped[["background"]], max_na_frac))
  }
  regions <- regions[keep]; vals <- vals[keep, , drop = FALSE]
  rn <- ifelse(regions$set == "target", regions$region_id, paste0(regions$region_id, "|draw", regions$draw))
  rownames(vals) <- rn
  colnames(vals) <- colnames(se)
  out <- SummarizedExperiment(
    assays = list(beta = vals), rowRanges = regions, colData = colData(se),
    metadata = c(metadata(se), list(sets = c(sets@params, mode = sets@mode),
                                    summarise = list(fun = fun, max_na_frac = max_na_frac),
                                    level = "region", n_dropped_na = dropped))
  )
  rownames(out) <- rn
  out
}

.summarise_block <- function(reg, cpg, b, fun, min_cpgs) {
  out <- matrix(NA_real_, length(reg), ncol(b))
  h <- findOverlaps(reg, cpg)
  if (!length(h)) return(out)
  q <- queryHits(h)
  sub <- b[subjectHits(h), , drop = FALSE]
  ok <- !is.na(sub)
  cnt <- rowsum(ok * 1, q, reorder = TRUE)
  if (fun == "mean") {
    sub0 <- sub; sub0[!ok] <- 0
    v <- rowsum(sub0, q, reorder = TRUE) / cnt
  } else {
    v <- do.call(rbind, lapply(split(seq_along(q), q), function(i)
      apply(sub[i, , drop = FALSE], 2, stats::median, na.rm = TRUE)))
  }
  v[cnt < min_cpgs] <- NA_real_
  out[as.integer(rownames(cnt)), ] <- v
  out
}
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `Rscript -e 'devtools::document(); devtools::test(filter = "summarise")'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/summarise.R tests/testthat/test-summarise.R NAMESPACE man/
git commit -m "Add summarise_regions(): NA-aware per-region beta summaries for target and background"
```

---

### Task 9: `test_dm()` at both levels, with provenance

**Files:**
- Create: `R/test_dm.R`, `R/provenance.R`, `tests/testthat/test-test_dm.R`

**Interfaces:**
- Consumes:
  - CoRSIV level: the region SE from `summarise_regions()`, using `metadata$level`, `$sets`, `$summarise` and `$manifest`.
  - CpG level: the annotated SE plus `CoRSIVSets`.
- Produces: `test_dm(x, design, coef, level = c("corsiv", "cpg"), sets = NULL)` → `DataFrame`.
  - Columns: `id, region_id, set, draw, chr, start, end, n_cpg, source_corsivs, delta_beta, logFC_M, AveM, t, P, FDR, sd_beta`.
  - `metadata(res)` holds every provenance field from spec §6.5: `CoRSIVtools_version, corsiv_version, manifest, background_mode, n_draws, seed, flank, min_cpgs, fun, max_na_frac, design, coef, level, n_samples, R_version, Bioc_version, timestamp`.

- [ ] **Step 1: Write the failing tests**

`tests/testthat/test-test_dm.R`:
```r
setup_sim <- function(n_dm = 10, n_corsivs = 40, background = "array") {
  s <- simulate_corsiv_data(n_dm = n_dm, n_corsivs = n_corsivs)
  a <- suppressMessages(annotate_probes(s$se, s$manifest))
  list(s = s, a = a, x = define_sets(a, s$corsivs, background = background))
}
user_bg <- function(r) r$s$controls

test_that("CoRSIV-level test recovers simulated DM regions with correct effect size", {
  r <- setup_sim()
  x <- define_sets(r$a, r$s$corsivs, background = user_bg(r))
  reg <- summarise_regions(r$a, x, max_na_frac = 1)
  res <- test_dm(reg, design = ~ group + age, coef = "groupB")
  expect_s4_class(res, "DataFrame")
  expect_identical(names(res), c("id", "region_id", "set", "draw", "chr", "start", "end", "n_cpg",
                                 "source_corsivs", "delta_beta", "logFC_M", "AveM", "t", "P", "FDR", "sd_beta"))
  tg <- res[res$set == "target", ]
  hit <- tg$id[tg$FDR < 0.05]
  expect_gte(sum(r$s$dm_ids %in% hit), 8L)
  expect_lt(abs(mean(tg$delta_beta[tg$id %in% r$s$dm_ids]) - 0.15), 0.03)
  expect_identical(tg$id, tg$region_id)
})

test_that("under the null, P < 0.05 occurs at about the nominal rate", {
  r <- setup_sim(n_dm = 0, n_corsivs = 200)
  x <- define_sets(r$a, r$s$corsivs, background = user_bg(r))
  reg <- summarise_regions(r$a, x, max_na_frac = 1)
  tg <- test_dm(reg, ~ group, "groupB"); tg <- tg[tg$set == "target", ]
  ci <- stats::qbinom(c(0.005, 0.995), nrow(tg), 0.05) / nrow(tg)
  expect_gte(mean(tg$P < 0.05), ci[1]); expect_lte(mean(tg$P < 0.05), ci[2])
})

test_that("target statistics are identical whatever the background", {
  r <- setup_sim()
  xa <- define_sets(r$a, r$s$corsivs, background = "array")
  xu <- define_sets(r$a, r$s$corsivs, background = user_bg(r))
  ra <- test_dm(r$a, ~ group, "groupB", level = "cpg", sets = xa)
  ru <- test_dm(r$a, ~ group, "groupB", level = "cpg", sets = xu)
  ta <- ra[ra$set == "target", ]; tu <- ru[ru$set == "target", ]
  expect_identical(ta$id, tu$id)
  expect_equal(ta$P, tu$P); expect_equal(ta$FDR, tu$FDR)
  expect_true(all(ta$region_id %in% targetRegions(xa)$region_id))
  expect_true(all(ta$n_cpg == 1L))
})

test_that("bad designs give clear errors; NA covariates drop samples", {
  r <- setup_sim()
  expect_error(test_dm(r$a, ~ group, "groupB"), "summarise_regions")
  expect_error(test_dm(r$a, ~ group, "groupB", level = "cpg"), "sets")
  expect_error(test_dm(r$a, ~ group, "HCCCase", level = "cpg", sets = r$x), "groupB")
  expect_error(test_dm(r$a, ~ nope, "x", level = "cpg", sets = r$x), "nope")
  a2 <- r$a; a2$dup <- a2$group                                            # Review Focus 3
  expect_error(test_dm(a2, ~ group + dup, "groupB", level = "cpg", sets = r$x), "full rank")
  a3 <- r$a; a3$age[1:2] <- NA
  expect_message(res <- test_dm(a3, ~ group + age, "groupB", level = "cpg", sets = r$x), "2 samples dropped")
  expect_identical(metadata(res)$n_samples, 18L)
})

test_that("an all-NA row returns NA statistics instead of erroring", {  # Review Focus 2
  r <- setup_sim()
  a <- r$a
  i <- which(IRanges::overlapsAny(rowRanges(a), targetRegions(r$x)))[1]
  assay(a)[i, ] <- NA
  res <- test_dm(a, ~ group, "groupB", level = "cpg", sets = r$x)
  row <- res[res$id == rownames(a)[i], ]
  expect_identical(nrow(row), 1L)
  expect_true(is.na(row$P)); expect_true(is.na(row$delta_beta))
  expect_false(all(is.na(res$P)))
})

test_that("provenance metadata records every setting", {
  r <- setup_sim()
  res <- test_dm(r$a, ~ group, "groupB", level = "cpg", sets = r$x)
  md <- metadata(res)
  expect_true(all(c("CoRSIVtools_version", "corsiv_version", "manifest", "background_mode", "n_draws",
                    "seed", "flank", "min_cpgs", "fun", "max_na_frac", "design", "coef", "level",
                    "n_samples", "R_version", "Bioc_version", "timestamp") %in% names(md)))
  expect_identical(md$background_mode, "array")
  expect_identical(md$design, "~group")
  expect_identical(md$manifest, "simulated")
  expect_identical(md$corsiv_version, "sim")
})
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `Rscript -e 'devtools::test(filter = "test_dm")'`
Expected: FAIL with `could not find function "test_dm"`.

- [ ] **Step 3: Implement `R/provenance.R`**

```r
.provenance <- function(md, set_params, design, coef, level, n_samples) {
  list(
    CoRSIVtools_version = as.character(utils::packageVersion("CoRSIVtools")),
    corsiv_version = set_params$corsiv_version,
    manifest = if (is.null(md$manifest$source)) NA_character_ else md$manifest$source,
    background_mode = set_params$mode, n_draws = set_params$n_draws, seed = set_params$seed,
    flank = set_params$flank, min_cpgs = set_params$min_cpgs,
    fun = if (is.null(md$summarise$fun)) NA_character_ else md$summarise$fun,
    max_na_frac = if (is.null(md$summarise$max_na_frac)) NA_real_ else md$summarise$max_na_frac,
    design = paste(deparse(design), collapse = ""), coef = coef, level = level, n_samples = n_samples,
    R_version = R.version.string,
    Bioc_version = tryCatch(as.character(utils::packageVersion("BiocVersion")),
                            error = function(e) NA_character_),
    timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
}
```

- [ ] **Step 4: Implement `R/test_dm.R`**

```r
#' Differential methylation of CoRSIVs (and background) with limma
#'
#' Fits limma on M-values separately for the target set and each background
#' draw, so target moderation and FDR never depend on the background.
#' `delta_beta` is the same design fitted on the beta scale.
#'
#' @param x Region-level SE from [summarise_regions()] (`level = "corsiv"`)
#'   or annotated SE from [annotate_probes()] (`level = "cpg"`).
#' @param design One-sided formula using `colData(x)` columns.
#' @param coef Name of the design-matrix column to test, e.g. `"groupB"`.
#' @param level `"corsiv"` or `"cpg"`.
#' @param sets [define_sets()] output; required for `level = "cpg"`.
#' @return `DataFrame`; see the spec for the column schema. Provenance in
#'   `metadata()`.
#' @export
test_dm <- function(x, design, coef, level = c("corsiv", "cpg"), sets = NULL) {
  level <- match.arg(level)
  if (level == "corsiv") {
    if (!identical(metadata(x)$level, "region")) {
      stop("level = 'corsiv' needs the output of summarise_regions()", call. = FALSE)
    }
    rr <- rowRanges(x)
    info <- data.frame(id = rr$region_id, region_id = rr$region_id, set = rr$set, draw = rr$draw,
                       chr = as.character(seqnames(rr)), start = start(rr), end = end(rr),
                       n_cpg = rr$n_cpg, source_corsivs = rr$source_corsivs)
    b <- assay(x, "beta")
    set_params <- metadata(x)$sets
  } else {
    if (!is(sets, "CoRSIVSets")) stop("level = 'cpg' needs sets = define_sets(...)", call. = FALSE)
    .check_annotated(x)
    built <- .cpg_rows(x, sets)
    info <- built$info; b <- built$b
    set_params <- c(sets@params, mode = sets@mode)
  }
  cd <- as.data.frame(colData(x))
  vars <- all.vars(design)
  miss <- setdiff(vars, names(cd))
  if (length(miss)) stop("design variables not in colData: ", paste(miss, collapse = ", "), call. = FALSE)
  keep <- stats::complete.cases(cd[, vars, drop = FALSE])
  if (any(!keep)) message(sprintf("%d samples dropped for missing design variables", sum(!keep)))
  cd <- cd[keep, , drop = FALSE]; b <- b[, keep, drop = FALSE]
  mm <- stats::model.matrix(design, data = cd)
  if (!coef %in% colnames(mm)) {
    stop(sprintf("coef '%s' is not a design-matrix column; available: %s", coef,
                 paste(colnames(mm), collapse = ", ")), call. = FALSE)
  }
  q <- qr(mm)
  if (q$rank < ncol(mm)) {
    stop("design matrix is not of full rank; aliased column(s): ",
         paste(colnames(mm)[q$pivot[(q$rank + 1):ncol(mm)]], collapse = ", "), call. = FALSE)
  }
  grp <- paste(info$set, info$draw)
  blocks <- split(seq_len(nrow(info)), factor(grp, levels = unique(grp)))
  res <- do.call(rbind, lapply(blocks, function(i) .fit_block(b[i, , drop = FALSE], mm, coef, info[i, , drop = FALSE])))
  rownames(res) <- NULL
  res <- DataFrame(res, check.names = FALSE)
  metadata(res) <- .provenance(metadata(x), set_params, design, coef, level, sum(keep))
  res
}

.fit_block <- function(b, mm, coef, info) {
  stat_cols <- c("delta_beta", "logFC_M", "AveM", "t", "P", "FDR")
  for (cc in stat_cols) info[[cc]] <- NA_real_
  info$sd_beta <- apply(b, 1, stats::sd, na.rm = TRUE)
  ok <- rowSums(!is.na(b)) > ncol(mm)          # rows with residual df > 0; others stay NA
  if (!any(ok)) return(info)
  bo <- b[ok, , drop = FALSE]
  bc <- pmin(pmax(bo, 1e-6), 1 - 1e-6)
  M <- log2(bc / (1 - bc))
  fit <- limma::eBayes(limma::lmFit(M, mm), trend = TRUE, robust = TRUE)
  tt <- limma::topTable(fit, coef = coef, number = Inf, sort.by = "none")
  fb <- limma::lmFit(bo, mm)
  info$delta_beta[ok] <- unname(fb$coefficients[, coef])
  info$logFC_M[ok] <- tt$logFC; info$AveM[ok] <- tt$AveExpr; info$t[ok] <- tt$t
  info$P[ok] <- tt$P.Value; info$FDR[ok] <- tt$adj.P.Val
  info
}

.cpg_rows <- function(x, sets) {
  cpg <- rowRanges(x); b <- assay(x, "beta")
  one <- function(regions) {
    h <- findOverlaps(cpg, regions)
    i <- queryHits(h); r <- regions[subjectHits(h)]
    list(b = b[i, , drop = FALSE],
         info = data.frame(id = rownames(x)[i], region_id = r$region_id, set = r$set, draw = r$draw,
                           chr = as.character(seqnames(cpg))[i], start = start(cpg)[i],
                           end = end(cpg)[i], n_cpg = 1L, source_corsivs = r$source_corsivs))
  }
  t <- one(sets@target); g <- one(sets@background)
  list(b = rbind(t$b, g$b), info = rbind(t$info, g$info))
}
```

- [ ] **Step 5: Run the tests and confirm they pass**

Run: `Rscript -e 'devtools::document(); devtools::test(filter = "test_dm")'`
Expected: PASS. Rows without enough non-NA samples to leave residual degrees of freedom are never sent to limma; they come back with NA statistics. FDR is computed over the testable rows only.

- [ ] **Step 6: Commit**

```bash
git add R/test_dm.R R/provenance.R tests/testthat/test-test_dm.R NAMESPACE man/
git commit -m "Add test_dm(): per-block limma at CoRSIV and CpG level with full provenance"
```

---

### Task 10: `compare_to_background()`

**Files:**
- Create: `R/compare.R`, `tests/testthat/test-compare.R`

**Interfaces:**
- Consumes: the `test_dm()` result (columns `set, draw, P, FDR, delta_beta, sd_beta`).
- Produces: `compare_to_background(res, p = 0.05, fdr = 0.05)` → a named list of `DataFrame`s:
  - `enrichment`: `threshold_P, target_n, target_n_sig, target_frac_sig, bg_mean_frac_sig, bg_sd_frac_sig, empirical_p, fisher_or, fisher_p, target_n_fdr, bg_mean_n_fdr, n_draws`;
  - `per_draw`: `draw, n, n_sig, n_fdr, median_sd, frac_sig`;
  - `distributions`: `quantile, target_abs_delta_beta, bg_abs_delta_beta, target_P, bg_P`;
  - `qq`: `set, expected, observed`;
  - `variability`: `target_median_sd, bg_mean_median_sd, empirical_p`.

  `empirical_p` is NA when there are fewer than 20 draws.

- [ ] **Step 1: Write the failing tests**

`tests/testthat/test-compare.R`:
```r
fake_res <- function(n_draws = 20) {
  tg <- data.frame(set = "target", draw = NA_integer_, P = c(rep(0.01, 30), rep(0.5, 70)),
                   FDR = c(rep(0.04, 10), rep(0.5, 90)), delta_beta = 0.1, sd_beta = 0.1)
  bg <- do.call(rbind, lapply(seq_len(n_draws), function(d)
    data.frame(set = "background", draw = d, P = c(rep(0.01, 5), rep(0.5, 95)),
               FDR = 0.9, delta_beta = 0.01, sd_beta = 0.03)))
  S4Vectors::DataFrame(rbind(tg, bg))
}

test_that("enrichment, empirical p and Fisher test are computed correctly", {
  cmp <- compare_to_background(fake_res())
  e <- cmp$enrichment
  expect_identical(e$target_n_sig, 30L)
  expect_equal(e$target_frac_sig, 0.30)
  expect_equal(e$bg_mean_frac_sig, 0.05)
  expect_equal(e$empirical_p, 1 / 21)
  ft <- stats::fisher.test(matrix(c(30, 70, 100, 1900), 2))
  expect_equal(e$fisher_p, ft$p.value)
  expect_identical(e$target_n_fdr, 10L)
  expect_identical(nrow(cmp$per_draw), 20L)
  expect_equal(cmp$variability$target_median_sd, 0.1)
  expect_equal(cmp$variability$empirical_p, 1 / 21)
  expect_identical(sort(unique(cmp$qq$set)), c("background", "target"))
  expect_identical(cmp$distributions$quantile, c(0.05, 0.25, 0.5, 0.75, 0.95))
})

test_that("empirical p is NA with fewer than 20 draws; missing background errors", {
  expect_true(is.na(compare_to_background(fake_res(1))$enrichment$empirical_p))
  r <- fake_res(); r <- r[r$set == "target", ]
  expect_error(compare_to_background(r), "no background rows")
})
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `Rscript -e 'devtools::test(filter = "compare")'`
Expected: FAIL with `could not find function "compare_to_background"`.

- [ ] **Step 3: Implement `R/compare.R`**

```r
#' Compare CoRSIV results with the background
#'
#' @param res Output of [test_dm()].
#' @param p Nominal P threshold for enrichment.
#' @param fdr FDR threshold reported alongside.
#' @return Named list of `DataFrame`s: `enrichment`, `per_draw`,
#'   `distributions`, `qq`, `variability`.
#' @export
compare_to_background <- function(res, p = 0.05, fdr = 0.05) {
  res <- as.data.frame(res)
  tg <- res[res$set == "target", , drop = FALSE]
  bg <- res[res$set == "background", , drop = FALSE]
  if (!nrow(bg)) stop("res has no background rows; run define_sets() with a background", call. = FALSE)
  per_draw <- do.call(rbind, lapply(split(bg, bg$draw), function(d) data.frame(
    draw = d$draw[1], n = sum(!is.na(d$P)), n_sig = sum(d$P < p, na.rm = TRUE),
    n_fdr = sum(d$FDR < fdr, na.rm = TRUE), median_sd = stats::median(d$sd_beta, na.rm = TRUE))))
  per_draw$frac_sig <- per_draw$n_sig / per_draw$n
  nd <- nrow(per_draw)
  emp <- function(stat, obs) if (nd >= 20) (1 + sum(stat >= obs)) / (1 + nd) else NA_real_
  t_n <- sum(!is.na(tg$P)); t_sig <- sum(tg$P < p, na.rm = TRUE); t_frac <- t_sig / t_n
  b_n <- sum(!is.na(bg$P)); b_sig <- sum(bg$P < p, na.rm = TRUE)
  ft <- stats::fisher.test(matrix(c(t_sig, t_n - t_sig, b_sig, b_n - b_sig), 2))
  enrichment <- DataFrame(threshold_P = p, target_n = t_n, target_n_sig = t_sig, target_frac_sig = t_frac,
                          bg_mean_frac_sig = mean(per_draw$frac_sig),
                          bg_sd_frac_sig = if (nd > 1) stats::sd(per_draw$frac_sig) else NA_real_,
                          empirical_p = emp(per_draw$frac_sig, t_frac),
                          fisher_or = unname(ft$estimate), fisher_p = ft$p.value,
                          target_n_fdr = sum(tg$FDR < fdr, na.rm = TRUE),
                          bg_mean_n_fdr = mean(per_draw$n_fdr), n_draws = nd)
  probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  qf <- function(v) unname(stats::quantile(v, probs, na.rm = TRUE))
  distributions <- DataFrame(quantile = probs, target_abs_delta_beta = qf(abs(tg$delta_beta)),
                             bg_abs_delta_beta = qf(abs(bg$delta_beta)), target_P = qf(tg$P), bg_P = qf(bg$P))
  t_sd <- stats::median(tg$sd_beta, na.rm = TRUE)
  variability <- DataFrame(target_median_sd = t_sd, bg_mean_median_sd = mean(per_draw$median_sd),
                           empirical_p = emp(per_draw$median_sd, t_sd))
  list(enrichment = enrichment, per_draw = DataFrame(per_draw),
       distributions = distributions,
       qq = DataFrame(rbind(.qq(tg$P, "target"), .qq(bg$P, "background"))),
       variability = variability)
}

.qq <- function(p, set) {
  p <- sort(p[!is.na(p)])
  data.frame(set = rep(set, length(p)), expected = -log10(stats::ppoints(length(p))), observed = -log10(p))
}
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `Rscript -e 'devtools::document(); devtools::test(filter = "compare")'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/compare.R tests/testthat/test-compare.R NAMESPACE man/
git commit -m "Add compare_to_background(): enrichment, distributions, QQ and variability vs background"
```

---

### Task 11: Plots

**Files:**
- Create: `R/plots.R`, `tests/testthat/test-plots.R`

**Interfaces:**
- Consumes: the `test_dm()` result, the region SE, and the `compare_to_background()` output.
- Produces `ggplot` objects:
  - `plot_manhattan(res, set = "target", fdr = 0.05)`
  - `plot_volcano(res, set = "target", fdr = 0.05, delta = NULL)`
  - `plot_corsiv(x, id, group)`
  - `plot_background(cmp, which = c("enrichment", "variability"))`
  - internal `.theme_corsiv()`

- [ ] **Step 1: Write the failing tests**

`tests/testthat/test-plots.R`:
```r
pl <- local({
  s <- simulate_corsiv_data()
  a <- suppressMessages(annotate_probes(s$se, s$manifest))
  x <- define_sets(a, s$corsivs, background = s$controls)
  reg <- summarise_regions(a, x, max_na_frac = 1)
  res <- test_dm(reg, ~ group, "groupB")
  list(reg = reg, res = res, cmp = compare_to_background(res), dm = s$dm_ids)
})

test_that("all plots return buildable ggplot objects", {
  for (g in list(plot_manhattan(pl$res), plot_volcano(pl$res, delta = 0.05),
                 plot_corsiv(pl$reg, pl$dm[1], "group"),
                 plot_background(pl$cmp, "enrichment"), plot_background(pl$cmp, "variability"))) {
    expect_s3_class(g, "ggplot")
    expect_silent(ggplot2::ggplot_build(g))
  }
})

test_that("plot_corsiv errors on an unknown region", {
  expect_error(plot_corsiv(pl$reg, "chr1:1-2", "group"), "not found")
})
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `Rscript -e 'devtools::test(filter = "plots")'`
Expected: FAIL with `could not find function "plot_manhattan"`.

- [ ] **Step 3: Implement `R/plots.R`**

```r
.col <- c(hyper = "#e74c3c", hypo = "#3498db", ns = "#bdc3c7", alt = "#7f8c8d")

.theme_corsiv <- function() {
  ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16),
                   axis.title = ggplot2::element_text(face = "bold", size = 15),
                   axis.text = ggplot2::element_text(size = 12),
                   legend.text = ggplot2::element_text(size = 11))
}

.subset_set <- function(res, set) {
  d <- as.data.frame(res); d <- d[d$set == set & !is.na(d$P), , drop = FALSE]
  if (!nrow(d)) stop("no rows with set == '", set, "'", call. = FALSE)
  d
}

#' Manhattan plot of CoRSIV results
#' @param res Output of [test_dm()].
#' @param set `"target"` or `"background"`.
#' @param fdr FDR threshold; significant points in red, dashed line at the
#'   largest P with FDR < `fdr`.
#' @return `ggplot`.
#' @export
plot_manhattan <- function(res, set = "target", fdr = 0.05) {
  d <- .subset_set(res, set)
  lv <- paste0("chr", c(1:22, "X", "Y"))
  d <- d[d$chr %in% lv, , drop = FALSE]
  d$chr <- factor(d$chr, levels = intersect(lv, unique(d$chr)))
  d <- d[order(d$chr, d$start), ]
  len <- tapply(as.numeric(d$end), d$chr, max)
  off <- stats::setNames(c(0, cumsum(len))[seq_along(len)], names(len))
  d$pos <- d$start + off[as.character(d$chr)]
  d$band <- factor(as.integer(d$chr) %% 2)
  d$sig <- d$FDR < fdr
  ax <- tapply(d$pos, d$chr, function(v) mean(range(v)))
  g <- ggplot2::ggplot(d, ggplot2::aes(.data$pos, -log10(.data$P))) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$band), size = 1.2, alpha = 0.6, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c(`0` = .col[["ns"]], `1` = .col[["alt"]])) +
    ggplot2::geom_point(data = d[d$sig, , drop = FALSE], colour = .col[["hyper"]], size = 2.5) +
    ggplot2::scale_x_continuous(breaks = ax, labels = sub("^chr", "", names(ax))) +
    ggplot2::labs(x = "Chromosome", y = expression(-log[10](P)),
                  title = sprintf("CoRSIV differential methylation (%s)", set)) +
    .theme_corsiv()
  if (any(d$sig)) {
    g <- g + ggplot2::geom_hline(yintercept = -log10(max(d$P[d$sig])), colour = "grey",
                                 linetype = "dashed", linewidth = 1, alpha = 0.7)
  }
  g
}

#' Volcano plot of CoRSIV results
#' @inheritParams plot_manhattan
#' @param delta Optional |delta_beta| reference lines.
#' @return `ggplot`.
#' @export
plot_volcano <- function(res, set = "target", fdr = 0.05, delta = NULL) {
  d <- .subset_set(res, set)
  d$cls <- ifelse(d$FDR < fdr & d$delta_beta > 0, "hyper",
                  ifelse(d$FDR < fdr & d$delta_beta < 0, "hypo", "ns"))
  d$cls <- factor(d$cls, levels = c("hyper", "hypo", "ns"))
  g <- ggplot2::ggplot(d, ggplot2::aes(.data$delta_beta, -log10(.data$P), colour = .data$cls)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$cls, alpha = .data$cls)) +
    ggplot2::scale_colour_manual(values = .col[c("hyper", "hypo", "ns")], drop = FALSE,
                                 labels = c("Hyper (FDR)", "Hypo (FDR)", "NS"), name = NULL) +
    ggplot2::scale_size_manual(values = c(hyper = 2.5, hypo = 2.5, ns = 1), guide = "none") +
    ggplot2::scale_alpha_manual(values = c(hyper = 0.85, hypo = 0.85, ns = 0.45), guide = "none") +
    ggplot2::geom_vline(xintercept = 0, colour = "black", linewidth = 0.5, alpha = 0.3) +
    ggplot2::labs(x = expression(Delta * beta), y = expression(-log[10](P)),
                  title = sprintf("CoRSIV volcano (%s)", set)) +
    .theme_corsiv()
  if (!is.null(delta)) {
    g <- g + ggplot2::geom_vline(xintercept = c(-delta, delta), colour = .col[["hyper"]],
                                 linetype = "dotted", alpha = 0.5)
  }
  g
}

#' Beta values of one CoRSIV by group
#' @param x Region-level SE from [summarise_regions()].
#' @param id Target `region_id`, e.g. `"chr1:629002-635000"`.
#' @param group `colData` column to group by.
#' @return `ggplot`.
#' @export
plot_corsiv <- function(x, id, group) {
  i <- which(rowRanges(x)$region_id == id & rowRanges(x)$set == "target")
  if (length(i) != 1L) stop("target region '", id, "' not found in x", call. = FALSE)
  d <- data.frame(beta = as.numeric(assay(x, "beta")[i, ]), group = factor(colData(x)[[group]]))
  ggplot2::ggplot(d, ggplot2::aes(.data$group, .data$beta)) +
    ggplot2::geom_boxplot(outlier.shape = NA, fill = "white") +
    ggplot2::geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.85, colour = .col[["hypo"]]) +
    ggplot2::labs(x = group, y = "Mean beta across CoRSIV CpGs", title = id) +
    .theme_corsiv()
}

#' Target vs background distribution across draws
#' @param cmp Output of [compare_to_background()].
#' @param which `"enrichment"` (fraction P < threshold) or `"variability"`
#'   (median SD of beta).
#' @return `ggplot`.
#' @export
plot_background <- function(cmp, which = c("enrichment", "variability")) {
  which <- match.arg(which)
  pd <- as.data.frame(cmp$per_draw)
  if (which == "enrichment") {
    v <- pd$frac_sig; obs <- cmp$enrichment$target_frac_sig
    xl <- sprintf("Fraction with P < %g", cmp$enrichment$threshold_P)
  } else {
    v <- pd$median_sd; obs <- cmp$variability$target_median_sd
    xl <- "Median SD of beta across samples"
  }
  ggplot2::ggplot(data.frame(v = v), ggplot2::aes(.data$v)) +
    ggplot2::geom_histogram(bins = max(5L, min(30L, length(v))), fill = .col[["ns"]], colour = "white") +
    ggplot2::geom_vline(xintercept = obs, colour = .col[["hyper"]], linewidth = 1.2) +
    ggplot2::labs(x = xl, y = "Background draws",
                  title = sprintf("CoRSIVs (red) vs %d background draw(s)", length(v))) +
    .theme_corsiv()
}
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `Rscript -e 'devtools::document(); devtools::test(filter = "plots")'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/plots.R tests/testthat/test-plots.R NAMESPACE man/
git commit -m "Add Manhattan, volcano, per-CoRSIV and background plots in lab figure style"
```

---

### Task 12: Vignettes, citation, CI, and a clean R CMD check plus BiocCheck

**Files:**
- Create: `vignettes/CoRSIVtools.Rmd`, `vignettes/methods.Rmd`, `inst/CITATION`, `NEWS.md`, `README.md`, `.github/workflows/R-CMD-check.yaml`, `.github/workflows/BiocCheck.yaml`, `.github/workflows/test-coverage.yaml`

**Interfaces:**
- Consumes: every exported function.
- Produces: a package that passes `R CMD check` with 0 errors and 0 warnings, and BiocCheck with 0 errors.

- [ ] **Step 1: Write `vignettes/CoRSIVtools.Rmd`**

````markdown
---
title: "CoRSIV-focused differential methylation with CoRSIVtools"
output: BiocStyle::html_document
vignette: >
  %\VignetteIndexEntry{CoRSIVtools}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, message = FALSE}
library(CoRSIVtools)
```

## Inputs

Real analyses supply a normalised beta matrix, a metadata table with one row per
array sample, and an hg38 probe manifest (BED or Illumina CSV, or `"EPICv1"` /
`"EPICv2"`). For EPIC-Xtra, request the manifest from the vendor. Here we use
simulated data on real chr21 sequence.

```{r data}
sim <- simulate_corsiv_data(n_samples = 20, n_corsivs = 40, n_dm = 10)
se  <- sim$se                                   # from_matrix(betas, meta, "sample_id")
se  <- annotate_probes(se, sim$manifest)        # read_manifest("path/to/manifest.bed")
```

## Target CoRSIVs and a matched background

```{r sets}
sets <- define_sets(se, corsivs = sim$corsivs, background = "matched", n_draws = 20, seed = 1)
sets
coverage_report(sets)
```

`background = "array"` uses every CpG outside the CoRSIVs (CpG level only);
a `GRanges` uses your own control regions.

## CoRSIV-level test

```{r test}
reg <- summarise_regions(se, sets)
res <- test_dm(reg, design = ~ group + age, coef = "groupB")
head(res[res$set == "target", ][order(res$P[res$set == "target"]), ])
```

## Compare with the background

```{r compare}
cmp <- compare_to_background(res)
cmp$enrichment
cmp$variability
```

## Plots

```{r plots, fig.width = 9, fig.height = 5}
plot_volcano(res)
plot_background(cmp, "enrichment")
plot_corsiv(reg, sim$dm_ids[1], "group")
```

## Exporting for bedtools

```{r bed}
f <- tempfile(fileext = ".bed")
write_bed(res[res$set == "target", ], f)
readLines(f, 3)
```

## Session info

```{r}
sessionInfo()
```
````

- [ ] **Step 2: Write `vignettes/methods.Rmd`**

````markdown
---
title: "CoRSIVtools statistical methods"
output: BiocStyle::html_document
vignette: >
  %\VignetteIndexEntry{CoRSIVtools methods}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

# Probe mapping and replicate probes

Probes are mapped to hg38 CpG positions with a user-supplied or built-in
manifest. Probes without coordinates are dropped and counted. EPICv2 and
EPIC-Xtra include replicate probes (suffixes such as `_TC11`, `_BC11`) that
target the same CpG; they are averaged per sample (ignoring NA) **only when
their coordinates are identical**. Probes sharing a `cg` base name but mapped to
different positions are kept separate, because a shared name does not
guarantee a shared locus.

# CoRSIV definition and summary

CoRSIVs are taken from CoRSIV 2.0 (merged population-specific calls, hg38;
Gunasekara et al. 2019 for the original atlas). Each CoRSIV is identified by its
1-based coordinates (`chr:start-end`); the original merged names are retained
but are not unique and are never used as keys. A CoRSIV is analysed if it
contains at least `min_cpgs` (default 3) distinct CpGs on the array. For sample
*j* and CoRSIV *r* with measured CpGs *C_rj*,

$$\bar\beta_{rj} = \frac{1}{|C_{rj}|} \sum_{c \in C_{rj}} \beta_{cj}, \quad |C_{rj}| \ge \texttt{min\_cpgs},$$

and NA otherwise. Regions with more than `max_na_frac` (default 0.2) NA samples
are removed. Because CpGs within a CoRSIV are, by definition, strongly
correlated across individuals, the mean is an efficient low-noise summary.

# Differential methylation model

For each region the summary is converted to an M-value,
$M = \log_2\{\beta/(1-\beta)\}$ with $\beta$ clamped to $[10^{-6}, 1-10^{-6}]$,
and modelled with a linear model $M_r = X\gamma_r + \varepsilon_r$, where $X$ is
the user's design (e.g. `~ group + age + sex + cell fractions`). Variances are
moderated with limma's empirical Bayes (`trend = TRUE` for the mean–variance
trend, `robust = TRUE` to down-weight hypervariable regions; Smyth 2004;
Phipson et al. 2016). M-values are used for inference because their variance is
approximately homogeneous across the methylation range (Du et al. 2010); the
reported effect size, $\Delta\beta$, is the same design coefficient fitted on
the beta scale.

**Multiple testing.** Benjamini–Hochberg FDR is computed within the target set
only (and separately within each background draw). The FDR therefore refers to
the family of CoRSIV tests, not to the whole array; this is a deliberate
restriction of the hypothesis family and must be reported as such.

**Assumptions and limitations.** Independent samples (no repeated measures in
v0.1); approximately normal residuals on the M scale; covariates enter only
through the design. CoRSIVs are enriched for methylation QTLs, so genotype can
confound group comparisons — include ancestry PCs where available. chrX
CoRSIVs are excluded by default because their methylation differs by sex and
requires sex-stratified models.

# Background comparison

**Matched control regions.** For each of `n_draws` draws (seed `seed + r`),
each analysed CoRSIV is paired with one random control region that has (1) the
same width, (2) the same number of array CpGs (±1 when ≥ 10), (3) the same
chromosome type, and (4) CpG density (CpGs per 100 bp, hg38) within ±20%
(minimum ±0.25). Controls are proposed by choosing a random array CpG outside
every CoRSIV ± `flank` bp and placing a window of the target's width at a
uniformly random offset containing it, then rejection-sampling until all four
criteria hold (up to `max_tries` proposals); controls within a draw do not
overlap. This follows the permutation-with-constraints approach of regioneR
(Gel et al. 2016), adapted to array coverage, which is the dominant determinant
of measurement precision for region summaries.

**Statistics.** The enrichment statistic is the fraction of regions with
P < 0.05. The empirical p-value is $(1 + \#\{r : f_r \ge f_{target}\})/(1 + R)$
over $R$ draws (Phipson & Smyth 2010), reported when $R \ge 20$. Fisher's exact
test of target versus pooled background is also given and is the primary test
for `background = "array"` (all non-CoRSIV CpGs). The same empirical procedure
is applied to the median between-individual SD of beta, which provides a
data-quality check: CoRSIVs are expected to be more variable between individuals
than matched control regions.

# References

- Du P, et al. (2010) Comparison of Beta-value and M-value methods for quantifying methylation levels by microarray analysis. *BMC Bioinformatics* 11:587.
- Gel B, et al. (2016) regioneR: an R/Bioconductor package for the association analysis of genomic regions based on permutation tests. *Bioinformatics* 32:289–291.
- Gunasekara CJ, et al. (2019) A genomic atlas of systemic interindividual epigenetic variation in humans. *Genome Biology* 20:105.
- Phipson B, Smyth GK (2010) Permutation P-values should never be zero. *Stat Appl Genet Mol Biol* 9:39.
- Phipson B, et al. (2016) Robust hyperparameter estimation protects against hypervariable genes and improves power to detect differential expression. *Ann Appl Stat* 10:946–963.
- Smyth GK (2004) Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. *Stat Appl Genet Mol Biol* 3:3.
````

- [ ] **Step 3: Write `inst/CITATION`, `NEWS.md` and `README.md`**

`inst/CITATION`:
```r
bibentry(bibtype = "Manual",
         title = "CoRSIVtools: CoRSIV-Focused DNA Methylation Analysis",
         author = c(person("Chathura J.", "Gunasekara")),
         organization = "Baylor College of Medicine",
         year = "2026", note = "R package version 0.1.0",
         url = "https://github.com/waterlandlab/CoRSIVtools")
bibentry(bibtype = "Article",
         title = "A genomic atlas of systemic interindividual epigenetic variation in humans",
         author = c(person("Chathura J.", "Gunasekara"), person("others")),
         journal = "Genome Biology", year = "2019", volume = "20", pages = "105",
         doi = "10.1186/s13059-019-1708-1",
         textVersion = paste("Gunasekara CJ, et al. (2019) A genomic atlas of systemic interindividual",
                             "epigenetic variation in humans. Genome Biology 20:105.",
                             "doi:10.1186/s13059-019-1708-1. (CoRSIV 2.0 citation forthcoming.)"))
```

`NEWS.md`:
```markdown
# CoRSIVtools 0.1.0

* First release: from_matrix(), read_manifest(), annotate_probes(),
  corsiv_regions() (CoRSIV 2.0, hg38), define_sets() with matched / array /
  user backgrounds, coverage_report(), summarise_regions(), test_dm() at
  CoRSIV and CpG level, compare_to_background(), plots, write_bed().
```

`README.md`:
````markdown
# CoRSIVtools

CoRSIV-focused DNA methylation analysis for any Illumina array: map probes to
CoRSIVs, summarise per CoRSIV, test differential methylation with limma, and
compare against matched control regions or the rest of the array.

```r
# install.packages("BiocManager")
BiocManager::install("waterlandlab/CoRSIVtools")
```

```r
library(CoRSIVtools)
se   <- from_matrix(betas, meta, sample_col = "Sample_ID")
se   <- annotate_probes(se, read_manifest("EPICv2"))      # or a manifest file path
sets <- define_sets(se, background = "matched", n_draws = 100)
res  <- test_dm(summarise_regions(se, sets), ~ group + age + sex, coef = "groupCase")
cmp  <- compare_to_background(res)
```

See `vignette("CoRSIVtools")` and `vignette("methods", package = "CoRSIVtools")`.
Licence: MIT © Baylor College of Medicine.
````

- [ ] **Step 4: Add the CI workflows**

`.github/workflows/R-CMD-check.yaml`:
```yaml
on:
  push: {branches: [main]}
  pull_request: {branches: [main]}
name: R-CMD-check
jobs:
  R-CMD-check:
    runs-on: ${{ matrix.os }}
    strategy:
      fail-fast: false
      matrix:
        os: [ubuntu-latest, macos-latest]
    env:
      GITHUB_PAT: ${{ secrets.GITHUB_TOKEN }}
      R_KEEP_PKG_SOURCE: yes
    steps:
      - uses: actions/checkout@v4
      - uses: r-lib/actions/setup-pandoc@v2
      - uses: r-lib/actions/setup-r@v2
        with: {r-version: release, use-public-rspm: true}
      - uses: r-lib/actions/setup-r-dependencies@v2
        with: {extra-packages: any::rcmdcheck, needs: check}
      - uses: r-lib/actions/check-r-package@v2
        with: {upload-snapshots: true, error-on: '"warning"'}
```

`.github/workflows/BiocCheck.yaml`:
```yaml
on:
  push: {branches: [main]}
  pull_request: {branches: [main]}
name: BiocCheck
jobs:
  BiocCheck:
    runs-on: ubuntu-latest
    env: {GITHUB_PAT: "${{ secrets.GITHUB_TOKEN }}"}
    steps:
      - uses: actions/checkout@v4
      - uses: r-lib/actions/setup-pandoc@v2
      - uses: r-lib/actions/setup-r@v2
        with: {r-version: release, use-public-rspm: true}
      - uses: r-lib/actions/setup-r-dependencies@v2
        with: {extra-packages: bioc::BiocCheck}
      - run: Rscript -e 'res <- BiocCheck::BiocCheck(".", `quit-with-status` = TRUE, `no-check-bioc-help` = TRUE)'
```

`.github/workflows/test-coverage.yaml`:
```yaml
on:
  push: {branches: [main]}
  pull_request: {branches: [main]}
name: test-coverage
jobs:
  test-coverage:
    runs-on: ubuntu-latest
    env: {GITHUB_PAT: "${{ secrets.GITHUB_TOKEN }}"}
    steps:
      - uses: actions/checkout@v4
      - uses: r-lib/actions/setup-r@v2
        with: {use-public-rspm: true}
      - uses: r-lib/actions/setup-r-dependencies@v2
        with: {extra-packages: any::covr, needs: coverage}
      - run: Rscript -e 'cov <- covr::package_coverage(); print(cov); stopifnot(covr::percent_coverage(cov) >= 80)'
```

`README.md` and `NEWS.md` ship with the package, so neither goes in `.Rbuildignore`.

- [ ] **Step 5: Run the full check locally**

```bash
Rscript -e 'devtools::document(); devtools::test()'
Rscript -e 'devtools::check(error_on = "warning")'
Rscript -e 'BiocCheck::BiocCheck(".", `no-check-bioc-help` = TRUE)'
Rscript -e 'cov <- covr::package_coverage(); print(covr::percent_coverage(cov))'
```
Expected:
- all tests pass;
- `R CMD check` gives 0 errors and 0 warnings (NOTEs allowed; record each one in the commit message);
- BiocCheck gives 0 errors (its WARNINGs and NOTEs are recorded for the 0.99.0 submission step);
- coverage is at least 80%.

If `check` warns about undocumented arguments or a missing `@return`, fix the roxygen comment and re-run. Don't suppress the warning.

- [ ] **Step 6: Commit**

```bash
git add vignettes/ inst/CITATION NEWS.md README.md .github/ DESCRIPTION NAMESPACE man/
git commit -m "Add vignettes (workflow + formal methods), citation, CI; package passes R CMD check"
```

---

### Task 13: Real-data smoke run on Addanki (no inference)

This task checks that v0.1 runs at the scale of a real array: 1M probes, 18,649 CoRSIV probes and 990,766 background probes. It runs **no** group comparison; the Addanki analysis plan defines that separately. Nothing from this task is committed to the package. The script lives in the Addanki repo's `scratch/` folder, which is gitignored.

**Files:**
- Create: `<Addanki repo>/scratch/smoke_corsivtools.R`, where the Addanki repo is `~/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/WORK/WaterlandLab/XTraCoRSIV_Addanki/Baylor_Addanki_10022026/Normalized_Data`.

**Interfaces:**
- Consumes: the installed CoRSIVtools 0.1.0, `data/raw/betas_Xtra_filtered.RData` and `data/raw/XTRACoRSIV1_SS_20156551_A1.bed`.

- [ ] **Step 1: Install the package locally**

Run: `Rscript -e 'devtools::install("~/code/CoRSIVtools", upgrade = "never")'`
Expected: `* DONE (CoRSIVtools)`.

- [ ] **Step 2: Write and run the smoke script**

```r
library(CoRSIVtools)
e <- new.env(); load("data/raw/betas_Xtra_filtered.RData", envir = e)
b <- e$betas_Xtra_filtered
meta <- data.frame(array_id = colnames(b), dummy = rep(c("A", "B"), length.out = ncol(b)))
t0 <- Sys.time()
se <- from_matrix(b, meta, "array_id"); rm(b, e); gc()
se <- annotate_probes(se, read_manifest("data/raw/XTRACoRSIV1_SS_20156551_A1.bed"))
cat("annotated:", nrow(se), "CpGs\n")
sets <- define_sets(se, background = "matched", n_draws = 2, seed = 1)
print(sets); print(coverage_report(sets)); print(matchDiagnostics(sets))
reg <- summarise_regions(se, sets)
cat("regions:", nrow(reg), " elapsed:", format(Sys.time() - t0), "\n")
xa <- define_sets(se, background = "array")
print(coverage_report(xa))
```
Run from the Addanki repo root: `Rscript scratch/smoke_corsivtools.R 2>&1 | tee scratch/smoke_corsivtools.log`

Expected, from the facts recorded in the spec:
- `annotated:` is about 1,009,415 minus the replicates collapsed by coordinate (rows = distinct coordinates);
- the corsiv row of the coverage report has `n_regions = 6349` (autosomes) and `n_ge3` close to 3,083 minus the chrX contribution;
- both draws have ≤5% unmatched;
- the background row of the `xa` coverage report has `n_regions` close to 990,766 minus the CpGs within 1 kb of a CoRSIV.

Record the elapsed time per matched draw. If one draw takes more than 5 minutes, open a NEWS item proposing a speed-up (cached CpG densities per anchor) for v0.2. Don't optimise in v0.1.

- [ ] **Step 3: Report**

Write a summary of the log in the hand-off message: counts, unmatched fraction and timing. Don't commit anything in either repo for this task.

---

## Release gate (not a task)

Pushing to the public `waterlandlab/CoRSIVtools` repo needs the user's explicit go-ahead, after the open item in spec §12.4 (open-source release sign-off) is resolved. When approved:

```bash
gh repo create waterlandlab/CoRSIVtools --public --source . --remote origin --push
```
