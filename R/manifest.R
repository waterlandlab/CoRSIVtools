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
#' @examples
#' f <- tempfile(fileext = ".bed")
#' writeLines(c("chr21\t30000100\t30000101\tcg00000001",
#'              "chr21\t30000200\t30000201\tcg00000002"), f)
#' read_manifest(f)
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
