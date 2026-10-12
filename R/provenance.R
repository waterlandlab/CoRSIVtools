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
