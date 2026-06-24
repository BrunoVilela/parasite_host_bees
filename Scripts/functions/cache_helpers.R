# Modular cache helpers for the parasite-host niche workflow.
#
# The workflow prefers qs when the package is available because it is fast for R
# objects. The current R installation does not provide qs, so every helper falls
# back to base RDS files without changing calling code.

cache_safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  gsub("^_+|_+$", "", x)
}

cache_backend_extension <- function() {
  if (requireNamespace("qs", quietly = TRUE)) ".qs" else ".rds"
}

pack_cache_value <- function(x) {
  if (requireNamespace("terra", quietly = TRUE) && inherits(x, "SpatRaster")) {
    return(terra::wrap(x))
  }
  if (is.list(x) && !is.data.frame(x)) {
    x <- lapply(x, pack_cache_value)
  }
  x
}

unpack_cache_value <- function(x) {
  if (requireNamespace("terra", quietly = TRUE) && inherits(x, "PackedSpatRaster")) {
    return(terra::unwrap(x))
  }
  if (is.list(x) && !is.data.frame(x)) {
    x <- lapply(x, unpack_cache_value)
  }
  x
}

cache_file_path <- function(cache_dir, name) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  safe_name <- cache_safe_name(name)
  preferred <- file.path(cache_dir, paste0(safe_name, cache_backend_extension()))
  alternatives <- file.path(cache_dir, paste0(safe_name, c(".qs", ".rds")))

  existing <- alternatives[file.exists(alternatives)]
  if (length(existing) > 0) {
    return(existing[[1]])
  }
  preferred
}

cache_read <- function(path) {
  if (grepl("\\.qs$", path) && requireNamespace("qs", quietly = TRUE)) {
    return(unpack_cache_value(qs::qread(path)))
  }
  unpack_cache_value(readRDS(path))
}

cache_write <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  value <- pack_cache_value(value)
  if (grepl("\\.qs$", path) && requireNamespace("qs", quietly = TRUE)) {
    qs::qsave(value, path)
  } else {
    saveRDS(value, path)
  }
  invisible(path)
}

input_fingerprint <- function(paths) {
  paths <- normalizePath(paths, mustWork = FALSE)
  exists <- file.exists(paths)
  info <- file.info(paths)
  md5 <- rep(NA_character_, length(paths))
  md5[exists] <- unname(tools::md5sum(paths[exists]))

  data.frame(
    path = paths,
    exists = exists,
    size = unname(info$size),
    mtime = as.character(info$mtime),
    md5 = md5,
    stringsAsFactors = FALSE
  )
}

cache_eval <- function(cache_dir, name, key, code, use_cache = TRUE) {
  path <- cache_file_path(cache_dir, name)

  if (isTRUE(use_cache) && file.exists(path)) {
    cached <- cache_read(path)
    if (is.list(cached) && identical(cached$key, key)) {
      return(list(value = cached$value, path = path, status = "hit"))
    }
  }

  value <- eval(substitute(code), envir = parent.frame())
  cache_write(
    list(
      key = key,
      value = value,
      created_at = as.character(Sys.time())
    ),
    path
  )
  list(value = value, path = path, status = "miss")
}
