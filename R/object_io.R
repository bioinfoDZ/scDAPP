# RDS and qs2 object IO. qs2 compression level is ignored for RDS.

#' Normalize a saving algorithm name.
#' @param saving_algorithm \code{"RDS"} or \code{"qs2"}.
#' @return \code{"RDS"} or \code{"qs2"}.
#' @keywords internal
normalize_saving_algorithm <- function(saving_algorithm = "RDS") {
  if (is.null(saving_algorithm) || length(saving_algorithm) < 1L ||
      !nzchar(as.character(saving_algorithm)[1])) {
    return("RDS")
  }
  x <- tolower(trimws(as.character(saving_algorithm)[1]))
  if (x %in% c("rds", "saverds")) {
    return("RDS")
  }
  if (x %in% c("qs2", "qs")) {
    return("qs2")
  }
  stop("saving_algorithm must be \"RDS\" or \"qs2\".", call. = FALSE)
}

#' File extension for a saving algorithm.
#' @keywords internal
scdapp_object_ext <- function(saving_algorithm = "RDS") {
  if (identical(normalize_saving_algorithm(saving_algorithm), "qs2")) ".qs2" else ".rds"
}

#' Path with the extension for \code{saving_algorithm}.
#'
#' An existing \code{.rds} or \code{.qs2} suffix is replaced.
#'
#' @param path File path, with or without an object extension.
#' @param saving_algorithm \code{"RDS"} or \code{"qs2"}.
#' @return Character path.
#' @export
scdapp_object_path <- function(path, saving_algorithm = "RDS") {
  path <- as.character(path)[1]
  if (is.na(path) || !nzchar(path)) {
    stop("path is empty.", call. = FALSE)
  }
  path <- sub("\\.(rds|qs2)$", "", path, ignore.case = TRUE)
  paste0(path, scdapp_object_ext(saving_algorithm))
}

#' Locate an object file that may be \code{.rds} or \code{.qs2}.
#'
#' An existing path is returned as given. Otherwise the stem is tried with
#' both extensions, preferring \code{saving_algorithm} when both exist.
#'
#' @param path Path or stem.
#' @param saving_algorithm Preferred format when both files exist.
#' @return Existing path.
#' @export
find_scdapp_object <- function(path, saving_algorithm = "RDS") {
  path <- as.character(path)[1]
  if (!is.na(path) && nzchar(path) && file.exists(path)) {
    return(path)
  }
  stem <- sub("\\.(rds|qs2)$", "", path, ignore.case = TRUE)
  prefer <- scdapp_object_path(stem, saving_algorithm)
  other <- if (grepl("\\.qs2$", prefer, ignore.case = TRUE)) {
    paste0(stem, ".rds")
  } else {
    paste0(stem, ".qs2")
  }
  if (file.exists(prefer)) {
    return(prefer)
  }
  if (file.exists(other)) {
    return(other)
  }
  stop(
    "Object file not found. Looked for: ",
    paste(unique(c(path, prefer, other)), collapse = ", "),
    call. = FALSE
  )
}

.require_qs2 <- function() {
  if (!requireNamespace("qs2", quietly = TRUE)) {
    stop(
      "qs2 is required for this file. Install it with install.packages(\"qs2\").",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Write an R object as RDS or qs2.
#'
#' @param object Object to save.
#' @param path Destination path. The extension follows \code{saving_algorithm}.
#' @param saving_algorithm \code{"RDS"} (default) or \code{"qs2"}.
#' @param compress_level qs2 zstd level. Ignored for RDS. Default 3.
#' @param nthreads Threads for qs2 compression. Use 1 inside an existing
#'   parallel worker. Ignored for RDS.
#' @return The path written, invisibly.
#' @export
write_scdapp_object <- function(
    object,
    path,
    saving_algorithm = "RDS",
    compress_level = 3L,
    nthreads = 1L
) {
  alg <- normalize_saving_algorithm(saving_algorithm)
  path <- scdapp_object_path(path, alg)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  nthreads <- max(1L, as.integer(nthreads)[1])
  if (is.na(nthreads)) {
    nthreads <- 1L
  }
  if (identical(alg, "qs2")) {
    .require_qs2()
    level <- as.integer(compress_level)[1]
    if (is.na(level)) {
      level <- 3L
    }
    qs2::qs_save(
      object,
      file = path,
      compress_level = level,
      nthreads = nthreads
    )
  } else {
    saveRDS(object, file = path)
  }
  invisible(path)
}

#' Read an RDS or qs2 object.
#'
#' The format is chosen from the file extension. A missing path is retried
#' with the other extension.
#'
#' @param path File path.
#' @param nthreads Threads for qs2 decompression. Ignored for RDS.
#' @return The stored object.
#' @export
read_scdapp_object <- function(path, nthreads = 1L) {
  path <- find_scdapp_object(path)
  nthreads <- max(1L, as.integer(nthreads)[1])
  if (is.na(nthreads)) {
    nthreads <- 1L
  }
  if (grepl("\\.qs2$", path, ignore.case = TRUE)) {
    .require_qs2()
    qs2::qs_read(path, nthreads = nthreads)
  } else {
    readRDS(path)
  }
}
