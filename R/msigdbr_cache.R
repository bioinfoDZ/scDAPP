#' MSigDB gs_subcat used for cluster-marker cell-type ORA (currently C8)
#'
#' @return Single normalized `gs_subcat` string.
#' @keywords internal
.msigdb_celltype_gs_subcat <- function() {
  "C8"
}

#' Default MSigDB pathway categories used by scDAPP
#'
#' @return Character vector of `gs_subcat` values after normalization.
#' @keywords internal
.default_msigdbr_pwaycats <- function() {
  pwaycats <- c(
    "HALLMARK", "GO_BP", "GO_MF", "GO_CC",
    "CP_REACTOME", "CP_KEGG", "TFT_GTRD", "TFT_TFT_Legacy",
    .msigdb_celltype_gs_subcat()
  )
  pwaycats <- gsub(":", "_", pwaycats)
  names(pwaycats) <- pwaycats
  pwaycats
}

#' @keywords internal
.msigdbr_species_slug <- function(species) {
  slug <- gsub("[^A-Za-z0-9._-]+", "_", species)
  slug <- gsub("_+", "_", slug)
  slug <- gsub("^_|_$", "", slug)
  if (!nzchar(slug)) {
    stop("Invalid species string for msigdbr cache naming.", call. = FALSE)
  }
  slug
}

#' @keywords internal
.msigdbr_cache_dir_writable <- function(dir) {
  if (is.null(dir) || !nzchar(dir)) return(FALSE)
  ok <- tryCatch({
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    testfile <- file.path(dir, paste0(".scdapp_write_test_", Sys.getpid()))
    writeLines("ok", testfile)
    unlink(testfile)
    TRUE
  }, error = function(e) FALSE)
  isTRUE(ok)
}

#' Resolve a writable directory for cached msigdbr pathway tables
#'
#' Uses, in order: user `cache_dir`; `XDG_CACHE_HOME/scDAPP` when set;
#' [tools::R_user_dir()] cache for scDAPP; then `fallback_dir` with a warning.
#'
#' @param cache_dir Optional user-provided cache root (typically .../scDAPP).
#' @param fallback_dir Directory used when the preferred cache is not writable
#'   (pipeline `outdir_int`).
#' @return Character path to a writable cache directory.
#' @export
resolve_msigdbr_cache_dir <- function(cache_dir = NULL, fallback_dir = NULL) {
  candidates <- list()

  if (!is.null(cache_dir) && nzchar(cache_dir)) {
    candidates <- c(candidates, list(list(path = cache_dir, label = "user cache_dir")))
  }

  xdg <- Sys.getenv("XDG_CACHE_HOME", unset = NA_character_)
  if (!is.na(xdg) && nzchar(xdg)) {
    candidates <- c(candidates, list(list(
      path = file.path(xdg, "scDAPP"),
      label = "XDG_CACHE_HOME/scDAPP"
    )))
  }

  candidates <- c(candidates, list(list(
    path = tools::R_user_dir("scDAPP", which = "cache"),
    label = "R_user_dir cache"
  )))

  for (cand in candidates) {
    if (.msigdbr_cache_dir_writable(cand$path)) {
      return(normalizePath(cand$path, winslash = "/", mustWork = FALSE))
    }
  }

  if (!is.null(fallback_dir) && nzchar(fallback_dir)) {
    fb <- file.path(fallback_dir, "msigdbr_cache")
    if (.msigdbr_cache_dir_writable(fb)) {
      warning(
        "Could not write msigdbr cache to user/XDG/R cache locations; ",
        "using fallback directory: ", fb,
        call. = FALSE
      )
      return(normalizePath(fb, winslash = "/", mustWork = FALSE))
    }
  }

  stop(
    "Could not find a writable msigdbr cache directory. ",
    "Set msigdbr_cache_dir to a writable path or ensure outdir_int is writable.",
    call. = FALSE
  )
}

#' Build the cache file path for prepared msigdbr pathways
#'
#' @param species Species passed to [msigdbr::msigdbr()].
#' @param cache_dir Optional cache root; resolved via [resolve_msigdbr_cache_dir()]
#'   when NULL.
#' @param fallback_dir Passed to [resolve_msigdbr_cache_dir()].
#' @param msigdbr_version msigdbr package version string. Defaults to installed version.
#' @param db_version MSigDB version from `msigdbr` output (required for exact path).
#' @return Character path to the cache RDS file.
#' @export
msigdbr_cache_path <- function(species,
                                 cache_dir = NULL,
                                 fallback_dir = NULL,
                                 msigdbr_version = NULL,
                                 db_version = NULL) {
  if (is.null(msigdbr_version)) {
    msigdbr_version <- as.character(utils::packageVersion("msigdbr"))
  }
  dir <- resolve_msigdbr_cache_dir(cache_dir = cache_dir, fallback_dir = fallback_dir)
  slug <- .msigdbr_species_slug(species)
  ver <- gsub("[^0-9A-Za-z._-]+", "_", msigdbr_version)
  if (is.null(db_version) || !nzchar(db_version)) {
    pattern <- paste0("^msigdbr_", slug, "_", ver, "_.*_prepared\\.rds$")
    hits <- list.files(dir, pattern = pattern, full.names = TRUE)
    if (length(hits) == 0) return(NA_character_)
    hits[order(file.info(hits)$mtime, decreasing = TRUE)[1]]
  } else {
    dbv <- gsub("[^0-9A-Za-z._-]+", "_", db_version)
    file.path(dir, paste0("msigdbr_", slug, "_", ver, "_", dbv, "_prepared.rds"))
  }
}

#' @keywords internal
.validate_msigdbr_raw <- function(pathways) {
  cols <- colnames(pathways)
  has_legacy <- all(c("gs_cat", "gs_subcat") %in% cols)
  has_modern <- all(c("gs_collection", "gs_subcollection") %in% cols)
  if (!all(c("gene_symbol", "gs_name") %in% cols)) {
    stop(
      "msigdbr output is missing required columns gene_symbol and/or gs_name. ",
      "Found: ", paste(cols, collapse = ", "),
      call. = FALSE
    )
  }
  if (!has_legacy && !has_modern) {
    stop(
      "msigdbr output is missing collection columns (gs_cat/gs_subcat or ",
      "gs_collection/gs_subcollection). Found: ", paste(cols, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(pathways)
}

#' @keywords internal
.normalize_msigdbr_pathways <- function(pathways, pwaycats = NULL) {
  if (is.null(pwaycats)) pwaycats <- .default_msigdbr_pwaycats()

  msigdbrcolnames <- colnames(pathways)
  if (any(!c("gs_subcat", "gs_cat") %in% msigdbrcolnames)) {
    pathways$gs_cat <- pathways$gs_collection
    pathways$gs_subcat <- pathways$gs_subcollection
    pathways$gs_subcat <- gsub(":", "_", pathways$gs_subcat)
    pathways[pathways$gs_subcat == "TFT_TFT_LEGACY", "gs_subcat"] <- "TFT_TFT_Legacy"
    pathways[pathways$gs_subcat == "CP_KEGG_LEGACY", "gs_subcat"] <- "CP_KEGG"
  }

  pathways$gs_subcat <- gsub(":", "_", pathways$gs_subcat)
  pathways[pathways$gs_cat == "H", "gs_subcat"] <- "HALLMARK"
  celltype_subcat <- .msigdb_celltype_gs_subcat()
  pathways[pathways$gs_cat == celltype_subcat, "gs_subcat"] <- celltype_subcat
  if ("gs_collection" %in% colnames(pathways)) {
    pathways[pathways$gs_collection == celltype_subcat, "gs_subcat"] <- celltype_subcat
  }
  pathways <- as.data.frame(pathways[pathways$gs_subcat %in% pwaycats, , drop = FALSE])
  gs_sizes <- table(pathways$gs_name)
  is_celltype <- pathways$gs_subcat == celltype_subcat
  keep <- rep(TRUE, nrow(pathways))
  keep[!is_celltype] <- gs_sizes[pathways$gs_name[!is_celltype]] <= 500L
  keep[is_celltype] <- gs_sizes[pathways$gs_name[is_celltype]] <= 1000L
  pathways <- pathways[keep, , drop = FALSE]
  pathways <- pathways[gs_sizes[pathways$gs_name] >= 3L, , drop = FALSE]
  pathways
}

#' @keywords internal
.validate_prepared_pathways <- function(pathways, pwaycats, species) {
  req <- c("gene_symbol", "gs_name", "gs_subcat")
  miss <- setdiff(req, colnames(pathways))
  if (length(miss)) {
    stop(
      "Prepared msigdbr pathways missing columns: ", paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  cached_species <- attr(pathways, "species")
  if (!is.null(cached_species) && !identical(cached_species, species)) {
    stop("Cached msigdbr pathways species mismatch.", call. = FALSE)
  }
  present <- unique(pathways$gs_subcat)
  empty <- setdiff(pwaycats, present)
  if (length(empty)) {
    warning(
      "Prepared msigdbr pathways have no gene sets for: ",
      paste(empty, collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(pathways) == 0L) {
    stop("Prepared msigdbr pathways table is empty.", call. = FALSE)
  }
  invisible(pathways)
}

#' @keywords internal
.load_msigdbr_cache <- function(species, pwaycats, cache_dir, fallback_dir) {
  cache_file <- msigdbr_cache_path(
    species = species,
    cache_dir = cache_dir,
    fallback_dir = fallback_dir
  )
  if (is.na(cache_file) || !file.exists(cache_file)) return(NULL)

  pathways <- tryCatch(readRDS(cache_file), error = function(e) NULL)
  if (is.null(pathways)) return(NULL)

  cached_pwaycats <- attr(pathways, "pwaycats")
  if (!is.null(cached_pwaycats)) {
    if (!identical(sort(unname(cached_pwaycats)), sort(unname(pwaycats)))) {
      return(NULL)
    }
  }

  ok <- tryCatch({
    .validate_prepared_pathways(pathways, pwaycats = pwaycats, species = species)
    TRUE
  }, error = function(e) FALSE)
  if (!ok) return(NULL)

  cached_ver <- attr(pathways, "msigdbr_version")
  cur_ver <- as.character(utils::packageVersion("msigdbr"))
  if (!is.null(cached_ver) && !identical(cached_ver, cur_ver)) {
    message(
      "Using cached msigdbr pathways from msigdbr ", cached_ver,
      " (installed: ", cur_ver, "). Delete cache or set refresh_msigdbr_cache = TRUE to rebuild."
    )
  } else {
    message("Reading cached msigdbr pathways: ", cache_file)
  }
  pathways
}

#' @keywords internal
.save_msigdbr_cache <- function(pathways, species, pwaycats, cache_dir, fallback_dir) {
  db_version <- unique(pathways$db_version)
  if (length(db_version) != 1L || is.na(db_version) || !nzchar(db_version)) {
    db_version <- "unknown"
  }
  cache_file <- msigdbr_cache_path(
    species = species,
    cache_dir = cache_dir,
    fallback_dir = fallback_dir,
    db_version = db_version
  )
  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  attr(pathways, "species") <- species
  attr(pathways, "msigdbr_version") <- as.character(utils::packageVersion("msigdbr"))
  attr(pathways, "db_version") <- db_version
  attr(pathways, "scdapp_version") <- as.character(utils::packageVersion("scDAPP"))
  attr(pathways, "prepared_at") <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  attr(pathways, "pwaycats") <- pwaycats
  saveRDS(pathways, cache_file)
  message("Saved msigdbr pathways cache: ", cache_file)
  invisible(cache_file)
}
