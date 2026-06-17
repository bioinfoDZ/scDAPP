#' @keywords internal
.pipeline_library_pkgs <- function(load_risc = FALSE) {
  pkgs <- c(
    "tidyverse",
    "patchwork",
    if (isTRUE(load_risc)) "RISC",
    "Seurat",
    "scDAPP",
    "DoubletFinder",
    "future",
    "parallel",
    "foreach",
    "glmGamPoi",
    "ComplexHeatmap",
    "ggdendro",
    "ggridges",
    "edgeR",
    "msigdbr",
    "hdf5r",
    "ggalluvial",
    "ggfittext",
    "ggrepel"
  )
  unique(pkgs)
}

#' @keywords internal
.pipeline_version_check_pkgs <- function() {
  c(
    "tidyverse",
    "Seurat",
    "patchwork",
    "ggdendro",
    "foreach",
    "msigdbr",
    "ggalluvial",
    "ggfittext",
    "ggrepel",
    "hdf5r",
    "edgeR",
    "glmGamPoi",
    "fgsea",
    "ComplexHeatmap",
    "DoubletFinder",
    "RISC",
    "scDAPP",
    "ggridges"
  )
}

#' @keywords internal
.lib_pkg <- function(pkg, quietly = FALSE) {
  if (quietly) {
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  } else {
    library(pkg, character.only = TRUE)
  }
  invisible(pkg)
}

#' @keywords internal
.check_pipeline_namespace_masks <- function(load_risc = FALSE) {
  if (!isTRUE(load_risc)) {
    return(invisible(NULL))
  }
  if (!requireNamespace("RISC", quietly = TRUE)) {
    return(invisible(NULL))
  }
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    return(invisible(NULL))
  }
  if (!"DimPlot" %in% ls(envir = .GlobalEnv, all.names = TRUE)) {
    return(invisible(NULL))
  }
  seurat_dimplot <- getFromNamespace("DimPlot", "Seurat")
  resolved <- get("DimPlot", envir = .GlobalEnv)
  if (!identical(resolved, seurat_dimplot)) {
    warning(
      "DimPlot in the global environment is not Seurat::DimPlot after attach. ",
      "Use Seurat::DimPlot() in package code and check library() order.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Attach scDAPP pipeline dependencies in a safe order
#'
#' Loads packages in the order used by
#' [inst/rmd/scRNAseq_clustering_integration.Rmd](inst/rmd/scRNAseq_clustering_integration.Rmd):
#' tidyverse and patchwork first, optional RISC, then Seurat (so Seurat wins
#' namespace conflicts such as `DimPlot`), then scDAPP and remaining dependencies.
#'
#' Call this before running the pipeline Rmd, smoke scripts, or interactive tests
#' that rely on Seurat/tidyverse. `devtools::load_all()` alone is not sufficient.
#'
#' @param load_risc logical; attach RISC when TRUE (pipeline default when
#'   `integration_method` is RISC).
#' @param set_seed optional seed passed to [set.seed()]; `NULL` skips seeding.
#' @param configure_parallel if TRUE, call [set_parallel_blas_threads()] and set
#'   `options(future.globals.maxSize = future_globals_max_size)`.
#' @param future_globals_max_size value for `future.globals.maxSize` when
#'   `configure_parallel` is TRUE.
#' @param quietly if TRUE, suppress package startup messages.
#'
#' @return Invisibly, a named character vector of package versions for packages
#'   that were requested to attach.
#' @export
#'
#' @seealso [r_package_test()] for install verification with the same attach order.
#'
#' @examples
#' \dontrun{
#' scDAPP::attach_scDAPP_pipeline_libraries(load_risc = FALSE)
#' }
attach_scDAPP_pipeline_libraries <- function(
    load_risc = FALSE,
    set_seed = 2022L,
    configure_parallel = TRUE,
    future_globals_max_size = 15000 * 1024^2,
    quietly = FALSE
) {
  if (isTRUE(load_risc) && !requireNamespace("RISC", quietly = TRUE)) {
    warning(
      'load_risc = TRUE but package "RISC" is not installed.',
      call. = FALSE
    )
  }

  pkgs <- .pipeline_library_pkgs(load_risc = load_risc)

  for (pkg in pkgs) {
    if (pkg == "RISC" && !requireNamespace("RISC", quietly = TRUE)) {
      next
    }
    .lib_pkg(pkg, quietly = quietly)
  }

  if (isTRUE(configure_parallel)) {
    scDAPP::set_parallel_blas_threads()
    options(future.globals.maxSize = future_globals_max_size)
  }

  if (!is.null(set_seed)) {
    set.seed(set_seed)
  }

  .check_pipeline_namespace_masks(load_risc = load_risc)

  vers <- vapply(
    pkgs,
    function(pkg) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        return(NA_character_)
      }
      tryCatch(
        as.character(utils::packageVersion(pkg)),
        error = function(e) NA_character_
      )
    },
    character(1)
  )

  invisible(vers)
}
