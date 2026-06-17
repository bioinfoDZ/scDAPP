#' Limit BLAS/OpenMP threads for fork-based parallel jobs
#'
#' Sets thread environment variables (and optionally \code{RhpcBLASctl}) so
#' \code{foreach}/\code{mclapply} workers do not each spawn multi-threaded BLAS.
#' Call once at the start of a pipeline session, before parallel or \code{irlba}
#' work (e.g. \code{RISC::InPlot} with \code{ncore > 1} on macOS).
#'
#' @return Invisible \code{NULL}.
#' @export
set_parallel_blas_threads <- function() {
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1"
  )
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L)
  }
  invisible(NULL)
}
