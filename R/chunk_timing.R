# Wall-clock seconds for each evaluated knitr chunk, in run order.
.scdapp_chunk_timing <- new.env(parent = emptyenv())

#' Time each evaluated knitr chunk.
#'
#' Registers a chunk hook that records elapsed wall-clock time. The hook adds
#' nothing to chunk output. Chunks with \code{eval = FALSE} are skipped.
#' Times are kept in a package environment so a later \code{rm(list = ls())}
#' in the report does not drop them.
#'
#' @return Invisible \code{NULL}.
#' @keywords internal
install_chunk_timing_hook <- function() {
  rm(list = ls(envir = .scdapp_chunk_timing, all.names = TRUE), envir = .scdapp_chunk_timing)
  .scdapp_chunk_timing$starts <- list()
  .scdapp_chunk_timing$elapsed <- list()
  .scdapp_chunk_timing$order <- character()

  knitr::knit_hooks$set(scdapp_chunk_time = function(before, options) {
    if (isFALSE(options$eval)) {
      return(NULL)
    }
    label <- options$label
    if (isTRUE(before)) {
      .scdapp_chunk_timing$starts[[label]] <- proc.time()[["elapsed"]]
    } else {
      start <- .scdapp_chunk_timing$starts[[label]]
      if (is.null(start)) {
        return(NULL)
      }
      .scdapp_chunk_timing$elapsed[[label]] <- proc.time()[["elapsed"]] - start
      .scdapp_chunk_timing$order <- c(.scdapp_chunk_timing$order, label)
    }
    NULL
  })
  knitr::opts_chunk$set(scdapp_chunk_time = TRUE)
  invisible(NULL)
}

#' Chunk timing table in the order chunks ran.
#'
#' @return A data frame with \code{step}, \code{chunk}, \code{seconds},
#'   \code{minutes}, \code{hours}, and \code{percent} of the summed chunk time.
#'   Zero rows when the hook has not recorded anything.
#' @keywords internal
chunk_timing_table <- function() {
  labels <- .scdapp_chunk_timing$order
  empty <- data.frame(
    step = integer(),
    chunk = character(),
    seconds = numeric(),
    minutes = numeric(),
    hours = numeric(),
    percent = numeric(),
    stringsAsFactors = FALSE
  )
  if (is.null(labels) || !length(labels)) {
    return(empty)
  }
  sec <- vapply(labels, function(lab) {
    val <- .scdapp_chunk_timing$elapsed[[lab]]
    if (is.null(val) || !is.finite(val)) NA_real_ else as.numeric(val)
  }, numeric(1), USE.NAMES = FALSE)
  total <- sum(sec, na.rm = TRUE)
  pct <- if (is.finite(total) && total > 0) 100 * sec / total else rep(NA_real_, length(sec))
  data.frame(
    step = seq_along(labels),
    chunk = labels,
    seconds = round(sec, 1),
    minutes = round(sec / 60, 2),
    hours = round(sec / 3600, 3),
    percent = round(pct, 1),
    stringsAsFactors = FALSE
  )
}
