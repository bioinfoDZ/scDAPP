# Read user-supplied CSVs (sample_metadata, comps), tolerating Excel export quirks.

#' Read a user CSV, allowing missing final newline and UTF-8 BOM.
#'
#' Excel "Save As CSV" often omits the last newline (`incomplete final line
#' found by readTableHeader`) and may write a UTF-8 BOM. This helper parses
#' anyway, warns about the missing newline, strips a BOM, and drops trailing
#' all-empty rows.
#'
#' @param path Path to a CSV file.
#' @param ... Passed to [utils::read.csv()] (except `text`, which is set here).
#' @return A data.frame.
#' @keywords internal
.read_user_csv <- function(path, ...) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("CSV path must be a single non-empty string.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("CSV file not found: ", path, call. = FALSE)
  }
  size <- file.info(path)$size
  if (is.na(size) || size < 1) {
    stop("CSV file is empty: ", path, call. = FALSE)
  }

  raw <- readBin(path, what = "raw", n = size)
  bom <- as.raw(c(0xEF, 0xBB, 0xBF))
  if (length(raw) >= 3L && identical(raw[seq_len(3L)], bom)) {
    raw <- raw[-seq_len(3L)]
  }
  if (!length(raw)) {
    stop("CSV file is empty after removing UTF-8 BOM: ", path, call. = FALSE)
  }

  missing_nl <- !(raw[length(raw)] %in% as.raw(c(0x0A, 0x0D)))
  if (isTRUE(missing_nl)) {
    warning(
      "CSV does not end with a newline (common when exporting from Excel): ",
      path,
      ". Parsed the last row anyway. Re-save as CSV (UTF-8) so the file ends ",
      "with a newline if you want to silence this warning.",
      call. = FALSE
    )
    raw <- c(raw, as.raw(0x0A))
  }

  txt <- rawToChar(raw)
  Encoding(txt) <- "UTF-8"

  dots <- list(...)
  if (!("stringsAsFactors" %in% names(dots))) {
    dots$stringsAsFactors <- FALSE
  }
  dots$text <- txt
  df <- do.call(utils::read.csv, dots)

  if (nrow(df) > 0L) {
    empty <- vapply(
      seq_len(nrow(df)),
      function(i) {
        r <- unlist(df[i, , drop = TRUE], use.names = FALSE)
        all(is.na(r) | !nzchar(trimws(as.character(r))))
      },
      logical(1)
    )
    while (length(empty) && isTRUE(empty[length(empty)])) {
      df <- df[-nrow(df), , drop = FALSE]
      empty <- empty[-length(empty)]
    }
  }
  df
}
