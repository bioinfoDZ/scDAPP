#' 10X version 3 multiplet rate
#'
#' Data.frame with doublet rates for 10X v3 data.
#'
#' @format A data frame with 11 rows and 3 variables:
#' \describe{
#'   \item{MultipletRate}{Rate of multiplets}
#'   \item{CellsLoaded_100.Viability}{Number of cells loaded, assuming 100% viability}
#'   \item{CellsRecovered}{Number of non-doublet cells recovered after doublet removal}
#' }
#' @source \url{https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-}
"dratedf"
