% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aPEAR_wrapper_functions.R
\name{apear_update_and_save_input}
\alias{apear_update_and_save_input}
\title{Update and Save aPEAR Input Data}
\usage{
apear_update_and_save_input(
  apear_input_list,
  findPathClusterres,
  output_dir = NULL,
  save_apear_input = TRUE,
  ...
)
}
\arguments{
\item{apear_input_list}{A list of sub-lists, which is the output from the \code{apear_data_prep} function.}

\item{findPathClusterres}{A list of sub-lists, which is the output from the \code{apear_find_clusters} function.}

\item{output_dir}{Directory where the updated \code{apear_input_list} and its corresponding CSV files are saved.\code{NULL} by default.}

\item{save_apear_input}{A boolean that controls whether the updated apear input list is saved. If \code{TRUE},
directories are created for each comparison group (\link{condition A}\emph{vs}\link{condition B}) with a subfolder "inputs" containing CSV files
for each cell cluster dataframe, which now include "number_of_LE_genes", "Select_for_aPEAR", and "aPEAR_Cluster" columns.}

\item{...}{Additional arguments passed to the underlying clustering functions.}
}
\value{
An updated \code{apear_input_list} containing sub-lists named by comparison groups (\link{condition A}\emph{vs}\link{condition B}).
Each sub-list includes dataframes corresponding to unique cell clusters, now with 12 columns: "Description",
"Category", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "core_enrichment", "Des_Full",
"number_of_LE_genes", "Select_for_aPEAR", and "aPEAR_Cluster".
}
\description{
This function integrates the outputs from \code{apear_data_prep} and \code{apear_find_clusters} functions.
It updates the \code{apear_input_list} with clustering information and saves the result if specified.
}
\examples{
# Assuming `apear_input_list` and `findPathClusterres` have been defined:
apear_input_list_updated <- apear_update_and_save_input(apear_input_list,
                                                        findPathClusterres,
                                                        output_dir = "path/to/save",
                                                        save_apear_input = TRUE)

}
