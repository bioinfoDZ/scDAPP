% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sc_processing.R
\name{doubletfinderwrapper}
\alias{doubletfinderwrapper}
\title{DoubletFinder Wrapper}
\usage{
doubletfinderwrapper(seuratobject, clusters, autofilterres, num.cores, sct)
}
\arguments{
\item{seuratobject}{A Seurat object, pre-proc with \code{Seurat::SCTransform()}}

\item{clusters}{string, should match a column name of seuratobject metadata with clusters or other cell group annotations. default = 'seurat_clusters'}

\item{autofilterres}{optional, output of \code{scDAPP::autofilter()} whether to return the the autofilter result with updated \code{autofilterres$cellstatus} taking into account doublet status}

\item{num.cores}{integer, num.cores to use for \code{DoubletFinder::paramSweep}}
}
\value{
if \code{autofilterres} is provided, it will return \code{autofilterres} with updated \code{autofilterres$cellstatus}, if not it will return a data.frame with doublet information and score.
}
\description{
Please see https://github.com/chris-mcginnis-ucsf/DoubletFinder
}
\details{
Models homotypic doublets using the following table:
https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
}
\examples{
# With autofiler res, will add in the result to af$cellstatus
af <- doubletfinderwrapper(sobj, autofilterres = af, num.cores = 5)

# remove doublets
goodcells <- af$cellstatus[af$cellstatus$filteredout==F,"barcodes"]
sobj <- sobj[,goodcells]

# w/o autofilter res:
ddf <- doubletfinderwrapper(sobj, num.cores = 5)

# remove doublets
singlets <- dfdf[dfdf$DoubletFinderClassification=='Singlet','cells']
sobj <- sobj[,singlets]
}
