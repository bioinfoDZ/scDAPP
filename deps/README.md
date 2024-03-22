The [RISC github repo](https://github.com/bioinfoDZ/RISC) will eventually be updated to 1.7, but for now, we host the tar.gz file here.


First, explicitly install key dependencies (which may not happen when installing RISC from a tar.gz file as we do here): 
```
install.packages(c("Matrix", "sparseMatrixStats", "Matrix.utils", "irlba", "doParallel", "foreach", "data.table", "Rtsne", "umap", "MASS", "pbapply", "Rcpp", "RcppArmadillo", "densityClust", "FNN", "igraph", "RColorBrewer", "ggplot2", "gridExtra", "pheatmap", "hdf5r"))
```

To download the file RISC 1.7 tar.gz file, click on the tar.gz file above, and then click the downward facing arrow on the top right.

Download the file, then run the line below, make sure to adjust path.
```
install.packages('~/Downloads/RISC_1.7.tar.gz', repos = NULL, type = "source")
```
