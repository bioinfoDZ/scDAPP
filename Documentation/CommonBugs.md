# List of known bugs and fixes

<br />
<br />

## 1 "function 'as_cholmod_sparse' not provided by package 'Matrix'"

2023.11.20


This can be fixed by re-installing irlba. In general, installing irlba after Matrix solves this issue.

```
install.packages("irlba", type = "source")
```

If that does not work, then it should be fixed by removing and reinstalling the packages in exactly this order:

```
remove.packages('Matrix')
remove.packages('irlba')
install.packages("Matrix", type = "source")
install.packages("irlba", type = "source")
```

See here: 
https://github.com/bwlewis/irlba/issues/70

<br />
<br />



## 2 Error in SCTransform (useNames = NA is defunct. Instead, specify either useNames = TRUE or useNames = FALSE)

2023.12.18

please see: 
https://github.com/satijalab/seurat/issues/8183
https://github.com/satijalab/seurat/issues/7501#issuecomment-1854571904
```
# restart your session and run previous scripts
remotes::install_version("matrixStats", version="1.1.0") 

```
<br />
<br />



## 3 Error: 'paramSweep_v3' is not an exported object from 'namespace:DoubletFinder'

2023.12.18

DoubletFinder removed the functions with "_v3" suffix. Please force reinstall DF as below

```
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force = T)
```
And make sure to use latest version of scDAPP.

<br />
<br />



## 4 Error: At vendor/cigraph/src/community/louvain.c:601 : Weight vector length must agree with number of edges. Invalid value

2024.03.17: install RISC 1.7.

2024.02.21

If using RISC 1.6.0, you need to downgrade the version of the C igraph package, version 0.10.7 should work. (Note there are three versions of igraph, including C igraph, R igraph package, and python-igraph, but our testing indicates this error is caused by C igraph version higher than 0.10.8 or 0.10.9 and above)

RISC will be updated with a fix soon.

<br />
<br />




## 5 RISC install error with matrix.utils --> no longer a dependency of RISC 1.7 (as of 03/17/2024)

The dependency was removed from CRAN, so we installed it from an archive. Try using one of the following:
```
install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix.utils/Matrix.utils_0.9.8.tar.gz", type = "source", repos = NULL)
```


```
remotes::install_github("cvarrichio/Matrix.utils")
```

<br />
<br />





## 6 configure: error: geos-config not found or not executable.
ERROR: configuration failed for package ‘rgeos’
* removing ‘/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/rgeos’
* restoring previous ‘/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/rgeos’


Fix: In the Mac terminal, use brew to install

```
brew install geos gdal
```


https://stackoverflow.com/questions/50997636/problems-installing-rgeos-and-rgdal-on-mac-os-x-high-sierra



## 7 "‘==’ only defined for equally-sized data frames"

This error may arise early in the pipeline, for example during the "filtering_main" block.

It may be caused by an issue with dependency packages ggplot2 and patchwork. Solution: update to at least ggplot2 >= 3.5.1 and patchwork >= 1.2.0.

Relevant github issues:
- https://github.com/satijalab/seurat/issues/8170
- https://github.com/thomasp85/patchwork/issues/342



