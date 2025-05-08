# Version changelog

For all changes, please update changelog and use Year-Month-Day

## 1.3.0
2025.03.31
- add a new feature: check if Seurat objects already have "percent.mito", "percent.hemoglobin" or "Phase" in the metadata before computing these, thus allowing user to pass these if using seurat objects as input. Useful for when running organisms besides human/mouse.
- fix a bug to fully utilize assay and slot (layer) in the de cross conditions module; allows use of the function outside of the pipeline, ie even if no RISC assay is present in the object
- aPEAR removed from CRAN on 2025.01.10; suggest a method to install from CRAN archive.
- add a function `twt_colored_heatmap()` to easily visualize a two-way-table of categorical variables with many levels, a nice alternative to alluvial plots
- add a fix for the MSIGDBR v10 update: adjust column names and soem subcategory ("subcollection") names to match old formats


## 1.2.2
2024.10.29
- In Description file, add BiocManager as an Imports dependency, hopefully to fix auto-install of sparseMatrixStats from fresh
- In Description file, add clusterProfiler to biocViews and Suggests
- Fix a bug in the "de_preplots" block related to adding missing samples as zeros, missing zero df was not being given the colnames of the missing samples
- update docs


## 1.2.1
2024.09.22
- fix a bug related to parallelization memory reservation caused by external dependency update
- add presto package to the description file remotes for auto-install
- place presto to new suggests field of description file and also move aPEAR from imports to suggests
- add an exception to make sure all condition levels in comps.csv are present in sample_metadata$Condition (second column)
- fix an error with DEG heatmap plotting when pseudobulk mode = T
- add a new vignette for subclustering and running comparative modules


## 1.2.0

2024.06.08
- `DE_test` parameter: add support for DESeq2, DESeq2-LRT and non-LRT EdgeR test in pseudobulk DE comparisons
- `run_ORA` parameter: add optional support for ORA test of GO, KEGG and other pathways including those also tested via GSEA
- add support for all test.use options from Seurat::FindMarkers as possible DE tests in non-pseudobulk DE comparisons (though so far only MAST has been tested - note some of the tests from Seurat return a very different output format and will likely not work, such as "roc")
- remove saving of unnecessary .rds file with list of DE tables in DE module, because it is already saved later on in pathway analysis module
- fix minor issue causing confusing warning messages related to theme_dimplot docs during install
- preppathways module now saves the subsetted and and reformatted (Hallmark entered as both category and subcategory) pathways object rather than the raw pathway tibble from msigdbr




## 1.1.0

2024.04.07
- fix a bug in DEG weighting for GSEA caused by incorrect gene matching
- fix a bug in pathway analysis module which sometimes caused errors due to active devices


## 1.0.2

2024.04.30
- fix typos in report
- update default title and author fields


## 1.0.1

2024.04.19
- fix issue with risc_reference when set as auto instead of NULL or left empty


## 1.0.0

2024.02.29 - leap day.
- Finalized pipeline for publication.



For older development versions, see [SDAP](https://github.com/FerrenaAlexander/SDAP) or [devel scDAPP](https://github.com/FerrenaAlexander/scDAPP/).
