# Version changelog

For all changes, please update changelog and use Year-Month-Day


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
