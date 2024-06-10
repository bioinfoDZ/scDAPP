# Version changelog

For all changes, please update changelog and use Year-Month-Day


## 1.2.0

2024.06.XX
- add support for DESeq2, DESeq2-LRT and non-LRT EdgeR test in pseudobulk DE comparisons
- add support for all test.use options from Seurat::FindMarkers as possible DE tests in non-pseudobulk DE comparisons
- remove saving of unnecessary .rds file with list of DE tables in DE module, because it is already saved later on in pathway analysis module



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
