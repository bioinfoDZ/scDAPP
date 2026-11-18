# Version changelog

For all changes, please update changelog and use Year-Month-Day

## 2.0.0
2026.11.18

#### Comparative analysis v2

- Pseudobulk DE now supports more complex contrasts, ie adjusting for confounders, studying interactions, or handling paired designs: `EdgeR`, `EdgeR-LRT`, `EdgeR-QLF`, `DESeq2`, `DESeq2-LRT`, and `Dream` (`(1|Patient)`); EdgeR-style vs DESeq2-style contrasts; blank contrast defaults to Condition c1 vs c0; early preflight (comps/design and integration packages) and soft-skip for sparse clusters.
- Cross-condition API: `comps` uses `c0`/`c1` (legacy `c2` → `c0`) with optional `formula`, `contrast`, and `label`; DE / compositional / pathway / ORA modules return flat tables.
- Propeller compositional analysis follows the same formulas/contrasts (covariates, interactions, blocked pairing).
- DE now performs joint normalization and gene dispersion estimation of all samples (instead of just per-A vs B comparison)
- New comparative designs cookbook ([Comparative_Designs.md](Comparative_Designs.md)) and Usage updates.
- DE visalization for every design, including plain `~ Condition`: one pseudobulk heatmap per cluster (columns are samples, annotations are the formula variables), and a second heatmap only when that cluster has more than 50 significant genes.
- Pathway report: one summary dot plot with positive and negative NES together, plus a leading-edge mean heatmap per cluster using the same sample annotations. `pathwayanalysis_crosscondition_legacy()` keeps the previous up/down plots.
- Optional Milo neighbourhood differential abundance (`run_MILO`, `milo_k`, `milo_k_max`). Saved under `multisample_integration/differential_abundance_Milo/`. The report shows the top 20 significant neighbourhoods. Paired formulas use Milo's mixed model with `force = TRUE`, and drop `(1|var)` with a warning if that fit fails or the contrast is not the last fixed coefficient.

#### Integration methods and stability selection v2

- Multi-backend integration via `integration_method` (RISC default, plus Seurat CCA/RPCA/Harmony and SCT variants); integrated RISC assay renamed to `Integrated_RISC`.
- Cluster-stability auto-tuning: `pcs_int` / `res_int` can be `"auto"` (bootstrap ARI + Jaccard) before the final integration.
- RISC reference default `"autoV2"` (product of InPlot cluster × Stv / mean KS, after a KS outlier veto); chosen sample is frozen for the stability sweep and final integration (`"auto"` / `"autoV1"` keeps the legacy heuristic).


#### Additional pipeline updates
- Combined per-sample cluster-marker and MSigDB cell-type ORA CSVs (`all_samples_clustermarkers.csv`, `all_samples_ora_results.csv`). Based on cluster marker genes, find associations with published cell type signatures; aids with celltype annotation of numeric clusters.
- `subset_module()` subsets an integrated object, reclusters, then runs composition, DE, GSEA, optional Milo, and writes `subset_module_report.html`. RISC re-embeds `Integrated_RISC` and can bootstrap PCs and resolution. Harmony, CCA, RPCA, and SCT re-integrate.
- `comparative_module()` runs composition, DE, GSEA, and optional Milo on an existing integrated object grouped by a metadata column (for example called cell types), without subsetting or reclustering. Writes `comparative_module_report.html`.

#### Performance and efficiency updates
- Safe pipeline resume: fingerprint-gated caches under `{outdir}/.scdapp_resume/` for per-sample QC/processing, stability auto-tuning, and integration (including saved InPlot); comps/DE always re-run. Force a full recompute with a new outdir or by deleting `.scdapp_resume/`.
- `saving_algorithm` selects `qs2` (default) or `RDS` for pipeline object files. `qs2` saves large objects faster. Sample, reference, and marker inputs may be either `.rds` or `.qs2`.
- The HTML report Run time section lists wall-clock time for each code chunk in the order it ran, then the total hours for the whole report.
- Add `conda_env/2026_scdappv2.yml`: R >= 4.5, Seurat >= 5, and the v2 feature packages that are on conda (`harmony`, `mclust`, `clusterProfiler`, `variancePartition`, `miloR`, `DESeq2`, `presto`). `2025scdapp.yml` remains for existing environments. `RISC` and `DoubletFinder` are still installed from GitHub.
- MSigDB pathway cache; cell-type marker ORA; DEG bar/dot plots; per-sample QC modules and `attach_scDAPP_pipeline_libraries()`.



<br />
<br />


## 1.3.1
2025.05.14
- add a fix for the MSIGDBR v10 to v24 update: remove "msigdbdf" from remotes imports
- add new Conda environment .yml file "2025scdapp.yml" (including Conda-based Seurat installation and many Bioconductor dependencies)
- update docs: account for new conda env, remove references to mamba since conda seems to use it by default now, other small changes

## 1.3.0
2025.05.09
- add a new feature: check if Seurat objects already have "percent.mito", "percent.hemoglobin" or "Phase" in the metadata before computing these, thus allowing user to pass these if using seurat objects as input. Useful for when running organisms besides human/mouse.
- fix a bug to fully utilize assay and slot (layer) in the de cross conditions module; allows use of the function outside of the pipeline, ie even if no RISC assay is present in the object
- aPEAR was removed from CRAN on 2025.01.10; suggest a method to install from CRAN archive in the docs.
- add a function `twt_colored_heatmap()` to easily visualize a two-way-table of categorical variables with many levels, a nice alternative to alluvial plots
- add a fix for the MSIGDBR v10 update: adjust column names and some subcategory ("subcollection") names to match old formats


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
