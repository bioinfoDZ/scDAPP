# Version changelog

For all changes, please update changelog and use Year-Month-Day

## 2.0.0
2026.06.09

Major release consolidating integration backends, cross-condition API refactor, MSigDB caching, and cluster-stability bootstrapping.

### Breaking changes
- **`comps.csv`:** uses `c0` (reference) and `c1` (test); legacy `c2` is auto-renamed to `c0`. Optional columns: `formula` (default `~ Condition`), `contrast`, `label`.
- **Module return types:** `de_across_conditions_module()` returns flat `de_results`; `compositional_analysis_module()` returns flat `composition_results`; `pathwayanalysis_crosscondition_module()` returns flat `pathway_results`; `ORA_crosscondition_module()` returns flat `ora_results`. Use `de_results_by_cluster()`, `composition_results_by_comparison()`, `pathway_results_by_cluster()`, and `ora_results_by_cluster()` for per-cluster or per-comparison subsets.
- **RISC assay rename:** batch-corrected expression is stored in Seurat assay `Integrated_RISC` (formerly `RISC`). Update downstream code that referenced assay `"RISC"`.
- **MSigDB cache location:** prepared pathways are cached under user/XDG/R cache (or `outdir/msigdbr_cache` fallback), not in `multisample_integration/pathwayanalysis_crosscondition/msigdb_pathways.rds`. Use `preppathways_pathwayanalysis_crosscondition_module()` or `msigdbr_cache_path()` for downstream work.
- **`integration_method`:** choices are `RISC` (default), `CCAIntegration`, `RPCAIntegration`, `CCAIntegration_SCT`, `RPCAIntegration_SCT`, and `HarmonyIntegration`.
- **Pre-filter clustering:** SCT+Louvain on the full unfiltered matrix is opt-in via `cluster_unfiltered` (default `FALSE`). `Unfiltered-SeuratObject-*.rds` outputs are written only when `cluster_unfiltered = TRUE`.

### Integration
- New `integration_method` pipeline parameter and unified `run_integration()` entry point with Seurat v5 `IntegrateLayers` backends (`prepare_integration_inputs()`, `integration_method_choices()`, `resolve_integration_config()`, `check_integration_dependencies()`).
- Integrated clustering, plotting, and DE assays are chosen per method via shared config. `HarmonyIntegration` uses CRAN `harmony` via Seurat's `HarmonyIntegration`, not SeuratWrappers.
- Seurat 5.5 compatibility: assay-aware `merge()`/`split()`; capped integration `dims`; `FindVariableFeatures` before SCT `IntegrateLayers`; `IntegrateLayers` passes `k.filter = 200` for CCA/RPCA.
- `CCAIntegration_SCT` / `RPCAIntegration_SCT` follow Seurat v5 workflow (merge, split RNA by sample, `SCTransform`, `IntegrateLayers(..., normalization.method = "SCT")`). Integrated cluster markers call `PrepSCTFindMarkers()` for SCT-based integration.
- Integration report (Rmd): method-conditional prose; RISC reference-selection section only for RISC; Seurat backends get a dedicated integration summary note.
- RISC count matrices (`input_seurat_obj = TRUE`): per-sample data read once at pipeline start; `prepare_integration_inputs()` builds and saves `.concatmatrix.rds` for integration and stability sweeps.

### Cluster-stability auto-tuning
- `pcs_int` and `res_int` accept `"auto"` for bootstrap cluster-stability selection (ARI + Jaccard) before the final integration run.
- New pipeline parameters: `stability_outdir`, `stability_numreps` (default 50), `stability_sweep_maxPCs`, `stability_sweep_res`, `stability_propcells.perrep` (default 0.8), `stability_workernum`.
- New exports: `cluster_stability_sweep()`, `auto_integration_cluster_params()`, `resolve_integration_cluster_params()`, `integration_stability_plots_module()`.
- Stability CSVs, checkpoint resume, and PDF plots under `{stability_outdir}/plots/`; figures rendered in the HTML report when auto tuning runs. Requires CRAN package `mclust` (added to Suggests).

### Cross-condition analysis
- Centralized DEG thresholds and counting: `crosscondition_de_threshold_defaults()`, `count_crosscondition_degs()`, `filter_significant_de_results()`. Mode-specific defaults when pipeline thresholds are omitted (pseudobulk padj 0.1 / LFC 0 / min.pct 0.1; Wilcox padj 0.05 / LFC 0.25 / min.pct 0).
- New pipeline parameter `crossconditionDE_min.pct` (runner + Rmd); HTML report reuses precomputed DEG count tables from `de_prepplots` instead of recomputing.
- `count_crosscondition_degs()` validated against numDEGs summary CSV logic (same filtering rules as `*_numDEGs_summary.csv`).
- Pseudobulk DE fits joint models per `comps$formula` group (e.g. `~ Condition + Batch`). Compositional analysis (`propeller` mode) respects `comps$formula` via `propeller.ttest` when covariates are present.
- New DE report plots: `plot_crosscondition_deg_barplot()` (signed Up/Down counts per cluster) and `plot_crosscondition_deg_dotplot()` (comparison summary and per-cluster Seurat DotPlots of top DEGs).
- New exports: `normalize_comps()`, `select_top_crosscondition_deg_genes()`, `select_top_crosscondition_deg_genes_pooled()`.

### MSigDB and ORA
- Central MSigDB cache keyed by species, msigdbr version, and MSigDB `db_version`; validated for msigdbr v10+ column names. New pipeline parameter `msigdbr_cache_dir`; module parameter `refresh_msigdbr_cache` on `preppathways_pathwayanalysis_crosscondition_module()`.
- New exports: `msigdbr_cache_path()`, `resolve_msigdbr_cache_dir()`.
- Cell-type signature ORA: `ORA_cluster_markers_msigdb_celltype_module()` and `ORA_cluster_markers_msigdb_celltype_batch_module()` on individual and integrated marker tables, with top-5-per-cluster summary dotplots in the HTML report. Outputs under `celltype_marker_prediction/` (CSV tables and dotplot PDF). MSigDB cache retains C8 cell-type gene sets.
- New pipeline parameters: `run_msigdb_celltype_ora` (default TRUE), `msigdb_celltype_ora_marker_padj_thres`, `msigdb_celltype_ora_top_markers`.

### QC, reporting, and infrastructure
- Per-sample QC extracted to module functions: `per_sample_qc_module()` and `add_per_sample_qc_metadata()`.
- `attach_scDAPP_pipeline_libraries()` centralizes pipeline `library()` order; `r_package_test()` uses it.
- `set_parallel_blas_threads()` limits BLAS threads during pipeline runs (uses `RhpcBLASctl` when installed; added to Suggests).
- Pipeline report and QC PDFs show `twt_colored_heatmap()` alongside alluvial plots. `twt_colored_heatmap()` gains `color_by` (`row_prop`, `count`, `row_scaled`; default `row_prop`) and optional `title`.
- `apear_data_prep()` accepts flat `pathway_results` (legacy nested RDS still supported).
- Docs: `Documentation/Usage.md` and integration runtime expectations.

## 1.3.3
2025.05.18
- apptainer-related update: copy the .rmd file from the package to outdir so tmp file creation and rendering happens in a executable dir (since apptainer is strict with read-only once made)

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
