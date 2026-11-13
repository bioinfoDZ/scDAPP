# scDAPP downstream: subclustering with `subset_module()`

`subset_module()` takes an integrated Seurat object from a finished scDAPP run, keeps the clusters or cells you pass, reclusters that subset, and runs composition, differential expression, pathway GSEA, and a short HTML report.

Comparisons use the same `comps` table as the original run. See the [comparative designs cookbook](../Comparative_Designs.md).

## Load the object

```r
library(scDAPP)
scDAPP::attach_scDAPP_pipeline_libraries(load_risc = TRUE)
# scDAPP::r_package_test() does the same load and prints a version table

sobjint <- scDAPP::read_scdapp_object(
  "PATH/TO/PIPELINE/OUTPUT/multisample_integration/data_objects/Seurat-object_integrated.qs2"
)
sample_metadata <- read.csv("PATH/TO/sample_metadata.csv")
comps <- read.csv("PATH/TO/comps.csv")
```

`read_scdapp_object()` reads `.qs2` or `.rds`. The default pipeline save is `.qs2`. Pass `load_risc = FALSE` when the parent integration was not RISC.

The object needs `Code` and `Condition`. If your columns have other names, pass `sample_column` and `condition_column`. `cluster_column` defaults to `seurat_clusters`. Use the integrated column from the parent run when that is the grouping you want to subset, for example `RISC_Louvain_npc40_res0.1`.

## Run it

```r
sub <- scDAPP::subset_module(
  sobjint = sobjint,
  sample_metadata = sample_metadata,
  comps = comps,
  outdir = "subset_cluster1",
  clusters = "1",
  cluster_column = "seurat_clusters",
  Pseudobulk_mode = TRUE,
  pcs_int = "auto",
  res_int = "auto",
  species = "Homo sapiens"
)
sub$sobj
sub$cluster_column
```

`clusters` is matched against `cluster_column`. Pass `cells` instead when you already have a vector of cell names. Do not pass both.

`Pseudobulk_mode = TRUE` uses propeller, pseudobulk DE (EdgeR-LRT unless you set `DE_test`), and Milo. `FALSE` uses a chi-square composition test and Wilcoxon, and skips Milo. GSEA runs by default. ORA stays off unless `run_ora = TRUE`.

## How the subset is reclustered

The parent method is read from `sobjint@misc$scDAPP$integration_method`. You can pass `integration_method` if that slot is missing.

**RISC.** The default keeps the parent `Integrated_RISC` assay. Variable genes, scaling, and PCA are computed on that assay, then neighbors and Louvain use the new PCA. This does not call RISC again. The new metadata column is `RISC_subset_npc{pcs}_res{res}`.

`pcs_int` or `res_int = "auto"` selects those values with the same bootstrap used by the pipeline: subsample cells, refit PCA, cluster, and score ARI and Jaccard. Each replicate refits that PCA. It does not re-integrate. `stability_numreps` is the number of replicates (default 50). The sweep is saved under `outdir/cluster_stability/`.

Pass `reintegrate = TRUE` to split the subset back into samples and run RISC again instead of re-embedding.

**Harmony, CCA, RPCA, and both SCT methods.** The subset is split by sample and passed through `run_integration()`. `pcs_int` or `res_int = "auto"` uses that method's integration stability sweep. `reintegrate = FALSE` is only allowed for RISC.

## What is written

Files go in `outdir`, using the same folder names as `multisample_integration/` in a pipeline run:

- `subset_module_report.html`
- `data_objects/subset_object.qs2` (`saving_algorithm = "RDS"` writes `.rds`)
- `cluster_stability/` when PC or resolution is `"auto"`
- `compositional_proportion_analysis/`
- `differentialexpression_crosscondition/`
- `pathwayanalysis_crosscondition/`
- `differential_abundance_Milo/` when Milo runs
- `overrepresentation_pathway_analysis/` when `run_ora = TRUE`

The return value also holds the subset object, the new cluster column, and the module results (`sub$de_results`, `sub$composition`, `sub$pathways`, `sub$milo`).

## Cell types without reclustering

To compare called cell types on the clusters you already have, add a metadata column and pass it as `grouping_variable` to the comparative modules. You do not need `subset_module()`.

```r
sobjint$Celltype <- plyr::mapvalues(
  sobjint$seurat_clusters,
  from = levels(sobjint$seurat_clusters),
  to = c("NK", "B", "Myeloid", "T")
)

de_results <- scDAPP::de_across_conditions_module(
  sobjint = sobjint,
  grouping_variable = "Celltype",
  sample_metadata = sample_metadata,
  comps = comps,
  Pseudobulk_mode = TRUE,
  assay = "RNA",
  slot = "counts",
  outdir_int = "celltype_comparisons"
)
```

`compositional_analysis_module()`, `pathwayanalysis_crosscondition_module()`, and `milo_neighbourhood_da_module()` take the same `grouping_variable`. Details are in `?de_across_conditions_module` and the [cookbook](../Comparative_Designs.md).

## Earlier guide

The manual recluster-then-call-each-module walkthrough from before v2 is [here](Subclustering_and_ComparativeModule_Vignette_prior_to_v2.md).
