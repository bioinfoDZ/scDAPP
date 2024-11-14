# scDAPP downstream: sub-clustering and comparative analysis modules


Sometimes, comparative analyses such as compositional, differential expression, or pathway analyses need to be re-run. For example, we may subset the object to remove cells or samples, re-cluster the object, call cell types, or even run these analyses on cell type sub-clusters.

Using the scDAPP pipeline module functions, it is straightforward to run these types of analyses.


## Step 1: Load the scDAPP library and dependency libraries in a specific order

First, we load the scDAPP library and then run the `r_package_test()` function, which itself loads important libraries in the correct order. If the libraries are not loaded or loaded in an incorrect order, errors will occur and the module functions will not work.
```
library(scDAPP)
scDAPP::r_package_test() #this function will load a bunch of relevant packages in the correct order
```

## Step 2: Prepare inputs: read in the Seurat object, sample_metadata.csv, and comps.csv; and define an outdir_int


First, read in the sample_metadata.csv and comps.csv files you used to run the analysis:
```
sample_metadata <- read.csv('PATH/TO/sample_metadata.csv')
comps <- read.csv('PATH/TO/comps.csv')
```

Next, we'll define a directory to save outputs to:

```
outdir_int <- 'PATH/TO/Downstream_Comparisons/'
```

We will read in the saved integrated Seurat object from the pipeline output:

```
sobjint <- readRDS('PATH/TO/PIPELINE/OUTPUT/multisample_integration/data_objects/Seurat-object_integrated.rds')
```



The object has some special metadata including the sample Code, Condition, and clustering:

Minimally, the module functions require Code (designating each sample / replicate) and Condition (designating experimental groups/conditions).
The clustering is an example of a cell-level grouping which will be used as a "grouping_variable". See below.
```
head(sobjint@meta.data, 3)[,c('Code', 'Condition', 'seurat_clusters')]

                                Code        Condition             seurat_clusters
Set1_AAACCTGTCAGGCAAG-1-gPlexJ1 Healthy_2   Healthy               3
Set1_AAAGATGGTTCCCTTG-1-gPlexJ1 Healthy_2   Healthy               2
Set1_AACACGTGTAGCGTGA-1-gPlexJ1 Healthy_2   Healthy               2

```


## Step 3: Prepare a `grouping_variable`

The `grouping_variable` refers to a grouping of cells, such as clusters. For example, across each level of `grouping_variable` (ie, across each cluster), you can compare Condition A vs B. You will likely be interested in modifying this.

Below are three examples:
- Example A: define a new metadata column called "Celltype", and we will pass the string "Celltype" as the `grouping_variable` argument in the comparison module functions.
- Example B: remove some cells, and pass "seurat_clusters" as the `grouping_variable` argument.
- Example C: subset for a single cluster, perform sub-clustering, and pass the new sub-clustering metadata column name as a string to `grouping_variable` in the module functions.


### Example A: defining a "Celltype" metadata column to label clusters

For example, we can combine cluster levels together as celltype calls, and then run the comparative analysis at the celltype level rather than the cluster level. Below, we will create a new column called "Celltype" for this:


```

#check the levels of the cluster variable.
levels(sobjint$seurat_clusters)

# "1" "2" "3" "4" "5" "6" "7"

# Let's say we called a celltype for each of these clusters.
# We'll define a vector where each. Make sure each element below matches the desired cluster above.
# IE, cluster "1" above will be marked as "NK" below.

celltypes <- c('NK', 'B', 'MNP', 'T', 'MNP', 'T', 'MNP')

#create a new metadata column, mapping the cluster levels to the new celltypes.
sobjint$Celltype <- plyr::mapvalues(sobjint$seurat_clusters,
from = levels(sobjint$seurat_clusters),
to = celltypes

)


# plot it to check:
DimPlot(sobjint, group.by = 'seurat_clusters', label = T) + DimPlot(sobjint, group.by = 'Celltype', label = T) 

```

In the above code, we defined a new column called Celltype. We can use this as our `grouping_variable` in the comparative analysis functions.

### Example B: subsetting to remove cells


Here is how you could subset the object, IE, we'll remove cluster 7 in this example.
```
#get metadata
md <- sobjint@meta.data

#subset metadata; in this example, we'll remove cluster 7
md <- md[md$seurat_clusters != 7,]

#the rownames of md are the cell barcodes. We can use this to subset the seurat object:
sobjint_SUBSETTED_noC7 <- sobjint[,rownames(md)]

```

This could be done to remove cells for any reason and re-run the analysis using the same or any different cell grouping as the `grouping_variable`. For example, it can be used for removing weird outlier clusters or samples. Note: be prepared to justify your reasoning whenever removing what you designate as outliers.


### Example C: sub-clustering

You may also be interested in finding sub-clusters for a certain cell grouping or cluster. Here is an example of how you could subset for cluster 1 cells and perform a sub-clustering:

```
#get metadata
md <- sobjint@meta.data

#subset metadata; in this example, we'll select only cluster 1
md <- md[md$seurat_clusters == 1,]

#subset seurat object
sobjint_SUBSETTED_C1only <- sobjint[,rownames(md)]

## SUBCLUSTER - this is how we subcluster with RISC.
# no renormalizing (basically it's just a log norm anyway, plus some special RISC normalization)
# just find HVGs, scale, PCA, graph, and cluster


#set RISC assay as default.
DefaultAssay(sobjint_SUBSETTED_C1only) <- 'RISC'


# find new HVGs of subsetted object
sobjint_SUBSETTED_C1only <- Seurat::FindVariableFeatures(sobjint_SUBSETTED_C1only)

#scale these HVGs and run PCA
sobjint_SUBSETTED_C1only <- Seurat::ScaleData(sobjint_SUBSETTED_C1only)
sobjint_SUBSETTED_C1only <- Seurat::RunPCA(object = sobjint_SUBSETTED_C1only, verbose = F)

#check elbow plot
ElbowPlot(sobjint_SUBSETTED_C1only)

#based on elbow plot, choose PCS to use for graph and UMAP
sobjint_SUBSETTED_C1only <- Seurat::FindNeighbors(object = sobjint_SUBSETTED_C1only, dims = 1:10, verbose = F)
sobjint_SUBSETTED_C1only <- Seurat::RunUMAP(sobjint_SUBSETTED_C1only, dims = 1:10)

#perform louvain sub-clustering with desired resolution
sobjint_SUBSETTED_C1only <- Seurat::FindClusters(object = sobjint_SUBSETTED_C1only, resolution = 0.1, verbose = F)

#inspect
DimPlot(sobjint_SUBSETTED_C1only, label=T)
DimPlot(sobjint_SUBSETTED_C1only, split.by = 'Condition', group.by = 'Code')

#note that this overwrites the "seurat_clusters" variable.


# a note about the clustering nomenclatures you may find after scDAPP: 

#pipeline post-intergration clustering is named like: 
# RISC_Louvain_npc30_res0.5

# pipeline individual sample analysis clustering (before integration) is named like:
# SCT_snn_res.0.5

# the sub-clustering ran here will be named like:
# RISC_snn_res.0.1
# additionally, anytime FindClusters is run, Seurat will create or overwrite the 'seurat_clusters" variable using the latest clustering result

```

After sub-clustering, you can do other types of analysis, such as marker analysis, using Seurat functions:
```
#check assay; set to RISC if not RISC
DefaultAssay(sobjint_SUBSETTED_C1only)
DefaultAssay(sobjint_SUBSETTED_C1only) <- 'RISC'

#find markers
m <- FindAllMarkers(sobjint_SUBSETTED_C1only, only.pos = T)
m$score <- (m$pct.1 - m$pct.2) * m$avg_log2FC #this is how we rank marker genes in scDAPP

#prep genes for heatmap
n <- 5
top <- m %>% group_by(cluster) %>% top_n(n = n, wt = score)

#make sure markers are in the scale.data slot
sobjint_SUBSETTED_C1only <- ScaleData(sobjint_SUBSETTED_C1only, features = top$gene)

#plot heatmap
DoHeatmap(sobjint_SUBSETTED_C1only, features = top$gene)

#if there are no clear markers, you can try to change your number of PCs (dims) or your louvain resolution; or, there may not be real clear sub-clusters in this context

```

After sub-clustering, you can also pass the sub-clustering metadata column name as the `grouping_variable` in the modules as shown below. To be explicit we will pass "RISC_snn_res.0.1" in this example.

<br />
<br />

## Step 4: Compositional analysis module

Using any of the examples above, we can first perform compositional analysis. First, you should check the documentation of the module to see the different options:

```
?scDAPP::compositional_analysis_module()
```

To apply it, you can use something like below. The most important things you need to decide are the `grouping_variable`, and the `compositional_test`. For the latter, we use the [Propeller](https://academic.oup.com/bioinformatics/article/38/20/4720/6675456) test for pseudobulk, and the chi-square test for non-pseudobulk comparison.



Below is how we can run it for Example A, where we defined the Celltype variable:


```
comp_result <- scDAPP::compositional_analysis_module(sobjint,
                                                     grouping_variable = 'Celltype',
                                                     comps=comps, sample_metadata=sample_metadata,
                                                     compositional_test = 'propeller', #'propeller' if pseudobulk; 'chisq' if non-pseudobulk
                                                     outdir_int = outdir_int
)
```


Or we can run it for Example C, where we subsetted for one cluster and then sub-clustered:

```
comp_result <- scDAPP::compositional_analysis_module(sobjint_SUBSETTED_C1only,
                                                     grouping_variable = 'RISC_snn_res.0.1',
                                                     comps=comps, sample_metadata=sample_metadata,
                                                     compositional_test = 'propeller', #'propeller' if pseudobulk; 'chisq' if non-pseudobulk
                                                     outdir_int = outdir_int
)

```

## Step 5: Differential expression and pathway analysis

Next, we will use Example A (Celltype re-mapping) for differential expression (DE) analysis and pathway analysis.


The documentation for these modules is extensive and explains how to use the modules in the correct order. The order is also below:
```
# check DE module
?scDAPP::de_across_conditions_module()

# check pathway preparation module
?scDAPP::preppathways_pathwayanalysis_crosscondition_module()

# check pathway analysis module (uses GSEA)
?scDAPP::pathwayanalysis_crosscondition_module()

# check ORA pathway analysis module
# ?scDAPP::ORA_crosscondition_module()
```



### 5.1 - DE analysis


DE analysis be run with something like below. The most important things to decide are the `grouping_variable` and `Pseudobulk_mode`.

Below is how we can run DE analysis for Example A (Celltype re-mapping):

```
m_bycluster_crosscondition_de_comps <- scDAPP::de_across_conditions_module(sobjint,
                                                                           grouping_variable = 'Celltype', ## grouping variable is the most important thing to change.
                                                                           sample_metadata=sample_metadata, comps=comps,
                                                                           Pseudobulk_mode = T, #set to F if no replicates
                                                                           assay = 'RNA', slot = 'counts', #we use raw counts for pseudobulk analysis; set to normalized RISC values (assay = 'RISC' and slot = 'data') for non-pseudobulk analysis
                                                                           cluster_prefix = F, #if set to T, will append the prefix "cluster_" to each `grouping_variable` level
                                                                           outdir_int = outdir_int
)
```

For Example C, we can change "sobjint" to "sobjint_SUBSETTED_C1only"; and we can change "`grouping_variable` = 'Celltype'" to "`grouping_variable` = 'RISC_snn_res.0.1'", or whatever suits your data.


### 5.2 - Prep pathways

The scDAPP pipeline downloads a pathways database (Molecular Signatures Database, MSIGDB) using the "msigdbr" package. You can re-use the downloaded copy from your prior scDAPP pipeline run (preferred), or just re-download from scratch.

To re-use the downloaded one: read it in from the pipeline output, as shown below. However, note that you may need to do a quick and easy pre-processing if the pipeline was run on older versions. This involves sub-selecting some pathways and changing the pathway table format slightly. Prior to scDAPP v1.2.0 we used to save the whole raw pathway table (but in new versions, the subsetted/reformatted table is saved).

```
#read in pathway table
pathways <- readRDS("PATH/TO/PIPELINE/OUTPUT/multisample_integration/pathwayanalysis_crosscondition/msigdb_pathways.rds")
pathways <- as.data.frame(pathways)


#test if need to pre-process:
'HALLMARK' %in% pathways$gs_subcat

#if TRUE, you can proceed to the next step.

#if FALSE, you need to run the lines below:


#if so, just copy and run the following lines:
#it needs to be subsetted and formatted a bit
pwaycats <- c("HALLMARK", "GO_BP", "GO_MF", "GO_CC", "CP_REACTOME",
              "CP_KEGG", "TFT_GTRD", "TFT_TFT_Legacy")
pwaycats <- gsub(":", "_", pwaycats)
names(pwaycats) <- pwaycats
pathways$gs_subcat <- gsub(":", "_", pathways$gs_subcat)
pathways[pathways$gs_cat == "H", "gs_subcat"] <- "HALLMARK"
pathways <- as.data.frame(pathways[pathways$gs_subcat %in%
                                     pwaycats, ])
pathways <- pathways[table(pathways$gs_name) <= 500, ]
pathways <- pathways[table(pathways$gs_name) >= 3, ]
invisible(gc(full = T, reset = F, verbose = F))

```



Or, you can just re-download the pathways (not preferred, since the pathways database changes over time):

```
#Set species
# see `msigdbr::msigdbr_species()` for a list of usable species
species = 'Homo sapiens'

pathways <- scDAPP::preppathways_pathwayanalysis_crosscondition_module(species = species, outdir_int = outdir_int)
```





### 5.3 - Run pathway analysis (GSEA)

Above, we ran DE analysis and prepared the pathways from a database. We will use both of these to run pathway analysis using GSEA.

```
pathway_analysis_mainlist_comps <- pathwayanalysis_crosscondition_module(
  m_bycluster_crosscondition_de_comps = m_bycluster_crosscondition_de_comps, #output of DE analysis
  pathways = pathways, #the prepped pathways
  sample_metadata = sample_metadata,
  comps = comps,
  workernum = 1, # number of CPUs you want to use.
  outdir_int = outdir_int
)

```


It will save a bunch of tables and plots for each comparison, across each level of the grouping_variable, and for each pathway category. 


At this point, you are pretty much done, and can start just looking at the outputs!


Note that the pathway module function will also save a .rds object of a list containing the DEG object and the pathway results. This can be handy for quickly reading in the data, making new plots, etc.

```
#read in list of DE and pathway results
DE_pathways_plot_objects_list <- readRDS('PATH/TO/SAVED/OUTPUTS/pathwayanalysis_crosscondition/DE_pathways_plot_objects_list.rds')

#this is a complex nested list object:


# DE tables are here:
#check the names here, it will show comparisons
names( DE_pathways_plot_objects_list$m_bycluster_crosscondition_de_comps ) 

#nesting levels are:
#1. comparisons; 2. grouping_variable levels.



# DE pathways detailed analysis is here:
#check the names here, it will show comparisons
names( DE_pathways_plot_objects_list$pathway_analysis_mainlist_comps  )

#nesting levels are:
# 1. comparisons; 2. pathway database categories; 3. grouping_variable levels; 4. pathway table, and pathway plot (two elements)



# DE summary plots are here (these plots summarize GSEA results for a condition and a pathway category, ie, shows pathways from GO_BP upregulated in condition A vs B across all clusters)

#check the names here, it will show comparisons
names(DE_pathways_plot_objects_list$pathwaysummplots_comps)

#nesting levels are:
#1. comparisons; 2. pathway database categories 3. comparison "direction" (Condition A, or Condition B)

```




If you are interested, you can use the new pathway results for aPEAR analysis. Read in the pathways object like below:

```
#read in all DE and pathway results:
DE_pathways_plot_objects_list <- readRDS('PATH/TO/SAVED/OUTPUTS/pathwayanalysis_crosscondition/DE_pathways_plot_objects_list.rds')

#get pathways result list:
pathway_analysis_mainlist_comps <- DE_pathways_plot_objects_list$pathway_analysis_mainlist_comps


```

This object will be the input for aPEAR. We have prepared a [detailed vignette on applying aPEAR](https://github.com/bioinfoDZ/scDAPP/blob/main/Documentation/downstream_postpipeline/aPEAR_wrapper.md) with this.


### 5.4 - Run pathway analysis (ORA)


At a reviewer's request, we added functionality for over-representation analysis (ORA) to scDAPP. This is a fisher-test based analysis, where we actually make a cutoff of DEGs (i.e. select a threshold for pvalue / padj and LFC). It is different from GSEA, the default pathway analysis method shown above, which does not require sometimes arbitrary cutoffs.

The ORA module is similar to the GSEA module.

You must define some cutoffs:
- `crossconditionDE_padj_thres`: adjusted pvalue cutoff. In the pipeline, we set this to 0.1 for pseudobulk, and 0.05 for non-pseudobulk.
- `crossconditionDE_lfc_thres`: log fold change cutoff, provided as an absolute value. In the pipeline, we use no cutoff (set to 0) for pseudobulk, and set this to 0.25 (it will select +/-0.25 l2fc) for non-pseudobulk.
- `crossconditionDE_min.pct`: numeric; Minimum percent of cells expressing gene required to count as a DEG. For positive LFC genes (up in condition A); pct.1 must be at least this value (percent of cells in A must be at least this value); for negative LFC genes, pct.2 cells must be at least this value. Only used for counting DEGs. In the pipeline, we set this to 0.1 if pseudobulk is used; and 0 if wilcox is used.

To run the module, you can use something like below:

```
ORA_analysis_mainlist_comps <- scDAPP::ORA_crosscondition_module(
m_bycluster_crosscondition_de_comps = m_bycluster_crosscondition_de_comps,
pathways = pathways,
sample_metadata = sample_metadata,
comps = comps,
workernum = 1, #num CPUs
outdir_int = outdir_int,
crossconditionDE_padj_thres = 0.1, #we suggest 0.1 for pseudobulk, 0.05 for non-pseudobulk
crossconditionDE_lfc_thres = 0, #we suggest 0 for pseudobulk, 0.25 for non-pseudobulk
crossconditionDE_min.pct = 0.1 #we suggest 0.1 for pseudobulk, 0 for non-pseudobulk
)
```



This vignette was prepared on Sep 22, 2024 by Alexander Ferrena using scDAPP v1.2.1.

