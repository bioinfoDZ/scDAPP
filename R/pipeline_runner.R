# https://r-pkgs.org/whole-game.html

#' Test installation of the packages
#'
#' This will load the key dependencies and report the versions of those dependencies.
#'
#' @return a data.frame with packages and version numbers
#' @export
#'
#' @examples
#' \dontrun{
#' SpatialSingulomics::r_package_test()
#' }
r_package_test <- function(){


  packages <- c(
    #CRAN
    "tidyverse",  # general data wrangling
    "Seurat",     # spatial analysis
    "patchwork",  # combine plots
    "ggdendro",       #for clustering dendrograms
    "foreach",    # parallelization
    "msigdbr",          #get pathways (cross species", from msigdb
    "ggalluvial", # part of alluvial plot
    "ggfittext", # part of alluvial plot
    "ggrepel", # part of alluvial plot
    "hdf5r", # generally a hard oen to install, seurat dep


    #Bioconductor
    "edgeR",     # optional, for edgeR pseudobulk DE
    "glmGamPoi",  # for faster  SCT
    "fgsea",              #GSEA / pathway analysis
    "ComplexHeatmap", # for heatmaps

    #Github
    "DoubletFinder",
    "RISC",
    "scDAPP"

  )

  packages <- data.frame(pkg = packages)

  packages$vers <- sapply(packages$pkg, function(pkg){
    tryCatch({as.character(packageVersion(pkg))},
             error=function(cond) {
               warning('package "',pkg, '" not detected!!')
               # Choose a return value in case of error
               return(NA)
             })
  }, simplify = T)



  library(tidyverse)
  library(patchwork)  # combine plots
  library(RISC)
  library(Seurat)
  library(scDAPP)
  library(DoubletFinder)
  library(future)
  library(parallel)
  library(foreach)
  library(glmGamPoi)  # for faster SCT
  library(ComplexHeatmap) # for heatmaps
  library(ggdendro)       #for clustering dendrograms
  library(ggridges) # qc ridgeplots
  library(edgeR)
  library(msigdbr)          #get pathways (cross species) from msigdb
  library(hdf5r) # HARD TO INSTALL: installed thru mamba
  library(ggalluvial) # part of alluvial plot
  library(ggfittext) # part of alluvial plot
  library(ggrepel) # part of alluvial plot

  return(packages)
}



#' scRNA-seq analysis, integration, and comparative DE pipeline
#'
#' This will run a pipeline of Seurat individual sample analysis, RISC integration, and comparative DE. Multiple conditions (A vs B vs C) are supported. Also, pseudobulk DE or Wilcox are available for comparative DE. Finally, there is an option to use label transfer with a reference single-cell RNAseq dataset.
#'
#' @param datadir string, path to folder containing Cellranger output folders for each sample
#' @param outdir string, path to output folder, will be created if doesn't already exist
#' @param use_labeltransfer T/F, whether to use labeltransfer for cell type prediction, default = F
#' @param refdatapath string, path to a Seurat object .rds file for labe latransfer, pre-processed with `Seurat::SCTransform()`, with a column called "Celltype" in its meta.data. Ignored if `use_labeltransfer` = F.
#' @param m_reference string, path to .rds file containing output of `Seurat::FindAllMarkers` run on the reference object specified above. Ignored if `use_labeltransfer` = F
#' @param sample_metadata string, path to a .csv file containing at least two columns: "Sample", matching exactly the sample names in `datadir`, and "Condition", giving the experiment status of that sample, such as WT or KO, Case vs Control, etc. Optionally, can provide a third column "Code" giving a nickname for each sample; this is set to "Sample_Condition" for each sample if not.
#' @param comps string, path to a .csv file containing two columns called "c1" and "c2". Each row will be used to compare conditions from the "sample_metadata" csv; multiple comparisons are supported.
#' @param risc_reference string, name of sample to use as RISC reference sample, if not provided will automate the choice
#' @param min_num_UMI numeric, default is 500, if no filter is desired set to -Inf
#' @param min_num_Feature numeric, default is 200, if no filter is desired set to -Inf
#' @param max_perc_mito numeric, default is 25, if no filter is desired set to Inf
#' @param max_perc_hemoglobin numeric, default is 25, if no filter is desired set to Inf
#' @param autofilter_complexity T/F, default T, whether to filter cells with lower than expected number of genes given number of UMIs
#' @param autofilter_mito T/F, default T, whether to filter cells with higher than normal mito content
#' @param autofilter_nUMI T/F, default T, whether to filter cells with lower than normal UMI content
#' @param autofilter_medianabsolutedev_threshold numeric, default is 3, threshold for median abs deviation thresholding, ie cutoffs set to `median +/- mad * threshold`
#' @param autofilter_loess_negative_residual_threshold numeric, cutoff for loess residuals applied in complexity filtering, default is -5, if you set it high (ie any higher than -2) you will probably remove many good cells.
#' @param doubletFinder T/F, default is T, whether to filter doublets with `DoubletFinder`
#' @param pcs_indi integer, default = 30; number of PCs to use in individual sample processing / clustering
#' @param res_indi numeric, default = 0.5; Louvain resolution for individual sample clustering
#' @param pcs_int integer, default = 30; number of PCs to use in integrated data processing / clustering
#' @param res_int numeric, default = 0.5 ; Louvain resolution for louvain clustring of RISC integrated dataset; see `scDAPP::scCluster_louvain_res()`
#' @param RISC_louvain_neighbors integer, default = 10; number of nearest neighbors to consider during clustering; see `RISC::scCluster()` or `scDAPP::scCluster_louvain_res()` where implementation of this is unchanged
#' @param Pseudobulk_mode T/F. Sets the cross-conditional analysis mode. TRUE uses pseudobulk EdgeR for DE testing and propeller for compositional analysis. FALSE uses single-cell wilcox test within Seurat for DE testing and 2-prop Z test within the `prop.test()` function for compositional analysis.
#' @param DE_test a string, default is 'EdgeR-LRT' when Pseudobulk_mode is set to True, or 'wilcox' when Pseudobulk_mode is False. Can be either "DESeq2", "DESeq2-LRT", "EdgeR", "EdgeR-LRT" for pseudobulk, or any of the tests supported by the "test.use" argument in the FindMarkers function in Seurat; see `?Seurat::FindMarkers` for more.
#' @param crossconditionDE_padj_thres numeric, numeric; adjusted p value threshold for significant DE genes in cross condition DE; if `Pseudobulk_mode` is set to T default is 0.1; if `Pseudobulk_mode` is F default is 0.05
#' @param crossconditionDE_lfc_thres numeric, absolute value of LFC threshold for significant DE genes in cross condition DE; if `Pseudobulk_mode` is set to F default is 0 (no minimum lFC); if `Pseudobulk_mode` is F default is 0.25
#' @param pathway_padj_thres numeric, threshold for significant DE pathways via GSEA test; default is 0.1
#' @param species string, for example 'Homo sapiens' or 'Mus musculus', default = 'Homo sapiens'; this is for pathway analysis, see `msigdbr::msigdbr_species()`
#' @param workernum integer, number of CPU threads, default = 1
#' @param input_seurat_obj T/F. If true, will read in Seurat objects from `datadir` with names matching the sample column of `sample_metadata`. Ie, if datadir contains objects called "Sample1.rds", "Sample2.rds", and psuedobulk_metadata has "Sample1" in the Sample column, only Sample1.rds will be read in. Useful for data with some preprocessing or hashed data input.
#' @param title string, title of HTML report. Default is "10X analysis - clustering and integration".
#' @param author string, name of authors which will be shown on HTML report. We recommend passing a comma separated string. Default is "Alexander Ferrena, Deyou Zheng".
#' @param pseudobulk_metadata FOR BACKWARDS COMPATIBILITY ONLY. Will be set to `sample_metadata`. string, path to a .csv file containing at least two columns: "Sample", matching exactly the sample names in `datadir`, and "Condition", giving the experiment status of that sample, such as WT or KO, Case vs Control, etc. Optionally, can provide a third column "Code" giving a nickname for each sample; this is set to "Sample_Condition" for each sample if not.
#' @param de.test.use FOR BACKWARDS COMPATIBILITY ONLY. String, will set `Pseudobulk_mode` to TRUE if "pseudobulk_edgeR" is passed, or will set `Pseudobulk_mode` to FALSE if "wilcox" is passed.
#'
#' @return Returns the warnings from the pipeline
#' @export
#'
#' @examples
#' \dontrun{
#' # `sample_metadata` is a string path to a .csv file
#' # that looks like this, with "samples" matching
#' # the Spaceranger output folder names in `datadir`
#'
#' #sample_metadata.csv file
#' Sample,Condition
#' SampleXYZ1,Control
#' SampleXYZ2,Control
#' SampleABC1,KO1
#' SampleABC2,KO1
#' SampleJKL1,KO2
#' SampleJKL2,KO2
#'
#' # Optionally, can provide a third column
#' # to `sample_metadata` csv file like below:
#'
#' #sample_metadata.csv file with optional Code column
#' Sample,Condition,Code
#' SampleXYZ1,Control,Control1
#' SampleXYZ2,Control,Control2
#' SampleABC1,KO1,KO1_1
#' SampleABC2,KO1,KO1_2
#' SampleJKL1,KO2,KO2_1
#' SampleJKL2,KO2,KO2_1
#'
#'
#' `comps` leads to a .csv file that looks like this:
#'
#' #comps.csv file
#' c1,c2
#' KO1,Control
#' KO2,Control
#' KO1,K2
#'
#'
#' # Run pipeline like so:
#'
#' scRNAseq_pipeline_runner(datadir = 'path/to/cellranger/output/folder',
#'                outdir = 'path/to/output/folder',
#'                sample_metadata = 'path/to/sample_metadata.csv'
#'                comps = 'path/to/comps.csv',
#'                Pseudobulk_mode = T
#'                )
#'
#'
#' # Or with options, like so:
#'
#' scRNAseq_pipeline_runner(datadir = 'path/to/cellranger/output/folder',
#'                outdir = 'path/to/output/folder',
#'                sample_metadata = 'path/to/sample_metadata.csv'
#'                comps = 'path/to/comps.csv',
#'                Pseudobulk_mode = T,
#'                use_labeltransfer = T,
#'                refdatapath = 'path/to/referenceSeuratObject.rds',
#'                m_reference = 'path/to/reference_FindAllMarkers.rds',
#'                species = 'Mus musculus',
#'                risc_reference = 'SampleXYZ1',
#'                workernum = 4
#'                )
#'
#'
#' }
scRNAseq_pipeline_runner <- function(  datadir,
                                       outdir,
                                       use_labeltransfer,
                                       refdatapath,
                                       m_reference,
                                       
                                       sample_metadata,
                                       comps,
                                       
                                       Pseudobulk_mode,

                                       risc_reference,

                                       min_num_UMI,
                                       min_num_Feature,
                                       max_perc_mito,
                                       max_perc_hemoglobin,
                                       autofilter_complexity,
                                       autofilter_mito,
                                       autofilter_nUMI,
                                       autofilter_medianabsolutedev_threshold,
                                       autofilter_loess_negative_residual_threshold,
                                       doubletFinder,

                                       pcs_indi,
                                       res_indi,
                                       pcs_int,
                                       res_int,
                                       RISC_louvain_neighbors,

                                       
                                       DE_test,
                                       crossconditionDE_padj_thres,
                                       crossconditionDE_lfc_thres,
                                       pathway_padj_thres,
                                       species,
                                       workernum,

                                       input_seurat_obj,

                                       title,
                                       author,
                                       
                                       de.test.use,
                                       pseudobulk_metadata
                                       ){


  message('\n\nBegin pipeline\n\n')



  ### note: make sure to all params to here and to the rmd params section
  if(!missing(pseudobulk_metadata) & missing(sample_metadata) ){sample_metadata <- pseudobulk_metadata}
  if(!missing(de.test.use) & missing(Pseudobulk_mode) ){Pseudobulk_mode <- ifelse(de.test.use == 'pseudobulk_edgeR', yes= T, no = F) }
  
  if(missing(datadir)){datadir = NULL}
  if(missing(outdir)){outdir = NULL}
  if(missing(use_labeltransfer)){use_labeltransfer = FALSE}
  if(missing(refdatapath)){ refdatapath = NULL}
  if(missing(m_reference)){m_reference = NULL}
  # if(missing(SeuratLabelTransfer.normalization.method)){SeuratLabelTransfer.normalization.method = 'auto'}
  # SeuratLabelTransfer.normalization.method  string, either "auto", "SCT", or "LogNormalize". Default is "auto". This is passed to `Seurat::FindTransferAnchors()`. "SCT" is ideal; "auto" searches for "SCT" assay in reference and uses if detected. "LogNormalize" can be used if SCT is not possible, for example if raw counts are hard to get for a published dataset.



  if(missing(sample_metadata)){ sample_metadata = NULL}
  if(missing(comps)){ comps = NULL}

  if(missing(risc_reference)){ risc_reference =  NULL}

  if(missing(min_num_UMI)){ min_num_UMI =  500}
  if(missing(min_num_Feature)){ min_num_Feature =  200}
  if(missing(max_perc_mito)){ max_perc_mito =  25}
  if(missing(max_perc_hemoglobin)){ max_perc_hemoglobin =  25}
  if(missing(autofilter_complexity)){ autofilter_complexity =  TRUE}
  if(missing(autofilter_mito)){ autofilter_mito =  TRUE}
  if(missing(autofilter_nUMI)){ autofilter_nUMI =  TRUE}
  if(missing(autofilter_medianabsolutedev_threshold)){ autofilter_medianabsolutedev_threshold =  3}
  if(missing(autofilter_loess_negative_residual_threshold)){ autofilter_loess_negative_residual_threshold =  -5}
  if(missing(doubletFinder)){ doubletFinder =  TRUE}


  if(missing(pcs_indi)){pcs_indi =  30}
  if(missing(res_indi)){res_indi = 0.5}
  if(missing(pcs_int)){ pcs_int = 30}
  if(missing(res_int)){ res_int = 0.5}
  if(missing(RISC_louvain_neighbors)){ RISC_louvain_neighbors = 10}


  if(missing(DE_test)){
    if(Pseudobulk_mode == T){DE_test = 'EdgeR-LRT'}
    if(Pseudobulk_mode == F){DE_test = 'wilcox'}
  }
  if(missing(crossconditionDE_padj_thres)){ crossconditionDE_padj_thres = NULL}
  if(missing(crossconditionDE_lfc_thres)){ crossconditionDE_lfc_thres = NULL}
  if(missing(pathway_padj_thres)){ pathway_padj_thres = 0.1}
  if(missing(species)){ species = 'Homo sapiens'}
  if(missing(workernum)){ workernum = 1}

  if(missing(input_seurat_obj)){ input_seurat_obj = FALSE}

  if(missing(title)){title = 'scDAPP Report'}
  if(missing(author)){author = 'Pipeline prepared by Alexander Ferrena, Deyou Zheng, and colleagues'}



  #locate the pipeline file
  rmdfile <- system.file("rmd", "scRNAseq_clustering_integration.Rmd", package = "scDAPP")

  message('Will run rmd file at:\n',
          rmdfile,
          '\n\n')


  rmarkdown::render(rmdfile,
                    params=list(
                      datadir = datadir,
                      outdir = outdir,
                      use_labeltransfer = use_labeltransfer,
                      refdatapath = refdatapath,
                      m_reference = m_reference,
                      
                      sample_metadata = sample_metadata,
                      comps = comps,

                      risc_reference = risc_reference,

                      min_num_UMI = min_num_UMI,
                      min_num_Feature = min_num_Feature,
                      max_perc_mito = max_perc_mito,
                      max_perc_hemoglobin = max_perc_hemoglobin,
                      autofilter_complexity = autofilter_complexity,
                      autofilter_mito = autofilter_mito,
                      autofilter_nUMI = autofilter_nUMI,
                      autofilter_medianabsolutedev_threshold = autofilter_medianabsolutedev_threshold,
                      autofilter_loess_negative_residual_threshold = autofilter_loess_negative_residual_threshold,
                      doubletFinder = doubletFinder,

                      pcs_indi = pcs_indi,
                      res_indi = res_indi,
                      pcs_int = pcs_int,
                      res_int = res_int,
                      RISC_louvain_neighbors = RISC_louvain_neighbors,

                      Pseudobulk_mode = Pseudobulk_mode,
                      DE_test = DE_test,
                      crossconditionDE_padj_thres = crossconditionDE_padj_thres,
                      crossconditionDE_lfc_thres = crossconditionDE_lfc_thres,
                      pathway_padj_thres = pathway_padj_thres,
                      species = species,
                      workernum = workernum,
                      input_seurat_obj = input_seurat_obj,

                      title = title,
                      author = author,

                      force_redo = FALSE #maybe in future...
                    ),

                    #this line ensures html prints to outdir folder
                    output_dir = outdir

  )



  message('\n\nPipeline completed!\n\n')

  return(warnings())


}






