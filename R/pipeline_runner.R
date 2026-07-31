# https://r-pkgs.org/whole-game.html

#' Test installation of pipeline packages
#'
#' Attaches dependencies via [attach_scDAPP_pipeline_libraries()] (including RISC)
#' and returns a data frame of package versions. Use
#' [attach_scDAPP_pipeline_libraries()] directly in scripts and the Rmd when you
#' only need to load libraries without printing versions.
#'
#' @return A `data.frame` with columns `pkg` and `vers`.
#' @export
#'
#' @examples
#' \dontrun{
#' scDAPP::r_package_test()
#' }
r_package_test <- function() {
  pkg_names <- .pipeline_version_check_pkgs()

  packages <- data.frame(pkg = pkg_names, stringsAsFactors = FALSE)

  packages$vers <- vapply(
    packages$pkg,
    function(pkg) {
      tryCatch(
        as.character(utils::packageVersion(pkg)),
        error = function(cond) {
          warning('package "', pkg, '" not detected!!', call. = FALSE)
          NA_character_
        }
      )
    },
    character(1)
  )

  scDAPP::attach_scDAPP_pipeline_libraries(load_risc = TRUE, quietly = TRUE)

  packages
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
#' @param comps string, path to a .csv file with columns c0 (reference), c1 (test), and optional formula, contrast, label. Legacy c2 is accepted as alias for c0. Multiple comparisons are supported.
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
#' @param cluster_unfiltered logical, default FALSE. If TRUE, run SCT and Louvain clustering
#'   on the full unfiltered matrix before autofilter (for cluster-level QC diagnostics).
#' @param pcs_indi integer, default = 30; number of PCs to use in individual sample processing / clustering
#' @param res_indi numeric, default = 0.5; Louvain resolution for individual sample clustering
#' @param pcs_int integer or `"auto"`, default = 30; number of PCs for integrated clustering, or automatic selection via bootstrap cluster stability
#' @param res_int numeric or `"auto"`, default = 0.5; Louvain resolution for integrated clustering, or `"auto"` (see `auto_integration_cluster_params()`)
#' @param stability_outdir optional path for stability sweep outputs when `pcs_int` or `res_int` is `"auto"` (default: `outdir/integrated_analysis/cluster_stability`)
#' @param stability_numreps bootstrap replicates for auto tuning (default 50)
#' @param stability_sweep_maxPCs PC grid when `pcs_int = "auto"`
#' @param stability_sweep_res resolution grid when `res_int = "auto"`
#' @param stability_propcells.perrep cell fraction per bootstrap replicate (default 0.8)
#' @param RISC_louvain_neighbors integer, default = 10; number of nearest neighbors to consider during clustering; see `RISC::scCluster()` or `scDAPP::scCluster_louvain_res()` where implementation of this is unchanged
#' @param integration_method character; sample integration backend. One of `RISC` (default) or Seurat v5 `IntegrateLayers` methods (`CCAIntegration`, `RPCAIntegration`, `CCAIntegration_SCT`, `RPCAIntegration_SCT`, `HarmonyIntegration`). See `integration_method_choices()`.
#' @param Pseudobulk_mode T/F. Sets the cross-conditional analysis mode. TRUE uses pseudobulk EdgeR (or Dream) for DE testing and propeller for compositional analysis. FALSE uses single-cell wilcox test within Seurat for DE testing and 2-prop Z test within the `prop.test()` function for compositional analysis.
#' @param DE_test a string, default is 'EdgeR-LRT' when Pseudobulk_mode is set to True, or 'wilcox' when Pseudobulk_mode is False. Can be "DESeq2", "DESeq2-LRT", "EdgeR", "EdgeR-LRT", or "Dream" for pseudobulk (Dream requires Bioconductor `variancePartition` and a `(1|var)` term in `comps$formula`, e.g. `~ Condition + (1|Patient)`), or any of the tests supported by the "test.use" argument in the FindMarkers function in Seurat; see `?Seurat::FindMarkers` for more. Note the Seurat "roc" test is not included, and some additional packages like DESeq2 may require installation.
#' @param crossconditionDE_padj_thres numeric, numeric; adjusted p value threshold for significant DE genes in cross condition DE; if `Pseudobulk_mode` is set to T default is 0.1; if `Pseudobulk_mode` is F default is 0.05
#' @param crossconditionDE_lfc_thres numeric, absolute value of LFC threshold for significant DE genes in cross condition DE; if `Pseudobulk_mode` is T default is 0 (no minimum LFC); if `Pseudobulk_mode` is F default is 0.25
#' @param crossconditionDE_min.pct numeric, minimum expression fraction for DEG counting and ORA (`pct.1` if up, `pct.2` if down); if `Pseudobulk_mode` is T default is 0.1; if F default is 0. Pass `NULL` to use mode defaults via `crosscondition_de_threshold_defaults()`.
#' @param pathway_padj_thres numeric, threshold for significant DE pathways via GSEA test; default is 0.1
#' @param species string, for example 'Homo sapiens' or 'Mus musculus', default = 'Homo sapiens'; this is for pathway analysis, see `msigdbr::msigdbr_species()`
#' @param msigdbr_cache_dir optional string, directory for cached MSigDB pathway tables prepared by `preppathways_pathwayanalysis_crosscondition_module()`. When NULL, uses `XDG_CACHE_HOME/scDAPP` if set, else `tools::R_user_dir("scDAPP", "cache")`. Falls back to a subfolder of the pipeline output directory with a warning if those locations are not writable.
#' @param workernum integer, number of CPU threads for final RISC integration and Dream pseudobulk DE (`BiocParallel`), default = 1
#' @param stability_workernum integer or NULL; parallel workers for auto stability sweep (NULL uses workernum)
#' @param run_ORA T/F, default is F. Whether to run OverRepresentation Analysis (ORA) using fisher exact tests as implemented in `clusterProfiler::enricher()`. clusterProfiler must be installed for this. Will save table outputs.
#' @param run_msigdb_celltype_ora T/F, default is TRUE. Whether to run ORA of cluster markers against MSigDB cell-type signature gene sets on individual and integrated marker tables. clusterProfiler must be installed.
#' @param msigdb_celltype_ora_marker_padj_thres numeric, adjusted p-value cutoff for cluster marker genes included in MSigDB cell-type ORA (default 0.05).
#' @param msigdb_celltype_ora_top_markers integer, max markers per cluster after padj filter, ranked by score (default 100).
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
                                       cluster_unfiltered,
                                       
                                       pcs_indi,
                                       res_indi,
                                       pcs_int,
                                       res_int,
                                       stability_outdir,
                                       stability_numreps,
                                       stability_sweep_maxPCs,
                                       stability_sweep_res,
                                       stability_propcells.perrep,
                                       RISC_louvain_neighbors,
                                       integration_method,
                                       
                                       DE_test,
                                       crossconditionDE_padj_thres,
                                       crossconditionDE_lfc_thres,
                                       crossconditionDE_min.pct,
                                       pathway_padj_thres,
                                       species,
                                       msigdbr_cache_dir,
                                       workernum,
                                       stability_workernum,
                                       run_ORA,
                                       run_msigdb_celltype_ora,
                                       msigdb_celltype_ora_marker_padj_thres,
                                       msigdb_celltype_ora_top_markers,
                                       
                                       input_seurat_obj,
                                       
                                       title,
                                       author,
                                       
                                       de.test.use,
                                       pseudobulk_metadata
){
  
  
  message('\n\nBegin pipeline\n\n')

  scDAPP::set_parallel_blas_threads()
  
  
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
  if(missing(cluster_unfiltered)){ cluster_unfiltered = FALSE}
  
  
  if(missing(pcs_indi)){pcs_indi =  30}
  if(missing(res_indi)){res_indi = 0.5}
  if(missing(pcs_int)){ pcs_int = 30}
  if(missing(res_int)){ res_int = 0.5}
  if(missing(stability_outdir)){ stability_outdir <- NULL}
  if(missing(stability_numreps)){ stability_numreps <- 50L}
  if(missing(stability_sweep_maxPCs)){
    stability_sweep_maxPCs <- c(5, 10, 15, 20, 25, 30, 40, 50)
  }
  if(missing(stability_sweep_res)){ stability_sweep_res <- seq(0.1, 1.5, by = 0.2)}
  if(missing(stability_propcells.perrep)){ stability_propcells.perrep <- 0.8}
  if(missing(RISC_louvain_neighbors)){ RISC_louvain_neighbors = 10}
  if(missing(integration_method)){ integration_method = "RISC"}
  integration_method <- match.arg(integration_method, scDAPP::integration_method_choices())
  
  
  if(missing(DE_test)){
    if(Pseudobulk_mode == T){DE_test = 'EdgeR-LRT'}
    if(Pseudobulk_mode == F){DE_test = 'wilcox'}
  }
  if(missing(crossconditionDE_padj_thres)){ crossconditionDE_padj_thres = NULL}
  if(missing(crossconditionDE_lfc_thres)){ crossconditionDE_lfc_thres = NULL}
  if(missing(crossconditionDE_min.pct)){ crossconditionDE_min.pct = NULL}
  if(missing(pathway_padj_thres)){ pathway_padj_thres = 0.1}
  if(missing(species)){ species = 'Homo sapiens'}
  if(missing(msigdbr_cache_dir)){ msigdbr_cache_dir = NULL}
  if(missing(workernum)){ workernum = 1}
  if(missing(stability_workernum)){ stability_workernum <- NULL}
  if(missing(run_ORA)){ run_ORA = F }
  if(missing(run_msigdb_celltype_ora)){ run_msigdb_celltype_ora = TRUE }
  if(missing(msigdb_celltype_ora_marker_padj_thres)){ msigdb_celltype_ora_marker_padj_thres = 0.05 }
  if(missing(msigdb_celltype_ora_top_markers)){ msigdb_celltype_ora_top_markers = 100L }
  
  if(missing(input_seurat_obj)){ input_seurat_obj = FALSE}
  
  if(missing(title)){title = 'scDAPP Report'}
  if(missing(author)){author = 'Pipeline prepared by Alexander Ferrena, Deyou Zheng, and colleagues'}
  
  
  

  
  
  
  #### write some exceptions 
  
  #DE tests must be in a set of tests
  if(Pseudobulk_mode == T){
    
    if(!DE_test %in% c('EdgeR', 'EdgeR-LRT', 'DESeq2', 'DESeq2-LRT', 'Dream')){
      stop("With 'Pseudobulk_mode' set to T, DE_test must be one of: 'EdgeR', 'EdgeR-LRT', 'DESeq2', 'DESeq2-LRT', 'Dream'; value ", DE_test, " was passed")
    }
    
  }
  
  #DE tests must be in a set of tests
  if(Pseudobulk_mode == F){
    
    if(!DE_test %in% c('wilcox' , 'wilcox_limma', 'bimod', 't', 'negbinom', 'poisson', 'LR', 'MAST', 'DESeq2' )){
      stop("With 'Pseudobulk_mode' set to F, DE_test must be one of: 'wilcox' , 'wilcox_limma', 'bimod', 't', 'negbinom', 'poisson', 'LR', 'MAST', 'DESeq2'; value ", DE_test, " was passed")
    }
    
  }
  
  
  
  
  # make sure DE packages are installed
  if(DE_test == 'MAST'){
    if(!('MAST' %in% rownames(installed.packages()))){
      stop('DE_test is set to "MAST". Please install MAST first.')
    }
  }
  
  if(DE_test == 'DESeq2' | DE_test == 'DESeq2-LRT'){
    if(!('DESeq2' %in% rownames(installed.packages()))){
      stop('DE_test is set to "DESeq2". Please install DESeq2 first')
    }
  }

  if (identical(DE_test, "Dream")) {
    if (!requireNamespace("variancePartition", quietly = TRUE)) {
      stop(
        'DE_test is set to "Dream". Please install variancePartition first ',
        "(BiocManager::install(\"variancePartition\"))."
      )
    }
  }

  # Paired formulas (1|var) must match DE_test
  if (isTRUE(Pseudobulk_mode) && !is.null(comps)) {
    comps_df <- if (is.character(comps) && length(comps) == 1L && file.exists(comps)) {
      utils::read.csv(comps, stringsAsFactors = FALSE)
    } else {
      comps
    }
    .validate_comps_de_formulas(comps_df, DE_test)
  }
  
  
  #if using run_ORA make sure clusterProfiler is installed
  if(run_ORA == T){
    if(!('clusterProfiler' %in% rownames(installed.packages()))){
      stop('run_ORA is set to T. Please install clusterProfiler first')
    }
  }

  if(run_msigdb_celltype_ora == TRUE){
    if(!('clusterProfiler' %in% rownames(installed.packages()))){
      stop('run_msigdb_celltype_ora is set to TRUE. Please install clusterProfiler first')
    }
  }
  
  ####
  
  
  
  #locate the pipeline file
  rmdfile <- system.file("rmd", "scRNAseq_clustering_integration.Rmd", package = "scDAPP")
  
  message('Found rmd file at:\n',
          rmdfile,
          '\n\n')
  
  ## update 2025.03.18; try copying the rmd file to outdir first and rendering there
  # this is to make things work nicely with singularity
  # tmp_rmd <- tempfile(tmpdir = outdir, fileext = ".Rmd")
  tmp_rmd <- paste0(outdir, '/', basename(rmdfile))

  message('Placing runnable copy of .rmd file in outdir:\n',
          tmp_rmd,
          '\n\n')
  dir.create(outdir, recursive = T)
  file.copy(rmdfile, tmp_rmd, overwrite = TRUE)
  
  
  
  rmarkdown::render(tmp_rmd,
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
                      cluster_unfiltered = cluster_unfiltered,
                      
                      pcs_indi = pcs_indi,
                      res_indi = res_indi,
                      pcs_int = pcs_int,
                      res_int = res_int,
                      stability_outdir = stability_outdir,
                      stability_numreps = stability_numreps,
                      stability_sweep_maxPCs = stability_sweep_maxPCs,
                      stability_sweep_res = stability_sweep_res,
                      stability_propcells.perrep = stability_propcells.perrep,
                      RISC_louvain_neighbors = RISC_louvain_neighbors,
                      integration_method = integration_method,
                      
                      Pseudobulk_mode = Pseudobulk_mode,
                      DE_test = DE_test,
                      crossconditionDE_padj_thres = crossconditionDE_padj_thres,
                      crossconditionDE_lfc_thres = crossconditionDE_lfc_thres,
                      crossconditionDE_min.pct = crossconditionDE_min.pct,
                      pathway_padj_thres = pathway_padj_thres,
                      species = species,
                      msigdbr_cache_dir = msigdbr_cache_dir,
                      workernum = workernum,
                      stability_workernum = stability_workernum,
                      run_ORA = run_ORA,
                      run_msigdb_celltype_ora = run_msigdb_celltype_ora,
                      msigdb_celltype_ora_marker_padj_thres = msigdb_celltype_ora_marker_padj_thres,
                      msigdb_celltype_ora_top_markers = msigdb_celltype_ora_top_markers,
                      
                      input_seurat_obj = input_seurat_obj,
                      
                      title = title,
                      author = author,
                      
                      force_redo = FALSE #maybe in future...
                    ),
                    
                    #this line ensures html prints to outdir folder
                    output_dir = outdir,
                    
                    #update 2025.03.18 below lines will hopefully write tmp files to outdir to help with singularity read-only restrictions
                    intermediates_dir = outdir,
                    knit_root_dir = outdir
                    
  )
  
  
  
  message('\n\nPipeline completed!\n\n')
  
  return(warnings())
  
  
}






