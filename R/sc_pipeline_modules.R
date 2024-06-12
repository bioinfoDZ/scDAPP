# https://r-pkgs.org/whole-game.html

#' Differential expression analysis across conditions for integrated Seurat objects
#'
#' This is a modular component of the scDAPP scRNAseq pipeline. Perform DE analysis across conditions. Supports A vs B vs C pairwise (multiple conditions) comparisons. Options for Pseudobulk DE via EdgeR - LRT, or old-school scRNAseq DE via wilcoxon test.
#'
#' @param sobjint integrated Seurat object. metadata needs two special columns: one called "Condition" that contains the A vs B conditions, and a second that matches the  `grouping_variable` parameter of this function.
#' @param sample_metadata data.frame with three columns called Sample, Condition, Code.
#' @param comps data.frame with two columns called c1, c2; these are the conditions you want to set up. Should be present in psuedobulk_metadata Condition column.
#' @param grouping_variable string, column name of identity in Seurat object meta.data to stratify DE by. For example, clusters or celltype. Will perform A vs B DE in each of these groupings. Default is "seurat_clusters". It is not mandatory, but will use factor level ordering of this variable in the meta.data to control analysis order, and if not will sort by alphanumeric order (cluster 1, then 2, cluster A, then B, etc)
#' @param Pseudobulk_mode T/F. Sets the cross-conditional analysis mode. TRUE uses pseudobulk EdgeR for DE testing and propeller for compositional analysis. FALSE uses single-cell wilcox test within Seurat for DE testing and 2-prop Z test within the `prop.test()` function for compositional analysis.
#' @param DE_test a string, default is 'EdgeR-LRT' when Pseudobulk_mode is set to True, or 'wilcox' when Pseudobulk_mode is False. Can be either "DESeq2", "DESeq2-LRT", "EdgeR", "EdgeR-LRT" for pseudobulk, or any of the tests supported by the "test.use" argument in the FindMarkers function in Seurat; see `?Seurat::FindMarkers` for more. Note the Seurat "roc" test is not included, and some additional packages like DESeq2 may require installation.
#' @param outdir_int string, path to save results to, if not provided will not save. Will create a sub-directory called "differentialexpression_crosscondition" and save inside of there.
#' @param assay string, name of Seurat assay to use, default is DefaultAssay(sobjint)
#' @param slot string, name of Seurat assay slot to use, default is 'data'
#' @param cluster_prefix T/F, whether to append prefix "cluster_" to grouping_variable levels, useful if grouping_variable is a cluster. Rather than saving results with names such as "1", "2", "3", will save as "cluster_1", and so on.
#' @param crossconditionDE_padj_thres numeric; adjusted P value maximum threshold for calling DEGs. Only used in counting number of DEGs. Default is 0.1
#' @param crossconditionDE_lfc_thres numeric; Log Fold Change magnitude (absolute value) minimum threshold for calling DEGs. Only used in counting number of DEGs. Default is 0.
#' @param crossconditionDE_min.pct numeric; Minimum percent of cells expressing gene required to count as a DEG. For positive LFC genes (up in condition A); pct.1 must be at least this value (percent of cells in A must be at least this value); for negative LFC genes, pct.2 cells must be at least this value. Only used for counting DEGs. Default is 0.1 if pseudobulk_edgeR is used; 0 if wilcox is used.
#'
#' @return a complex nested list with several layers. Nest level 1: comparisons performed, as given by comps. Level 2: data.frames with A vs B results stratified by cluster. Note that "Weight" will be -log10(pval) * sign of L2FC.
#' @export
#'
#' @examples
#' \dontrun{
#' # `sample_metadata` looks like this:
#' Sample,Condition,Code
#' SampleXYZ1,Control,Control1
#' SampleXYZ2,Control,Control2
#' SampleABC1,KO1,KO1_1
#' SampleABC2,KO1,KO1_2
#' SampleJKL1,KO2,KO2_1
#' SampleJKL2,KO2,KO2_1
#'
#'
#' # `comps` looks like this:
#' c1,c2
#' KO1,Control
#' KO2,Control
#' KO1,K2
#'
#' # Run DE analysis:
#' m_bycluster_crosscondition_de_comps <- de_across_conditions_module(
#'  sobjint = sobjint,
#'  sample_metadata = sample_metadata,
#'  comps = comps,
#'  outdir_int = outdir_int,
#'  grouping_variable = 'seurat_clusters',
#'  Pseudobulk_mode = T
#'  )
#'
#'  # check outputs:
#'  # Level 1 will be comparisons (A vs B; B vs C etc)
#'  names(m_bycluster_crosscondition_de_comps)
#'  "KO_vs_WT"
#'
#'  # Level 2 will be DE results of each grouping_variable level stratified by cluster
#'  m_bycluster_crosscondition_de_comps$KO_vs_WT
#'  "cluster_1"  "cluster_2"  "cluster_3"  "cluster_4"  "cluster_5"
#'
#'  #each of these is a data.frame with the DE results for this comparison in this cluster
#'  head(m_bycluster_crosscondition_de_comps$KO_vs_WT$cluster_1)
#'
#' )
#'
#' }
de_across_conditions_module <- function(sobjint,
                                        sample_metadata,
                                        comps,
                                        grouping_variable,
                                        Pseudobulk_mode,
                                        DE_test,
                                        outdir_int,
                                        assay,
                                        slot,
                                        cluster_prefix,
                                        crossconditionDE_padj_thres,
                                        crossconditionDE_lfc_thres,
                                        crossconditionDE_min.pct
                                        
){
  
  
  if( missing(sobjint)) { stop('Provide Seurat object') }
  if( missing(grouping_variable)) { stop(grouping_variable <- 'seurat_clusters') }
  # if( missing(min.pct)) { stop(min.pct <- 0.1) } # for GSEA, don't do this
  
  if(missing(assay)){assay <- DefaultAssay(sobjint)}
  if(missing(slot)){slot <- 'data'}
  if(missing(cluster_prefix)){cluster_prefix <- NULL}
  
  if(missing(crossconditionDE_padj_thres)){ crossconditionDE_padj_thres = 0.1 }
  if(missing(crossconditionDE_lfc_thres)){crossconditionDE_lfc_thres = 0}
  if(missing(crossconditionDE_min.pct)){
    if(Pseudobulk_mode == T){crossconditionDE_min.pct = 0.1} else{
      crossconditionDE_min.pct = 0
    }
    
  }
  
  if(missing(DE_test)){
    if(Pseudobulk_mode == T){DE_test = 'EdgeR-LRT'}
    if(Pseudobulk_mode == F){DE_test = 'wilcox'}
  }
  
  
  
  
  
  
  #prep names
  comps$labels <- paste0(comps$c1, '_vs_', comps$c2)
  
  #read sobjlist back in? keep it in?
  # will need to optimize memory
  
  #get cluster object name
  # grouping_variable <- risc_clust_lab
  
  
  #check if factor
  groupingvec <- sobjint@meta.data[,grouping_variable]
  if(!is.factor(groupingvec)){
    warning('The grouping variable (ie clusters within which to perform A vs B DE) is not a factor.\nThe the comparison order will default to alphanumeric order.')
    
    groupinglevs <- unique(groupingvec)
    groupinglevs <- stringr::str_sort(groupinglevs, numeric=T)
    groupingvec <- factor(groupingvec, levels = groupinglevs)
    
  }
  
  
  #get the actual clusters (grouping levels)
  groupinglevs <- levels( groupingvec )
  
  #remove empty levels
  groupinglevs <- groupinglevs[groupinglevs %in% groupingvec]
  
  # for later, if they are clusters, we want better labels than just numerics
  # try to check if they are clusters; max str len will probably be 3 (in huge datasets...)
  if( is.null(cluster_prefix) ) {
    
    if( max(str_length(groupinglevs)) <= 3 ){cluster_prefix <- T}
    
  }
  
  if(cluster_prefix==T){
    groupinglev_nicelabs <- paste0('cluster_', groupinglevs)
  } else{groupinglev_nicelabs <- groupinglevs}
  
  
  ### DE: pseudobulk or single-cell
  
  if(Pseudobulk_mode == T){
    
    
    # how to deal with clusters that are missing from samples?
    # pseudobulk all even if too few cells
    # later if all 0 just remove the column from analysis
    # if this removes too many sampless, make sure we can still run at least 2 v 2 by condition
    
    
    
    #for sample, pseudobulk by celltype
    samples <- sample_metadata$Code
    samp <- samples[1]
    
    pblist <- lapply(samples, function(samp){
      
      # message('\n\n',samp, '\n')
      
      md <- sobjint@meta.data
      md <- md[md$Code == samp,]
      cells <- rownames(md)
      sobjint_ct <- sobjint[,cells]
      
      #try to unlog the RISC counts...
      sobjint_ct@assays$RISC@data <- expm1(sobjint_ct@assays$RISC@data)
      
      
      
      suppressMessages(
        pb <- scDAPP::pseudobulk(sobjint_ct,
                                 assay = assay,
                                 slot = slot,
                                 grouping_colname_in_md = grouping_variable,
                                 min_cells = 0)
      )
      
      
      ### round it for edgeR
      pb <- round(pb)
      
      pb
      
    })
    
    
    names(pblist) <- samples
    
    
    ### make sure cluster is in each sample, by adding fake column if needed...
    
    pblist <- lapply(pblist, function(pb){
      
      
      #if any cluster is missing,
      # loop thru missing clusters and create columns of 0s
      if(any(!(groupinglevs %in% colnames(pb)))){
        fakectcols <- lapply(groupinglevs, function(ct){
          if(!(ct %in% colnames(pb))){
            fakectcol <- data.frame(ct = rep(0, nrow(pb)))
            colnames(fakectcol) <- ct
            fakectcol
          }
        })
        fakectcolsdf <- dplyr::bind_cols(fakectcols)
        
        #add the columns of zeros to the gem
        pb <- cbind(pb, fakectcolsdf)
      }
      
      
      
      #make sure the clusters are ordered properly
      # (ie if we ahve 8 clusters and cluster 5 was missing)
      pb <- pb[,match(groupinglevs,colnames(pb))]
      
      pb
      
    })
    
    
    # save it as pblist overall, since it gets subsetted in the lapply
    pblist_overall <- pblist
    
    
    ### loop thru comparisons ###
    # DE in each comparison
    # make sure to select clusters shared by the two conditons
    # for those clusters, use all samples for edgeR, and use c1 vs c2 for contrast functon
    
    
    compidx = 1 # for testing
    
    compslen <- 1:nrow(comps)
    m_bycluster_crosscondition_de_comps <- lapply(compslen, function(compidx){
      
      
      
      #get comparison condition levels
      c1 <- comps[compidx,1]
      c2 <- comps[compidx,2]
      
      
      #get lab
      lab <- comps[compidx, 3]
      
      message('\n', lab)
      
      
      ## get only this comp samples
      #subset MD
      comp_pseudobulk_md <- sample_metadata[sample_metadata$Condition %in% c(c1,c2),]
      
      #for EDGER, c1 needs to be level 2, c2 needs to be level 1
      comp_pseudobulk_md$Condition <- factor(comp_pseudobulk_md$Condition, levels = c(c2, c1))
      
      #subet pblist for this comp
      pblist <- pblist_overall
      pblist <- pblist[ match(comp_pseudobulk_md$Code, names(pblist) )]
      
      
      
      
      
      
      #### for each int cluster, loop thru and compare condition 1 vs condition 2 ####
      
      ## loop thru SHARED clusters ##
      clusters <- groupinglevs
      names(clusters) <- clusters
      
      clust = clusters[1] # for testing
      
      
      m_bycluster_crosscondition_de <- lapply(clusters, function(clust){
        
        
        message(clust)
        
        #get the pseudobulks of each cluster
        
        gemlist <- lapply( names(pblist) ,  function(samp){
          
          # message(samp)
          pb <- pblist[[samp]]
          
          pbcol <- pb[,colnames(pb)==clust, drop=F]
          colnames(pbcol) <- samp
          
          pbcol
        })
        
        #combine the cluster pseudobulks to one gene exp matrix
        gem <- dplyr::bind_cols(gemlist)
        
        
        
        #prep for edgeR
        
        # remove low exp genes
        gem <- gem[Matrix::rowSums(gem) > 3,]
        
        # remove empty columns (samples with all 0 for this cell type)
        gem <- gem[,Matrix::colSums(gem) > 10,drop=F]
        
        #make sure enough samples from each condition remain... if not need to remove this one...
        
        comp_pseudobulk_md <- comp_pseudobulk_md[match(colnames(gem), comp_pseudobulk_md$Code),]
        
        #check or skip
        condtab <- table(comp_pseudobulk_md$Condition)
        if(condtab[c1] < 2 | condtab[c2] < 2){
          return()
        }
        
        
        # COUNT CELLS
        #subset object, to count num cells c1 vs c2
        md <- sobjint@meta.data
        md <- md[md[,grouping_variable] == clust,]
        md <- md[md$Condition %in% comp_pseudobulk_md$Condition,]
        sobjint_comp_ct <- sobjint[,rownames(md)]
        mat <- sobjint_comp_ct@assays$RISC@data
        mat <- expm1(mat)
        
        
        #split to c1 and c2, and count
        cells_c1 <- rownames( md[md$Condition == c1,] )
        mat_c1 <- mat[,colnames(mat) %in% cells_c1]
        cellsexp_c1 <- tabulate(mat_c1@i + 1L, nrow(mat_c1) )
        
        # pct.1 <- rowSums(x = mat_c1 > 0) / length(x = cells_c1)
        
        
        cells_c2 <- rownames( md[md$Condition == c2,] )
        mat_c2 <- mat[,colnames(mat) %in% cells_c2]
        cellsexp_c2 <- tabulate(mat_c2@i + 1L, nrow(mat_c2) )
        
        #  pct.2 <- rowSums(x = mat_c2 > 0) / length(x = cells_c2)
        
        
        #convert to proportions
        cellsexp_c1 <- cellsexp_c1 / ncol(mat_c1)
        cellsexp_c2 <- cellsexp_c2 / ncol(mat_c2)
        
        cellsexp <- data.frame(gene = rownames(mat),
                               pct.1 = cellsexp_c1,
                               pct.2 = cellsexp_c2)
        
        #get difference
        cellsexp$pct.diff <- cellsexp$pct.1 - cellsexp$pct.2
        
        ### select only genes exp in at least min.pct in either c1 or c2 ###
        # actually this removes a lot of genes and may violate GSEA, so we'll calculate without it and leave it up to user.
        # will also implement filtering and ORA later.
        # cellsexp <- cellsexp[cellsexp$pct.1 > min.pct | cellsexp$pct.2 > min.pct,]
        cellsexp <- cellsexp[cellsexp$gene %in% rownames(gem),]
        gem <- gem[match(cellsexp$gene, rownames(gem)),]
        
        
        
        
        
        
        
        
        
        
        ### EdgeR
        
        if(DE_test == 'EdgeR' | DE_test == 'EdgeR-LRT'){
          
          require(edgeR)
          
          
          #counts and "group"
          eobj <- DGEList(counts = gem, group = comp_pseudobulk_md$Condition)
          
          #size factors
          eobj <- calcNormFactors(eobj)
          
          #design, using group variable, factor levels are important
          
          ## Note, this may be passed as a user parameter...
          #FORMULA PARAMETER MAY BE IMPLEMENTED LATER
          design <- model.matrix(~comp_pseudobulk_md$Condition)
          
          
          
          
          #dispersion
          eobj <- estimateDisp(eobj, design)
          
          
          ### run the test
          
          #normal edgeR, added Jun 5 2024 as per reviewer request
          if(DE_test == 'EdgeR'){
            
            
            #run exact test as per qCML method, "vanilla" EdgeR
            et <- exactTest(eobj)
            
            #get res
            res <- as.data.frame ( topTags(et, n = Inf) )
            
            
          }
          
          #edgeR-LRT, only original test for pseudobulk comparisons based on Squair et al pseudobulk benchmark
          if(DE_test == 'EdgeR-LRT'){
            
            # coef = 2 is higher factor level, which should be c1
            fit <- glmFit(eobj,design)
            
            ## when using with user-supplied formula, we need to find the column name with Conditon... this can get complicated...
            # coef = grep(colnames(design), pattern = 'Condition') #not fully implemented currently
            #FORMULA PARAMETER MAY BE IMPLEMENTED LATER
            lrt <- glmLRT(fit,coef=2)
            
            #get res
            res <- as.data.frame ( topTags(lrt, n = Inf) )
            
          }
          
          
          
          
          #cbind pct.1 and pct.2
          cellsexp <- cellsexp[match(rownames(res), cellsexp$gene),]
          res <- cbind(res, cellsexp[,-1])
          
          
          
          #attach gene name as a column
          res <- cbind(rownames(res), res)
          colnames(res)[1] <- 'gene_symbol'
          
          
          
          
          
          # #add weight:
          # # -log10 pvalue * sign LFC * abs value of percent difference
          # # ie, significe of DE * sign of DE * num cells exp gene in that direction
          # ## downweight mismatch sign vs pct diff genes... do this by squring them, which makes decimal smaller ##
          # # add + 1 to abs pct diff --> do this to keep FGSEA scores high, otherwise we are actually dividing them by pct diff
          # pctdiff <- res$pct.diff
          # pctdiff[sign(pctdiff) != sign(res$logFC)] <- pctdiff[sign(pctdiff) != sign(res$logFC)] ^ 2
          # pctdiff <- abs(pctdiff) + 1
          # res$weight <- -log10(res$PValue) * res$logFC * pctdiff
          
          ## UPDATE DEC 7 2023, WEIGHT BY -LOG10(PVAL) * SIGN OF LFC
          
          ## prep weighted list ##
          
          #we have to deal with underflow...
          # get -log10 pvalues, sort
          scores <- -log10(res$PValue)
          scores <- scores * sign(res$logFC)
          names(scores) <- res$gene_symbol
          
          #sort by log pval with names
          scores <- sort(scores,decreasing = T)
          
          
          #also get logFC vector; for INf, we will sort them by LFC...
          logFC_vec <- res$logFC; names(logFC_vec) <- res$gene_symbol
          
          
          # fix the underflow...
          scores <- scDAPP::fix_underflow(scores, logFC_vec)
          
          #make sure scores is in order of genes...
          ### update May 7 2024 --> there was an error in this before
          # it caused scrambling of genes and incorrect pathway analysis :/
          # scores <- scores[match(res$gene_symbol, res$gene_symbol)]
          scores <- scores[match(res$gene_symbol, names(scores))]
          
          
          #put in res
          res$weight <- scores
          rm(scores, logFC_vec)
          
          
          #order by weight
          res <- res[order(res$weight, decreasing = T),]
          
          
          
          #cbind normcounts
          nc <- cpm(eobj)
          nc <- nc[match(rownames(res), rownames(nc)),]
          
          #add norm counts and raw counts
          rc <- gem
          rc <- rc[match(rownames(res), rownames(rc)),]
          
          
          #cbind norm counts and rawcounts
          colnames(nc) <- paste0('normcounts_', colnames(nc))
          colnames(rc) <- paste0('rawcounts_', colnames(rc))
          res <- cbind(res, nc, rc)
          
          
          
        }
        
        
        ### DESeq2
        
        if(DE_test == 'DESeq2' | DE_test == 'DESeq2-LRT'){
          
          require(DESeq2)
          
          
          #create dds obj
          dds <- DESeqDataSetFromMatrix(gem, comp_pseudobulk_md,
                                        design = ~ Condition)
          
          
          if(DE_test == 'DESeq2'){
            
            #run DESeq2 
            dds <- DESeq(dds)
            
          }
          
          
          if(DE_test == 'DESeq2-LRT'){
            
            #run DESeq2 
            dds <- DESeq(dds, test = 'LRT', reduced = ~1)
            
            
          }
          
          
          
          #get res
          res <- as.data.frame(results(dds))
          
          #set NA padj values to 1 as per guide
          # https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#i-want-to-benchmark-deseq2-comparing-to-other-de-tools.
          res[is.na(res$padj), "padj"] <- 1
          
          
          
          #add gene_symbol
          # attach gene name as a column
          res <- cbind(rownames(res), res)
          colnames(res)[1] <- 'gene_symbol'
          
          
          
          #cbind pct.1 and pct.2
          cellsexp <- cellsexp[match(rownames(res), cellsexp$gene),]
          res <- cbind(res, cellsexp[,-1])
          
          
          
          ## prep weighted list ##
          
          #we have to deal with underflow...
          # get -log10 pvalues, sort
          scores <- -log10(res$pvalue)
          scores <- scores * sign(res$log2FoldChange)
          names(scores) <- res$gene_symbol
          
          #sort by log pval with names
          scores <- sort(scores,decreasing = T)
          
          
          #also get logFC vector; for INf, we will sort them by LFC...
          logFC_vec <- res$log2FoldChange; names(logFC_vec) <- res$gene_symbol
          
          
          # fix the underflow...
          scores <- scDAPP::fix_underflow(scores, logFC_vec)
          
          #make sure scores is in order of genes...
          ### update May 7 2024 --> there was an error in this before
          # it caused scrambling of genes and incorrect pathway analysis :/
          # scores <- scores[match(res$gene_symbol, res$gene_symbol)]
          scores <- scores[match(res$gene_symbol, names(scores))]
          
          
          #put in res
          res$weight <- scores
          rm(scores, logFC_vec)
          
          
          #add norm counts and raw counts
          nc <- counts(dds, normalized = T)
          rc <- counts(dds, normalized = F)
          
          
          #cbind norm counts and rawcounts
          colnames(nc) <- paste0('normcounts_', colnames(nc))
          colnames(rc) <- paste0('rawcounts_', colnames(rc))
          res <- cbind(res, nc, rc)
          
          
          
          #order by weight
          res <- res[order(res$weight, decreasing = T),]
          
          
          
        }
        
        
        
        
        
        
        
        
        return(res)
        
        
      })
      
      #name them by cluster
      # use the nicelabs defined above
      names( m_bycluster_crosscondition_de ) <- groupinglev_nicelabs
      
      
      #remove empty clusters
      m_bycluster_crosscondition_de <- m_bycluster_crosscondition_de[lengths(m_bycluster_crosscondition_de) > 0]
      
      
      
      return(m_bycluster_crosscondition_de)
      
    })
    
    
    
    names(m_bycluster_crosscondition_de_comps) <- comps$labels
    
  }
  
  
  
  
  
  ### DE: pseudobulk or single-cell
  
  if(Pseudobulk_mode == F){
    
    ### loop thru comparisons ###
    # DE in each comparison
    # make sure to select clusters shared by the two conditons
    # for those clusters, use all samples for edgeR, and use c1 vs c2 for contrast functon
    
    
    
    
    compslen <- 1:nrow(comps)
    m_bycluster_crosscondition_de_comps <- lapply(compslen, function(compidx){
      
      
      
      
      
      #get comparison condition levels
      c1 <- comps[compidx,1]
      c2 <- comps[compidx,2]
      
      
      #get lab
      lab <- comps[compidx, 3]
      
      message('\n', lab)
      
      
      ## get only this comp samples
      #subset MD
      comp_pseudobulk_md <- sample_metadata[sample_metadata$Condition %in% c(c1,c2),]
      
      #adjust levels, for wilcox c1 = first, c2 = second
      comp_pseudobulk_md$Condition <- factor(comp_pseudobulk_md$Condition, levels = c(c1, c2))
      
      
      
      
      
      
      #### for each int cluster, loop thru and compare condition 1 vs condition 2 ####
      
      ## loop thru SHARED clusters ##
      clusters <- groupinglevs
      names(clusters) <- clusters
      
      clust = clusters[1] #for testing
      
      
      
      m_bycluster_crosscondition_de <- lapply(clusters, function(clust){
        
        
        message(clust)
        
        
        #subset for each cluster
        bigmd <- sobjint@meta.data
        bigmd <- bigmd[bigmd$Condition %in% c(c1,c2),]
        clustmd <- bigmd[bigmd[,grouping_variable] == clust,]
        
        #check if minimum num cells
        clustmd$Condition <- factor(clustmd$Condition, levels = c(c1, c2))
        cellnums <- table(clustmd$Condition)
        
        if( (cellnums[c1] < 5 | cellnums[c2] < 5) ){
          return()
        }
        
        
        
        #if all good then subset and run
        sobjsub <- sobjint[,rownames(clustmd)]
        
        #rerun sct adjustment?
        # sobjsub <- PrepSCTFindMarkers(sobjsub)
        
        
        #do DE
        
        #turn off parallelization, this step caused memory leak even on tiny datasets
        #future::plan('multisession', workers=workernum)
        
        res <- FindMarkers(sobjsub, logfc.threshold = 0, min.pct = 0,
                           ident.1 = c1, ident.2 = c2,
                           assay = assay, slot = slot,
                           group.by = 'Condition',
                           test.use = DE_test)
        
        # future::plan(strategy = 'sequential')
        
        
        # rm(sobjsub)
        
        #reformat table
        #gene symbol
        res <- cbind(rownames(res), res)
        colnames(res)[1] <- 'gene_symbol'
        rownames(res) <- NULL
        
        #add pct.diff
        res$pct.diff <- res$pct.1 - res$pct.2
        
        
        # #add weight:
        # # -log10 pvalue * sign LFC * abs value of percent difference
        # # ie, significe of DE * sign of DE * num cells exp gene in that direction
        # ## downweight mismatch sign vs pct diff genes... do this by squring them, which makes decimal smaller ##
        # # add + 1 to abs pct diff --> do this to keep FGSEA scores high, otherwise we are actually dividing them by pct diff
        # pctdiff <- res$pct.diff
        # pctdiff[sign(pctdiff) != sign(res$avg_log2FC)] <- pctdiff[sign(pctdiff) != sign(res$avg_log2FC)] ^ 2
        # pctdiff <- abs(pctdiff) + 1
        # res$weight <- -log10(res$p_val) * res$avg_log2FC * pctdiff
        
        ## UPDATE DEC 7 2023, WEIGHT BY -LOG10(PVAL) * SIGN OF LFC
        
        ## prep weighted list ##
        
        #we have to deal with underflow...
        # get -log10 pvalues, sort
        scores <- -log10(res$p_val)
        scores <- scores * sign(res$avg_log2FC)
        names(scores) <- res$gene_symbol
        
        #sort by log pval with names
        scores <- sort(scores,decreasing = T)
        
        
        #also get logFC vector; for INf, we will sort them by LFC...
        logFC_vec <- res$avg_log2FC; names(logFC_vec) <- res$gene_symbol
        
        
        # fix the underflow...
        scores <- scDAPP::fix_underflow(scores, logFC_vec)
        
        
        #make sure scores is in order of genes...
        ### update May 7 2024 --> there was an error in this before
        # it caused scrambling of genes and incorrect pathway analysis :/
        # scores <- scores[match(res$gene_symbol, res$gene_symbol)]
        scores <- scores[match(res$gene_symbol, names(scores))]
        
        #put in res
        res$weight <- scores
        rm(scores, logFC_vec)
        
        #order by weight
        res <- res[order(res$weight, decreasing = T),]
        
        res
        
        
      })
      
      #name them by cluster
      # use the nicelabs defined above
      names( m_bycluster_crosscondition_de ) <- groupinglev_nicelabs
      
      
      #remove empty clusters
      m_bycluster_crosscondition_de <- m_bycluster_crosscondition_de[lengths(m_bycluster_crosscondition_de) > 0]
      
      
      
      
      return(m_bycluster_crosscondition_de)
      
    })
    
    
    
    names(m_bycluster_crosscondition_de_comps) <- comps$labels
    
    
    
  }
  
  
  
  
  
  ## write out the DE results for each comparison, cross condition for each cluster
  # we will then reformat the results, and then last we will write out a table of num DEGs
  compslen <- 1:nrow(comps)
  
  invisible(
    lapply(compslen, function(compidx){
      
      
      #get comparison condition levels
      c1 <- comps[compidx,1]
      c2 <- comps[compidx,2]
      
      #get comp lab
      lab <- comps[compidx,3]
      
      #get cross conditions res per cluster list
      m_bycluster_crosscondition_de <- m_bycluster_crosscondition_de_comps[[compidx]]
      
      
      
      de_cross_conditions_dir <- paste0(outdir_int, '/differentialexpression_crosscondition/',c1,'_vs_', c2, '/')
      
      dir.create(de_cross_conditions_dir, recursive = T)
      
      
      #for each cluster, write out to csv files
      invisible(
        lapply( 1:length(m_bycluster_crosscondition_de), function(i){
          
          #get cluster name
          clustername <- names(m_bycluster_crosscondition_de)[i]
          
          #get cluster DE result
          m <- m_bycluster_crosscondition_de[[i]]
          
          #save file
          de_cross_conditions_file <- paste0(de_cross_conditions_dir, '/', clustername, '.csv')
          write.csv(m, de_cross_conditions_file, quote = F, row.names = F)
          
        })
      )
      
      
      
      # ### write out num DEGs summary table ###
      # #prep NUMDEGS object using THRESHOLD OBJECTS
      # numdegs <- sapply(m_bycluster_crosscondition_de, function(m){
      #   
      #   #normal fdr and padj thresholds
      #   m <- m[m$FDR < crossconditionDE_padj_thres,, drop=F]
      #   m <- m[abs(m$logFC) > crossconditionDE_lfc_thres,, drop=F]
      #   
      #   #pct thresholds: +FC, pct1 > 0.1; -FC, pct2 > 0.1
      #   
      #   upm <- m[m$logFC > 0,,drop=F]
      #   upm <- upm[upm$pct.1 > crossconditionDE_min.pct,, drop=F]
      #   
      #   dnm <- m[m$logFC < 0,,drop=F]
      #   dnm <- dnm[dnm$pct.2 > crossconditionDE_min.pct,, drop=F]
      #   
      #   m <- rbind(upm,dnm)
      #   
      #   try( table( factor(sign(m$logFC), levels=c(-1,1)) ) )
      # })
      # 
      # numdegs <- t(numdegs)
      # colnames(numdegs) <- c(c2, c1)
      # 
      # #make sure all clusters are shown
      # # make a fake df and replace fake with real res
      # numdegs_all <- data.frame(Cluster = groupinglev_nicelabs,
      #                           c1 = 0, c2 = 0)
      # colnames(numdegs_all) <- c('Cluster', c1, c2)
      # rownames(numdegs_all) <- numdegs_all$Cluster
      # 
      # numdegs_all[rownames(numdegs), c1] <- numdegs[,c1]
      # numdegs_all[rownames(numdegs), c2] <- numdegs[,c2]
      # rownames(numdegs_all) <- NULL
      # 
      # 
      # write.csv(numdegs_all, paste0(outdir_int, '/differentialexpression_crosscondition/',c1,'_vs_', c2, '_numDEGs_summary.csv'), quote = F, row.names = T)
      # 
      
      
    })
  )
  
  
  
  
  #if wilcoxon is used, reformat the dataframe to match edgeR, to ease the pathway analysis
  
  if( Pseudobulk_mode == F | DE_test == 'DESeq2' | DE_test == 'DESeq2-LRT' ){
    
    
    
    #for each comparison, loop thru each cluster, and reformat the table
    
    ##  REFORMAT FOR SINGLE CELL WILCOX TEST OUTPUT
    if(Pseudobulk_mode == F){
      
      
      
      m_bycluster_crosscondition_de_comps <- lapply(m_bycluster_crosscondition_de_comps, function(m_bycluster_crosscondition_de){
        
        
        m_bycluster_crosscondition_de <- lapply(m_bycluster_crosscondition_de, function(res){
          
          #reformat all
          res <- res[,c("gene_symbol", "avg_log2FC","p_val", "p_val_adj", "pct.1", "pct.2", "pct.diff", "weight" )]
          
          #rename to match edgeR
          colnames(res) <- c("gene_symbol", "logFC","PValue", "qvalue_bonferroni", "pct.1", "pct.2","pct.diff", "weight"  )
          
          res$FDR <- p.adjust(res$PValue, method = 'fdr')
          
          res
          
        })
        
      })
      
    }
    
    
    
    ##  REFORMAT FOR PSEUDOBULK DESEQ2 / DESEQ2-LRT CELL WILCOX TEST OUTPUT
    if( Pseudobulk_mode == T & (DE_test == 'DESeq2' | DE_test == 'DESeq2-LRT') ){
      
      
      m_bycluster_crosscondition_de_comps <- lapply(m_bycluster_crosscondition_de_comps, function(m_bycluster_crosscondition_de){
        
        
        m_bycluster_crosscondition_de <- lapply(m_bycluster_crosscondition_de, function(res){
          
          #select columns
          # we need to select the columns with raw and norm counts
          counts_colnames <- colnames(res)[(grep(pattern = 'weight', x = colnames(res)) + 1):ncol(res)]
          
          #subset cols
          # res <- res[,c("gene_symbol", "avg_log2FC","p_val", "p_val_adj", "pct.1", "pct.2", "pct.diff", "weight" , counts_colnames)]
          res <- res[,c("gene_symbol", "log2FoldChange","pvalue", "padj", "pct.1", "pct.2", "pct.diff", "weight" , counts_colnames)]
          
          #rename to match edgeR
          colnames(res) <- c("gene_symbol", "logFC","PValue", "qvalue_bonferroni", "pct.1", "pct.2","pct.diff", "weight" , counts_colnames  )
          
          res$FDR <- p.adjust(res$PValue, method = 'fdr')
          
          res
          
        })
        
      })
      
    }
    
    
  }
  
  
  
  
  
  ## write a table with num DEGs
  compslen <- 1:nrow(comps)
  invisible(
    lapply(compslen, function(compidx){
      
      
      #get comparison condition levels
      c1 <- comps[compidx,1]
      c2 <- comps[compidx,2]
      
      #get comp lab
      lab <- comps[compidx,3]
      
      #get cross conditions res per cluster list
      m_bycluster_crosscondition_de <- m_bycluster_crosscondition_de_comps[[compidx]]
      
      
      ### write out num DEGs summary table ###
      #prep NUMDEGS object using THRESHOLD OBJECTS
      numdegs <- sapply(m_bycluster_crosscondition_de, function(m){
        
        #normal fdr and padj thresholds
        m <- m[m$FDR < crossconditionDE_padj_thres,, drop=F]
        m <- m[abs(m$logFC) > crossconditionDE_lfc_thres,, drop=F]
        
        #pct thresholds: +FC, pct1 > 0.1; -FC, pct2 > 0.1
        
        upm <- m[m$logFC > 0,,drop=F]
        upm <- upm[upm$pct.1 > crossconditionDE_min.pct,, drop=F]
        
        dnm <- m[m$logFC < 0,,drop=F]
        dnm <- dnm[dnm$pct.2 > crossconditionDE_min.pct,, drop=F]
        
        m <- rbind(upm,dnm)
        
        try( table( factor(sign(m$logFC), levels=c(-1,1)) ) )
      })
      
      numdegs <- t(numdegs)
      colnames(numdegs) <- c(c2, c1)
      
      #make sure all clusters are shown
      # make a fake df and replace fake with real res
      numdegs_all <- data.frame(Cluster = groupinglev_nicelabs,
                                c1 = 0, c2 = 0)
      colnames(numdegs_all) <- c('Cluster', c1, c2)
      rownames(numdegs_all) <- numdegs_all$Cluster
      
      numdegs_all[rownames(numdegs), c1] <- numdegs[,c1]
      numdegs_all[rownames(numdegs), c2] <- numdegs[,c2]
      rownames(numdegs_all) <- NULL
      
      
      write.csv(numdegs_all, paste0(outdir_int, '/differentialexpression_crosscondition/',c1,'_vs_', c2, '_numDEGs_summary.csv'), quote = F, row.names = T)
      
      
      
    })
  )
  
  
  
  
  #save RDS file of DE res list -- > we save this later anyway, just wastes memory
  # deres_rds_file <- paste0(outdir_int, '/differentialexpression_crosscondition/m_bycluster_crosscondition_de_comps.rds')
  # saveRDS(m_bycluster_crosscondition_de_comps, deres_rds_file)
  
  
  return(m_bycluster_crosscondition_de_comps)
  
  
}









#' Prep MSIGDB pathways for pathway analysis
#'
#' This is a modular component of the scRNAseq analysis pipeline. Prep pathways from MSIGDB via the msigdbr package. We include the Hallmarks category; Gene Ontology BP, MF and CC; Reactome; KEGG; transcription factor CHIP-seq targets in the Gene Transcription Regulation Database (TFT_GTRD); inferred transcription factor targets via motif analysis from Xie et al Nature 2005 (TFT_Legacy). Only gene sets with < 500 genes are included.
#'
#' @param species string, species such as "Homo sapeins" or "Mus musculus"
#' @param outdir_int string, directory to save pathways to. Will create a sub-directory called "pathwayanalysis_crosscondition" and save inside of there. We save pathways since the database updates over time.
#'
#' @return a data.frame similar to the output of `msigdbr::msigdbr`, but filtering for some specific categories / subcategories.
#' @export
#'
#' @examples
#' \dontrun{
#' pathways <- preppathways_pathwayanalysis_crosscondition_module(
#' species = 'Mus musculus',
#' outdir_int = 'path/to/directory')
#' }
preppathways_pathwayanalysis_crosscondition_module <- function(species,
                                                               outdir_int)
{


  require(msigdbr)


  if( !missing(outdir_int) ) {

    #prep the pathways
    # make sure to save it. database can update over time
    pwayoutdir <- paste0(outdir_int, '/pathwayanalysis_crosscondition/')
    dir.create(pwayoutdir, recursive = T)
    
    
    if( !file.exists( paste0(pwayoutdir, '/msigdb_pathways.rds') ) ){

      message('Accessing MSIGDBR database')
      pathways <- msigdbr::msigdbr(species = species)

      #replace : with _ in actual pathway names:
      pathways$gs_subcat <- gsub(':', '_', pathways$gs_subcat)

      saveRDS(pathways, paste0(pwayoutdir, '/msigdb_pathways.rds') )

    } else{
      message('Reading cached msigdbr pathways')
      pathways <- readRDS(paste0(pwayoutdir, '/msigdb_pathways.rds') )
    }


  }
  else{

    message('Accessing MSIGDBR database')
    #read pathways
    pathways <- msigdbr::msigdbr(species = species)

    #replace : with _ in actual pathway names:
    pathways$gs_subcat <- gsub(':', '_', pathways$gs_subcat)
  }


  #### picking default categories
  # because hallmark is a "category" and rest are "subcategories", it is hard to make this automated
  # guess it may be possible if we set missing subcat as cat...
  # table( pathways[pathways$gs_subcat=='',"gs_cat"] )
  # for now hardcode these
  pwaycats <- c("HALLMARK", "GO_BP", "GO_MF", "GO_CC", "CP_REACTOME", "CP_KEGG", "TFT_GTRD", "TFT_TFT_Legacy")


  ### try to replace ':' with "_"
  # in pwaycats, user provided subcategories:
  pwaycats <- gsub(':', '_', pwaycats)
  names(pwaycats) <- pwaycats

  #in actual pathway names:
  pathways$gs_subcat <- gsub(':', '_', pathways$gs_subcat)

  #also get hallmarks...
  pathways[pathways$gs_cat == 'H', 'gs_subcat'] <- "HALLMARK"

  #prep pathways using categories defined by user
  pathways <- as.data.frame( pathways[pathways$gs_subcat %in% pwaycats,] )

  #fgsea recommends no pathways over 500 genes
  pathways <- pathways[table(pathways$gs_name) <= 500,]

  #let's also remove pathways with less than 3 genes
  pathways <- pathways[table(pathways$gs_name) >= 3,]


  #purge mem
  invisible(gc(full = T, reset = F, verbose = F))


  return(pathways)


}





#' GSEA analysis for cross condition DE analysis
#'
#' This function is a modular component of the scRNAseq pipeline. Perform GSEA analysis via the FGSEA package on the results of differential expression (DE) analysis for cross-condition comparison. Multiple conditions are supported. Plots and tables are saved.
#'
#' @param m_bycluster_crosscondition_de_comps the output of `scDAPP::de_across_conditions_module()`.
#' @param pathways data.frame, the output of `scDAPP::preppathways_pathwayanalysis_crosscondition_module()`
#' @param sample_metadata data.frame with sample names and conditions, same as in `scDAPP::de_across_conditions_module()`, see that function's documentation for description.
#' @param comps data.frame with conditions to test in GSEA, same as in `scDAPP::de_across_conditions_module()`, see that function for description
#' @param pathway_padj_thres numeric, threshold for significance of pathway enrichment after multiple test correction
#' @param pwaycats UNTESTED CURRENTLY. character vector of msigdb pathways data.frame in gs_subcat to run.
#' @param workernum integer. number of CPUs. default = 1.
#' @param outdir_int string, directory to save pathways to. Will create a sub-directory called "pathwayanalysis_crosscondition" and save inside of there.
#' @param cp.font.size numeric. size of pathway name test in summary plots. larger than 5 will likely result in name overlap and unreadable plots, currently not easy to solve
#'
#' @return a list with two elements. First contains the raw results and plots as a nested loop: first level is A vs B comparison, then category-by-category of MSIGDB database, then cluster-by-cluster. The second list element is similar but just has a summary plot showing top pathways across clusters.
#' @export
#'
#' @examples
#' \dontrun{
#'
#'
#' # FIRST: run `scDAPP::de_across_conditions_module()`. the output object of that is used as the main input for this pathway analysis function.
#'
#' # SECOND: prep pathways before running
#' pathways <- scDAPP::preppathways_pathwayanalysis_crosscondition_module(species = species,
#' outdir_int = outdir_int)
#'
#' # THIRD: run the pathway analysis.
#' # see `scDAPP::de_across_conditions_module()` for a description of the sample_metadata and comps files.
#' pways_output_list <- pathwayanalysis_crosscondition_module(
#' m_bycluster_crosscondition_de_comps = m_bycluster_crosscondition_de_comps,
#' pathways = pathways,
#' sample_metadata = sample_metadata,
#' comps = comps,
#' workernum = 6,
#' outdir_int = outdir_int
#' )
#'
#' # Access the output files
#' # note that all results will be saved to outdir_int.
#' pathway_analysis_mainlist_comps <- pways_output_list$pathway_analysis_mainlist_comps
#' pathwaysummplots_comps <- pways_output_list$pathwaysummplots_comps
#'
#' #these are both lists that contain the raw tables and per-cluster plots (first object);
#' # and summaryplots that show pathway enrichment across clusters (second object).
#'
#' }
pathwayanalysis_crosscondition_module <- function(m_bycluster_crosscondition_de_comps,
                                                  pathways,
                                                  sample_metadata,
                                                  comps,
                                                  pathway_padj_thres,
                                                  pwaycats,
                                                  workernum,
                                                  outdir_int,
                                                  cp.font.size
){
  
  
  require(fgsea)
  require(foreach)
  require(doParallel)
  require(parallel)
  
  
  #   UPDATE DECEMBER 7 2023 deg.weight has been deprecated, we will stick with -log10pval * sign FC, but weight value can be modified for res before this
  #   if( missing(deg.weight) ){deg.weight <- 'pval'}
  #   if( !(deg.weight %in% c('auto', 'pval') ) ){
  #     stop('deg.weight must be either "auto" or "pval" ')
  #   }
  
  if( missing(cp.font.size) ) {
    
    #set font size
    # this is the best size, any bigger there will be overlap...
    cp.font.size <- 5
  }
  
  
  
  
  #### picking default categories
  # because hallmark is a "category" and rest are "subcategories", it is hard to make this automated
  # guess it may be possible if we set missing subcat as cat...
  # table( pathways[pathways$gs_subcat=='',"gs_cat"] )
  # for now hardcode these
  if(missing(pwaycats)){
    pwaycats <- c("HALLMARK", "GO_BP", "GO_MF", "GO_CC", "CP_REACTOME", "CP_KEGG", "TFT_GTRD", "TFT_TFT_Legacy")
  }
  
  if(missing(pathway_padj_thres)){
    pathway_padj_thres <- 0.1
  }
  if(missing(workernum)){
    workernum <- 1
  }
  
  
  ### pathway analysis for each condition comparison
  
  # for each condition comparison,
  # for each cluster
  # do pathway analysis and make plot
  
  
  
  
  
  #prep names
  comps$labels <- paste0(comps$c1, '_vs_', comps$c2)
  
  compslen <- 1:nrow(comps)
  compidx = 1 #for testing
  
  
  
  pathway_analysis_mainlist_comps <- lapply(compslen, function(compidx){
    
    
    #get comparison condition levels
    c1 <- comps[compidx,1]
    c2 <- comps[compidx,2]
    
    #get comp lab
    lab <- comps[compidx,3]
    
    message(lab)
    
    #get cross conditions res per cluster list
    m_bycluster_crosscondition_de <- m_bycluster_crosscondition_de_comps[[compidx]]
    
    
    
    pwayoutdir <- paste0(outdir_int, '/pathwayanalysis_crosscondition/',c1,'_vs_', c2, '/')
    if( !dir.exists(pwayoutdir) ){ dir.create(pwayoutdir, recursive = T) }
    
    
    ### loop thru pathway categories
    names(pwaycats) <- pwaycats
    
    
    #get clust / grouping names
    clusters <- names(m_bycluster_crosscondition_de)
    
    
    
    
    #set gene universe
    pwaycat <- pwaycats[1] #for testing
    
    pathway_analysis_mainlist <- lapply(pwaycats, function(pwaycat){
      
      message('\n\n', pwaycat, '\n\n')
      
      
      #get pways and genes in this category
      term2gene <- pathways[pathways$gs_subcat == pwaycat,c('gs_name', 'gene_symbol')]
      
      
      #pways as list for gsea
      pwayl = split(term2gene$gene_symbol, term2gene$gs_name)
      rm(term2gene)
      
      
      #get list of pathways upreg in each cluster
      
      cl <- parallel::makeCluster(workernum, rscript_args = c("--no-init-file", "--no-site-file", "--no-environ"))
      doParallel::registerDoParallel(cl)
      
      
      
      #pwayres_DE_across_conditions_per_cluster <- lapply(clusters, function(clust){
      pwayres_DE_across_conditions_per_cluster <- foreach(clust = clusters,
                                                          .packages = c('fgsea', 'ggplot2'),
                                                          .export = c( 'm_bycluster_crosscondition_de', 'pathway_padj_thres', 'cp.font.size'),
                                                          .noexport = c('pathways'),
                                                          .verbose = T) %dopar%
        {
          
          
          
          invisible(gc(full = T, reset = F, verbose = F))
          
          
          #get DEG res
          res <- m_bycluster_crosscondition_de[[clust]]
          
          
          
          
          # UPDATE DEC 7 2023 WE DEPRECATED DEG.WEIGHT
          # AND MADE DEFAULT WEIGHT FROM DE MODULE AS -LOG10(PVAL ) * SIGN(L2FC)
          # if(deg.weight == 'pval'){
          #
          #   ## prep weighted list ##
          #
          #   #we have to deal with underflow...
          #   # get -log10 pvalues, sort
          #   scores <- -log10(res$PValue)
          #   scores <- scores * sign(res$logFC)
          #   names(scores) <- res$gene_symbol
          #
          #   #sort by log pval with names
          #   scores <- sort(scores,decreasing = T)
          #
          #
          #   #also get logFC vector; for INf, we will sort them by LFC...
          #   logFC_vec <- res$logFC; names(logFC_vec) <- res$gene_symbol
          #
          #
          #   # fix the underflow...
          #   scores <- fix_underflow(scores, logFC_vec)
          #
          #
          #   #put in res
          #   res$weight <- scores
          #   rm(scores, logFC_vec)
          #
          #
          # }
          
          #input for gsea is weight named by gene
          res <- res[order(res$weight, decreasing = T),]
          gl <- res$weight; names(gl) <- res$gene_symbol
          
          #clean env
          rm(res)
          
          ## run GSEA ##
          # first run multilevel, then try npermsimple = 1000
          gseares <- fgsea::fgsea(pathways=pwayl, stats=gl, nproc = 1)
          
          invisible(gc(full = T, reset = F, verbose = F))
          
          
          #sometimes there are NAs due to "severely unbalanced pathways", try to fix
          if( any(is.na(gseares$NES)) ){
            rm(gseares)
            
            gseares <- fgsea::fgsea(pathways=pwayl, stats=gl, nPermSimple=10000, nproc = 1)
            
            invisible(gc(full = T, reset = F, verbose = F))
            
          }
          
          #clean env
          rm(gl)
          
          
          #format as data.frame instead of data.table
          gseares <- as.data.frame(gseares)
          
          #ensure no NAs are kept
          gseares <- gseares[complete.cases(gseares[,1:7]),,drop=F]
          
          #order by NES
          gseares <- gseares[order(gseares$NES, decreasing = T),]
          
          #apply cutoff of pathway_padj_thres
          gseares <- gseares[gseares$padj < pathway_padj_thres, ,drop=F]
          
          #select pathways with more than just 1 gene in the list
          gseares <- gseares[gseares$size > 2,,drop=F]
          
          #skip if no significant results
          if( nrow(gseares)==0){ return() }
          
          
          #prep for plot, leave out leading edge
          gseares_plot <- gseares[,-8]
          
          #if more than 20 ,select just 20
          gseares_plot <- rbind( head( gseares_plot[gseares_plot$NES>0,,drop=F], 10) ,
                                 tail( gseares_plot[gseares_plot$NES<0,,drop=F], 10) )
          
          #make pathway names more readable by using spaces instead of underscores
          gseares_plot$pathway <- gsub(gseares_plot$pathway, pattern = '_', replacement = ' ')
          
          #make pathway names more readable by splitting long ones to multiple lines
          gseares_plot$pathway <- stringr::str_wrap(gseares_plot$pathway, width = 35)
          
          #make sure order is by -log(padj) * NES
          gseares_plot$weight <- -log(gseares_plot$padj) * sign(gseares_plot$NES)
          gseares_plot <- gseares_plot[order(gseares_plot$weight, decreasing = T),]
          gseares_plot$pathway <- factor(gseares_plot$pathway, levels = rev(gseares_plot$pathway)  )
          
          
          #plot it
          
          #fix color issue when just 1 obs
          if(nrow(gseares_plot) == 1){
            dp_single_col <- sign(gseares_plot$NES)
            dp_single_col <- ifelse(dp_single_col==1, yes = 'red', no = 'steelblue')
            
            
            dp <- ggplot(gseares_plot, aes(-log10(padj), pathway, col=NES, size = size))+
              geom_point()+
              theme_linedraw()+
              theme(axis.text=element_text(size=cp.font.size) )+
              scale_color_gradientn(colors = dp_single_col) +
              scale_size(range=c(2,6))
            
          } else{
            
            
            dp <- ggplot(gseares_plot, aes(-log10(padj), pathway, col=NES, size = size))+
              geom_point()+
              theme_linedraw()+
              theme(axis.text=element_text(size=cp.font.size) )+
              scale_color_gradient2(low = 'steelblue', high = 'red', mid = 'white', midpoint = 0, name = 'Normalized\nEnrichment\nScore')+
              scale_size(range=c(2,6))
            
          }
          
          #return the result table and the plot
          return(list(gseares = gseares, dp = dp))
          
          
          
          
        } #per-cluster loop for this pathway category loop end
      
      parallel::stopCluster(cl)
      
      
      names(pwayres_DE_across_conditions_per_cluster) <- clusters
      
      
      return(pwayres_DE_across_conditions_per_cluster)
      
    }) #per-pathway category loop end
    
    
    
    
    ### remove all NULLS (clusters with no pathways)
    
    # remove null categories, ie entire category had no significant pathways
    
    #recursively set all missing to 0
    pathway_analysis_mainlist = lapply(pathway_analysis_mainlist, function(pwayres_cats){
      
      pwayres_cats <- lapply(pwayres_cats, function(pwayres_clusts){
        pwayres_clusts[lengths(pwayres_clusts) > 0]
      })
      
      pwayres_cats[lengths(pwayres_cats) > 0]
      
      
    })
    
    #remove any missing categories
    pathway_analysis_mainlist <- pathway_analysis_mainlist[lengths(pathway_analysis_mainlist)>0]
    
    
    invisible(gc(full = T, reset = F, verbose = F))
    
    
    #save pway analysis for this comparison
    
    # loop over this new pwayruns, since some categories theoretically don't have any enriched though unlikely
    pwayruns <- names(pathway_analysis_mainlist)
    
    
    #for each category, get clster res in that cateogry,
    # for each cluster, save the up/down csv and plots
    invisible(
      finalpwayouts <- lapply(pwayruns, function(pwaycat){
        
        
        
        # message(pwaycat)
        
        
        subcatout <- paste0(pwayoutdir, '/', pwaycat, '/')
        
        # dir.create(subcatout) --> do this with recursive later, maybe prevent even making it if all don't work
        
        clustres <- pathway_analysis_mainlist[[pwaycat]]
        
        clusters <- names(clustres)
        
        
        #for each cluster, get up/down csv, up/dwon plot, and save
        numpways <- lapply(clusters, function(clust){
          
          
          
          
          pwayres_DE_across_conditions_per_cluster <- clustres[[clust]]
          
          
          
          if(is.null(pwayres_DE_across_conditions_per_cluster)){return()}
          
          gseares <- pwayres_DE_across_conditions_per_cluster$gseares
          dp <- pwayres_DE_across_conditions_per_cluster$dp
          
          #save PDFs and CSVs
          
          subcatout_clustdir <- paste0(subcatout, '/', clust, '/')
          
          suppressWarnings(dir.create(subcatout_clustdir, recursive = T))
          
          
          #upcsv
          subcatout_clustdir_gseares <- paste0(subcatout_clustdir, '/pathwaytable.csv')
          
          #gseares, leading edge needs to be adjusted...
          gseares$leadingEdge <- sapply(gseares$leadingEdge, function(x){ paste(x, collapse = '/') })
          
          write.csv(gseares, subcatout_clustdir_gseares, quote = F, row.names = F)
          
          
          subcatout_clustdir_dp <- paste0(subcatout_clustdir, '/dotplot_toppathways.pdf')
          
          pdf(subcatout_clustdir_dp)
          print( dp )
          
          while (!is.null(dev.list()))  dev.off()
          
          
          
          nrow(gseares)
          
          
          
          
          
          
          
        } ) # close clusters lapply
        
        
      }) # close saving loop for all categories
      
    ) # close invisible wrap around lapply
    
    
    
    
    
    
    
    
    
    
    ### close any open devices
    while (!is.null(dev.list()))  dev.off()
    
    
    
    
    return(pathway_analysis_mainlist)
    
    
  }) # close cross condition lapply
  
  
  names(pathway_analysis_mainlist_comps) <- comps$labels
  
  
  
  
  #remove big objects
  rm(pathways)
  invisible(gc(full = T, reset = F, verbose = F))
  
  
  
  
  ### prep summary plots for each category
  
  compslen <- 1:nrow(comps)
  pathwaysummplots_comps <- lapply(compslen, function(compidx){
    
    #get pway analysis
    pathway_analysis_mainlist <- pathway_analysis_mainlist_comps[[compidx]]
    
    #get comparison condition levels
    c1 <- comps[compidx,1]
    c2 <- comps[compidx,2]
    
    #get comp lab
    lab <- comps[compidx,3]
    
    
    ### extract the table from all categories
    
    cat_cpres_list <- lapply(pathway_analysis_mainlist, function(pwaycatlist){
      
      #in each cluster:
      # get the tables from each up/dn
      
      # use cluster index, we need the cluster name
      
      clust_cpres <- lapply(1:length(pwaycatlist), function(clustidx){
        
        clustname <- names(pwaycatlist)[clustidx]
        pwayres_DE_across_conditions_per_cluster <- pwaycatlist[[clustidx]]
        
        #use the dotplot data for the table
        gseares_plot <- pwayres_DE_across_conditions_per_cluster$dp$data
        
        gseares_plot$cluster = clustname
        gseares_plot$condition = c1
        gseares_plot[sign(gseares_plot$NES) == -1, "condition"] = c2
        
        return(gseares_plot)
        
        
      })
      
      dplyr::bind_rows(clust_cpres)
      
    })
    
    
    #loop thru each category's result data.frame, splitting c1 and c2, and plotting
    
    summplots_cats <- lapply( 1:length(cat_cpres_list) , function(catdex){
      
      cpres_cat <- cat_cpres_list[[catdex]]
      catname <- names(cat_cpres_list)[catdex]
      
      
      #make plots for c1 and c2 direction
      # some categories have no pathways significant for condition, just return null
      
      summplots_conds <- lapply( c(c1,c2) , function(cond){
        
        #get result tbale for this condition
        cpres_cat_cond <- cpres_cat[cpres_cat$condition == cond,,drop=F]
        
        #if no conditions, it will haev nrow=0 so just return null
        if(nrow(cpres_cat_cond) == 0){ return() }
        
        # subselect categories if more than 30 total, use just top 5 per pathway
        if(nrow(cpres_cat_cond) > 30){
          #select the ones to pick
          cpres_cat_cond_sub <- cpres_cat_cond %>%
            group_by(cluster) %>%
            top_n(n=5, wt = -log10(padj)) %>%
            top_n(n=5, wt = abs(NES)) %>%
            as.data.frame()
          
          #get them, doing it this way allows viewing shared pathways
          cpres_cat_cond <- cpres_cat_cond[cpres_cat_cond$pathway %in% cpres_cat_cond_sub$pathway,]
          
        }
        
        
        #make sure orders are proper
        # for clusters:
        cpres_cat_cond$cluster <- factor(cpres_cat_cond$cluster, levels = unique(cpres_cat_cond$cluster))
        
        
        #for pathways
        cpres_cat_cond$pathway <- factor(cpres_cat_cond$pathway, levels = rev(unique(cpres_cat_cond$pathway)))
        
        
        ggplot(cpres_cat_cond, aes(x=cluster, y=pathway ,size = -log10(padj), col = NES))+
          geom_point()+
          theme_linedraw()+
          theme(axis.text=element_text(size=cp.font.size),
                axis.text.x = element_text(angle = 45, vjust = 1, hjust=1) )+
          scale_color_gradient2(low = 'steelblue', high = 'red', mid = 'white', midpoint = 0, name = 'Normalized\nEnrichment\nScore')+
          scale_size(range=c(2,6), name = '-log10(padj)')+
          xlab('Cluster')+ylab('')+
          ggtitle(catname, subtitle = cond)
        
        
      }) # close cross-condition loop for summary plots
      
      names(summplots_conds) <- c(c1,c2)
      
      summplots_conds
      
      
    })
    
    
    names(summplots_cats) <- names(cat_cpres_list)
    
    
    
    #print them to pdfs...
    
    pwayoutdir <- paste0(outdir_int, '/pathwayanalysis_crosscondition/',c1,'_vs_', c2, '/')
    
    
    summarypdf <- paste0(pwayoutdir, '/SummaryDotPlots.pdf')
    pdf(summarypdf, width = 7, height = 7)
    
    print(summplots_cats)
    
    while (!is.null(dev.list()))  dev.off()
    
    return(summplots_cats)
    
    
    
  }) # close summary plot across conditons loop
  
  
  names(pathwaysummplots_comps) <- comps$labels
  
  
  
  
  ### for easily reproducing plots and etc, save them as R objects...
  pwayoutdir <- paste0(outdir_int, '/pathwayanalysis_crosscondition/')
  DE_pathways_plot_objects_list <- list(comps = comps,
                                        m_bycluster_crosscondition_de_comps = m_bycluster_crosscondition_de_comps,
                                        pathway_analysis_mainlist_comps = pathway_analysis_mainlist_comps,
                                        pathwaysummplots_comps = pathwaysummplots_comps
  )
  
  
  #save object sizes...
  # pwayobjsizedf <- data.frame(obj = names(DE_pathways_plot_objects_list))
  # pwayobjsizedf$size_bytes <- sapply(DE_pathways_plot_objects_list, object.size, simplify = T)
  #
  #
  # pwayoutdir <- paste0(outdir_int, '/pathwayanalysis_crosscondition/')
  #
  # objsizefile <- paste0(pwayoutdir, '/OBJSIZES_DE_pathways_plot_objects_list.csv')
  #
  # write.csv(pwayobjsizedf, objsizefile, quote = F, row.names = F)
  
  
  
  DE_pathways_plot_objects_list_file <- paste0(pwayoutdir, '/DE_pathways_plot_objects_list.rds')
  
  saveRDS(DE_pathways_plot_objects_list, DE_pathways_plot_objects_list_file)
  
  
  
  ### return just pathway outlist and pathway summplots
  pways_output_list <- list(
    pathway_analysis_mainlist_comps = pathway_analysis_mainlist_comps,
    pathwaysummplots_comps = pathwaysummplots_comps
  )
  
  return(pways_output_list)
  
  
}










#' Overrepresentation analysis for cross condition pathway analysis
#'
#' This function is a modular component of the scRNAseq pipeline. Perform OverRepresentation Analysis (ORA) via the ClusterProfiler package on the results of differential expression (DE) analysis for cross-condition comparison. Multiple conditions are supported. ClusterProfiler objects and tables are saved.
#'
#' @param m_bycluster_crosscondition_de_comps the output of `scDAPP::de_across_conditions_module()`.
#' @param pathways data.frame, the output of `scDAPP::preppathways_pathwayanalysis_crosscondition_module()`
#' @param sample_metadata data.frame with sample names and conditions, same as in `scDAPP::de_across_conditions_module()`, see that function's documentation for description.
#' @param comps data.frame with conditions to test in GSEA, same as in `scDAPP::de_across_conditions_module()`, see that function for description
#' @param pathway_padj_thres numeric, threshold for significance of pathway enrichment after multiple test correction, passed to qvalueCutoff in `clusterProfiler::enricher`
#' @param pwaycats UNTESTED CURRENTLY. character vector of msigdb pathways data.frame in gs_subcat to run.
#' @param workernum integer. number of CPUs. default = 1.
#' @param outdir_int string, directory to save pathways to. Will create a sub-directory called "pathwayanalysis_crosscondition_ORA" and save inside of there.
#'
#' @return a list that contains the raw results and plots as a nested loop: first level is A vs B comparison, then category-by-category of MSIGDB database, then a data.frame of cluster-by-cluster results
#' @export
#'
#' @examples
#' \dontrun{
#'
#'
#' # FIRST: run `scDAPP::de_across_conditions_module()`. the output object of that is used as the main input for this pathway analysis function.
#'
#' # SECOND: prep pathways before running
#' pathways <- scDAPP::preppathways_pathwayanalysis_crosscondition_module(species = species,
#' outdir_int = outdir_int)
#'
#' # THIRD: run the pathway analysis.
#' # see `scDAPP::de_across_conditions_module()` for a description of the sample_metadata and comps files.
#' pways_output_list <- scDAPP::ORA_crosscondition_module(
#' m_bycluster_crosscondition_de_comps = m_bycluster_crosscondition_de_comps,
#' pathways = pathways,
#' sample_metadata = sample_metadata,
#' comps = comps,
#' workernum = 6,
#' outdir_int = outdir_int
#' )
#'
#' # Access the output files
#' # note that all results will be saved to outdir_int.
#' pathway_analysis_mainlist_comps <- pways_output_list$pathway_analysis_mainlist_comps
#' pathwaysummplots_comps <- pways_output_list$pathwaysummplots_comps
#'
#' #these are both lists that contain the raw tables and per-cluster plots (first object);
#' # and summaryplots that show pathway enrichment across clusters (second object).
#'
#' }
ORA_crosscondition_module <- function(m_bycluster_crosscondition_de_comps,
                                      pathways,
                                      sample_metadata,
                                      comps,
                                      pathway_padj_thres,
                                      pwaycats,
                                      workernum,
                                      outdir_int
){
  
  
  require(clusterProfiler)
  require(foreach)
  require(doParallel)
  require(parallel)
  
  
  
  
  
  #### picking default categories
  # for now hardcode these
  if(missing(pwaycats)){
    # pwaycats <- c("HALLMARK", "GO_BP", "GO_MF", "GO_CC", "CP_REACTOME", "CP_KEGG", "TFT_GTRD", "TFT_TFT_Legacy")
    pwaycats <- c("HALLMARK", "GO_BP", "GO_MF", "GO_CC", "CP_REACTOME", "CP_KEGG", "TFT_GTRD", "TFT_TFT_Legacy")
  }
  
  if(missing(pathway_padj_thres)){
    pathway_padj_thres <- 0.1
  }
  if(missing(workernum)){
    workernum <- 1
  }
  
  
  ### pathway analysis for each condition comparison
  
  # for each condition comparison,
  # for each cluster
  # do pathway analysis and make plot
  
  
  
  
  
  #prep names
  comps$labels <- paste0(comps$c1, '_vs_', comps$c2)
  
  compslen <- 1:nrow(comps)
  compidx = 1 #for testing
  
  
  
  pathway_analysis_mainlist_comps <- lapply(compslen, function(compidx){
    
    
    #get comparison condition levels
    c1 <- comps[compidx,1]
    c2 <- comps[compidx,2]
    
    #get comp lab
    lab <- comps[compidx,3]
    
    message(lab)
    
    #get cross conditions res per cluster list
    m_bycluster_crosscondition_de <- m_bycluster_crosscondition_de_comps[[compidx]]
    
    
    
    pwayoutdir <- paste0(outdir_int, '/pathwayanalysis_crosscondition_ORA/',c1,'_vs_', c2, '/')
    if( !dir.exists(pwayoutdir) ){ dir.create(pwayoutdir, recursive = T) }
    
    
    ### loop thru pathway categories
    names(pwaycats) <- pwaycats
    
    
    #get clust / grouping names
    clusters <- names(m_bycluster_crosscondition_de)
    
    
    
    
    #set gene universe
    pwaycat <- pwaycats[1] #for testing
    
    pathway_analysis_mainlist <- lapply(pwaycats, function(pwaycat){
      
      message('\n\n', pwaycat, '\n\n')
      
      
      #get pways and genes in this category
      term2gene <- pathways[pathways$gs_subcat == pwaycat,c('gs_name', 'gene_symbol')]
      
      
      #get list of pathways upreg in each cluster
      
      cl <- parallel::makeCluster(workernum, rscript_args = c("--no-init-file", "--no-site-file", "--no-environ"))
      doParallel::registerDoParallel(cl)
      
      
      clust = clusters[1] #for test
      
      #pwayres_DE_across_conditions_per_cluster <- lapply(clusters, function(clust){
      pwayres_DE_across_conditions_per_cluster <- foreach(clust = clusters,
                                                          .packages = c('clusterProfiler'),
                                                          .export = c( 'm_bycluster_crosscondition_de', 'pathway_padj_thres', 'crossconditionDE_padj_thres', 'crossconditionDE_lfc_thres', 'crossconditionDE_min.pct', 'c1', 'c2'),
                                                          .noexport = c('pathways'),
                                                          .verbose = T) %dopar%
        {
          
          
          
          invisible(gc(full = T, reset = F, verbose = F))
          
          
          #get DEG res
          res <- m_bycluster_crosscondition_de[[clust]]
          
          
          
          ## split by lfc sign and run ##
          sign_nums <- c(1, -1)
          ts <- 1 #test
          signres_l <- lapply(sign_nums, function(ts){
            
            #for this sign (ts), subset res and run clusterProfiler
            subres <- res[sign(res$logFC) == ts,,drop=F]
            
            #subset by padj
            subres <- subres[subres$FDR < crossconditionDE_padj_thres,,drop = F]
            
            #subset by lfc thres, can use abs val
            subres <- subres[abs(subres$FDR) > crossconditionDE_lfc_thres,,drop = F]
            
            #subset by min.pct 1 for positive, min pct.2 for negative
            crossconditionDE_min.pct <- 0
            if(ts == 1){
              subres <- subres[subres$pct.1 > crossconditionDE_min.pct,,drop = F] 
            } else{
              subres <- subres[subres$pct.2 > crossconditionDE_min.pct,,drop = F] 
            }
            
            
            
            #important, select min num genes; let's say 7 genes min for good luck
            if(nrow(subres) < 7){
              return()
            }
            
            
            ## if not, we can proceed with clusterProfiler
            genenames <- subres$gene_symbol
            
            
            
            
            ## run cluster profiler; any error, let's just return a null, maybe dangerous
            tryCatch(
              
              expr = {
                
                ora_res <- enricher(genenames,
                                    TERM2GENE = term2gene,
                                    qvalueCutoff = pathway_padj_thres,
                                    pvalueCutoff = 1
                                    
                )
                
              },
              
              error = function(e){
                
                ora_res <- data.frame()
                
              }
              
            )
            
            
            
            #clusterProfiler has its own trycatch, which returns NULL... i guess we'll use ours though
            if(length(ora_res) == 0){ora_res <- data.frame()}
            
            
            #return null if none 
            if(nrow(ora_res) == 0){return()}
            
            #don't use their useless class, just a less useful data.frame
            ora_res <- as.data.frame(ora_res)
            
            
            # add a direction column
            ora_res$Direction <- ifelse(ts == 1, c1, c2)
            
            
            
            return(ora_res)
            
          })
          
          
          
          #remove empty res
          signres_l <- signres_l[lengths(signres_l) > 0]
          
          if(length(signres_l) == 0){return()} # if there are just no pathways, return null
          
          signres <- dplyr::bind_rows(signres_l)
          
          
          #return the result table and the plot
          return(signres)
          
          
          
          
        } #per-cluster loop for this pathway category loop end
      
      parallel::stopCluster(cl)
      
      
      names(pwayres_DE_across_conditions_per_cluster) <- clusters
      
      
      return(pwayres_DE_across_conditions_per_cluster)
      
    }) #per-pathway category loop end
    
    
    
    
    ### remove all NULLS (clusters with no pathways)
    
    # remove null categories, ie entire category had no significant pathways
    
    #recursively set all missing to 0
    pathway_analysis_mainlist = lapply(pathway_analysis_mainlist, function(pwayres_cats){
      
      pwayres_cats <- lapply(pwayres_cats, function(pwayres_clusts){
        pwayres_clusts[lengths(pwayres_clusts) > 0]
      })
      
      pwayres_cats[lengths(pwayres_cats) > 0]
      
      
    })
    
    #remove any missing categories
    pathway_analysis_mainlist <- pathway_analysis_mainlist[lengths(pathway_analysis_mainlist)>0]
    
    
    invisible(gc(full = T, reset = F, verbose = F))
    
    
    #save pway analysis for this comparison
    
    # loop over this new pwayruns, since some categories theoretically don't have any enriched though unlikely
    pwayruns <- names(pathway_analysis_mainlist)
    
    
    #for each category, get clster res in that cateogry,
    # for each cluster, save the up/down csv and plots
    invisible(
      finalpwayouts <- lapply(pwayruns, function(pwaycat){
        
        
        
        # message(pwaycat)
        
        
        subcatout <- paste0(pwayoutdir, '/', pwaycat, '/')
        
        # dir.create(subcatout) --> do this with recursive later, maybe prevent even making it if all don't work
        
        clustres <- pathway_analysis_mainlist[[pwaycat]]
        
        clusters <- names(clustres)
        
        
        #for each cluster, get csv and save
        numpways <- lapply(clusters, function(clust){
          
          
          
          
          pwayres_DE_across_conditions_per_cluster <- clustres[[clust]]
          
          
          
          if(is.null(pwayres_DE_across_conditions_per_cluster)){return()}
          
          #this is the table
          signres <- pwayres_DE_across_conditions_per_cluster
          
          #save csv
          
          subcatout_clustdir <- paste0(subcatout, '/', clust, '/')
          
          suppressWarnings(dir.create(subcatout_clustdir, recursive = T))
          
          
          #csv file
          subcatout_clustdir_gseares <- paste0(subcatout_clustdir, '/pathwaytable.csv')
          
          
          write.csv(signres, subcatout_clustdir_gseares, quote = F, row.names = F)
          
          
          
          
          nrow(signres)
          
          
          
          
          
          
          
        } ) # close clusters lapply
        
        
      }) # close saving loop for all categories
      
    ) # close invisible wrap around lapply
    
    
    
    
    
    
    
    
    
    
    ### close any open devices
    while (!is.null(dev.list()))  dev.off()
    
    
    
    
    return(pathway_analysis_mainlist)
    
    
  }) # close cross condition lapply
  
  
  names(pathway_analysis_mainlist_comps) <- comps$labels
  
  
  
  
  #remove big objects
  rm(pathways)
  invisible(gc(full = T, reset = F, verbose = F))
  
  
  
  
  pathway_analysis_mainlist_comps
  
  
  ### for easily reproducing plots and etc, save them as R objects...
  pwayoutdir <- paste0(outdir_int, '/pathwayanalysis_crosscondition_ORA/')
  DE_pathways_plot_objects_list_file <- paste0(pwayoutdir, '/DE_ORA_list_object.rds')
  
  saveRDS(pathway_analysis_mainlist_comps, DE_pathways_plot_objects_list_file)
  
  
  
  return(pathway_analysis_mainlist_comps)
  
  
}

















#' Compositional analysis comparing proportional abundnace across conditions for integrated Seurat objects
#'
#' This is a modular component of the scDAPP scRNAseq pipeline. Perform compositional analysis across conditions to compare the proportion of cell types. Supports multiple condtions (A vs B vs C). Supports "pseudobulk" replicate-aware analysis as implemented in the propeller test with arcsin transformation, as recommended by Simmons 2022 (https://doi.org/10.1101/2022.02.04.479123). Alternatively supports old-school, non-replicate aware chisq test as implemented by the `prop.test()` function.
#'
#' @param sobjint integrated Seurat object. Metadata should have two columns: "Condition" corresponding to the A vs B conditions to compare across, and "Code" corresponding to sample / replicate names. A third column for clusters or celltypes should also be in the metadata and the name of that column will be passed to the `grouping_variable` parameter.
#' @param comps data.frame with two columns called c1, c2; these are the conditions you want to set up. Should be present in psuedobulk_metadata Condition column.
#' @param sample_metadata data.frame with three columns called Sample, Condition, Code.
#' @param outdir_int string, path to save results to. Will create a sub-directory called "compositional_proportion_analysis" and save inside of there.
#' @param grouping_variable string, column name of identity in Seurat object meta.data to stratify DE by. For example, clusters or celltype. Will perform A vs B DE in each of these groupings. Default is "seurat_clusters"
#' @param compositional_test string, which test to use, either "propeller" or "chisq"; will use chisq if not set and issue a warning
#' @param fill_barplots T/F, whether to "fill" the annotation barplots on the side of the heatmap, normalizing to 1; default = T
#'
#' @return will return a list with two items: the A vs B compositional analysis in the "composition_comps" element, with a sub-element for each comparison, and some global compositional information in the "globalcomposition" element.
#' @export
#'
#' @examples
#' \dontrun{
#'
#'
#' # `sample_metadata` looks like this:
#' Sample,Condition,Code
#' SampleXYZ1,Control,Control1
#' SampleXYZ2,Control,Control2
#' SampleABC1,KO1,KO1_1
#' SampleABC2,KO1,KO1_2
#' SampleJKL1,KO2,KO2_1
#' SampleJKL2,KO2,KO2_1
#'
#'
#' # `comps` looks like this:
#' c1,c2
#' KO1,Control
#' KO2,Control
#' KO1,K2
#'
#'
#' ## Run the analysis ##
#' comp_result <- compositional_analysis_module(sobjint,
#' comps,
#' sample_metadata,
#' outdir_int,
#' grouping_variable,
#' compositional_test = 'propeller')
#'
#' #Get the output for KO1 vs Control; replace KO1_vs_Control to get the other comparisons
#' # differential composition results - heatmap
#' comp_result$composition_comps$KO1_vs_Control$hmprop_comp
#'
#' # differential composition results - table
#' comp_result$composition_comps$KO1_vs_Control$compres
#'
#' #Get some global information including cell numbers, proportions
#' # table of cell numbers
#' comp_result$globalcomposition$cellstab
#'
#' #table of cell proportions, ie the table above with each column divided by column sum
#' comp_result$globalcomposition$proptab
#'
#' #heatmap of cell proportions, without any statistical testing
#' comp_result$globalcomposition$hmprop
#' }
#'
compositional_analysis_module <- function(sobjint,
                                          comps,
                                          sample_metadata,
                                          outdir_int,
                                          grouping_variable,
                                          compositional_test,
                                          fill_barplots
){

  require(ComplexHeatmap)
  require(circlize)
  require(speckle) #package with propeller test

  if(missing(grouping_variable)){grouping_variable = "seurat_clusters"}
  if(missing(compositional_test)){warning("No compositional test selected, will use chisq test"); compositional_test = 'chisq'}
  if(missing(fill_barplots)){fill_barplots = T}


  #prep names
  comps$labels <- paste0(comps$c1, '_vs_', comps$c2)

  outdir_comp <- paste0(outdir_int, '/compositional_proportion_analysis/')
  dir.create(outdir_comp, recursive = T)

  ### get the composition table

  #num cells per cluster table
  md <- sobjint@meta.data
  cellstab <- table(md[,grouping_variable], md$Code)



  #prop table, divide num cells by total cells for each samp
  # ensure no div by zero
  cellstab2 <- cellstab[,!Matrix::colSums(cellstab) == 0]
  sample_metadata <- sample_metadata[sample_metadata$Code %in% colnames(cellstab2),]
  proptab <- t(t(cellstab2)/Matrix::colSums(cellstab2))

  #make a heatmap of the prop table
  col_fun = circlize::colorRamp2(c(0, max(proptab)), c( "white", "red"))

  hmprop <- ComplexHeatmap::Heatmap(proptab, name = "Proportion", col = col_fun,
                                    rect_gp = gpar(col = "black", lwd = 0.1),
                                    border_gp = gpar(col = "black", lwd = 1),
                                    column_split = sample_metadata$Condition,
                                    cluster_rows = T, cluster_columns = F,
                                    cell_fun = function(j, i, x, y, width, height, fill) {
                                      grid.text(sprintf("%.2f", proptab[i, j]), x, y, gp = gpar(fontsize = 10))
                                    })



  ## save global tables/heatmap
  cellnumfile <- paste0(outdir_comp, '/NumberCells.csv')
  write.csv(cellstab, cellnumfile, quote = F)

  cellpropfile <- paste0(outdir_comp, '/ProportionCells.csv')
  write.csv(proptab, cellpropfile, quote = F)

  hmprop_file <- paste0(outdir_comp, '/HeatmapProportions.pdf')
  pdf(hmprop_file, height = 7, width = 7)
  print(hmprop)
  dev.off()


  # proplong <- reshape2::melt(proptab)
  # proplong$Var1 <- factor(proplong$Var1, levels = str_sort(unique(proplong$Var1), numeric = T))

  # ggplot(proplong, aes(Var1, value, fill=Var2))+
  #   geom_col(position = 'fill')





  #for each comparison, do the compositional analysis

  compslen <- 1:nrow(comps)
  compidx = 1 #for testing





  composition_comps <- lapply(compslen, function(compidx){

    #get comparison condition levels
    c1 <- comps[compidx,1]
    c2 <- comps[compidx,2]

    #get comp lab
    lab <- comps[compidx,3]

    message(lab)


    #get just this comp samples
    subpmd <- sample_metadata[sample_metadata$Condition %in% c(c1,c2),]



    #order columns of heatmap by c1 vs c2
    code_comps_order <- subpmd[subpmd$Condition == c1,"Code"]
    code_comps_order <- c(code_comps_order, subpmd[subpmd$Condition == c2,"Code"])

    #also prepare factor for ordering of heatmap
    condition_vector_ordering <- factor(subpmd[match(code_comps_order, subpmd$Code), "Condition"], levels = c(c1,c2))



    ### use propeller if multiple samples ###

    if(compositional_test == 'propeller'){
      #subet seurat for metadata
      bigmd <- sobjint@meta.data
      md <- bigmd[bigmd$Condition %in% c(c1,c2),]

      subpmd <- sample_metadata[sample_metadata$Condition %in% c(c1,c2),]

      #remove empty levels from all three vars
      md[,grouping_variable] <- factor(md[,grouping_variable],
                                       levels = stringr::str_sort(unique(md[,grouping_variable]), numeric = T) )

      md$Code <- factor(md$Code,
                        levels = unique(subpmd$Code))

      md$Condition <- factor(md$Condition,
                             levels = c(c1,c2))

      #run propeller with ARCSIN (asin) transform, based on paper
      # https://www.biorxiv.org/content/10.1101/2022.02.04.479123v1
      pres <- speckle::propeller(clusters = md[,grouping_variable],
                                 sample = md$Code,
                                 group = md$Condition,
                                 transform = 'asin')


      #get rid of negatives... not sure why this happens...
      pres$PropRatio[sign(pres$PropRatio)==-1] <- pres$PropRatio[sign(pres$PropRatio)==-1] * -1

      #order by diff
      pres <- pres[order(pres$PropRatio, decreasing = T),]

      #subset table of proportions
      comp_proptab <- proptab[,subpmd$Code]


      #make sure columns are in right order
      comp_proptab <- comp_proptab[,code_comps_order]

      #match table of proporitons order with result
      comp_proptab <- comp_proptab[match(pres$BaselineProp.clusters, rownames(comp_proptab)),]

      #order by pval * sign of diff...
      # pres <- pres[ order( -log10(pres$P.Value) * log(pres$PropRatio) , decreasing = T) ,]
      # comp_proptab <- comp_proptab[match(pres$BaselineProp.clusters, rownames(comp_proptab)),]


      #add signficance marks
      pres$BaselineProp.clusters <- as.character(pres$BaselineProp.clusters)
      pres[pres$P.Value<0.05, "BaselineProp.clusters"] <- paste0(
        '* ',
        pres[pres$P.Value<0.05, "BaselineProp.clusters"],
        ' *'
      )

      rownames(comp_proptab) <- pres$BaselineProp.clusters


      # row annotation: pvalue (eprecated)
      #deal with underflow
      # neglog <- -log10(pres$P.Value)
      # if(any(neglog == Inf)){neglog[neglog==Inf] <- 210}
      # row_ha = rowAnnotation( "-log10P" = anno_barplot( neglog ) )

      #negative p value * propratio
      # directional_pval <- -log10(pres$P.Value) * sign(log(pres$PropRatio))
      # row_ha = rowAnnotation( "-log10P * PropRatio" = anno_barplot( directional_pval ) )

      #row annotation: PropRatio
      # row_ha = rowAnnotation( "LogPropRatio" = anno_lines( log(pres$PropRatio), smooth = T,
      #                                                      axis_param = list(direction = "reverse"))
      #                          )

      #row annotation: barplots of proportions
      propmat_annot <- pres[,3:4]
      if(fill_barplots==T){
        propmat_annot[sign(propmat_annot)==-1] <- propmat_annot[sign(propmat_annot)==-1] * -1 # deal with negatives; why does this happen? must be all zero... cannot reproduce...
        propmat_annot$sum <- rowSums(propmat_annot)
        propmat_annot[,1] <- propmat_annot[,1] / propmat_annot[,3]
        propmat_annot[,2] <- propmat_annot[,2] / propmat_annot[,3]
      }
      propmat_annot <- as.matrix(propmat_annot[,1:2])
      row_ha = rowAnnotation( "Proportions" = anno_barplot( propmat_annot,
                                                            # axis_param = list(direction = "reverse"),
                                                            gp = gpar(fill = c('firebrick', 'steelblue'), col = c('firebrick', 'steelblue') ))
      )


      hmprop_comp <- Heatmap(comp_proptab, name = "Proportion", col = col_fun,
                             rect_gp = gpar(col = "black", lwd = 0.1),
                             border_gp = gpar(col = "black", lwd = 1),
                             column_split = condition_vector_ordering,
                             cluster_rows = F, cluster_columns = F,
                             right_annotation = row_ha,
                             # width = ncol(comp_proptab)*unit(20, "mm"),
                             cell_fun = function(j, i, x, y, width, height, fill) {
                               grid.text(sprintf("%.2f", comp_proptab[i, j]), x, y, gp = gpar(fontsize = 10))
                             })


      compres <- pres

    } else{



      # if no replicates, use prop.test

      #calculate cells per condition w/o regard to sample replicates
      md <- sobjint@meta.data
      md <- md[md$Condition %in% c(c1,c2),]
      cells_condtab <- table(md[,grouping_variable], md$Condition)

      tots <- table(md$Condition)
      c1tot <- tots[c1]
      c2tot <- tots[c2]

      cells_condtab <- cells_condtab[,c(c1,c2)]

      #remove all zero clusters
      cells_condtab = cells_condtab[Matrix::rowSums(cells_condtab)>0,]


      #in each cluster, calculate proportion of c1 and c2 from overall c1 and c2, then compare with prop.test
      clust_proptest_res <- lapply(1:nrow(cells_condtab), function(i){

        clust <- rownames(cells_condtab)[i]

        #get this clusters' props
        vec <- cells_condtab[i,,drop=F]

        x = c(vec[1,1], vec[1,2])
        n = c(c1tot, c2tot)

        pt <- prop.test(x=x,n=n)

        ptdf <- data.frame(cluster = clust,
                           c1prop = pt$estimate[1],
                           c2prop = pt$estimate[2],
                           asin_ratio = plogis(pt$estimate[1]) / plogis(pt$estimate[2]),
                           difference = pt$estimate[1] - pt$estimate[2],
                           p = pt$p.value,
                           row.names = NULL)


      })




      #make results to data.frame
      clust_proptest_resdf <- dplyr::bind_rows(clust_proptest_res)

      # asin ratio = Inf, means c2 was zero..
      clust_proptest_resdf$asin_ratio[clust_proptest_resdf$asin_ratio == Inf] <- 1

      #add FDR
      clust_proptest_resdf$FDR <- p.adjust(clust_proptest_resdf$p)

      #sort by directional pvalue?
      # dirpval <- -log1p(clust_proptest_resdf$p) * sign(clust_proptest_resdf$difference)
      # dirpval <- order(dirpval, decreasing = T)
      # clust_proptest_resdf <- clust_proptest_resdf[order(clust_proptest_resdf$difference, decreasing = T),]

      #sort by difference, deprecated
      clust_proptest_resdf <- clust_proptest_resdf[order(clust_proptest_resdf$difference, decreasing = T),]

      #sort by "asin ratio"
      # clust_proptest_resdf <- clust_proptest_resdf[order(clust_proptest_resdf$asin_ratio, decreasing = T),]

      #subset table of proportions
      comp_proptab <- proptab[,subpmd$Code]

      #make sure columns are in right order
      comp_proptab <- comp_proptab[,code_comps_order]

      #match table of proporitons order with result
      comp_proptab <- comp_proptab[match(clust_proptest_resdf$cluster, rownames(comp_proptab)),]

      #put asterisk if significant
      clust_proptest_resdf[clust_proptest_resdf$p < 0.05, "cluster"] <- paste0( '* ', clust_proptest_resdf[clust_proptest_resdf$p < 0.05, "cluster"], ' *')

      rownames(comp_proptab) <- clust_proptest_resdf$cluster


      #annotate with pvalue, deprecated
      # #deal with underflow
      # neglog <- -log10(clust_proptest_resdf$p)
      # if(any(neglog == Inf)){neglog[neglog==Inf] <- 210}
      # row_ha = rowAnnotation( "-log10P" = anno_barplot( neglog ) )

      #annotate with pval * sign of difference
      # neglog <- -log10(clust_proptest_resdf$p)
      # if(any(neglog == Inf)){neglog[neglog==Inf] <- 210}
      # directional_pval <- neglog * sign(clust_proptest_resdf$difference)
      # row_ha = rowAnnotation( "-log10P * PropRatio" = anno_barplot( directional_pval ) )

      #row annotation: barplots of proportions
      propmat_annot <- clust_proptest_resdf[,2:3]
      if(fill_barplots==T){
        propmat_annot$sum <- rowSums(propmat_annot)
        propmat_annot[,1] <- propmat_annot[,1] / propmat_annot[,3]
        propmat_annot[,2] <- propmat_annot[,2] / propmat_annot[,3]
      }
      propmat_annot <- as.matrix(propmat_annot[,1:2])
      row_ha = rowAnnotation( "Proportions" = anno_barplot( propmat_annot,
                                                            # axis_param = list(direction = "reverse"),
                                                            gp = gpar(fill = c('firebrick', 'steelblue'), col = c('firebrick', 'steelblue') ))
      )


      hmprop_comp <- Heatmap(comp_proptab, name = "Proportion", col = col_fun,
                             rect_gp = gpar(col = "black", lwd = 0.1),
                             border_gp = gpar(col = "black", lwd = 1),
                             column_split = condition_vector_ordering,
                             cluster_rows = F, cluster_columns = F,
                             right_annotation = row_ha,
                             # width = ncol(comp_proptab)*unit(20, "mm"),
                             cell_fun = function(j, i, x, y, width, height, fill) {
                               grid.text(sprintf("%.2f", comp_proptab[i, j]), x, y, gp = gpar(fontsize = 10))
                             })



      compres <- clust_proptest_resdf

    }




    ### save it ###
    if(compositional_test == 'propeller'){
      outdir_comp_thiscomparison <- paste0(outdir_comp,
                                           '/Propellertest/',
                                           lab, '/'
      )
    } else{
      outdir_comp_thiscomparison <- paste0(outdir_comp,
                                           '/2PropZtest/',
                                           lab, '/'
      )
    }

    dir.create(outdir_comp_thiscomparison, recursive = T)

    compresfile <- paste0(outdir_comp_thiscomparison, '/CompositionAnalysis.csv')
    write.csv(compres, file = compresfile, quote = F)

    hmprop_comp_file <- paste0(outdir_comp_thiscomparison, '/HeatmapProportions.pdf')
    pdf(hmprop_comp_file, height = 7, width = 7)
    print(hmprop_comp)
    dev.off()

    return(list(compres = compres,
                hmprop_comp = hmprop_comp))

  })

  names(composition_comps) <- comps$labels

  #also keep other prop stuff in a nice list
  globalcomposition <- list(cellstab = cellstab,
                            proptab = proptab,
                            hmprop = hmprop)

  comp_out <- list(composition_comps = composition_comps,
                   globalcomposition = globalcomposition)


  return(comp_out)

}








# # ### testing ###
# #
# library(Seurat)
# library(tidyverse)
# library(edgeR)
#
# # rm(sobj)
# sobjint <- readRDS('~/Dropbox/Result_from_Alex/deyoudata/hto_prostate/outs/scRNAseqpipeline/multisample_integration/data_objects/Seurat-object_integrated.rds')
# sample_metadata <- read.csv('~/Dropbox/Result_from_Alex/deyoudata/hto_prostate/data/metadata/sample_metadata.csv')
# comps <- read.csv('~/Dropbox/Result_from_Alex/deyoudata/hto_prostate/data/metadata/comps.csv')
# grouping_variable <- 'seurat_clusters'
# Pseudobulk_mode <- T
# outdir_int <- 'tests/de_module_test' ; dir.create(outdir_int, recursive = T)
# cluster_prefix = T
# # min.pct <- 0
# assay = 'RISC'
# slot = 'data'
# crossconditionDE_padj_thres = 0.1
# crossconditionDE_lfc_thres = 0
# crossconditionDE_min.pct = 0.1
#
#
# set.seed(2022)
#
#
# m_bycluster_crosscondition_de_comps <- de_across_conditions_module(
#   sobjint = sobjint,
#   sample_metadata = sample_metadata,
#   comps = comps,
#   outdir_int = outdir_int,
#   grouping_variable = 'seurat_clusters',
#   Pseudobulk_mode = T
#
# )
#
#
# m_bycluster_crosscondition_de_comps_WILCOX <- de_across_conditions_module(
#   sobjint = sobjint,
#   sample_metadata = sample_metadata,
#   comps = comps,
#   grouping_variable = 'seurat_clusters',
#   Pseudobulk_mode = T
#
# )
#
#
#
# library(Seurat)
# library(tidyverse)
# library(fgsea)
# library(foreach)
# library(doParallel)
# library(parallel)
#
# rm(sobj)
# rm(sobjint)
#
# sample_metadata <- read.csv('~/Dropbox/Result_from_Alex/deyoudata/hto_prostate/data/metadata/sample_metadata.csv')
# comps <- read.csv('~/Dropbox/Result_from_Alex/deyoudata/hto_prostate/data/metadata/comps.csv')
#
#
# outdir_int <- 'tests/de_module_test' ; dir.create(outdir_int, recursive = T)
#
# deg.weight = 'auto'
# pathway_padj_thres <- 0.1
# workernum <- 4
# cp.font.size <- 5
# pwaycats <- c("HALLMARK", "GO_BP", "GO_MF", "GO_CC", "CP_REACTOME", "CP_KEGG", "TFT_GTRD", "TFT_TFT_Legacy")
#
#
#
#
#
#
# pathways <- preppathways_pathwayanalysis_crosscondition_module(species = species,
#                                                                outdir_int = outdir_int)
#
#
#
#
# pways_output_list <- pathwayanalysis_crosscondition_module(
#   m_bycluster_crosscondition_de_comps = m_bycluster_crosscondition_de_comps,
#   pathways = pathways,
#   sample_metadata = sample_metadata,
#   comps = comps,
#   workernum = 6,
#   outdir_int = outdir_int
# )
#
# beepr::beep()




# pways_output_list <- pathwayanalysis_crosscondition_module(
#   m_bycluster_crosscondition_de_comps = m_bycluster_crosscondition_de_comps,
#   pathways = pathways,
#   sample_metadata = sample_metadata,
#   deg.weight = "pval",
#   comps = comps,
#   workernum = workernum,
#   outdir_int = outdir_int
# )
