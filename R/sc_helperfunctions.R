# https://r-pkgs.org/whole-game.html


#' Create an alluvial plot from long categorical data
#'
#' Wrapper around ggalluvium package for ggplot based alluvial plot, for quickly making alluvial plot from "long", "raw" categorical data (such as Seurat object meta.data), rather than two-way counts of categories.
#'
#'
#'
#'
#'
#' @param labelsdf data.frame with two columns of raw categorical label: for example, each row is a cell (or other observation), and each column is metadata column 1 and metadata column 2
#' @param repel T/F, whether to repel the labels, default = T
#' @param nudge_x numeric, default nudge to the left and right of the repel labels, default 0.2
#' @param ggfittext T/F - whether to use ggfittext, to try to squeeze or remove tiny stratum labels
#' @param ...
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' #labelsdf can look like this, row.names of cells not needed,
#' # just two columns of categorical data as a data.frame:
#'                      From      To
#' AACCCAAGCATGCGA-1    Malignant  2
#' AAACCCAAGTAGGTTA-1    Malignant  2
#' AAACCCACAAAGCACG-1   Neutrophil  0
#' AAACCCACAGCAGTAG-1    Malignant  2
#' AAACCCACATACCGTA-1    Malignant  2
alluvialplot <- function(labelsdf, repel, nudge_x, ggfittext, ...){


  if( missing(repel) ){ repel = T}
  if( missing(nudge_x) ){nudge_x = 0.3}

  if( missing(ggfittext) ){ggfittext = F}


  #if levels not set, get them by ordering hi > lo

  labelsdf2 <- lapply(labelsdf, function(i){
    if( !is.factor(i) ){
      factor(i, levels = names(sort(table(i),decreasing = T)))
    } else{
      i
    }
  })

  labelsdf <- data.frame(labelsdf2, row.names = rownames(labelsdf))

  require(ggalluvial)

  #for ease, we'll set colnames to from and to
  # colnames(labelsdf)[1:2] <- c('From', 'To')
  # do this with .data trick now to keep colname!


  # turn it into a matrix
  mat <- table(labelsdf[,1], labelsdf[,2])


  #make the table long format
  longfreqs <- reshape2::melt(mat)
  colnames(longfreqs) <- c(colnames(labelsdf)[1:2], 'Freq')


  #factorize, using input levels or existing levels
  # inputting levels is mostly about order.


  longfreqs[,1] <- factor(longfreqs[,1], levels = levels(labelsdf[,1]) )
  longfreqs[,2] <- factor(longfreqs[,2], levels = levels(labelsdf[,2]))

  #if both are numerics, it seems to cause an issue, so convert to char vector...
  # if(is.numeric(longfreqs)[1] & is.numeric(longfreqs)[2])
  # can't reporduce that problem...


  if(ggfittext == T){

    require(ggfittext)
    ap <- ggplot(longfreqs, aes(y = Freq, axis1=.data[[colnames(longfreqs[1])]], axis2= .data[[colnames(longfreqs[2])]] ))+
      geom_alluvium(aes(fill= .data[[colnames(longfreqs[1])]] )) +
      geom_stratum()+
      #geom_label(stat = "stratum", aes(label = after_stat(stratum)))+
      ggfittext::geom_fit_text(stat = "stratum",aes(label = after_stat(stratum)), width = 1/4, min.size = 3) +
      theme_void()


  } else if(repel==T){

    require(ggrepel)

    ap <- ggplot(longfreqs, aes(y = Freq, axis1=.data[[colnames(longfreqs[1])]], axis2= .data[[colnames(longfreqs[2])]] ) )+
      scale_x_discrete(expand = c(.4, 0))+
      geom_alluvium( aes(fill= .data[[colnames(longfreqs[1])]] ), width = 1/4 ) +
      geom_stratum(width = 1/4) +
      scale_linetype_manual(values = c("blank", "solid")) +

      ggrepel::geom_label_repel(
        aes(label = .data[[colnames(longfreqs[1])]] ),
        stat = "stratum", nudge_x = nudge_x * -1, ...) +

      ggrepel::geom_label_repel(
        aes(label = .data[[colnames(longfreqs[2])]]),
        stat = "stratum", nudge_x = nudge_x, ...) +
      theme_void()

  } else{

    ap <- ggplot(longfreqs, aes(y = Freq, axis1=.data[[colnames(longfreqs[1])]], axis2= .data[[colnames(longfreqs[2])]] ) )+
      geom_alluvium(aes(fill= .data[[colnames(longfreqs[1])]] )) +
      geom_stratum()+
      geom_label(stat = "stratum", aes(label = after_stat(stratum)))+
      # ggfittext::geom_fit_text(stat = "stratum",aes(label = after_stat(stratum)), width = 1/4, min.size = 3) +
      theme_void()


  }
  ap

}











#' Pseudobulk Seurat objects at whole sample or celltype level
#'
#' Pseudobulking is performed adding up gene expression values for each cell, for either cell types or whole sample.
#'
#' @param obj - seurat object, or a matrix
#' @param grouping_colname_in_md - optional, a string, the column name of `sobj@meta.data` (if obj is a seurat object) or metadata (if using matrix and metadata input) to use as "cell types" or any other grouping for pseudobulking at group level. if not provided, pseudobulk the entire matrix. default, not used.
#' @param metadata - data.frame with cell metadata, similar to `seuratobject@meta.data`. pass this only if
#' @param rawh5_path - optional, a string, the path to a rawH5 file, if provided will use the raw matrix subsetted by cells in sobj; if not will just pseudobulk the seurat object. useful if some filtering was applied to seurat object but you want to pseudobulk the whole matrix without that filtering, but with only using cells in seurat object. Default is not to use this.
#' @param assay - optional, a string, the name of the Seurat object assay to pull matrix from if rawh5_path is not provided. Default is "RNA" assay
#' @param slot - optional, a string, the name of the Seurat object slot within the designated object assay to pull matrix from if rawh5_path is not provided. Default is "counts" slot
#' @param min_cells - numeric, minimum number of cells in the cluster to pseudobulk. default is 7 cells.
#'
#' @return a data.frame. if grouping_colname_in_md is provided, each celltype will have a pseudobulked column, if not the data.frame will just be one column for the whole sample matrix.
#' @export
#'
#' @examples
pseudobulk <- function(obj, grouping_colname_in_md, metadata, rawh5_path, assay, slot, min_cells){

  if(missing(assay)){assay = 'RNA'}
  if(missing(slot)){slot = 'counts'}
  if(missing(min_cells)){min_cells <- 7}

  require(Seurat)
  require(Matrix)

  #if rawh5_path is given, read in from raw data for all genes
  # if not, just use the seurat object as is

  #if grouping_colname_in_md is given, pseudobulk at celltype level
  # if not, pseudobulk whole object


  #get matrix and md


  if( any(grepl('Seurat', is(obj), ignore.case = T))  ){

    message('Seurat object detected')
    #md from seurat obj
    sobj <- obj
    md <- sobj@meta.data

    #mat: read in H5, or use
    if( !missing(rawh5_path) ){

      message('Reading raw matrix from:\n', rawh5_path)

      #read in raw mat
      mat  <- Read10X_h5(rawh5)

      mat <- mat[, match(colnames(sobj), colnames(mat)) ]

    } else{

      message('Using matrix from Seurat object:',
              '\n - Assay = ', assay,
              '\n - Slot = ', slot)

      mat <- Seurat::GetAssayData(sobj, assay=assay, slot=slot)
    }
  } else{
    message('Assuming input is matrix-like')

    mat <- obj

    if(!missing(grouping_colname_in_md)){
      if(missing(metadata)){stop('Please pass metadata dataframe to "metadata" argument')} else{md = metadata}
    }

  }




  #pseudobulk (at whole or celltype level)

  if( !missing(grouping_colname_in_md) ){
    message('For grouping, using metadata column: "', grouping_colname_in_md, '"')

    #get celltypes by order of number hi-->lo
    if( is.factor(md[,grouping_colname_in_md]) ){cts <- levels(md[,grouping_colname_in_md])} else{
      cts <- names( sort(table(as.vector(md[,grouping_colname_in_md])), decreasing = T) )
    }

    ## skip if zero cells...
    cttab <- table(md[,grouping_colname_in_md])
    cttab <- cttab[cts]
    cts <- cts[cttab>0]



    #skip if below min cells
    ## skip if zero cells...
    cttab <- table(md[,grouping_colname_in_md])
    cttab <- cttab[cts]
    cts <- cts[cttab>min_cells]

    #for each cell type, pseudobulk
    dflist <- lapply(cts, function(ct){

      #subset md
      md_ct <- md[md[,grouping_colname_in_md]==ct,]

      #subset mat; drop =F applies if only 1 cell is left...
      mat_ct <- mat[,match(rownames(md_ct), colnames(mat)), drop=F]

      df <- data.frame(Matrix::rowSums(mat_ct))
      colnames(df) <- ct
      df
    })

    df <- dplyr::bind_cols(dflist)

  } else{
    message('Groupings not provided, will pseudobulk whole dataset')


    df <- data.frame(Matrix::rowSums(mat))
    colnames(df) <- NULL

  }


  df


}









#' Calculate hemoglobin features in Seurat object
#'
#' @param sobj Seurat object
#' @param hemoglobin.features character vector, hemoglobin gene symbols to search in dataset, default will search all mouse and human hemoglobin gene symbols
#'
#' @return data.frame with percent.hemoglobin in all cells
#' @export
#'
#' @examples
#' sobj[['percent.hemoglobin']] <- Seurat::PercentageFeatureSet(sobj, features = hemoglobin.features)
calculate_percent.hemoglobin <- function(sobj, hemoglobin.features){


  # is missing, check all human / mouse hemoglobin genes
  if( missing(hemoglobin.features) ){
    hemoglobin.features <- c('HBA1', 'HBA2', 'HBB', 'HBD', 'HBE1', 'HBG1', 'HBG2', 'HBM', 'HBQ1', 'HBZ',
                             'Hba', 'Hba-a1', 'Hba-a2', 'Hba-ps3', 'Hba-ps4', 'Hba-x', 'Hbb', 'Hbb-ar', 'Hbb-b1', 'Hbb-b2', 'Hbb-bh0', 'Hbb-bh1', 'Hbb-bh2', 'Hbb-bh3', 'Hbb-bs', 'Hbb-bt', 'Hbb-y')
  }


  message('Calculating percent.hemoglobin')
  #hemoglobin content, add to metadata

  #hemoglobin features, all mouse and human hemoglobin genes are searched for, or hemoglobins are user-provided
  hemoglobin.features <- hemoglobin.features[hemoglobin.features %in% rownames(sobj)]

  hb <- Seurat::PercentageFeatureSet(sobj, features = hemoglobin.features)

  message('Returning percent.hemoglobin as data.frame column\nMake sure to add to sobj$percent.hemoglobin')

  return(hb)

}





#' Make a table that can print to a pdf page
#'
#' @param tabledf data.frame to put as table on pdf
#' @param title title for table
#' @param titlesize title font size, default is 15
#' @param padding whitespace between title and table, default=1
#'
#' @return
#' @export
#'
#' @examples
pdftable <- function(tabledf, title, titlesize, padding){
  
  require(grid)

  if(missing(titlesize)){titlesize <- 15}
  if(missing(padding)){padding <- 1}


  table <- gridExtra::tableGrob(tabledf)


  if(missing(title)){

    grid::grid.newpage()

    return(grid::grid.draw(table))

  } else{


    #set up title and table
    title <- grid::textGrob(label = title,
                            gp=grid::gpar(fontsize=titlesize) )


    #add padding and table
    # https://stackoverflow.com/a/33738678
    padding <- unit(padding,"line")

    table <- gtable::gtable_add_rows(
      table, heights = grid::grobHeight(title) + padding, pos = 0
    )
    table <- gtable::gtable_add_grob(
      table, list(title),
      t = 1, l = 1, r = ncol(table)
    )

    grid::grid.newpage()

    return(grid::grid.draw(table))

  }
}







#' a minimalistic theme for publication-quality nonlinear dimensional reduction plots like tSNE/UMAP
#'
#' numeric axes are supposed to be not really meaningful in tSNE/UMAP
#'
#' @param titlesize numeric, size of font
#'
#' @return
#' @export
#'
#' @examples
#' Seurat::DimPlot(object = seurat_object) + theme_dimplot()
theme_dimplot <- function(titlesize = 15) {

  base_size = 11
  base_family = ""
  base_line_size = base_size/22
  base_rect_size = base_size/22

  theme_grey(base_size = base_size, base_family = base_family,
             base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace%

    theme(panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour = "black"), panel.grid = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey85", colour = "grey20"),
          legend.key = element_rect(fill = "white", colour = NA), complete = TRUE,
          axis.ticks = element_blank(), axis.text = element_blank(),
          plot.title = element_text(hjust = 0.5, vjust = 1, face = 'bold', size = titlesize)
    )
}




#' Fix underflow in a sorted list of -log pvalues
#'
#' Given a list of P values from Seurat Wilcox or EdgeR test, one should multiply P values times -log10 times sign of lfc. With this input, sometimes underflow occurs. This function will fix this underflow (ie for input to GSEA) by adding +1 to all Inf values and -1 to all -Inf values.
#'
#' @param scores named vector of -log or -log10(pvalues) times sign(LFC), sorted hi to low. For the +inf and -inf values, one suggestion is to sort them by logFoldChange. names are gene sybols or IDs, values are -logP values * sign of LFC
#' @param logFC_vec named vector of some effect size, such as log2fc, by which the underflow P value DEGs may be sorted, or else the sorting may be arbitrary. names are gene symbols or IDs, values are L2FC (or some other effect size such as L2FC * pct.diff)
#'
#' @return a vector of P value scores, where underflow -logP values have been ranked appropriated
#' @export
#'
#' @examples
fix_underflow <- function(scores,
                          logFC_vec){


  if(missing(logFC_vec)){logFC_vec = NULL; warning('it would be optimal to provide a vector of logFC values or some other effect size to logFC_vec, or else sorting of top DEGs may be arbitrary')}

  bool = ( scores == Inf | scores == -Inf)
  tbl = table(factor(bool, levels = c(F,T)))

  message(tbl['TRUE'], ' underflow genes detected')


  #underflow, positive
  if(any(scores == Inf)){
    scores_uf <- scores[scores == 'Inf']

    #sort INF values by lfc
    if(!is.null(logFC_vec)){
      logFC_vec_uf <- logFC_vec[names(logFC_vec) %in% names(scores_uf)]
      logFC_vec_uf <- sort(abs(logFC_vec_uf), decreasing = T)
      scores_uf <- scores_uf[match(names(logFC_vec_uf), names(scores_uf))]
    }

    last <- scores[length(scores_uf) + 1]

    replacement <- c()
    for(i in 1:length(scores_uf)){
      to_replace = ifelse(i == 1, yes = last, no = replacement[i-1])
      to_replace = to_replace + 1 #for pos + 1
      replacement <- c(replacement, to_replace)
    }

    # reverse , lo to hi
    replacement <- rev(replacement)
    names(replacement) <- names(scores_uf)


    scores[names(replacement)] <- replacement

  }


  #underflow, negative
  if(any(scores == -Inf)){

    scores <- rev(scores)
    scores_uf <- scores[scores == '-Inf']

    #sort INF values by lfc
    if(!is.null(logFC_vec)){
      logFC_vec_uf <- logFC_vec[names(logFC_vec) %in% names(scores_uf)]
      logFC_vec_uf <- sort(abs(logFC_vec_uf), decreasing = T)
      scores_uf <- scores_uf[names(logFC_vec_uf)]
    }


    last <- scores[length(scores_uf) + 1]

    replacement <- c()
    for(i in 1:length(scores_uf)){
      to_replace = ifelse(i == 1, yes = last, no = replacement[i-1])
      to_replace = to_replace - 1 #for negative only, minus 1...
      replacement <- c(replacement, to_replace)
    }

    #for negative only
    replacement <- rev(replacement)
    names(replacement) <- names(scores_uf)



    scores[names(replacement)] <- replacement

    scores <- rev(scores)

  }

  scores <- sort(scores, decreasing = T)


  return(scores)
}





### try to collapse clusters in seurat object with some cutoff of correlation?
# what's the reference for this?
# "extend" the collapsing?
# A may be cor with B, and B with C, but A not with C...
# collapse A, B and C?

# #collapse clusters: make cor matrix
# avgs <- AverageExpression(sobj, assays = 'SCT', return.seurat = F, slot = 'data')$SCT
#
# #pairwise correlation
# cormat <- cor(as.matrix(avgs))
#
# #collapse them:
# mask <- cormat>correlation_threshold_collapse_hires_clusters
#
# #for each cluster, record which ones have cor > threshold
# clusts <- colnames(mask)
# clusts_to_collapse <- lapply(clusts, function(clust){
#
#   #get correlations of this cluster
#   maskcol <- mask[,clust]
#
#
#   #get any clusts above thres
#   if( any(maskcol) == T ){
#     return( names(maskcol[maskcol==T]) )
#   } else{
#     return()
#   }
#
# })
# names(clusts_to_collapse) <- clusts
#
#
# #check if list elements are null; if so they have no clusters above thres, so drop them
# clusts_to_collapse <- clusts_to_collapse[ sapply(clusts_to_collapse, function(x){length(x)>1}) ]
#
#
# # if there are clusters to collapse, we should get them...
# identdf <- data.frame(orig = names(clusts_to_collapse),
#                       collapse = sapply(clusts_to_collapse,function(x){paste(x,collapse = '_')}) )
#
# identdf
