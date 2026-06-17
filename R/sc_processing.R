# https://r-pkgs.org/whole-game.html


#' Perform data-driven filtering of scRNAseq data
#'
#' Apply simple cutoffs and discover data-driven thresholds for poor quality cells in scRNAseq.
#'
#' Simple cutoffs include minimum number of UMIs, minimum number of unique genes detected, maximum percent mito, and maximum percent hemoglobin.
#' More complex cutoffs are learnt for lower than expected complexity (defined for each cell as num unique genes / num UMIs). Additionally, median absolute deviation is used to exclude remaining cells with high mito content or low UMI content.
#'
#' Specifically, for complexity, a two-part model is used to model log(num Genes) ~ log(num UMIs) for each cell. A linear model and a Loess model are both set up in this way. Outliers with low complexity are called as cells with > 4/n cooks distance cells in the linear model, and low residuals in the loess model. The residual cutoff is set to -3 by default, capturing very low complexity outlier cells.
#'
#' @param sobj seurat object
#' @param min_num_UMI numeric, default is 500, if no filter is desired set to -Inf
#' @param min_num_Feature numeric, default is 200, if no filter is desired set to -Inf
#' @param max_perc_mito numeric, default is 25, if no filter is desired set to Inf
#' @param max_perc_hemoglobin numeric, default is 25, if no filter is desired set to Inf
#' @param loess_negative_residual_threshold numeric, cutoff for loess residuals applied in complexity filtering, default is -5, if you set it high (ie any higher than -2) you will probably remove many good cells.
#' @param mad.score.threshold numeric, default is 3, threshold for median abs deviation thresholding, ie cutoffs set to `median +/- mad * threshold`
#' @param globalfilter.complexity T/F, default T, whether to filter cells with lower than expected number of genes given number of UMIs
#' @param globalfilter.mito T/F, default T, whether to filter cells with higher than normal mito content
#' @param globalfilter.libsize T/F, default T, whether to filter cells with lower than normal UMI content
#'
#' @return a list object.
#'
#' 'cellstatus' = data.frame with cells, filtered out (T/F), filter reason, and other information.
#'
#' 'filtersummary' = small data.frame summarizing the `cellstatus$filterreason` information.
#'
#' 'allcommands' = commands passed to the autofilter function
#'
#' 'baseline_qc_summary' = summarizes distributions of key QC variables
#'
#' 'globalfilter.complexity' = summarizes the complexity filtering with plots and number cells removed
#'
#' 'globalfilter.libsize' = summarizes the libsize filtering with plots and number cells removed
#'
#' 'globalfilter.mito' = summarizes the mito filtering with plots and number cells removed
#' @export
#'
#' @examples
#' # identify outliers
#' af <- autofilter(sobj)
#'
#' # remove outliers
#' goodcells <- af$cellstatus[af$cellstatus$filteredout==F,"barcodes"]
#' sobj <- sobj[,goodcells]
autofilter <- function(
    sobj,

    min_num_UMI,
    min_num_Feature,
    max_perc_mito,
    max_perc_hemoglobin,

    loess_negative_residual_threshold,

    mad.score.threshold,

    globalfilter.complexity,
    globalfilter.mito,
    globalfilter.libsize

){

  #check if seurat object is there
  if( missing( sobj )) {stop('please input a Seurat object')}

  #set the basic cutoffs
  if( missing(min_num_UMI ) ){min_num_UMI <- 500}
  if( missing(min_num_Feature ) ){min_num_Feature <- 200}
  if( missing(max_perc_mito ) ){max_perc_mito <- 25}
  if( missing(max_perc_hemoglobin ) ){max_perc_hemoglobin <- 25}

  #set lowess cutoff
  if( missing( loess_negative_residual_threshold )) {loess_negative_residual_threshold <- -5}

  #set maximum distance of deviations from median tolerated before outlier classification
  if( missing( mad.score.threshold )) {mad.score.threshold <- 3}

  # baseline (global) filtration
  if( missing( globalfilter.complexity )) {globalfilter.complexity <- T}
  if( missing( globalfilter.mito )) {globalfilter.mito <- T}
  if( missing( globalfilter.libsize )) {globalfilter.libsize <- T}



  #output list object; add to this as needed for all filtering.
  reportlist <- list()

  #start keeping a cell status DF
  cellstatus <- data.frame(barcodes = colnames(sobj), filteredout = F, filterreason = 'Unfiltered')

  reportlist[['cellstatus']] <- cellstatus



  #start making summary of removed cells
  reportlist[['filtersummary']] <- list()


  #report commands used
  reportlist[['allcommands']] <- data.frame(
    Command = c("mad.score.threshold", "loess_negative_residual_threshold",
                'min_num_UMI', 'min_num_Feature', 'max_perc_mito', 'max_perc_hemoglobin',
                "globalfilter.complexity", "globalfilter.libsize", "globalfilter.mito"),

    Option = c(mad.score.threshold, loess_negative_residual_threshold,
               min_num_UMI, min_num_Feature, max_perc_mito, max_perc_hemoglobin,
               globalfilter.complexity, globalfilter.libsize, globalfilter.mito)

  )


  if( !('percent.mito' %in% colnames(sobj@meta.data)) ){

    message('Calculating percent.mito')

    #mito content, add to metadata
    mito.features <- grep(pattern = "^mt-", x = rownames(x = sobj), value = TRUE, ignore.case = T)
    sobj[["percent.mito"]] <- Seurat::PercentageFeatureSet(sobj, features = mito.features)

  }


  if( !('percent.hemoglobin' %in% colnames(sobj@meta.data)) ){


    sobj[["percent.hemoglobin"]] <- scDAPP::calculate_percent.hemoglobin(sobj)

  }

  #report baseline stuff
  md <- sobj@meta.data
  reportlist[['baseline_qc_summary']] <- data.frame(summary_nCount_RNA = as.vector(summary(md$nCount_RNA)),
                                                 summary_nFeature_RNA = as.vector(summary(md$nFeature_RNA)),
                                                 summary_perc.mito = as.vector(summary(md$percent.mito)),
                                                 summary_perc.hemoglobin = as.vector(summary(md$percent.hemoglobin)),
                                                 row.names = names(summary(md$nCount_RNA)))


  ### apply cutoffs based on min umi, min gene, max mito, max hemoglobin

  md <- sobj@meta.data

  bad <- md[md$nCount_RNA < min_num_UMI | md$nFeature_RNA < min_num_Feature | md$percent.mito > max_perc_mito | md$percent.hemoglobin > max_perc_hemoglobin,]

  cellstatus[cellstatus$barcodes %in% rownames(bad),"filteredout"] <- T
  cellstatus[cellstatus$barcodes %in% rownames(bad),"filterreason"] <- 'BasicFilter'


  cellstatus$BasicFilter <- 'No'
  cellstatus[cellstatus$filterreason == 'BasicFilter', 'BasicFilter'] <- 'BasicFilter'

  #add reasons for filtering
  #min UMI:
  explanation_string <- paste0('_nUMI-below-', min_num_UMI)
  bad_reason <- rownames( bad[bad$nCount_RNA < min_num_UMI,] )
  cellstatus[ cellstatus$barcodes %in% bad_reason , 'BasicFilter'] <- paste0(cellstatus[ cellstatus$barcodes %in% bad_reason , 'BasicFilter'],
                                                                             explanation_string)

  #min feature:
  explanation_string <- paste0('_nFeature-below-', min_num_Feature)
  bad_reason <- rownames( bad[bad$nFeature_RNA < min_num_Feature,] )
  cellstatus[ cellstatus$barcodes %in% bad_reason , 'BasicFilter'] <- paste0(cellstatus[ cellstatus$barcodes %in% bad_reason , 'BasicFilter'],
                                                                             explanation_string)

  #max mito:
  explanation_string <- paste0('_percent.mito-above-', max_perc_mito)
  bad_reason <- rownames( bad[bad$percent.mito > max_perc_mito,] )
  cellstatus[ cellstatus$barcodes %in% bad_reason , 'BasicFilter'] <- paste0(cellstatus[ cellstatus$barcodes %in% bad_reason , 'BasicFilter'],
                                                                             explanation_string)

  #max hemoglobin:
  explanation_string <- paste0('_percent.hemoglobin-above-', max_perc_hemoglobin)
  bad_reason <- rownames( bad[bad$percent.hemoglobin > max_perc_hemoglobin,] )
  cellstatus[ cellstatus$barcodes %in% bad_reason , 'BasicFilter'] <- paste0(cellstatus[ cellstatus$barcodes %in% bad_reason , 'BasicFilter'],
                                                                             explanation_string)



  #retain good cells
  good <- cellstatus[cellstatus$filteredout==F,"barcodes"]
  sobj <- sobj[,good]


  if(globalfilter.complexity == T){

    message(' - Initiating baseline complexity filtration')

    # complexity = num Genes given num UMIs.
    # nFeature_RNA / nCount_RNA

    # use a linear model to try to ID outliers with low-complexity


    #set up md
    md <- sobj@meta.data

    #calculate and order by complexity
    md$complexity <-   md$nFeature_RNA / md$nCount_RNA
    md <- md[order(md$complexity),]

    #two-step model:
    # loess: must have low negative residuals < -1
    # lm: must have cooks distance > 4 / n

    #lm
    m <- stats::lm(data = md, formula = log(nFeature_RNA) ~ log(nCount_RNA))

    #calculate cooks distance for each point
    cd <- cooks.distance(m)

    #get scaled residuals for each point
    resid_lm <- scale( residuals(m))

    # log transforms are used to flatten variance, make things more normal-looking, bring to similar scale
    m_loess <- stats::loess(data = md, formula = log(nFeature_RNA) ~ log(nCount_RNA))

    #get scaled residuals for each point from LOESS
    resid_loess <- scale( residuals(m_loess))

    out_decision <- data.frame(cd_lm = cd,
                               resid_lm  =resid_lm,
                               resid_loess = resid_loess)

    #to be an outlier, must have scaled loess residuals < -1 and lm cooks distance > 4/n
    outs <- rownames( out_decision[out_decision$cd_lm > (4 / nrow(out_decision)) & out_decision$resid_loess < loess_negative_residual_threshold,] )

    #plot outliers
    md$outlier <- 'nonoutlier'
    md[rownames(md) %in% outs, "outlier"] <- 'outlier'


    #plot relationship with outlier calls
    complexityplot.outliers <- ggplot2::ggplot(md, ggplot2::aes(log(nCount_RNA), log(nFeature_RNA), col = outlier)) +
      ggplot2::geom_point(alpha = 0.7, size = 0.7)+
      ggplot2::scale_color_brewer(palette = 'Set1', direction = -1)+
      ggplot2::labs( caption = paste0('Outliers = Cooks Distance > (4/n)\n& Loess residual < ',loess_negative_residual_threshold) )


    #plot the relationship; we must use log transforms to see things more clearly.
    complexityplot <- ggplot2::ggplot(md, ggplot2::aes(log(nCount_RNA), log(nFeature_RNA), col = complexity)) +
      ggplot2::geom_point(alpha = 0.7, size = 0.7)+
      ggplot2::labs(caption = 'Complexity = nFeature / nCount')


    #plot loess with residuals and lm with cd
    md <- cbind(md, out_decision)
    md$CooksDistance_lm <- md$cd_lm

    cp_lm_cd <- ggplot2::ggplot(md, ggplot2::aes(log(nCount_RNA), log(nFeature_RNA))) +
      ggplot2::geom_point(alpha = 0.7, size = 0.7, ggplot2::aes(col = CooksDistance_lm))+
      ggplot2::geom_smooth(method = 'lm')+
      ggplot2::scale_color_gradient2(high = 'red', low = 'blue', mid = 'grey'  )+
      ggplot2::labs(caption = paste0('Cutoff = cooks > 4 / N, here = ', round(4/nrow(md), 5) ) )

    cp_ls_rd <- ggplot2::ggplot(md, ggplot2::aes(log(nCount_RNA), log(nFeature_RNA))) +
      ggplot2::geom_point(alpha = 0.7, size = 0.7, ggplot2::aes(col = resid_loess))+
      ggplot2::geom_smooth(method = 'loess')+
      ggplot2::scale_color_gradient2(mid = 'grey', low = 'red', high = 'blue')+
      ggplot2::labs(caption = paste('Cutoff = loess residuals < ', loess_negative_residual_threshold ))


    finalplot <- patchwork::wrap_plots(complexityplot,complexityplot.outliers, cp_lm_cd, cp_ls_rd)+
      patchwork::plot_annotation(title = 'Complexity Filtering',
                                 subtitle = paste0('Total Cells: ', ncol(sobj),
                                                   '; Num Outlers: ', length(outs)) )


    #subset and save
    bad <- outs
    filteredcells <- rownames( md[md$outlier == 'nonoutlier',] )

    sobj <- sobj[,filteredcells]

    #report
    cellstatus[cellstatus$barcodes %in% bad,"filteredout"] <- T
    cellstatus[cellstatus$barcodes %in% bad,"filterreason"] <- 'globalfilter.complexity'


    reportlist[["globalfilter.complexity"]] <- list(plot = finalplot,
                                                    table = table(md$outlier) )





  }



  #find libsize cutoff low.
  if(globalfilter.libsize == T){

    message(' - Initiating global baseline lib size low filtration')


    #the use of median absolute deviation in outlier detection is inspired by the following:
    # https://www.sciencedirect.com/science/article/abs/pii/S0022103113000668
    # https://eurekastatistics.com/using-the-median-absolute-deviation-to-find-outliers/

    # don't use long-tailed hack.


    # get variable, in this case libsize
    # log transform to deal with extreme tails and negative values
    x <- log( sobj$nCount_RNA )


    # make a mad cutoff; similar to mean +/- sd*2.5, we use median +/- mad*2.5
    # get lib size cutoffs; use a lower cutoff only
    globalfilter.libsize.cutoff.lo <- median(x) - (mad.score.threshold  * mad(x) )

    #try multimode analysis...
    # mm <- multimode::locmodes(log(sobj$nCount_RNA), 2)
    #globalfilter.libsize.cutoff.lo <-mm$locations[1] - ( mm$locations[2] - mm$locations[1] )


    #get cells outside threshold
    bad <- x[x < globalfilter.libsize.cutoff.lo]
    #bad <- x[x < globalfilter.libsize.cutoff.lo | x > globalfilter.libsize.cutoff.hi]

    #exponentiate to get out of log space
    nonexp <- globalfilter.libsize.cutoff.lo
    globalfilter.libsize.cutoff.lo <- exp(globalfilter.libsize.cutoff.lo)


    #get cells within threshold
    filteredcells <- names(x[!names(x) %in% names(bad)] )


    #plot cutoffs for library size / UMIs
    g_lib_hist <- ggplot2::ggplot(sobj@meta.data) +
      ggplot2::geom_histogram(ggplot2::aes(nCount_RNA),
                              color="black", fill = "steelblue",
                              # binwidth = 0.05
      )+
      #ggplot2::geom_vline(xintercept = globalfilter.libsize.cutoff.hi, linetype = "dashed", colour = "red")+
      ggplot2::geom_vline(xintercept =  globalfilter.libsize.cutoff.lo, linetype = "dashed", colour = "red")+
      ggplot2::geom_vline(xintercept = median(sobj@meta.data$nCount_RNA),
                          linetype = "dotted", colour = "red", size = 1.2)+
      ggplot2::scale_x_log10(labels = scales::comma)+
      ggplot2::labs(x = "", y = "Number of cells" ,
                    caption = paste0('Num cells presubset = ', ncol(sobj),
                                     '\nNum cells remaining = ', length(filteredcells))  )+
      ggplot2::theme_linedraw()+
      ggplot2::coord_flip()

    g_lib_vln <- ggplot2::ggplot(sobj@meta.data, ggplot2::aes(x = 0, y = nCount_RNA))+
      ggplot2::geom_violin(fill='steelblue')+
      ggplot2::geom_jitter(height = 0, width = 0.25, size = 0.1)+
      ggplot2::geom_hline(yintercept =  globalfilter.libsize.cutoff.lo, col = 'red', linetype = 'dashed')+
      ggplot2::geom_hline(yintercept = median(sobj@meta.data$nCount_RNA),
                          linetype = "dotted", colour = "red", size = 1.2)+
      ggplot2::theme_linedraw()+
      ggplot2::scale_y_log10(labels = scales::comma)+
      ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                     axis.ticks.x=ggplot2::element_blank())+
      ggplot2::labs(y = "nUMI, log10 scale",
                    caption = paste0("Low UMI cutoff: ", as.character(round(globalfilter.libsize.cutoff.lo, digits = 3)),
                                     "\nMedian lib size = ", median(sobj@meta.data$nCount_RNA)) )+
      ggplot2::xlab('Cells')



    #report
    cellstatus[cellstatus$barcodes %in% names(bad),"filteredout"] <- T
    cellstatus[cellstatus$barcodes %in% names(bad),"filterreason"] <- 'globalfilter.libsize'


    #set up reportables
    finalplot <- patchwork::wrap_plots(g_lib_vln, g_lib_hist)+
      patchwork::plot_annotation(title = 'Num UMI Filtering',
                                 subtitle = paste0('Total Cells: ', ncol(sobj),
                                                   '; Num Outlers: ', length(bad)) )

    finaltable <- data.frame(var = x, var_cond = x < nonexp)
    finaltable$outlier <- 'nonoutlier'
    finaltable[finaltable$var_cond==T,'outlier'] <- 'outlier'
    finaltable <- table(finaltable$outlier)

    reportlist[['globalfilter.libsize']] <- list(plot = finalplot, table = finaltable)


    sobj <- sobj[,filteredcells]

  } #close baseline libsize loop



  #find mito cutoff high
  if(globalfilter.mito == T){

    message(' - Initiating global baseline mito hi filtration')



    # get variable, in this case mito
    x <- sobj$percent.mito

    # make a mad cutoff; similar to mean +/- sd*2.5, we use median +/- mad*2.5
    # get lib size cutoffs; use a high cutoff only
    globalfilter.mito.cutoff <- median(x) + (mad.score.threshold  * mad(x) )
    
    # UPDATE 2024 MARCH 6
    # apply an absolute minimum of 5%
    # if above line picks lower, just use 5%
    globalfilter.mito.cutoff <- ifelse(globalfilter.mito.cutoff <= 5, yes = 5, no = globalfilter.mito.cutoff)


    #get cells above the threshold
    bad <- x[x > globalfilter.mito.cutoff]

    #get cells within threshold
    filteredcells <- names(x[!names(x) %in% names(bad)] )



    g_mito_hist <- ggplot2::ggplot(sobj@meta.data) +
      ggplot2::geom_histogram(ggplot2::aes(percent.mito),
                              color="black", fill = "steelblue",
                              binwidth = 0.5)+
      ggplot2::geom_vline(xintercept =  globalfilter.mito.cutoff, linetype = "dashed", colour = "red")+
      ggplot2::geom_vline(xintercept = median(sobj@meta.data$percent.mito),
                          linetype = "dotted", colour = "red", size = 1.2)+
      ggplot2::labs(y = "Number of cells" , x = '',
                    caption = paste0('Num cells presubset = ', ncol(sobj),
                                     '\nNum cells remaining = ', length(filteredcells))  )+
      ggplot2::theme_linedraw()+
      ggplot2::coord_flip()

    g_mito_vln <- ggplot2::ggplot(sobj@meta.data, ggplot2::aes(x = 0, y = percent.mito))+
      ggplot2::geom_violin(fill='steelblue')+
      ggplot2::geom_jitter(height = 0, width = 0.25, size = 0.1)+
      ggplot2::geom_hline(yintercept = globalfilter.mito.cutoff, col = 'red', linetype = 'dashed')+
      ggplot2::geom_hline(yintercept = median(sobj@meta.data$percent.mito),
                          linetype = "dotted", colour = "red", size = 1.2)+
      ggplot2::scale_y_continuous(breaks = c(0, 25, 50, 75, 100, round(globalfilter.mito.cutoff, digits = 2)))+
      ggplot2::theme_linedraw()+
      ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                     axis.ticks.x=ggplot2::element_blank())+
      ggplot2::labs(
        caption = paste0('Percent.mito cutoff = ', round(globalfilter.mito.cutoff, 2) ,
                         "\nMedian percent.mito = ", round( median(sobj@meta.data$percent.mito), 2))
      )+
      ggplot2::xlab('Cells')




    #report
    cellstatus[cellstatus$barcodes %in% names(bad),"filteredout"] <- T
    cellstatus[cellstatus$barcodes %in% names(bad),"filterreason"] <- 'globalfilter.mito'



    #set up reportables
    finalplot <- patchwork::wrap_plots(g_mito_vln, g_mito_hist)+
      patchwork::plot_annotation(title = 'Percent Mito Filtering',
                                 subtitle = paste0('Total Cells: ', ncol(sobj),
                                                   '; Num Outlers: ', length(bad)) )

    finaltable <- data.frame(var = x, var_cond = x > globalfilter.mito.cutoff)
    finaltable$outlier <- 'nonoutlier'
    finaltable[finaltable$var_cond==T,'outlier'] <- 'outlier'
    finaltable <- table(finaltable$outlier)

    reportlist[['globalfilter.mito']] <- list(plot = finalplot, table = finaltable)


    sobj <- sobj[,filteredcells]


  } #end mito baseline block



  #update reportlist with finalized cellstatus:
  reportlist[["cellstatus"]] <- cellstatus


  # update reportlist with summary of filterng results
  filtersummary <- data.frame(table(cellstatus$filterreason))
  colnames(filtersummary) <- c('FilterReason', 'numCells')

  tot <- data.frame(FilterReason = 'Total', numCells = nrow(cellstatus))

  filtersummary <- rbind(filtersummary, tot)
  reportlist$filtersummary <- filtersummary


  #return whole result list.
  reportlist


}




#' DoubletFinder Wrapper
#'
#'
#' Please see https://github.com/chris-mcginnis-ucsf/DoubletFinder
#'
#' Models homotypic doublets using the following table:
#' https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
#'
#' @param seuratobject A Seurat object, pre-proc with `Seurat::SCTransform()`
#' @param clusters string, should match a column name of seuratobject metadata with clusters or other cell group annotations. default = 'seurat_clusters'
#' @param autofilterres optional, output of `scDAPP::autofilter()` whether to return the the autofilter result with updated `autofilterres$cellstatus` taking into account doublet status
#' @param num.cores integer, num.cores to use for `DoubletFinder::paramSweep`
#'
#' @return if `autofilterres` is provided, it will return `autofilterres` with updated `autofilterres$cellstatus`, if not it will return a data.frame with doublet information and score.
#' @export
#'
#' @examples \dontrun{
#' # With autofilter res, will add in the result to af$cellstatus
#' af <- doubletfinderwrapper(sobj, autofilterres = af, num.cores = 5)
#' goodcells <- af$cellstatus[af$cellstatus$filteredout == FALSE, "barcodes"]
#' sobj <- sobj[, goodcells]
#' }
doubletfinderwrapper <- function(seuratobject, clusters, autofilterres, num.cores, sct){

  if( missing( clusters )) {clusters <- "seurat_clusters"}
  if( !(clusters %in% colnames(seuratobject@meta.data)) ) {stop('No clusters detected in Seurat object')}
  if( missing(sct) ){  sct <- ifelse('SCT' %in% names(seuratobject@assays), yes = T, no = F) }

  if( missing(num.cores) ){num.cores <- 1}

  message('Running DoubletFinder paramsweep (may take a while)')

  #param sweep
  sweepres <- DoubletFinder::paramSweep(seu = seuratobject, PCs = 1:30, num.cores = num.cores, sct=sct)

  message('DF Parameter Sweep Completed')

  sweepstats <- DoubletFinder::summarizeSweep(sweepres)

  pdf(NULL) #prevent plotted output
  bcmvn <- DoubletFinder::find.pK(sweepstats)
  dev.off()

  maxscorepk <- bcmvn[bcmvn$BCmetric == max(bcmvn$BCmetric),2]
  maxscorepk <- as.numeric( levels(maxscorepk)[maxscorepk] )





  #homotypic doublet modelling

  ### using 10x table, use linear regression --> important for predicting homotypic / total doublet number
  #tamlabscpipeline::dratedf
  #dratedf <- read.delim('/Users/ferrenaa/Documents/tam/scripts/doublets/doubletrate.txt', header = T)
  #
  #   dratedf[,1] <- as.numeric(sub("%","",dratedf[,1]))/100
  #
  #   names(dratedf) <- c('MultipletRate', 'CellsLoaded_100%Viability', 'CellsRecovered')

  dratedf <- scDAPP::dratedf

  dbmodel <- lm(MultipletRate ~ CellsRecovered, data = dratedf)

  predicteddoubletrate <- as.numeric((dbmodel$coefficients['CellsRecovered'] * ncol(seuratobject)) + dbmodel$coefficients[1])

  #choose annotations to model homotypic doublets
  homotypicprop <- DoubletFinder::modelHomotypic(    as.vector(  seuratobject@meta.data[,clusters] )   )

  nexppoi <- round(predicteddoubletrate * length(rownames(seuratobject@meta.data)))
  nexppoiadj <- round(nexppoi * (1 - homotypicprop))

  #classify doublets
  message('Final DoubletFinder run:')

  seuratobject <- suppressWarnings(DoubletFinder::doubletFinder(seu = seuratobject, PCs = 1:30, pN = 0.25, pK = maxscorepk, nExp = nexppoi, sct = sct))

  ddf <- data.frame(cells = rownames(seuratobject@meta.data),
                    DoubletFinderClassification = seuratobject@meta.data[,ncol(seuratobject@meta.data)],
                    DoubletFinder_pANN_doublet_score = seuratobject@meta.data[,ncol(seuratobject@meta.data) - 1])

  #if autofilterres is there,add it, if not just return the ddf
  if( missing(autofilterres) ){
    return(ddf)
  } else{

    message('\nAdding DoubletFinder result to autofilterres$cellstatus')
    cellstatus <- autofilterres$cellstatus

    cellstatus$DoubletFinder_pANN_doublet_score <- NA

    cellstatus[match(ddf$cells, cellstatus$barcodes),'DoubletFinder_pANN_doublet_score'] <- ddf$DoubletFinder_pANN_doublet_score

    #mark doublets as filt
    ddf_d <- ddf[ddf$DoubletFinderClassification=='Doublet',]

    cellstatus[match(ddf_d$cells, cellstatus$barcodes),"filteredout"] <- T
    cellstatus[match(ddf_d$cells, cellstatus$barcodes),"filterreason"] <- 'DoubletFinder_doublet'


    autofilterres$cellstatus <- cellstatus

    # update reportlist with summary of filterng results
    filtersummary <- data.frame(table(cellstatus$filterreason))
    colnames(filtersummary) <- c('FilterReason', 'numCells')

    tot <- data.frame(FilterReason = 'Total', numCells = nrow(cellstatus))

    filtersummary <- rbind(filtersummary, tot)
    autofilterres$filtersummary <- filtersummary



    return(autofilterres)



  }

}














### testing below

# hemoglobin.features <- c('HBA1', 'HBA2', 'HBB', 'HBD', 'HBE1', 'HBG1', 'HBG2', 'HBM', 'HBQ1', 'HBZ',
#                          'Hba', 'Hba-a1', 'Hba-a2', 'Hba-ps3', 'Hba-ps4', 'Hba-x', 'Hbb', 'Hbb-ar', 'Hbb-b1', 'Hbb-b2', 'Hbb-bh0', 'Hbb-bh1', 'Hbb-bh2', 'Hbb-bh3', 'Hbb-bs', 'Hbb-bt', 'Hbb-y')
#
#
# min_num_UMI <- 1000
# min_num_Feature <- 200
# max_perc_mito <- 25
# max_perc_hemoglobin <- 25
#
# loess_negative_residual_threshold <- -3
# mad.score.threshold = 2.5
# globalfilter.complexity <- T
# globalfilter.mito <- T
# globalfilter.libsize <- T
#
# rawh5 <- '~/Dropbox/data/bangdata/scrnaseq-TKO-DKOAA-DKO/rawdata/Sample-06_TL494/filtered_feature_bc_matrix.h5'
# rawh5 <- '~/Dropbox/data/bangdata/scrnaseq-TKO-DKOAA-DKO/rawdata/Sample-04_DJ582M11/filtered_feature_bc_matrix.h5'
#
# library(tidyverse) ; library(Seurat) ; library(FerrenaSCRNAseq)
# sobj <- CreateSeuratObject(   Read10X_h5(rawh5), min.cells= 3)
#
# af <- autofilter(sobj)
#
#
# goodcells <- af$cellstatus[af$cellstatus$filteredout==F,"barcodes"]
#
# sobj <- sobj[,goodcells]
#
# #normalize and cluster
# suppressWarnings(sobj <- Seurat::SCTransform(sobj, verbose = T, method="glmGamPoi"))
#
# sobj <- Seurat::RunPCA(object = sobj, verbose = F)
#
# sobj <- Seurat::FindNeighbors(object = sobj, dims = 1:30, verbose = F)
# sobj <- Seurat::FindClusters(object = sobj, resolution = 0.1, verbose = F, algorithm = 1)
#
# sobj <- Seurat::RunUMAP(sobj, dims = 1:30)
#
# # get the doubletfinder  dataframe
# dfdf <- scDAPP::doubletfinderwrapper(sobj, #autofilterres = af,
#                              num.cores = 5)
#
# af <- doubletfinderwrapper(sobj, autofilterres = af,
#                            num.cores = 5)


#' Add per-sample QC metadata columns to a Seurat object
#'
#' Computes percent mitochondrial and hemoglobin content and cell-cycle phase
#' when those columns are not already present in object metadata.
#'
#' @param sobj Seurat object
#' @return Seurat object with QC metadata columns added when missing
#' @export
add_per_sample_qc_metadata <- function(sobj) {
  if (!("percent.mito" %in% colnames(sobj@meta.data))) {
    mito.features <- grep(
      pattern = "^mt-",
      x = rownames(x = sobj),
      value = TRUE,
      ignore.case = TRUE
    )
    sobj[["percent.mito"]] <- Seurat::PercentageFeatureSet(sobj, features = mito.features)
  }

  if (!("percent.hemoglobin" %in% colnames(sobj@meta.data))) {
    sobj$percent.hemoglobin <- scDAPP::calculate_percent.hemoglobin(sobj)
  }

  if (!("Phase" %in% colnames(sobj@meta.data))) {
    try(
      sobj <- Seurat::CellCycleScoring(
        sobj,
        s.features = Seurat::cc.genes.updated.2019$s.genes,
        g2m.features = Seurat::cc.genes.updated.2019$g2m.features
      ),
      silent = TRUE
    )
  }

  sobj
}


.prefilter_cluster_seurat <- function(sobj, prefilter_pcs, prefilter_resolution) {
  dims <- seq_len(prefilter_pcs)
  suppressWarnings(
    sobj <- Seurat::SCTransform(sobj, verbose = TRUE, method = "glmGamPoi")
  )
  sobj <- Seurat::RunPCA(object = sobj, verbose = FALSE)
  sobj <- Seurat::FindNeighbors(object = sobj, dims = dims, verbose = FALSE)
  sobj <- Seurat::FindClusters(
    object = sobj,
    resolution = prefilter_resolution,
    verbose = FALSE,
    algorithm = 1
  )
  sobj <- Seurat::RunUMAP(sobj, dims = dims)
  sobj
}


.has_dim_reduction <- function(sobj) {
  length(Seurat::Reductions(sobj)) > 0L
}


#' Per-sample QC filtering module for the scDAPP pipeline
#'
#' Orchestrates optional pre-filter clustering, automated filtering via
#' [autofilter()], optional DoubletFinder, diagnostic plots, and saving
#' filtered and unfiltered objects.
#'
#' @param sobj Seurat object for one sample (with QC metadata)
#' @param code sample code used in output file names
#' @param min_num_UMI,min_num_Feature,max_perc_mito,max_perc_hemoglobin passed to [autofilter()]
#' @param autofilter_complexity,autofilter_mito,autofilter_nUMI logical; map to
#'   `globalfilter.complexity`, `globalfilter.mito`, and `globalfilter.libsize` in [autofilter()]
#' @param autofilter_medianabsolutedev_threshold,autofilter_loess_negative_residual_threshold
#'   passed to [autofilter()] as `mad.score.threshold` and `loess_negative_residual_threshold`
#' @param doubletFinder logical; run DoubletFinder on autofiltered cells
#' @param cluster_unfiltered logical; if TRUE, SCT+cluster the full matrix
#'   before autofilter (default FALSE)
#' @param qctmp_path path to save filtered object RDS (overwrites input tmp)
#' @param rawobj_path optional path to save unfiltered annotated object RDS; written
#'   only when `cluster_unfiltered` is TRUE and this is non-NULL
#' @param qc_pdf_path path for per-sample QC summary PDF
#' @param prefilter_pcs number of PCs for pre-filter / DoubletFinder clustering
#' @param prefilter_resolution Louvain resolution for pre-filter clustering
#'
#' @return list with components `af` (autofilter result plus plots) and
#'   `rawmd` (metadata from the unfiltered object with filter annotations)
#' @export
per_sample_qc_module <- function(
    sobj,
    code,
    min_num_UMI,
    min_num_Feature,
    max_perc_mito,
    max_perc_hemoglobin,
    autofilter_complexity,
    autofilter_mito,
    autofilter_nUMI,
    autofilter_medianabsolutedev_threshold,
    autofilter_loess_negative_residual_threshold,
    doubletFinder = TRUE,
    cluster_unfiltered = FALSE,
    qctmp_path,
    rawobj_path = NULL,
    qc_pdf_path,
    prefilter_pcs = 30L,
    prefilter_resolution = 0.1
) {
  message(code)

  if (isTRUE(cluster_unfiltered)) {
    sobj <- .prefilter_cluster_seurat(sobj, prefilter_pcs, prefilter_resolution)
  }

  af <- scDAPP::autofilter(
    sobj,
    min_num_UMI = min_num_UMI,
    min_num_Feature = min_num_Feature,
    max_perc_mito = max_perc_mito,
    max_perc_hemoglobin = max_perc_hemoglobin,
    globalfilter.complexity = autofilter_complexity,
    globalfilter.mito = autofilter_mito,
    globalfilter.libsize = autofilter_nUMI,
    mad.score.threshold = autofilter_medianabsolutedev_threshold,
    loess_negative_residual_threshold = autofilter_loess_negative_residual_threshold
  )

  sobjraw <- sobj
  cellstatus <- af$cellstatus
  goodcells <- cellstatus[cellstatus$filteredout == FALSE, "barcodes"]
  sobj <- sobj[, goodcells]

  if (isTRUE(doubletFinder)) {
    if (
      (utils::packageVersion("DoubletFinder") == "2.0.3") &&
        (utils::packageVersion("Seurat") >= "5.0.0") &&
        "SCT" %in% names(sobj@assays)
    ) {
      sobj_df <- Seurat::GetAssayData(sobj, assay = "SCT", layer = "data")
      sobj_df <- Seurat::CreateAssayObject(sobj_df)
      sobj_df <- Seurat::CreateSeuratObject(sobj_df)
      warning(
        "Will attempt to coerce v5 Seurat object to work with DoubletFinder v2.0.3; ",
        "this is unstable and does not always work! If any errors arise, set ",
        "doubletFinder to FALSE in pipeline runner",
        call. = FALSE
      )
    } else {
      sobj_df <- sobj
    }

    sobj_df <- .prefilter_cluster_seurat(
      sobj_df,
      prefilter_pcs,
      prefilter_resolution
    )

    try(
      expr = {
        af <- scDAPP::doubletfinderwrapper(
          sobj_df,
          autofilterres = af,
          num.cores = 1
        )
        cellstatus <- af$cellstatus
        goodcells <- cellstatus[cellstatus$filteredout == FALSE, "barcodes"]
        sobj <- sobj[, goodcells]
      },
      silent = FALSE
    )

    rm(sobj_df)
    gc(full = TRUE)
  }

  prefilter_cols <- grep("SCT_snn_res", colnames(sobj@meta.data), value = TRUE)
  if (length(prefilter_cols) > 0L) {
    colnames(sobj@meta.data)[
      grepl("SCT_snn_res.0.1", colnames(sobj@meta.data))
    ] <- "PREFILTER_SCT_snn_res.0.1"
  }

  sobjsave <- sobj
  rm(sobj)

  sobjraw@meta.data <- cbind(sobjraw@meta.data, af$cellstatus[, -1, drop = FALSE])

  has_clusters <- isTRUE(cluster_unfiltered) &&
    "seurat_clusters" %in% colnames(sobjraw@meta.data)
  has_embedding <- .has_dim_reduction(sobjraw)

  af$d_rawclust <- NULL
  af$d_raw_filt <- NULL
  af$d_raw_filt_reason <- NULL
  af$fp_raw_qc <- NULL
  af$tab_filt_by_clust <- NULL
  af$hm_raw <- NULL
  af$ap_filt <- NULL
  af$ap_filt_reason <- NULL
  af$twt_filt <- NULL
  af$twt_filt_reason <- NULL

  if (has_clusters) {
    m <- Seurat::FindAllMarkers(sobjraw, only.pos = TRUE)
    m$score <- (m$pct.1 - m$pct.2) * m$avg_log2FC

    if (length(levels(sobjraw$seurat_clusters)) > 1L) {
      n <- 5L
      top <- dplyr::slice_max(
        dplyr::group_by(m, .data$cluster),
        order_by = .data$score,
        n = n,
        with_ties = FALSE
      )
    } else {
      top <- NULL
    }

    af$d_rawclust <- Seurat::DimPlot(
      sobjraw,
      group.by = "seurat_clusters",
      label = TRUE,
      repel = TRUE
    ) +
      ggplot2::ggtitle(
        "Unfiltered data clusters",
        subtitle = paste0("Louvain res = ", prefilter_resolution)
      )

    if (has_embedding) {
      af$d_raw_filt <- Seurat::DimPlot(
        sobjraw,
        group.by = "filteredout",
        label = FALSE,
        repel = TRUE
      )
      sobjraw$filterreason <- factor(
        sobjraw$filterreason,
        levels = names(sort(table(sobjraw$filterreason), decreasing = TRUE))
      )
      af$d_raw_filt_reason <- Seurat::DimPlot(
        sobjraw,
        group.by = "filterreason",
        label = FALSE,
        repel = TRUE
      )
      af$fp_raw_qc <- Seurat::FeaturePlot(
        sobjraw,
        c("nCount_RNA", "nFeature_RNA", "percent.mito", "percent.hemoglobin"),
        order = TRUE
      )
    }

    tab_filt_by_clust <- table(sobjraw$filterreason, sobjraw$seurat_clusters)
    tab_filt_by_clust <- t(tab_filt_by_clust)
    rownames(tab_filt_by_clust) <- paste0("cluster_", rownames(tab_filt_by_clust))
    colnames(tab_filt_by_clust) <- gsub(
      x = colnames(tab_filt_by_clust),
      pattern = "\\.",
      replacement = "\n"
    )
    af$tab_filt_by_clust <- tab_filt_by_clust

    if (!is.null(top) && length(levels(sobjraw$seurat_clusters)) > 1L) {
      af$hm_raw <- Seurat::DoHeatmap(sobjraw, top$gene, raster = FALSE) +
        Seurat::NoLegend() +
        ggplot2::labs(title = "Pre-filter cluster markers")
    }

    sobjraw$filterreason <- factor(
      sobjraw$filterreason,
      levels = names(sort(table(sobjraw$filterreason), decreasing = TRUE))
    )
    pal <- grDevices::colorRampPalette(
      RColorBrewer::brewer.pal("Dark2", n = 8)
    )(length(levels(sobjraw$seurat_clusters)))

    af$ap_filt <- scDAPP::alluvialplot(
      sobjraw@meta.data[, c("seurat_clusters", "filteredout"), drop = FALSE]
    ) +
      ggplot2::scale_fill_manual(values = pal) +
      ggplot2::labs(title = "Cluster filtering")

    af$twt_filt <- scDAPP::twt_colored_heatmap(
      sobjraw@meta.data[, c("seurat_clusters", "filteredout"), drop = FALSE],
      title = "Cluster filtering"
    )

    af$ap_filt_reason <- scDAPP::alluvialplot(
      sobjraw@meta.data[, c("seurat_clusters", "filterreason"), drop = FALSE]
    ) +
      ggplot2::scale_fill_manual(values = pal) +
      ggplot2::labs(title = "Cluster filtering reason")

    af$twt_filt_reason <- scDAPP::twt_colored_heatmap(
      sobjraw@meta.data[, c("seurat_clusters", "filterreason"), drop = FALSE],
      title = "Cluster filtering reason"
    )
  } else if (has_embedding) {
    af$d_raw_filt <- Seurat::DimPlot(
      sobjraw,
      group.by = "filteredout",
      label = FALSE,
      repel = TRUE
    )
    sobjraw$filterreason <- factor(
      sobjraw$filterreason,
      levels = names(sort(table(sobjraw$filterreason), decreasing = TRUE))
    )
    af$d_raw_filt_reason <- Seurat::DimPlot(
      sobjraw,
      group.by = "filterreason",
      label = FALSE,
      repel = TRUE
    )
    af$fp_raw_qc <- Seurat::FeaturePlot(
      sobjraw,
      c("nCount_RNA", "nFeature_RNA", "percent.mito", "percent.hemoglobin"),
      order = TRUE
    )
  }

  comm <- af$allcommands
  rownames(comm) <- comm$Command

  af$vln_umi <- Seurat::VlnPlot(sobjraw, "nCount_RNA", group.by = "orig.ident") +
    ggplot2::scale_y_log10(labels = scales::label_comma()) +
    ggplot2::geom_hline(
      yintercept = comm["min_num_UMI", 2],
      linetype = "dotted"
    ) +
    ggplot2::labs(caption = paste0("cutoff = ", comm["min_num_UMI", 2]))

  af$vln_feature <- Seurat::VlnPlot(sobjraw, "nFeature_RNA", group.by = "orig.ident") +
    ggplot2::scale_y_log10(labels = scales::label_comma()) +
    ggplot2::geom_hline(
      yintercept = comm["min_num_Feature", 2],
      linetype = "dotted"
    ) +
    ggplot2::labs(caption = paste0("cutoff = ", comm["min_num_Feature", 2]))

  af$vln_mito <- Seurat::VlnPlot(sobjraw, "percent.mito", group.by = "orig.ident") +
    ggplot2::geom_hline(
      yintercept = comm["max_perc_mito", 2],
      linetype = "dotted"
    ) +
    ggplot2::labs(caption = paste0("cutoff = ", comm["max_perc_mito", 2]))

  af$vln_hemo <- Seurat::VlnPlot(
    sobjraw,
    "percent.hemoglobin",
    group.by = "orig.ident"
  ) +
    ggplot2::geom_hline(
      yintercept = comm["max_perc_hemoglobin", 2],
      linetype = "dotted"
    ) +
    ggplot2::labs(caption = paste0("cutoff = ", comm["max_perc_hemoglobin", 2]))

  colnames(af$baseline_qc_summary) <- gsub(
    "summary_",
    "summary\n",
    colnames(af$baseline_qc_summary)
  )

  if (isTRUE(cluster_unfiltered) && !is.null(rawobj_path)) {
    saveRDS(sobjraw, rawobj_path)
  }

  grDevices::pdf(qc_pdf_path, height = 10, width = 10)
  scDAPP::pdftable(af$filtersummary, title = "Cell Filtering Summary")
  scDAPP::pdftable(af$allcommands, title = "Filter parameters")
  scDAPP::pdftable(round(af$baseline_qc_summary, 2), title = "QC summary stats")
  print(af$vln_umi)
  print(af$vln_feature)
  print(af$vln_mito)
  print(af$vln_hemo)
  if (!is.null(af$globalfilter.complexity)) print(af$globalfilter.complexity)
  if (!is.null(af$globalfilter.libsize)) print(af$globalfilter.libsize)
  if (!is.null(af$globalfilter.mito)) print(af$globalfilter.mito)
  if (!is.null(af$d_rawclust)) print(af$d_rawclust)
  if (!is.null(af$d_raw_filt)) print(af$d_raw_filt)
  if (!is.null(af$d_raw_filt_reason)) print(af$d_raw_filt_reason)
  if (!is.null(af$fp_raw_qc)) print(af$fp_raw_qc)
  if (!is.null(af$tab_filt_by_clust)) {
    scDAPP::pdftable(af$tab_filt_by_clust, title = "Cell filtering per cluster")
  }
  if (!is.null(af$hm_raw)) print(af$hm_raw)
  if (!is.null(af$ap_filt)) print(af$ap_filt)
  if (!is.null(af$twt_filt)) print(af$twt_filt)
  if (!is.null(af$ap_filt_reason)) print(af$ap_filt_reason)
  if (!is.null(af$twt_filt_reason)) print(af$twt_filt_reason)
  grDevices::dev.off()

  saveRDS(sobjsave, qctmp_path)

  rawmd <- sobjraw@meta.data
  rm(sobjsave, sobjraw)
  invisible(gc(full = TRUE, reset = FALSE, verbose = FALSE))

  list(af = af, rawmd = rawmd)
}


