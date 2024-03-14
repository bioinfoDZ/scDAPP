# https://r-pkgs.org/whole-game.html




# UPDATE MARCH 14 2024 - IGRAPH UPDATE BROKE RISC LOUVAIN, 1.7 WAS UPDATED ON THIS DATE, ALSO INTRODUCED RESOLUTION TO DEFAULT RISC FUNCTION


#' 
#' #' Cluster RISC object with Louvain Resolution option
#' #'
#' #' This function is almost identical to `RISC::scCluster()` function with method = "louvain", with the addition of louvain resolution as a parameter. Higher resolution = more clusters. The R package `igraph` from CRAN must be installed and up to date, tested with version 1.4.3. Additionally, clusters will now be sorted biggest to smallest by number of cells (1 is biggest). The Cluster variable added to the object coldata will include resolution in the column name. Density clustering based on low dimensional embedding (ie UMAP) is not allowed in this function, only Louvain is supported, just use the original RISC verison for that, though it is not recommended.
#' #'
#' #' @param object RISC object: a framework dataset.
#' #' @param slot The dimension_reduction slot for cell clustering. The default is "cell.umap" under RISC object "DimReduction" item for UMAP method, but the customer can add new dimension_reduction method under DimReduction and use it.
#' #' @param neighbor The neighbor cells for "igraph" method.
#' #' @param resolution Optional resolution parameter that allows the user to adjust the resolution parameter of the modularity function that the algorithm uses internally. Lower values typically yield fewer, larger clusters. The original definition of modularity is recovered when the resolution parameter is set to 1.
#' #' @param algorithm The algorithm for knn, the default is "kd_tree", all options: "kd_tree", "cover_tree", "CR", "brute".
#' #' @param npc The number of PCA or PLS used for cell clustering.
#' #' @param redo Whether re-cluster the cells.
#' #' @param random.seed The random seed, the default is 123
#' #'
#' #' @return RISC single cell dataset, the cluster slot; the column in coldata will be named as "Cluster_res.(resolution)"
#' #' @export
#' #'
#' #' @examples
#' #' \dontrun{
#' #' obj0 = scCluster(obj0, slot = 'cell.pls', resolution = 0.5)
#' #' DimPlot(obj0, slot = "cell.umap", colFactor = 'Cluster_res.0.5', size = 2)
#' #'}
#' scCluster_louvain_res <- function (object,
#'                                    slot = "cell.umap",
#'                                    neighbor = 10,
#'                                    resolution = 1,
#'                                    algorithm = "kd_tree",
#'                                    npc = 20,
#'                                    redo = TRUE,
#'                                    random.seed = 123
#'                                    )
#' {
#' 
#' 
#'   # require(FNN) -> get.knn function?
#'   # require(igraph) --> graph_from_data_frame , simplify, and cluster_louvain
#' 
#'   ## hard code use louvain not density
#'   # turn off k, remove ref to "dc", set method = 'louvain'
#'   set.seed(random.seed)
#'   # k = as.integer(k)
#'   neighbor = as.integer(neighbor)
#'   algorithm = as.character(algorithm)
#'   npc = as.integer(npc)
#'   random.seed = as.integer(random.seed)
#'   method = 'louvain'
#'   resolution = as.numeric(resolution)
#' 
#' 
#'   message("Clustering RISC object with:",
#'           "\n - slot = ", slot,
#'           "\n - neighbor = ", neighbor,
#'           "\n - resolution = ", resolution,
#'           "\n - knn algorithm = ", algorithm,
#'           "\n - npc = ", npc,
#'           "\n - random.seed = ", random.seed
#'           )
#' 
#'   # if (is.null(dc)) {
#'   #   dc = object@metadata$dcluster$dc
#'   #   if (isTRUE(redo)) {
#'   #     dc = NULL
#'   #   }
#'   #   else {
#'   #     dc = dc
#'   #   }
#'   # }
#'   # else {
#'   #   dc = dc
#'   # }
#' 
#'   if (!is.null(object@vargene) & !is.null(object@DimReduction)) {
#' 
#' 
#'     slot0 = as.character(slot)
#'     dimReduce0 = object@DimReduction[[slot0]]
#'     if (is.null(dimReduce0)) {
#'       stop("Do not include this dimention_reduction slot, try another one")
#'     }
#'     else if (ncol(dimReduce0) >= npc) {
#'       count = as.matrix(dimReduce0[, 1:npc])
#'     }
#'     else {
#'       count = as.matrix(dimReduce0)
#'     }
#'     # if (method == "density") {
#'     #   dist0 = dist(count)
#'     #   if (is.null(dc)) {
#'     #     dataClust = densityClust(dist0, gaussian = TRUE)
#'     #   }
#'     #   else {
#'     #     dataClust = densityClust(dist0, dc = dc, gaussian = TRUE)
#'     #   }
#'     #   delta.rho = data.frame(rho = dataClust$rho, delta = dataClust$delta,
#'     #                          stringsAsFactors = FALSE)
#'     #   delta.rho = delta.rho[order(delta.rho$delta, decreasing = TRUE),
#'     #   ]
#'     #   delta.cut = delta.rho$delta[k + 1L]
#'     #   clust0 = findClusters(dataClust, 0, delta.cut)
#'     #   object@cluster = object@coldata$Cluster = as.factor(clust0$clusters)
#'     #   names(object@cluster) = names(object@coldata$Cluster) = rownames(count)
#'     #   object@metadata[["clustering"]] = data.frame(Method = "densityClust",
#'     #                                                Distance = as.numeric(clust0$dc), rho = as.numeric(clust0$threshold[1]),
#'     #                                                delta = as.numeric(clust0$threshold[2]), stringsAsFactors = F)
#'     # }
#'     #
#'     #
#'     # else
#' 
#' 
#' 
#'     if (method == "louvain") {
#' 
#'       message('\nFinding neighbors')
#'       clust0 = FNN::get.knn(count, k = neighbor, algorithm = algorithm)
#' 
#'       message('\nClustering')
#'       clust1 = data.frame(NodStar = rep(1L:nrow(count),
#'                                         neighbor), NodEnd = as.vector(clust0$nn.index),
#'                           stringsAsFactors = FALSE)
#'       clust1 = igraph::graph_from_data_frame(clust1, directed = FALSE)
#'       igraph::E(clust1)$weight = 1/(1 + as.vector(clust0$nn.dist)) # ADDED IN RISC 1.7
#'       clust1 = igraph::simplify(clust1)
#'       clust1 = igraph::cluster_louvain(clust1, weights = 1/(1 +
#'                                                       as.vector(clust0$nn.dist)),
#'                                resolution = resolution
#'                                )
#' 
#'       ### remap clust, biggest to smallest...
#'       #get risc clust
#'       rc <- as.factor(clust1$membership)
#' 
#'       #sort by biggest to smallest
#'       bs <- sort(table(rc), decreasing = T)
#'       rc <- plyr::mapvalues(rc, from = names(bs), to = c(1:length(bs)) )
#'       rc <- factor(rc, levels = 1:length(bs))
#' 
#' 
#' 
#'       # in metadata name, name it as Cluster_res.[resolution]
#'       clustvarname <- paste0("Cluster_res.", resolution)
#' 
#'       #add size-sorted clusters to object
#'       object@coldata$tmpclust = rc
#' 
#'       #rename tmpclust to clustvarname
#'       colnames(object@coldata)[ncol(object@coldata)] <- clustvarname
#' 
#' 
#'       object@cluster = object@coldata[,clustvarname]
#'       names(object@cluster) = names(object@coldata[,clustvarname]) = rownames(count)
#'       object@metadata[["clustering"]] = data.frame(Method = "louvain",
#'                                                    PCs = npc, Neighbors = neighbor, louvain_resolution = resolution,
#'                                                    stringsAsFactors = F)
#'     }
#'     else {
#'       stop("A new method later")
#'     }
#'   }
#'   else {
#'     stop("Please calculate dispersion and dimention reduction first")
#'   }
#'   return(object)
#' }
#' 
#' 



### test ####


# object <- readRDS("~/Dropbox/Result_from_Alex/deyoudata/hto_prostate/outs/NEWDATA_scRNAseqpipeline/multisample_integration/data_objects/RISC-object_integrated.rds")
#
#
# slot = "cell.pls"
# neighbor = 10
# algorithm = "kd_tree"
# method = "louvain"
# npc = 30
# k = 10
# dc = NULL
# redo = TRUE
# random.seed = 123
# resolution = 0.5
#
# rm(k,dc,method)
#
#
# object <- scCluster_louvain_res(object, slot = 'cell.pls', resolution = 0.5, npc = 30)
# RISC::DimPlot(object, colFactor = "Cluster_res.0.5")






