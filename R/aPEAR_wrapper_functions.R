#' Convert scDAPP Pathway Analysis Results to aPEAR Input Format
#'
#' Converts gene set enrichment analysis (GSEA) results from scDAPP across 
#' different conditions into a format compatible with aPEAR. It organizes data 
#' from `pathway_analysis_mainlist_comp` by condition and cluster, arranging 
#' MSIGDB categories (e.g., Hallmark) into a designated column. This function 
#' requires `tidyverse` for data manipulation and string operations.
#'
#' @param object An object from scDAPP pathway analysis, expected to be an RDS 
#' file located at `multisample_integration/pathwayanalysis_crosscondition/DE_pathways_plot_objects_list.rds`. 
#' The object must be a list containing pathway analysis results structured 
#' specifically for compatibility with this function.
#' @param show_Category Logical; if TRUE (default), includes the MSIGDB category in 
#' the cluster description. This affects the transformation of pathway names based 
#' on the presence of certain prefixes (e.g., "HALLMARK_").
#' @param ... Additional arguments passed to the underlying clustering functions.
#'
#' @return Returns a nested list where each key represents a comparison condition 
#' (e.g., [condition A]_vs_[condition B]). Each corresponding value is a list of 
#' data frames, one for each unique cell cluster identified in the analysis. Data 
#' frames contain columns: "Description", "Category", "setSize", "enrichmentScore", 
#' "NES", "pvalue", "p.adjust", "core_enrichment", and "Des_Full", detailing the 
#' cluster's description, MSIGDB category, set size, enrichment score, normalized 
#' enrichment score (NES), p-value, adjusted p-value, core enrichment genes, and 
#' full description, respectively.
#'
#' @examples
#' # Assuming `pathway_analysis_result` is your scDAPP analysis result object:
#' apear_input_list <- apear_data_prep(pathway_analysis_result)
#'
#' @details The function performs several key operations: validating the input object's 
#' structure, creating a mapping of pathways to clusters, and transforming the pathway 
#' analysis results into the aPEAR input format. Error handling is robust, with 
#' specific messages provided for common issues such as missing or incorrectly 
#' structured input data. It is crucial that the input object adheres to the expected 
#' structure for successful transformation.
#'
#' @note This function assumes that `tidyverse` is installed and available 
#' for use. If not, please install these packages using `install.packages()`.
#' 
#' @export
apear_data_prep <- function(object,
                            show_Category=TRUE,
                            ...
){
  ######################check the object####################################### 
  # Check if object is provided and has the expected structure
  if (is.null(object)) {
    stop("An pathway list object is not provided. Please provide an pathway list object .")
  }
  
  if (!is.list(object) || is.null(object$pathway_analysis_mainlist_comps) || length(object$pathway_analysis_mainlist_comps) == 0) {
    stop("An pathway list object does not have the expected structure or is empty. Please check!")
  }
  ############################################################################ 
  
  
  
  
  ##############################Create a map###################################
  #create a map/dataframe to keep track of different layers of the object
  map_list <- vector("list", length = length(object$pathway_analysis_mainlist_comps))
  
  # loop through comparison conditions e.g. "YoungKO_vs_YoungWT"
  for(x in 1:length(object$pathway_analysis_mainlist_comps)){
    tmp_comps <- object$pathway_analysis_mainlist_comps[[x]]
    
    tmp_map <- data.frame() # create a tmp dataframe that will summarize the pathway layers and cluster sub-layers
    
    # loop through pathways e.g. "HALLMARK" and make a summary of pathway names and cluster names
    for (y in 1:length(tmp_comps)) {
      tmp_pathway <- data.frame(pathway = rep(names(tmp_comps)[y],length(tmp_comps[[y]])),
                                cluster = names(tmp_comps[[y]]))
      tmp_map <- rbind(tmp_map,tmp_pathway)
    }
    map_list[[x]] <- tmp_map # add the dataframe to a list for each comparison
  }
  names(map_list) <- names(object$pathway_analysis_mainlist_comps) # give the dataframe the same name as the comparison
  
  unique_clusters <- lapply(map_list, function(df) unique(df$cluster)) # go through the list to identify unique clusters
  
  
  
  # Create a function to add a column with the name of the list element
  add_name_column <- function(df, name) {
    df <- df %>%
      mutate(category = name)
    return(df)
  }
  ############################################################################
  
  
  
  
  #################Convert the data for apear input##########################
  
  # create an empty list for aPEAR inputs
  apear_input_list <- vector("list", length = length(map_list))
  
  
  # Auxiliary function to clean pathway names based on the show_Category parameter
  clean_pathway_names <- function(pathway) {
    string_removal <- c("HALLMARK_", "GOCC_", "GOBP_", "GOMF_", "KEGG_", "REACTOME_")
    if (show_Category) {
      return(str_replace_all(pathway, "_", " "))
    } else {
      return(pathway %>% 
               str_replace_all(paste(string_removal, collapse = "|"), "") %>%
               str_replace_all("_", " "))
    }
  }
  
  # Auxiliary function to prepare aPEAR input from result_combined
  prepare_apear_input <- function(result_combined) {
    result_combined %>%
      mutate(Description = clean_pathway_names(pathway),
             Des_Full = pathway,
             Category = category,
             setSize = size,
             enrichmentScore = ES,
             NES = NES,
             pvalue = pval,
             p.adjust = padj,
             core_enrichment = sapply(leadingEdge, function(les) paste(les, collapse = "/"))) %>%
      select(Description, Category, setSize, enrichmentScore, NES, pvalue, p.adjust, core_enrichment, Des_Full)
  }
  
  
  # loop through the map to pullout the information
  for (x in 1:length(map_list)) {
    tmp_map <-  map_list[[x]] # pull out pathway and cluster info from each comparison 
    tmp_uniq_cluster <- unique_clusters[[x]] #pull out unique cluster info for each comparison 
    
    #create a tmp list for the result after we integrate the info from above
    tmp_result_list <- vector("list", length = length(tmp_uniq_cluster))
    
    #loop through each unique cluster
    for (y in 1:length(tmp_uniq_cluster)) {
      tmp_pathway <- tmp_map$pathway[tmp_map$cluster %in% tmp_uniq_cluster[y]] #pull out only the pathways that contain each unique cluster
      
      #for each of these pathways extract their gsea results 
      result_list <- lapply(tmp_pathway, function(pathway) {
        object$pathway_analysis_mainlist_comps[[x]][[pathway]][[tmp_uniq_cluster[y]]][["gseares"]]
      }) 
      names(result_list) <- tmp_pathway
      
      #Apply the function to each element in the list
      result_list <- Map(add_name_column, result_list, names(result_list))
      
      
      #combine their gsea results into a dataframe 
      result_combined <- do.call(rbind, result_list)
      rownames(result_combined) <- NULL
      
      # transform the result_combined dataframe into apear acceptable format
      aPEAR_input<- prepare_apear_input(result_combined)
      
      dup_description <- unique(aPEAR_input$Description[duplicated(aPEAR_input$Description)]) 
      aPEAR_input[aPEAR_input$Description %in% dup_description,]$Description <- aPEAR_input[aPEAR_input$Description %in% dup_description,]$Des_Full %>% gsub("_"," ",.) 
      
      
      #save the aPEAR_input to the a tmp_result_list
      tmp_result_list[[y]] <- aPEAR_input
      
    } 
    
    #name the tmp_result_list with same unique cluster 
    names(tmp_result_list) <- tmp_uniq_cluster
    
    #store tmp_result_list into apear_input_list
    apear_input_list[[x]] <- tmp_result_list
  }
  
  #name the tmp_result_list with same comparison
  names(apear_input_list) <- names(map_list)
  
  ############################################################################
  
  
  
  return(apear_input_list)
}







#' Select Pathways and Calculate aPEAR Clusters
#'
#' This function serves as a wrapper for the `findPathCluster` function in the aPEAR package. 
#' It processes the output from the `apear_data_prep` function, filters pathways based on the number 
#' of leading edge genes, and selects pathways with the most significant p-values. It employs dynamic 
#' clustering, primarily using the Markov clustering method, adjusting the minimum cluster size as needed 
#' to optimize the number of labels displayed on the final graph. If Markov clustering fails, it 
#' automatically attempts hierarchical clustering as a fallback.
#'
#' @param apear_input_list The output list from the `apear_data_prep` function, structured by comparison groups.
#' @param initial_minCS The initial minimum cluster size for `aPEAR::findPathCluster`. Defaults to 4.
#' @param max_label_threshold The maximum number of labels allowed on the final graph. Defaults to 50.
#' @param max_sign_pathway The maximum number of pathways considered, based on p-value significance. Defaults to 500.
#' @param min_leading_edge_threshold The minimum number of leading edge genes required for a pathway. Defaults to 2.
#' @param output_dir The directory where the `findPathCluster` results will be saved. `NULL` by default.
#' @param save_cluster_parameters_and_results Logical; if `TRUE`, the cluster parameters and results are saved. Defaults to TRUE.
#' @param ... Additional arguments passed to the underlying clustering functions.
#'
#' @return Returns a list of two sub-lists named `clusterres` and `params`. `clusterres` contains nested list with the comparison groups (e.g., [condition A]_vs_[condition B]). 
#' Each sub-list within `clusterres` contains two nested lists: `findPathCluster_opt_inputs`, which includes the `clustMethod` and `minCS` parameters, 
#' and `findPathClusterres`, which contains cluster information and similarity values. `params` contains `max_sign_pathway` and `min_leading_edge_threshold` information.
#'
#' @examples
#' # Assuming `apear_data_prep_output` is your prepared input list:
#' findPathClusterres <- apear_find_clusters(apear_input_list,
#'                                           max_sign_pathway = 500,
#'                                           min_leading_edge_threshold = 2,
#'                                           output_dir = "path/to/save",
#'                                           save_cluster_parameters_and_results=TRUE)
#'
#' @note Ensure that the `output_dir` is specified if `save_cluster_parameters_and_results` is set to TRUE. 
#' The parameter `save_cluster_parameters_and_results` is strongly encouraged to be set as `TRUE` to ensure the reproducibility.
#' If the clusters for the final graph are too crowded, you may consider to decrease `max_sign_pathway`. 
#'
#' @export
apear_find_clusters <- function(apear_input_list,
                                initial_minCS=4,
                                max_label_threshold = 50,
                                max_sign_pathway= 500,
                                min_leading_edge_threshold = 2,
                                output_dir = NULL,
                                save_cluster_parameters_and_results=TRUE, 
                                ...
                                
){
  
  
  ############### Check before find path calculation ########################### 
  if(save_cluster_parameters_and_results) {
    # Check if output_dir is NULL when save_cluster_parameters_and_results is TRUE
    if(is.null(output_dir)) {
      stop("save_cluster_parameters_and_results is TRUE but output_dir is not set. Please provide a valid output directory before proceeding.")
    } else{
      # Ensure the output directory exists or create it
      if (!dir.exists(output_dir)) {
        message(paste0("Output directory:",output_dir, " does NOT exist! Will create one!") )
        dir.create(output_dir, recursive = TRUE)
      }
    } 
  } else {
    # Handle case for save_cluster_parameters_and_results being FALSE
    # You can add logic here if there's anything specific to do when not saving results
    if(!is.null(output_dir)) {
      message("Note: Output directory is provided but save_cluster_parameters_and_results is FALSE, so no results will be saved.")
    }
    # Else, nothing needs to be saved, so you can optionally include logic for handling this scenario.
  }
  ############################################################################ 
  
  
  
  
  ########################Count Leading Edge genes########################
  #for each comparison in apear_input_list, count Leading Edge genes in core enrichment. number of "/" +1
  for (x in 1:length(apear_input_list)) {
    tmp_apear_input <- apear_input_list[[x]]
    for (y in 1:length(tmp_apear_input)) {
      tmp_str <- apear_input_list[[x]][[y]]$core_enrichment
      apear_input_list[[x]][[y]]$number_of_LE_genes <- sapply(tmp_str, function(s) sum(gregexpr("/", s)[[1]] >= 0)+1)
      apear_input_list[[x]][[y]] <- as.data.frame(apear_input_list[[x]][[y]])
    }
    
  }  
  
  ############################################################################
  
  
  
  
  
  
  ######################Create an output result list##########################
  
  #define unique_clusters
  unique_clusters <- lapply(apear_input_list, function(df) names(df))
  
  # create a list to store results for findPathCluster results of each cluster
  findPathClusterres <- vector("list",length = length(apear_input_list))
  names(findPathClusterres) <- names(apear_input_list)
  
  
  # loop through each comparison 
  for (x in 1:length(apear_input_list)) {
    # create a tmp cluster to work with in this iteration
    tmp_cluster <- unique_clusters[[names(apear_input_list)[x]]]
    
    # create a list to store result for each cluster
    findPathClusterres[[x]] <- vector("list",length = length(tmp_cluster))
    names(findPathClusterres[[x]]) <- tmp_cluster
    
    
    # loop through each cluster
    for (y in 1:length(tmp_cluster)) {
      
      findPathClusterres[[x]][[y]] <- vector("list",length = 2)
      names(findPathClusterres[[x]][[y]])<- c("findPathCluster_opt_inputs","findPathClusterres") #store both parameters and results
      
    }
  }
  ############################################################################
  
  
  
  
  
  
  ######################Calculate findPathCluster##############################
  # Loop through each comparison in apear_input_list
  for (x in 1:length(apear_input_list)) {
    tmp_apear_input <- apear_input_list[[x]]
    
    # Loop through each tmp_apear_input and pull out 1 cluster at a time
    for (y in 1:length(tmp_apear_input)) {
      tmp_cluster <- tmp_apear_input[[y]]
      
      # Filter out pathways with too few leading edge genes that cause troubles for markov clustering
      tmp_cluster <- tmp_cluster[which(tmp_cluster$number_of_LE_genes>min_leading_edge_threshold),]
      
      #filter number of pathways with too many pvalues
      if(nrow(tmp_cluster) > max_sign_pathway) {
        tmp_cluster <- tmp_cluster %>% arrange(pvalue)
        tmp_cluster <- tmp_cluster %>% top_n(., max_sign_pathway, -pvalue)
      }
      
      msg_name <- paste0(names(apear_input_list)[x], "(", names(tmp_apear_input)[y], ")")
      
      
      # Function to create enrichment network; w/ dynamic minCS cutoff
      cal_minCS <- function(clustMethod) {
        tmp_minCS <- initial_minCS
        while (TRUE) {
          apear_clusters <- aPEAR::findPathClusters(tmp_cluster,cluster = clustMethod,minClusterSize = tmp_minCS) 
          num_clusters <- length(unique(apear_clusters$clusters$Cluster))
          if (num_clusters < max_label_threshold) {
            print(paste0(msg_name, ": ",clustMethod," clustering Done! ","minClusterSize=",tmp_minCS))
            return(tmp_minCS) 
          } else {
            tmp_minCS <- tmp_minCS + 1
          }
        }
      }
      
      # In the original dataframe indicate whether a pathway is used in the findPathClusters; In other words, it passes the min_leading_edge_threshold and max_sign_pathway
      tmp_df <- apear_input_list[[x]][[y]]
      Select_for_aPEAR_vec <- rep("No", times = nrow(tmp_df))
      Select_for_aPEAR_vec[tmp_df$Des_Full %in% tmp_cluster$Des_Full] <- "Yes"
      apear_input_list[[x]][[y]]$Select_for_aPEAR <- Select_for_aPEAR_vec
      
      
      # Attempt 1: Markov clustering
      result_markov <- tryCatch({
        cal_minCS("markov")
      }, error = function(e1) {
        cat(msg_name, ":", conditionMessage(e1), "Will try hierarchical clustering!\n")
        NULL  # Return NULL to trigger the second attempt
      })
      
      # Attempt 2: Hierarchical clustering
      if (is.null(result_markov)) {
        result_hier <- tryCatch({
          cal_minCS("hier")
        }, error = function(e2) {
          cat(msg_name, ":", conditionMessage(e2), msg_name, "cannot be processed via aPEAR!\n")
          NULL  # Return NULL if hierarchical clustering also fails
        })
        
        if (!is.null(result_hier)) {
          # store both parameters and results 
          findPathClusterres[[x]][[y]][["findPathCluster_opt_inputs"]] <- list("hier",result_hier) 
          names(findPathClusterres[[x]][[y]][["findPathCluster_opt_inputs"]]) <- c("clustMethod","minCS")
          findPathClusterres[[x]][[y]][["findPathClusterres"]] <- aPEAR::findPathClusters(tmp_cluster,cluster = "hier",minClusterSize = result_hier)
        }
      } else {
        # store both parameters and results
        findPathClusterres[[x]][[y]][["findPathCluster_opt_inputs"]] <- list("markov",result_markov)
        names(findPathClusterres[[x]][[y]][["findPathCluster_opt_inputs"]]) <- c("clustMethod","minCS")
        findPathClusterres[[x]][[y]][["findPathClusterres"]] <- aPEAR::findPathClusters(tmp_cluster,cluster = "markov",minClusterSize = result_markov)
      }
      
    }
  }
  ############################################################################
  
  
  ######################Packing findPathCluster and parameters##################
  # Package the results and the used parameters
  findPathClusterres <- list(
    clusterres = findPathClusterres,  
    params = list(
      max_sign_pathway = max_sign_pathway,
      min_leading_edge_threshold = min_leading_edge_threshold
    )
  )
  ############################################################################
  
  
  ######################Save findPathCluster##################################
  #save findPathCluster results and parameters and results to avoid recalculation
  if(save_cluster_parameters_and_results){
    if(!is.null(output_dir)){
      saveRDS(findPathClusterres,paste0(output_dir,"findPathClusterres.rds"))
    }
  }
  ############################################################################
  
  return(findPathClusterres)
}






#' Update and Save aPEAR Input Data
#'
#' This function integrates the outputs from `apear_data_prep` and `apear_find_clusters` functions.
#' It updates the `apear_input_list` with clustering information and saves the result if specified.
#'
#' @param apear_input_list A list of sub-lists, which is the output from the `apear_data_prep` function.
#' @param findPathClusterres A list of sub-lists, which is the output from the `apear_find_clusters` function.
#' @param output_dir Directory where the updated `apear_input_list` and its corresponding CSV files are saved.`NULL` by default.
#' @param save_apear_input A boolean that controls whether the updated apear input list is saved. If `TRUE`,
#' directories are created for each comparison group ([condition A]_vs_[condition B]) with a subfolder "inputs" containing CSV files
#' for each cell cluster dataframe, which now include "number_of_LE_genes", "Select_for_aPEAR", and "aPEAR_Cluster" columns.
#' @param ... Additional arguments passed to the underlying clustering functions.
#'
#' @return An updated `apear_input_list` containing sub-lists named by comparison groups ([condition A]_vs_[condition B]). 
#' Each sub-list includes dataframes corresponding to unique cell clusters, now with 12 columns: "Description",
#' "Category", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "core_enrichment", "Des_Full",
#' "number_of_LE_genes", "Select_for_aPEAR", and "aPEAR_Cluster".
#'
#' @examples
#' # Assuming `apear_input_list` and `findPathClusterres` have been defined:
#' apear_input_list_updated <- apear_update_and_save_input(apear_input_list,
#'                                                         findPathClusterres,
#'                                                         output_dir = "path/to/save",
#'                                                         save_apear_input = TRUE)
#' 
#' @export
apear_update_and_save_input <- function(apear_input_list,
                                        findPathClusterres,
                                        output_dir = NULL,
                                        save_apear_input = TRUE,
                                        ...){
  
  ############### Check before find path calculation ########################### 
  if(save_apear_input) {
    # Check if output_dir is NULL when save_apear_input is TRUE
    if(is.null(output_dir)) {
      stop("save_apear_input is TRUE but output_dir is not set. Please provide a valid output directory before proceeding.")
    }else{
      # Ensure the output directory exists or create it
      if (!dir.exists(output_dir)) {
        message(paste0("Output directory:",output_dir, " does NOT exist! Will create one!") )
        dir.create(output_dir, recursive = TRUE)
      }
    }  
  } else {
    # Handle case for save_cluster_parameters_and_results being FALSE
    # You can add logic here if there's anything specific to do when not saving results
    if(!is.null(output_dir)) {
      message("Note: Output directory is provided but save_cluster_parameters_and_results is FALSE, so no results will be saved.")
    }
    # Else, nothing needs to be saved, so you can optionally include logic for handling this scenario.
  }
  ############################################################################ 
  
  
  
  
  ########################Define and create output dir########################
  #for each comparison in apear_input_list, if the dir does not exist, create one
  if(save_apear_input) {
    for(x in 1:length(apear_input_list)){
      tmp_dir <- paste0(output_dir,names(apear_input_list)[x],"/")
      if(!dir.exists(paste0(tmp_dir,"inputs/"))){
        dir.create(paste0(tmp_dir,"inputs/"),recursive = T)
      }
    }
  }
  ############################################################################
  
  
  ######################Unpacking findPathCluster and parameters##################
  # Extract parameters from the findPathClusterres object
  max_sign_pathway <- findPathClusterres$params$max_sign_pathway
  min_leading_edge_threshold <- findPathClusterres$params$min_leading_edge_threshold
  
  findPathClusterres <- findPathClusterres$clusterres
  ############################################################################
  
  
  
  
  ########################Count Leading Edge genes########################
  #for each comparison in apear_input_list, count Leading Edge genes in core enrichment. number of "/" +1
  for (x in 1:length(apear_input_list)) {
    tmp_apear_input <- apear_input_list[[x]]
    for (y in 1:length(tmp_apear_input)) {
      tmp_str <- apear_input_list[[x]][[y]]$core_enrichment
      apear_input_list[[x]][[y]]$number_of_LE_genes <- sapply(tmp_str, function(s) sum(gregexpr("/", s)[[1]] >= 0)+1)
      apear_input_list[[x]][[y]] <- as.data.frame(apear_input_list[[x]][[y]])
    }
    
  }  
  
  ############################################################################
  
  
  ######################Find selected pathway##############################
  # Loop through each comparison in apear_input_list
  for (x in 1:length(apear_input_list)) {
    tmp_apear_input <- apear_input_list[[x]]
    
    # Loop through each tmp_apear_input and pull out 1 cluster at a time
    for (y in 1:length(tmp_apear_input)) {
      tmp_cluster <- tmp_apear_input[[y]]
      
      # Filter out pathways with too few leading edge genes that cause troubles for markov clustering
      tmp_cluster <- tmp_cluster[which(tmp_cluster$number_of_LE_genes>min_leading_edge_threshold),]
      
      #filter number of pathways with too many pvalues
      if(nrow(tmp_cluster) > max_sign_pathway) {
        tmp_cluster <- tmp_cluster %>% arrange(pvalue)
        tmp_cluster <- tmp_cluster %>% top_n(., max_sign_pathway, -pvalue)
      }
      
      # In the original dataframe indicate whether a pathway is used in the findPathClusters; In other words, it passes the min_leading_edge_threshold and max_sign_pathway
      tmp_df <- apear_input_list[[x]][[y]]
      Select_for_aPEAR_vec <- rep("No", times = nrow(tmp_df))
      Select_for_aPEAR_vec[tmp_df$Des_Full %in% tmp_cluster$Des_Full] <- "Yes"
      apear_input_list[[x]][[y]]$Select_for_aPEAR <- Select_for_aPEAR_vec
      
    }
  }
  
  ############################################################################
  
  
  
  
  #######################Update apear_input obj #############################
  # Loop through each comparison in apear_input_list
  for (x in 1:length(apear_input_list)) {
    tmp_apear_input <- apear_input_list[[x]]
    
    # Loop through each tmp_apear_input and pull out 1 cluster at a time
    for (y in 1:length(tmp_apear_input)) {
      tmp_cluster <- tmp_apear_input[[y]]
      tmp_cluster$Pathway <- tmp_cluster$Description
      tmp_clusterres <- as.data.frame(findPathClusterres[[x]][[y]][["findPathClusterres"]][["clusters"]])
      
      # combine findPathClusterres and apear_input dataframes based on their pathway to pass the cluster info
      merged_df <- left_join(tmp_cluster, tmp_clusterres, by = "Pathway")
      merged_df$Pathway <- NULL
      merged_df$Cluster[is.na(merged_df$Cluster)] <- c("NA")
      
      #store the final cluster results 
      apear_input_list[[x]][[y]]$aPEAR_Cluster <- merged_df$Cluster
    }
  }
  ############################################################################
  
  
  
  
  #######################Save the outputs#####################################
  if(save_apear_input) {
    if(!is.null(output_dir)){
      # Loop through each comparison in apear_input_list
      for (x in 1:length(apear_input_list)) {
        tmp_output_dir <- paste0(output_dir, names(apear_input_list)[x], "/")
        tmp_apear_input <- apear_input_list[[x]]
        # Loop through each tmp_apear_input and pull out 1 cluster at a time
        for (y in 1:length(tmp_apear_input)) {
          
          apear_input_tmp <-  apear_input_list[[x]][[y]]
          #save the result inputs as csv for future usage 
          write.csv(apear_input_tmp, paste0(tmp_output_dir,"inputs/",names(apear_input_list)[x],"_",names(apear_input_list[[x]])[y],"_aPEAR_input.csv"),row.names = F,quote = F)
        }
      }
    }else{
      message("Please provide an output directory for the output_dir paramater")
    }  
  } 
  
  ############################################################################
  
  
  return(apear_input_list)
}






#' Plot aPEAR Clusters
#'
#' Visualizes clusters generated by aPEAR analysis using the updated aPEAR input list and the
#' clustering results. The function allows customization of the plot appearance and the option
#' to save the plot and aPEAR results object.
#'
#' @param apear_input_list The updated list output from the `apear_update_and_save_input` function.
#' @param findPathClusterres The list output from the `apear_find_clusters` function containing
#' clustering results.
#' @param color_mapping A named vector of colors to be assigned to comparison groups for plotting. `NULL` by default.
#' @param repelLabels Logical; if `TRUE`, labels are adjusted to minimize overlap in the plot. Defaults to `FALSE`.
#' @param drawEllipses Logical; if `TRUE`, ellipses are drawn around the clusters. Defaults to `FALSE`.
#' @param fontSize Numeric; the font size used for cluster labels in the plot. Defaults to 2.5.
#' @param show_parameters Logical; if `TRUE`, plot includes clustering parameters. Defaults to `TRUE`. 
#' @param output_dir A string specifying the directory where output files will be saved. Defaults to `TRUE`.
#' @param output_width Numeric; the width of the output plot. Defaults to 10.
#' @param output_height Numeric; the height of the output plot. Defaults to 10.
#' @param save_network_plot Logical; if `TRUE`, the network plot is saved to the specified `output_dir`. Defaults to `TRUE`.
#' @param save_apearres_obj Logical; if `TRUE`, the aPEAR results object is saved to the specified `output_dir`. Defaults to `TRUE`.
#' @param ... Additional arguments passed to the underlying clustering functions.
#'
#' @return Generates and optionally saves a visual representation of the aPEAR clusters and the aPEAR
#' results object.
#'
#' @examples
#' # Assuming `apear_input_list` and `findPathClusterres` have been defined, and color mapping is set:
#' apear_network_list <-  apear_plot_cluster(apear_input_list_updated,
#'                                           findPathClusterres,
#'                                           color_mapping = c(Group1 = "blue", Group2 = "red"),
#'                                           repelLabels = FALSE,
#'                                           drawEllipses = FALSE,
#'                                           fontSize = 2.5,
#'                                           show_parameters = TRUE,
#'                                           output_dir = "path/to/save",
#'                                           output_width = 10,
#'                                           output_height = 10,
#'                                           save_network_plot = TRUE,
#'                                           save_apearres_obj = TRUE)
#'
#' @export
apear_plot_cluster <- function(apear_input_list,
                               findPathClusterres,
                               color_mapping = NULL,
                               repelLabels = FALSE,
                               drawEllipses = FALSE,
                               fontSize = 2.5,
                               show_parameters=TRUE,
                               output_dir = NULL,
                               output_width = 10,
                               output_height = 10,
                               save_network_plot=TRUE, 
                               save_apearres_obj=TRUE, 
                               ...){
  
  
  ############### Check before plot cluster ########################### 
  # Check for `output_dir` only if needed to save
  if((save_network_plot || save_apearres_obj) && is.null(output_dir)) {
    stop("Please provide an output directory for the output_dir parameter")
  }else{
    if(!is.null(output_dir)) {
      if(!(save_network_plot || save_apearres_obj)){
        message("Note: Output directory is provided but save_network_plot and/or save_apearres_obj is FALSE, some results will be saved.")
      }
    }else{
      # Ensure the output directory exists or create it
      if (!dir.exists(output_dir)) {
        message(paste0("Output directory:",output_dir, " does NOT exist! Will create one!") )
        dir.create(output_dir, recursive = TRUE)
      }
    } 
  }
  
  # Check for color mapping
  
  check_color_mapping_format <- function(color_mapping) {
    # Check if color_mapping is a named vector
    if (!is.vector(color_mapping) || is.null(names(color_mapping))) {
      return(FALSE)
    }
    
    # Check if all values are valid colors
    valid_colors <- colors()
    all_valid_colors <- all(color_mapping %in% valid_colors | grepl("^#[0-9A-Fa-f]{6}$", color_mapping))
    
    return(all_valid_colors)
  }
  
  
  if(is.null(color_mapping)) {
    message("Note: Color key is not provided. Will use the default color scheme!")
  }else{
    if(!check_color_mapping_format(color_mapping)){
      stop("Note: Color key provided is in the wrong format! Please double check")
    }
  }
  
  #########################################################################################
  
  
  ######################Unpacking findPathCluster and parameters##################
  # Extract parameters from the findPathClusterres object
  max_sign_pathway <- findPathClusterres$params$max_sign_pathway
  min_leading_edge_threshold <- findPathClusterres$params$min_leading_edge_threshold
  
  findPathClusterres <- findPathClusterres$clusterres
  ############################################################################
  
  
  
  ########################Create apearres##############################################
  #define unique_clusters
  unique_clusters <- lapply(apear_input_list, function(df) names(df))
  
  # create a list to store results for apear outputs of each cluster
  apearres <- vector("list",length = length(apear_input_list))
  names(apearres) <- names(apear_input_list)
  
  
  # loop through each comparison 
  for (x in 1:length(apear_input_list)) {
    # create a tmp cluster to work with in this iteration
    tmp_cluster <- unique_clusters[[names(apear_input_list)[x]]]
    
    # create a list to store result for each cluster
    apearres[[x]] <- vector("list",length = length(tmp_cluster))
    names(apearres[[x]]) <- tmp_cluster
    
    
    # loop through each cluster
    for (y in 1:length(tmp_cluster)) {
      apearres[[x]][[y]] <- vector("list",length = 2)
      names(apearres[[x]][[y]])<- c("apear_input","netp") #store both inputs and networkplot
      
      
      apearres[[x]][[y]][["apear_input"]] <- apear_input_list[[x]][[y]]
    }
  }
  #########################################################################################
  
  
  
  ########################Define color usage##############################################
  get_color <- function(key) {
    
    generate_random_color <- function() {
      red <- runif(1, 0, 1)
      green <- runif(1, 0, 1)
      blue <- runif(1, 0, 1)
      random_color <- rgb(red, green, blue)
      return(random_color)
    }
    
    if(key %in% names(color_mapping)){
      return(color_mapping[key])
    }else{
      random_color <- generate_random_color()
      # Update the color_mapping with the new key and generated random color
      color_mapping <<- c(color_mapping, setNames(random_color, key))
      message(paste0(key," is NOT found in color maps! Will use ",random_color," instead. Please double check!"))
      return(random_color)
    }
    
    #return(ifelse(key %in% names(color_mapping), color_mapping[key], generate_random_color()))
  }
  
  
  #########################################################################################
  
  
  ########################Output and Plotting##############################################
  
  # Loop through each comparison in apear_input_list
  for (x in 1:length(apear_input_list)) {
    tmp_apear_input <- apear_input_list[[x]]
    tmp_output_dir <- paste0(output_dir, names(apear_input_list)[x], "/")
    
    if (!dir.exists(tmp_output_dir)) {
      dir.create(tmp_output_dir, recursive = TRUE)
    }
    
    if(!is.null(color_mapping)){
      tmp_Hgroup <- sapply(strsplit(names(apear_input_list),"_vs_"),`[`,1)[x]
      tmp_Lgroup <- sapply(strsplit(names(apear_input_list),"_vs_"),`[`,2)[x]
      Hcolor <- get_color(tmp_Hgroup)
      Lcolor <- get_color(tmp_Lgroup)
    }else{
      Hcolor <- "firebrick"
      Lcolor <- "steelblue"
    }
    
    # Loop through each tmp_apear_input and pull out 1 cluster at a time
    for (y in 1:length(tmp_apear_input)) {
      
      # Obtain the parameters and results from findPathCluster_obj
      tmp_cluster <- tmp_apear_input[[y]]  
      tmp_parameters <- findPathClusterres[[x]][[y]][["findPathCluster_opt_inputs"]]
      clustMethod <- tmp_parameters[[1]]
      ClusterSize <- tmp_parameters[[2]]
      
      # Calculate positive and negative NES pathways
      pos_NES_numb <- sum(tmp_cluster$NES[!c(tmp_cluster$aPEAR_Cluster %in% c("NA"))]>0) 
      neg_NES_numb <- sum(tmp_cluster$NES[!c(tmp_cluster$aPEAR_Cluster %in% c("NA"))]<0)
      
      tmp_res<- findPathClusterres[[x]][[y]][["findPathClusterres"]]
      
      msg_name <- paste0(names(apear_input_list)[x], "(", names(tmp_apear_input)[y], ")")
      
      
      # whether to show the parameter in the final output
      if(show_parameters){
        parameters_string <- paste0(clustMethod," clustering; Minimum cluster size: ",ClusterSize, "\nMaxmium # of significant pathways: ", max_sign_pathway, "\nMinimum # of leading edge (LE) genes per pathway: ", min_leading_edge_threshold,"\n# of pos NES pathways: ",pos_NES_numb,"; # of neg NES pathways: ",neg_NES_numb)
      }else{
        parameters_string <-""
      }
      
      
      # Plot the network plot
      apearres_output <- aPEAR::plotPathClusters(tmp_cluster,
                                                 sim = tmp_res$similarity,
                                                 clusters = tmp_res$clusters,
                                                 fontSize = fontSize,
                                                 nodeSize = "number_of_LE_genes",
                                                 repelLabels = repelLabels, drawEllipses = drawEllipses)+
        ggtitle(paste0(msg_name),
                subtitle = parameters_string)+
        labs(size = "# of LE genes")+
        scale_colour_gradient2(low = Lcolor,high = Hcolor, mid = "white", midpoint = 0, name = "Normalized\nEnrichment\nScore") 
      
      if(save_network_plot){
        if(dir.exists(tmp_output_dir)){
          #save the network plot as pdfs 
          pdf(paste0(tmp_output_dir, names(tmp_apear_input)[y], "_aPEAR_", clustMethod, "_network.pdf"),width = output_width, height = output_height)
          print(apearres_output)
          dev.off()
        }
      }
      
      #update the network slots
      apearres[[x]][[y]][["netp"]] <- apearres_output 
      
    }
  }
  #########################################################################################
  
  
  
  
  
  ##################################save the apearres obj#######################################
  
  if(save_apearres_obj){
    if(!is.null(output_dir)){
      saveRDS(apearres,paste0(output_dir,"aPEARres.rds"))
    }else{
      message("Please provide an output directory for the output_dir paramater")
    }  
  }
  #########################################################################################
  
  
  
  return(apearres)
  
}


