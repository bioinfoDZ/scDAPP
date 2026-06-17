# ---- Integration dispatcher and input prep ----

#' Run multi-sample integration (dispatcher)
#'
#' @param integration_method One of \code{integration_method_choices()}; default \code{"RISC"}.
#' @param sample_metadata data.frame with Sample, Condition, Code.
#' @param sobjlist Named list of per-sample Seurat objects (for prep), or NULL if using \code{tmpobjdir}.
#' @param mdlist Named list of metadata per sample; required for RISC.
#' @param tmpobjdir Path to saved per-sample objects (\code{Code.rds}).
#' @param datadir Cell Ranger folder; used only when \code{input_seurat_obj = FALSE} (h5 at prep).
#' @param outdir_int Integrated analysis output folder.
#' @param outdir_indi Individual-sample output folder.
#' @param pcs_int PCs for integration, or \code{"auto"} for stability-based selection.
#' @param res_int Louvain resolution, or \code{"auto"}.
#' @param risc_reference RISC reference sample (Code or Sample); NULL = auto.
#' @param RISC_louvain_neighbors Louvain neighbors (RISC only).
#' @param workernum Parallel workers for final RISC integration (\code{InPlot}, \code{scMultiIntegrate}).
#' @param stability_workernum Parallel workers for auto stability sweep; \code{NULL} uses \code{workernum}.
#' @param input_seurat_obj If \code{TRUE}, RISC counts come from \code{sobjlist} at prep (not datadir).
#' @param stability_outdir Output directory for stability sweep (default \code{outdir_int/cluster_stability}).
#' @param stability_numreps Bootstrap replicates when \code{pcs_int} or \code{res_int} is \code{"auto"}.
#' @param stability_sweep_maxPCs PC grid for auto \code{pcs_int}.
#' @param stability_sweep_res Resolution grid for auto \code{res_int}.
#' @param stability_propcells.perrep Cell fraction per bootstrap replicate.
#' @return List with \code{sobjint}, \code{int_clust_lab}, \code{ip}, \code{config}.
#'   For RISC, also \code{refscore} (named by sample Code) and
#'   \code{selected_risc_reference}; both are \code{NULL} for Seurat methods.
#'   When auto tuning was used, includes \code{stability} and \code{auto_selection}.
#' @export
run_integration <- function(
    integration_method = "RISC",
    sample_metadata,
    sobjlist = NULL,
    mdlist = NULL,
    tmpobjdir = NULL,
    datadir = NULL,
    outdir_int,
    outdir_indi = NULL,
    pcs_int = 30L,
    res_int = 0.5,
    risc_reference = NULL,
    RISC_louvain_neighbors = 10L,
    workernum = 1L,
    stability_workernum = NULL,
    input_seurat_obj = FALSE,
    stability_outdir = NULL,
    stability_numreps = 50L,
    stability_sweep_maxPCs = c(5, 10, 15, 20, 25, 30, 40, 50),
    stability_sweep_res = seq(0.1, 1.5, by = 0.2),
    stability_propcells.perrep = 0.8
) {
  integration_method <- match.arg(integration_method, integration_method_choices())
  auto_requested <- .integration_param_needs_auto(pcs_int, res_int)

  auto_result <- NULL
  if (auto_requested) {
    if (is.null(stability_outdir)) {
      stability_outdir <- file.path(outdir_int, "cluster_stability")
    }
    if (is.null(tmpobjdir) && !is.null(outdir_indi)) {
      tmpobjdir <- file.path(outdir_indi, ".tmp_Seurat_objects")
    }
    matlist_path <- file.path(outdir_int, "data_objects", ".concatmatrix.rds")
    auto_result <- auto_integration_cluster_params(
      integration_method = integration_method,
      sample_metadata = sample_metadata,
      pcs_int = pcs_int,
      res_int = res_int,
      tmpobjdir = tmpobjdir,
      mdlist = mdlist,
      matlist_path = matlist_path,
      datadir = datadir,
      outdir = stability_outdir,
      numreps = stability_numreps,
      propcells.perrep = stability_propcells.perrep,
      sweep_maxPCs = stability_sweep_maxPCs,
      sweep_res = stability_sweep_res,
      workernum = workernum,
      stability_workernum = stability_workernum,
      risc_reference = risc_reference,
      RISC_louvain_neighbors = RISC_louvain_neighbors,
      input_seurat_obj = input_seurat_obj,
      verbose = TRUE
    )
    resolved <- resolve_integration_cluster_params(pcs_int, res_int, auto_result)
    pcs_int <- resolved$pcs_int
    res_int <- resolved$res_int
  }

  config <- resolve_integration_config(
    integration_method,
    pcs_int,
    res_int
  )
  check_integration_dependencies(config, pcs_int = pcs_int, res_int = res_int)

  if (config$engine == "RISC") {
    if (is.null(mdlist)) {
      stop("mdlist is required for RISC integration.", call. = FALSE)
    }
    matlist_path <- file.path(outdir_int, "data_objects", ".concatmatrix.rds")
    if (!file.exists(matlist_path) && !isTRUE(input_seurat_obj) && is.null(datadir)) {
      stop(
        "RISC requires prepare_integration_inputs() or datadir when input_seurat_obj is FALSE.",
        call. = FALSE
      )
    }
    res <- run_risc_integration(
      mdlist = mdlist,
      sample_metadata = sample_metadata,
      outdir_int = outdir_int,
      config = config,
      matlist_path = matlist_path,
      risc_reference = risc_reference,
      workernum = workernum,
      RISC_louvain_neighbors = RISC_louvain_neighbors
    )
  } else {
    if (is.null(tmpobjdir)) {
      if (!is.null(outdir_indi)) {
        tmpobjdir <- file.path(outdir_indi, ".tmp_Seurat_objects")
      } else {
        stop("tmpobjdir or outdir_indi required for Seurat integration.", call. = FALSE)
      }
    }
    res <- run_seurat_integration(
      tmpobjdir = tmpobjdir,
      sample_metadata = sample_metadata,
      outdir_int = outdir_int,
      config = config,
      workernum = workernum
    )
  }

  outdir_int_objects <- file.path(outdir_int, "data_objects")
  dir.create(outdir_int_objects, recursive = TRUE, showWarnings = FALSE)
  saveRDS(res$sobjint, file.path(outdir_int_objects, "Seurat-object_integrated.rds"))

  out <- c(res, list(config = config))
  if (!is.null(auto_result)) {
    out$stability <- auto_result$stability
    out$auto_selection <- list(
      pcs_int = auto_result$pcs_int,
      res_int = auto_result$res_int,
      selected_params_i = auto_result$selected_params_i,
      sweep_dir = auto_result$sweep_dir
    )
  }
  out
}

#' Prepare integration inputs: annotate objects, save tmp, build concat matrix
#'
#' Call before \code{run_integration()} from the pipeline Rmd when \code{sobjlist} is in memory.
#'
#' @param sobjlist Named list of Seurat objects.
#' @param sample_metadata data.frame with Code, Sample, Condition.
#' @param outdir_indi Individual-sample output directory.
#' @param datadir Raw h5 root; required only when \code{input_seurat_obj = FALSE}.
#' @param outdir_int Integrated output directory.
#' @param input_seurat_obj Logical.
#' @return List with \code{mdlist}, \code{tmpobjdir}, and \code{matlist_path}.
#' @export
prepare_integration_inputs <- function(
    sobjlist,
    sample_metadata,
    outdir_indi,
    datadir = NULL,
    outdir_int,
    input_seurat_obj = FALSE
) {
  sobjlist <- .annotate_sobjlist_for_integration(sobjlist, sample_metadata)
  tmpobjdir <- file.path(outdir_indi, ".tmp_Seurat_objects")
  dir.create(tmpobjdir, recursive = TRUE, showWarnings = FALSE)

  outdir_int_objects <- file.path(outdir_int, "data_objects")
  dir.create(outdir_int_objects, recursive = TRUE, showWarnings = FALSE)
  matlist_path <- file.path(outdir_int_objects, ".concatmatrix.rds")

  if (isTRUE(input_seurat_obj)) {
    matlist <- .build_matlist_from_sobjlist(sobjlist, sample_metadata)
  } else {
    if (is.null(datadir) || !nzchar(datadir)) {
      stop("datadir is required when input_seurat_obj is FALSE.", call. = FALSE)
    }
  }

  mdlist <- lapply(sample_metadata$Code, function(code) {
    sobj <- sobjlist[[code]]
    saveRDS(sobj, file.path(tmpobjdir, paste0(code, ".rds")))
    sobj@meta.data
  })
  names(mdlist) <- sample_metadata$Code

  if (!isTRUE(input_seurat_obj)) {
    matlist <- .build_matlist_from_h5_mdlist(mdlist, sample_metadata, datadir)
  }
  saveRDS(matlist, matlist_path)
  rm(matlist, sobjlist)
  invisible(gc(full = TRUE, reset = FALSE, verbose = FALSE))

  list(mdlist = mdlist, tmpobjdir = tmpobjdir, matlist_path = matlist_path)
}

# ---- RISC integration ----

#' Run RISC scMultiIntegrate workflow and build Integrated_RISC Seurat object.
#'
#' @param mdlist Named list of per-sample metadata data.frames (rownames = cell barcodes).
#' @param sample_metadata data.frame with Sample, Condition, Code.
#' @param outdir_int Integrated analysis output directory.
#' @param config List from \code{resolve_integration_config()}.
#' @param matlist_path Path to \code{.concatmatrix.rds} from \code{prepare_integration_inputs()}.
#' @param risc_reference Optional Code or Sample name for reference; NULL = auto.
#' @param workernum Parallel workers.
#' @param RISC_louvain_neighbors Nearest neighbors for Louvain clustering.
#' @return List with \code{sobjint}, \code{int_clust_lab}, \code{ip} (InPlot patchwork),
#'   \code{refscore} (named numeric vector by sample Code), and
#'   \code{selected_risc_reference} (Code of the reference sample used).
#' @keywords internal
run_risc_integration <- function(
    mdlist,
    sample_metadata,
    outdir_int,
    config,
    matlist_path = NULL,
    risc_reference = NULL,
    workernum = 1L,
    RISC_louvain_neighbors = 10L
) {
  outdir_int_objects <- file.path(outdir_int, "data_objects")
  dir.create(outdir_int_objects, recursive = TRUE, showWarnings = FALSE)

  if (is.null(matlist_path) || !nzchar(matlist_path)) {
    matlist_path <- file.path(outdir_int_objects, ".concatmatrix.rds")
  }
  if (!file.exists(matlist_path)) {
    stop(
      "Missing ", matlist_path, "; call prepare_integration_inputs() first.",
      call. = FALSE
    )
  }
  matlist <- readRDS(matlist_path)

  bigmat <- do.call(cbind, matlist)
  num_nonzeros <- tabulate(bigmat@i + 1L, nbins = nrow(bigmat))
  joint_filt_genes <- rownames(bigmat)[num_nonzeros >= 3L]

  risclist <- lapply(sample_metadata$Code, function(code) {
    md <- mdlist[[code]]
    mat0 <- matlist[[code]]
    mat0 <- mat0[rownames(mat0) %in% joint_filt_genes, , drop = FALSE]
    mat0 <- mat0[match(joint_filt_genes, rownames(mat0)), , drop = FALSE]
    coldata0 <- md
    barcodes <- stringr::str_split_fixed(rownames(coldata0), "-", 2)[, 1]
    barcodes <- paste0(coldata0$orig.ident, ".", barcodes)
    coldata0 <- cbind(barcodes, coldata0)
    rowdata0 <- data.frame(Symbol = rownames(mat0), row.names = rownames(mat0))
    RISC::readsc(mat0, coldata0, rowdata0, is.filter = FALSE)
  })
  names(risclist) <- sample_metadata$Code

  rm(bigmat, num_nonzeros, joint_filt_genes)
  invisible(gc(full = TRUE, reset = FALSE, verbose = FALSE))

  process0 <- function(obj0) {
    obj0 <- RISC::scFilter(obj0, min.UMI = 0, max.UMI = Inf, min.gene = 0, min.cell = 0, is.filter = FALSE)
    obj0 <- RISC::scNormalize(obj0, ncore = 1)
    obj0 <- RISC::scDisperse(obj0)
    obj0
  }

  cl <- parallel::makeCluster(workernum)
  doParallel::registerDoParallel(cl)
  risclist <- foreach::foreach(
    dat0 = risclist,
    .packages = c("RISC")
  ) %dopar% process0(dat0)
  parallel::stopCluster(cl)
  invisible(gc(full = TRUE, reset = FALSE, verbose = FALSE))
  names(risclist) <- sample_metadata$Code

  var0 <- Reduce(intersect, lapply(risclist, function(x) x@rowdata$Symbol))

  grDevices::pdf(NULL)
  ip <- RISC::InPlot(risclist, var.gene = var0, Std.cut = 0.95, ncore = workernum)
  ip <- patchwork::wrap_plots(ip)
  grDevices::dev.off()

  ref <- NULL
  if (!is.null(risc_reference)) {
    if (any(risc_reference %in% sample_metadata$Code)) {
      ref <- which(sample_metadata$Code == risc_reference)[1]
    } else if (any(risc_reference %in% sample_metadata$Sample)) {
      ref <- which(sample_metadata$Sample == risc_reference)[1]
    }
  }

  numclusts <- vapply(risclist, function(dat0) length(unique(dat0@coldata$seurat_clusters)), integer(1))
  numcells_per_sample <- vapply(risclist, function(dat0) nrow(dat0@coldata), numeric(1))
  numcells_per_sample <- numcells_per_sample / max(numcells_per_sample)
  numclusts <- numclusts * numcells_per_sample

  pbvar <- vapply(risclist, function(dat0) {
    mat <- dat0@assay$logcount
    md <- dat0@coldata
    pb <- scDAPP::pseudobulk(obj = mat, metadata = md, grouping_colname_in_md = "seurat_clusters")
    numcells <- table(md$seurat_clusters)
    pb <- sweep(pb, 2, numcells, FUN = "/")
    clustervar <- apply(pb, 2, var)
    mean(clustervar)
  }, numeric(1))

  refscore <- numclusts * pbvar
  names(refscore) <- sample_metadata$Code
  if (is.null(ref)) {
    ref <- which.max(refscore)
  }
  selected_risc_reference <- sample_metadata$Code[ref]

  if (ref != 1L) {
    data0 <- list(risclist[[ref]])
    names(data0) <- names(risclist)[ref]
    for (i in seq_along(risclist)) {
      if (i != ref) {
        data0[[names(risclist)[i]]] <- risclist[[i]]
      }
    }
  } else {
    data0 <- risclist
  }
  rm(risclist)
  invisible(gc(full = TRUE, reset = FALSE, verbose = FALSE))

  eigens <- config$pcs_int
  data0 <- RISC::scMultiIntegrate(
    objects = data0,
    eigens = eigens,
    add.Id = NULL,
    var.gene = var0,
    align = "OLS",
    npc = 50,
    adjust = TRUE,
    ncore = workernum
  )
  rm(var0)
  invisible(gc(full = TRUE, reset = FALSE, verbose = FALSE))

  data0 <- RISC::scUMAP(data0, npc = eigens, use = "PLS")
  data0 <- RISC::scCluster(
    data0,
    slot = "cell.pls",
    method = "louvain",
    npc = eigens,
    res = config$res_int,
    neighbor = RISC_louvain_neighbors
  )

  data0@coldata$Cluster <- .remap_clusters_by_size(data0@coldata$Cluster)

  saveRDS(data0, file.path(outdir_int_objects, "RISC-object_integrated.rds"))

  mat <- do.call(cbind, data0@assay$logcount)
  md <- data0@coldata
  umap <- data0@DimReduction$cell.umap
  pca <- data0@DimReduction$cell.pls

  if ("scBarcode" %in% colnames(md) && !"Barcode" %in% colnames(md)) {
    md$Barcode <- md$scBarcode
  }

  int_clust_lab <- .integration_cluster_column(config)
  colnames(md)[ncol(md)] <- int_clust_lab

  risc_assay <- config$integrated_assay
  sobjint <- Seurat::CreateSeuratObject(
    counts = Seurat::CreateAssayObject(data = mat),
    assay = risc_assay,
    project = "Integrated",
    meta.data = md
  )
  sobjint[["umap"]] <- Seurat::CreateDimReducObject(umap, assay = risc_assay, key = "UMAP_")
  sobjint[["pca"]] <- Seurat::CreateDimReducObject(pca, assay = risc_assay, key = "PCA_")

  rm(umap, pca, mat, md, data0)
  invisible(gc(full = TRUE, reset = FALSE, verbose = FALSE))

  sobjint <- .append_concat_rna_assay(sobjint, outdir_int_objects, sample_metadata)
  sobjint <- .finalize_sobjint_clusters(sobjint, config, int_clust_lab)

  sobjint$Code <- factor(sobjint$Code, levels = sample_metadata$Code)
  sobjint$Condition <- factor(sobjint$Condition, levels = levels(sample_metadata$Condition))

  list(
    sobjint = sobjint,
    int_clust_lab = int_clust_lab,
    ip = ip,
    refscore = refscore,
    selected_risc_reference = selected_risc_reference
  )
}

# ---- Seurat v5 integration ----

#' Join Seurat v5 split RNA layers before merge/split integration workflow.
#' @keywords internal
.join_assay_layers_if_needed <- function(obj, assay_name) {
  if (!assay_name %in% Seurat::Assays(obj)) {
    return(obj)
  }
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    return(obj)
  }
  layers <- tryCatch(
    as.character(SeuratObject::Layers(obj, assay = assay_name)),
    error = function(e) character(0)
  )
  if (length(layers) > 1L) {
    obj <- SeuratObject::JoinLayers(obj, assay = assay_name)
  }
  obj
}

#' Map integration_method string to Seurat::IntegrateLayers function.
#' @keywords internal
.get_integrate_layers_fn <- function(method_name) {
  switch(
    method_name,
    CCAIntegration = Seurat::CCAIntegration,
    RPCAIntegration = Seurat::RPCAIntegration,
    HarmonyIntegration = Seurat::HarmonyIntegration,
    stop("Unknown IntegrateLayers method: ", method_name, call. = FALSE)
  )
}

#' SCTransform merged object with split RNA layers for SCT IntegrateLayers path.
#' @keywords internal
.sctransform_merged_split_rna <- function(obj) {
  sct_args <- list(
    object = obj,
    assay = "RNA",
    verbose = FALSE,
    vst.flavor = "v2"
  )
  if (requireNamespace("glmGamPoi", quietly = TRUE)) {
    sct_args$method <- "glmGamPoi"
  }
  suppressWarnings(do.call(Seurat::SCTransform, sct_args))
}

#' Merge per-sample RDS from tmpobjdir into one object for Seurat integration.
#' @keywords internal
.seurat_merge_from_tmpobjdir <- function(tmpobjdir, sample_metadata) {
  codes <- sample_metadata$Code
  sobjlist <- .read_sobjlist_from_dir(tmpobjdir, codes)
  obj <- if (length(sobjlist) == 1L) {
    sobjlist[[1]]
  } else {
    merge(x = sobjlist[[1]], y = sobjlist[-1])
  }
  rm(sobjlist)
  invisible(gc(full = TRUE, reset = FALSE, verbose = FALSE))
  obj$Code <- factor(obj$Code, levels = sample_metadata$Code)
  if (requireNamespace("SeuratObject", quietly = TRUE) && "RNA" %in% Seurat::Assays(obj)) {
    obj <- .join_assay_layers_if_needed(obj, "RNA")
  }
  obj
}

#' Split RNA assay by sample Code and cap PCs by smallest layer cell count.
#' @return list(obj, split_by, npcs, dims)
#' @keywords internal
.seurat_prepare_split_pcs <- function(obj, npcs) {
  split_by <- obj$Code
  min_layer_cells <- min(as.integer(table(split_by)))
  npcs <- .stability_cap_pcs(npcs, min_layer_cells)
  dims <- seq_len(npcs)
  obj[["RNA"]] <- split(obj[["RNA"]], f = split_by)
  list(obj = obj, split_by = split_by, npcs = npcs, dims = dims)
}

#' Run normalization, PCA, and IntegrateLayers at npcs without clustering.
#' @keywords internal
.seurat_integrate_layers_at_pcs <- function(obj, config, npcs) {
  prep <- .seurat_prepare_split_pcs(obj, npcs)
  obj <- prep$obj
  npcs <- prep$npcs
  dims <- prep$dims
  integrate_fn <- .get_integrate_layers_fn(config$integrate_layers_method)

  if (config$normalization_method == "SCT") {
    obj <- .sctransform_merged_split_rna(obj)
    obj <- Seurat::RunPCA(obj, assay = "SCT", npcs = npcs, verbose = FALSE)
    il_args <- list(
      object = obj,
      method = integrate_fn,
      normalization.method = "SCT",
      orig.reduction = "pca",
      new.reduction = config$cluster_reduction,
      dims = dims,
      k.filter = 200L,
      verbose = FALSE
    )
  } else {
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    obj <- Seurat::FindVariableFeatures(obj, verbose = FALSE)
    obj <- Seurat::ScaleData(obj, verbose = FALSE)
    obj <- Seurat::RunPCA(obj, npcs = npcs, verbose = FALSE)
    il_args <- list(
      object = obj,
      method = integrate_fn,
      orig.reduction = "pca",
      new.reduction = config$cluster_reduction,
      dims = dims,
      k.filter = 200L,
      verbose = FALSE
    )
  }
  do.call(Seurat::IntegrateLayers, il_args)
}

#' Louvain cluster and UMAP on batch-corrected reduction after IntegrateLayers.
#' @keywords internal
.seurat_cluster_umap_on_reduction <- function(
    obj,
    config,
    npcs,
    res_int,
    cluster_name = NULL
) {
  emb_dims <- tryCatch(
    seq_len(min(
      npcs,
      ncol(Seurat::Embeddings(obj, config$cluster_reduction))
    )),
    error = function(e) seq_len(npcs)
  )
  obj <- Seurat::FindNeighbors(
    obj,
    reduction = config$cluster_reduction,
    dims = emb_dims,
    verbose = FALSE
  )
  if (is.null(cluster_name)) {
    cluster_name <- .integration_cluster_column(config)
  }
  obj <- Seurat::FindClusters(
    obj,
    resolution = res_int,
    cluster.name = cluster_name,
    verbose = FALSE
  )
  obj@meta.data[[cluster_name]] <- .remap_clusters_by_size(obj@meta.data[[cluster_name]])
  obj <- Seurat::RunUMAP(
    obj,
    reduction = config$cluster_reduction,
    dims = emb_dims,
    reduction.name = config$umap_reduction,
    verbose = FALSE
  )
  list(obj = obj, int_clust_lab = cluster_name, emb_dims = emb_dims)
}

#' Stability sweep path: full IntegrateLayers + cluster at one PC/res grid point.
#' @keywords internal
.seurat_stability_integrate_cluster_at <- function(sobj_merged, config, pcs_int, res_int) {
  obj <- .seurat_integrate_layers_at_pcs(sobj_merged, config, pcs_int)
  .seurat_stability_cluster_at(obj, config, pcs_int, res_int)
}

#' Stability sweep: re-cluster an already-integrated object at one grid point.
#' @keywords internal
.seurat_stability_cluster_at <- function(obj, config, maxdim, res) {
  dims <- seq_len(as.integer(maxdim)[1])
  params_i <- paste0("PCs_1-", maxdim, ".res_", res)
  cluster_name <- .stability_sweep_cluster_name(params_i)
  obj <- Seurat::FindNeighbors(
    obj,
    reduction = config$cluster_reduction,
    dims = dims,
    verbose = FALSE
  )
  obj <- Seurat::FindClusters(
    obj,
    resolution = res,
    cluster.name = cluster_name,
    verbose = FALSE
  )
  .stability_extract_cluster_df(obj@meta.data, params_i)
}
#' Run Seurat v5 IntegrateLayers integration and attach RNA counts for pseudobulk DE.
#'
#' @param tmpobjdir Directory with per-sample .rds objects (Code.rds).
#' @param sample_metadata data.frame with Code, Condition, Sample.
#' @param outdir_int Integrated output directory.
#' @param config List from \code{resolve_integration_config()}.
#' @param workernum Not used for Seurat integration (reserved).
#' @return List with \code{sobjint}, \code{int_clust_lab}, \code{ip} NULL.
#' @keywords internal
run_seurat_integration <- function(
    tmpobjdir,
    sample_metadata,
    outdir_int,
    config,
    workernum = 1L
) {
  outdir_int_objects <- file.path(outdir_int, "data_objects")
  dir.create(outdir_int_objects, recursive = TRUE, showWarnings = FALSE)

  obj <- .seurat_merge_from_tmpobjdir(tmpobjdir, sample_metadata)
  npcs <- as.integer(config$pcs_int)[1]
  obj <- .seurat_integrate_layers_at_pcs(obj, config, npcs)
  clust <- .seurat_cluster_umap_on_reduction(
    obj,
    config,
    npcs,
    config$res_int
  )
  obj <- clust$obj
  int_clust_lab <- clust$int_clust_lab

  if (requireNamespace("SeuratObject", quietly = TRUE) && "RNA" %in% Seurat::Assays(obj)) {
    obj <- .join_assay_layers_if_needed(obj, "RNA")
  }

  sobjint <- obj
  rm(obj)
  invisible(gc(full = TRUE, reset = FALSE, verbose = FALSE))

  if (file.exists(file.path(outdir_int_objects, ".concatmatrix.rds"))) {
    sobjint <- .append_concat_rna_assay(sobjint, outdir_int_objects, sample_metadata)
  }

  sobjint <- .finalize_sobjint_clusters(sobjint, config, int_clust_lab)
  sobjint$Code <- factor(sobjint$Code, levels = sample_metadata$Code)
  sobjint$Condition <- factor(sobjint$Condition, levels = levels(sample_metadata$Condition))

  list(
    sobjint = sobjint,
    int_clust_lab = int_clust_lab,
    ip = NULL,
    refscore = NULL,
    selected_risc_reference = NULL
  )
}
