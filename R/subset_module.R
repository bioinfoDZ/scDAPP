# Subset an integrated object, recluster, and run comparative modules.
#
# RISC keeps the parent Integrated_RISC correction and reclusters a new PCA of
# that assay. Harmony, CCA, RPCA, and SCT re-integrate the subset. Auto PC or
# resolution on the RISC path bootstraps that frozen PCA. It does not call
# cluster_stability_sweep(), which re-integrates at every PC count.

.subset_require_columns <- function(df, cols, what) {
  missing <- setdiff(cols, colnames(df))
  if (length(missing)) {
    stop(what, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(df)
}

.subset_standardize_sample_metadata <- function(sample_metadata, sample_column, condition_column) {
  md <- as.data.frame(sample_metadata, stringsAsFactors = FALSE)
  if (!"Code" %in% colnames(md)) {
    if (!sample_column %in% colnames(md)) {
      stop(
        "sample_metadata needs Code or a column named '", sample_column, "'.",
        call. = FALSE
      )
    }
    md$Code <- md[[sample_column]]
  }
  if (!"Condition" %in% colnames(md)) {
    if (!condition_column %in% colnames(md)) {
      stop(
        "sample_metadata needs Condition or a column named '", condition_column, "'.",
        call. = FALSE
      )
    }
    md$Condition <- md[[condition_column]]
  }
  if (!"Sample" %in% colnames(md)) md$Sample <- md$Code
  md$Code <- as.character(md$Code)
  md$Sample <- as.character(md$Sample)
  md$Condition <- as.character(md$Condition)
  md
}

.subset_standardize_object_metadata <- function(obj, sample_column, condition_column) {
  md <- obj@meta.data
  if (!sample_column %in% colnames(md)) {
    stop("Object metadata has no sample column '", sample_column, "'.", call. = FALSE)
  }
  if (!condition_column %in% colnames(md)) {
    stop("Object metadata has no condition column '", condition_column, "'.", call. = FALSE)
  }
  obj$Code <- as.character(md[[sample_column]])
  obj$Condition <- as.character(md[[condition_column]])
  obj
}

# Integrated objects can carry Condition as a factor whose values were dropped
# when sample_metadata$Condition had no factor levels. Restore it from the table.
.subset_apply_condition_from_metadata <- function(obj, sample_metadata) {
  code <- as.character(obj$Code)
  cond <- as.character(sample_metadata$Condition)[
    match(code, as.character(sample_metadata$Code))
  ]
  if (anyNA(cond)) {
    stop("sample_metadata is missing Condition for some cells.", call. = FALSE)
  }
  obj$Condition <- cond
  obj
}

.subset_select_cells <- function(obj, clusters, cells, cluster_column) {
  if (!is.null(cells) && !is.null(clusters)) {
    stop("Pass clusters or cells, not both.", call. = FALSE)
  }
  if (is.null(cells) && is.null(clusters)) {
    stop("Pass clusters or cells to keep.", call. = FALSE)
  }
  if (!is.null(cells)) {
    cells <- as.character(cells)
    if (!any(cells %in% colnames(obj)) && "Barcode" %in% colnames(obj@meta.data)) {
      hit <- match(cells, as.character(obj$Barcode))
      if (anyNA(hit)) {
        stop(
          sum(is.na(hit)), " cells were not found in colnames or Barcode.",
          call. = FALSE
        )
      }
      cells <- colnames(obj)[hit]
    }
    missing <- setdiff(cells, colnames(obj))
    if (length(missing)) {
      stop(
        "Some cells are not in the object: ",
        paste(utils::head(missing, 5L), collapse = ", "),
        call. = FALSE
      )
    }
    return(unique(cells))
  }
  if (!cluster_column %in% colnames(obj@meta.data)) {
    stop("cluster_column not in metadata: ", cluster_column, call. = FALSE)
  }
  clusters <- as.character(clusters)
  have <- as.character(obj@meta.data[[cluster_column]])
  if (!any(have %in% clusters)) {
    stop("None of the requested clusters are in ", cluster_column, ".", call. = FALSE)
  }
  colnames(obj)[have %in% clusters]
}

.subset_parent_method <- function(obj, integration_method) {
  if (!is.null(integration_method) && nzchar(as.character(integration_method)[1])) {
    return(match.arg(as.character(integration_method)[1], integration_method_choices()))
  }
  stored <- obj@misc$scDAPP$integration_method
  if (is.null(stored) || !nzchar(as.character(stored)[1])) {
    stop(
      "Could not read integration_method from sobjint@misc$scDAPP. ",
      "Pass integration_method.",
      call. = FALSE
    )
  }
  match.arg(as.character(stored)[1], integration_method_choices())
}

.subset_filter_metadata <- function(obj, sample_metadata) {
  codes <- unique(as.character(obj$Code))
  md <- sample_metadata[as.character(sample_metadata$Code) %in% codes, , drop = FALSE]
  if (!nrow(md)) {
    stop("No sample_metadata rows match the subset Code values.", call. = FALSE)
  }
  missing_codes <- setdiff(codes, md$Code)
  if (length(missing_codes)) {
    stop(
      "sample_metadata is missing Code values present in the subset: ",
      paste(missing_codes, collapse = ", "),
      call. = FALSE
    )
  }
  obj$Code <- factor(as.character(obj$Code), levels = as.character(md$Code))
  cond_levels <- unique(as.character(md$Condition))
  obj$Condition <- factor(as.character(obj$Condition), levels = cond_levels)
  md$Code <- as.character(md$Code)
  md$Condition <- as.character(md$Condition)
  list(obj = obj, sample_metadata = md)
}

.subset_check_comps <- function(sample_metadata, comps) {
  comps <- normalize_comps(comps)
  need <- unique(c(as.character(comps$c0), as.character(comps$c1)))
  have <- unique(as.character(sample_metadata$Condition))
  miss <- setdiff(need, have)
  if (length(miss)) {
    stop(
      "Subset does not contain comps condition levels: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  comps
}

.subset_pc_cap <- function(n_cells, n_features = NULL) {
  cap <- as.integer(n_cells)[1] - 1L
  if (!is.null(n_features)) cap <- min(cap, as.integer(n_features)[1] - 1L)
  cap
}

.subset_cap_pc_values <- function(pcs, cap) {
  pcs <- unique(as.integer(pcs))
  pcs <- sort(pcs[pcs >= 2L & pcs <= cap])
  if (!length(pcs)) {
    if (cap < 2L) {
      stop("Not enough cells or genes for PCA (cap ", cap, ").", call. = FALSE)
    }
    pcs <- cap
  }
  pcs
}

.subset_obj_by_stability_ids <- function(obj, ids) {
  ids <- as.character(ids)
  if ("Barcode" %in% colnames(obj@meta.data)) {
    hit <- match(ids, as.character(obj@meta.data$Barcode))
    if (anyNA(hit)) {
      stop(
        "Stability ids missing from Barcode: ", sum(is.na(hit)),
        call. = FALSE
      )
    }
    return(obj[, hit, drop = FALSE])
  }
  missing <- setdiff(ids, colnames(obj))
  if (length(missing)) {
    stop("Stability ids missing from cell names.", call. = FALSE)
  }
  obj[, ids, drop = FALSE]
}

.subset_stability_barcodes <- function(obj) {
  if ("Barcode" %in% colnames(obj@meta.data)) {
    return(as.character(obj$Barcode))
  }
  colnames(obj)
}

.subset_reembed_prepare_pca <- function(obj, npcs) {
  if (!"Integrated_RISC" %in% Seurat::Assays(obj)) {
    stop("RISC subset re-embed needs an Integrated_RISC assay.", call. = FALSE)
  }
  # Vignette: variable genes, scale, and PCA all use the corrected assay.
  Seurat::DefaultAssay(obj) <- "Integrated_RISC"
  obj <- Seurat::FindVariableFeatures(obj, verbose = FALSE)
  feats <- Seurat::VariableFeatures(obj)
  if (length(feats) < 3L) {
    stop("Too few variable genes in Integrated_RISC.", call. = FALSE)
  }
  npcs <- min(
    as.integer(npcs)[1],
    .subset_pc_cap(ncol(obj), length(feats))
  )
  if (npcs < 2L) {
    stop("Subset has too few cells or genes for PCA.", call. = FALSE)
  }
  obj <- Seurat::ScaleData(obj, features = feats, verbose = FALSE)
  scale_feats <- rownames(Seurat::GetAssayData(
    obj, assay = "Integrated_RISC", layer = "scale.data"
  ))
  if (length(scale_feats) < 3L) {
    stop("ScaleData left too few Integrated_RISC features.", call. = FALSE)
  }
  npcs <- min(npcs, .subset_pc_cap(ncol(obj), length(scale_feats)))
  obj <- Seurat::RunPCA(
    obj,
    assay = "Integrated_RISC",
    features = scale_feats,
    npcs = npcs,
    verbose = FALSE,
    reduction.name = "pca"
  )
  obj
}

.subset_reembed_graph_cluster <- function(obj, pcs, res, cluster_name, run_umap = FALSE) {
  dims <- seq_len(as.integer(pcs)[1])
  have <- ncol(Seurat::Embeddings(obj, "pca"))
  if (max(dims) > have) {
    stop("PCA has ", have, " components; requested ", max(dims), ".", call. = FALSE)
  }
  k_use <- min(20L, ncol(obj) - 1L)
  if (k_use < 2L) {
    stop("Subset has too few cells to build a neighbor graph.", call. = FALSE)
  }
  obj <- Seurat::FindNeighbors(
    obj,
    reduction = "pca",
    dims = dims,
    k.param = k_use,
    graph.name = c("subset_nn", "subset_snn"),
    verbose = FALSE
  )
  obj <- Seurat::FindClusters(
    obj,
    graph.name = "subset_snn",
    resolution = as.numeric(res)[1],
    cluster.name = cluster_name,
    verbose = FALSE
  )
  if (!cluster_name %in% colnames(obj@meta.data)) {
    stop("FindClusters did not create column ", cluster_name, ".", call. = FALSE)
  }
  if (isTRUE(run_umap)) {
    obj <- Seurat::RunUMAP(
      obj,
      reduction = "pca",
      dims = dims,
      reduction.name = "umap",
      verbose = FALSE
    )
  }
  obj
}

.subset_reembed_cluster_df <- function(obj, pcs, res) {
  clustered <- .subset_reembed_graph_cluster(
    obj,
    pcs = pcs,
    res = res,
    cluster_name = "reembed_tmp",
    run_umap = FALSE
  )
  data.frame(
    barcode = .subset_stability_barcodes(clustered),
    cluster = as.character(clustered@meta.data$reembed_tmp),
    stringsAsFactors = FALSE
  )
}

.subset_reembed_backend <- function(obj, saving_algorithm = "qs2") {
  list(
    engine = "reembed",
    saving_algorithm = saving_algorithm,
    integrate_at_pcs = function(pcs) {
      have <- ncol(Seurat::Embeddings(obj, "pca"))
      if (as.integer(pcs)[1] > have) return(NULL)
      obj
    },
    cluster_at_integrated = function(pca_obj, pcs, res) {
      .subset_reembed_cluster_df(pca_obj, pcs, res)
    }
  )
}

.subset_reembed_bootstrap_worker <- function(
    repi,
    paramfield,
    path_obj,
    max_pcs,
    propcells.perrep,
    bootstrapresdir,
    verbose = FALSE,
    saving_algorithm = "qs2"
) {
  set_parallel_blas_threads()
  set.seed(repi)
  tryCatch({
    obj <- read_scdapp_object(path_obj, nthreads = 1L)
    sub_ids <- .stability_subsample_cells(obj@meta.data, "Code", propcells.perrep)
    obj <- .subset_obj_by_stability_ids(obj, sub_ids)
    obj <- .subset_reembed_prepare_pca(obj, max_pcs)
    backend <- .subset_reembed_backend(obj, saving_algorithm = saving_algorithm)
    .stability_sweep_save_params(
      backend,
      paramfield,
      bootstrapresdir,
      filename_prefix = paste0("Bootstrap-", repi, "."),
      verbose = verbose
    )
  }, error = function(e) {
    warning(
      "Bootstrap-", repi, " failed: ", conditionMessage(e),
      call. = FALSE
    )
  })
  invisible(NULL)
}

.subset_assign_cluster_labels <- function(obj, barcodes, clusters, cluster_name) {
  ids <- .subset_stability_barcodes(obj)
  lab <- clusters[match(ids, barcodes)]
  if (anyNA(lab)) {
    stop("Reference clustering is missing subset cells.", call. = FALSE)
  }
  lab <- .remap_clusters_by_size(lab)
  obj@meta.data[[cluster_name]] <- lab
  obj$seurat_clusters <- lab
  Seurat::SetIdent(obj, value = lab)
}

.subset_reembed_stability_sweep <- function(
    obj,
    pcs_int,
    res_int,
    sweep_maxPCs,
    sweep_res,
    numreps,
    propcells.perrep,
    stability_workernum,
    outdir,
    verbose = TRUE,
    saving_algorithm = "qs2"
) {
  .check_stability_dependencies()
  saving_algorithm <- normalize_saving_algorithm(saving_algorithm)
  n_ok <- sum(table(obj$Code) >= 2L)
  if (n_ok < 2L) {
    stop(
      "Auto PC/resolution needs at least two samples with >= 2 cells in the subset.",
      call. = FALSE
    )
  }

  cap <- .subset_pc_cap(ncol(obj))
  if (.is_auto_integration_param(pcs_int)) {
    sweep_maxPCs <- .subset_cap_pc_values(sweep_maxPCs, cap)
  } else {
    requested <- as.integer(pcs_int)[1]
    pcs_int <- .subset_cap_pc_values(requested, cap)[1]
    if (!identical(pcs_int, requested)) {
      warning(
        "pcs_int ", requested, " exceeds the subset PCA cap ", cap,
        ". Using ", pcs_int, ".",
        call. = FALSE
      )
    }
  }

  paramfield <- .build_stability_paramfield(
    sweep_maxPCs,
    sweep_res,
    pcs_int = pcs_int,
    res_int = res_int
  )
  max_pcs <- max(.stability_pcs_values_from_paramfield(paramfield))

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  rawoutsdir <- file.path(outdir, "RawOuts")
  globrefdir <- file.path(rawoutsdir, "GlobalReferenceClustering")
  bootstrapresdir <- file.path(rawoutsdir, "PerBootstrapClustering")
  dir.create(globrefdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(bootstrapresdir, recursive = TRUE, showWarnings = FALSE)

  if (verbose) message("Subset re-embed: PCA on frozen Integrated_RISC")
  obj_pca <- .subset_reembed_prepare_pca(obj, max_pcs)
  backend <- .subset_reembed_backend(obj_pca, saving_algorithm = saving_algorithm)
  .stability_sweep_save_params(
    backend,
    paramfield,
    globrefdir,
    verbose = verbose
  )

  path_obj <- scdapp_object_path(file.path(rawoutsdir, "TEMP_SUBSET"), saving_algorithm)
  write_scdapp_object(obj, path_obj, saving_algorithm = saving_algorithm, nthreads = 1L)

  reps <- seq_len(as.integer(numreps)[1])
  n_workers <- .stability_effective_workers(length(reps), stability_workernum)
  if (verbose) {
    message(
      "Subset re-embed bootstrap: ", length(reps),
      " replicates, ", n_workers, " workers"
    )
  }
  if (n_workers <= 1L) {
    for (repi in reps) {
      .subset_reembed_bootstrap_worker(
        repi,
        paramfield,
        path_obj,
        max_pcs,
        propcells.perrep,
        bootstrapresdir,
        verbose = verbose,
        saving_algorithm = saving_algorithm
      )
    }
  } else {
    cl <- parallel::makeCluster(n_workers, outfile = "")
    doParallel::registerDoParallel(cl)
    .stability_cluster_register_workers(cl)
    parallel::clusterExport(
      cl,
      c(
        "paramfield", "path_obj", "max_pcs", "propcells.perrep",
        "bootstrapresdir", "verbose", "saving_algorithm"
      ),
      envir = environment()
    )
    foreach::foreach(
      repi = reps,
      .packages = .stability_foreach_packages(c("Seurat", "scDAPP")),
      .verbose = verbose
    ) %dopar% {
      scDAPP:::.subset_reembed_bootstrap_worker(
        repi,
        paramfield,
        path_obj,
        max_pcs,
        propcells.perrep,
        bootstrapresdir,
        verbose,
        saving_algorithm
      )
    }
    parallel::stopCluster(cl)
  }

  globref_l <- lapply(paramfield, function(params_i) {
    f <- scdapp_object_path(file.path(globrefdir, params_i), saving_algorithm)
    if (file.exists(f)) read_scdapp_object(f, nthreads = 1L) else NULL
  })
  names(globref_l) <- paramfield
  globref_l <- globref_l[lengths(globref_l) > 0L]
  paramfield <- paramfield[paramfield %in% names(globref_l)]
  if (!length(paramfield)) {
    stop("All subset re-embed reference clusterings failed.", call. = FALSE)
  }

  bootstrap_index <- do.call(rbind, lapply(reps, function(repi) {
    data.frame(
      Filename = paste0("Bootstrap-", repi, ".", paramfield),
      BootstrapID = paste0("Bootstrap-", repi),
      ParamID = paramfield,
      stringsAsFactors = FALSE
    )
  }))
  repres_flat <- lapply(bootstrap_index$Filename, function(fn) {
    f <- scdapp_object_path(file.path(bootstrapresdir, fn), saving_algorithm)
    if (file.exists(f)) read_scdapp_object(f, nthreads = 1L) else NULL
  })
  names(repres_flat) <- bootstrap_index$Filename
  repres_flat <- repres_flat[lengths(repres_flat) > 0L]
  if (!length(repres_flat)) {
    stop(
      "No subset bootstrap clusterings were saved. The subset may be too small for this PC grid.",
      call. = FALSE
    )
  }
  bsindex <- bootstrap_index[bootstrap_index$Filename %in% names(repres_flat), , drop = FALSE]

  if (verbose) message("Scoring subset re-embed stability")
  fulloutlist <- .compute_stability_metrics(
    globref_l,
    repres_flat,
    bsindex,
    paramfield,
    verbose = verbose,
    stability_workernum = stability_workernum
  )
  table_names <- c(
    "perparam_meanscores",
    "perbootstrap_perparam_scores",
    "perclust_jaccard_mean_acrossbootstraps",
    "perclust_jaccard_per_bootstrap_long"
  )
  for (nm in table_names) {
    utils::write.csv(
      fulloutlist[[nm]],
      file.path(outdir, paste0(nm, ".csv")),
      row.names = FALSE
    )
  }
  best <- .pick_best_stability_params(fulloutlist$perparam_meanscores)
  resolved_pcs <- if (.is_auto_integration_param(pcs_int)) {
    best$pcs_int
  } else {
    as.integer(pcs_int)[1]
  }
  resolved_res <- if (.is_auto_integration_param(res_int)) {
    best$res_int
  } else {
    as.numeric(res_int)[1]
  }
  selected_params_i <- fulloutlist$perparam_meanscores$params_i[
    fulloutlist$perparam_meanscores$MAX == "*"
  ][1]
  # When one of pcs/res was fixed, the starred row is still the best joint score.
  # Prefer the row that matches the resolved pair when it exists.
  resolved_id <- paste0("PCs_1-", resolved_pcs, ".res_", resolved_res)
  if (resolved_id %in% fulloutlist$perparam_meanscores$params_i) {
    selected_params_i <- resolved_id
  }
  ref <- globref_l[[selected_params_i]]
  if (is.null(ref)) {
    stop("Missing reference clustering for ", selected_params_i, ".", call. = FALSE)
  }

  plot_dir <- file.path(outdir, "plots")
  fulloutlist$plots <- tryCatch(
    integration_stability_plots_module(
      stability = fulloutlist,
      outdir = plot_dir,
      selected_params_i = selected_params_i,
      save = TRUE
    ),
    error = function(e) {
      warning("Subset stability plots failed: ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
  if (dir.exists(rawoutsdir)) unlink(rawoutsdir, recursive = TRUE)

  list(
    obj_pca = obj_pca,
    pcs_int = resolved_pcs,
    res_int = resolved_res,
    selected_params_i = selected_params_i,
    reference_labels = ref,
    stability = fulloutlist,
    stability_dir = outdir
  )
}

.subset_reembed_finish <- function(
    obj,
    pcs,
    res,
    cluster_name,
    reference_labels = NULL,
    pca_ready = FALSE
) {
  if (!isTRUE(pca_ready)) obj <- .subset_reembed_prepare_pca(obj, pcs)
  if (is.null(reference_labels)) {
    obj <- .subset_reembed_graph_cluster(
      obj,
      pcs = pcs,
      res = res,
      cluster_name = cluster_name,
      run_umap = FALSE
    )
    obj@meta.data[[cluster_name]] <- .remap_clusters_by_size(obj@meta.data[[cluster_name]])
    obj$seurat_clusters <- obj@meta.data[[cluster_name]]
    obj <- Seurat::SetIdent(obj, value = obj@meta.data[[cluster_name]])
  } else {
    obj <- .subset_assign_cluster_labels(
      obj,
      barcodes = reference_labels$barcode,
      clusters = reference_labels$cluster,
      cluster_name = cluster_name
    )
  }
  n_pca <- ncol(Seurat::Embeddings(obj, "pca"))
  dims <- seq_len(min(as.integer(pcs)[1], n_pca))
  n_neighbors <- min(30L, ncol(obj) - 1L)
  obj <- Seurat::RunUMAP(
    obj,
    reduction = "pca",
    dims = dims,
    n.neighbors = n_neighbors,
    reduction.name = "umap",
    verbose = FALSE
  )
  Seurat::DefaultAssay(obj) <- "Integrated_RISC"
  obj
}

.subset_rna_counts <- function(obj) {
  mat <- tryCatch(
    Seurat::GetAssayData(obj, assay = "RNA", layer = "counts"),
    error = function(e) NULL
  )
  if (is.null(mat) || nrow(mat) == 0L || ncol(mat) == 0L) {
    stop("Subset re-integration needs RNA counts.", call. = FALSE)
  }
  mat
}

.subset_split_samples <- function(obj, sample_metadata) {
  codes <- as.character(sample_metadata$Code)
  sobjlist <- lapply(codes, function(code) {
    cells <- colnames(obj)[as.character(obj$Code) == code]
    if (length(cells) < 2L) return(NULL)
    sub <- obj[, cells, drop = FALSE]
    counts <- .subset_rna_counts(sub)
    md <- sub@meta.data
    rownames(md) <- colnames(counts)
    Seurat::CreateSeuratObject(counts = counts, meta.data = md, assay = "RNA")
  })
  names(sobjlist) <- codes
  sobjlist <- sobjlist[!vapply(sobjlist, is.null, logical(1))]
  md <- sample_metadata[as.character(sample_metadata$Code) %in% names(sobjlist), , drop = FALSE]
  if (nrow(md) < 2L) {
    stop(
      "Re-integration needs at least two samples with >= 2 cells in the subset.",
      call. = FALSE
    )
  }
  list(sobjlist = sobjlist, sample_metadata = md)
}

.subset_reintegrate_pc_grid <- function(obj, pcs_int, sweep_maxPCs) {
  min_cells <- min(as.integer(table(obj$Code)))
  cap <- .stability_cap_pcs(max(c(as.integer(sweep_maxPCs), 2L)), min_cells)
  if (.is_auto_integration_param(pcs_int)) {
    pcs <- .subset_cap_pc_values(sweep_maxPCs, cap)
    return(list(pcs_int = "auto", sweep_maxPCs = pcs, cap = cap))
  }
  requested <- as.integer(pcs_int)[1]
  used <- .subset_cap_pc_values(requested, .stability_cap_pcs(requested, min_cells))[1]
  if (!identical(used, requested)) {
    warning(
      "pcs_int ", requested, " exceeds the smallest subset sample. Using ", used, ".",
      call. = FALSE
    )
  }
  list(pcs_int = used, sweep_maxPCs = used, cap = cap)
}

.subset_label_transfer_summary <- function(obj, cluster_column) {
  if (!"predictions" %in% Seurat::Assays(obj)) return(NULL)
  pred <- tryCatch(
    Seurat::GetAssayData(obj, assay = "predictions", layer = "data"),
    error = function(e) NULL
  )
  if (is.null(pred) || !ncol(pred)) return(NULL)
  drop <- tolower(rownames(pred)) %in% c("max", "prediction.score.max")
  if (any(drop)) pred <- pred[!drop, , drop = FALSE]
  if (!nrow(pred)) return(NULL)
  md <- obj@meta.data
  cl <- as.character(md[[cluster_column]])
  cells <- intersect(colnames(pred), rownames(md))
  if (!length(cells)) return(NULL)
  pred <- pred[, cells, drop = FALSE]
  cl <- cl[match(cells, rownames(md))]
  clusters <- unique(cl)
  means <- sapply(clusters, function(cv) {
    Matrix::rowMeans(pred[, cl == cv, drop = FALSE])
  })
  if (is.null(dim(means))) {
    means <- matrix(means, ncol = 1L, dimnames = list(rownames(pred), clusters))
  }
  colnames(means) <- clusters
  top <- apply(means, 2, function(x) {
    data.frame(
      top_label = names(x)[which.max(x)],
      score = max(x),
      stringsAsFactors = FALSE
    )
  })
  top <- dplyr::bind_rows(top)
  top$cluster <- colnames(means)
  top$label <- top$top_label
  top$label[top$score < 0.3] <- NA_character_
  mid <- !is.na(top$score) & top$score >= 0.3 & top$score < 0.5
  top$label[mid] <- paste0(top$top_label[mid], "_putative")
  top <- top[, c("cluster", "top_label", "score", "label"), drop = FALSE]
  list(means = means, summary = top)
}

.subset_de_top_table <- function(de_results, n = 10L) {
  if (is.null(de_results) || !is.data.frame(de_results) || !nrow(de_results)) {
    return(de_results)
  }
  lfc <- if ("logFC" %in% colnames(de_results)) {
    "logFC"
  } else if ("avg_log2FC" %in% colnames(de_results)) {
    "avg_log2FC"
  } else {
    NULL
  }
  if (is.null(lfc) || !"cluster" %in% colnames(de_results)) {
    return(utils::head(de_results, n))
  }
  key <- if ("label" %in% colnames(de_results)) {
    paste(de_results$label, de_results$cluster, sep = "||")
  } else {
    as.character(de_results$cluster)
  }
  pieces <- lapply(split(de_results, key), function(df) {
    df <- df[order(-abs(df[[lfc]])), , drop = FALSE]
    utils::head(df, n)
  })
  out <- dplyr::bind_rows(pieces)
  rownames(out) <- NULL
  out
}

.subset_dimplots <- function(obj, cluster_column, reduction) {
  list(
    clusters = Seurat::DimPlot(
      obj, reduction = reduction, group.by = cluster_column, label = TRUE, repel = TRUE
    ),
    condition = Seurat::DimPlot(
      obj, reduction = reduction, group.by = "Condition"
    ),
    code = Seurat::DimPlot(
      obj, reduction = reduction, group.by = "Code"
    )
  )
}

.subset_record_misc <- function(obj, info) {
  prev <- obj@misc$scDAPP
  if (is.null(prev)) prev <- list()
  prev$subset_module <- info
  obj@misc$scDAPP <- prev
  obj
}

#' Subset an integrated object, recluster, and run comparative analysis
#'
#' Keep selected clusters or cells from an integrated Seurat object, recluster
#' that subset, and run the cross-condition analyses on the new clusters.
#' Composition uses propeller when \code{Pseudobulk_mode} is \code{TRUE} and a
#' chi-square test otherwise. Differential expression uses \code{DE_test}.
#' Pathway GSEA runs by default. Milo runs when pseudobulk mode is on. ORA runs
#' only when \code{run_ora = TRUE}. If a \code{predictions} assay is present,
#' label-transfer scores are summarized on the new clusters. An HTML report is
#' written to \code{outdir} unless \code{render_html = FALSE}.
#'
#' Reclustering follows the parent integration method. RISC, by default, keeps
#' the existing \code{Integrated_RISC} assay and recomputes variable genes,
#' scaling, and PCA on that assay. Harmony, CCA, RPCA, and SCT re-integrate the
#' subset with \code{run_integration()}. \code{reintegrate = TRUE} re-integrates
#' a RISC subset as well. \code{pcs_int} or \code{res_int = "auto"} selects the
#' PC count and Louvain resolution by bootstrap.
#'
#' @param sobjint Integrated Seurat object.
#' @param sample_metadata Sample table. \code{Code} and \code{Condition} are
#'   used when present. Otherwise \code{sample_column} and \code{condition_column}
#'   are copied onto those names.
#' @param comps Comparison table with \code{c0} and \code{c1}. See \code{normalize_comps()}.
#' @param outdir Output directory for the subset run and HTML report.
#' @param clusters Cluster ids to keep, matched to \code{cluster_column}.
#' @param cells Cell names or \code{Barcode} values to keep. Do not pass with \code{clusters}.
#' @param cluster_column Metadata column used to interpret \code{clusters}.
#'   Default \code{"seurat_clusters"}.
#' @param sample_column Sample id column on the object. Default \code{"Code"}.
#' @param condition_column Condition column on the object. Default \code{"Condition"}.
#' @param integration_method Parent integration method. Default reads
#'   \code{sobjint@misc$scDAPP$integration_method}.
#' @param Pseudobulk_mode If \code{TRUE}, composition uses propeller, DE uses a
#'   pseudobulk test on \code{RNA}/\code{counts}, and Milo runs. If \code{FALSE},
#'   composition is a chi-square test and DE is a single-cell test on the
#'   integration method's normalized assay.
#' @param pcs_int PCs for the subset clustering, or \code{"auto"}.
#' @param res_int Louvain resolution, or \code{"auto"}.
#' @param reintegrate If \code{NULL} (default), RISC re-embeds and every other
#'   method re-integrates. \code{TRUE} re-integrates RISC as well.
#' @param species Species for \code{msigdbr} when \code{run_pathways = TRUE}.
#' @param run_pathways Run GSEA. Default \code{TRUE}.
#' @param run_ora Run ORA. Default \code{FALSE}.
#' @param run_milo Run Milo. Default is \code{Pseudobulk_mode}.
#' @param DE_test Differential expression test, passed to
#'   \code{de_across_conditions_module()}. \code{NULL} uses that module's
#'   default: \code{EdgeR-LRT} when \code{Pseudobulk_mode} is \code{TRUE},
#'   \code{wilcox} when \code{FALSE}. Pseudobulk choices are \code{DESeq2},
#'   \code{DESeq2-LRT}, \code{EdgeR}, \code{EdgeR-LRT}, \code{EdgeR-QLF}, and
#'   \code{Dream} (\code{Dream} requires \code{variancePartition} and a
#'   \code{(1|var)} term in \code{comps$formula}). Single-cell choices are any
#'   \code{test.use} value from \code{Seurat::FindMarkers} except \code{roc}.
#' @param stability_numreps Bootstrap replicates when PC or resolution is \code{"auto"}.
#' @param stability_sweep_maxPCs PC grid for auto selection.
#' @param stability_sweep_res Resolution grid for auto selection.
#' @param stability_propcells.perrep Cell fraction per bootstrap replicate.
#' @param workernum Workers for integration, DE, and pathways.
#' @param stability_workernum Workers across bootstrap replicates. \code{NULL} uses \code{workernum}.
#' @param render_html Write \code{subset_module_report.html} in \code{outdir}.
#' @param verbose Print progress.
#' @param saving_algorithm \code{"qs2"} or \code{"rds"}.
#' @return List with the subset object, new cluster column, stability results,
#'   and comparative module results.
#' @export
subset_module <- function(
    sobjint,
    sample_metadata,
    comps,
    outdir,
    clusters = NULL,
    cells = NULL,
    cluster_column = "seurat_clusters",
    sample_column = "Code",
    condition_column = "Condition",
    integration_method = NULL,
    Pseudobulk_mode = TRUE,
    pcs_int = "auto",
    res_int = "auto",
    reintegrate = NULL,
    species = "Homo sapiens",
    run_pathways = TRUE,
    run_ora = FALSE,
    run_milo = NULL,
    DE_test = NULL,
    stability_numreps = 50L,
    stability_sweep_maxPCs = c(5, 10, 15, 20, 25, 30),
    stability_sweep_res = seq(0.1, 1.5, by = 0.2),
    stability_propcells.perrep = 0.8,
    workernum = 1L,
    stability_workernum = NULL,
    render_html = TRUE,
    verbose = TRUE,
    saving_algorithm = "qs2"
) {
  if (missing(sobjint)) stop("Provide an integrated Seurat object.", call. = FALSE)
  if (missing(outdir) || !nzchar(as.character(outdir)[1])) {
    stop("Provide outdir.", call. = FALSE)
  }
  saving_algorithm <- normalize_saving_algorithm(saving_algorithm)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  workernum <- as.integer(workernum)[1]
  if (is.na(workernum) || workernum < 1L) workernum <- 1L
  stability_workernum <- resolve_stability_workernum(workernum, stability_workernum)

  method <- .subset_parent_method(sobjint, integration_method)
  if (is.null(reintegrate)) reintegrate <- !identical(method, "RISC")
  reintegrate <- isTRUE(reintegrate)
  if (!reintegrate && !identical(method, "RISC")) {
    stop(
      "reintegrate = FALSE is only supported for RISC, which stores Integrated_RISC.",
      call. = FALSE
    )
  }
  if (is.null(run_milo)) run_milo <- isTRUE(Pseudobulk_mode)
  if (isTRUE(run_milo) && !isTRUE(Pseudobulk_mode)) {
    warning("Milo is skipped when Pseudobulk_mode is FALSE.", call. = FALSE)
    run_milo <- FALSE
  }

  sample_metadata <- .subset_standardize_sample_metadata(
    sample_metadata, sample_column, condition_column
  )
  obj <- .subset_standardize_object_metadata(sobjint, sample_column, condition_column)
  obj <- .subset_apply_condition_from_metadata(obj, sample_metadata)
  if (
    !is.null(cluster_column) &&
    cluster_column %in% colnames(obj@meta.data)
  ) {
    obj@meta.data[[paste0("parent_", cluster_column)]] <- as.character(
      obj@meta.data[[cluster_column]]
    )
  }
  keep <- .subset_select_cells(obj, clusters, cells, cluster_column)
  obj <- obj[, keep, drop = FALSE]
  filtered <- .subset_filter_metadata(obj, sample_metadata)
  obj <- filtered$obj
  sample_metadata <- filtered$sample_metadata
  comps <- .subset_check_comps(sample_metadata, comps)
  if (verbose) {
    message(
      "Subset: ", ncol(obj), " cells, ",
      length(unique(as.character(obj$Code))), " samples, method ", method,
      if (reintegrate) " (re-integrate)" else " (re-embed Integrated_RISC)"
    )
  }

  stability <- NULL
  stability_dir <- NULL
  selected_params_i <- NULL
  pcs_use <- pcs_int
  res_use <- res_int
  auto <- .integration_param_needs_auto(pcs_int, res_int)

  if (!reintegrate) {
    if (auto) {
      stability_dir <- file.path(outdir, "cluster_stability")
      swept <- .subset_reembed_stability_sweep(
        obj,
        pcs_int = pcs_int,
        res_int = res_int,
        sweep_maxPCs = stability_sweep_maxPCs,
        sweep_res = stability_sweep_res,
        numreps = stability_numreps,
        propcells.perrep = stability_propcells.perrep,
        stability_workernum = stability_workernum,
        outdir = stability_dir,
        verbose = verbose,
        saving_algorithm = saving_algorithm
      )
      pcs_use <- swept$pcs_int
      res_use <- swept$res_int
      stability <- swept$stability
      selected_params_i <- swept$selected_params_i
      cluster_col <- paste0("RISC_subset_npc", pcs_use, "_res", res_use)
      obj <- .subset_reembed_finish(
        swept$obj_pca,
        pcs = pcs_use,
        res = res_use,
        cluster_name = cluster_col,
        reference_labels = swept$reference_labels,
        pca_ready = TRUE
      )
    } else {
      pcs_use <- .subset_cap_pc_values(as.integer(pcs_int)[1], .subset_pc_cap(ncol(obj)))[1]
      res_use <- as.numeric(res_int)[1]
      cluster_col <- paste0("RISC_subset_npc", pcs_use, "_res", res_use)
      obj <- .subset_reembed_finish(
        obj,
        pcs = pcs_use,
        res = res_use,
        cluster_name = cluster_col
      )
    }
    cfg <- resolve_integration_config(method, pcs_use, res_use)
  } else {
    split <- .subset_split_samples(obj, sample_metadata)
    sample_metadata <- split$sample_metadata
    comps <- .subset_check_comps(sample_metadata, comps)
    pc_grid <- .subset_reintegrate_pc_grid(obj, pcs_int, stability_sweep_maxPCs)
    outdir_int <- file.path(outdir, "reintegration")
    outdir_indi <- file.path(outdir, "individual")
    dir.create(outdir_indi, recursive = TRUE, showWarnings = FALSE)
    if (verbose) message("Re-integrating subset with ", method)
    prep <- prepare_integration_inputs(
      sobjlist = split$sobjlist,
      sample_metadata = sample_metadata,
      outdir_indi = outdir_indi,
      outdir_int = outdir_int,
      input_seurat_obj = TRUE,
      saving_algorithm = saving_algorithm,
      nthreads = workernum
    )
    stability_dir <- if (auto) file.path(outdir_int, "cluster_stability") else NULL
    int_res <- run_integration(
      integration_method = method,
      sample_metadata = sample_metadata,
      mdlist = prep$mdlist,
      tmpobjdir = prep$tmpobjdir,
      outdir_int = outdir_int,
      outdir_indi = outdir_indi,
      pcs_int = pc_grid$pcs_int,
      res_int = res_int,
      workernum = workernum,
      stability_workernum = stability_workernum,
      input_seurat_obj = TRUE,
      stability_outdir = stability_dir,
      stability_numreps = stability_numreps,
      stability_sweep_maxPCs = pc_grid$sweep_maxPCs,
      stability_sweep_res = stability_sweep_res,
      stability_propcells.perrep = stability_propcells.perrep,
      pipeline_outdir = outdir,
      saving_algorithm = saving_algorithm
    )
    obj <- int_res$sobjint
    cfg <- int_res$config
    cluster_col <- int_res$int_clust_lab
    pcs_use <- cfg$pcs_int
    res_use <- cfg$res_int
    stability <- int_res$stability
    selected_params_i <- if (!is.null(int_res$auto_selection)) {
      int_res$auto_selection$selected_params_i
    } else {
      NULL
    }
  }

  obj <- .subset_record_misc(obj, list(
    parent_method = method,
    reintegrate = reintegrate,
    embedding = if (reintegrate) "reintegrate" else "reembed",
    parent_cluster_column = cluster_column,
    clusters = if (is.null(clusters)) character() else as.character(clusters),
    pcs_int = pcs_use,
    res_int = res_use,
    cluster_column = cluster_col
  ))
  obj_dir <- file.path(outdir, "data_objects")
  dir.create(obj_dir, recursive = TRUE, showWarnings = FALSE)
  write_scdapp_object(
    obj,
    file.path(obj_dir, "subset_object"),
    saving_algorithm = saving_algorithm,
    nthreads = workernum
  )

  label_transfer <- .subset_label_transfer_summary(obj, cluster_col)
  if (!is.null(label_transfer)) {
    utils::write.csv(
      label_transfer$summary,
      file.path(outdir, "label_transfer_by_subset_cluster.csv"),
      row.names = FALSE
    )
  }

  comparative <- .run_comparative_stack(
    obj = obj,
    sample_metadata = sample_metadata,
    comps = comps,
    outdir = outdir,
    grouping_variable = cluster_col,
    cfg = cfg,
    Pseudobulk_mode = Pseudobulk_mode,
    DE_test = DE_test,
    species = species,
    run_pathways = run_pathways,
    run_ora = run_ora,
    run_milo = run_milo,
    milo_d = as.integer(pcs_use)[1],
    workernum = workernum,
    verbose = verbose,
    saving_algorithm = saving_algorithm
  )
  composition <- comparative$composition
  de_results <- comparative$de_results
  pathways <- comparative$pathways
  ora <- comparative$ora
  milo <- comparative$milo

  reduction <- if (!reintegrate) "umap" else cfg$umap_reduction
  dimplots <- tryCatch(
    .subset_dimplots(obj, cluster_col, reduction),
    error = function(e) {
      warning("DimPlots failed: ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )

  html_path <- NULL
  if (isTRUE(render_html)) {
    if (verbose) message("Rendering subset HTML report")
    report <- list(
      method = method,
      reintegrate = reintegrate,
      embedding = if (reintegrate) "reintegrate" else "reembed",
      n_cells = ncol(obj),
      n_samples = nrow(sample_metadata),
      parent_cluster_column = cluster_column,
      clusters = if (is.null(clusters)) character() else as.character(clusters),
      n_cells_requested = if (is.null(cells)) NA_integer_ else length(cells),
      pcs_int = pcs_use,
      res_int = res_use,
      cluster_column = cluster_col,
      pseudobulk = isTRUE(Pseudobulk_mode),
      selected_params_i = selected_params_i,
      dimplots = dimplots,
      stability = stability,
      label_transfer = label_transfer,
      composition = composition,
      de_top = .subset_de_top_table(de_results, 10L),
      pathways = pathways,
      ora = ora,
      milo = milo
    )
    report_rds <- file.path(outdir, "subset_module_report_data.rds")
    saveRDS(report, report_rds)
    rmd <- system.file("rmd", "subset_module_report.Rmd", package = "scDAPP")
    if (!nzchar(rmd)) {
      stop("Could not find inst/rmd/subset_module_report.Rmd.", call. = FALSE)
    }
    html_path <- file.path(outdir, "subset_module_report.html")
    rmarkdown::render(
      input = rmd,
      output_file = basename(html_path),
      output_dir = outdir,
      intermediates_dir = file.path(outdir, ".knit"),
      params = list(report_rds = normalizePath(report_rds, mustWork = TRUE)),
      envir = new.env(parent = globalenv()),
      quiet = !isTRUE(verbose)
    )
    unlink(file.path(outdir, ".knit"), recursive = TRUE)
  }

  list(
    sobj = obj,
    cluster_column = cluster_col,
    integration_method = method,
    reintegrate = reintegrate,
    embedding = if (reintegrate) "reintegrate" else "reembed",
    pcs_int = pcs_use,
    res_int = res_use,
    stability = stability,
    stability_dir = stability_dir,
    selected_params_i = selected_params_i,
    composition = composition,
    de_results = de_results,
    pathways = pathways,
    ora = ora,
    milo = milo,
    label_transfer = label_transfer,
    html = html_path,
    outdir = outdir
  )
}
