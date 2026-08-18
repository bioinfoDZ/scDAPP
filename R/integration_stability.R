# ---- Cluster stability: helpers ----

#' Build PC x resolution grid labels for stability sweep (params_i ids).
#' @keywords internal
.build_stability_paramfield <- function(
    sweep_maxPCs,
    sweep_res,
    pcs_int = NULL,
    res_int = NULL
) {
  if (.is_auto_integration_param(pcs_int)) {
    pc_vals <- as.integer(sweep_maxPCs)
  } else {
    pc_vals <- as.integer(pcs_int)[1]
  }
  if (.is_auto_integration_param(res_int)) {
    res_vals <- as.numeric(sweep_res)
  } else {
    res_vals <- as.numeric(res_int)[1]
  }
  pc_vals <- unique(pc_vals)
  res_vals <- unique(res_vals)
  paramfield <- unlist(lapply(pc_vals, function(pci) {
    paste0("PCs_1-", pci, ".res_", res_vals)
  }), use.names = FALSE)
  names(paramfield) <- paramfield
  paramfield
}

#' Unique PC values from stability param ids
#' @keywords internal
.stability_pcs_values_from_paramfield <- function(paramfield) {
  sort(unique(vapply(
    paramfield,
    function(p) .parse_stability_param_id(p)$pcs_int,
    integer(1)
  )))
}

#' Param ids in \code{paramfield} matching a PC count
#' @keywords internal
.stability_params_with_pcs <- function(paramfield, pcs_int) {
  pcs_int <- as.integer(pcs_int)[1]
  paramfield[vapply(
    paramfield,
    function(p) .parse_stability_param_id(p)$pcs_int == pcs_int,
    logical(1)
  )]
}

#' Save integrate+cluster results grouped by PC (integrate once per PC, cluster all res)
#' @keywords internal
.stability_sweep_save_params <- function(
    backend,
    params_to_run,
    savedir,
    filename_prefix = "",
    verbose = FALSE
) {
  for (pcs_k in .stability_pcs_values_from_paramfield(params_to_run)) {
    params_k <- .stability_params_with_pcs(params_to_run, pcs_k)
    int_state <- NULL
    try({
      int_state <- backend$integrate_at_pcs(pcs_k)
    }, silent = !verbose)
    if (is.null(int_state)) next
    for (params_i in params_k) {
      out_path <- file.path(savedir, paste0(filename_prefix, params_i, ".rds"))
      if (file.exists(out_path)) next
      parsed <- .parse_stability_param_id(params_i)
      outdf <- NULL
      try({
        outdf <- backend$cluster_at_integrated(
          int_state,
          parsed$pcs_int,
          parsed$res_int
        )
      }, silent = !verbose)
      if (!is.null(outdf)) {
        saveRDS(outdf, out_path)
      }
    }
    rm(int_state)
    invisible(gc(full = TRUE, reset = FALSE, verbose = FALSE))
  }
  invisible(NULL)
}

#' Parse PCs_1-N.res_R stability grid id into pcs_int and res_int integers.
#' @keywords internal
.parse_stability_param_id <- function(params_i) {
  m <- regexpr("^PCs_1-([0-9]+)[.]res_([0-9.]+)$", params_i, perl = TRUE)
  if (m < 0L) {
    stop("Invalid stability param id: ", params_i, call. = FALSE)
  }
  g <- regmatches(params_i, regexec("^PCs_1-([0-9]+)[.]res_([0-9.]+)$", params_i, perl = TRUE))[[1]]
  list(
    pcs_int = as.integer(g[2]),
    res_int = as.numeric(g[3])
  )
}

#' Temporary cluster column name during stability sweep for one grid point.
#' @keywords internal
.stability_sweep_cluster_name <- function(params_i) {
  paste0("stability_", params_i)
}

#' Extract barcode/cluster data.frame from stability cluster column for metrics.
#' @keywords internal
.stability_extract_cluster_df <- function(meta, params_i) {
  col_pat <- .stability_sweep_cluster_name(params_i)
  hit <- grep(col_pat, colnames(meta), fixed = TRUE)
  if (!length(hit)) {
  stop("Missing stability cluster column matching ", col_pat, call. = FALSE)
  }
  data.frame(
    barcode = rownames(meta),
    cluster = as.vector(meta[, hit[1], drop = TRUE]),
    stringsAsFactors = FALSE
  )
}

#' Select best PC/res from stability summary table
#' @keywords internal
.pick_best_stability_params <- function(perparam_meanscores) {
  if (!nrow(perparam_meanscores)) {
    stop("Stability sweep produced no parameter scores.", call. = FALSE)
  }
  perparam_meanscores$combinedscore <-
    0.5 * perparam_meanscores$ARI_mean +
    0.5 * perparam_meanscores$Jaccard_mean_of_clustermeans
  best_row <- perparam_meanscores[which.max(perparam_meanscores$combinedscore), , drop = FALSE]
  parsed <- .parse_stability_param_id(best_row$params_i)
  c(parsed, list(perparam_meanscores = perparam_meanscores))
}

#' Fail fast if mclust is missing when stability auto-tuning is requested.
#' @keywords internal
.check_stability_dependencies <- function() {
  if (!requireNamespace("mclust", quietly = TRUE)) {
    stop(
      "Cluster stability (pcs_int/res_int = 'auto') requires the mclust package. ",
      "Install with install.packages('mclust').",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Cap PC count by smallest sample layer size minus one (Seurat PCA limit).
#' @keywords internal
.stability_cap_pcs <- function(npcs, n_cells_per_layer_min) {
  max(1L, min(as.integer(npcs)[1], as.integer(n_cells_per_layer_min) - 1L))
}

#' Strip integrated assays/reductions; keep counts-only object for stability workers.
#' @keywords internal
.stability_downsize_seurat_counts <- function(sobj, assay = "RNA") {
  if (!assay %in% Seurat::Assays(sobj)) {
    assay <- Seurat::DefaultAssay(sobj)
  }
  sobj <- .join_assay_layers_if_needed(sobj, assay)
  counts <- Seurat::GetAssayData(sobj, assay = assay, layer = "counts")
  Seurat::CreateSeuratObject(
    counts = counts,
    assay = assay,
    meta.data = sobj@meta.data
  )
}

#' Subsample cells for bootstrap (stratified by sample); keep samples with >= 2 cells
#' @keywords internal
.stability_subsample_cells <- function(meta, samplecolumn, propcells.perrep) {
  use_barcode <- "Barcode" %in% colnames(meta)
  cell_ids <- function(rows) {
    if (use_barcode) as.character(meta$Barcode[rows]) else rownames(meta)[rows]
  }
  codes <- unique(as.character(meta[[samplecolumn]]))
  subcells <- unlist(lapply(codes, function(code) {
    idx <- which(as.character(meta[[samplecolumn]]) == code)
    n <- length(idx)
    if (n < 2L) {
      return(character(0))
    }
    size <- max(2L, round(n * propcells.perrep))
    size <- min(size, n)
    cell_ids(sample(idx, size = size, replace = FALSE))
  }), use.names = FALSE)
  subcells <- unique(subcells)
  if (length(subcells) < 2L) {
    stop(
      "Bootstrap subsample has too few cells; increase propcells.perrep or num cells.",
      call. = FALSE
    )
  }
  if (use_barcode) {
    submd <- meta[as.character(meta$Barcode) %in% subcells, , drop = FALSE]
  } else {
    submd <- meta[subcells, , drop = FALSE]
  }
  codetab <- table(submd[[samplecolumn]])
  if (any(codetab < 2L)) {
    drop_codes <- names(codetab)[codetab < 2L]
    submd <- submd[!(submd[[samplecolumn]] %in% drop_codes), , drop = FALSE]
  }
  if (length(unique(submd[[samplecolumn]])) < 2L) {
    stop(
      "Bootstrap subsample must retain at least two samples with >= 2 cells each.",
      call. = FALSE
    )
  }
  if (use_barcode) {
    return(as.character(submd$Barcode))
  }
  rownames(submd)
}

#' Append worker stdout/stderr to stability sweep log file during foreach jobs.
#' @keywords internal
.stability_write_worker_log <- function(logfile, expr) {
  con <- file(logfile, open = "a")
  sink(con)
  sink(con, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(con)
  }, add = TRUE)
  force(expr)
}

#' Limit BLAS threads inside a stability sweep worker process.
#' @keywords internal
.stability_limit_threads_worker <- function() {
  scDAPP::set_parallel_blas_threads()
}

#' Register doParallel cluster workers with scDAPP (dev load_all or library)
#' @keywords internal
.stability_cluster_register_workers <- function(cl) {
  pkg_root <- getOption("scDAPP.dev_pkg_root", NULL)
  if (!is.null(pkg_root) && nzchar(pkg_root)) {
    parallel::clusterExport(cl, "pkg_root", envir = environment())
    parallel::clusterEvalQ(cl, {
      loaded <- FALSE
      if (requireNamespace("pkgload", quietly = TRUE)) {
        pkgload::load_all(pkg_root, quiet = TRUE)
        loaded <- TRUE
      } else if (requireNamespace("devtools", quietly = TRUE)) {
        devtools::load_all(pkg_root, quiet = TRUE)
        loaded <- TRUE
      }
      if (!loaded) {
        stop(
          "Option scDAPP.dev_pkg_root is set but pkgload/devtools is unavailable on a worker.",
          call. = FALSE
        )
      }
    })
  }
  parallel::clusterEvalQ(cl, {
    if (!"package:scDAPP" %in% search() &&
        requireNamespace("scDAPP", quietly = TRUE)) {
      suppressPackageStartupMessages(library(scDAPP))
    }
    scDAPP::set_parallel_blas_threads()
  })
  invisible(cl)
}

#' Foreach .packages vector; omit scDAPP when dev pkg root is set (already loaded)
#' @keywords internal
.stability_foreach_packages <- function(extra = "scDAPP") {
  pkgs <- unique(c(extra, "scDAPP"))
  pkg_root <- getOption("scDAPP.dev_pkg_root", NULL)
  if (!is.null(pkg_root) && nzchar(pkg_root)) {
    pkgs <- setdiff(pkgs, "scDAPP")
  }
  pkgs
}

#' Resolve parallel workers for the stability sweep
#' @param workernum Pipeline worker count (final RISC integration).
#' @param stability_workernum Sweep-only workers; \code{NULL} uses \code{workernum}.
#' @keywords internal
resolve_stability_workernum <- function(workernum, stability_workernum = NULL) {
  if (is.null(stability_workernum)) {
    return(as.integer(workernum)[1])
  }
  as.integer(stability_workernum)[1]
}

#' RISC \code{ncore} during sweep: avoid oversubscription when sweep uses multiple workers
#' @keywords internal
risc_ncore_for_sweep <- function(workernum, stability_workernum = NULL) {
  sw <- resolve_stability_workernum(workernum, stability_workernum)
  if (sw > 1L) {
    return(1L)
  }
  as.integer(workernum)[1]
}

#' Cap foreach worker count by number of stability tasks to run.
#' @keywords internal
.stability_effective_workers <- function(n_tasks, stability_workernum) {
  max(1L, min(as.integer(stability_workernum)[1], as.integer(n_tasks)[1]))
}

# ---- Cluster stability: RISC/Seurat backends ----

#' Build RISC object list from mdlist (optionally cell subset)
#' @param risc_ncore Passed to \code{scMultiIntegrate}; per-sample parallel uses min(risc_ncore, n samples).
#' @keywords internal
.risc_stability_build_risclist <- function(
    mdlist,
    sample_metadata,
    matlist = NULL,
    tmpobjdir = NULL,
    datadir = NULL,
    matlist_path = NULL,
    cell_rownames = NULL,
    input_seurat_obj = FALSE,
    risc_ncore = 1L
) {
  risc_ncore <- as.integer(risc_ncore)[1]
  if (is.null(matlist)) {
    matlist <- .stability_load_matlist(
      matlist_path = matlist_path,
      tmpobjdir = tmpobjdir,
      mdlist = mdlist,
      sample_metadata = sample_metadata,
      datadir = datadir,
      input_seurat_obj = input_seurat_obj
    )
  }
  bigmat <- do.call(cbind, matlist)
  num_nonzeros <- tabulate(bigmat@i + 1L, nbins = nrow(bigmat))
  joint_filt_genes <- rownames(bigmat)[num_nonzeros >= 3L]

  risclist <- lapply(sample_metadata$Code, function(code) {
    md <- mdlist[[code]]
    if (!is.null(cell_rownames)) {
      if ("Barcode" %in% colnames(md)) {
        md <- md[as.character(md$Barcode) %in% cell_rownames, , drop = FALSE]
      } else {
        md <- md[rownames(md) %in% cell_rownames, , drop = FALSE]
      }
    }
    mat0 <- matlist[[code]]
    mat0 <- mat0[rownames(mat0) %in% joint_filt_genes, , drop = FALSE]
    mat0 <- mat0[match(joint_filt_genes, rownames(mat0)), , drop = FALSE]
    if (nrow(md) == 0L || ncol(mat0) < 1L) {
      return(NULL)
    }
    cell_order <- rownames(md)
    cell_order <- cell_order[cell_order %in% colnames(mat0)]
    if (length(cell_order) < 2L) {
      return(NULL)
    }
    md <- md[cell_order, , drop = FALSE]
    mat0 <- mat0[, cell_order, drop = FALSE]
    coldata0 <- md
    barcodes <- stringr::str_split_fixed(rownames(coldata0), "-", 2)[, 1]
    barcodes <- paste0(coldata0$orig.ident, ".", barcodes)
    coldata0 <- cbind(barcodes, coldata0)
    rowdata0 <- data.frame(Symbol = rownames(mat0), row.names = rownames(mat0))
    RISC::readsc(mat0, coldata0, rowdata0, is.filter = FALSE)
  })
  keep <- !vapply(risclist, is.null, logical(1))
  if (sum(keep) < 2L) {
    stop("RISC stability requires at least two samples with cells.", call. = FALSE)
  }
  risclist <- risclist[keep]
  sample_metadata <- sample_metadata[keep, , drop = FALSE]

  process0 <- function(obj0) {
    obj0 <- RISC::scFilter(obj0, min.UMI = 0, max.UMI = Inf, min.gene = 0, min.cell = 0, is.filter = FALSE)
    obj0 <- RISC::scNormalize(obj0, ncore = 1)
    obj0 <- RISC::scDisperse(obj0)
    obj0
  }

  if (risc_ncore > 1L && length(risclist) > 1L) {
    cl <- parallel::makeCluster(min(risc_ncore, length(risclist)))
    doParallel::registerDoParallel(cl)
    risclist <- foreach::foreach(
      dat0 = risclist,
      .packages = c("RISC")
    ) %dopar% process0(dat0)
    parallel::stopCluster(cl)
  } else {
    risclist <- lapply(risclist, process0)
  }
  names(risclist) <- sample_metadata$Code
  risclist
}

#' Look up a frozen RISC reference in a (possibly subset) risclist.
#'
#' Does not auto-select. Returns \code{NA} if the sample is absent (e.g. bootstrap dropout).
#' @keywords internal
.risc_stability_select_reference <- function(
    risclist,
    sample_metadata,
    risc_reference = NULL
) {
  ref <- .risc_reference_index_in_list(risclist, sample_metadata, risc_reference)
  list(ref = ref, refscore = NULL)
}

#' Integrate RISC objects at given PC count
#' @keywords internal
.risc_stability_integrate <- function(
    risclist,
    sample_metadata,
    pcs_int,
    risc_reference = NULL,
    risc_ncore = 1L
) {
  var0 <- Reduce(intersect, lapply(risclist, function(x) x@rowdata$Symbol))
  ref_info <- .risc_stability_select_reference(risclist, sample_metadata, risc_reference)
  ref <- ref_info$ref
  if (is.na(ref) || is.null(ref)) {
    stop(
      "RISC reference '", risc_reference,
      "' is not present in this RISC object list.",
      call. = FALSE
    )
  }
  data0 <- .risc_put_reference_first(risclist, ref)
  data0 <- RISC::scMultiIntegrate(
    objects = data0,
    eigens = as.integer(pcs_int)[1],
    add.Id = NULL,
    var.gene = var0,
    align = "OLS",
    npc = 50,
    adjust = TRUE,
    ncore = as.integer(risc_ncore)[1]
  )
  list(risc = data0, refscore = ref_info$refscore, var0 = var0)
}

#' Global reference: integrate+cluster for one PC/res (Seurat; merged counts)
#' @keywords internal
.stability_global_ref_worker_seurat <- function(
    params_i,
    path_sobj_merged_temp,
    sweep_config,
    globrefdir,
    verbose = FALSE
) {
  scDAPP::set_parallel_blas_threads()
  outdf <- NULL
  try({
    sobj <- readRDS(path_sobj_merged_temp)
    parsed <- scDAPP:::.parse_stability_param_id(params_i)
    outdf <- scDAPP:::.seurat_stability_integrate_cluster_at(
      sobj,
      sweep_config,
      parsed$pcs_int,
      parsed$res_int
    )
    saveRDS(outdf, file.path(globrefdir, paste0(params_i, ".rds")))
  }, silent = !verbose)
  outdf
}

#' Global reference: integrate+cluster for one PC/res (RISC; worker loads risclist)
#' @keywords internal
.stability_global_ref_worker_risc <- function(
    params_i,
    path_risclist_temp,
    sample_metadata,
    risc_reference,
    risc_ncore,
    risc_louvain_neighbors,
    globrefdir,
    verbose = FALSE
) {
  scDAPP::set_parallel_blas_threads()
  outdf <- NULL
  try({
    risclist <- readRDS(path_risclist_temp)
    parsed <- scDAPP:::.parse_stability_param_id(params_i)
    outdf <- scDAPP:::.risc_stability_integrate_cluster_at(
      risclist,
      sample_metadata,
      parsed$pcs_int,
      parsed$res_int,
      risc_reference = risc_reference,
      risc_ncore = risc_ncore,
      RISC_louvain_neighbors = risc_louvain_neighbors
    )
    saveRDS(outdf, file.path(globrefdir, paste0(params_i, ".rds")))
  }, silent = !verbose)
  outdf
}

#' Bootstrap replicate worker (Seurat): integrate per PC, cluster all res at that PC
#' @keywords internal
.stability_bootstrap_worker_seurat <- function(
    repi,
    paramfield,
    path_sobj_merged_temp,
    sweep_config,
    samplecolumn,
    propcells.perrep,
    bootstrapresdir,
    downsize_seurat,
    verbose = FALSE
) {
  scDAPP::set_parallel_blas_threads()
  set.seed(repi)
  bs_id <- paste0("Bootstrap-", repi)
  try({
    sobj <- readRDS(path_sobj_merged_temp)
    if (downsize_seurat) {
      sobj <- scDAPP:::.stability_downsize_seurat_counts(sobj)
    }
    subcells <- scDAPP:::.stability_subsample_cells(
      sobj@meta.data,
      samplecolumn,
      propcells.perrep
    )
    sobj <- sobj[, subcells]
    backend <- list(
      integrate_at_pcs = function(pcs) {
        scDAPP:::.seurat_integrate_layers_at_pcs(sobj, sweep_config, pcs)
      },
      cluster_at_integrated = function(obj, pcs, res) {
        scDAPP:::.seurat_stability_cluster_at(obj, sweep_config, pcs, res)
      }
    )
    scDAPP:::.stability_sweep_save_params(
      backend,
      paramfield,
      bootstrapresdir,
      filename_prefix = paste0(bs_id, "."),
      verbose = verbose
    )
    invisible(gc(full = TRUE, reset = FALSE, verbose = FALSE))
  }, silent = !verbose)
  invisible(NULL)
}

#' Bootstrap replicate worker (RISC): build risclist once, integrate per PC, cluster all res
#' @keywords internal
.stability_bootstrap_worker_risc <- function(
    repi,
    paramfield,
    path_globmd_temp,
    mdlist,
    sample_metadata,
    path_matlist_temp,
    risc_reference,
    risc_ncore,
    RISC_louvain_neighbors,
    samplecolumn,
    propcells.perrep,
    bootstrapresdir,
    verbose = FALSE
) {
  scDAPP::set_parallel_blas_threads()
  set.seed(repi)
  bs_id <- paste0("Bootstrap-", repi)
  try({
    globmd <- readRDS(path_globmd_temp)
    subcells <- scDAPP:::.stability_subsample_cells(
      globmd,
      samplecolumn,
      propcells.perrep
    )
    matlist <- readRDS(path_matlist_temp)
    risclist <- scDAPP:::.risc_stability_build_risclist(
      mdlist,
      sample_metadata,
      matlist = matlist,
      cell_rownames = subcells,
      risc_ncore = risc_ncore
    )
    ref_idx <- scDAPP:::.risc_reference_index_in_list(
      risclist,
      sample_metadata,
      risc_reference
    )
    if (is.na(ref_idx)) {
      warning(
        bs_id, ": RISC reference '", risc_reference,
        "' missing from subsample; skipping replicate.",
        call. = FALSE
      )
    } else {
      backend <- list(
        integrate_at_pcs = function(pcs) {
          scDAPP:::.risc_stability_integrate(
            risclist,
            sample_metadata,
            pcs,
            risc_reference = risc_reference,
            risc_ncore = risc_ncore
          )
        },
        cluster_at_integrated = function(risc_state, pcs, res) {
          scDAPP:::.risc_stability_cluster_at(
            risc_state,
            pcs,
            res,
            RISC_louvain_neighbors = RISC_louvain_neighbors
          )
        }
      )
      scDAPP:::.stability_sweep_save_params(
        backend,
        paramfield,
        bootstrapresdir,
        filename_prefix = paste0(bs_id, "."),
        verbose = verbose
      )
      invisible(gc(full = TRUE, reset = FALSE, verbose = FALSE))
    }
  }, silent = !verbose)
  invisible(NULL)
}

#' Louvain cluster RISC integrated object at one stability grid point.
#' @keywords internal
.risc_stability_cluster_at <- function(risc_state, pcs_int, res, RISC_louvain_neighbors = 10L) {
  data0 <- risc_state$risc
  data0 <- RISC::scCluster(
    data0,
    slot = "cell.pls",
    method = "louvain",
    npc = as.integer(pcs_int)[1],
    res = as.numeric(res)[1],
    neighbor = RISC_louvain_neighbors
  )
  md <- data0@coldata
  clusters <- as.character(md$Cluster)
  if ("scBarcode" %in% colnames(md)) {
    barcodes <- as.character(md$scBarcode)
  } else if ("Barcode" %in% colnames(md)) {
    barcodes <- as.character(md$Barcode)
  } else {
    barcodes <- rownames(md)
  }
  data.frame(
    barcode = barcodes,
    cluster = clusters,
    stringsAsFactors = FALSE
  )
}

#' Integrate RISC at pcs_int then cluster at res_int for stability reference run.
#' @keywords internal
.risc_stability_integrate_cluster_at <- function(
    risclist,
    sample_metadata,
    pcs_int,
    res_int,
    risc_reference = NULL,
    risc_ncore = 1L,
    RISC_louvain_neighbors = 10L
) {
  risc_state <- .risc_stability_integrate(
    risclist,
    sample_metadata,
    pcs_int,
    risc_reference = risc_reference,
    risc_ncore = risc_ncore
  )
  .risc_stability_cluster_at(
    risc_state,
    pcs_int,
    res_int,
    RISC_louvain_neighbors = RISC_louvain_neighbors
  )
}

#' RISC stability backend closures
#' @keywords internal
.stability_risc_backend <- function(
    config,
    mdlist,
    sample_metadata,
    tmpobjdir = NULL,
    matlist_path = NULL,
    datadir = NULL,
    risc_reference = NULL,
    risc_ncore = 1L,
    input_seurat_obj = FALSE,
    RISC_louvain_neighbors = 10L
) {
  state <- new.env(parent = emptyenv())
  state$mdlist <- mdlist
  state$sample_metadata <- sample_metadata
  state$tmpobjdir <- tmpobjdir
  state$matlist_path <- matlist_path
  state$datadir <- datadir
  state$input_seurat_obj <- input_seurat_obj
  state$matlist <- NULL
  state$risc_ncore <- as.integer(risc_ncore)[1]
  state$risc_reference <- risc_reference
  state$RISC_louvain_neighbors <- RISC_louvain_neighbors
  state$risclist <- NULL
  state$globmd <- NULL

  list(
    engine = "RISC",
    samplecolumn = "Code",
    config = config,
    mdlist = mdlist,
    tmpobjdir = tmpobjdir,
    matlist_path = matlist_path,
    datadir = datadir,
    sample_metadata = sample_metadata,
    risc_reference = risc_reference,
    risc_ncore = as.integer(risc_ncore)[1],
    input_seurat_obj = input_seurat_obj,
    RISC_louvain_neighbors = RISC_louvain_neighbors,
    get_globmd = function() state$globmd,
    get_risclist = function() state$risclist,
    get_matlist = function() state$matlist,
    set_risc_reference = function(code) {
      state$risc_reference <<- code
    },
    prep_reference = function() {
      state$matlist <<- .stability_load_matlist(
        matlist_path = state$matlist_path,
        tmpobjdir = state$tmpobjdir,
        mdlist = state$mdlist,
        sample_metadata = state$sample_metadata,
        datadir = state$datadir,
        input_seurat_obj = state$input_seurat_obj
      )
      state$risclist <<- .risc_stability_build_risclist(
        state$mdlist,
        state$sample_metadata,
        matlist = state$matlist,
        cell_rownames = NULL,
        risc_ncore = state$risc_ncore
      )
      state$globmd <<- do.call(rbind, state$mdlist)
      invisible(state$risclist)
    },
    integrate_at_pcs = function(pcs) {
      if (is.null(state$risclist)) {
        stop("Call prep_reference() before integrate_at_pcs() for RISC.", call. = FALSE)
      }
      .risc_stability_integrate(
        state$risclist,
        state$sample_metadata,
        pcs,
        risc_reference = state$risc_reference,
        risc_ncore = state$risc_ncore
      )
    },
    cluster_at_integrated = function(risc_state, pcs, res) {
      .risc_stability_cluster_at(
        risc_state,
        pcs,
        res,
        RISC_louvain_neighbors = state$RISC_louvain_neighbors
      )
    },
    integrate_cluster_at = function(pcs, res) {
      if (is.null(state$risclist)) {
        stop("Call prep_reference() before integrate_cluster_at() for RISC.", call. = FALSE)
      }
      .risc_stability_integrate_cluster_at(
        state$risclist,
        state$sample_metadata,
        pcs,
        res,
        risc_reference = state$risc_reference,
        risc_ncore = state$risc_ncore,
        RISC_louvain_neighbors = state$RISC_louvain_neighbors
      )
    }
  )
}

#' Seurat object backend from merged tmp objects
#' @keywords internal
.stability_backend_for_config <- function(
    config,
    sample_metadata,
    tmpobjdir = NULL,
    mdlist = NULL,
    datadir = NULL,
    matlist_path = NULL,
    risc_reference = NULL,
    risc_ncore = 1L,
    input_seurat_obj = FALSE,
    RISC_louvain_neighbors = 10L
) {
  if (config$engine == "RISC") {
    if (is.null(mdlist)) {
      stop("mdlist is required for RISC stability sweep.", call. = FALSE)
    }
    if (isTRUE(input_seurat_obj)) {
      if (is.null(matlist_path) && is.null(tmpobjdir)) {
        stop(
          "RISC stability requires matlist_path or tmpobjdir when input_seurat_obj is TRUE.",
          call. = FALSE
        )
      }
    } else if (is.null(datadir)) {
      stop("datadir is required for RISC stability when input_seurat_obj is FALSE.", call. = FALSE)
    }
    return(.stability_risc_backend(
      config,
      mdlist,
      sample_metadata,
      tmpobjdir = tmpobjdir,
      matlist_path = matlist_path,
      datadir = datadir,
      risc_reference = risc_reference,
      risc_ncore = risc_ncore,
      input_seurat_obj = input_seurat_obj,
      RISC_louvain_neighbors = RISC_louvain_neighbors
    ))
  }
  if (is.null(tmpobjdir)) {
    stop("tmpobjdir is required for Seurat stability sweep.", call. = FALSE)
  }
  sweep_config <- config
  state <- new.env(parent = emptyenv())
  state$sobj <- .seurat_merge_from_tmpobjdir(tmpobjdir, sample_metadata)
  state$sweep_config <- sweep_config
  list(
    engine = "Seurat",
    samplecolumn = "Code",
    sweep_config = sweep_config,
    sobj = state$sobj,
    sobj_full = state$sobj,
    set_sobj = function(x) {
      state$sobj <- x
    },
    get_sobj = function() state$sobj,
    prep_reference = function() {
      invisible(NULL)
    },
    integrate_at_pcs = function(pcs) {
      .seurat_integrate_layers_at_pcs(state$sobj, state$sweep_config, pcs)
    },
    cluster_at_integrated = function(obj, pcs, res) {
      .seurat_stability_cluster_at(obj, state$sweep_config, pcs, res)
    },
    integrate_cluster_at = function(pcs, res) {
      .seurat_stability_integrate_cluster_at(state$sobj, state$sweep_config, pcs, res)
    }
  )
}

# ---- Cluster stability: sweep API ----

#' Compute stability metrics (ARI / Jaccard) for one parameter combo
#' @keywords internal
.stability_metrics_one_param <- function(
    params_i,
    globref,
    bs_filenames,
    repres_flat,
    verbose = FALSE
) {
  bsindex_params_i <- bs_filenames
  thisparam_bsres_l <- lapply(bsindex_params_i, function(bsidx_i) {
    if (verbose) message(" - ", bsidx_i)
    bsidx_i_df <- repres_flat[[bsidx_i]]
    globref_bsidx_i <- globref
    intbarcodes <- intersect(globref_bsidx_i$barcode, bsidx_i_df$barcode)
    bsidx_i_df <- bsidx_i_df[match(intbarcodes, bsidx_i_df$barcode), , drop = FALSE]
    globref_bsidx_i <- globref_bsidx_i[match(intbarcodes, globref_bsidx_i$barcode), , drop = FALSE]

    ari_bsidx_i <- mclust::adjustedRandIndex(
      as.vector(globref_bsidx_i$cluster),
      as.vector(bsidx_i_df$cluster)
    )

    fullclusts <- stringr::str_sort(unique(globref$cluster), numeric = TRUE)
    bsclusts <- stringr::str_sort(unique(bsidx_i_df$cluster), numeric = TRUE)

    jaccard_l <- lapply(fullclusts, function(fci) {
      jaccardscore <- NA_real_
      try({
        globref_bsidx_i_fci <- globref_bsidx_i[globref_bsidx_i$cluster == fci, , drop = FALSE]
        bstabmatch_l <- lapply(bsclusts, function(bci) {
          bsidx_i_df_bci <- bsidx_i_df[bsidx_i_df$cluster == bci, , drop = FALSE]
          tabmatch <- table(factor(
            bsidx_i_df_bci$barcode %in% globref_bsidx_i_fci$barcode,
            levels = c("FALSE", "TRUE")
          ))
          data.frame(
            BootstrapClust = bci,
            Match = tabmatch["TRUE"],
            NoMatch = tabmatch["FALSE"],
            row.names = paste0("BSclust.", bci)
          )
        })
        bstabmatch <- dplyr::bind_rows(bstabmatch_l)
        bsmatch <- bstabmatch[which.max(bstabmatch$Match), , drop = FALSE]
        fci_barcodes <- globref_bsidx_i_fci$barcode
        bci_barcodes <- bsidx_i_df[bsidx_i_df$cluster == bsmatch$BootstrapClust, "barcode", drop = TRUE]
        jaccardscore <- length(intersect(fci_barcodes, bci_barcodes)) /
          length(union(fci_barcodes, bci_barcodes))
      }, silent = TRUE)
      jaccardscore
    })

    jaccard_vec <- unlist(jaccard_l, use.names = FALSE)
    jaccard_vec_nonna <- stats::na.omit(jaccard_vec)
    mean_jaccard_score <- if (length(jaccard_vec_nonna) == 0) {
      0
    } else {
      mean(jaccard_vec_nonna)
    }

    perclust_jaccard_df <- data.frame(
      Jaccard = jaccard_vec,
      row.names = fullclusts,
      stringsAsFactors = FALSE
    )
    colnames(perclust_jaccard_df) <- bsidx_i

    list(
      ari_meanjaccard_df = data.frame(
        params_i = params_i,
        bsidx_i = bsidx_i,
        ARI = ari_bsidx_i,
        MeanJaccard_PerClust = mean_jaccard_score,
        nClust_Ref = length(fullclusts),
        nClust_Bootstrap = length(bsclusts),
        stringsAsFactors = FALSE
      ),
      perclust_jaccard_df = perclust_jaccard_df
    )
  })

  ari_meanjaccard_df_l <- dplyr::bind_rows(lapply(thisparam_bsres_l, `[[`, "ari_meanjaccard_df"))
  rownames(ari_meanjaccard_df_l) <- NULL

  parammeandf <- data.frame(
    params_i = params_i,
    nClust_Ref = ari_meanjaccard_df_l$nClust_Ref[1],
    nClust_Bootstrap_Mean = mean(ari_meanjaccard_df_l$nClust_Bootstrap),
    nClust_Bootstrap_SD = stats::sd(ari_meanjaccard_df_l$nClust_Bootstrap),
    ARI_mean = mean(ari_meanjaccard_df_l$ARI),
    ARI_sd = stats::sd(ari_meanjaccard_df_l$ARI),
    Jaccard_mean_of_clustermeans = mean(ari_meanjaccard_df_l$MeanJaccard_PerClust),
    Jaccard_sd_of_clustermeans = stats::sd(ari_meanjaccard_df_l$MeanJaccard_PerClust),
    stringsAsFactors = FALSE
  )

  perclust_jaccard_df <- dplyr::bind_cols(lapply(thisparam_bsres_l, `[[`, "perclust_jaccard_df"))
  perclust_jaccard_long <- perclust_jaccard_df %>%
    tibble::rownames_to_column("cluster") %>%
    tidyr::pivot_longer(cols = dplyr::starts_with("Bootstrap"))
  perclust_jaccard_long$Bootstrap <- stringr::str_split_fixed(perclust_jaccard_long$name, "\\.", 2)[, 1]
  perclust_jaccard_long$Param <- stringr::str_split_fixed(perclust_jaccard_long$name, "\\.", 2)[, 2]

  perclust_jaccard_mean <- data.frame(
    cparams_i = params_i,
    cluster = rownames(perclust_jaccard_df),
    perclust_jaccard_mean = rowMeans(perclust_jaccard_df),
    perclust_jaccard_sd = apply(perclust_jaccard_df, 1, stats::sd),
    stringsAsFactors = FALSE
  )
  rownames(perclust_jaccard_mean) <- NULL

  list(
    ari_meanjaccard_df_l = ari_meanjaccard_df_l,
    parammeandf = parammeandf,
    perclust_jaccard_long = perclust_jaccard_long,
    perclust_jaccard_mean = perclust_jaccard_mean
  )
}

#' Aggregate stability metrics across the parameter field
#' @keywords internal
.compute_stability_metrics <- function(
    globref_l,
    repres_flat,
    bsindex,
    paramfield,
    verbose = FALSE,
    stability_workernum = 1L
) {
  metrics_worker <- function(params_i) {
    if (verbose) message(params_i)
    bs_fns <- bsindex$Filename[bsindex$ParamID == params_i]
    .stability_metrics_one_param(
      params_i,
      globref_l[[params_i]],
      bs_fns,
      repres_flat,
      verbose = verbose
    )
  }

  n_workers <- .stability_effective_workers(length(paramfield), stability_workernum)
  if (n_workers <= 1L) {
    param_outlist_l <- lapply(paramfield, metrics_worker)
  } else {
    cl <- parallel::makeCluster(n_workers, outfile = "")
    doParallel::registerDoParallel(cl)
    .stability_cluster_register_workers(cl)
    parallel::clusterExport(
      cl,
      c("globref_l", "repres_flat", "bsindex", "paramfield", "verbose"),
      envir = environment()
    )
    param_outlist_l <- foreach::foreach(
      params_i = paramfield,
      .packages = .stability_foreach_packages(
        c("mclust", "dplyr", "stringr", "tibble", "tidyr")
      ),
      .verbose = verbose
    ) %dopar% {
      scDAPP::set_parallel_blas_threads()
      bs_fns <- bsindex$Filename[bsindex$ParamID == params_i]
      scDAPP:::.stability_metrics_one_param(
        params_i,
        globref_l[[params_i]],
        bs_fns,
        repres_flat,
        verbose = verbose
      )
    }
    names(param_outlist_l) <- paramfield
    parallel::stopCluster(cl)
  }

  ari_meanjaccard_df <- dplyr::bind_rows(lapply(param_outlist_l, `[[`, "ari_meanjaccard_df_l"))
  parammeandf <- dplyr::bind_rows(lapply(param_outlist_l, `[[`, "parammeandf"))
  perclust_jaccard_long <- dplyr::bind_rows(lapply(param_outlist_l, `[[`, "perclust_jaccard_long"))
  perclust_jaccard_mean <- dplyr::bind_rows(lapply(param_outlist_l, `[[`, "perclust_jaccard_mean"))
  rownames(perclust_jaccard_mean) <- NULL

  parammeandf$combinedscore <-
    0.5 * parammeandf$ARI_mean + 0.5 * parammeandf$Jaccard_mean_of_clustermeans
  parammeandf$MAX <- ""
  parammeandf[which.max(parammeandf$combinedscore), "MAX"] <- "*"

  list(
    perparam_meanscores = parammeandf,
    perbootstrap_perparam_scores = ari_meanjaccard_df,
    perclust_jaccard_mean_acrossbootstraps = perclust_jaccard_mean,
    perclust_jaccard_per_bootstrap_long = perclust_jaccard_long
  )
}

#' Run global + bootstrap clustering over a PC/res grid and score stability
#'
#' Uses resampling stability (ARI and cluster-matched Jaccard) to compare
#' bootstrap clusterings to full-data reference clusterings.
#'
#' @param backend List from \code{.stability_backend_for_config()} with
#'   \code{prep_reference}, \code{integrate_at_pcs}, \code{cluster_at_integrated}, and \code{samplecolumn}.
#' @param paramfield Character vector of parameter ids (see \code{.build_stability_paramfield()}).
#' @param numreps Number of bootstrap replicates.
#' @param propcells.perrep Fraction of cells per bootstrap.
#' @param stability_workernum Parallel workers for sweep stages (global ref, bootstrap, metrics).
#' @param outdir Output directory (checkpoint/resume supported).
#' @param save_plots If \code{TRUE}, run \code{integration_stability_plots_module()}
#'   after metrics are computed (PDFs under \code{outdir/plots/}).
#' @param label_top_k Top-scoring combinations to label on the ARI vs Jaccard scatter.
#' @param remove_rawouts If \code{TRUE}, delete \code{outdir/RawOuts/} after a
#'   successful run (CSVs and plots retained). Failed runs keep \code{RawOuts} for resume.
#' @param verbose Logical.
#' @param downsize_seurat For Seurat backends, strip to counts assay before sweep.
#' @return List with stability tables (see \code{perparam_meanscores}) and optional
#'   \code{plots} from \code{integration_stability_plots_module()}.
#' @export
cluster_stability_sweep <- function(
    backend,
    paramfield,
    numreps = 50L,
    propcells.perrep = 0.8,
    stability_workernum = 1L,
    outdir = "./cluster_stability",
    save_plots = TRUE,
    label_top_k = 10L,
    remove_rawouts = TRUE,
    verbose = TRUE,
    downsize_seurat = TRUE
) {
  .check_stability_dependencies()
  set_parallel_blas_threads()
  set.seed(54321)

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  logdir <- file.path(outdir, "logs")
  dir.create(logdir, recursive = TRUE, showWarnings = FALSE)
  rawoutsdir <- file.path(outdir, "RawOuts")
  dir.create(rawoutsdir, recursive = TRUE, showWarnings = FALSE)
  globrefdir <- file.path(rawoutsdir, "GlobalReferenceClustering")
  dir.create(globrefdir, recursive = TRUE, showWarnings = FALSE)
  bootstrapresdir <- file.path(rawoutsdir, "PerBootstrapClustering")
  dir.create(bootstrapresdir, recursive = TRUE, showWarnings = FALSE)

  samplecolumn <- backend$samplecolumn
  is_risc <- identical(backend$engine, "RISC")
  stability_workernum <- as.integer(stability_workernum)[1]
  risc_reference_selection <- NULL

  if (is_risc) {
    existing_sel <- .risc_reference_selection_read(outdir)
    spec_now <- .parse_risc_reference_arg(backend$risc_reference)
    reuse_sel <- !is.null(existing_sel) &&
      !is.null(existing_sel$selected_code) &&
      identical(existing_sel$user_arg, spec_now$raw) &&
      existing_sel$selected_code %in% backend$sample_metadata$Code
    if (is.null(backend$get_risclist())) {
      if (verbose) message("Preparing RISC objects for stability sweep")
      backend$prep_reference()
    }
    if (isTRUE(reuse_sel)) {
      risc_reference_selection <- existing_sel
      if (is.function(backend$set_risc_reference)) {
        backend$set_risc_reference(existing_sel$selected_code)
      }
      backend$risc_reference <- existing_sel$selected_code
      if (verbose) {
        message(
          "RISC reference (", existing_sel$mode, "): ",
          existing_sel$selected_code, " (from previous selection)"
        )
      }
    } else {
      risc_reference_selection <- .stability_risc_freeze_reference(
        backend, outdir, verbose
      )
      backend$risc_reference <- risc_reference_selection$selected_code
    }
  }

  # --- Global reference clustering ---
  saved_globrefs <- gsub("\\.rds$", "", list.files(globrefdir))
  globparamfield_run <- paramfield[!(paramfield %in% saved_globrefs)]

  if (length(globparamfield_run) > 0L) {
    if (verbose) message("Computing reference clusters across parameter field")
    if (!is_risc && downsize_seurat && !is.null(backend$get_sobj())) {
      if (verbose) message("Downsizing Seurat object to counts")
      backend$set_sobj(.stability_downsize_seurat_counts(backend$get_sobj()))
    }
    if (!is_risc) {
      if (verbose) message("Preparing reference data (integrate per PC during sweep)")
      backend$prep_reference()
    }

    ref_workers <- .stability_effective_workers(
      length(globparamfield_run),
      stability_workernum
    )
    if (ref_workers <= 1L) {
      .stability_sweep_save_params(
        backend,
        globparamfield_run,
        globrefdir,
        verbose = verbose
      )
    } else if (is_risc) {
      path_risclist_temp <- file.path(rawoutsdir, "TEMP_RISCLIST.RDS")
      saveRDS(backend$get_risclist(), path_risclist_temp)
      risc_louvain_neighbors <- backend$RISC_louvain_neighbors
      risc_ncore <- backend$risc_ncore
      risc_reference <- backend$risc_reference
      sample_metadata <- backend$sample_metadata
      cl <- parallel::makeCluster(ref_workers, outfile = "")
      doParallel::registerDoParallel(cl)
      .stability_cluster_register_workers(cl)
      parallel::clusterExport(
        cl,
        c(
          "path_risclist_temp",
          "sample_metadata",
          "risc_reference",
          "risc_ncore",
          "risc_louvain_neighbors",
          "globrefdir",
          "verbose"
        ),
        envir = environment()
      )
      foreach::foreach(
        params_i = globparamfield_run,
        .packages = .stability_foreach_packages(),
        .verbose = verbose
      ) %dopar% {
        scDAPP:::.stability_global_ref_worker_risc(
          params_i,
          path_risclist_temp,
          sample_metadata,
          risc_reference,
          risc_ncore,
          risc_louvain_neighbors,
          globrefdir,
          verbose
        )
      }
      parallel::stopCluster(cl)
    } else {
      path_sobj_merged_temp <- file.path(rawoutsdir, "TEMP_SOBJ_MERGED.RDS")
      saveRDS(backend$get_sobj(), path_sobj_merged_temp)
      sweep_config <- backend$sweep_config
      cl <- parallel::makeCluster(ref_workers, outfile = "")
      doParallel::registerDoParallel(cl)
      .stability_cluster_register_workers(cl)
      parallel::clusterEvalQ(cl, Sys.sleep(stats::runif(1, 1, 3)))
      parallel::clusterExport(
        cl,
        c("path_sobj_merged_temp", "sweep_config", "globrefdir", "verbose"),
        envir = environment()
      )
      foreach::foreach(
        params_i = globparamfield_run,
        .packages = .stability_foreach_packages(),
        .verbose = verbose
      ) %dopar% {
        scDAPP:::.stability_global_ref_worker_seurat(
          params_i,
          path_sobj_merged_temp,
          sweep_config,
          globrefdir,
          verbose
        )
      }
      parallel::stopCluster(cl)
    }
    invisible(gc(full = TRUE, reset = FALSE, verbose = FALSE))
  }

  globref_l <- lapply(paramfield, function(params_i) {
    f <- file.path(globrefdir, paste0(params_i, ".rds"))
    if (file.exists(f)) readRDS(f) else NULL
  })
  names(globref_l) <- paramfield
  globref_l <- globref_l[lengths(globref_l) > 0L]
  paramfield <- paramfield[paramfield %in% names(globref_l)]
  if (!length(paramfield)) {
    stop("All global reference parameter combinations failed.", call. = FALSE)
  }

  # --- Bootstrap clustering ---
  bootstrapreps_run_full_df <- do.call(
    rbind,
    lapply(seq_len(numreps), function(repi) {
      data.frame(
        Filename = paste0("Bootstrap-", repi, ".", paramfield),
        BootstrapID = paste0("Bootstrap-", repi),
        ParamID = paramfield,
        stringsAsFactors = FALSE
      )
    })
  )
  saved_bs <- gsub("\\.rds$", "", list.files(bootstrapresdir))
  bootstrapreps_run <- bootstrapreps_run_full_df$Filename
  bootstrapreps_run <- bootstrapreps_run[!(bootstrapreps_run %in% saved_bs)]
  bootstrapreps_run_unique <- unique(
    stringr::str_split_fixed(bootstrapreps_run, "\\.", 2)[, 1]
  )

  if (length(bootstrapreps_run) > 0L) {
    if (verbose) message("Bootstrap clustering across parameter field")

    bootstrap_rep_ids <- vapply(
      bootstrapreps_run_unique,
      function(bs_id) as.integer(sub("^Bootstrap-", "", bs_id)),
      integer(1)
    )

    if (is_risc) {
      if (is.null(backend$get_globmd())) {
        backend$prep_reference()
      }
      path_globmd_temp <- file.path(rawoutsdir, "TEMP_GLOBMD.RDS")
      saveRDS(backend$get_globmd(), path_globmd_temp)
      path_mdlist_temp <- file.path(rawoutsdir, "TEMP_STABILITY_MDLIST.RDS")
      saveRDS(backend$mdlist, path_mdlist_temp)
      if (is.null(backend$get_matlist())) {
        backend$prep_reference()
      }
      path_matlist_temp <- file.path(rawoutsdir, "TEMP_MATLIST.RDS")
      saveRDS(backend$get_matlist(), path_matlist_temp)
      bs_workers <- .stability_effective_workers(
        length(bootstrap_rep_ids),
        stability_workernum
      )
      risc_louvain_neighbors <- backend$RISC_louvain_neighbors
      risc_ncore <- backend$risc_ncore
      risc_reference <- backend$risc_reference
      sample_metadata <- backend$sample_metadata
      if (bs_workers <= 1L) {
        for (repi in bootstrap_rep_ids) {
          .stability_bootstrap_worker_risc(
            repi,
            paramfield,
            path_globmd_temp,
            readRDS(path_mdlist_temp),
            sample_metadata,
            path_matlist_temp,
            risc_reference,
            risc_ncore,
            risc_louvain_neighbors,
            samplecolumn,
            propcells.perrep,
            bootstrapresdir,
            verbose
          )
        }
      } else {
        cl <- parallel::makeCluster(bs_workers, outfile = "")
        doParallel::registerDoParallel(cl)
        .stability_cluster_register_workers(cl)
        mdlist <- readRDS(path_mdlist_temp)
        parallel::clusterExport(
          cl,
          c(
            "paramfield",
            "path_globmd_temp",
            "mdlist",
            "sample_metadata",
            "path_matlist_temp",
            "risc_reference",
            "risc_ncore",
            "risc_louvain_neighbors",
            "samplecolumn",
            "propcells.perrep",
            "bootstrapresdir",
            "verbose"
          ),
          envir = environment()
        )
        foreach::foreach(
          repi = bootstrap_rep_ids,
          .packages = .stability_foreach_packages(),
          .verbose = verbose
        ) %dopar% {
          scDAPP:::.stability_bootstrap_worker_risc(
            repi,
            paramfield,
            path_globmd_temp,
            mdlist,
            sample_metadata,
            path_matlist_temp,
            risc_reference,
            risc_ncore,
            risc_louvain_neighbors,
            samplecolumn,
            propcells.perrep,
            bootstrapresdir,
            verbose
          )
        }
        parallel::stopCluster(cl)
      }
    } else {
      path_sobj_merged_temp <- file.path(rawoutsdir, "TEMP_SOBJ_MERGED.RDS")
      if (downsize_seurat) {
        saveRDS(
          .stability_downsize_seurat_counts(backend$get_sobj()),
          path_sobj_merged_temp
        )
      } else {
        saveRDS(backend$get_sobj(), path_sobj_merged_temp)
      }
      sweep_config <- backend$sweep_config
      bs_workers <- .stability_effective_workers(
        length(bootstrap_rep_ids),
        stability_workernum
      )
      if (bs_workers <= 1L) {
        for (repi in bootstrap_rep_ids) {
          .stability_bootstrap_worker_seurat(
            repi,
            paramfield,
            path_sobj_merged_temp,
            sweep_config,
            samplecolumn,
            propcells.perrep,
            bootstrapresdir,
            downsize_seurat = FALSE,
            verbose
          )
        }
      } else {
        cl <- parallel::makeCluster(bs_workers, outfile = "")
        doParallel::registerDoParallel(cl)
        .stability_cluster_register_workers(cl)
        parallel::clusterExport(
          cl,
          c(
            "paramfield",
            "path_sobj_merged_temp",
            "sweep_config",
            "samplecolumn",
            "propcells.perrep",
            "bootstrapresdir",
            "verbose"
          ),
          envir = environment()
        )
        foreach::foreach(
          repi = bootstrap_rep_ids,
          .packages = .stability_foreach_packages(),
          .verbose = verbose
        ) %dopar% {
          scDAPP:::.stability_bootstrap_worker_seurat(
            repi,
            paramfield,
            path_sobj_merged_temp,
            sweep_config,
            samplecolumn,
            propcells.perrep,
            bootstrapresdir,
            downsize_seurat = FALSE,
            verbose
          )
        }
        parallel::stopCluster(cl)
      }
    }
    invisible(gc(full = TRUE, reset = FALSE, verbose = FALSE))
  }

  repres_flat <- lapply(bootstrapreps_run_full_df$Filename, function(fn) {
    f <- file.path(bootstrapresdir, paste0(fn, ".rds"))
    if (file.exists(f)) readRDS(f) else NULL
  })
  names(repres_flat) <- bootstrapreps_run_full_df$Filename
  repres_flat <- repres_flat[lengths(repres_flat) > 0L]

  bsindex <- bootstrapreps_run_full_df[
    bootstrapreps_run_full_df$Filename %in% names(repres_flat),
    ,
    drop = FALSE
  ]

  if (verbose) message("Computing stability metrics")
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
    write.csv(
      fulloutlist[[nm]],
      file.path(outdir, paste0(nm, ".csv")),
      row.names = FALSE
    )
  }

  if (isTRUE(save_plots)) {
    plot_dir <- file.path(outdir, "plots")
    fulloutlist$plots <- tryCatch(
      integration_stability_plots_module(
        stability = fulloutlist,
        outdir = plot_dir,
        label_top_k = label_top_k,
        save = TRUE
      ),
      error = function(e) {
        warning(
          "Stability plots failed: ",
          conditionMessage(e),
          call. = FALSE
        )
        NULL
      }
    )
  }

  if (isTRUE(remove_rawouts) && dir.exists(rawoutsdir)) {
    unlink(rawoutsdir, recursive = TRUE)
    if (verbose) {
      message("Removed RawOuts checkpoint directory: ", rawoutsdir)
    }
  }

  if (!is.null(risc_reference_selection)) {
    fulloutlist$risc_reference_selection <- risc_reference_selection
  }

  fulloutlist
}

#' Select integration PCs and resolution via cluster stability
#'
#' @param integration_method One of \code{integration_method_choices()}.
#' @param sample_metadata data.frame with Sample, Condition, Code.
#' @param pcs_int \code{"auto"} or fixed integer (used when not auto).
#' @param res_int \code{"auto"} or fixed numeric.
#' @param tmpobjdir Per-sample Seurat RDS directory (Seurat methods).
#' @param mdlist Named metadata lists (RISC).
#' @param datadir Raw data path (RISC).
#' @param outdir Stability output directory.
#' @param numreps Bootstrap replicates (default 50).
#' @param propcells.perrep Cell fraction per bootstrap (default 0.8).
#' @param sweep_maxPCs PC grid when \code{pcs_int = "auto"}.
#' @param sweep_res Resolution grid when \code{res_int = "auto"}.
#' @param workernum Pipeline workers (used for final RISC integration and RISC ncore when sweep is serial).
#' @param stability_workernum Parallel workers for sweep \code{foreach} stages; \code{NULL} uses \code{workernum}.
#' @param risc_reference Optional RISC reference: \code{"autoV2"} (default),
#'   \code{"auto"} for the legacy heuristic, or a sample Code/Sample name.
#' @param RISC_louvain_neighbors RISC Louvain neighbors.
#' @param matlist_path \code{.concatmatrix.rds} from prep (RISC when \code{input_seurat_obj = TRUE}).
#' @param datadir Raw h5 path (RISC when \code{input_seurat_obj = FALSE}).
#' @param input_seurat_obj If \code{TRUE}, use \code{matlist_path}/\code{tmpobjdir} not \code{datadir}.
#' @param verbose Logical.
#' @return List with \code{pcs_int}, \code{res_int}, \code{stability}, \code{sweep_dir},
#'   and for RISC \code{selected_risc_reference} plus \code{risc_reference_selection}.
#' @export
auto_integration_cluster_params <- function(
    integration_method,
    sample_metadata,
    pcs_int = "auto",
    res_int = "auto",
    tmpobjdir = NULL,
    mdlist = NULL,
    matlist_path = NULL,
    datadir = NULL,
    outdir,
    numreps = 50L,
    propcells.perrep = 0.8,
    sweep_maxPCs = c(5, 10, 15, 20, 25, 30, 40, 50),
    sweep_res = seq(0.1, 1.5, by = 0.2),
    workernum = 1L,
    stability_workernum = NULL,
    risc_reference = "autoV2",
    RISC_louvain_neighbors = 10L,
    input_seurat_obj = FALSE,
    verbose = TRUE
) {
  if (!.integration_param_needs_auto(pcs_int, res_int)) {
    stop(
      "auto_integration_cluster_params() requires pcs_int and/or res_int = 'auto'.",
      call. = FALSE
    )
  }
  .check_stability_dependencies()

  method <- match.arg(integration_method, integration_method_choices())
  sweep_config <- resolve_integration_config(method, pcs_int = 30L, res_int = 0.5)
  check_integration_dependencies(sweep_config)

  paramfield <- .build_stability_paramfield(
    sweep_maxPCs,
    sweep_res,
    pcs_int = pcs_int,
    res_int = res_int
  )
  sweep_workers <- resolve_stability_workernum(workernum, stability_workernum)
  risc_ncore <- risc_ncore_for_sweep(workernum, stability_workernum)

  backend <- .stability_backend_for_config(
    sweep_config,
    sample_metadata,
    tmpobjdir = tmpobjdir,
    mdlist = mdlist,
    matlist_path = matlist_path,
    datadir = datadir,
    risc_reference = risc_reference,
    risc_ncore = risc_ncore,
    input_seurat_obj = input_seurat_obj,
    RISC_louvain_neighbors = RISC_louvain_neighbors
  )
  if (sweep_config$engine != "RISC") {
    backend$tmpobjdir <- tmpobjdir
    backend$sample_metadata <- sample_metadata
  }

  stability <- cluster_stability_sweep(
    backend = backend,
    paramfield = paramfield,
    numreps = numreps,
    propcells.perrep = propcells.perrep,
    stability_workernum = sweep_workers,
    outdir = outdir,
    verbose = verbose
  )

  best <- .pick_best_stability_params(stability$perparam_meanscores)
  resolved_pcs <- if (.is_auto_integration_param(pcs_int)) best$pcs_int else as.integer(pcs_int)[1]
  resolved_res <- if (.is_auto_integration_param(res_int)) best$res_int else as.numeric(res_int)[1]

  list(
    pcs_int = resolved_pcs,
    res_int = resolved_res,
    stability = stability,
    sweep_dir = outdir,
    selected_params_i = stability$perparam_meanscores$params_i[
      stability$perparam_meanscores$MAX == "*"
    ][1],
    selected_risc_reference = if (!is.null(stability$risc_reference_selection)) {
      stability$risc_reference_selection$selected_code
    } else {
      NULL
    },
    risc_reference_selection = stability$risc_reference_selection
  )
}

#' Resolve numeric pcs_int / res_int after optional auto selection
#'
#' @param pcs_int Integer or \code{"auto"}.
#' @param res_int Numeric or \code{"auto"}.
#' @param auto_result Optional list from \code{auto_integration_cluster_params()}.
#' @return Named list with numeric \code{pcs_int} and \code{res_int}.
#' @export
resolve_integration_cluster_params <- function(
    pcs_int,
    res_int,
    auto_result = NULL
) {
  if (!.integration_param_needs_auto(pcs_int, res_int)) {
    return(list(
      pcs_int = as.integer(pcs_int)[1],
      res_int = as.numeric(res_int)[1]
    ))
  }
  if (is.null(auto_result)) {
    stop(
      "pcs_int and/or res_int is 'auto' but auto_result was not provided.",
      call. = FALSE
    )
  }
  resolved_pcs <- if (.is_auto_integration_param(pcs_int)) {
    as.integer(auto_result$pcs_int)[1]
  } else {
    as.integer(pcs_int)[1]
  }
  resolved_res <- if (.is_auto_integration_param(res_int)) {
    as.numeric(auto_result$res_int)[1]
  } else {
    as.numeric(res_int)[1]
  }
  list(pcs_int = resolved_pcs, res_int = resolved_res)
}

# ---- Cluster stability: plots ----

#' Load stability sweep tables from a list or output directory
#' @keywords internal
.load_stability_tables <- function(stability = NULL, stability_dir = NULL) {
  if (!is.null(stability)) {
    need <- c(
      "perparam_meanscores",
      "perbootstrap_perparam_scores",
      "perclust_jaccard_mean_acrossbootstraps"
    )
    miss <- setdiff(need, names(stability))
    if (length(miss)) {
      stop(
        "stability must include: ",
        paste(need, collapse = ", "),
        " (missing: ",
        paste(miss, collapse = ", "),
        ")",
        call. = FALSE
      )
    }
    return(stability[need])
  }
  if (is.null(stability_dir) || !nzchar(stability_dir)) {
    stop("Provide stability (list) or stability_dir.", call. = FALSE)
  }
  stability_dir <- normalizePath(stability_dir, mustWork = FALSE)
  if (!dir.exists(stability_dir)) {
    stop("stability_dir does not exist: ", stability_dir, call. = FALSE)
  }
  read_one <- function(nm) {
    f <- file.path(stability_dir, paste0(nm, ".csv"))
    if (!file.exists(f)) {
      stop("Missing stability table: ", f, call. = FALSE)
    }
    utils::read.csv(f, stringsAsFactors = FALSE)
  }
  list(
    perparam_meanscores = read_one("perparam_meanscores"),
    perbootstrap_perparam_scores = read_one("perbootstrap_perparam_scores"),
    perclust_jaccard_mean_acrossbootstraps = read_one("perclust_jaccard_mean_acrossbootstraps")
  )
}

#' Add combinedscore and factor order for stability summary plotting.
#' @keywords internal
.prepare_stability_perparam <- function(perparam) {
  if (!"combinedscore" %in% colnames(perparam)) {
    perparam$combinedscore <-
      0.5 * perparam$ARI_mean + 0.5 * perparam$Jaccard_mean_of_clustermeans
  }
  parsed <- lapply(perparam$params_i, .parse_stability_param_id)
  perparam$pcs_int <- vapply(parsed, function(x) x$pcs_int, integer(1L))
  perparam$res_int <- vapply(parsed, function(x) x$res_int, numeric(1))
  perparam <- perparam[order(-perparam$combinedscore), , drop = FALSE]
  id_levels <- rev(unique(as.character(perparam$params_i)))
  perparam$params_i <- factor(perparam$params_i, levels = id_levels)
  lab_map <- unique(data.frame(
    params_i = as.character(perparam$params_i),
    params_lab = sprintf("PC%d / res %s", perparam$pcs_int, perparam$res_int),
    stringsAsFactors = FALSE
  ))
  lab_levels <- lab_map$params_lab[match(id_levels, lab_map$params_i)]
  perparam$params_lab <- factor(
    lab_map$params_lab[match(as.character(perparam$params_i), lab_map$params_i)],
    levels = lab_levels
  )
  perparam
}

#' Resolve winning params_i from explicit arg or MAX flag or top combined score.
#' @keywords internal
.stability_selected_params_i <- function(perparam, selected_params_i = NULL) {
  if (!is.null(selected_params_i) && nzchar(selected_params_i)) {
    return(selected_params_i)
  }
  if ("MAX" %in% colnames(perparam)) {
    hit <- perparam$params_i[perparam$MAX == "*"]
    if (length(hit)) {
      return(as.character(hit[1]))
    }
  }
  as.character(perparam$params_i[which.max(perparam$combinedscore)])
}

#' Shared ggplot2 theme for stability sweep PDF figures.
#' @keywords internal
.stability_plot_theme <- function() {
  ggplot2::theme_linedraw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

#' Adaptive axis text size for dense params_i labels on the Y axis.
#' @keywords internal
.stability_param_axis_text_size <- function(n_param) {
  n_param <- as.integer(n_param)[1]
  if (is.na(n_param) || n_param <= 16L) {
    return(10)
  }
  if (n_param <= 32L) {
    return(9)
  }
  8
}

#' Save a ggplot to PDF with explicit width and height (non-interactive devices).
#' @keywords internal
.save_stability_ggplot <- function(p, path, width, height) {
  grDevices::pdf(path, width = width, height = height)
  print(p)
  grDevices::dev.off()
  invisible(path)
}

#' PDF/HTML width/height per plot type for a given sweep grid.
#' Ranked Y-axis plots (bar, bootstrap ARI, nclust) cap near 8-10 in so HTML
#' does not downscale labels off the page.
#' @keywords internal
.stability_plot_save_dims <- function(perparam, perclust_df, selected_params_i) {
  perparam <- .prepare_stability_perparam(perparam)
  n_param <- nrow(perparam)
  n_pc <- length(unique(perparam$pcs_int))
  n_res <- length(unique(perparam$res_int))

  # Ranked Y-axis plots: keep height near the HTML column so labels are not
  # downscaled; width stays close to typical report width (~8-10 in).
  bar_w <- min(10, max(8, 0.05 * n_param + 7.5))
  bar_h <- min(10, max(6, 0.12 * n_param + 4))

  perclust_sel <- perclust_df[
    perclust_df$cparams_i == selected_params_i,
    ,
    drop = FALSE
  ]
  n_clust <- max(1L, nrow(perclust_sel))

  list(
    combinedscore_bar = list(width = bar_w, height = bar_h),
    combinedscore_heatmap = list(
      width = min(16, max(7, 1 * n_res + 3)),
      height = min(14, max(5, 0.55 * n_pc + 2))
    ),
    ari_jaccard_scatter = list(width = 10, height = 7),
    bootstrap_ari = list(width = bar_w, height = bar_h),
    nclust_ref = list(width = bar_w, height = bar_h),
    perclust_jaccard = list(
      width = min(20, max(7, 0.35 * n_clust + 2)),
      height = 4
    )
  )
}

#' Bar chart of combined stability scores ranked by params_i.
#' @keywords internal
.plot_stability_combinedscore_bar <- function(perparam, selected_params_i) {
  perparam <- .prepare_stability_perparam(perparam)
  sel <- .stability_selected_params_i(perparam, selected_params_i)
  plot_df <- perparam
  plot_df$is_selected <- as.character(plot_df$params_i) == sel
  y_size <- .stability_param_axis_text_size(nrow(plot_df))

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$combinedscore,
      y = .data$params_lab,
      fill = .data$is_selected
    )
  ) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::scale_fill_manual(
      values = c("FALSE" = "grey70", "TRUE" = "#2166AC"),
      labels = c("FALSE" = "grid point", "TRUE" = "selected"),
      name = NULL
    ) +
    ggplot2::labs(
      title = "Cluster stability by parameter combination",
      subtitle = "Combined score = 0.5 * mean ARI + 0.5 * mean per-cluster Jaccard",
      x = "Combined stability score",
      y = NULL
    ) +
    .stability_plot_theme() +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      axis.text.y = ggplot2::element_text(size = y_size)
    )
}

#' Heatmap of combined stability score over PC x resolution grid.
#' @keywords internal
.plot_stability_combinedscore_heatmap <- function(perparam, selected_params_i) {
  perparam <- .prepare_stability_perparam(perparam)
  sel <- .stability_selected_params_i(perparam, selected_params_i)
  plot_df <- perparam
  plot_df$is_selected <- as.character(plot_df$params_i) == sel
  plot_df$pcs_lab <- factor(plot_df$pcs_int, levels = sort(unique(plot_df$pcs_int)))
  plot_df$res_lab <- factor(plot_df$res_int, levels = sort(unique(plot_df$res_int)))

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$res_lab,
      y = .data$pcs_lab,
      fill = .data$combinedscore
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::geom_point(
      data = plot_df[plot_df$is_selected, , drop = FALSE],
      ggplot2::aes(shape = "selected"),
      color = "black",
      size = 3
    ) +
    ggplot2::scale_shape_manual(values = c("selected" = 4), name = NULL) +
    ggplot2::scale_fill_viridis_c(option = "C", name = "Combined\nscore", limits = c(0, 1)) +
    ggplot2::labs(
      title = "Stability score across PC and resolution grid",
      x = "Louvain resolution",
      y = "Integration PCs (1:N)"
    ) +
    .stability_plot_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      legend.position = "right"
    )
}

#' Scatter of mean ARI vs mean Jaccard with top grid points labeled.
#' @keywords internal
.plot_stability_ari_jaccard_scatter <- function(
    perparam,
    selected_params_i,
    label_top_k = 10L
) {
  perparam <- .prepare_stability_perparam(perparam)
  sel <- .stability_selected_params_i(perparam, selected_params_i)
  plot_df <- perparam
  plot_df$is_selected <- as.character(plot_df$params_i) == sel
  plot_df$label_me <- FALSE
  top_k <- min(as.integer(label_top_k)[1], nrow(plot_df))
  if (top_k > 0L) {
    ord <- order(-plot_df$combinedscore)
    plot_df$label_me[ord[seq_len(top_k)]] <- TRUE
  }

  other_df <- plot_df[!plot_df$is_selected, , drop = FALSE]
  sel_df <- plot_df[plot_df$is_selected, , drop = FALSE]

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$ARI_mean,
      y = .data$Jaccard_mean_of_clustermeans,
      size = .data$combinedscore
    )
  ) +
    ggplot2::geom_point(
      data = other_df,
      color = "grey50",
      alpha = 0.75
    ) +
    ggplot2::geom_point(
      data = sel_df,
      color = "#B2182B",
      alpha = 1,
      stroke = 0.6
    ) +
    ggplot2::scale_size_continuous(range = c(2, 6), name = "Combined\nscore") +
    ggplot2::labs(
      title = "Mean ARI vs mean per-cluster Jaccard",
      subtitle = sprintf(
        "Top %d combinations labeled; selected point drawn on top in red",
        top_k
      ),
      x = "Mean ARI (bootstrap vs reference)",
      y = "Mean per-cluster Jaccard",
      caption = "Red = selected parameter combination; grey = other grid points"
    ) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    .stability_plot_theme() +
    ggplot2::theme(
      legend.position = "right",
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5)
    )

  if (any(plot_df$label_me)) {
    p <- p + ggrepel::geom_text_repel(
      data = plot_df[plot_df$label_me, , drop = FALSE],
      ggplot2::aes(label = .data$params_i),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.35,
      show.legend = FALSE,
      color = "grey20"
    )
  }
  p
}

#' Boxplot of bootstrap ARI distributions per params_i vs reference clustering.
#' @keywords internal
.plot_stability_bootstrap_ari <- function(perbootstrap, perparam) {
  perparam <- .prepare_stability_perparam(perparam)
  plot_df <- perbootstrap
  plot_df$params_lab <- factor(
    perparam$params_lab[match(as.character(plot_df$params_i), as.character(perparam$params_i))],
    levels = levels(perparam$params_lab)
  )
  y_size <- .stability_param_axis_text_size(nlevels(plot_df$params_lab))

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$ARI,
      y = .data$params_lab,
      group = .data$params_lab
    )
  ) +
    ggplot2::geom_boxplot(fill = "grey85", outlier.size = 0.8, width = 0.6) +
    ggplot2::geom_jitter(height = 0.12, alpha = 0.5, size = 0.9) +
    ggplot2::labs(
      title = "Bootstrap ARI by parameter combination",
      x = "ARI (bootstrap vs reference)",
      y = NULL
    ) +
    .stability_plot_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      axis.text.y = ggplot2::element_text(size = y_size)
    )
}

#' Bar chart of reference cluster counts per params_i (full-data integration).
#' @keywords internal
.plot_stability_nclust_ref <- function(perparam, selected_params_i) {
  perparam <- .prepare_stability_perparam(perparam)
  sel <- .stability_selected_params_i(perparam, selected_params_i)
  plot_df <- perparam
  plot_df$is_selected <- as.character(plot_df$params_i) == sel
  y_size <- .stability_param_axis_text_size(nrow(plot_df))

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$nClust_Ref,
      y = .data$params_lab,
      fill = .data$is_selected
    )
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(
      values = c("FALSE" = "grey70", "TRUE" = "#2166AC"),
      guide = "none"
    ) +
    ggplot2::labs(
      title = "Reference cluster count by parameter combination",
      x = "Number of clusters (full data)",
      y = NULL
    ) +
    .stability_plot_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      axis.text.y = ggplot2::element_text(size = y_size)
    )
}

#' Heatmap of per-cluster Jaccard scores for the selected winning params_i.
#' @keywords internal
.plot_stability_perclust_jaccard <- function(perclust_mean, selected_params_i) {
  plot_df <- perclust_mean[
    perclust_mean$cparams_i == selected_params_i,
    ,
    drop = FALSE
  ]
  if (!nrow(plot_df)) {
    stop(
      "No per-cluster Jaccard rows for selected_params_i: ",
      selected_params_i,
      call. = FALSE
    )
  }
  plot_df$cluster <- factor(
    plot_df$cluster,
    levels = plot_df$cluster[order(-plot_df$perclust_jaccard_mean)]
  )
  plot_df$jaccard_lab <- sprintf("%.3f", plot_df$perclust_jaccard_mean)

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$cluster,
      y = 1,
      fill = .data$perclust_jaccard_mean
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$jaccard_lab),
      color = "black",
      size = 3.2
    ) +
    ggplot2::scale_fill_viridis_c(
      option = "C",
      name = "Mean\nJaccard",
      limits = c(0, 1)
    ) +
    ggplot2::labs(
      title = "Per-cluster Jaccard (selected parameter combination)",
      subtitle = selected_params_i,
      x = "Reference cluster",
      y = NULL,
      caption = paste0(
        "Mean best-matched membership overlap vs reference clustering ",
        "across bootstrap replicates (0-1 scale; near 1 = highly reproducible)"
      )
    ) +
    .stability_plot_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )
}

#' Integration cluster-stability visualization module
#'
#' Builds diagnostic plots from \code{cluster_stability_sweep()} tables:
#' ranked combined-score bar chart (params on Y), PC x resolution heatmap,
#' ARI vs Jaccard scatter (selected point drawn last), bootstrap ARI boxplots
#' and reference cluster-count bars (params on Y), and a per-cluster Jaccard
#' heatmap for the selected parameter combination only (tile values annotated).
#'
#' @param stability Optional list returned by \code{cluster_stability_sweep()}
#'   (or a subset with the three main tables).
#' @param stability_dir Directory containing stability CSV exports; used when
#'   \code{stability} is \code{NULL}.
#' @param outdir If \code{save = TRUE}, PDFs are written here (created if needed).
#' @param selected_params_i \code{params_i} id for the winner (e.g.
#'   \code{"PCs_1-10.res_0.3"}). When \code{NULL}, uses the row with
#'   \code{MAX == "*"} or the highest \code{combinedscore}.
#' @param label_top_k Number of top-scoring combinations to label on the ARI vs
#'   Jaccard scatter (default 10).
#' @param save Logical; write PDFs to \code{outdir}. When \code{TRUE}, PDF width and
#'   height are computed from the sweep grid (see \code{?cluster_stability_sweep} docs).
#' @param width,height Ignored when \code{save = TRUE} (reserved for compatibility).
#' @return List with \code{plots} (named ggplot objects), \code{plot_dims}
#'   (named width/height inches for HTML or PDF), \code{plot_paths}
#'   (when saved), and \code{selected_params_i}.
#' @export
integration_stability_plots_module <- function(
    stability = NULL,
    stability_dir = NULL,
    outdir = NULL,
    selected_params_i = NULL,
    label_top_k = 10L,
    save = !is.null(outdir),
    width = 10,
    height = 7
) {
  tabs <- .load_stability_tables(stability, stability_dir)
  perparam <- tabs$perparam_meanscores
  if (!nrow(perparam)) {
    stop("perparam_meanscores has no rows.", call. = FALSE)
  }
  perparam <- .prepare_stability_perparam(perparam)
  sel <- .stability_selected_params_i(perparam, selected_params_i)

  plots <- list(
    combinedscore_bar = .plot_stability_combinedscore_bar(perparam, sel),
    combinedscore_heatmap = .plot_stability_combinedscore_heatmap(perparam, sel),
    ari_jaccard_scatter = .plot_stability_ari_jaccard_scatter(
      perparam,
      sel,
      label_top_k = label_top_k
    ),
    bootstrap_ari = .plot_stability_bootstrap_ari(
      tabs$perbootstrap_perparam_scores,
      perparam
    ),
    nclust_ref = .plot_stability_nclust_ref(perparam, sel),
    perclust_jaccard = .plot_stability_perclust_jaccard(
      tabs$perclust_jaccard_mean_acrossbootstraps,
      sel
    )
  )

  plot_dims <- .stability_plot_save_dims(
    tabs$perparam_meanscores,
    tabs$perclust_jaccard_mean_acrossbootstraps,
    sel
  )

  plot_paths <- NULL
  if (isTRUE(save)) {
    if (is.null(outdir) || !nzchar(outdir)) {
      stop("outdir is required when save = TRUE.", call. = FALSE)
    }
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    plot_paths <- vapply(names(plots), function(nm) {
      path <- file.path(outdir, paste0("stability_", nm, ".pdf"))
      dims <- plot_dims[[nm]]
      .save_stability_ggplot(
        plots[[nm]],
        path,
        width = dims$width,
        height = dims$height
      )
      path
    }, character(1))
  }

  list(
    plots = plots,
    plot_dims = plot_dims,
    plot_paths = plot_paths,
    selected_params_i = sel
  )
}
