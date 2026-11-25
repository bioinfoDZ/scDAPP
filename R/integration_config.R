#' @title Integration method registry for scDAPP
#' @name integration_config
#' @description Configuration and validation for multi-method sample integration.
#'
#' Downstream DE policy: pseudobulk DE always uses \code{RNA}/\code{counts} (raw UMI
#' counts from the concatenated per-sample matrix). Wilcox DE uses integration-normalized
#' expression per method (\code{Integrated_RISC}/\code{data}, \code{SCT}/\code{data}, or
#' \code{RNA}/\code{data}). Seurat v5 \code{IntegrateLayers} methods store batch-corrected
#' signal in reductions, not a merged \code{integrated} expression assay.
NULL

#' Suffix appended to check_integration_dependencies() errors pointing to Usage.md.
#' @keywords internal
.INTEGRATION_DEPS_DOC <- "See Documentation/Usage.md (Integration method dependencies)."

#' Return allowed integration_method strings for match.arg and pipeline validation.
#' @export
integration_method_choices <- function() {
  c(
    "RISC",
    "CCAIntegration",
    "RPCAIntegration",
    "CCAIntegration_SCT",
    "RPCAIntegration_SCT",
    "HarmonyIntegration"
  )
}

# ---- Auto PC/resolution detection ----

#' Detect pcs_int or res_int == "auto" so config and run can branch into stability sweep.
#' @keywords internal
.is_auto_integration_param <- function(x) {
  is.character(x) && length(x) == 1L && tolower(x) == "auto"
}

#' TRUE when either pcs_int or res_int requests bootstrap auto-tuning.
#' @keywords internal
.integration_param_needs_auto <- function(pcs_int, res_int) {
  .is_auto_integration_param(pcs_int) || .is_auto_integration_param(res_int)
}


#' Build config list for SCT-based IntegrateLayers methods (shared by CCA_SCT and RPCA_SCT).
#' @param method Value stored in \code{config$method}.
#' @param integrate_layers_method \code{CCAIntegration} or \code{RPCAIntegration}.
#' @keywords internal
.sct_integrate_config <- function(
    method,
    integrate_layers_method,
    pcs_int,
    res_int,
    pcs_int_auto = FALSE,
    res_int_auto = FALSE
) {
  is_rpca <- integrate_layers_method == "RPCAIntegration"
  list(
    method = method,
    engine = "Seurat",
    normalization_method = "SCT",
    integrate_layers_method = integrate_layers_method,
    cluster_reduction = if (is_rpca) "integrated.dr" else "integrated.cca",
    umap_reduction = if (is_rpca) "umap.dr" else "umap.cca",
    integrated_assay = NULL,
    plot_assay = "SCT",
    de_pseudobulk_assay = "RNA",
    de_pseudobulk_slot = "counts",
    de_wilcox_assay = "SCT",
    de_wilcox_slot = "data",
    default_assay = "SCT",
    cluster_column_prefix = if (is_rpca) {
      "Integrated_SCT_RPCA"
    } else {
      "Integrated_SCT_CCA"
    },
    report_label = paste0(
      "Seurat v5 ",
      integrate_layers_method,
      " on merged split-RNA SCT (IntegrateLayers, normalization.method = SCT)"
    ),
    pcs_int = pcs_int,
    res_int = res_int,
    pcs_int_auto = pcs_int_auto,
    res_int_auto = res_int_auto
  )
}

#' Resolve per-method integration settings
#'
#' @param method Character; one of \code{integration_method_choices()}.
#' @param pcs_int Number of PCs for integration / clustering, or \code{"auto"}.
#' @param res_int Louvain resolution for integrated clustering, or \code{"auto"}.
#' @return Named list used by \code{run_integration()} and the pipeline Rmd.
#' @export
resolve_integration_config <- function(
    method,
    pcs_int = 30L,
    res_int = 0.5
) {
  method <- match.arg(method, integration_method_choices())
  pcs_auto <- .is_auto_integration_param(pcs_int)
  res_auto <- .is_auto_integration_param(res_int)
  if (pcs_auto || res_auto) {
    pcs_store <- if (pcs_auto) NA_integer_ else as.integer(pcs_int)[1]
    res_store <- if (res_auto) NA_real_ else as.numeric(res_int)[1]
  } else {
    pcs_store <- as.integer(pcs_int)[1]
    res_store <- as.numeric(res_int)[1]
  }
  pcs_int <- pcs_store
  res_int <- res_store

  base_seurat <- function(
    integrate_layers_method,
    cluster_reduction,
    umap_reduction,
    normalization_method = "RNA",
    plot_assay = "RNA",
    de_wilcox_assay = "RNA",
    default_assay = "RNA",
    report_label
  ) {
    list(
      method = method,
      engine = "Seurat",
      normalization_method = normalization_method,
      integrate_layers_method = integrate_layers_method,
      cluster_reduction = cluster_reduction,
      umap_reduction = umap_reduction,
      integrated_assay = NULL,
      plot_assay = plot_assay,
      de_pseudobulk_assay = "RNA",
      de_pseudobulk_slot = "counts",
      de_wilcox_assay = de_wilcox_assay,
      de_wilcox_slot = "data",
      default_assay = default_assay,
      cluster_column_prefix = paste0(
        "Integrated_",
        gsub("Integration$", "", integrate_layers_method)
      ),
      report_label = report_label,
      pcs_int = pcs_int,
      res_int = res_int,
      pcs_int_auto = pcs_auto,
      res_int_auto = res_auto
    )
  }

  if (method == "RISC") {
    return(list(
      method = method,
      engine = "RISC",
      normalization_method = NA_character_,
      integrate_layers_method = NA_character_,
      cluster_reduction = "pca",
      umap_reduction = "umap",
      integrated_assay = "Integrated_RISC",
      plot_assay = "Integrated_RISC",
      de_pseudobulk_assay = "RNA",
      de_pseudobulk_slot = "counts",
      de_wilcox_assay = "Integrated_RISC",
      de_wilcox_slot = "data",
      default_assay = "Integrated_RISC",
      cluster_column_prefix = "RISC_Louvain",
      report_label = "Robust Integration of Single Cell RNA-seq (RISC)",
      pcs_int = pcs_int,
      res_int = res_int,
      pcs_int_auto = pcs_auto,
      res_int_auto = res_auto
    ))
  }

  if (method == "RPCAIntegration_SCT") {
    return(.sct_integrate_config(
      method, "RPCAIntegration", pcs_int, res_int, pcs_auto, res_auto
    ))
  }
  if (method == "CCAIntegration_SCT") {
    return(.sct_integrate_config(
      method, "CCAIntegration", pcs_int, res_int, pcs_auto, res_auto
    ))
  }

  reduction_suffix <- gsub("Integration$", "", method)
  reduction_suffix <- tolower(reduction_suffix)
  if (method == "HarmonyIntegration") {
    reduction_suffix <- "harmony"
  }

  base_seurat(
    integrate_layers_method = method,
    cluster_reduction = if (method == "HarmonyIntegration") {
      "harmony"
    } else {
      paste0("integrated.", reduction_suffix)
    },
    umap_reduction = if (method == "HarmonyIntegration") {
      "umap.harmony"
    } else {
      paste0("umap.", reduction_suffix)
    },
    report_label = paste0("Seurat v5 ", method, " (IntegrateLayers, RNA normalization)")
  )
}

#' Metadata column name for integrated clusters, e.g. RISC_Louvain_npc30_res0.5.
#' @param config List from \code{resolve_integration_config()}.
#' @keywords internal
.integration_cluster_column <- function(config) {
  pcs <- as.integer(config$pcs_int)[1]
  res <- as.numeric(config$res_int)[1]
  paste0(config$cluster_column_prefix, "_npc", pcs, "_res", res)
}

#' Store integration_method and reduction names in sobjint@misc$scDAPP for downstream RDS.
#' @param sobjint Seurat object.
#' @param config List from \code{resolve_integration_config()}.
#' @keywords internal
.attach_integration_misc <- function(sobjint, config) {
  prev <- sobjint@misc$scDAPP
  if (is.null(prev)) prev <- list()
  sobjint@misc$scDAPP <- c(
    prev,
    list(
      integration_method = config$method,
      integration_reduction = config$cluster_reduction,
      normalization_method = config$normalization_method
    )
  )
  sobjint
}

#' Check optional packages for the chosen integration method
#' @param config List from \code{resolve_integration_config()}.
#' @export
check_integration_dependencies <- function(config, pcs_int = NULL, res_int = NULL) {
  missing_pkgs <- character(0)
  method <- config$method
  doc <- .INTEGRATION_DEPS_DOC

  needs_stability <- FALSE
  if (!is.null(pcs_int) && .is_auto_integration_param(pcs_int)) needs_stability <- TRUE
  if (!is.null(res_int) && .is_auto_integration_param(res_int)) needs_stability <- TRUE
  if (isTRUE(config$pcs_int_auto) || isTRUE(config$res_int_auto)) needs_stability <- TRUE
  if (needs_stability && !requireNamespace("mclust", quietly = TRUE)) {
    missing_pkgs <- c(missing_pkgs, "mclust (CRAN; required for pcs_int/res_int = 'auto')")
  }

  if (method == "HarmonyIntegration") {
    if (!requireNamespace("harmony", quietly = TRUE)) {
      missing_pkgs <- c(missing_pkgs, "harmony (CRAN; https://github.com/immunogenomics/harmony)")
    }
  }
  if (config$engine == "RISC") {
    if (!requireNamespace("RISC", quietly = TRUE)) {
      missing_pkgs <- c(missing_pkgs, "RISC")
    }
  }

  if (length(missing_pkgs)) {
    stop(
      "integration_method '", method, "' requires: ",
      paste(unique(missing_pkgs), collapse = ", "), ". ",
      doc,
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' Relabel clusters so ID 1 is the largest cluster for consistent plots and markers.
#' @param clusters Factor or vector of cluster labels.
#' @return Factor with numeric levels 1..K by size.
#' @keywords internal
.remap_clusters_by_size <- function(clusters) {
  clusters <- as.character(clusters)
  bs <- sort(table(clusters), decreasing = TRUE)
  out <- plyr::mapvalues(clusters, from = names(bs), to = seq_along(bs))
  factor(out, levels = seq_along(bs))
}

# ---- Input preparation helpers ----

#' Intersect genes across per-sample count matrices so cbind/matlist columns align.
#' @keywords internal
.intersect_matlist_genes <- function(matlist) {
  var0 <- Reduce(intersect, lapply(matlist, rownames))
  lapply(matlist, function(mat) mat[match(var0, rownames(mat)), , drop = FALSE])
}

#' Extract RNA counts per sample from sobjlist for RISC concat matrix (input_seurat_obj path).
#' @keywords internal
.build_matlist_from_sobjlist <- function(sobjlist, sample_metadata) {
  matlist <- lapply(sample_metadata$Code, function(code) {
    sobj <- sobjlist[[code]]
    if (is.null(sobj)) {
      stop("sobjlist missing Code: ", code, call. = FALSE)
    }
    Seurat::GetAssayData(sobj, assay = "RNA", layer = "counts")
  })
  names(matlist) <- sample_metadata$Code
  .intersect_matlist_genes(matlist)
}

#' Reload per-sample RDS from tmpobjdir and build matlist (stability / late RISC path).
#' @keywords internal
.build_matlist_from_tmpobjdir <- function(tmpobjdir, sample_metadata) {
  codes <- sample_metadata$Code
  sobjlist <- .read_sobjlist_from_dir(tmpobjdir, codes)
  .build_matlist_from_sobjlist(sobjlist, sample_metadata)
}

#' Read h5 counts per sample once at pipeline start; avoids repeated Cell Ranger IO.
#' @keywords internal
.build_matlist_from_h5_mdlist <- function(mdlist, sample_metadata, datadir) {
  if (is.null(datadir) || !nzchar(datadir)) {
    stop("datadir is required for h5-based count matrices.", call. = FALSE)
  }
  matlist <- lapply(sample_metadata$Code, function(code) {
    md <- mdlist[[code]]
    samp <- sample_metadata[sample_metadata$Code == code, "Sample", drop = TRUE]
    datafp <- paste0(datadir, "/", samp)
    h5_filename <- grep(
      pattern = "filtered_feature_bc_matrix.h5",
      list.files(datafp, recursive = TRUE, full.names = TRUE),
      value = TRUE
    )
    if (!length(h5_filename)) {
      stop("No filtered_feature_bc_matrix.h5 found under ", datafp, call. = FALSE)
    }
    mat0 <- Seurat::Read10X_h5(h5_filename[1])
    mat0[, match(rownames(md), colnames(mat0)), drop = FALSE]
  })
  names(matlist) <- sample_metadata$Code
  .intersect_matlist_genes(matlist)
}

#' Load cached .concatmatrix.rds or build matlist without re-reading datadir.
#' @keywords internal
.stability_load_matlist <- function(
    matlist_path = NULL,
    tmpobjdir = NULL,
    mdlist = NULL,
    sample_metadata = NULL,
    datadir = NULL,
    input_seurat_obj = FALSE
) {
  if (!is.null(matlist_path) && nzchar(matlist_path)) {
    found <- tryCatch(
      find_scdapp_object(matlist_path),
      error = function(e) NA_character_
    )
    if (!is.na(found) && file.exists(found)) {
      return(read_scdapp_object(found))
    }
  }
  if (isTRUE(input_seurat_obj)) {
    if (is.null(tmpobjdir) || !nzchar(tmpobjdir)) {
      stop(
        "tmpobjdir is required when input_seurat_obj is TRUE and matlist_path is missing.",
        call. = FALSE
      )
    }
    return(.build_matlist_from_tmpobjdir(tmpobjdir, sample_metadata))
  }
  .build_matlist_from_h5_mdlist(mdlist, sample_metadata, datadir)
}

#' Wrapper around .stability_load_matlist() for legacy call sites.
#' @keywords internal
.build_matlist_from_mdlist <- function(
    mdlist,
    sample_metadata,
    datadir = NULL,
    input_seurat_obj = FALSE,
    tmpobjdir = NULL,
    matlist_path = NULL
) {
  .stability_load_matlist(
    matlist_path = matlist_path,
    tmpobjdir = tmpobjdir,
    mdlist = mdlist,
    sample_metadata = sample_metadata,
    datadir = datadir,
    input_seurat_obj = input_seurat_obj
  )
}

#' Add Code, Sample, Condition, and Barcode columns before integration prep.
#' @keywords internal
.annotate_sobjlist_for_integration <- function(sobjlist, sample_metadata) {
  codes <- sample_metadata$Code
  sobjlist <- lapply(codes, function(code) {
    sample_metadata_code <- sample_metadata[sample_metadata$Code == code, , drop = FALSE]
    sobj <- sobjlist[[code]]
    sobj$orig.ident <- code
    sobj$Code <- code
    sobj$Sample <- sample_metadata_code$Sample
    sobj
  })
  names(sobjlist) <- codes

  sobjlist <- lapply(names(sobjlist), function(sampname) {
    sobj <- sobjlist[[sampname]]
    sobj$Condition <- sample_metadata[sample_metadata$Code == sampname, "Condition"]
    sobj
  })
  names(sobjlist) <- codes

  out <- lapply(names(sobjlist), function(sampname) {
    sobj <- sobjlist[[sampname]]
    md <- sobj@meta.data
    md <- cbind(Barcode = rownames(md), md)
    sobj@meta.data <- md
    sobj
  })
  names(out) <- codes
  out
}

#' Original per-sample cell ids.
#'
#' RISC stores these in \code{scBarcode}. Seurat \code{merge()} leaves them in
#' \code{Barcode} and renames colliding cell names (\code{AAAC..._1}).
#' @keywords internal
.original_cell_ids <- function(md) {
  if ("scBarcode" %in% colnames(md)) {
    return(as.character(md[["scBarcode"]]))
  }
  if ("Barcode" %in% colnames(md)) {
    return(as.character(md[["Barcode"]]))
  }
  rownames(md)
}

#' Per-sample query ids, before integration renamed cells.
#' @keywords internal
.sample_query_ids <- function(md) {
  if ("Barcode" %in% colnames(md)) {
    return(as.character(md[["Barcode"]]))
  }
  rownames(md)
}

.barcode_sample_label <- function(sample_code) {
  if (is.null(sample_code) || !nzchar(as.character(sample_code)[1])) {
    return("")
  }
  paste0(" (sample ", as.character(sample_code)[1], ")")
}

#' Integrated metadata rows in per-sample cell order.
#'
#' Matches each per-sample \code{Barcode} to integrated \code{scBarcode} when
#' that column exists, otherwise to \code{Barcode}. Row names of the result are
#' the integrated cell names. Stops when a hit is missing or duplicated.
#'
#' @param int_md Integrated metadata already restricted to one sample.
#' @param sample_md Per-sample metadata in object cell order.
#' @param sample_code Optional sample id used in error messages.
#' @keywords internal
.align_integrated_metadata_to_sample <- function(int_md, sample_md, sample_code = NULL) {
  where <- .barcode_sample_label(sample_code)
  query <- .sample_query_ids(sample_md)
  keys <- .original_cell_ids(int_md)
  if (length(query) != nrow(sample_md)) {
    stop("Per-sample cell ids do not match the metadata rows", where, ".", call. = FALSE)
  }
  if (anyDuplicated(keys)) {
    stop(
      "Integrated metadata has duplicate original cell ids",
      where,
      ".",
      call. = FALSE
    )
  }
  hit <- match(query, keys)
  if (anyNA(hit)) {
    stop(
      "Integrated metadata is missing ",
      sum(is.na(hit)),
      " of ",
      length(query),
      " per-sample cells",
      where,
      ".",
      call. = FALSE
    )
  }
  if (anyDuplicated(hit)) {
    stop(
      "Per-sample cell ids match the same integrated cell more than once",
      where,
      ".",
      call. = FALSE
    )
  }
  int_md[hit, , drop = FALSE]
}

#' Set prediction-assay column names to integrated cell names.
#' @keywords internal
.label_prediction_matrix <- function(predictionmat, integrated_cells, sample_code = NULL) {
  where <- .barcode_sample_label(sample_code)
  integrated_cells <- as.character(integrated_cells)
  if (ncol(predictionmat) != length(integrated_cells)) {
    stop(
      "Prediction assay has ",
      ncol(predictionmat),
      " columns and ",
      length(integrated_cells),
      " integrated cells",
      where,
      ".",
      call. = FALSE
    )
  }
  if (anyDuplicated(integrated_cells)) {
    stop("Integrated cell names are not unique", where, ".", call. = FALSE)
  }
  colnames(predictionmat) <- integrated_cells
  predictionmat
}

#' Column-bind per-sample prediction assays in integrated cell order.
#'
#' Rows are aligned to the first assay. Column names must already be integrated
#' cell names. Stops when an integrated cell is missing.
#' @keywords internal
.combine_prediction_matrices <- function(mats, integrated_cells) {
  if (!length(mats)) {
    stop("No prediction matrices to combine.", call. = FALSE)
  }
  integrated_cells <- as.character(integrated_cells)
  if (anyDuplicated(integrated_cells)) {
    stop("Integrated cell names are not unique.", call. = FALSE)
  }
  ref_rows <- rownames(mats[[1]])
  if (is.null(ref_rows) || !length(ref_rows)) {
    stop("Prediction assay is missing row names.", call. = FALSE)
  }
  mats <- lapply(mats, function(mat) {
    rows <- rownames(mat)
    if (is.null(rows) || !all(ref_rows %in% rows)) {
      stop("Prediction assays do not share the same rows.", call. = FALSE)
    }
    if (!identical(rows, ref_rows)) {
      mat <- mat[ref_rows, , drop = FALSE]
    }
    methods::as(mat, "CsparseMatrix")
  })
  predmat <- do.call(cbind, mats)
  rownames(predmat) <- ref_rows
  if (anyDuplicated(colnames(predmat))) {
    stop("Prediction assay column names are not unique.", call. = FALSE)
  }
  idx <- match(integrated_cells, colnames(predmat))
  if (length(colnames(predmat)) != length(integrated_cells) || anyNA(idx)) {
    n_miss <- if (anyNA(idx)) sum(is.na(idx)) else length(integrated_cells)
    stop(
      "Prediction assay is missing ",
      n_miss,
      " integrated cells.",
      call. = FALSE
    )
  }
  predmat[, idx, drop = FALSE]
}

#' Place each sample's count matrix onto integrated cell names.
#'
#' Matches original ids (\code{scBarcode}, else \code{Barcode}) within
#' \code{Code}. A shared 10x barcode is not looked up in the concatenated matrix.
#' @keywords internal
.align_concat_matlist_to_object <- function(matlist, sobjint) {
  if (!length(matlist)) {
    stop("concat matrix is empty.", call. = FALSE)
  }
  bigmat <- do.call(cbind, matlist)
  num_nonzeros <- tabulate(bigmat@i + 1L, nbins = nrow(bigmat))
  min_samples <- min(3L, length(matlist))
  genes <- rownames(bigmat)[num_nonzeros >= min_samples]
  rm(bigmat)
  if (!length(genes)) {
    stop("concat matrix has no genes expressed in enough cells.", call. = FALSE)
  }

  md <- sobjint@meta.data
  if (!"Code" %in% colnames(md)) {
    stop(
      "Integrated object is missing Code; cannot align concat counts by sample.",
      call. = FALSE
    )
  }
  integrated_cells <- colnames(sobjint)
  if (anyDuplicated(integrated_cells)) {
    stop("Integrated cell names are not unique.", call. = FALSE)
  }
  orig <- .original_cell_ids(md)
  if (length(orig) != length(integrated_cells)) {
    stop("Original cell ids do not match integrated cells.", call. = FALSE)
  }
  codes <- as.character(md[["Code"]])
  if (anyNA(codes) || any(!nzchar(codes))) {
    stop("Integrated object has cells with no Code.", call. = FALSE)
  }

  pieces <- lapply(unique(codes), function(code) {
    cells <- which(codes == code)
    mat <- matlist[[code]]
    if (is.null(mat)) {
      stop("concat matrix is missing sample ", code, ".", call. = FALSE)
    }
    gene_hit <- match(genes, rownames(mat))
    if (anyNA(gene_hit)) {
      stop(
        "concat matrix is missing ",
        sum(is.na(gene_hit)),
        " genes for sample ",
        code,
        ".",
        call. = FALSE
      )
    }
    mat <- mat[gene_hit, , drop = FALSE]
    keys <- colnames(mat)
    if (anyDuplicated(keys)) {
      stop("concat matrix has duplicate cell names in sample ", code, ".", call. = FALSE)
    }
    query <- orig[cells]
    if (anyDuplicated(query)) {
      stop(
        "Integrated metadata has duplicate original cell ids (sample ",
        code,
        ").",
        call. = FALSE
      )
    }
    hit <- match(query, keys)
    if (anyNA(hit)) {
      stop(
        "concat matrix missing ",
        sum(is.na(hit)),
        " of ",
        length(cells),
        " cells in sample ",
        code,
        ". Check Barcode metadata alignment.",
        call. = FALSE
      )
    }
    mat <- mat[, hit, drop = FALSE]
    colnames(mat) <- integrated_cells[cells]
    mat
  })
  bigmat <- do.call(cbind, pieces)
  bigmat[, integrated_cells, drop = FALSE]
}

#' Attach joined raw RNA counts assay to integrated object for pseudobulk DE.
#' @keywords internal
.append_concat_rna_assay <- function(
    sobjint,
    outdir_int_objects,
    sample_metadata,
    saving_algorithm = "qs2"
) {
  matlist_path <- find_scdapp_object(
    file.path(outdir_int_objects, ".concatmatrix"),
    saving_algorithm
  )
  matlist <- read_scdapp_object(matlist_path)
  bigmat <- .align_concat_matlist_to_object(matlist, sobjint)
  rm(matlist)

  rnaassay <- Seurat::CreateAssayObject(counts = bigmat)
  sobjint[["RNA"]] <- rnaassay
  sobjint <- Seurat::NormalizeData(sobjint, assay = "RNA", verbose = FALSE)

  risc_assay <- "Integrated_RISC"
  if (risc_assay %in% Seurat::Assays(sobjint)) {
    sobjint@assays[[risc_assay]]@counts <- expm1(sobjint@assays[[risc_assay]]@data)
  }

  unlink(matlist_path)
  sobjint
}

#' Set active ident, default assay, and misc metadata after integration clustering.
#' @keywords internal
.finalize_sobjint_clusters <- function(sobjint, config, int_clust_lab) {
  sobjint <- Seurat::SetIdent(sobjint, value = sobjint@meta.data[[int_clust_lab]])
  sobjint$seurat_clusters <- sobjint@meta.data[[int_clust_lab]]
  Seurat::DefaultAssay(sobjint) <- config$default_assay
  sobjint <- .attach_integration_misc(sobjint, config)
  sobjint
}

#' Load Code.rds objects from pipeline tmpobjdir for Seurat merge integration.
#' @keywords internal
.read_sobjlist_from_dir <- function(objdir, codes, saving_algorithm = "qs2") {
  out <- lapply(codes, function(code) {
    read_scdapp_object(
      find_scdapp_object(file.path(objdir, code), saving_algorithm)
    )
  })
  names(out) <- codes
  out
}
