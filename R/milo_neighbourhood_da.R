# Neighbourhood differential abundance (Milo)
#
# Graph and neighbourhoods are built once on the integrated clustering embedding.
# Each comps row is one contrast, using the same formula and contrast strings as
# EdgeR-QLF. Gene expression is not used.

.milo_require <- function() {
  if (!requireNamespace("miloR", quietly = TRUE)) {
    stop(
      "miloR is required for milo_neighbourhood_da_module(). ",
      "Install it with BiocManager::install(\"miloR\").",
      call. = FALSE
    )
  }
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop(
      "SingleCellExperiment is required for milo_neighbourhood_da_module(). ",
      "It is installed with miloR.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.milo_file_stub <- function(label) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(label)[1])
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if (!nzchar(x)) x <- "comparison"
  x
}

.milo_resolve_k <- function(milo_k, milo_k_max, n_samples, n_cells) {
  n_cap <- max(1L, as.integer(n_cells) - 1L)
  milo_k_max <- as.integer(milo_k_max)[1]
  if (length(milo_k_max) != 1L || is.na(milo_k_max) || milo_k_max < 1L) {
    stop("milo_k_max must be a positive integer.", call. = FALSE)
  }
  auto <- is.character(milo_k) && length(milo_k) == 1L &&
    tolower(trimws(milo_k)) == "auto"
  if (auto) {
    uncapped <- 5L * as.integer(n_samples)
    k <- min(uncapped, milo_k_max, n_cap)
    return(list(k = k, uncapped = uncapped, capped = unname(k < uncapped), auto = TRUE))
  }
  k <- suppressWarnings(as.integer(milo_k)[1])
  if (length(k) != 1L || is.na(k) || k < 1L) {
    stop("milo_k must be \"auto\" or a positive integer.", call. = FALSE)
  }
  if (k > n_cap) {
    stop(
      "milo_k (", k, ") must be smaller than the number of cells (", n_cells, ").",
      call. = FALSE
    )
  }
  list(k = k, uncapped = k, capped = FALSE, auto = FALSE)
}

# Numeric contrast aligned to a treatment-coded model matrix.
# DESeq2 triples and name= strings are translated onto EdgeR column names.
.milo_contrast_vector <- function(design, parsed) {
  cn <- colnames(design)
  if (identical(parsed$style, "ignore")) {
    stop(
      "Milo uses the pseudobulk contrast grammar. DE_test is a single-cell test, ",
      "which ignores comps$contrast.",
      call. = FALSE
    )
  }
  if (identical(parsed$style, "deseq2") && !is.null(parsed$deseq_name)) {
    nm <- parsed$deseq_name
    candidates <- unique(c(
      nm,
      gsub(".", ":", nm, fixed = TRUE),
      gsub(":", ".", nm, fixed = TRUE)
    ))
    hit <- candidates[candidates %in% cn]
    if (!length(hit)) {
      stop(
        "Could not map DESeq2 name='", nm, "' to design columns: ",
        paste(cn, collapse = ", "),
        call. = FALSE
      )
    }
    cv <- rep(0, length(cn))
    cv[match(hit[1], cn)] <- 1
    return(cv)
  }
  expr <- if (identical(parsed$style, "deseq2") && !is.null(parsed$deseq_triple)) {
    trip <- parsed$deseq_triple
    paste0(trip[[1]], trip[[2]], " - ", trip[[1]], trip[[3]])
  } else {
    parsed$edger_expr
  }
  if (is.null(expr) || !nzchar(as.character(expr)[1])) {
    stop("Empty contrast expression.", call. = FALSE)
  }
  .edger_expr_to_vector(design, expr)
}

# makeContrasts string for the non-zero entries of a contrast vector.
# The pipeline parser only emits a single coefficient or an A - B pair.
.milo_contrast_string <- function(cv, cn) {
  nz <- which(abs(cv) > 1e-8)
  if (!length(nz)) stop("Contrast is all zeros.", call. = FALSE)
  if (!all(abs(abs(cv[nz]) - 1) < 1e-8)) {
    stop(
      "Milo contrasts support a single coefficient or an A - B difference. Got coefficients: ",
      paste(paste0(cn[nz], "=", cv[nz]), collapse = ", "),
      call. = FALSE
    )
  }
  pos <- cn[nz][cv[nz] > 0]
  neg <- cn[nz][cv[nz] < 0]
  if (!length(neg)) return(paste(pos, collapse = " + "))
  paste(c(paste(pos, collapse = " + "), paste("-", neg)), collapse = " ")
}

.milo_tests_last_coef <- function(cv) {
  p <- length(cv)
  p >= 1L && abs(cv[p] - 1) < 1e-8 && sum(abs(cv) > 1e-8) == 1L
}

.milo_design_df <- function(sample_metadata, formula, sample_ids) {
  md <- .prepare_pseudobulk_coldata(sample_metadata, formula)
  md$Code <- as.character(sample_metadata$Code)
  if (anyDuplicated(md$Code)) {
    stop("sample_metadata$Code values must be unique for Milo.", call. = FALSE)
  }
  rownames(md) <- md$Code
  miss <- setdiff(sample_ids, rownames(md))
  if (length(miss)) {
    stop(
      "Neighbourhood counts include samples missing from sample_metadata$Code: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  md <- md[sample_ids, , drop = FALSE]
  vars <- all.vars(formula)
  for (v in vars) {
    if (is.factor(md[[v]])) md[[v]] <- droplevels(md[[v]])
  }
  md
}

.milo_standardize_da <- function(da, cluster_col) {
  da <- as.data.frame(da, stringsAsFactors = FALSE)
  lfc_col <- grep("^logFC", colnames(da), value = TRUE)
  if (!length(lfc_col)) stop("Milo result is missing logFC.", call. = FALSE)
  if (!"SpatialFDR" %in% colnames(da)) {
    stop("Milo result is missing SpatialFDR.", call. = FALSE)
  }
  if (!"PValue" %in% colnames(da)) da$PValue <- NA_real_
  if (!"FDR" %in% colnames(da)) {
    da$FDR <- stats::p.adjust(da$PValue, method = "BH")
  }
  if (!"F" %in% colnames(da)) da$F <- NA_real_
  if (!"logCPM" %in% colnames(da)) da$logCPM <- NA_real_
  if (!"Nhood" %in% colnames(da)) da$Nhood <- seq_len(nrow(da))
  frac_col <- paste0(cluster_col, "_fraction")
  cluster <- if (cluster_col %in% colnames(da)) as.character(da[[cluster_col]]) else rep(NA_character_, nrow(da))
  fraction <- if (frac_col %in% colnames(da)) as.numeric(da[[frac_col]]) else rep(NA_real_, nrow(da))
  data.frame(
    Nhood = as.integer(da$Nhood),
    logFC = as.numeric(da[[lfc_col[1]]]),
    logCPM = as.numeric(da$logCPM),
    F = as.numeric(da$F),
    PValue = as.numeric(da$PValue),
    FDR = as.numeric(da$FDR),
    SpatialFDR = as.numeric(da$SpatialFDR),
    cluster = cluster,
    cluster_fraction = fraction,
    stringsAsFactors = FALSE
  )
}

.milo_lfc_color_scale <- function() {
  ggplot2::scale_color_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    na.value = "grey80",
    name = "Log fold change"
  )
}

.milo_ggsave <- function(plot, path) {
  ggplot2::ggsave(filename = path, plot = plot, width = 9, height = 7, units = "in")
  invisible(path)
}

.milo_save_heatmap <- function(ht, path) {
  grDevices::pdf(path, width = 11, height = 8)
  on.exit(grDevices::dev.off(), add = TRUE)
  ComplexHeatmap::draw(ht)
  invisible(path)
}

.milo_overlap_long <- function(milo, grouping_variable) {
  nh <- miloR::nhoods(milo)
  cl <- as.character(SummarizedExperiment::colData(milo)[[grouping_variable]])
  parts <- vector("list", ncol(nh))
  for (j in seq_len(ncol(nh))) {
    cells_j <- which(nh[, j] != 0)
    if (!length(cells_j)) next
    tab <- table(cl[cells_j])
    parts[[j]] <- data.frame(
      Nhood = j,
      cluster = names(tab),
      n_cells = as.integer(tab),
      stringsAsFactors = FALSE
    )
  }
  parts <- parts[lengths(parts) > 0]
  if (!length(parts)) {
    return(data.frame(Nhood = integer(), cluster = character(), n_cells = integer()))
  }
  do.call(rbind, parts)
}

.milo_top_labels <- function(milo, grouping_variable, nhood_ids) {
  if (!length(nhood_ids)) {
    return(data.frame(Cluster = character(), Nhood = character()))
  }
  nh <- miloR::nhoods(milo)
  cl <- as.character(SummarizedExperiment::colData(milo)[[grouping_variable]])
  parts <- vector("list", length(nhood_ids))
  for (i in seq_along(nhood_ids)) {
    j <- nhood_ids[i]
    cells_j <- which(nh[, j] != 0)
    if (!length(cells_j)) next
    parts[[i]] <- data.frame(
      Cluster = cl[cells_j],
      Nhood = as.character(j),
      stringsAsFactors = FALSE
    )
  }
  parts <- parts[lengths(parts) > 0]
  if (!length(parts)) {
    return(data.frame(Cluster = character(), Nhood = character()))
  }
  labelsdf <- do.call(rbind, parts)
  labelsdf$Nhood <- factor(labelsdf$Nhood, levels = as.character(nhood_ids))
  labelsdf$Cluster <- factor(
    labelsdf$Cluster,
    levels = .milo_cluster_level_order(labelsdf$Cluster)
  )
  labelsdf
}

# Numeric cluster ids sort as 1, 2, 3. Other labels sort alphabetically.
.milo_cluster_level_order <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(x)
  num <- suppressWarnings(as.numeric(x))
  if (all(!is.na(num))) x[order(num)] else x[order(x)]
}

# ggplot builds some warnings only when drawn. Keep them in the console log.
.milo_muffle_to_log <- function(expr) {
  withCallingHandlers(
    force(expr),
    warning = function(w) {
      cat(conditionMessage(w), "\n", file = stderr())
      invokeRestart("muffleWarning")
    }
  )
}

.milo_cluster_summary <- function(da, spatial_fdr) {
  sig <- da[!is.na(da$SpatialFDR) & da$SpatialFDR < spatial_fdr, , drop = FALSE]
  clusters <- .milo_cluster_level_order(da$cluster)
  if (!length(clusters)) {
    return(data.frame(
      cluster = character(),
      n_significant = integer(),
      n_enriched = integer(),
      n_depleted = integer(),
      stringsAsFactors = FALSE
    ))
  }
  rows <- lapply(clusters, function(cl) {
    sub <- sig[sig$cluster == cl, , drop = FALSE]
    data.frame(
      cluster = cl,
      n_significant = nrow(sub),
      n_enriched = sum(sub$logFC > 0, na.rm = TRUE),
      n_depleted = sum(sub$logFC < 0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.milo_call_test <- function(milo, design, design.df, model_contrasts, force) {
  args <- list(
    x = milo,
    design = design,
    design.df = design.df,
    fdr.weighting = "graph-overlap",
    norm.method = "TMM",
    robust = TRUE,
    reduced.dim = "PCA"
  )
  if (!is.null(model_contrasts)) args$model.contrasts <- model_contrasts
  if (isTRUE(force)) args$force <- TRUE
  do.call(miloR::testNhoods, args)
}

.milo_one_comparison <- function(milo, design.df, formula, parsed, spatial_fdr) {
  fixed <- .fixed_effects_formula(formula)
  re_vars <- .formula_random_effect_vars(formula)
  has_re <- length(re_vars) > 0L
  design <- stats::model.matrix(fixed, data = design.df)
  if (qr(design)$rank < ncol(design)) {
    stop("Design is not full rank. Columns: ", paste(colnames(design), collapse = ", "), call. = FALSE)
  }
  resid_df <- nrow(design) - qr(design)$rank
  if (resid_df < 1L) {
    stop(
      "No residual degrees of freedom (", nrow(design), " samples, ",
      qr(design)$rank, " coefficients). Milo needs biological replicates.",
      call. = FALSE
    )
  }
  cv <- .milo_contrast_vector(design, parsed)
  # Interaction columns contain ':'. limma::makeContrasts rejects those names.
  # Rename after the contrast is matched so the tested columns stay the same.
  contrast_label <- .milo_contrast_string(cv, colnames(design))
  colnames(design) <- make.names(colnames(design), unique = TRUE)
  rownames(design) <- rownames(design.df)
  contrast_string <- .milo_contrast_string(cv, colnames(design))
  pairing_status <- "none"
  pairing_note <- ""
  model_contrast <- contrast_string

  if (!has_re) {
    raw <- .milo_call_test(milo, design, design.df, contrast_string, force = FALSE)
  } else if (!.milo_tests_last_coef(cv)) {
    pairing_status <- "fixed_fallback"
    pairing_note <- paste0(
      "Milo's mixed model only tests the last fixed coefficient and ignores contrasts. ",
      "This contrast (", contrast_label, ") was fit without the random effect (",
      paste(paste0("(1|", re_vars, ")"), collapse = " + "),
      ")."
    )
    warning(pairing_note, call. = FALSE)
    raw <- .milo_call_test(milo, design, design.df, contrast_string, force = FALSE)
  } else {
    n_s <- nrow(design.df)
    if (n_s < 60L) {
      warning(
        "Milo mixed model with N=", n_s,
        " and force=TRUE. Estimates can be unstable below 60 samples.",
        call. = FALSE
      )
    }
    glmm <- tryCatch(
      .milo_call_test(milo, formula, design.df, model_contrasts = NULL, force = TRUE),
      error = function(e) e
    )
    if (inherits(glmm, "error")) {
      pairing_status <- "fixed_fallback"
      pairing_note <- paste0(
        "Milo mixed model failed (", conditionMessage(glmm),
        "). Refit without ",
        paste(paste0("(1|", re_vars, ")"), collapse = " + "),
        ". Pairing was dropped."
      )
      warning(pairing_note, call. = FALSE)
      raw <- .milo_call_test(milo, design, design.df, contrast_string, force = FALSE)
    } else {
      pairing_status <- "glmm"
      pairing_note <- paste0(
        "Fit with Milo's negative-binomial mixed model (force=TRUE) using ",
        paste(paste0("(1|", re_vars, ")"), collapse = " + "),
        ". Mixed-model estimates with fewer than 60 samples can be unstable."
      )
      model_contrast <- "glmm_last_coef"
      raw <- glmm
    }
  }

  list(
    raw = raw,
    pairing_status = pairing_status,
    pairing_note = pairing_note,
    model_contrast = model_contrast,
    contrast_string = contrast_string
  )
}

.milo_build_object <- function(sobjint, grouping_variable, cluster_reduction,
                               umap_reduction, d, k, seed) {
  if (!cluster_reduction %in% names(sobjint@reductions)) {
    stop(
      "Clustering reduction '", cluster_reduction, "' was not found. Available: ",
      paste(names(sobjint@reductions), collapse = ", "),
      call. = FALSE
    )
  }
  if (!umap_reduction %in% names(sobjint@reductions)) {
    stop(
      "UMAP reduction '", umap_reduction, "' was not found. Available: ",
      paste(names(sobjint@reductions), collapse = ", "),
      call. = FALSE
    )
  }
  cells <- colnames(sobjint)
  emb <- Seurat::Embeddings(sobjint, cluster_reduction)
  if (!all(cells %in% rownames(emb))) {
    stop("Clustering embedding is missing cells from the integrated object.", call. = FALSE)
  }
  d_req <- suppressWarnings(as.integer(d)[1])
  if (length(d_req) != 1L || is.na(d_req) || d_req < 1L) {
    d_use <- ncol(emb)
  } else {
    d_use <- min(d_req, ncol(emb))
  }
  if (d_use < 2L) {
    stop("Milo needs at least 2 embedding dimensions.", call. = FALSE)
  }
  if (k >= nrow(emb)) {
    stop("k must be smaller than the number of cells.", call. = FALSE)
  }
  pca <- as.matrix(emb[cells, seq_len(d_use), drop = FALSE])
  umap <- as.matrix(Seurat::Embeddings(sobjint, umap_reduction)[cells, , drop = FALSE])
  if (ncol(umap) > 2L) umap <- umap[, seq_len(2L), drop = FALSE]
  md <- as.data.frame(sobjint@meta.data[cells, , drop = FALSE])
  keep <- vapply(md, function(col) is.atomic(col) || is.factor(col), logical(1))
  md <- md[, keep, drop = FALSE]
  if (!grouping_variable %in% colnames(md)) {
    stop("Cluster column '", grouping_variable, "' is not in the integrated metadata.", call. = FALSE)
  }
  if (!"Code" %in% colnames(md)) {
    stop("Integrated metadata must include Code (sample id) for Milo.", call. = FALSE)
  }
  md$Code <- as.character(md$Code)
  rownames(md) <- cells

  dummy <- Matrix::Matrix(
    0, nrow = 1L, ncol = length(cells), sparse = TRUE,
    dimnames = list("milo_placeholder", cells)
  )
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = dummy, counts = dummy)
  )
  SingleCellExperiment::reducedDim(sce, "PCA") <- pca
  SingleCellExperiment::reducedDim(sce, "UMAP") <- umap
  SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(md, check.names = FALSE)

  set.seed(as.integer(seed)[1])
  milo <- miloR::Milo(sce)
  milo <- miloR::buildGraph(milo, k = k, d = d_use, reduced.dim = "PCA")
  milo <- miloR::makeNhoods(
    milo, prop = 0.1, k = k, d = d_use, refined = TRUE,
    reduced_dims = "PCA", refinement_scheme = "reduced_dim"
  )
  meta <- as.data.frame(SummarizedExperiment::colData(milo))
  rownames(meta) <- colnames(milo)
  milo <- miloR::countCells(milo, samples = "Code", meta.data = meta)
  dist_try <- tryCatch(
    miloR::calcNhoodDistance(milo, d = d_use, reduced.dim = "PCA"),
    error = function(e) e
  )
  if (inherits(dist_try, "error")) {
    warning(
      "calcNhoodDistance failed (", conditionMessage(dist_try),
      "). Spatial FDR uses graph-overlap, which does not need those distances.",
      call. = FALSE
    )
  } else {
    milo <- dist_try
  }
  milo <- miloR::buildNhoodGraph(milo)
  list(milo = milo, d_use = d_use)
}

#' Neighbourhood differential abundance with Milo
#'
#' Builds a kNN graph on the integrated clustering embedding, defines overlapping
#' neighbourhoods, and tests sample-level neighbourhood counts with Milo's
#' edgeR quasi-likelihood model. One comps row is one contrast. The graph is
#' built once and reused for every comparison.
#'
#' Automatic \code{milo_k} is \code{min(5 * n_samples, milo_k_max, n_cells - 1)}.
#' An integer \code{milo_k} replaces that value. Paired formulas are fit with
#' Milo's mixed model and \code{force = TRUE}. If that model errors, or if the
#' contrast is not the last fixed coefficient (the only term the mixed model
#' returns), the fit is repeated without \code{(1|var)} and a warning is issued.
#'
#' @param sobjint Integrated Seurat object with \code{Code} and the cluster column.
#' @param comps Comparison table (\code{c0}, \code{c1}, optional \code{formula},
#'   \code{contrast}, \code{label}). See \code{normalize_comps()}.
#' @param sample_metadata Sample table with \code{Code}, \code{Condition}, and
#'   any covariates named in \code{comps$formula}.
#' @param outdir_int Integrated output directory. Results are written to
#'   \code{differential_abundance_Milo/} inside it.
#' @param grouping_variable Metadata column of integrated clusters.
#' @param cluster_reduction Name of the Seurat reduction used for clustering.
#' @param umap_reduction Name of the Seurat UMAP reduction.
#' @param d Number of embedding dimensions. Defaults to all columns of
#'   \code{cluster_reduction} when missing or non-numeric. Capped at the number
#'   of columns in that reduction.
#' @param DE_test Contrast grammar. EdgeR / Dream use coefficient expressions.
#'   DESeq2 triples and \code{name=} are translated onto those coefficients.
#' @param milo_k \code{"auto"} or a positive integer.
#' @param milo_k_max Cap used when \code{milo_k = "auto"}. Default 100.
#' @param spatial_fdr SpatialFDR cutoff for significant neighbourhoods. Default 0.1.
#' @param n_report How many significant neighbourhoods to flag for the report
#'   (tables, alluvial, two-way heatmap). Default 20. Ranked by SpatialFDR,
#'   then absolute log fold change.
#' @param seed Random seed for neighbourhood sampling. Default 1.
#'
#' @return A list with \code{run_info}, \code{size_plot}, \code{nhood_counts},
#'   \code{overlap}, \code{outdir}, and \code{comparisons} (named by comps label).
#'   Each comparison has the full DA table, the top-report subset, cluster
#'   summary, plots, and \code{pairing_status} (\code{none}, \code{glmm}, or
#'   \code{fixed_fallback}).
#' @export
milo_neighbourhood_da_module <- function(
    sobjint,
    comps,
    sample_metadata,
    outdir_int,
    grouping_variable = "seurat_clusters",
    cluster_reduction,
    umap_reduction,
    d = NULL,
    DE_test = "EdgeR-LRT",
    milo_k = "auto",
    milo_k_max = 100L,
    spatial_fdr = 0.1,
    n_report = 20L,
    seed = 1L
) {
  .milo_require()
  if (missing(cluster_reduction) || missing(umap_reduction)) {
    stop("cluster_reduction and umap_reduction are required.", call. = FALSE)
  }
  comps <- .apply_contrast_defaults(.normalize_comps(comps), DE_test)
  style <- .contrast_style(DE_test)
  sample_metadata <- as.data.frame(sample_metadata, stringsAsFactors = FALSE)
  if (!"Code" %in% colnames(sample_metadata)) {
    stop("sample_metadata must include Code.", call. = FALSE)
  }
  sample_metadata$Code <- as.character(sample_metadata$Code)

  sample_n <- table(as.character(sobjint$Code))
  n_samples <- length(sample_n)
  n_cells <- ncol(sobjint)
  if (n_samples < 1L || n_cells < 3L) {
    stop("Milo needs cells from at least one sample.", call. = FALSE)
  }
  k_info <- .milo_resolve_k(milo_k, milo_k_max, n_samples, n_cells)
  n_small <- sum(sample_n < 5L)
  if (n_small > 0L) {
    warning(
      n_small, " sample(s) have fewer than 5 cells. ",
      "Neighbourhoods may contain fewer than 5 cells from those samples. ",
      "k is not reduced to the smallest sample.",
      call. = FALSE
    )
  }
  if (isTRUE(k_info$capped)) {
    warning(
      "Automatic Milo k was capped at ", k_info$k,
      " (5 * n_samples = ", k_info$uncapped,
      ", milo_k_max = ", as.integer(milo_k_max), "). ",
      "Pass an integer milo_k to use a larger neighbourhood.",
      call. = FALSE
    )
  }

  built <- .milo_build_object(
    sobjint, grouping_variable, cluster_reduction, umap_reduction,
    d = d, k = k_info$k, seed = seed
  )
  milo <- built$milo
  outdir <- file.path(outdir_int, "differential_abundance_Milo")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  n_nhoods <- ncol(miloR::nhoods(milo))
  run_info <- data.frame(
    k = k_info$k,
    k_uncapped = k_info$uncapped,
    k_capped = k_info$capped,
    k_auto = k_info$auto,
    milo_k_max = as.integer(milo_k_max),
    d = built$d_use,
    n_cells = n_cells,
    n_samples = n_samples,
    n_nhoods = n_nhoods,
    n_samples_under_5_cells = n_small,
    spatial_fdr = spatial_fdr,
    n_report = as.integer(n_report),
    seed = as.integer(seed),
    stringsAsFactors = FALSE
  )
  utils::write.csv(run_info, file.path(outdir, "milo_run_info.csv"), row.names = FALSE)

  size_plot <- miloR::plotNhoodSizeHist(milo) +
    ggplot2::ggtitle(paste0("Neighborhood size (k = ", k_info$k, ", d = ", built$d_use, ")"))
  .milo_ggsave(size_plot, file.path(outdir, "nhood_size_histogram.pdf"))

  nhood_counts <- as.data.frame(as.matrix(miloR::nhoodCounts(milo)))
  nhood_counts <- cbind(Nhood = seq_len(nrow(nhood_counts)), nhood_counts)
  utils::write.csv(nhood_counts, file.path(outdir, "nhood_counts.csv"), row.names = FALSE)

  overlap <- .milo_overlap_long(milo, grouping_variable)
  utils::write.csv(overlap, file.path(outdir, "nhood_cluster_overlap.csv"), row.names = FALSE)

  counts_mat <- miloR::nhoodCounts(milo)
  lib_size <- Matrix::colSums(counts_mat)
  sample_ids <- names(lib_size)[lib_size > 0]
  if (!length(sample_ids)) {
    stop("Every sample has zero cells in the neighbourhood counts.", call. = FALSE)
  }

  comparisons <- setNames(vector("list", nrow(comps)), comps$label)
  for (i in seq_len(nrow(comps))) {
    lab <- comps$label[i]
    message("Milo comparison: ", lab)
    stub <- .milo_file_stub(lab)
    comp_dir <- file.path(outdir, stub)
    dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
    one <- tryCatch({
      formula <- .parse_comp_formula(comps$formula[i])
      parsed <- .parse_contrast(
        comps$contrast[i], comps$c0[i], comps$c1[i], style = style
      )
      design.df <- .milo_design_df(sample_metadata, formula, sample_ids)
      fit <- .milo_one_comparison(milo, design.df, formula, parsed, spatial_fdr)
      raw <- as.data.frame(fit$raw, stringsAsFactors = FALSE)
      if ("Nhood" %in% colnames(raw)) {
        raw <- raw[order(as.integer(raw$Nhood)), , drop = FALSE]
        rownames(raw) <- NULL
      }
      raw <- miloR::annotateNhoods(milo, raw, coldata_col = grouping_variable)
      da <- .milo_standardize_da(raw, grouping_variable)
      sizes <- Matrix::colSums(miloR::nhoods(milo))
      da$n_cells <- as.integer(sizes[da$Nhood])
      da$significant <- !is.na(da$SpatialFDR) & da$SpatialFDR < spatial_fdr
      da$label <- lab
      da$formula <- comps$formula[i]
      da$contrast <- parsed$display
      da$model_contrast <- fit$model_contrast
      da$pairing_status <- fit$pairing_status
      da$c0 <- as.character(comps$c0[i])
      da$c1 <- as.character(comps$c1[i])

      sig <- da[da$significant, , drop = FALSE]
      if (nrow(sig)) {
        sig <- sig[order(sig$SpatialFDR, -abs(sig$logFC), sig$Nhood), , drop = FALSE]
      }
      n_keep <- min(as.integer(n_report), nrow(sig))
      da_top <- if (n_keep > 0L) sig[seq_len(n_keep), , drop = FALSE] else sig
      summary_df <- .milo_cluster_summary(da, spatial_fdr)

      # Milo's beeswarm drops rows with NA logFC, then mis-aligns positions.
      # Plot non-converged neighbourhoods as logFC 0 / not significant.
      da_plot <- da
      da_plot$logFC[is.na(da_plot$logFC)] <- 0
      da_plot$logCPM[is.na(da_plot$logCPM)] <- 0
      da_plot$SpatialFDR[is.na(da_plot$SpatialFDR)] <- 1
      da_plot$cluster[is.na(da_plot$cluster) | !nzchar(da_plot$cluster)] <- "unassigned"
      da_plot$cluster <- factor(
        da_plot$cluster,
        levels = .milo_cluster_level_order(da_plot$cluster)
      )
      beeswarm <- tryCatch(
        miloR::plotDAbeeswarm(da_plot, group.by = "cluster", alpha = spatial_fdr) +
          .milo_lfc_color_scale() +
          ggplot2::labs(
            x = "Integrated cluster",
            y = "Log fold change",
            title = lab
          ) +
          ggplot2::guides(
            color = ggplot2::guide_colorbar(title = "Log fold change")
          ),
        error = function(e) {
          warning("Beeswarm failed for ", lab, ": ", conditionMessage(e), call. = FALSE)
          ggplot2::ggplot(
            da_plot,
            ggplot2::aes(x = logFC, y = cluster, color = logFC)
          ) +
            ggplot2::geom_jitter(width = 0, height = 0.15, alpha = 0.7) +
            .milo_lfc_color_scale() +
            ggplot2::labs(
              x = "Log fold change",
              y = "Integrated cluster",
              title = lab
            ) +
            ggplot2::theme_classic()
        }
      )
      milo_plot <- miloR::buildNhoodGraph(milo)
      graph <- miloR::plotNhoodGraphDA(milo_plot, da_plot, alpha = spatial_fdr) +
        ggplot2::ggtitle(lab)
      ma <- miloR::plotNhoodMA(da_plot, alpha = spatial_fdr) +
        ggplot2::ggtitle(lab)

      labelsdf <- .milo_top_labels(milo, grouping_variable, da_top$Nhood)
      alluvial <- NULL
      twt <- NULL
      if (nrow(labelsdf) > 0L && length(unique(labelsdf$Nhood)) >= 1L) {
        alluvial <- alluvialplot(labelsdf) +
          ggplot2::ggtitle(paste0(lab, " (top neighborhoods)"))
        twt <- twt_colored_heatmap(labelsdf, title = paste0(lab, " (top neighborhoods)"))
      }

      utils::write.csv(da, file.path(comp_dir, "da_results.csv"), row.names = FALSE)
      utils::write.csv(da_top, file.path(comp_dir, "da_results_top.csv"), row.names = FALSE)
      utils::write.csv(summary_df, file.path(comp_dir, "cluster_summary.csv"), row.names = FALSE)
      utils::write.csv(labelsdf, file.path(comp_dir, "top_nhood_cell_overlap.csv"), row.names = FALSE)
      .milo_ggsave(beeswarm, file.path(comp_dir, "beeswarm.pdf"))
      .milo_ggsave(graph, file.path(comp_dir, "nhood_graph_da.pdf"))
      .milo_ggsave(ma, file.path(comp_dir, "ma_plot.pdf"))
      if (!is.null(alluvial)) {
        .milo_muffle_to_log(
          .milo_ggsave(alluvial, file.path(comp_dir, "alluvial_top.pdf"))
        )
      }
      if (!is.null(twt)) .milo_save_heatmap(twt, file.path(comp_dir, "twt_top.pdf"))

      list(
        skipped = FALSE,
        skip_reason = NA_character_,
        label = lab,
        c0 = as.character(comps$c0[i]),
        c1 = as.character(comps$c1[i]),
        formula = comps$formula[i],
        contrast = parsed$display,
        model_contrast = fit$model_contrast,
        pairing_status = fit$pairing_status,
        pairing_note = fit$pairing_note,
        da = da,
        da_top = da_top,
        cluster_summary = summary_df,
        n_significant = nrow(sig),
        beeswarm = beeswarm,
        graph = graph,
        ma = ma,
        alluvial = alluvial,
        twt = twt,
        outdir = comp_dir
      )
    }, error = function(e) {
      warning("Skipping Milo comparison '", lab, "': ", conditionMessage(e), call. = FALSE)
      list(
        skipped = TRUE,
        skip_reason = conditionMessage(e),
        label = lab,
        c0 = as.character(comps$c0[i]),
        c1 = as.character(comps$c1[i]),
        formula = comps$formula[i],
        contrast = comps$contrast[i],
        model_contrast = NA_character_,
        pairing_status = NA_character_,
        pairing_note = conditionMessage(e),
        da = NULL,
        da_top = NULL,
        cluster_summary = NULL,
        n_significant = NA_integer_,
        beeswarm = NULL,
        graph = NULL,
        ma = NULL,
        alluvial = NULL,
        twt = NULL,
        outdir = comp_dir
      )
    })
    comparisons[[lab]] <- one
  }

  list(
    run_info = run_info,
    size_plot = size_plot,
    nhood_counts = nhood_counts,
    overlap = overlap,
    outdir = outdir,
    comparisons = comparisons,
    k = k_info$k,
    d = built$d_use
  )
}
