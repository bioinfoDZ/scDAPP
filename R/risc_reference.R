# ---- RISC reference selection (autoV2 / auto / manual) ----

#' Parse \code{risc_reference} into a method spec.
#'
#' @param risc_reference \code{NULL}, \code{"autoV2"}, \code{"auto"} / \code{"autoV1"},
#'   or a sample Code/Sample name.
#' @return List with \code{mode} (\code{"autoV2"}, \code{"auto"}, \code{"manual"}),
#'   \code{value} (manual name or NULL), and \code{raw} (fingerprint string).
#' @keywords internal
.parse_risc_reference_arg <- function(risc_reference = NULL) {
  if (is.null(risc_reference) || (length(risc_reference) == 1L && is.na(risc_reference))) {
    return(list(mode = "autoV2", value = NULL, raw = "autoV2"))
  }
  x <- trimws(as.character(risc_reference)[1])
  if (!nzchar(x)) {
    return(list(mode = "autoV2", value = NULL, raw = "autoV2"))
  }
  xl <- tolower(x)
  if (xl %in% c("autov2", "auto_v2", "auto-v2")) {
    return(list(mode = "autoV2", value = NULL, raw = "autoV2"))
  }
  if (xl %in% c("auto", "autov1", "auto_v1", "auto-v1")) {
    return(list(mode = "auto", value = NULL, raw = "auto"))
  }
  list(mode = "manual", value = x, raw = x)
}

#' Canonical \code{risc_reference} string for resume fingerprints.
#' @keywords internal
.risc_reference_fingerprint_value <- function(risc_reference = NULL) {
  .parse_risc_reference_arg(risc_reference)$raw
}

#' Shared gene symbols across a RISC object list (same set InPlot uses here).
#' @keywords internal
.risc_shared_var_genes <- function(risclist) {
  Reduce(intersect, lapply(risclist, function(x) as.character(x@rowdata$Symbol)))
}

#' Align \code{sample_metadata} rows to \code{names(risclist)}.
#' @keywords internal
.risc_align_sample_metadata <- function(risclist, sample_metadata) {
  codes <- names(risclist)
  if (is.null(codes) || any(!nzchar(codes))) {
    if (length(risclist) != nrow(sample_metadata)) {
      stop("risclist is unnamed and does not match sample_metadata rows.", call. = FALSE)
    }
    names(risclist) <- sample_metadata$Code
    codes <- names(risclist)
  }
  idx <- match(codes, sample_metadata$Code)
  if (anyNA(idx)) {
    stop(
      "RISC list names are not a subset of sample_metadata$Code: ",
      paste(codes[is.na(idx)], collapse = ", "),
      call. = FALSE
    )
  }
  list(risclist = risclist, sample_metadata = sample_metadata[idx, , drop = FALSE])
}

#' Index of a Code/Sample in sample_metadata.
#' @keywords internal
.risc_lookup_reference_index <- function(sample_metadata, risc_reference, required = TRUE) {
  if (is.null(risc_reference) || !nzchar(as.character(risc_reference)[1])) {
    if (isTRUE(required)) {
      stop("risc_reference is empty.", call. = FALSE)
    }
    return(NULL)
  }
  x <- as.character(risc_reference)[1]
  if (any(sample_metadata$Code == x)) {
    return(which(sample_metadata$Code == x)[1])
  }
  if ("Sample" %in% names(sample_metadata) && any(sample_metadata$Sample == x)) {
    return(which(sample_metadata$Sample == x)[1])
  }
  if (isTRUE(required)) {
    stop(
      "risc_reference '", x,
      "' does not match any sample Code or Sample in sample_metadata.",
      call. = FALSE
    )
  }
  NULL
}

#' Index of the frozen reference inside a (possibly subset) risclist.
#' @return 1-based index, or \code{NA_integer_} if the sample is absent.
#' @keywords internal
.risc_reference_index_in_list <- function(risclist, sample_metadata, risc_reference) {
  if (is.null(risc_reference) || !nzchar(as.character(risc_reference)[1])) {
    return(NA_integer_)
  }
  x <- as.character(risc_reference)[1]
  nms <- names(risclist)
  if (is.null(nms)) {
    nms <- sample_metadata$Code[seq_along(risclist)]
  }
  if (x %in% nms) {
    return(which(nms == x)[1])
  }
  if ("Sample" %in% names(sample_metadata) && x %in% sample_metadata$Sample) {
    code <- sample_metadata$Code[match(x, sample_metadata$Sample)][1]
    if (!is.na(code) && code %in% nms) {
      return(which(nms == code)[1])
    }
  }
  NA_integer_
}

#' Put the reference RISC object first (RISC \code{scMultiIntegrate} convention).
#' @keywords internal
.risc_put_reference_first <- function(risclist, ref) {
  ref <- as.integer(ref)[1]
  if (ref == 1L) {
    return(risclist)
  }
  data0 <- list(risclist[[ref]])
  names(data0) <- names(risclist)[ref]
  for (i in seq_along(risclist)) {
    if (i != ref) {
      data0[[names(risclist)[i]]] <- risclist[[i]]
    }
  }
  data0
}

#' Legacy (auto / autoV1) reference score: size-weighted cluster count times pb variance.
#' @keywords internal
.risc_reference_score_v1 <- function(risclist, sample_metadata) {
  aligned <- .risc_align_sample_metadata(risclist, sample_metadata)
  risclist <- aligned$risclist
  sample_metadata <- aligned$sample_metadata
  numclusts <- vapply(
    risclist,
    function(dat0) length(unique(dat0@coldata$seurat_clusters)),
    integer(1)
  )
  numcells_per_sample <- vapply(risclist, function(dat0) nrow(dat0@coldata), numeric(1))
  numcells_per_sample <- numcells_per_sample / max(numcells_per_sample)
  numclusts <- numclusts * numcells_per_sample
  pbvar <- vapply(risclist, function(dat0) {
    mat <- dat0@assay$logcount
    md <- dat0@coldata
    pb <- suppressMessages(
      scDAPP::pseudobulk(
        obj = mat,
        metadata = md,
        grouping_colname_in_md = "seurat_clusters"
      )
    )
    numcells <- table(md$seurat_clusters)
    pb <- sweep(pb, 2, numcells, FUN = "/")
    clustervar <- apply(pb, 2, stats::var)
    mean(clustervar)
  }, numeric(1))
  refscore <- numclusts * pbvar
  names(refscore) <- sample_metadata$Code
  refscore
}

#' InPlot-equivalent per-sample cluster, Stv, and KS scores.
#'
#' Mirrors \code{RISC::InPlot} (irlba on scaled logcounts of shared genes; Louvain
#' cluster counts across PC bins; Stv slope; mean KS of gene loadings vs Normal).
#'
#' @keywords internal
.risc_inplot_reference_scores <- function(
    risclist,
    var.gene = NULL,
    nPC = 20L,
    minPC = 11L,
    neighbor = 30L,
    res = 1,
    Std.cut = 0.95,
    bin = 5L,
    ncore = 1L,
    algorithm = "kd_tree"
) {
  set.seed(123)
  if (is.null(names(risclist)) || any(!nzchar(names(risclist)))) {
    stop("risclist must be named by sample Code.", call. = FALSE)
  }
  nset <- length(risclist)
  if (nset < 1L) {
    stop("risclist is empty.", call. = FALSE)
  }
  gene0 <- .risc_shared_var_genes(risclist)
  if (is.null(var.gene)) {
    var0 <- gene0
  } else {
    var0 <- intersect(gene0, as.character(var.gene))
  }
  if (length(var0) < 3L) {
    stop("Too few shared genes for RISC InPlot-equivalent scoring.", call. = FALSE)
  }
  npc <- as.integer(nPC)[1]
  minpc <- as.integer(minPC)[1]
  neighbor <- as.integer(neighbor)[1]
  res <- as.numeric(res)[1]
  bin <- as.integer(bin)[1]
  ncore <- max(1L, as.integer(ncore)[1])
  # Std.cut is InPlot's cumulative-variance cutoff for the figure; ranking uses the Stv slope.
  force(Std.cut)
  n_cells <- vapply(risclist, function(x) ncol(x@assay$logcount), integer(1))
  npc <- min(npc, length(var0) - 1L, max(2L, min(n_cells) - 1L))
  if (minpc >= npc) {
    minpc <- npc
  }
  bin0 <- as.integer(seq(minpc, npc, length.out = min(bin, npc - minpc + 1L)))
  bin0 <- unique(bin0[bin0 >= 1L & bin0 <= npc])
  if (!length(bin0)) {
    bin0 <- npc
  }

  pca_one <- function(obj) {
    vari <- obj@assay$logcount[var0, , drop = FALSE]
    vari <- scale(as.matrix(vari), center = TRUE, scale = TRUE)
    vari[!is.finite(vari)] <- 0
    irlba::irlba(vari, nv = npc)
  }
  # One-time full-data scoring; sequential to avoid foreach export issues.
  invisible(ncore)
  PC0 <- lapply(risclist, pca_one)

  ks_mean <- numeric(nset)
  stv_score <- numeric(nset)
  for (i in seq_len(nset)) {
    ka <- apply(PC0[[i]]$u, 2, function(x) {
      suppressWarnings(
        stats::ks.test(x, "pnorm", mean = mean(x), sd = stats::sd(x))$statistic
      )
    })
    ks_mean[i] <- mean(ka)
    pvar <- PC0[[i]]$d^2 / sum(PC0[[i]]$d^2)
    pvar <- cumsum(pvar)
    stv_coef <- stats::glm(pvar ~ seq_len(npc), family = stats::poisson)$coefficients[2]
    stv_score[i] <- 1 - (1 / npc - stv_coef) / (1 / npc)
  }

  nclust_one <- function(i, j) {
    v <- PC0[[i]]$v
    nn <- min(neighbor, nrow(v) - 1L)
    if (nn < 1L) {
      return(1L)
    }
    k0 <- FNN::get.knn(v[, seq_len(j), drop = FALSE], k = nn, algorithm = algorithm)
    ki <- data.frame(
      NodStar = rep(seq_len(nrow(v)), nn),
      NodEnd = as.vector(k0$nn.index),
      stringsAsFactors = FALSE
    )
    gi <- igraph::graph_from_data_frame(ki, directed = FALSE)
    igraph::E(gi)$weight <- 1 / (1 + as.vector(k0$nn.dist))
    gi <- igraph::simplify(gi)
    clu <- igraph::cluster_louvain(gi, resolution = res)
    as.integer(length(unique(clu$membership)))
  }

  clust_mat <- matrix(NA_integer_, nset, length(bin0))
  for (jj in seq_along(bin0)) {
    for (i in seq_len(nset)) {
      clust_mat[i, jj] <- nclust_one(i, bin0[jj])
    }
  }

  cluster_median <- as.numeric(sparseMatrixStats::rowMedians(clust_mat))
  cluster_score <- cluster_median / max(cluster_median)
  data.frame(
    Code = names(risclist),
    Set = paste0("Set-", seq_len(nset)),
    cluster_median = cluster_median,
    cluster_score = cluster_score,
    stv_score = as.numeric(stv_score),
    ks_mean = as.numeric(ks_mean),
    stringsAsFactors = FALSE
  )
}

#' Rank InPlot scores: cluster > stv > KS, with a KS outlier veto.
#'
#' @param scores data.frame from \code{.risc_inplot_reference_scores()}.
#' @return List with \code{selected_code}, \code{vetoed} (logical), \code{rank}.
#' @keywords internal
.pick_risc_reference_autov2 <- function(scores) {
  n <- nrow(scores)
  if (n < 1L) {
    stop("No samples to rank for RISC autoV2 reference selection.", call. = FALSE)
  }
  veto <- rep(FALSE, n)
  ks <- as.numeric(scores$ks_mean)
  if (n >= 4L) {
    fence <- stats::median(ks) + 2 * stats::mad(ks)
    veto <- ks > fence
    if (all(veto)) {
      veto[] <- FALSE
    }
  } else if (n >= 2L) {
    ord <- order(ks, decreasing = TRUE)
    second <- ks[ord[2]]
    if (is.finite(second) && second > 0 && ks[ord[1]] > 1.5 * second &&
        ks[ord[1]] > ks[ord[2]]) {
      veto[ord[1]] <- TRUE
    }
  }
  keep <- which(!veto)
  if (!length(keep)) {
    keep <- seq_len(n)
    veto[] <- FALSE
  }
  sub <- scores[keep, , drop = FALSE]
  o <- order(
    -as.numeric(sub$cluster_score),
    -as.numeric(sub$stv_score),
    as.numeric(sub$ks_mean),
    as.character(sub$Code)
  )
  selected <- as.character(sub$Code)[o[1]]
  rank_keep <- integer(length(keep))
  rank_keep[o] <- seq_along(o)
  rank_all <- rep(NA_integer_, n)
  rank_all[keep] <- rank_keep
  list(
    selected_code = selected,
    vetoed = veto,
    rank = rank_all
  )
}

#' Build the autoV2 HTML/CSV table.
#' @keywords internal
.risc_reference_table_autov2 <- function(scores, pick, selected_code = NULL) {
  if (is.null(selected_code)) {
    selected_code <- pick$selected_code
  }
  data.frame(
    Set = scores$Set,
    Code = scores$Code,
    cluster = round(as.numeric(scores$cluster_score), 4),
    cluster_median = round(as.numeric(scores$cluster_median), 2),
    stv = round(as.numeric(scores$stv_score), 4),
    KS = round(as.numeric(scores$ks_mean), 4),
    Vetoed = ifelse(pick$vetoed, "yes", ""),
    Rank = pick$rank,
    Selected = ifelse(as.character(scores$Code) == as.character(selected_code), "*", ""),
    stringsAsFactors = FALSE
  )
}

#' Build the legacy auto HTML/CSV table.
#' @keywords internal
.risc_reference_table_v1 <- function(refscore, selected_code) {
  data.frame(
    Code = names(refscore),
    Set = paste0("Set-", seq_along(refscore)),
    RefScore = as.numeric(refscore),
    Max = ifelse(names(refscore) == selected_code, "*", ""),
    stringsAsFactors = FALSE
  )
}

#' Select or validate the RISC reference sample.
#'
#' @param ncore Workers for InPlot-equivalent PCA/Louvain (autoV2 / manual comparison).
#' @return List with \code{selected_code}, \code{would_select}, \code{mode}, \code{table},
#'   \code{refscore} (V1 named vector or NULL), \code{inplot_scores} (or NULL).
#' @keywords internal
.risc_resolve_reference <- function(
    risclist,
    sample_metadata,
    risc_reference = NULL,
    ncore = 1L,
    var.gene = NULL
) {
  spec <- .parse_risc_reference_arg(risc_reference)
  aligned <- .risc_align_sample_metadata(risclist, sample_metadata)
  risclist <- aligned$risclist
  sm <- aligned$sample_metadata
  if (is.null(var.gene)) {
    var.gene <- .risc_shared_var_genes(risclist)
  }
  ncore <- max(1L, as.integer(ncore)[1])

  table <- NULL
  refscore <- NULL
  inplot_scores <- NULL
  pick <- NULL

  if (identical(spec$mode, "manual")) {
    idx <- .risc_lookup_reference_index(sm, spec$value, required = TRUE)
    selected_code <- as.character(sm$Code)[idx]
    inplot_scores <- .risc_inplot_reference_scores(
      risclist,
      var.gene = var.gene,
      ncore = ncore
    )
    pick <- .pick_risc_reference_autov2(inplot_scores)
    would_select <- pick$selected_code
    table <- .risc_reference_table_autov2(inplot_scores, pick, selected_code)
  } else if (identical(spec$mode, "auto")) {
    refscore <- .risc_reference_score_v1(risclist, sm)
    selected_code <- names(refscore)[which.max(refscore)]
    would_select <- selected_code
    table <- .risc_reference_table_v1(refscore, selected_code)
  } else {
    inplot_scores <- .risc_inplot_reference_scores(
      risclist,
      var.gene = var.gene,
      ncore = ncore
    )
    pick <- .pick_risc_reference_autov2(inplot_scores)
    selected_code <- pick$selected_code
    would_select <- selected_code
    table <- .risc_reference_table_autov2(inplot_scores, pick, selected_code)
  }

  vetoed_codes <- character(0)
  if (!is.null(pick) && any(pick$vetoed)) {
    vetoed_codes <- as.character(inplot_scores$Code)[pick$vetoed]
  }

  list(
    selected_code = selected_code,
    would_select = would_select,
    mode = spec$mode,
    user_arg = spec$raw,
    table = table,
    refscore = refscore,
    inplot_scores = inplot_scores,
    vetoed_codes = vetoed_codes,
    var.gene_n = length(var.gene)
  )
}

#' Write RISC reference selection artifacts.
#' @keywords internal
.risc_reference_selection_write <- function(selection, dir) {
  if (is.null(dir) || !nzchar(as.character(dir)[1])) {
    return(invisible(NULL))
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(selection, file.path(dir, "risc_reference_selection.rds"))
  if (!is.null(selection$table)) {
    utils::write.csv(
      selection$table,
      file.path(dir, "risc_reference_selection.csv"),
      row.names = FALSE
    )
  }
  invisible(file.path(dir, "risc_reference_selection.rds"))
}

#' Read RISC reference selection artifacts if present.
#' @keywords internal
.risc_reference_selection_read <- function(dir) {
  fp <- file.path(dir, "risc_reference_selection.rds")
  if (!file.exists(fp)) {
    return(NULL)
  }
  tryCatch(readRDS(fp), error = function(e) NULL)
}

#' Freeze RISC reference on a stability backend after \code{prep_reference()}.
#' @keywords internal
.stability_risc_freeze_reference <- function(backend, outdir, verbose = TRUE) {
  if (!identical(backend$engine, "RISC")) {
    return(NULL)
  }
  if (is.null(backend$get_risclist())) {
    backend$prep_reference()
  }
  sel <- .risc_resolve_reference(
    backend$get_risclist(),
    backend$sample_metadata,
    risc_reference = backend$risc_reference,
    ncore = backend$risc_ncore
  )
  if (is.function(backend$set_risc_reference)) {
    backend$set_risc_reference(sel$selected_code)
  }
  .risc_reference_selection_write(sel, outdir)
  if (isTRUE(verbose)) {
    message(
      "RISC reference (", sel$mode, "): ", sel$selected_code,
      if (!identical(sel$would_select, sel$selected_code)) {
        paste0(" (autoV2 would select ", sel$would_select, ")")
      } else {
        ""
      }
    )
    if (length(sel$vetoed_codes)) {
      message(
        "RISC autoV2 KS veto: ",
        paste(sel$vetoed_codes, collapse = ", ")
      )
    }
  }
  sel
}
