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

#' Palette used by RISC::InPlot for up to 35 samples.
#' @keywords internal
.risc_inplot_colors <- function(n) {
  qual <- RColorBrewer::brewer.pal.info[
    rownames(RColorBrewer::brewer.pal.info) %in% c("Set1", "Set3", "Dark2", "Paired"),
    ,
    drop = FALSE
  ]
  cols <- unique(unlist(mapply(RColorBrewer::brewer.pal, qual$maxcolors, rownames(qual))))
  cols <- cols[order(cols, decreasing = FALSE)]
  cols <- cols[-c(37:40)]
  rgb <- grDevices::col2rgb(cols)
  dist_color <- as.matrix(stats::dist(t(rgb)))
  diag(dist_color) <- 1e10
  while (length(cols) > n) {
    min_col <- apply(dist_color, 1, min)
    ids <- which(min_col == min(min_col))[1]
    dist_color <- dist_color[-ids, -ids, drop = FALSE]
    cols <- cols[-ids]
  }
  cols
}

#' Draw the three InPlot panels from one diagnostics run.
#' @keywords internal
.risc_inplot_draw <- function(pc3, pc4, nset) {
  cols <- NULL
  if (nset <= 35L) {
    cols <- .risc_inplot_colors(nset)
  }
  plot_stv <- ggplot2::ggplot(pc3, ggplot2::aes(x = Num, y = PC2, g = Group)) +
    ggplot2::geom_point(ggplot2::aes(fill = Group2), size = 3, shape = 21, alpha = 0.8) +
    ggplot2::geom_smooth(
      ggplot2::aes(color = Group),
      se = FALSE,
      method = "loess",
      formula = "y~x",
      span = 1,
      linetype = "twodash",
      show.legend = FALSE
    ) +
    ggplot2::labs(x = "PC Number", y = "Score (Stv by PCs)", fill = "Scores") +
    ggplot2::theme_bw(base_line_size = 0, base_size = 12) +
    ggplot2::theme(axis.title = ggplot2::element_text(size = 12, face = "bold")) +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 4)))
  plot_ks <- ggplot2::ggplot(pc3, ggplot2::aes(x = Group, y = PC1)) +
    ggplot2::geom_violin(ggplot2::aes(color = Group), show.legend = FALSE) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = Group1), width = 0.1, show.legend = FALSE) +
    ggplot2::geom_jitter(
      ggplot2::aes(fill = Group1),
      size = 3,
      pch = 21,
      alpha = 0.8,
      width = 0.25
    ) +
    ggplot2::geom_label(
      ggplot2::aes(x = Group, y = mean0, label = mean0),
      color = "firebrick",
      nudge_x = 0.25
    ) +
    ggplot2::labs(x = "", y = "Score (Kolmogorov-Smirnov)", fill = "Scores") +
    ggplot2::theme_bw(base_line_size = 0, base_size = 12) +
    ggplot2::theme(axis.title = ggplot2::element_text(size = 12, face = "bold")) +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 4)))
  plot_clust <- ggplot2::ggplot(pc4, ggplot2::aes(x = PCs, y = Cluster, group = Group)) +
    ggplot2::geom_bar(
      ggplot2::aes(fill = Group),
      stat = "identity",
      width = 0.8,
      position = "dodge"
    ) +
    ggplot2::labs(x = "PCs", y = "Cluster Num.", fill = "Scores") +
    ggplot2::theme_bw(base_line_size = 0, base_size = 12) +
    ggplot2::theme(axis.title = ggplot2::element_text(size = 12, face = "bold")) +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 4)))
  if (!is.null(cols)) {
    plot_stv <- plot_stv +
      ggplot2::scale_fill_manual(values = cols) +
      ggplot2::scale_color_manual(values = cols)
    plot_ks <- plot_ks +
      ggplot2::scale_fill_manual(values = cols) +
      ggplot2::scale_color_manual(values = cols)
    plot_clust <- plot_clust + ggplot2::scale_fill_manual(values = cols)
  }
  # Same gtable grid.arrange() returns, without drawing during the sweep.
  gridExtra::arrangeGrob(plot_clust, plot_stv, plot_ks, ncol = 1)
}

#' InPlot calculation used for autoV2: figure, ranking scores, and panel table.
#'
#' One pass over the samples. Scaling, SVD width, PC bins, KS, Stv, and Louvain
#' follow \code{RISC::InPlot}. Each task sets the seed so the result does not
#' depend on worker order. The cluster is stopped before return.
#'
#' @param make_plot Logical; when \code{FALSE}, \code{plot} is \code{NULL}.
#' @return List with \code{plot}, \code{scores}, and \code{panel}.
#' @keywords internal
.risc_inplot_diagnostics <- function(
    risclist,
    var.gene = NULL,
    nPC = 20L,
    minPC = 11L,
    neighbor = 30L,
    res = 1,
    Std.cut = 0.95,
    bin = 5L,
    ncore = 1L,
    algorithm = "kd_tree",
    make_plot = TRUE,
    seed = 123L
) {
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
    stop("Too few shared genes for RISC InPlot diagnostics.", call. = FALSE)
  }
  npc <- as.integer(nPC)[1]
  minpc <- as.integer(minPC)[1]
  neighbor <- as.integer(neighbor)[1]
  res <- as.numeric(res)[1]
  bin <- as.integer(bin)[1]
  ncore <- max(1L, as.integer(ncore)[1])
  seed <- as.integer(seed)[1]
  force(Std.cut)
  if (minpc >= npc) {
    minpc <- npc
  }
  bin0 <- as.integer(seq(minpc, npc, length.out = bin))

  use_par <- ncore > 1L && nset > 1L
  if (use_par) {
    cl <- parallel::makeCluster(ncore)
    doParallel::registerDoParallel(cl)
    on.exit({
      try(parallel::stopCluster(cl), silent = TRUE)
      foreach::registerDoSEQ()
    }, add = TRUE)
  }

  if (use_par) {
    PC0 <- foreach::`%dopar%`(
      foreach::foreach(obj = risclist, .packages = c("irlba", "Matrix")),
      {
        Sys.setenv(
          OMP_NUM_THREADS = "1",
          OPENBLAS_NUM_THREADS = "1",
          MKL_NUM_THREADS = "1",
          VECLIB_MAXIMUM_THREADS = "1"
        )
        set.seed(seed)
        vari <- obj@assay$logcount[var0, ]
        vari <- scale(vari, center = TRUE, scale = TRUE)
        irlba::irlba(vari, nv = npc)
      }
    )
  } else {
    PC0 <- lapply(risclist, function(obj) {
      set.seed(seed)
      vari <- obj@assay$logcount[var0, ]
      vari <- scale(vari, center = TRUE, scale = TRUE)
      irlba::irlba(vari, nv = npc)
    })
  }

  ks_stat <- vector("list", nset)
  ks_mean <- numeric(nset)
  ks_dev <- numeric(nset)
  stv_score <- numeric(nset)
  cumvar <- vector("list", nset)
  for (i in seq_len(nset)) {
    ka <- apply(PC0[[i]]$u, 2, function(x) {
      suppressWarnings(
        stats::ks.test(x, "pnorm", mean = mean(x), sd = stats::sd(x))$statistic
      )
    })
    ka <- as.numeric(ka)
    ks_stat[[i]] <- ka
    ks_mean[i] <- mean(ka)
    ks_dev[i] <- stats::glm(ka ~ seq_len(npc), family = stats::quasipoisson)$deviance
    pvar <- PC0[[i]]$d^2 / sum(PC0[[i]]$d^2)
    pvar <- cumsum(pvar)
    cumvar[[i]] <- as.numeric(pvar)
    stv_coef <- stats::glm(pvar ~ seq_len(npc), family = stats::poisson)$coefficients[2]
    stv_score[i] <- as.numeric(1 - (1 / npc - stv_coef) / (1 / npc))
  }
  ks_dev_scaled <- (1 - ks_dev) * (1 / max(1 - ks_dev))

  v_list <- lapply(PC0, `[[`, "v")
  if (use_par) {
    clust_vec <- foreach::`%dopar%`(
      foreach::foreach(
        i = rep(seq_len(nset), length(bin0)),
        j = rep(bin0, each = nset),
        .packages = c("FNN", "igraph")
      ),
      {
        set.seed(seed)
        v <- v_list[[i]]
        k0 <- FNN::get.knn(v[, 1:j, drop = FALSE], k = neighbor, algorithm = algorithm)
        ki <- data.frame(
          NodStar = rep(seq_len(nrow(v)), neighbor),
          NodEnd = as.vector(k0$nn.index),
          stringsAsFactors = FALSE
        )
        gi <- igraph::graph_from_data_frame(ki, directed = FALSE)
        igraph::E(gi)$weight <- 1 / (1 + as.vector(k0$nn.dist))
        gi <- igraph::simplify(gi)
        clu <- igraph::cluster_louvain(gi, resolution = res)
        as.integer(length(unique(clu$membership)))
      }
    )
  } else {
    clust_vec <- lapply(seq_along(rep(bin0, each = nset)), function(idx) {
      i <- rep(seq_len(nset), length(bin0))[idx]
      j <- rep(bin0, each = nset)[idx]
      set.seed(seed)
      v <- v_list[[i]]
      k0 <- FNN::get.knn(v[, 1:j, drop = FALSE], k = neighbor, algorithm = algorithm)
      ki <- data.frame(
        NodStar = rep(seq_len(nrow(v)), neighbor),
        NodEnd = as.vector(k0$nn.index),
        stringsAsFactors = FALSE
      )
      gi <- igraph::graph_from_data_frame(ki, directed = FALSE)
      igraph::E(gi)$weight <- 1 / (1 + as.vector(k0$nn.dist))
      gi <- igraph::simplify(gi)
      clu <- igraph::cluster_louvain(gi, resolution = res)
      as.integer(length(unique(clu$membership)))
    })
  }
  clust_mat <- matrix(as.integer(unlist(clust_vec)), nrow = nset, ncol = length(bin0))

  cluster_median <- as.numeric(sparseMatrixStats::rowMedians(clust_mat))
  cluster_score <- cluster_median / max(cluster_median)
  scores <- data.frame(
    Code = names(risclist),
    Set = paste0("Set-", seq_len(nset)),
    cluster_median = cluster_median,
    cluster_score = cluster_score,
    stv_score = stv_score,
    ks_mean = ks_mean,
    stringsAsFactors = FALSE
  )

  pc_labels <- paste0("PC-", bin0)
  panel <- data.frame(Code = names(risclist), clust_mat, stringsAsFactors = FALSE)
  colnames(panel)[-1] <- pc_labels
  panel$Stv <- round(stv_score, 2)
  panel$KS <- round(ks_mean, 2)

  plot <- NULL
  if (isTRUE(make_plot)) {
    set_ids <- paste0("Set-", seq_len(nset))
    pc3 <- data.frame(
      PC1 = unlist(ks_stat),
      PC2 = unlist(cumvar),
      Group = rep(set_ids, each = npc),
      Num = rep(seq_len(npc), nset),
      stringsAsFactors = FALSE
    )
    pc3$Group <- factor(pc3$Group, levels = set_ids)
    pc3$mean0 <- rep(round(ks_mean, 2), each = npc)
    pc3$Group1 <- factor(
      pc3$Group,
      labels = paste0(levels(pc3$Group), ": ", round(ks_dev_scaled, 2))
    )
    pc3$Group2 <- factor(
      pc3$Group,
      labels = paste0(levels(pc3$Group), ": ", round(stv_score, 2))
    )
    pc4 <- do.call(rbind, lapply(seq_len(nset), function(i) {
      data.frame(
        Set = set_ids[i],
        PCs = pc_labels,
        Cluster = as.integer(clust_mat[i, ]),
        stringsAsFactors = FALSE
      )
    }))
    cmean <- round(cluster_score, 2)
    pc4$Group <- factor(pc4$Set, levels = set_ids, labels = paste0(set_ids, ": ", cmean))
    pc4$PCs <- factor(pc4$PCs, levels = pc_labels)
    plot <- .risc_inplot_draw(pc3, pc4, nset)
  }

  list(plot = plot, scores = scores, panel = panel)
}

#' InPlot-equivalent per-sample cluster, Stv, and KS scores.
#'
#' Sequential scorer kept for manual reference comparison. autoV2 uses
#' \code{.risc_inplot_diagnostics()} instead.
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
#' @param make_plot Logical. autoV2 stores the diagnostics figure when \code{TRUE}.
#' @return List with \code{selected_code}, \code{would_select}, \code{mode}, \code{table},
#'   \code{refscore} (V1 named vector or NULL), \code{inplot_scores} (or NULL),
#'   and for autoV2 \code{plot} plus \code{panel}.
#' @keywords internal
.risc_resolve_reference <- function(
    risclist,
    sample_metadata,
    risc_reference = NULL,
    ncore = 1L,
    var.gene = NULL,
    make_plot = TRUE
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
  plot <- NULL
  panel <- NULL

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
    diag <- .risc_inplot_diagnostics(
      risclist,
      var.gene = var.gene,
      ncore = ncore,
      make_plot = isTRUE(make_plot)
    )
    inplot_scores <- diag$scores
    pick <- .pick_risc_reference_autov2(inplot_scores)
    selected_code <- pick$selected_code
    would_select <- selected_code
    table <- .risc_reference_table_autov2(inplot_scores, pick, selected_code)
    plot <- diag$plot
    panel <- diag$panel
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
    var.gene_n = length(var.gene),
    plot = plot,
    panel = panel
  )
}

#' Write RISC reference selection artifacts.
#' @keywords internal
.risc_reference_selection_write <- function(
    selection,
    dir,
    saving_algorithm = "RDS",
    nthreads = 1L
) {
  if (is.null(dir) || !nzchar(as.character(dir)[1])) {
    return(invisible(NULL))
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  write_scdapp_object(
    selection,
    file.path(dir, "risc_reference_selection"),
    saving_algorithm = saving_algorithm,
    nthreads = nthreads
  )
  if (!is.null(selection$table)) {
    utils::write.csv(
      selection$table,
      file.path(dir, "risc_reference_selection.csv"),
      row.names = FALSE
    )
  }
  if (!is.null(selection$panel)) {
    utils::write.csv(
      selection$panel,
      file.path(dir, "risc_inplot_scores.csv"),
      row.names = FALSE
    )
  }
  invisible(file.path(dir, "risc_reference_selection.rds"))
}

#' Read RISC reference selection artifacts if present.
#' @keywords internal
.risc_reference_selection_read <- function(dir, saving_algorithm = "RDS") {
  fp <- tryCatch(
    find_scdapp_object(file.path(dir, "risc_reference_selection"), saving_algorithm),
    error = function(e) NA_character_
  )
  if (is.na(fp) || !file.exists(fp)) {
    return(NULL)
  }
  tryCatch(read_scdapp_object(fp), error = function(e) NULL)
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
    ncore = backend$risc_ncore,
    make_plot = TRUE
  )
  if (is.function(backend$set_risc_reference)) {
    backend$set_risc_reference(sel$selected_code)
  }
  .risc_reference_selection_write(
    sel,
    outdir,
    saving_algorithm = if (is.null(backend$saving_algorithm)) "RDS" else backend$saving_algorithm
  )
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
