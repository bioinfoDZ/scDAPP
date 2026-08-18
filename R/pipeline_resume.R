# Pipeline resume: fingerprints, markers, and durable stage caches
# Layout under {outdir}/.scdapp_resume/

.resume_root <- function(outdir) {
  file.path(outdir, ".scdapp_resume")
}

.resume_marker_path <- function(outdir, stage) {
  file.path(.resume_root(outdir), paste0(stage, ".ok"))
}

.resume_fingerprint_path <- function(outdir) {
  file.path(.resume_root(outdir), "fingerprints.rds")
}

.resume_samples_dir <- function(outdir) {
  file.path(.resume_root(outdir), "samples")
}

.resume_qc_report_path <- function(outdir) {
  file.path(.resume_root(outdir), "qc_report_state.rds")
}

.resume_stages <- c("samples", "stability", "integration")

.resume_downstream_stages <- function(stage) {
  idx <- match(stage, .resume_stages)
  if (is.na(idx)) {
    stop("Unknown resume stage: ", stage, call. = FALSE)
  }
  .resume_stages[seq.int(idx, length(.resume_stages))]
}

#' Canonical digest of a named list for fingerprinting.
#' @keywords internal
.resume_digest <- function(x) {
  payload <- .resume_canonicalize(x)
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(payload, algo = "xxhash64"))
  }
  raw <- serialize(payload, connection = NULL, xdr = TRUE)
  # Stable enough fallback without Suggests packages
  sprintf(
    "s%x_l%d",
    sum(as.numeric(raw)) %% (2^31 - 1),
    length(raw)
  )
}

#' @keywords internal
.resume_canonicalize <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.data.frame(x)) {
    x <- as.list(x)
  }
  if (is.list(x) && !is.data.frame(x)) {
    nms <- names(x)
    if (is.null(nms)) {
      return(lapply(x, .resume_canonicalize))
    }
    nms <- sort(nms)
    out <- lapply(nms, function(nm) .resume_canonicalize(x[[nm]]))
    names(out) <- nms
    return(out)
  }
  if (is.factor(x)) {
    return(as.character(x))
  }
  if (is.atomic(x)) {
    if (is.numeric(x)) {
      return(as.numeric(x))
    }
    return(as.character(x))
  }
  as.character(x)
}

#' Sample metadata columns used in fingerprints.
#' @keywords internal
.resume_sample_key <- function(sample_metadata) {
  need <- c("Sample", "Code", "Condition")
  miss <- setdiff(need, colnames(sample_metadata))
  if (length(miss)) {
    stop(
      "sample_metadata missing columns for resume fingerprint: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  df <- sample_metadata[, need, drop = FALSE]
  df[order(as.character(df$Code)), , drop = FALSE]
}

#' @keywords internal
.resume_fingerprint_samples <- function(
    sample_metadata,
    datadir,
    input_seurat_obj,
    min_num_UMI,
    min_num_Feature,
    max_perc_mito,
    max_perc_hemoglobin,
    autofilter_complexity,
    autofilter_mito,
    autofilter_nUMI,
    autofilter_medianabsolutedev_threshold,
    autofilter_loess_negative_residual_threshold,
    doubletFinder,
    cluster_unfiltered,
    pcs_indi,
    res_indi,
    use_labeltransfer,
    refdatapath = NULL,
    m_reference = NULL
) {
  .resume_digest(list(
    stage = "samples",
    samples = .resume_sample_key(sample_metadata),
    datadir = as.character(datadir),
    input_seurat_obj = isTRUE(input_seurat_obj),
    min_num_UMI = min_num_UMI,
    min_num_Feature = min_num_Feature,
    max_perc_mito = max_perc_mito,
    max_perc_hemoglobin = max_perc_hemoglobin,
    autofilter_complexity = autofilter_complexity,
    autofilter_mito = autofilter_mito,
    autofilter_nUMI = autofilter_nUMI,
    autofilter_medianabsolutedev_threshold = autofilter_medianabsolutedev_threshold,
    autofilter_loess_negative_residual_threshold = autofilter_loess_negative_residual_threshold,
    doubletFinder = isTRUE(doubletFinder),
    cluster_unfiltered = isTRUE(cluster_unfiltered),
    pcs_indi = pcs_indi,
    res_indi = res_indi,
    use_labeltransfer = isTRUE(use_labeltransfer),
    refdatapath = if (isTRUE(use_labeltransfer)) as.character(refdatapath) else NA_character_,
    m_reference = if (isTRUE(use_labeltransfer)) as.character(m_reference) else NA_character_
  ))
}

#' @keywords internal
.resume_fingerprint_stability <- function(
    sample_metadata,
    integration_method,
    pcs_int,
    res_int,
    stability_numreps,
    stability_sweep_maxPCs,
    stability_sweep_res,
    stability_propcells.perrep,
    risc_reference,
    RISC_louvain_neighbors,
    input_seurat_obj,
    stability_outdir = NULL
) {
  .resume_digest(list(
    stage = "stability",
    samples = .resume_sample_key(sample_metadata),
    integration_method = as.character(integration_method),
    pcs_int = as.character(pcs_int),
    res_int = as.character(res_int),
    stability_numreps = as.integer(stability_numreps)[1],
    stability_sweep_maxPCs = as.numeric(stability_sweep_maxPCs),
    stability_sweep_res = as.numeric(stability_sweep_res),
    stability_propcells.perrep = as.numeric(stability_propcells.perrep)[1],
    risc_reference = .risc_reference_fingerprint_value(risc_reference),
    RISC_louvain_neighbors = as.integer(RISC_louvain_neighbors)[1],
    input_seurat_obj = isTRUE(input_seurat_obj),
    stability_outdir = if (is.null(stability_outdir)) NA_character_ else as.character(stability_outdir)
  ))
}

#' @keywords internal
.resume_fingerprint_integration <- function(
    sample_metadata,
    integration_method,
    pcs_int,
    res_int,
    risc_reference,
    RISC_louvain_neighbors,
    input_seurat_obj
) {
  .resume_digest(list(
    stage = "integration",
    samples = .resume_sample_key(sample_metadata),
    integration_method = as.character(integration_method),
    pcs_int = as.integer(pcs_int)[1],
    res_int = as.numeric(res_int)[1],
    risc_reference = .risc_reference_fingerprint_value(risc_reference),
    RISC_louvain_neighbors = as.integer(RISC_louvain_neighbors)[1],
    input_seurat_obj = isTRUE(input_seurat_obj)
  ))
}

#' @keywords internal
.resume_read_fingerprints <- function(outdir) {
  fp <- .resume_fingerprint_path(outdir)
  if (!file.exists(fp)) {
    return(list())
  }
  tryCatch(readRDS(fp), error = function(e) list())
}

#' @keywords internal
.resume_write_fingerprints <- function(outdir, fingerprints) {
  dir.create(.resume_root(outdir), recursive = TRUE, showWarnings = FALSE)
  saveRDS(fingerprints, .resume_fingerprint_path(outdir))
  invisible(.resume_fingerprint_path(outdir))
}

#' Clear success markers from a stage through integration.
#' @keywords internal
.resume_invalidate_from <- function(outdir, stage) {
  for (st in .resume_downstream_stages(stage)) {
    mk <- .resume_marker_path(outdir, st)
    if (file.exists(mk)) {
      unlink(mk)
    }
  }
  fps <- .resume_read_fingerprints(outdir)
  for (st in .resume_downstream_stages(stage)) {
    fps[[st]] <- NULL
  }
  .resume_write_fingerprints(outdir, fps)
  invisible(TRUE)
}

#' @keywords internal
.resume_mark_ok <- function(outdir, stage, fingerprint) {
  dir.create(.resume_root(outdir), recursive = TRUE, showWarnings = FALSE)
  mk <- .resume_marker_path(outdir, stage)
  writeLines(fingerprint, mk)
  fps <- .resume_read_fingerprints(outdir)
  fps[[stage]] <- fingerprint
  .resume_write_fingerprints(outdir, fps)
  invisible(TRUE)
}

#' Required stability CSV basenames for a completed sweep.
#' @keywords internal
.resume_stability_csv_names <- function() {
  c(
    "perparam_meanscores.csv",
    "perbootstrap_perparam_scores.csv",
    "perclust_jaccard_mean_acrossbootstraps.csv"
  )
}

#' @keywords internal
.resume_stability_csvs_ok <- function(stability_outdir) {
  if (is.null(stability_outdir) || !dir.exists(stability_outdir)) {
    return(FALSE)
  }
  files <- file.path(stability_outdir, .resume_stability_csv_names())
  all(file.exists(files))
}

#' @keywords internal
.resume_samples_artifacts_ok <- function(outdir, sample_codes) {
  sdir <- .resume_samples_dir(outdir)
  if (!dir.exists(sdir)) {
    return(FALSE)
  }
  if (!file.exists(.resume_qc_report_path(outdir))) {
    return(FALSE)
  }
  needed <- file.path(sdir, paste0(sample_codes, ".rds"))
  all(file.exists(needed))
}

#' @keywords internal
.resume_integration_artifacts_ok <- function(outdir_int, engine = "RISC") {
  objdir <- file.path(outdir_int, "data_objects")
  seu <- file.path(objdir, "Seurat-object_integrated.rds")
  if (!file.exists(seu)) {
    return(FALSE)
  }
  if (identical(engine, "RISC")) {
    risc <- file.path(objdir, "RISC-object_integrated.rds")
    if (!file.exists(risc)) {
      return(FALSE)
    }
  }
  meta <- file.path(objdir, "integration_meta.rds")
  file.exists(meta)
}

#' Whether a stage can be loaded from cache.
#' @keywords internal
.resume_can_load <- function(
    stage,
    outdir,
    fingerprint,
    sample_codes = NULL,
    stability_outdir = NULL,
    outdir_int = NULL,
    engine = "RISC"
) {
  mk <- .resume_marker_path(outdir, stage)
  if (!file.exists(mk)) {
    return(FALSE)
  }
  saved <- tryCatch(readLines(mk, warn = FALSE)[1], error = function(e) NA_character_)
  if (is.na(saved) || !identical(saved, fingerprint)) {
    return(FALSE)
  }
  fps <- .resume_read_fingerprints(outdir)
  if (!identical(fps[[stage]], fingerprint)) {
    return(FALSE)
  }
  if (identical(stage, "samples")) {
    return(.resume_samples_artifacts_ok(outdir, sample_codes))
  }
  if (identical(stage, "stability")) {
    return(.resume_stability_csvs_ok(stability_outdir))
  }
  if (identical(stage, "integration")) {
    return(.resume_integration_artifacts_ok(outdir_int, engine = engine))
  }
  FALSE
}

#' Save per-sample Seurat objects + QC report state; mark samples OK.
#' @keywords internal
.resume_save_samples <- function(outdir, sobjlist, aflist, mdlist, fingerprint) {
  sdir <- .resume_samples_dir(outdir)
  dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
  codes <- names(sobjlist)
  if (is.null(codes)) {
    stop("sobjlist must be a named list of Seurat objects.", call. = FALSE)
  }
  for (code in codes) {
    saveRDS(sobjlist[[code]], file.path(sdir, paste0(code, ".rds")))
  }
  saveRDS(
    list(aflist = aflist, mdlist = mdlist),
    .resume_qc_report_path(outdir)
  )
  .resume_invalidate_from(outdir, "stability")
  .resume_mark_ok(outdir, "samples", fingerprint)
  message("Pipeline resume: saved samples cache under ", .resume_root(outdir))
  invisible(TRUE)
}

#' @keywords internal
.resume_load_samples <- function(outdir, sample_codes) {
  sdir <- .resume_samples_dir(outdir)
  sobjlist <- lapply(sample_codes, function(code) {
    readRDS(file.path(sdir, paste0(code, ".rds")))
  })
  names(sobjlist) <- sample_codes
  qc <- readRDS(.resume_qc_report_path(outdir))
  list(sobjlist = sobjlist, aflist = qc$aflist, mdlist = qc$mdlist)
}

#' Rebuild auto_selection from a completed stability_outdir.
#' @keywords internal
.resume_load_stability_auto_selection <- function(stability_outdir) {
  means <- utils::read.csv(
    file.path(stability_outdir, "perparam_meanscores.csv"),
    stringsAsFactors = FALSE
  )
  best <- .pick_best_stability_params(means)
  selected <- means$params_i[means$MAX == "*"]
  if (!length(selected) || is.na(selected[1])) {
    selected <- means$params_i[
      which.max(0.5 * means$ARI_mean + 0.5 * means$Jaccard_mean_of_clustermeans)
    ]
  }
  sel <- .risc_reference_selection_read(stability_outdir)
  list(
    pcs_int = best$pcs_int,
    res_int = best$res_int,
    selected_params_i = selected[1],
    sweep_dir = stability_outdir,
    stability = list(
      perparam_meanscores = means
    ),
    selected_risc_reference = if (!is.null(sel)) sel$selected_code else NULL,
    risc_reference_selection = sel
  )
}

#' Unlink RawOuts when stability fingerprint mismatches.
#' @keywords internal
.resume_invalidate_stability_rawouts <- function(stability_outdir) {
  raw <- file.path(stability_outdir, "RawOuts")
  if (dir.exists(raw)) {
    unlink(raw, recursive = TRUE)
    message("Pipeline resume: removed stale RawOuts under ", raw)
  }
  invisible(TRUE)
}

#' Copy resume sample RDS into tmpobjdir for update_indi_objects / Seurat int.
#' @keywords internal
.resume_restore_tmpobjdir <- function(outdir, tmpobjdir, sample_codes) {
  sdir <- .resume_samples_dir(outdir)
  dir.create(tmpobjdir, recursive = TRUE, showWarnings = FALSE)
  for (code in sample_codes) {
    src <- file.path(sdir, paste0(code, ".rds"))
    dst <- file.path(tmpobjdir, paste0(code, ".rds"))
    if (!file.exists(src)) {
      stop("Missing resume sample RDS: ", src, call. = FALSE)
    }
    file.copy(src, dst, overwrite = TRUE)
  }
  invisible(tmpobjdir)
}
