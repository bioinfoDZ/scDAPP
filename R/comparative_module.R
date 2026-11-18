# Run composition, DE, pathways, optional ORA, and optional Milo.
# Shared by comparative_module() and subset_module().

.run_comparative_stack <- function(
    obj,
    sample_metadata,
    comps,
    outdir,
    grouping_variable,
    cfg,
    Pseudobulk_mode,
    DE_test = NULL,
    species = "Homo sapiens",
    run_pathways = TRUE,
    run_ora = FALSE,
    run_milo = FALSE,
    milo_d = NULL,
    workernum = 1L,
    verbose = TRUE,
    saving_algorithm = "qs2"
) {
  if (verbose) message("Composition")
  comp_test <- if (isTRUE(Pseudobulk_mode)) "propeller" else "chisq"
  comp_de_test <- if (!is.null(DE_test)) {
    DE_test
  } else if (isTRUE(Pseudobulk_mode)) {
    "EdgeR-LRT"
  } else {
    "wilcox"
  }
  composition <- compositional_analysis_module(
    sobjint = obj,
    comps = comps,
    sample_metadata = sample_metadata,
    outdir_int = outdir,
    grouping_variable = grouping_variable,
    compositional_test = comp_test,
    DE_test = comp_de_test
  )

  if (verbose) message("Differential expression")
  if (isTRUE(Pseudobulk_mode)) {
    de_assay <- cfg$de_pseudobulk_assay
    de_slot <- cfg$de_pseudobulk_slot
  } else {
    de_assay <- cfg$de_wilcox_assay
    de_slot <- cfg$de_wilcox_slot
  }
  de_args <- list(
    sobjint = obj,
    sample_metadata = sample_metadata,
    comps = comps,
    grouping_variable = grouping_variable,
    Pseudobulk_mode = Pseudobulk_mode,
    outdir_int = outdir,
    assay = de_assay,
    slot = de_slot,
    cluster_prefix = FALSE,
    workernum = workernum
  )
  if (!is.null(DE_test)) de_args$DE_test <- DE_test
  de_results <- do.call(de_across_conditions_module, de_args)

  pathways <- NULL
  if (isTRUE(run_pathways)) {
    if (is.null(species) || !nzchar(as.character(species)[1])) {
      stop("species is required when run_pathways = TRUE.", call. = FALSE)
    }
    if (!is.data.frame(de_results) || !nrow(de_results)) {
      warning("Skipping pathways because DE returned no genes.", call. = FALSE)
    } else {
      if (verbose) message("Pathways")
      pathway_sets <- preppathways_pathwayanalysis_crosscondition_module(
        species = species,
        outdir_int = outdir,
        nthreads = workernum,
        saving_algorithm = saving_algorithm
      )
      pathways <- pathwayanalysis_crosscondition_module(
        de_results = de_results,
        pathways = pathway_sets,
        sample_metadata = sample_metadata,
        comps = comps,
        workernum = workernum,
        outdir_int = outdir,
        saving_algorithm = saving_algorithm,
        nthreads = workernum
      )
    }
  }

  ora <- NULL
  if (isTRUE(run_ora)) {
    if (is.null(pathways) && isTRUE(run_pathways)) {
      warning("Skipping ORA because pathway preparation did not run.", call. = FALSE)
    } else {
      if (verbose) message("ORA")
      if (!exists("pathway_sets", inherits = FALSE)) {
        pathway_sets <- preppathways_pathwayanalysis_crosscondition_module(
          species = species,
          outdir_int = outdir,
          nthreads = workernum,
          saving_algorithm = saving_algorithm
        )
      }
      ora <- ORA_crosscondition_module(
        de_results = de_results,
        pathways = pathway_sets,
        sample_metadata = sample_metadata,
        comps = comps,
        workernum = workernum,
        outdir_int = outdir
      )
    }
  }

  milo <- NULL
  if (isTRUE(run_milo)) {
    if (verbose) message("Milo")
    milo_args <- list(
      sobjint = obj,
      comps = comps,
      sample_metadata = sample_metadata,
      outdir_int = outdir,
      grouping_variable = grouping_variable,
      cluster_reduction = cfg$cluster_reduction,
      umap_reduction = cfg$umap_reduction,
      DE_test = comp_de_test
    )
    if (!is.null(milo_d) && !is.na(suppressWarnings(as.integer(milo_d)[1]))) {
      milo_args$d <- as.integer(milo_d)[1]
    }
    milo <- do.call(milo_neighbourhood_da_module, milo_args)
  }

  list(
    composition = composition,
    de_results = de_results,
    pathways = pathways,
    ora = ora,
    milo = milo
  )
}

#' Compare conditions within an existing grouping
#'
#' Run composition, differential expression, pathway GSEA, and optional Milo on
#' an integrated Seurat object, without subsetting or reclustering. Cells are
#' stratified by \code{grouping_variable}, for example a metadata column of
#' called cell types.
#'
#' Composition uses propeller when \code{Pseudobulk_mode} is \code{TRUE} and a
#' chi-square test otherwise. Differential expression uses \code{DE_test}.
#' Pathway GSEA runs by default. Milo runs when pseudobulk mode is on. ORA runs
#' only when \code{run_ora = TRUE}. An HTML report is written to \code{outdir}
#' unless \code{render_html = FALSE}.
#'
#' Assay, slot, and Milo reductions follow the parent integration method, read
#' from \code{sobjint@misc$scDAPP$integration_method} unless
#' \code{integration_method} is passed. Milo uses every column of the clustering
#' reduction.
#'
#' @param sobjint Integrated Seurat object.
#' @param sample_metadata Sample table. \code{Code} and \code{Condition} are
#'   used when present. Otherwise \code{sample_column} and \code{condition_column}
#'   are copied onto those names.
#' @param comps Comparison table with \code{c0} and \code{c1}. See \code{normalize_comps()}.
#' @param outdir Output directory for the comparative run and HTML report.
#' @param grouping_variable Metadata column that defines the groups, such as
#'   called cell types. Results are computed within each level of this column.
#' @param sample_column Sample id column on the object. Default \code{"Code"}.
#' @param condition_column Condition column on the object. Default \code{"Condition"}.
#' @param integration_method Parent integration method. Default reads
#'   \code{sobjint@misc$scDAPP$integration_method}.
#' @param Pseudobulk_mode If \code{TRUE}, composition uses propeller, DE uses a
#'   pseudobulk test on \code{RNA}/\code{counts}, and Milo runs. If \code{FALSE},
#'   composition is a chi-square test and DE is a single-cell test on the
#'   integration method's normalized assay.
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
#' @param workernum Workers for DE and pathways.
#' @param render_html Write \code{comparative_module_report.html} in \code{outdir}.
#' @param verbose Print progress.
#' @param saving_algorithm \code{"qs2"} or \code{"rds"}.
#' @return List with the object used, \code{grouping_variable}, comparative
#'   module results, the HTML path, and \code{outdir}.
#' @export
comparative_module <- function(
    sobjint,
    sample_metadata,
    comps,
    outdir,
    grouping_variable,
    sample_column = "Code",
    condition_column = "Condition",
    integration_method = NULL,
    Pseudobulk_mode = TRUE,
    species = "Homo sapiens",
    run_pathways = TRUE,
    run_ora = FALSE,
    run_milo = NULL,
    DE_test = NULL,
    workernum = 1L,
    render_html = TRUE,
    verbose = TRUE,
    saving_algorithm = "qs2"
) {
  if (missing(sobjint)) stop("Provide an integrated Seurat object.", call. = FALSE)
  if (missing(outdir) || !nzchar(as.character(outdir)[1])) {
    stop("Provide outdir.", call. = FALSE)
  }
  if (missing(grouping_variable)) {
    stop(
      "Provide grouping_variable, a metadata column such as called cell types.",
      call. = FALSE
    )
  }
  grouping_variable <- as.character(grouping_variable)[1]
  if (!nzchar(grouping_variable)) {
    stop(
      "Provide grouping_variable, a metadata column such as called cell types.",
      call. = FALSE
    )
  }

  saving_algorithm <- normalize_saving_algorithm(saving_algorithm)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  workernum <- as.integer(workernum)[1]
  if (is.na(workernum) || workernum < 1L) workernum <- 1L

  method <- .subset_parent_method(sobjint, integration_method)
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
  if (!grouping_variable %in% colnames(obj@meta.data)) {
    stop("grouping_variable not in metadata: ", grouping_variable, call. = FALSE)
  }
  filtered <- .subset_filter_metadata(obj, sample_metadata)
  obj <- filtered$obj
  sample_metadata <- filtered$sample_metadata
  comps <- .subset_check_comps(sample_metadata, comps)
  cfg <- resolve_integration_config(method)
  if (verbose) {
    message(
      "Comparative analysis: ", ncol(obj), " cells, ",
      length(unique(as.character(obj$Code))), " samples, grouped by ",
      grouping_variable, ", method ", method
    )
  }

  comparative <- .run_comparative_stack(
    obj = obj,
    sample_metadata = sample_metadata,
    comps = comps,
    outdir = outdir,
    grouping_variable = grouping_variable,
    cfg = cfg,
    Pseudobulk_mode = Pseudobulk_mode,
    DE_test = DE_test,
    species = species,
    run_pathways = run_pathways,
    run_ora = run_ora,
    run_milo = run_milo,
    milo_d = NULL,
    workernum = workernum,
    verbose = verbose,
    saving_algorithm = saving_algorithm
  )

  dimplots <- tryCatch(
    .subset_dimplots(obj, grouping_variable, cfg$umap_reduction),
    error = function(e) {
      warning("DimPlots failed: ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )

  html_path <- NULL
  if (isTRUE(render_html)) {
    if (verbose) message("Rendering comparative HTML report")
    n_groups <- length(unique(as.character(obj@meta.data[[grouping_variable]])))
    report <- list(
      method = method,
      grouping_variable = grouping_variable,
      n_cells = ncol(obj),
      n_samples = nrow(sample_metadata),
      n_groups = n_groups,
      pseudobulk = isTRUE(Pseudobulk_mode),
      dimplots = dimplots,
      composition = comparative$composition,
      de_top = .subset_de_top_table(comparative$de_results, 10L),
      pathways = comparative$pathways,
      ora = comparative$ora,
      milo = comparative$milo
    )
    report_rds <- file.path(outdir, "comparative_module_report_data.rds")
    saveRDS(report, report_rds)
    rmd <- system.file("rmd", "comparative_module_report.Rmd", package = "scDAPP")
    if (!nzchar(rmd)) {
      stop("Could not find inst/rmd/comparative_module_report.Rmd.", call. = FALSE)
    }
    html_path <- file.path(outdir, "comparative_module_report.html")
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
    grouping_variable = grouping_variable,
    composition = comparative$composition,
    de_results = comparative$de_results,
    pathways = comparative$pathways,
    ora = comparative$ora,
    milo = comparative$milo,
    html = html_path,
    outdir = outdir
  )
}
