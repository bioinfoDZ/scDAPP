# https://r-pkgs.org/whole-game.html

# =============================================================================
# Internal helpers (not exported)
# Used for comps parsing, design matrices, contrasts, and result formatting.
# =============================================================================

# Standardize comps columns (c0/c1, formula, contrast, label) and backward-compat c2 -> c0.
.normalize_comps <- function(comps) {
  comps <- as.data.frame(comps, stringsAsFactors = FALSE)
  cn <- colnames(comps)
  if ("c2" %in% cn && !"c0" %in% cn) {
    warning("comps column 'c2' is deprecated; renamed to 'c0' (reference level).", call. = FALSE)
    colnames(comps)[cn == "c2"] <- "c0"
  }
  if (!all(c("c0", "c1") %in% colnames(comps))) {
    stop("comps must include columns c0 (reference) and c1 (test).", call. = FALSE)
  }
  if (!"formula" %in% colnames(comps)) comps$formula <- NA_character_
  if (!"contrast" %in% colnames(comps)) comps$contrast <- NA_character_
  if (!"label" %in% colnames(comps)) comps$label <- NA_character_
  for (i in seq_len(nrow(comps))) {
    if (is.na(comps$formula[i]) || !nzchar(trimws(comps$formula[i]))) {
      comps$formula[i] <- "~ Condition"
    } else {
      f <- trimws(as.character(comps$formula[i]))
      if (!startsWith(f, "~")) f <- paste("~", f)
      comps$formula[i] <- f
    }
    if (is.na(comps$label[i]) || !nzchar(trimws(comps$label[i]))) {
      comps$label[i] <- paste0(comps$c1[i], "_vs_", comps$c0[i])
    }
  }
  comps$labels <- comps$label
  comps
}

# Parse a per-row formula string; default is ~ Condition.
.parse_comp_formula <- function(formula_entry) {
  if (is.null(formula_entry) || length(formula_entry) == 0 ||
      (length(formula_entry) == 1 && (is.na(formula_entry) || !nzchar(trimws(as.character(formula_entry)))))) {
    return(as.formula("~ Condition"))
  }
  f <- trimws(as.character(formula_entry))
  if (!startsWith(f, "~")) f <- paste("~", f)
  as.formula(f)
}

# Character form of a formula for regex parsing of (1|var) terms.
.formula_chr <- function(formula) {
  paste(deparse(as.formula(formula), width.cutoff = 500L), collapse = " ")
}

# Extract intercept random-effect variables from (1|Var) terms.
.formula_random_effect_vars <- function(formula) {
  f_chr <- .formula_chr(formula)
  m <- gregexpr("\\(\\s*1\\s*\\|\\s*([A-Za-z.][A-Za-z0-9_.]*)\\s*\\)", f_chr, perl = TRUE)
  starts <- as.integer(m[[1]])
  if (length(starts) == 1L && starts[1] == -1L) return(character())
  lens <- attr(m[[1]], "match.length")
  caps <- attr(m[[1]], "capture.start")
  cap_lens <- attr(m[[1]], "capture.length")
  vapply(seq_along(starts), function(i) {
    substr(f_chr, caps[i], caps[i] + cap_lens[i] - 1L)
  }, character(1))
}

.formula_has_random_effect <- function(formula) {
  length(.formula_random_effect_vars(formula)) > 0L
}

# Strip (1|var) terms so model.matrix / propeller fixed designs stay valid.
.fixed_effects_formula <- function(formula) {
  f_chr <- .formula_chr(formula)
  f_chr <- gsub("\\(\\s*1\\s*\\|\\s*[A-Za-z.][A-Za-z0-9_.]*\\s*\\)", "", f_chr)
  f_chr <- gsub("\\+\\s*\\+", "+", f_chr)
  f_chr <- gsub("~\\s*\\+", "~ ", f_chr)
  f_chr <- gsub("\\+\\s*$", "", f_chr)
  f_chr <- trimws(gsub("\\s+", " ", f_chr))
  if (!nzchar(f_chr) || grepl("^~\\s*$", f_chr)) f_chr <- "~ 1"
  as.formula(f_chr)
}

# Enforce Dream vs edgeR/DESeq2 pairing rules for a comps formula.
.validate_de_formula <- function(DE_test, formula) {
  has_re <- .formula_has_random_effect(formula)
  if (identical(DE_test, "Dream")) {
    if (!has_re) {
      stop(
        "DE_test = 'Dream' requires a random effect in comps$formula, ",
        "e.g. '~ Condition + (1|Patient)'. Got: ", .formula_chr(formula),
        call. = FALSE
      )
    }
    return(invisible(TRUE))
  }
  if (DE_test %in% c("EdgeR", "EdgeR-LRT", "EdgeR-QLF", "DESeq2", "DESeq2-LRT") && has_re) {
    stop(
      "comps$formula includes a random effect (", .formula_chr(formula), "). ",
      "Use DE_test = 'Dream' for paired mixed-model DE, or remove (1|var) for edgeR/DESeq2.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Scan all comps formulas for Dream / edgeR-DESeq2 consistency (pipeline preflight).
.validate_comps_de_formulas <- function(comps, DE_test) {
  comps <- .normalize_comps(comps)
  forms <- unique(comps$formula)
  has_any_re <- any(vapply(forms, function(f) {
    .formula_has_random_effect(.parse_comp_formula(f))
  }, logical(1)))
  if (identical(DE_test, "Dream")) {
    if (!has_any_re) {
      stop(
        "DE_test = 'Dream' requires at least one comps$row with a random effect, ",
        "e.g. formula = '~ Condition + (1|Patient)'.",
        call. = FALSE
      )
    }
  } else if (DE_test %in% c("EdgeR", "EdgeR-LRT", "EdgeR-QLF", "DESeq2", "DESeq2-LRT") && has_any_re) {
    stop(
      "At least one comps$formula includes (1|var). ",
      "Use DE_test = 'Dream' for paired mixed-model DE, or remove random effects for edgeR/DESeq2.",
      call. = FALSE
    )
  }
  for (f in forms) .validate_de_formula(DE_test, .parse_comp_formula(f))
  invisible(TRUE)
}

# Soft check that a parsed contrast is estimable under a design matrix / resultsNames.
.validate_contrast_for_design <- function(parsed, design_or_names, style,
                                          context = "") {
  cn <- if (is.matrix(design_or_names) || is.data.frame(design_or_names)) {
    colnames(design_or_names)
  } else {
    as.character(design_or_names)
  }
  prefix <- if (nzchar(context)) paste0(context, ": ") else ""
  if (identical(style, "edger_expr")) {
    invisible(.edger_expr_to_vector(
      matrix(0, nrow = 1, ncol = length(cn), dimnames = list(NULL, cn)),
      parsed$edger_expr
    ))
    return(invisible(TRUE))
  }
  if (identical(style, "deseq2")) {
    if (!is.null(parsed$deseq_name)) {
      nm <- parsed$deseq_name
      candidates <- unique(c(
        nm,
        gsub("\\.", ":", nm, fixed = FALSE),
        gsub(":", ".", nm, fixed = TRUE)
      ))
      if (!any(candidates %in% cn)) {
        stop(
          prefix, "contrast name='", nm,
          "' not in design/resultsNames: ", paste(cn, collapse = ", "),
          call. = FALSE
        )
      }
      return(invisible(TRUE))
    }
    trip <- parsed$deseq_triple
    num_col <- paste0(trip[1], trip[2])
    denom_col <- paste0(trip[1], trip[3])
    # Treatment coding often drops the reference level column; presence of either
    # side (or both) is enough for preflight on full metadata. Per-cluster may still fail.
    if (!(num_col %in% cn) && !(denom_col %in% cn) &&
        !(trip[1] %in% cn)) {
      # Also accept factor name alone when contrast uses results() contrast= API
      # which does not require both columns in resultsNames.
      warning(
        prefix, "contrast ", paste(trip, collapse = ";"),
        " may not map cleanly to columns ", paste(cn, collapse = ", "),
        " (often OK for DESeq2 results(contrast=); check per-cluster).",
        call. = FALSE
      )
    }
    return(invisible(TRUE))
  }
  invisible(TRUE)
}

#' Preflight comps formulas/contrasts against sample metadata and DE_test
#'
#' Validates that design variables exist, Condition covers c0/c1, contrasts parse
#' for the active DE_test style, Dream/(1|var) rules hold, and (when possible)
#' the fixed design is full rank on the full metadata. Per-cluster failures are
#' still soft-skipped later in DE. Propeller uses `.preflight_propeller_comps()`
#' on its cell-means design (not this treatment-coded DE design).
#'
#' @param sample_metadata data.frame with Condition/Code and design columns
#' @param comps comps data.frame (will be normalized)
#' @param DE_test character DE backend
#' @param Pseudobulk_mode logical
#' @param check_de_formula_rules logical; if FALSE, skip Dream vs edgeR/DESeq2
#'   exclusivity (compositional propeller may use `(1|var)` with any DE_test).
#' @return invisible TRUE
#' @keywords internal
.preflight_comps_design <- function(sample_metadata, comps, DE_test,
                                    Pseudobulk_mode = TRUE,
                                    check_de_formula_rules = TRUE) {
  comps <- .normalize_comps(comps)
  style <- .contrast_style(DE_test)
  md <- as.data.frame(sample_metadata, stringsAsFactors = FALSE)

  if (!"Condition" %in% colnames(md)) {
    stop("sample_metadata must include a Condition column.", call. = FALSE)
  }
  need_lv <- unique(c(as.character(comps$c0), as.character(comps$c1)))
  miss_lv <- setdiff(need_lv, as.character(md$Condition))
  if (length(miss_lv)) {
    stop(
      "comps c0/c1 levels missing from sample_metadata$Condition: ",
      paste(miss_lv, collapse = ", "),
      call. = FALSE
    )
  }

  if (!isTRUE(Pseudobulk_mode) || identical(style, "ignore")) {
    message(
      "Preflight: Pseudobulk_mode=", Pseudobulk_mode,
      " / DE_test=", DE_test,
      " — formula/contrast not used for single-cell DE (c0/c1 only)."
    )
    return(invisible(TRUE))
  }

  if (isTRUE(check_de_formula_rules)) {
    .validate_comps_de_formulas(comps, DE_test)
  }
  comps <- .apply_contrast_defaults(comps, DE_test)

  if (identical(DE_test, "EdgeR") && isTRUE(check_de_formula_rules)) {
    for (i in seq_len(nrow(comps))) {
      f <- .parse_comp_formula(comps$formula[i])
      if (!.is_simple_condition_design(f)) {
        stop(
          "DE_test='EdgeR' (exactTest) requires formula ~ Condition only ",
          "(row ", i, " label=", comps$label[i],
          "). Use EdgeR-LRT or EdgeR-QLF for covariates.",
          call. = FALSE
        )
      }
    }
    n_cond <- length(unique(md$Condition))
    if (n_cond > 2L) {
      warning(
        "DE_test='EdgeR' with >2 Condition levels: pipeline will use glmLRT ",
        "for multi-level designs (exactTest is 2-group only).",
        call. = FALSE
      )
    }
  }

  for (i in seq_len(nrow(comps))) {
    f <- .parse_comp_formula(comps$formula[i])
    fixed <- .fixed_effects_formula(f)
    vars <- all.vars(f)
    miss <- setdiff(vars, colnames(md))
    if (length(miss)) {
      stop(
        "sample_metadata missing design variables for comps row ", i,
        " (", comps$label[i], "): ", paste(miss, collapse = ", "),
        call. = FALSE
      )
    }
    parsed <- .parse_contrast(comps$contrast[i], comps$c0[i], comps$c1[i], style)
    coldata <- .prepare_pseudobulk_coldata(md, f)
    design <- tryCatch(
      model.matrix(fixed, data = coldata),
      error = function(e) {
        stop(
          "Could not build design for comps row ", i, " (", comps$label[i], "): ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
    if (qr(design)$rank != ncol(design)) {
      warning(
        "Design not full rank on full sample_metadata for comps row ", i,
        " (", comps$label[i], ", formula=", comps$formula[i],
        "). Some clusters may be skipped.",
        call. = FALSE
      )
    } else if (identical(style, "edger_expr")) {
      tryCatch(
        .validate_contrast_for_design(
          parsed, design, style,
          context = paste0("comps row ", i, " (", comps$label[i], ")")
        ),
        error = function(e) {
          warning(
            conditionMessage(e),
            " — may still succeed under Dream cell-means coding or fail only in sparse clusters.",
            call. = FALSE
          )
        }
      )
    } else if (identical(style, "deseq2") && !is.null(parsed$deseq_name)) {
      tryCatch(
        .validate_contrast_for_design(
          parsed, design, style,
          context = paste0("comps row ", i, " (", comps$label[i], ")")
        ),
        error = function(e) {
          warning(
            conditionMessage(e),
            " — DESeq2 resultsNames may use '.' instead of ':'; check after fit.",
            call. = FALSE
          )
        }
      )
    }
  }
  invisible(TRUE)
}

# Fail-fast check that each comps row maps onto propeller's cell-means design.
# Does not apply to 2-prop z-test (chisq); that path uses c0/c1 only.
.preflight_propeller_comps <- function(sample_metadata, comps, DE_test) {
  comps <- .apply_contrast_defaults(.normalize_comps(comps), DE_test)
  style <- .contrast_style(DE_test)
  md <- as.data.frame(sample_metadata, stringsAsFactors = FALSE)

  if (!"Condition" %in% colnames(md)) {
    stop("sample_metadata must include a Condition column.", call. = FALSE)
  }
  if (!"Code" %in% colnames(md)) {
    stop("sample_metadata must include a Code column for propeller.", call. = FALSE)
  }
  need_lv <- unique(c(as.character(comps$c0), as.character(comps$c1)))
  miss_lv <- setdiff(need_lv, as.character(md$Condition))
  if (length(miss_lv)) {
    stop(
      "comps c0/c1 levels missing from sample_metadata$Condition: ",
      paste(miss_lv, collapse = ", "),
      call. = FALSE
    )
  }

  for (i in seq_len(nrow(comps))) {
    lab <- comps$label[i]
    f <- .parse_comp_formula(comps$formula[i])
    vars <- all.vars(f)
    miss <- setdiff(vars, colnames(md))
    if (length(miss)) {
      stop(
        "sample_metadata missing design variables for propeller comps row ", i,
        " (", lab, "): ", paste(miss, collapse = ", "),
        call. = FALSE
      )
    }
    re_vars <- .formula_random_effect_vars(f)
    if (length(re_vars) > 1L) {
      stop(
        "Propeller paired mode supports one (1|var) block; got: ",
        paste(re_vars, collapse = ", "),
        ". Use a single random-effect term in comps$formula for comparison '",
        lab, "'.",
        call. = FALSE
      )
    }
    pd <- tryCatch(
      .propeller_design(md, f),
      error = function(e) {
        stop(
          "Could not build propeller design for comps row ", i, " (", lab, "): ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
    if (length(re_vars) == 1L) {
      block_var <- re_vars[1]
      if (!block_var %in% colnames(pd$coldata)) {
        stop(
          "Block variable '", block_var,
          "' is missing from sample_metadata for propeller comparison '", lab, "'.",
          call. = FALSE
        )
      }
      block <- pd$coldata[[block_var]][match(rownames(pd$design), pd$coldata$Code)]
      if (anyNA(block)) {
        stop(
          "Missing values in block variable '", block_var,
          "' for propeller comparison '", lab, "'.",
          call. = FALSE
        )
      }
      code_block <- unique(pd$coldata[, c("Code", block_var), drop = FALSE])
      if (any(duplicated(code_block$Code))) {
        stop(
          "Block variable '", block_var,
          "' is not unique per Code for comparison '", lab, "'.",
          call. = FALSE
        )
      }
    }

    use_simple <- length(re_vars) == 0L && .propeller_use_simple_wrapper(f, md)
    if (identical(style, "ignore")) {
      if (!use_simple) {
        stop(
          "Propeller with covariates, interactions, pairing, or >2 Condition levels ",
          "requires DE_test with EdgeR or DESeq2 contrast grammar. Got DE_test='",
          DE_test, "' for comparison '", lab, "'.",
          call. = FALSE
        )
      }
      next
    }

    map_contr <- function(contrast_entry) {
      tryCatch(
        .propeller_contrast_matrix(
          pd$design, contrast_entry, comps$c1[i], comps$c0[i], style = style
        ),
        error = function(e) e
      )
    }
    usr <- map_contr(comps$contrast[i])
    if (inherits(usr, "error")) {
      stop(
        "Propeller contrast for comps row ", i, " (", lab, "): ",
        conditionMessage(usr),
        call. = FALSE
      )
    }
    if (use_simple) {
      default_c <- .default_contrast(comps$c0[i], comps$c1[i], style)
      def <- map_contr(default_c)
      same_default <- !inherits(def, "error") &&
        isTRUE(all.equal(
          as.numeric(usr), as.numeric(def),
          check.attributes = FALSE
        ))
      if (!same_default) {
        stop(
          "Propeller simple two-group path (formula ~ Condition, two Condition levels) ",
          "ignores comps$contrast and always tests Condition c1 vs c0. ",
          "Comparison '", lab, "' has contrast '",
          as.character(comps$contrast[i])[1],
          "' which is not that default. Expand the formula to include the extra ",
          "terms, or omit contrast to use '", default_c, "'.",
          call. = FALSE
        )
      }
    }
  }
  invisible(TRUE)
}

# Contrast style by DE_test: EdgeR/Dream use design-coef expressions; DESeq2 uses
# Factor;num;denom or name=resultsName; Wilcox ignores contrast.
.contrast_style <- function(DE_test) {
  if (identical(DE_test, "Dream") ||
      DE_test %in% c("EdgeR", "EdgeR-LRT", "EdgeR-QLF")) {
    return("edger_expr")
  }
  if (DE_test %in% c("DESeq2", "DESeq2-LRT")) {
    return("deseq2")
  }
  "ignore"
}

.default_contrast <- function(c0, c1, style, condition_col = "Condition") {
  c0 <- as.character(c0)[1]
  c1 <- as.character(c1)[1]
  if (identical(style, "deseq2")) {
    return(paste(condition_col, c1, c0, sep = ";"))
  }
  if (identical(style, "edger_expr")) {
    return(paste0(condition_col, c1, " - ", condition_col, c0))
  }
  NA_character_
}

.is_blank_contrast <- function(contrast_entry) {
  is.null(contrast_entry) || length(contrast_entry) == 0L ||
    (length(contrast_entry) == 1L &&
       (is.na(contrast_entry) || !nzchar(trimws(as.character(contrast_entry)))))
}

# Parse comps$contrast into a structured list for the active style.
# Returns list(style, display, edger_expr=NULL, deseq_triple=NULL, deseq_name=NULL).
.parse_contrast <- function(contrast_entry, c0, c1, style,
                            condition_col = "Condition") {
  if (identical(style, "ignore")) {
    return(list(
      style = "ignore",
      display = NA_character_,
      edger_expr = NULL,
      deseq_triple = NULL,
      deseq_name = NULL
    ))
  }
  if (.is_blank_contrast(contrast_entry)) {
    contrast_entry <- .default_contrast(c0, c1, style, condition_col)
  }
  raw <- trimws(as.character(contrast_entry)[1])

  if (identical(style, "deseq2")) {
    if (grepl("^name\\s*=", raw, ignore.case = TRUE)) {
      nm <- sub("^name\\s*=\\s*", "", raw, ignore.case = TRUE)
      nm <- trimws(nm)
      if (!nzchar(nm)) {
        stop("DESeq2 contrast 'name=' is empty.", call. = FALSE)
      }
      return(list(
        style = "deseq2",
        display = paste0("name=", nm),
        edger_expr = NULL,
        deseq_triple = NULL,
        deseq_name = nm
      ))
    }
    parts <- if (length(contrast_entry) == 3L) {
      as.character(contrast_entry)
    } else {
      strsplit(raw, ";", fixed = TRUE)[[1]]
    }
    parts <- trimws(parts)
    if (length(parts) != 3L) {
      stop(
        "DESeq2 contrast must be 'Factor;numerator;denominator' or ",
        "'name=ResultsName'. Got: ", raw,
        call. = FALSE
      )
    }
    return(list(
      style = "deseq2",
      display = paste(parts, collapse = ";"),
      edger_expr = NULL,
      deseq_triple = parts,
      deseq_name = NULL
    ))
  }

  # edger_expr style — semicolon triples are DESeq2 grammar; reject here
  if (grepl(";", raw, fixed = TRUE)) {
    stop(
      "EdgeR/Dream contrast must be a design-coefficient expression ",
      "(e.g. 'ConditionKO - ConditionControl' or 'ConditionA:BatchB'), ",
      "not a DESeq2 semicolon triple. Got: ", raw,
      call. = FALSE
    )
  }
  list(
    style = "edger_expr",
    display = raw,
    edger_expr = raw,
    deseq_triple = NULL,
    deseq_name = NULL
  )
}

# Apply style-specific defaults onto comps$contrast (mutates display strings).
.apply_contrast_defaults <- function(comps, DE_test) {
  comps <- .normalize_comps(comps)
  style <- .contrast_style(DE_test)
  if (identical(style, "ignore")) return(comps)
  for (i in seq_len(nrow(comps))) {
    parsed <- .parse_contrast(comps$contrast[i], comps$c0[i], comps$c1[i], style)
    comps$contrast[i] <- parsed$display
  }
  comps
}

# Dream contrast expression string (EdgeR-style).
.dream_contrast_expr <- function(contrast_entry, c1, c0) {
  parsed <- .parse_contrast(contrast_entry, c0, c1, style = "edger_expr")
  parsed$edger_expr
}

# Rewrite comps formula for Dream: cell-means Condition so c1 - c0 contrasts exist.
# e.g. ~ Condition + Batch + (1|Patient) -> ~ 0 + Condition + Batch + (1|Patient)
.dream_fit_formula <- function(formula) {
  fixed <- .fixed_effects_formula(formula)
  re_vars <- .formula_random_effect_vars(formula)
  labels <- attr(terms(fixed), "term.labels")
  if (!"Condition" %in% labels) {
    stop("Dream formula requires Condition as a fixed effect.", call. = FALSE)
  }
  # Keep non-Condition fixed terms (including interactions) as written
  other <- labels[labels != "Condition"]
  rhs <- c("0", "Condition", other)
  f_chr <- paste("~", paste(rhs, collapse = " + "))
  if (length(re_vars)) {
    f_chr <- paste(f_chr, "+", paste0("(1|", re_vars, ")", collapse = " + "))
  }
  as.formula(f_chr)
}

.dream_bpparam <- function(workernum) {
  workernum <- as.integer(workernum)[1]
  if (is.na(workernum) || workernum <= 1L) {
    return(BiocParallel::SerialParam())
  }
  BiocParallel::SnowParam(workers = workernum)
}

# Backward-compatible: always returns DESeq2-style character triple.
.parse_comp_contrast <- function(contrast_entry, c1, c0, condition_col = "Condition") {
  if (.is_blank_contrast(contrast_entry)) {
    return(c(condition_col, as.character(c1), as.character(c0)))
  }
  if (length(contrast_entry) == 3L) return(as.character(contrast_entry))
  raw <- trimws(as.character(contrast_entry)[1])
  if (grepl("^name\\s*=", raw, ignore.case = TRUE)) {
    stop(
      "Contrast uses DESeq2 name= form; call .parse_contrast(..., style='deseq2') instead.",
      call. = FALSE
    )
  }
  parts <- strsplit(raw, ";", fixed = TRUE)[[1]]
  if (length(parts) == 3L) return(trimws(parts))
  # EdgeR expr: fall back to Condition;c1;c0 for callers that need a triple
  if (grepl("\\s*-\\s*", raw)) {
    return(c(condition_col, as.character(c1), as.character(c0)))
  }
  stop(
    "contrast must be semicolon-separated Factor;num;denom for this helper. Got: ",
    raw,
    call. = FALSE
  )
}

# Build sample-level colData for edgeR/DESeq2; coerce design variables to factors.
.prepare_pseudobulk_coldata <- function(sample_metadata, design_formula) {
  md <- as.data.frame(sample_metadata, stringsAsFactors = FALSE)
  vars <- all.vars(design_formula)
  miss <- setdiff(vars, colnames(md))
  if (length(miss)) {
    stop("sample_metadata missing design variables: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  for (v in vars) {
    if (!is.numeric(md[[v]])) md[[v]] <- factor(md[[v]])
  }
  md
}

# TRUE when fixed design is only ~ Condition (required for edgeR exactTest / simple propeller).
.is_simple_condition_design <- function(design_formula) {
  fixed <- .fixed_effects_formula(design_formula)
  tl <- attr(terms(fixed), "term.labels")
  length(tl) == 1L && tl[[1]] == "Condition"
}

# DESeq2-LRT reduced model: drop Condition terms, keep covariates (e.g. ~ Batch).
# Kept for backward-compatible callers; prefer .deseq2_reduced_for_contrast().
.make_reduced_formula <- function(design_formula) {
  fixed <- .fixed_effects_formula(design_formula)
  labels <- attr(terms(fixed), "term.labels")
  keep <- labels[!grepl("Condition", labels, fixed = TRUE)]
  if (!length(keep)) return(as.formula("~ 1"))
  as.formula(paste("~", paste(keep, collapse = " + ")))
}

# Nested reduced formula for DESeq2-LRT from a parsed DESeq2 contrast.
# Triple Factor;num;denom -> drop that factor (and interactions containing it).
# name= with '.' (interaction resultsName) -> drop all interaction (:) terms.
.deseq2_reduced_for_contrast <- function(full_formula, parsed) {
  fixed <- .fixed_effects_formula(full_formula)
  labels <- attr(terms(fixed), "term.labels")
  if (!length(labels)) return(as.formula("~ 1"))

  if (!is.null(parsed$deseq_name)) {
    # Interaction coefficient -> additive reduced model
    if (grepl("\\.", parsed$deseq_name, fixed = FALSE) ||
        grepl(":", parsed$deseq_name, fixed = TRUE)) {
      keep <- labels[!grepl(":", labels, fixed = TRUE)]
      if (!length(keep)) return(as.formula("~ 1"))
      return(as.formula(paste("~", paste(keep, collapse = " + "))))
    }
    # Single main-effect resultsName: drop the factor prefix before first level
    # Best-effort: drop terms that are pure interactions only if name has no '.'
    # For non-interaction name=, fall through to dropping Condition if present.
    fac <- sub("^([^A-Z]*[A-Za-z.][A-Za-z0-9_.]*).*", "\\1", parsed$deseq_name)
    # Prefer matching a known term label that prefixes the name
    hit <- labels[vapply(labels, function(lb) {
      !grepl(":", lb, fixed = TRUE) && startsWith(parsed$deseq_name, lb)
    }, logical(1))]
    if (length(hit)) {
      fac <- hit[[1]]
      keep <- labels[labels != fac & !grepl(paste0("(^|:)", fac, "(:|$)"), labels)]
      if (!length(keep)) return(as.formula("~ 1"))
      return(as.formula(paste("~", paste(keep, collapse = " + "))))
    }
  }

  if (!is.null(parsed$deseq_triple)) {
    fac <- as.character(parsed$deseq_triple[[1]])
    # Drop main effect of fac and any interaction terms involving fac
    keep <- labels[
      labels != fac &
        !grepl(paste0("(^|:)", fac, "(:|$)"), labels)
    ]
    if (!length(keep)) return(as.formula("~ 1"))
    return(as.formula(paste("~", paste(keep, collapse = " + "))))
  }

  stop(
    "Cannot build DESeq2-LRT reduced formula from contrast.",
    call. = FALSE
  )
}

# Shallow rebuild so nbinomLRT / DESeq(LRT) does not mutate the dispersion-fitted base.
.copy_deseq_dds_for_lrt <- function(dds) {
  out <- DESeq2::DESeqDataSetFromMatrix(
    countData = DESeq2::counts(dds),
    colData = as.data.frame(SummarizedExperiment::colData(dds)),
    design = DESeq2::design(dds)
  )
  sf <- DESeq2::sizeFactors(dds)
  if (!is.null(sf) && !anyNA(sf)) DESeq2::sizeFactors(out) <- sf
  disp <- DESeq2::dispersions(dds)
  if (!is.null(disp) && !anyNA(disp)) {
    DESeq2::dispersions(out) <- disp
  }
  out
}

# Run / cache nested LRT for one reduced formula on a cluster's DESeq2 base fit.
.deseq2_lrt_fit_cached <- function(fitobj, reduced_f) {
  key <- paste(deparse(reduced_f, width.cutoff = 500L), collapse = " ")
  if (is.null(fitobj$lrt_cache)) {
    stop("DESeq2-LRT fit object missing lrt_cache.", call. = FALSE)
  }
  if (!exists(key, envir = fitobj$lrt_cache, inherits = FALSE)) {
    dds_work <- .copy_deseq_dds_for_lrt(fitobj$dds)
    dds_lrt <- DESeq2::DESeq(
      dds_work, test = "LRT", reduced = reduced_f, quiet = TRUE
    )
    assign(key, dds_lrt, envir = fitobj$lrt_cache)
  }
  get(key, envir = fitobj$lrt_cache, inherits = FALSE)
}

# DESeq2 dislikes hyphens in factor level names.
.sanitize_deseq2_levels <- function(x) {
  gsub("-", "_", as.character(x), fixed = TRUE)
}

# EdgeR-style expression -> numeric contrast vector aligned to colnames(design).
# Supports a single coef name (test vs 0) or "A - B".
# Under treatment coding, reference levels are absorbed in the intercept, so
# "ConditionB - ConditionA" with only ConditionB present becomes +1 on ConditionB.
.edger_expr_to_vector <- function(design, expr) {
  cn <- colnames(design)
  expr <- trimws(as.character(expr)[1])
  if (!nzchar(expr)) {
    stop("Empty EdgeR-style contrast expression.", call. = FALSE)
  }
  cv <- setNames(rep(0, length(cn)), cn)
  if (expr %in% cn) {
    cv[expr] <- 1
    return(unname(as.numeric(cv)))
  }
  if (grepl("-", expr, fixed = TRUE)) {
    parts <- trimws(strsplit(expr, "\\s*-\\s*")[[1]])
    if (length(parts) == 2L) {
      a <- parts[1]
      b <- parts[2]
      a_ok <- a %in% cn
      b_ok <- b %in% cn
      if (a_ok && b_ok) {
        cv[a] <- 1
        cv[b] <- -1
        return(unname(as.numeric(cv)))
      }
      # Treatment coding: reference level missing from design -> single coef vs intercept
      if (a_ok && !b_ok && "(Intercept)" %in% cn) {
        cv[a] <- 1
        return(unname(as.numeric(cv)))
      }
      if (!a_ok && b_ok && "(Intercept)" %in% cn) {
        # Rare: user wrote ref - nonref; flip sign
        cv[b] <- -1
        return(unname(as.numeric(cv)))
      }
      miss <- setdiff(parts, cn)
      stop(
        "EdgeR contrast '", expr, "' references unknown design column(s): ",
        paste(miss, collapse = ", "),
        ". Available: ", paste(cn, collapse = ", "),
        call. = FALSE
      )
    }
  }
  stop(
    "Could not parse EdgeR-style contrast '", expr,
    "'. Use a design coefficient name (e.g. 'ConditionKO1' or ",
    "'ConditionCOVID_SEV:BatchJa005E') or 'A - B'. Available columns: ",
    paste(cn, collapse = ", "),
    call. = FALSE
  )
}

# Map EdgeR-style contrast expression string to numeric glmLRT contrast.
.edgR_contrast_vector <- function(design, contrast) {
  if (is.character(contrast) && length(contrast) == 1L) {
    return(.edger_expr_to_vector(design, contrast))
  }
  stop(
    "EdgeR contrast must be a single design-coefficient expression string.",
    call. = FALSE
  )
}

# Fraction of cells expressing each gene (handles sparse or dense assay matrices).
.pct_cells_expressing <- function(mat, cells) {
  if (length(cells) == 0) return(rep(0, nrow(mat)))
  sub <- mat[, colnames(mat) %in% cells, drop = FALSE]
  if (ncol(sub) == 0) return(rep(0, nrow(mat)))
  if (inherits(sub, "dgCMatrix")) {
    tabulate(sub@i + 1L, nrow(sub)) / ncol(sub)
  } else {
    rowSums(sub > 0) / ncol(sub)
  }
}

# Single-cell pct.1 / pct.2 for the test (c1) and reference (c0) conditions.
# Descriptive only — not adjusted for covariates in the pseudobulk model.
.compute_cellsexp_c1_c0 <- function(sobjint, clust, grouping_variable, c0, c1, assay, slot, genes) {
  md <- sobjint@meta.data
  md <- md[md[, grouping_variable] == clust, , drop = FALSE]
  md <- md[md$Condition %in% c(c0, c1), , drop = FALSE]
  if (!nrow(md)) return(NULL)
  mat <- Seurat::GetAssayData(sobjint[, rownames(md)], assay = assay, layer = slot)
  cells_c1 <- rownames(md[md$Condition == c1, , drop = FALSE])
  cells_c0 <- rownames(md[md$Condition == c0, , drop = FALSE])
  cellsexp <- data.frame(
    gene = rownames(mat),
    pct.1 = .pct_cells_expressing(mat, cells_c1),
    pct.2 = .pct_cells_expressing(mat, cells_c0),
    stringsAsFactors = FALSE
  )
  cellsexp$pct.diff <- cellsexp$pct.1 - cellsexp$pct.2
  cellsexp <- cellsexp[cellsexp$gene %in% genes, , drop = FALSE]
  cellsexp
}

# Require min_n pseudobulk samples per condition for a given contrast.
.has_min_replicates <- function(coldata, c0, c1, min_n = 2L) {
  ct <- table(coldata$Condition)
  v0 <- as.character(c0)
  v1 <- as.character(c1)
  if (!all(c(v0, v1) %in% names(ct))) return(FALSE)
  as.integer(ct[[v0]]) >= min_n && as.integer(ct[[v1]]) >= min_n
}

# Ranked gene weights for downstream GSEA: -log10(p) * sign(LFC), with underflow fix.
.add_pseudobulk_gene_weights <- function(res, lfc_col, pval_col) {
  scores <- -log10(res[[pval_col]])
  scores <- scores * sign(res[[lfc_col]])
  names(scores) <- res$gene_symbol
  scores <- sort(scores, decreasing = TRUE)
  logFC_vec <- res[[lfc_col]]
  names(logFC_vec) <- res$gene_symbol
  scores <- scDAPP::fix_underflow(scores, logFC_vec)
  scores <- scores[match(res$gene_symbol, names(scores))]
  res$weight <- scores
  res[order(res$weight, decreasing = TRUE), , drop = FALSE]
}

# Resolved contrast as a single string for output columns.
.format_comp_contrast_string <- function(contrast_entry, c1, c0, DE_test = NULL) {
  if (is.null(DE_test)) {
    # Prefer storing whatever was already resolved / written on comps
    if (!.is_blank_contrast(contrast_entry)) {
      return(trimws(as.character(contrast_entry)[1]))
    }
    return(paste(.parse_comp_contrast(contrast_entry, c1, c0), collapse = ";"))
  }
  style <- .contrast_style(DE_test)
  parsed <- .parse_contrast(contrast_entry, c0, c1, style)
  if (identical(style, "ignore") || is.na(parsed$display)) {
    return(paste0("c1=", c1, ";c0=", c0))
  }
  parsed$display
}

# Bind per-cluster DE tables into one long data.frame.
.flatten_de_results <- function(comps, res_list) {
  pieces <- list()
  for (i in seq_along(res_list)) {
    clust_list <- res_list[[i]]
    if (is.null(clust_list) || !length(clust_list)) next
    contrast_str <- .format_comp_contrast_string(comps$contrast[i], comps$c1[i], comps$c0[i])
    meta <- data.frame(
      label = comps$label[i],
      formula = comps$formula[i],
      c0 = comps$c0[i],
      c1 = comps$c1[i],
      contrast = contrast_str,
      stringsAsFactors = FALSE
    )
    for (clust_name in names(clust_list)) {
      res <- clust_list[[clust_name]]
      if (is.null(res) || !is.data.frame(res) || !nrow(res)) next
      meta_rep <- meta[rep(1L, nrow(res)), , drop = FALSE]
      row_df <- cbind(meta_rep, cluster = clust_name, res, stringsAsFactors = FALSE)
      pieces[[length(pieces) + 1L]] <- row_df
    }
  }
  if (!length(pieces)) {
    return(data.frame(
      label = character(), formula = character(), c0 = character(), c1 = character(),
      contrast = character(), cluster = character(),
      stringsAsFactors = FALSE
    ))
  }
  out <- dplyr::bind_rows(pieces)
  meta_cols <- c("label", "formula", "c0", "c1", "contrast", "cluster")
  other_cols <- setdiff(colnames(out), meta_cols)
  out[, c(meta_cols, other_cols), drop = FALSE]
}

# Significant DEGs using the same rules as the numDEGs summary (FDR, |logFC|, pct).
.significant_de_index <- function(df, padj_thres, lfc_thres, min_pct) {
  if (!nrow(df)) return(logical(0))
  base <- df$FDR < padj_thres & abs(df$logFC) > lfc_thres
  up <- df$logFC > 0 & df$pct.1 > min_pct
  down <- df$logFC < 0 & df$pct.2 > min_pct
  base & (up | down)
}

.filter_significant_de_results <- function(df, padj_thres, lfc_thres, min_pct) {
  df[.significant_de_index(df, padj_thres, lfc_thres, min_pct), , drop = FALSE]
}

#' Default cross-condition DE thresholds for DEG counting (and ORA gene sets)
#'
#' @param Pseudobulk_mode logical; \code{TRUE} uses lenient pseudobulk defaults,
#'   \code{FALSE} uses stricter Wilcox defaults.
#' @param padj optional adjusted P value threshold; \code{NULL} uses the mode default.
#' @param lfc optional absolute logFC threshold; \code{NULL} uses the mode default.
#' @param min_pct optional minimum expression fraction; \code{NULL} uses the mode default.
#' @return Named list with \code{padj}, \code{lfc}, and \code{min_pct}.
#' @export
crosscondition_de_threshold_defaults <- function(Pseudobulk_mode,
                                                   padj = NULL,
                                                   lfc = NULL,
                                                   min_pct = NULL) {
  if (isTRUE(Pseudobulk_mode)) {
    def_padj <- 0.1
    def_lfc <- 0
    def_min_pct <- 0.1
  } else {
    def_padj <- 0.05
    def_lfc <- 0.25
    def_min_pct <- 0
  }
  list(
    padj = if (is.null(padj)) def_padj else padj,
    lfc = if (is.null(lfc)) def_lfc else lfc,
    min_pct = if (is.null(min_pct)) def_min_pct else min_pct
  )
}

#' Filter cross-condition DE results to significant DEGs
#'
#' Uses the same rules as \code{count_crosscondition_degs()} and
#' \code{de_across_conditions_module()} significant-gene output.
#'
#' @param df data.frame of DE results for one or more genes (harmonized columns).
#' @param padj_thres adjusted P value threshold (strict \code{<}).
#' @param lfc_thres absolute logFC threshold (strict \code{>}).
#' @param min_pct minimum \code{pct.1} (up) or \code{pct.2} (down) threshold (strict \code{>}).
#' @return Subset of \code{df} with only significant DEG rows.
#' @export
filter_significant_de_results <- function(df, padj_thres, lfc_thres, min_pct) {
  .filter_significant_de_results(df, padj_thres, lfc_thres, min_pct)
}

#' Count significant DEGs per cluster for a cross-condition comparison
#'
#' @param de_by_cluster Named list of per-cluster DE tables (e.g. from
#'   \code{de_results_by_cluster()}) or a single cluster table.
#' @param padj_thres adjusted P value threshold.
#' @param lfc_thres absolute logFC threshold.
#' @param min_pct minimum expression fraction threshold.
#' @param c0 reference condition label (column for genes higher in \code{c0}, negative LFC).
#' @param c1 test condition label (column for genes higher in \code{c1}, positive LFC).
#' @return Matrix with rownames = cluster names and columns \code{c0}, \code{c1} (counts).
#' @export
count_crosscondition_degs <- function(de_by_cluster,
                                    padj_thres,
                                    lfc_thres,
                                    min_pct,
                                    c0,
                                    c1) {
  if (is.data.frame(de_by_cluster)) {
    de_by_cluster <- list(cluster = de_by_cluster)
  }
  clusters <- names(de_by_cluster)
  out <- matrix(0L, nrow = length(clusters), ncol = 2L,
                dimnames = list(clusters, c(c0, c1)))
  for (cl in clusters) {
    m <- de_by_cluster[[cl]]
    if (is.null(m) || !nrow(m)) next
    sig <- .filter_significant_de_results(m, padj_thres, lfc_thres, min_pct)
    if (!nrow(sig)) next
    out[cl, c0] <- sum(sig$logFC < 0, na.rm = TRUE)
    out[cl, c1] <- sum(sig$logFC > 0, na.rm = TRUE)
  }
  out
}

#' Diverging barplot of significant DEG counts per cluster
#'
#' Visualizes output from \code{count_crosscondition_degs()} as a horizontal
#' diverging barplot: clusters on the Y axis, signed DEG counts on the X axis
#' (negative = higher in \code{c0}, positive = higher in \code{c1}).
#'
#' @param numdegs matrix from \code{count_crosscondition_degs()} or a data.frame
#'   with a cluster column and two count columns for \code{c0} and \code{c1}.
#' @param c0 reference condition label.
#' @param c1 test condition label.
#' @param cluster_levels optional character vector for Y-axis cluster order.
#' @param title optional plot title (e.g. comparison label).
#' @return A \code{ggplot} object.
#' @export
plot_crosscondition_deg_barplot <- function(numdegs,
                                            c0,
                                            c1,
                                            cluster_levels = NULL,
                                            title = NULL) {
  if (is.matrix(numdegs)) {
    counts_df <- data.frame(
      Cluster = rownames(numdegs),
      numdegs,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    counts_df <- numdegs
  }

  cluster_col <- if ("Cluster" %in% colnames(counts_df)) {
    "Cluster"
  } else if ("cluster" %in% colnames(counts_df)) {
    "cluster"
  } else {
    stop("numdegs must include a Cluster (or cluster) column when passed as a data.frame.")
  }

  c0_col <- if (c0 %in% colnames(counts_df)) {
    c0
  } else {
    grep(paste0("^", c0), colnames(counts_df), value = TRUE)[1]
  }
  c1_col <- if (c1 %in% colnames(counts_df)) {
    c1
  } else {
    grep(paste0("^", c1), colnames(counts_df), value = TRUE)[1]
  }
  if (is.na(c0_col) || is.na(c1_col)) {
    stop("Could not find count columns for c0 and c1 in numdegs.")
  }

  if (is.null(cluster_levels)) {
    cluster_levels <- counts_df[[cluster_col]]
  }

  counts_df[[cluster_col]] <- factor(counts_df[[cluster_col]], levels = cluster_levels)
  counts_df[[c0_col]] <- as.numeric(counts_df[[c0_col]])
  counts_df[[c1_col]] <- as.numeric(counts_df[[c1_col]])
  counts_df[[c0_col]][is.na(counts_df[[c0_col]])] <- 0
  counts_df[[c1_col]][is.na(counts_df[[c1_col]])] <- 0

  plot_df <- rbind(
    data.frame(
      Cluster = counts_df[[cluster_col]],
      direction = c0,
      signed_count = -counts_df[[c0_col]],
      stringsAsFactors = FALSE
    ),
    data.frame(
      Cluster = counts_df[[cluster_col]],
      direction = c1,
      signed_count = counts_df[[c1_col]],
      stringsAsFactors = FALSE
    )
  )
  plot_df$direction <- factor(plot_df$direction, levels = c(c0, c1))

  fill_vals <- stats::setNames(
    c("#2166AC", "#B2182B"),
    c(c0, c1)
  )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$signed_count,
      y = .data$Cluster,
      fill = .data$direction
    )
  ) +
    ggplot2::geom_col(position = "identity", width = 0.75) +
    ggplot2::geom_vline(xintercept = 0, color = "grey30", linewidth = 0.4) +
    ggplot2::scale_x_continuous(
      labels = function(x) abs(x),
      expand = ggplot2::expansion(mult = c(0.05, 0.05))
    ) +
    ggplot2::scale_y_discrete(limits = rev(cluster_levels)) +
    ggplot2::scale_fill_manual(
      values = fill_vals,
      name = "Higher in"
    ) +
    ggplot2::labs(
      title = title,
      x = "Number of significant DEGs",
      y = NULL,
      caption = paste0("Left (negative): higher in ", c0, "  |  Right (positive): higher in ", c1)
    ) +
    ggplot2::theme_linedraw() +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  p
}

#' Select top up/down cross-condition DEG genes for one cluster
#'
#' @param de_table data.frame of DE results for one cluster.
#' @param n_top maximum genes per direction (default 10).
#' @param padj_thres adjusted P value threshold.
#' @param lfc_thres absolute logFC threshold.
#' @param min_pct minimum expression fraction threshold.
#' @return List with \code{up}, \code{down}, and \code{all} gene symbols (up block first).
#' @export
select_top_crosscondition_deg_genes <- function(de_table,
                                                n_top = 10L,
                                                padj_thres,
                                                lfc_thres,
                                                min_pct) {
  n_top <- as.integer(n_top)[1]
  if (is.null(de_table) || !nrow(de_table)) {
    return(list(up = character(), down = character(), all = character()))
  }
  sig <- .filter_significant_de_results(de_table, padj_thres, lfc_thres, min_pct)
  .select_top_crosscondition_deg_genes_core(sig, n_top)
}

#' Select top up/down DEG genes pooled across clusters for one comparison
#'
#' Significant genes are deduplicated by \code{gene_symbol}, keeping the row with
#' the largest absolute logFC.
#'
#' @param de_by_cluster named list of per-cluster DE tables.
#' @param n_top maximum genes per direction (default 10).
#' @param padj_thres adjusted P value threshold.
#' @param lfc_thres absolute logFC threshold.
#' @param min_pct minimum expression fraction threshold.
#' @return List with \code{up}, \code{down}, and \code{all} gene symbols.
#' @export
select_top_crosscondition_deg_genes_pooled <- function(de_by_cluster,
                                                       n_top = 10L,
                                                       padj_thres,
                                                       lfc_thres,
                                                       min_pct) {
  n_top <- as.integer(n_top)[1]
  if (is.data.frame(de_by_cluster)) {
    de_by_cluster <- list(cluster = de_by_cluster)
  }
  pieces <- lapply(de_by_cluster, function(m) {
    if (is.null(m) || !nrow(m)) {
      return(NULL)
    }
    .filter_significant_de_results(m, padj_thres, lfc_thres, min_pct)
  })
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  if (!length(pieces)) {
    return(list(up = character(), down = character(), all = character()))
  }
  sig <- dplyr::bind_rows(pieces)
  if (!nrow(sig)) {
    return(list(up = character(), down = character(), all = character()))
  }
  sig <- sig[order(-abs(sig$logFC), sig$gene_symbol), , drop = FALSE]
  sig <- sig[!duplicated(sig$gene_symbol), , drop = FALSE]
  .select_top_crosscondition_deg_genes_core(sig, n_top)
}

.select_top_crosscondition_deg_genes_core <- function(sig, n_top) {
  if (!nrow(sig)) {
    return(list(up = character(), down = character(), all = character()))
  }

  has_weight <- "weight" %in% colnames(sig)

  up_df <- sig[sig$logFC > 0, , drop = FALSE]
  if (nrow(up_df)) {
    if (has_weight) {
      up_df <- up_df[order(-up_df$logFC, -up_df$weight, up_df$gene_symbol), , drop = FALSE]
    } else {
      up_df <- up_df[order(-up_df$logFC, up_df$gene_symbol), , drop = FALSE]
    }
    up_genes <- head(up_df$gene_symbol, n_top)
  } else {
    up_genes <- character()
  }

  dn_df <- sig[sig$logFC < 0, , drop = FALSE]
  if (nrow(dn_df)) {
    if (has_weight) {
      dn_df <- dn_df[order(dn_df$logFC, dn_df$weight, dn_df$gene_symbol), , drop = FALSE]
    } else {
      dn_df <- dn_df[order(dn_df$logFC, dn_df$gene_symbol), , drop = FALSE]
    }
    dn_genes <- head(dn_df$gene_symbol, n_top)
  } else {
    dn_genes <- character()
  }

  list(up = up_genes, down = dn_genes, all = c(up_genes, dn_genes))
}

#' Seurat DotPlot of top cross-condition DEGs by sample Code
#'
#' @param sobj Seurat object (integrated).
#' @param genes character vector of gene symbols to plot.
#' @param sample_codes sample \code{Code} values to include (typically from the comparison).
#' @param c0 reference condition label.
#' @param c1 test condition label.
#' @param assay assay name for expression.
#' @param slot assay layer/slot (default \code{"data"}).
#' @param title optional plot title.
#' @param cluster_id optional cluster id (without \code{cluster_} prefix); \code{NULL} uses all cells in the comparison.
#' @param sample_metadata optional metadata with \code{Code} and \code{Condition} for sample ordering.
#' @return A \code{ggplot} object, or \code{NULL} if nothing to plot.
#' @export
plot_crosscondition_deg_dotplot <- function(sobj,
                                            genes,
                                            sample_codes,
                                            c0,
                                            c1,
                                            assay,
                                            slot = "data",
                                            title = NULL,
                                            cluster_id = NULL,
                                            sample_metadata = NULL) {
  if (is.null(genes) || !length(genes)) {
    return(NULL)
  }
  genes <- unique(as.character(genes))
  genes <- genes[nzchar(genes)]
  if (!length(genes)) {
    return(NULL)
  }

  md <- sobj@meta.data
  cells <- rownames(md)
  if (!is.null(cluster_id)) {
    keep <- md$seurat_clusters == cluster_id
    cells <- cells[keep]
    md <- md[cells, , drop = FALSE]
  }
  keep <- md$Condition %in% c(c0, c1)
  cells <- cells[keep]
  md <- md[cells, , drop = FALSE]
  if (!length(cells)) {
    return(NULL)
  }
  if (!is.null(sample_codes)) {
    keep <- md$Code %in% sample_codes
    cells <- cells[keep]
    md <- md[cells, , drop = FALSE]
  }
  if (!length(cells)) {
    return(NULL)
  }

  if (!is.null(sample_metadata)) {
    subpmd <- sample_metadata[sample_metadata$Condition %in% c(c1, c0), , drop = FALSE]
    code_order <- c(
      subpmd$Code[subpmd$Condition == c1],
      subpmd$Code[subpmd$Condition == c0]
    )
    code_order <- code_order[code_order %in% unique(md$Code)]
  } else {
    code_order <- c(
      unique(md$Code[md$Condition == c1]),
      unique(md$Code[md$Condition == c0])
    )
  }

  sub <- sobj[, cells]
  sub$Code <- factor(as.character(sub$Code), levels = code_order)

  if (!(assay %in% Seurat::Assays(sub))) {
    return(NULL)
  }
  genes_in <- genes[genes %in% rownames(sub[[assay]])]
  if (!length(genes_in)) {
    return(NULL)
  }

  dot_args <- list(
    object = sub,
    features = rev(genes_in),
    group.by = "Code",
    assay = assay
  )
  dot_formals <- names(formals(Seurat::DotPlot))
  if ("layer" %in% dot_formals) {
    dot_args$layer <- slot
  } else if ("slot" %in% dot_formals) {
    dot_args$slot <- slot
  }

  p <- do.call(Seurat::DotPlot, dot_args) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(
      title = title,
      caption = paste0("Positive logFC genes higher in ", c1, "; negative logFC higher in ", c0)
    )

  p
}

# Safe filename stem encoding run parameters.
.de_output_basename <- function(DE_test, grouping_variable, assay, slot, pseudobulk_mode) {
  mode_chr <- if (isTRUE(pseudobulk_mode)) "pseudobulk" else "singlecell"
  sanitize <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  paste(
    "crosscondition_DE",
    mode_chr,
    paste0("test_", sanitize(DE_test)),
    paste0("grp_", sanitize(grouping_variable)),
    paste0("assay_", sanitize(assay)),
    paste0("layer_", sanitize(slot)),
    sep = "_"
  )
}

.de_thresholds_suffix <- function(padj_thres, lfc_thres, min_pct) {
  sanitize <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  paste0(
    "_padj", sanitize(padj_thres),
    "_lfc", sanitize(lfc_thres),
    "_minpct", sanitize(min_pct)
  )
}

# Per-comparison x cluster DEG counts (significant genes only).
.write_numdegs_summary <- function(de_sig, comps, groupinglev_nicelabs, filepath) {
  rows <- lapply(seq_len(nrow(comps)), function(i) {
    lab <- comps$label[i]
    c0 <- comps$c0[i]
    c1 <- comps$c1[i]
    sub <- de_sig[de_sig$label == lab, , drop = FALSE]
    lapply(groupinglev_nicelabs, function(cl) {
      m <- sub[sub$cluster == cl, , drop = FALSE]
      data.frame(
        label = lab,
        cluster = cl,
        c0 = c0,
        c1 = c1,
        n_down = sum(m$logFC < 0, na.rm = TRUE),
        n_up = sum(m$logFC > 0, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
  })
  numdegs_all <- dplyr::bind_rows(unlist(rows, recursive = FALSE))
  write.csv(numdegs_all, filepath, quote = FALSE, row.names = FALSE)
  invisible(numdegs_all)
}

#' Differential expression analysis across conditions for integrated Seurat objects
#'
#' This is a modular component of the scDAPP scRNAseq pipeline. Perform DE analysis across conditions. Supports A vs B vs C pairwise (multiple conditions) comparisons. Options for Pseudobulk DE via EdgeR - LRT, DREAM (paired / mixed models), or old-school scRNAseq DE via wilcoxon test.
#'
#' @param sobjint integrated Seurat object. metadata needs two special columns: one called "Condition" that contains the A vs B conditions, and a second that matches the  `grouping_variable` parameter of this function.
#' @param sample_metadata data.frame with three columns called Sample, Condition, Code.
#' @param comps data.frame defining comparisons. Required: \code{c0} (reference) and
#'   \code{c1} (test; log2FC is c1 vs c0). Optional: \code{formula} (default
#'   \code{~ Condition}), \code{contrast}, and \code{label} (default
#'   \code{c1_vs_c0}). For paired pseudobulk DE with \code{DE_test = "Dream"}, include
#'   an intercept random effect such as \code{~ Condition + (1|Patient)}. Deprecated
#'   \code{c2} is renamed to \code{c0}.
#' @param grouping_variable string, column name of identity in Seurat object meta.data to stratify DE by. For example, clusters or celltype. Will perform A vs B DE in each of these groupings. Default is "seurat_clusters". It is not mandatory, but will use factor level ordering of this variable in the meta.data to control analysis order, and if not will sort by alphanumeric order (cluster 1, then 2, cluster A, then B, etc)
#' @param Pseudobulk_mode T/F. Sets the cross-conditional analysis mode. TRUE uses pseudobulk EdgeR for DE testing and propeller for compositional analysis. FALSE uses single-cell wilcox test within Seurat for DE testing and 2-prop Z test within the `prop.test()` function for compositional analysis.
#' @param DE_test a string, default is 'EdgeR-LRT' when Pseudobulk_mode is set to True, or 'wilcox' when Pseudobulk_mode is False. Can be "DESeq2", "DESeq2-LRT", "EdgeR", "EdgeR-LRT", "EdgeR-QLF", or "Dream" for pseudobulk (Dream requires \code{variancePartition} and a \code{(1|var)} term in \code{comps$formula}), or any of the tests supported by the "test.use" argument in the FindMarkers function in Seurat; see `?Seurat::FindMarkers` for more. Note the Seurat "roc" test is not included, and some additional packages like DESeq2 may require installation.
#' @param outdir_int optional path to save CSVs under
#'   \code{differentialexpression_crosscondition/}. Writes
#'   natural method output column names to \code{{params}_all.csv},
#'   \code{{params}_{thresholds}_significant.csv}, and
#'   \code{{params}_numDEGs_summary.csv} where \code{params} encodes
#'   \code{DE_test}, pseudobulk mode, \code{grouping_variable}, \code{assay}, and
#'   \code{slot} (layer). The returned data.frame is harmonized to edgeR-style
#'   column names for downstream code.
#' @param assay string, name of Seurat assay to use, default is DefaultAssay(sobjint)
#' @param slot string, name of Seurat assay slot to use, default is 'data'
#' @param cluster_prefix T/F, whether to append prefix "cluster_" to grouping_variable levels, useful if grouping_variable is a cluster. Rather than saving results with names such as "1", "2", "3", will save as "cluster_1", and so on.
#' @param crossconditionDE_padj_thres numeric; adjusted P value maximum threshold for calling DEGs (counting and ORA). If missing, 0.1 when \code{Pseudobulk_mode} is \code{TRUE}, 0.05 when \code{FALSE}.
#' @param crossconditionDE_lfc_thres numeric; absolute logFC minimum for calling DEGs. If missing, 0 when pseudobulk, 0.25 when Wilcox.
#' @param crossconditionDE_min.pct numeric; minimum \code{pct.1} (up) or \code{pct.2} (down) to count as a DEG. If missing, 0.1 when pseudobulk, 0 when Wilcox. See \code{crosscondition_de_threshold_defaults()}.
#' @param workernum integer; parallel workers for Dream (\code{BiocParallel}); default 1.
#'
#' @return A single harmonized data.frame with one row per gene x cluster x comparison. Metadata
#'   columns: \code{label}, \code{formula}, \code{c0}, \code{c1}, \code{contrast},
#'   \code{cluster}; then DE statistics (\code{gene_symbol}, \code{logFC}, \code{PValue},
#'   \code{FDR}, \code{pct.1}, \code{pct.2}, \code{pct.diff}, \code{weight}, and optional
#'   count columns). \code{weight} is \eqn{-\log_{10}(p) \times \mathrm{sign}(log2FC)}.
#' @export
#'
#' @examples
#' \dontrun{
#' de_results <- de_across_conditions_module(
#'  sobjint = sobjint,
#'  sample_metadata = sample_metadata,
#'  comps = comps,
#'  outdir_int = outdir_int,
#'  grouping_variable = 'seurat_clusters',
#'  Pseudobulk_mode = T
#'  )
#' }
de_across_conditions_module <- function(sobjint,
                                        sample_metadata,
                                        comps,
                                        grouping_variable,
                                        Pseudobulk_mode,
                                        DE_test,
                                        outdir_int,
                                        assay,
                                        slot,
                                        cluster_prefix,
                                        crossconditionDE_padj_thres,
                                        crossconditionDE_lfc_thres,
                                        crossconditionDE_min.pct,
                                        workernum = 1L
                                        
){
  
  # ---------------------------------------------------------------------------
  # Argument defaults
  # ---------------------------------------------------------------------------
  if( missing(sobjint)) { stop('Provide Seurat object') }
  if( missing(grouping_variable)) { grouping_variable <- 'seurat_clusters' }
  # if( missing(min.pct)) { stop(min.pct <- 0.1) } # for GSEA, don't do this
  
  if(missing(assay)){assay <- Seurat::DefaultAssay(sobjint)}
  if(missing(slot)){slot <- 'data'}
  if(missing(cluster_prefix)){cluster_prefix <- NULL}
  
  th <- crosscondition_de_threshold_defaults(
    Pseudobulk_mode,
    padj = if (missing(crossconditionDE_padj_thres)) NULL else crossconditionDE_padj_thres,
    lfc = if (missing(crossconditionDE_lfc_thres)) NULL else crossconditionDE_lfc_thres,
    min_pct = if (missing(crossconditionDE_min.pct)) NULL else crossconditionDE_min.pct
  )
  crossconditionDE_padj_thres <- th$padj
  crossconditionDE_lfc_thres <- th$lfc
  crossconditionDE_min.pct <- th$min_pct
  
  if(missing(DE_test)){
    if(Pseudobulk_mode == T){DE_test = 'EdgeR-LRT'}
    if(Pseudobulk_mode == F){DE_test = 'wilcox'}
  }
  if (missing(workernum) || is.null(workernum)) workernum <- 1L
  workernum <- as.integer(workernum)[1]
  if (is.na(workernum) || workernum < 1L) workernum <- 1L
  
  # Standardize comps (c0=reference, c1=test, formula, contrast, label).
  comps <- .apply_contrast_defaults(.normalize_comps(comps), DE_test)
  if (isTRUE(Pseudobulk_mode)) {
    .preflight_comps_design(sample_metadata, comps, DE_test, Pseudobulk_mode = TRUE)
    .validate_comps_de_formulas(comps, DE_test)
    if (identical(DE_test, "Dream") &&
        !requireNamespace("variancePartition", quietly = TRUE)) {
      stop(
        "DE_test = 'Dream' requires the Bioconductor package variancePartition. ",
        "Install with BiocManager::install('variancePartition').",
        call. = FALSE
      )
    }
  } else {
    .preflight_comps_design(sample_metadata, comps, DE_test, Pseudobulk_mode = FALSE)
  }
  
  # ---------------------------------------------------------------------------
  # Grouping variable (clusters / cell types to stratify DE)
  # ---------------------------------------------------------------------------
  groupingvec <- sobjint@meta.data[,grouping_variable]
  if(!is.factor(groupingvec)){
    warning('The grouping variable (ie clusters within which to perform A vs B DE) is not a factor.\nThe the comparison order will default to alphanumeric order.')
    
    groupinglevs <- unique(groupingvec)
    groupinglevs <- stringr::str_sort(groupinglevs, numeric=T)
    groupingvec <- factor(groupingvec, levels = groupinglevs)
    
  }
  
  
  #get the actual clusters (grouping levels)
  groupinglevs <- levels( groupingvec )
  
  #remove empty levels
  groupinglevs <- groupinglevs[groupinglevs %in% groupingvec]
  
  # for later, if they are clusters, we want better labels than just numerics
  # try to check if they are clusters; max str len will probably be 3 (in huge datasets...)
  if( is.null(cluster_prefix) ) {
    
    if( max(stringr::str_length(groupinglevs)) <= 3 ){cluster_prefix <- T} else{cluster_prefix = F}
    
  }
  
  if(cluster_prefix==T){
    groupinglev_nicelabs <- paste0('cluster_', groupinglevs)
  } else{groupinglev_nicelabs <- groupinglevs}
  
  
  # ===========================================================================
  # PSEUDOBULK MODE (edgeR / DESeq2)
  # ===========================================================================
  if(Pseudobulk_mode == T){
    
    # --- Step 1: pseudobulk counts per sample x cluster ---
    # One matrix per sample (rows=genes, cols=grouping levels).
    # min_cells=0 so rare clusters are kept; zeros added later if missing.
    samples <- sample_metadata$Code
    samp <- samples[1]
    
    pblist <- lapply(samples, function(samp){
      
      # message('\n\n',samp, '\n')
      
      md <- sobjint@meta.data
      md <- md[md$Code == samp,]
      cells <- rownames(md)
      sobjint_ct <- sobjint[,cells]
      
      # legacy note: previous code experimented with unlogging integrated assay data.
      # UPDATE feb 18 2025 - i think the line above is irrelevant, we unlog RISC values way earlier
      
      
      suppressMessages(
        pb <- scDAPP::pseudobulk(sobjint_ct,
                                 assay = assay,
                                 slot = slot,
                                 grouping_colname_in_md = grouping_variable,
                                 min_cells = 0)
      )
      
      
      ### round it for edgeR
      pb <- round(pb)
      
      pb
      
    })
    
    
    names(pblist) <- samples
    
    
    # --- Step 2: pad missing cluster columns with zeros (per sample) ---
    # Ensures every sample has the same cluster columns for binding.
    pblist <- lapply(pblist, function(pb){
      if(any(!(groupinglevs %in% colnames(pb)))){
        fakectcols <- lapply(groupinglevs, function(ct){
          if(!(ct %in% colnames(pb))){
            fakectcol <- data.frame(ct = rep(0, nrow(pb)))
            colnames(fakectcol) <- ct
            fakectcol
          }
        })
        fakectcols <- Filter(Negate(is.null), fakectcols)
        if (length(fakectcols)) {
          fakectcolsdf <- dplyr::bind_cols(fakectcols)
          pb <- cbind(pb, fakectcolsdf)
        }
      }
      
      
      
      #make sure the clusters are ordered properly
      # (ie if we ahve 8 clusters and cluster 5 was missing)
      pb <- pb[,match(groupinglevs,colnames(pb))]
      
      pb
      
    })
    
    
    pblist_overall <- pblist
    
    # --- Step 3: joint DE per design formula ---
    # Group comps rows that share the same formula (e.g. ~ Condition + Batch).
    # For each formula: fit once per cluster on ALL samples, then extract each
    # comps contrast (c1 vs c0) without re-subsetting counts.
    is_edger <- DE_test %in% c('EdgeR', 'EdgeR-LRT', 'EdgeR-QLF')
    is_deseq <- DE_test %in% c('DESeq2', 'DESeq2-LRT')
    is_dream <- identical(DE_test, "Dream")
    if (is_edger) require(edgeR)
    if (is_deseq) require(DESeq2)
    
    formula_groups <- split(seq_len(nrow(comps)), comps$formula)
    m_bycluster_crosscondition_de_comps <- vector("list", nrow(comps))
    
    for (formula_chr in names(formula_groups)) {
      row_idx <- formula_groups[[formula_chr]]
      design_formula <- .parse_comp_formula(formula_chr)
      .validate_de_formula(DE_test, design_formula)
      coldata <- .prepare_pseudobulk_coldata(sample_metadata, design_formula)
      pblist <- pblist_overall[match(coldata$Code, names(pblist_overall))]
      
      # Dream contrast expressions for all comps sharing this formula
      dream_contrast_exprs <- NULL
      dream_coef_names <- character()
      if (is_dream) {
        dream_contrast_exprs <- vapply(row_idx, function(compidx) {
          .dream_contrast_expr(
            comps$contrast[compidx], comps$c1[compidx], comps$c0[compidx]
          )
        }, character(1))
        dream_coef_names <- make.names(comps$label[row_idx], unique = TRUE)
        names(dream_contrast_exprs) <- dream_coef_names
      }
      
      # --- Step 3a: one edgeR/DESeq2/Dream fit per cluster (all samples) ---
      cluster_fits <- setNames(
        lapply(seq_along(groupinglevs), function(i) {
          clust <- groupinglevs[i]
          message(clust)
          
          # genes x samples count matrix for this cluster
          gemlist <- lapply(names(pblist), function(samp) {
            pb <- pblist[[samp]]
            pbcol <- pb[, colnames(pb) == clust, drop = FALSE]
            colnames(pbcol) <- samp
            pbcol
          })
          gem <- dplyr::bind_cols(gemlist)
          gem <- gem[Matrix::rowSums(gem) > 3, , drop = FALSE]      # low-expression genes
          gem <- gem[, Matrix::colSums(gem) > 10, drop = FALSE]     # empty pseudobulk samples
          if (!ncol(gem) || !nrow(gem)) return(NULL)
          
          fit_coldata <- coldata[match(colnames(gem), coldata$Code), , drop = FALSE]
          rownames(fit_coldata) <- fit_coldata$Code
          
          if (is_dream) {
            out_dream <- tryCatch({
              gem_mat <- as.matrix(gem)
              bp <- .dream_bpparam(workernum)
              dream_formula <- .dream_fit_formula(design_formula)
              L_fit <- variancePartition::makeContrastsDream(
                dream_formula, fit_coldata, contrasts = dream_contrast_exprs
              )
              vobj <- variancePartition::voomWithDreamWeights(
                gem_mat, dream_formula, fit_coldata, BPPARAM = bp
              )
              fit_mm <- variancePartition::dream(
                vobj, dream_formula, fit_coldata, L = L_fit, BPPARAM = bp
              )
              fit_mm <- variancePartition::eBayes(fit_mm)
              y_norm <- edgeR::DGEList(counts = gem_mat)
              y_norm <- edgeR::calcNormFactors(y_norm)
              list(
                type = "dream", fit = fit_mm, gem = gem_mat, coldata = fit_coldata,
                y = y_norm, coef_names = dream_coef_names
              )
            }, error = function(e) {
              warning(
                "Dream fit failed for cluster ", clust, " (", formula_chr, "): ",
                conditionMessage(e), "; skipping.",
                call. = FALSE
              )
              NULL
            })
            out_dream
          } else {
            fixed_formula <- .fixed_effects_formula(design_formula)
            design <- model.matrix(fixed_formula, data = fit_coldata)
            if (qr(design)$rank != ncol(design)) {
              warning("Design not full rank for cluster ", clust, " (", formula_chr, "); skipping.", call. = FALSE)
              return(NULL)
            }
            
            if (is_edger) {
              y <- DGEList(counts = gem, samples = fit_coldata,
                           group = fit_coldata$Condition)
              y <- calcNormFactors(y)
              y <- estimateDisp(y, design)
              if (identical(DE_test, "EdgeR-QLF")) {
                fit_glm <- glmQLFit(y, design)
                list(
                  type = "edger_qlf", y = y, design = design, fit = fit_glm,
                  gem = gem, coldata = fit_coldata
                )
              } else {
                fit_glm <- glmFit(y, design)
                list(
                  type = "edger", y = y, design = design, fit = fit_glm,
                  gem = gem, coldata = fit_coldata
                )
              }
            } else {
              col_dds <- fit_coldata
              vars <- all.vars(fixed_formula)
              for (v in vars) {
                if (is.factor(col_dds[[v]])) {
                  levels(col_dds[[v]]) <- .sanitize_deseq2_levels(levels(col_dds[[v]]))
                }
              }
              dds <- DESeqDataSetFromMatrix(gem, col_dds, design = fixed_formula)
              if (DE_test == 'DESeq2-LRT') {
                # Size factors + dispersions only; nested LRT runs per comps row at extract.
                dds <- DESeq2::estimateSizeFactors(dds)
                dds <- DESeq2::estimateDispersions(dds, quiet = TRUE)
                list(
                  type = "deseq2", dds = dds, gem = gem, coldata = col_dds,
                  use_lrt = TRUE, full_formula = fixed_formula,
                  lrt_cache = new.env(parent = emptyenv())
                )
              } else {
                dds <- DESeq(dds)
                list(
                  type = "deseq2", dds = dds, gem = gem, coldata = col_dds,
                  use_lrt = FALSE, full_formula = fixed_formula,
                  lrt_cache = NULL
                )
              }
            }
          }
        }),
        groupinglev_nicelabs
      )
      
      # --- Step 3b: extract each comparison contrast from the shared fits ---
      contrast_style <- .contrast_style(DE_test)
      for (j in seq_along(row_idx)) {
        compidx <- row_idx[j]
        c0 <- comps$c0[compidx]   # reference
        c1 <- comps$c1[compidx]   # test (log2FC = c1 vs c0)
        lab <- comps$label[compidx]
        message('\n', lab)
        parsed_contrast <- .parse_contrast(
          comps$contrast[compidx], c0, c1, style = contrast_style
        )
        comps$contrast[compidx] <- parsed_contrast$display
        dream_coef <- if (is_dream) dream_coef_names[j] else NULL
        
        m_bycluster_crosscondition_de <- lapply(seq_along(groupinglevs), function(i) {
          clust <- groupinglevs[i]
          fitobj <- cluster_fits[[i]]
          if (is.null(fitobj)) return(NULL)
          if (!.has_min_replicates(fitobj$coldata, c0, c1)) return(NULL)
          
          cellsexp <- .compute_cellsexp_c1_c0(
            sobjint, clust, grouping_variable, c0, c1, assay, slot, rownames(fitobj$gem)
          )
          if (is.null(cellsexp)) return(NULL)
          
          out <- tryCatch({
          # --- Dream: topTable for this contrast coefficient ---
          if (fitobj$type == "dream") {
            coef_use <- dream_coef
            if (is.null(coef_use) || !(coef_use %in% colnames(fitobj$fit$coefficients))) {
              # fall back to first contrast column if naming drifted
              coef_use <- colnames(fitobj$fit$coefficients)[1]
            }
            res <- as.data.frame(
              variancePartition::topTable(
                fitobj$fit, coef = coef_use, number = Inf, sort.by = "none"
              )
            )
            res <- cbind(rownames(res), res)
            colnames(res)[1] <- 'gene_symbol'
            # Harmonize limma column names toward edgeR-style
            if ("P.Value" %in% colnames(res)) {
              colnames(res)[colnames(res) == "P.Value"] <- "PValue"
            }
            if ("adj.P.Val" %in% colnames(res)) {
              colnames(res)[colnames(res) == "adj.P.Val"] <- "FDR"
            }
            cellsexp <- cellsexp[match(rownames(res), cellsexp$gene), , drop = FALSE]
            res <- cbind(res, cellsexp[, -1, drop = FALSE])
            res <- .add_pseudobulk_gene_weights(res, "logFC", "PValue")
            nc <- edgeR::cpm(fitobj$y)
            nc <- nc[match(rownames(res), rownames(nc)), , drop = FALSE]
            rc <- fitobj$gem[match(rownames(res), rownames(fitobj$gem)), , drop = FALSE]
            colnames(nc) <- paste0('normcounts_', colnames(nc))
            colnames(rc) <- paste0('rawcounts_', colnames(rc))
            cbind(res, nc, rc)
          } else if (fitobj$type %in% c("edger", "edger_qlf")) {
            # --- edgeR: glmLRT / glmQLFTest with contrast (or exactTest for simple 2-group) ---
            simple <- .is_simple_condition_design(design_formula)
            n_cond <- length(unique(fitobj$coldata$Condition))
            
            if (DE_test == 'EdgeR') {
              if (!simple) {
                warning(
                  "EdgeR exactTest requires formula ~ Condition only; skipping cluster ",
                  clust, " for ", lab, ". Use EdgeR-LRT or EdgeR-QLF with covariates.",
                  call. = FALSE
                )
                return(NULL)
              }
              if (n_cond > 2L) {
                warning("More than two Condition levels; using glmLRT (not exactTest) for ", lab, call. = FALSE)
                cv <- .edger_expr_to_vector(fitobj$design, parsed_contrast$edger_expr)
                lrt <- glmLRT(fitobj$fit, contrast = cv)
                res <- as.data.frame(topTags(lrt, n = Inf))
              } else {
                et <- exactTest(fitobj$y, pair = c(as.character(c0), as.character(c1)))
                res <- as.data.frame(topTags(et, n = Inf))
              }
            } else if (identical(DE_test, "EdgeR-QLF") || identical(fitobj$type, "edger_qlf")) {
              cv <- .edger_expr_to_vector(fitobj$design, parsed_contrast$edger_expr)
              qlf <- glmQLFTest(fitobj$fit, contrast = cv)
              res <- as.data.frame(topTags(qlf, n = Inf))
            } else {
              cv <- .edger_expr_to_vector(fitobj$design, parsed_contrast$edger_expr)
              lrt <- glmLRT(fitobj$fit, contrast = cv)
              res <- as.data.frame(topTags(lrt, n = Inf))
            }
            
            res <- cbind(rownames(res), res)
            colnames(res)[1] <- 'gene_symbol'
            cellsexp <- cellsexp[match(rownames(res), cellsexp$gene), , drop = FALSE]
            res <- cbind(res, cellsexp[, -1, drop = FALSE])
            res <- .add_pseudobulk_gene_weights(res, "logFC", "PValue")
            nc <- cpm(fitobj$y)
            nc <- nc[match(rownames(res), rownames(nc)), , drop = FALSE]
            rc <- fitobj$gem[match(rownames(res), rownames(fitobj$gem)), , drop = FALSE]
            colnames(nc) <- paste0('normcounts_', colnames(nc))
            colnames(rc) <- paste0('rawcounts_', colnames(rc))
            cbind(res, nc, rc)
          } else {
            # --- DESeq2: Wald, or per-comps-row nested LRT ---
            dds_use <- fitobj$dds
            if (isTRUE(fitobj$use_lrt) || identical(DE_test, "DESeq2-LRT")) {
              full_f <- if (!is.null(fitobj$full_formula)) {
                fitobj$full_formula
              } else {
                design_formula
              }
              reduced_f <- .deseq2_reduced_for_contrast(full_f, parsed_contrast)
              dds_use <- .deseq2_lrt_fit_cached(fitobj, reduced_f)
            }
            res_args <- list(object = dds_use)
            if (!is.null(parsed_contrast$deseq_name)) {
              res_args$name <- parsed_contrast$deseq_name
            } else {
              contrast_ds <- parsed_contrast$deseq_triple
              contrast_ds[2:3] <- .sanitize_deseq2_levels(contrast_ds[2:3])
              res_args$contrast <- contrast_ds
            }
            # After nested LRT, results() keeps LRT p-values; LFC follows contrast/name.
            res <- as.data.frame(do.call(DESeq2::results, res_args))
            res[is.na(res$padj), "padj"] <- 1
            res <- cbind(rownames(res), res)
            colnames(res)[1] <- 'gene_symbol'
            cellsexp <- cellsexp[match(rownames(res), cellsexp$gene), , drop = FALSE]
            res <- cbind(res, cellsexp[, -1, drop = FALSE])
            res <- .add_pseudobulk_gene_weights(res, "log2FoldChange", "pvalue")
            nc <- counts(dds_use, normalized = TRUE)
            rc <- counts(dds_use, normalized = FALSE)
            nc <- nc[match(rownames(res), rownames(nc)), , drop = FALSE]
            rc <- rc[match(rownames(res), rownames(rc)), , drop = FALSE]
            colnames(nc) <- paste0('normcounts_', colnames(nc))
            colnames(rc) <- paste0('rawcounts_', colnames(rc))
            cbind(res, nc, rc)
          }
          }, error = function(e) {
            warning(
              "Skipping cluster ", clust, " for comparison '", lab, "': ",
              conditionMessage(e),
              call. = FALSE
            )
            NULL
          })
          out
        })
        
        names(m_bycluster_crosscondition_de) <- groupinglev_nicelabs
        m_bycluster_crosscondition_de <- m_bycluster_crosscondition_de[lengths(m_bycluster_crosscondition_de) > 0]
        m_bycluster_crosscondition_de_comps[[compidx]] <- m_bycluster_crosscondition_de
      }
    }
    
    names(m_bycluster_crosscondition_de_comps) <- comps$labels
    
  } # end Pseudobulk_mode
  
  
  # ===========================================================================
  # SINGLE-CELL MODE (Seurat FindMarkers / Wilcoxon)
  # Pairwise only: subset to c0 and c1 cells per cluster (not joint modeling).
  # ===========================================================================
  if(Pseudobulk_mode == F){
    
    compslen <- 1:nrow(comps)
    m_bycluster_crosscondition_de_comps <- lapply(compslen, function(compidx){
      
      c0 <- comps$c0[compidx]
      c1 <- comps$c1[compidx]
      lab <- comps$label[compidx]
      
      message('\n', lab)
      
      comp_pseudobulk_md <- sample_metadata[sample_metadata$Condition %in% c(c0, c1),]
      comp_pseudobulk_md$Condition <- factor(comp_pseudobulk_md$Condition, levels = c(c0, c1))
      
      
      
      
      
      
      clusters <- groupinglevs
      names(clusters) <- clusters
      
      m_bycluster_crosscondition_de <- lapply(clusters, function(clust){
        
        
        message(clust)
        
        
        # cells in this cluster for c0 and c1 only
        bigmd <- sobjint@meta.data
        bigmd <- bigmd[bigmd$Condition %in% c(c0, c1),]
        clustmd <- bigmd[bigmd[,grouping_variable] == clust,]
        
        clustmd$Condition <- factor(clustmd$Condition, levels = c(c0, c1))
        cellnums <- table(clustmd$Condition)
        
        if ((cellnums[c1] < 5 | cellnums[c0] < 5)) {
          return()
        }
        
        
        
        #if all good then subset and run
        sobjsub <- sobjint[,rownames(clustmd)]
        
        #rerun sct adjustment?
        # sobjsub <- PrepSCTFindMarkers(sobjsub)
        
        
        #do DE
        
        #turn off parallelization, this step caused memory leak even on tiny datasets
        #future::plan('multisession', workers=workernum)
        
        res <- FindMarkers(sobjsub, logfc.threshold = 0, min.pct = 0,
                           ident.1 = c1, ident.2 = c0,
                           assay = assay, slot = slot,
                           group.by = 'Condition',
                           test.use = DE_test)
        
        # future::plan(strategy = 'sequential')
        
        
        # rm(sobjsub)
        
        #reformat table
        #gene symbol
        res <- cbind(rownames(res), res)
        colnames(res)[1] <- 'gene_symbol'
        rownames(res) <- NULL
        
        #add pct.diff
        res$pct.diff <- res$pct.1 - res$pct.2
        
        
        # #add weight:
        # # -log10 pvalue * sign LFC * abs value of percent difference
        # # ie, significe of DE * sign of DE * num cells exp gene in that direction
        # ## downweight mismatch sign vs pct diff genes... do this by squring them, which makes decimal smaller ##
        # # add + 1 to abs pct diff --> do this to keep FGSEA scores high, otherwise we are actually dividing them by pct diff
        # pctdiff <- res$pct.diff
        # pctdiff[sign(pctdiff) != sign(res$avg_log2FC)] <- pctdiff[sign(pctdiff) != sign(res$avg_log2FC)] ^ 2
        # pctdiff <- abs(pctdiff) + 1
        # res$weight <- -log10(res$p_val) * res$avg_log2FC * pctdiff
        
        ## UPDATE DEC 7 2023, WEIGHT BY -LOG10(PVAL) * SIGN OF LFC
        
        ## prep weighted list ##
        
        #we have to deal with underflow...
        # get -log10 pvalues, sort
        scores <- -log10(res$p_val)
        scores <- scores * sign(res$avg_log2FC)
        names(scores) <- res$gene_symbol
        
        #sort by log pval with names
        scores <- sort(scores,decreasing = T)
        
        
        #also get logFC vector; for INf, we will sort them by LFC...
        logFC_vec <- res$avg_log2FC; names(logFC_vec) <- res$gene_symbol
        
        
        # fix the underflow...
        scores <- scDAPP::fix_underflow(scores, logFC_vec)
        
        
        #make sure scores is in order of genes...
        ### update May 7 2024 --> there was an error in this before
        # it caused scrambling of genes and incorrect pathway analysis :/
        # scores <- scores[match(res$gene_symbol, res$gene_symbol)]
        scores <- scores[match(res$gene_symbol, names(scores))]
        
        #put in res
        res$weight <- scores
        rm(scores, logFC_vec)
        
        #order by weight
        res <- res[order(res$weight, decreasing = T),]
        
        res
        
        
      })
      
      #name them by cluster
      # use the nicelabs defined above
      names( m_bycluster_crosscondition_de ) <- groupinglev_nicelabs
      
      
      #remove empty clusters
      m_bycluster_crosscondition_de <- m_bycluster_crosscondition_de[lengths(m_bycluster_crosscondition_de) > 0]
      
      
      
      
      return(m_bycluster_crosscondition_de)
      
    })
    
    
    
    names(m_bycluster_crosscondition_de_comps) <- comps$labels
    
    
    
  }
  
  
  
  
  
  # Keep the method-native result tables for CSV output. The object below is
  # harmonized in-place afterward for the returned data.frame and filtering.
  m_bycluster_crosscondition_de_comps_natural <- m_bycluster_crosscondition_de_comps

  # ---------------------------------------------------------------------------
  # Harmonize column names to edgeR-style for downstream pathway analysis
  # ---------------------------------------------------------------------------
  if( Pseudobulk_mode == F | DE_test == 'DESeq2' | DE_test == 'DESeq2-LRT' ){
    
    
    
    #for each comparison, loop thru each cluster, and reformat the table
    
    ##  REFORMAT FOR SINGLE CELL WILCOX TEST OUTPUT
    if(Pseudobulk_mode == F){
      
      
      
      m_bycluster_crosscondition_de_comps <- lapply(m_bycluster_crosscondition_de_comps, function(m_bycluster_crosscondition_de){
        
        
        m_bycluster_crosscondition_de <- lapply(m_bycluster_crosscondition_de, function(res){
          
          #reformat all
          res <- res[,c("gene_symbol", "avg_log2FC","p_val", "p_val_adj", "pct.1", "pct.2", "pct.diff", "weight" )]
          
          #rename to match edgeR
          colnames(res) <- c("gene_symbol", "logFC","PValue", "qvalue_bonferroni", "pct.1", "pct.2","pct.diff", "weight"  )
          
          res$FDR <- p.adjust(res$PValue, method = 'fdr')
          
          res
          
        })
        
      })
      
    }
    
    
    
    ##  REFORMAT FOR PSEUDOBULK DESEQ2 / DESEQ2-LRT CELL WILCOX TEST OUTPUT
    if( Pseudobulk_mode == T & (DE_test == 'DESeq2' | DE_test == 'DESeq2-LRT') ){
      
      
      m_bycluster_crosscondition_de_comps <- lapply(m_bycluster_crosscondition_de_comps, function(m_bycluster_crosscondition_de){
        
        
        m_bycluster_crosscondition_de <- lapply(m_bycluster_crosscondition_de, function(res){
          
          #select columns
          # we need to select the columns with raw and norm counts
          counts_colnames <- colnames(res)[(grep(pattern = 'weight', x = colnames(res)) + 1):ncol(res)]
          
          #subset cols
          # res <- res[,c("gene_symbol", "avg_log2FC","p_val", "p_val_adj", "pct.1", "pct.2", "pct.diff", "weight" , counts_colnames)]
          res <- res[,c("gene_symbol", "log2FoldChange","pvalue", "padj", "pct.1", "pct.2", "pct.diff", "weight" , counts_colnames)]
          
          #rename to match edgeR
          colnames(res) <- c("gene_symbol", "logFC","PValue", "qvalue_bonferroni", "pct.1", "pct.2","pct.diff", "weight" , counts_colnames  )
          
          res$FDR <- p.adjust(res$PValue, method = 'fdr')
          
          res
          
        })
        
      })
      
    }
    
    
  }
  
  
  
  
  
  # ---------------------------------------------------------------------------
  # Flatten nested results to one data.frame; filter significant genes
  # ---------------------------------------------------------------------------
  de_results_all_natural <- .flatten_de_results(comps, m_bycluster_crosscondition_de_comps_natural)
  de_results_all <- .flatten_de_results(comps, m_bycluster_crosscondition_de_comps)
  
  de_results_significant <- .filter_significant_de_results(
    de_results_all,
    crossconditionDE_padj_thres,
    crossconditionDE_lfc_thres,
    crossconditionDE_min.pct
  )
  de_sig_idx <- .significant_de_index(
    de_results_all,
    crossconditionDE_padj_thres,
    crossconditionDE_lfc_thres,
    crossconditionDE_min.pct
  )
  de_results_significant_natural <- de_results_all_natural[de_sig_idx, , drop = FALSE]
  
  # ---------------------------------------------------------------------------
  # Optional file output (skipped when outdir_int is missing)
  # ---------------------------------------------------------------------------
  save_results <- !missing(outdir_int) && !is.null(outdir_int) && nzchar(as.character(outdir_int)[1])
  
  if (save_results) {
    out_dir <- file.path(outdir_int, "differentialexpression_crosscondition")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    file_stem <- .de_output_basename(
      DE_test = DE_test,
      grouping_variable = grouping_variable,
      assay = assay,
      slot = slot,
      pseudobulk_mode = Pseudobulk_mode
    )
    write.csv(
      de_results_all_natural,
      file.path(out_dir, paste0(file_stem, "_all.csv")),
      quote = FALSE,
      row.names = FALSE
    )
    write.csv(
      de_results_significant_natural,
      file.path(
        out_dir,
        paste0(
          file_stem,
          .de_thresholds_suffix(
            crossconditionDE_padj_thres,
            crossconditionDE_lfc_thres,
            crossconditionDE_min.pct
          ),
          "_significant.csv"
        )
      ),
      quote = FALSE,
      row.names = FALSE
    )
    .write_numdegs_summary(
      de_results_significant,
      comps,
      groupinglev_nicelabs,
      file.path(out_dir, paste0(file_stem, "_numDEGs_summary.csv"))
    )
  }
  
  return(de_results_all)
  
  
}
# Propeller-specific helpers (sourced into sc_pipeline_modules.R build)

.propeller_design <- function(sample_metadata, design_formula) {
  # Always use fixed-effects-only formula; (1|var) is invalid in model.matrix.
  fixed_formula <- .fixed_effects_formula(design_formula)
  labels <- attr(terms(fixed_formula), "term.labels")
  if (!"Condition" %in% labels) {
    stop("Propeller design requires Condition in formula.", call. = FALSE)
  }
  other <- labels[labels != "Condition"]
  if (length(other)) {
    f <- as.formula(paste("~ 0 + Condition +", paste(other, collapse = " + ")))
  } else {
    f <- as.formula("~ 0 + Condition")
  }
  # coldata needs all variables including random-effect block columns
  md <- .prepare_pseudobulk_coldata(sample_metadata, design_formula)
  rownames(md) <- md$Code
  design <- model.matrix(f, data = md)
  list(design = design, coldata = md)
}

# Build speckle-style propeller result rows: PropMean.<coef> columns, one row
# per cluster. Passing the coefficient matrix (not as.numeric()) lets
# data.frame() expand columns; flattening would recycle p-values onto phantom
# cluster IDs.
.propeller_ttest_result_table <- function(fit.prop, fit.cont, RR) {
  coefs <- fit.prop$coefficients
  if (is.null(dim(coefs))) {
    coefs <- matrix(
      coefs,
      ncol = 1L,
      dimnames = list(names(fit.prop$coefficients), NULL)
    )
  }
  fdr <- p.adjust(fit.cont$p.value[, 1], method = "BH")
  out <- data.frame(
    PropMean = coefs,
    PropRatio = as.numeric(RR),
    Tstatistic = fit.cont$t[, 1],
    P.Value = fit.cont$p.value[, 1],
    FDR = fdr,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (is.null(rownames(out)) || !length(rownames(out))) {
    rn <- rownames(coefs)
    if (!is.null(rn)) rownames(out) <- rn
  }
  out
}

.propeller_prop_ratio <- function(fit.prop, contrasts) {
  n_nz <- sum(contrasts != 0)
  coefs <- fit.prop$coefficients
  if (n_nz == 2L && length(as.numeric(contrasts)) == 2L) {
    z <- apply(coefs, 1, function(x) x^contrasts)
    return(apply(z, 2, prod))
  }
  new.cont <- as.numeric(contrasts[contrasts != 0])
  if (is.null(dim(coefs)) || ncol(coefs) == 1L) {
    return(as.numeric(coefs)^new.cont[1])
  }
  z <- apply(coefs, 1, function(x) x^new.cont)
  apply(z, 2, prod)
}

# Local copy of speckle::propeller.ttest with drop=FALSE so single-coefficient
# contrasts (e.g. interaction terms) do not collapse to a vector.
.propeller_ttest <- function(prop.list, design, contrasts,
                             robust = TRUE, trend = FALSE, sort = TRUE) {
  prop.trans <- prop.list$TransformedProps
  prop <- prop.list$Proportions
  if (nrow(prop.trans) <= 2) {
    message("Setting robust to FALSE for eBayes for less than 3 cell types")
    robust <- FALSE
  }
  fit <- limma::lmFit(prop.trans, design)
  fit.cont <- limma::contrasts.fit(fit, contrasts = contrasts)
  fit.cont <- limma::eBayes(fit.cont, robust = robust, trend = trend)
  n_nz <- sum(contrasts != 0)
  if (n_nz == 2L && length(as.numeric(contrasts)) == 2L) {
    fit.prop <- limma::lmFit(prop, design)
  } else {
    new.des <- design[, contrasts != 0, drop = FALSE]
    fit.prop <- limma::lmFit(prop, new.des)
  }
  RR <- .propeller_prop_ratio(fit.prop, contrasts)
  out <- .propeller_ttest_result_table(fit.prop, fit.cont, RR)
  if (sort) {
    o <- order(out$P.Value)
    return(out[o, , drop = FALSE])
  }
  out
}

# Paired propeller via limma duplicateCorrelation + blocked lmFit (speckle vignette).
.propeller_ttest_blocked <- function(prop.list, design, contrasts, block,
                                     robust = TRUE, trend = FALSE, sort = TRUE) {
  prop.trans <- prop.list$TransformedProps
  prop <- prop.list$Proportions
  if (nrow(prop.trans) <= 2) {
    message("Setting robust to FALSE for eBayes for less than 3 cell types")
    robust <- FALSE
  }
  block <- as.factor(block)
  if (length(block) != ncol(prop.trans)) {
    stop(
      "block length (", length(block), ") must match number of samples in ",
      "transformed proportions (", ncol(prop.trans), ").",
      call. = FALSE
    )
  }
  corfit <- limma::duplicateCorrelation(prop.trans, design, block = block)
  fit <- limma::lmFit(
    prop.trans, design,
    block = block,
    correlation = corfit$consensus
  )
  fit.cont <- limma::contrasts.fit(fit, contrasts = contrasts)
  fit.cont <- limma::eBayes(fit.cont, robust = robust, trend = trend)
  n_nz <- sum(contrasts != 0)
  if (n_nz == 2L && length(as.numeric(contrasts)) == 2L) {
    fit.prop <- limma::lmFit(prop, design)
  } else {
    new.des <- design[, contrasts != 0, drop = FALSE]
    fit.prop <- limma::lmFit(prop, new.des)
  }
  RR <- .propeller_prop_ratio(fit.prop, contrasts)
  out <- .propeller_ttest_result_table(fit.prop, fit.cont, RR)
  if (sort) {
    o <- order(out$P.Value)
    return(out[o, , drop = FALSE])
  }
  out
}

.propeller_contrast_matrix <- function(design, contrast_entry, c1, c0,
                                       style = "edger_expr") {
  parsed <- .parse_contrast(contrast_entry, c0, c1, style = style)
  cn <- colnames(design)
  contr <- matrix(0, nrow = length(cn), ncol = 1)
  rownames(contr) <- cn

  if (identical(style, "deseq2") && !is.null(parsed$deseq_name)) {
    nm <- parsed$deseq_name
    candidates <- unique(c(nm, gsub("\\.", ":", nm, fixed = FALSE),
                           gsub(":", ".", nm, fixed = TRUE)))
    hit <- candidates[candidates %in% cn]
    if (!length(hit)) {
      stop(
        "Could not map DESeq2 name='", nm,
        "' to propeller design columns: ", paste(cn, collapse = ", "),
        call. = FALSE
      )
    }
    colnames(contr) <- hit[1]
    contr[hit[1], 1] <- 1
    return(contr)
  }

  if (identical(style, "deseq2") && !is.null(parsed$deseq_triple)) {
    contrast <- parsed$deseq_triple
    num_col <- paste0(contrast[1], contrast[2])
    denom_col <- paste0(contrast[1], contrast[3])
    if (!(num_col %in% cn) || !(denom_col %in% cn)) {
      stop(
        "Could not map contrast ", num_col, " vs ", denom_col,
        " to propeller design columns: ", paste(cn, collapse = ", "),
        call. = FALSE
      )
    }
    colnames(contr) <- paste0(num_col, "_vs_", denom_col)
    contr[num_col, 1] <- 1
    contr[denom_col, 1] <- -1
    return(contr)
  }

  # EdgeR-style expression
  cv <- .edger_expr_to_vector(design, parsed$edger_expr)
  colnames(contr) <- make.names(parsed$display)
  contr[, 1] <- cv
  contr
}

.propeller_use_simple_wrapper <- function(design_formula, sample_metadata) {
  .is_simple_condition_design(design_formula) &&
    length(unique(sample_metadata$Condition)) == 2L
}

# Descriptive c1 vs c0 mean proportions from the displayed sample-by-cluster table.
.composition_sample_prop_means <- function(comp_proptab, sample_metadata, c1, c0) {
  codes <- colnames(comp_proptab)
  cond <- as.character(sample_metadata$Condition[match(codes, sample_metadata$Code)])
  c1_idx <- which(cond == as.character(c1))
  c0_idx <- which(cond == as.character(c0))
  if (!length(c1_idx) || !length(c0_idx)) {
    stop(
      "Could not match c1/c0 sample columns for composition heatmap annotation.",
      call. = FALSE
    )
  }
  cbind(
    rowMeans(comp_proptab[, c1_idx, drop = FALSE], na.rm = TRUE),
    rowMeans(comp_proptab[, c0_idx, drop = FALSE], na.rm = TRUE)
  )
}

.de_results_cluster_table <- function(de_results, label, cluster) {
  sub <- de_results[de_results$label == label & de_results$cluster == cluster, , drop = FALSE]
  if (!nrow(sub)) return(NULL)
  meta <- c("label", "formula", "c0", "c1", "contrast", "cluster")
  stat_cols <- setdiff(colnames(sub), meta)
  sub[, stat_cols, drop = FALSE]
}

.propeller_format_pres <- function(pres, c1, c0) {
  pres$PropRatio[sign(pres$PropRatio) == -1] <- pres$PropRatio[sign(pres$PropRatio) == -1] * -1
  pres <- pres[order(pres$PropRatio, decreasing = TRUE), , drop = FALSE]
  row_cn <- if ("BaselineProp.clusters" %in% colnames(pres)) {
    "BaselineProp.clusters"
  } else {
    rownames(pres)
  }
  if (is.character(row_cn) && length(row_cn) == 1L && row_cn %in% colnames(pres)) {
    pres$BaselineProp.clusters <- as.character(pres[[row_cn]])
  } else {
    pres$BaselineProp.clusters <- rownames(pres)
  }
  pres
}

#' Normalize comps for cross-condition DE and compositional modules
#'
#' @param comps data.frame with c0, c1, and optional formula, contrast, label
#' @return normalized comps data.frame
#' @export
normalize_comps <- function(comps) {
  .normalize_comps(comps)
}

#' Split flat DE results into a named list of per-cluster tables (one comparison)
#'
#' @param de_results output of `de_across_conditions_module()`
#' @param label comparison label from `comps$label`
#' @return named list of data.frames keyed by cluster
#' @export
de_results_by_cluster <- function(de_results, label) {
  de_sub <- de_results[de_results$label == label, , drop = FALSE]
  clusters <- unique(de_sub$cluster)
  out <- lapply(clusters, function(cl) {
    .de_results_cluster_table(de_sub, label, cl)
  })
  names(out) <- clusters
  out[lengths(out) > 0]
}

# Resolve propeller PropMean columns for c1/c0 (cell-means or speckle group names).
.propeller_propmean_pair_cols <- function(cn, c1, c0) {
  c1_cands <- c(paste0("PropMean.Condition", c1), paste0("PropMean.", c1))
  c0_cands <- c(paste0("PropMean.Condition", c0), paste0("PropMean.", c0))
  c1hit <- c1_cands[c1_cands %in% cn]
  c0hit <- c0_cands[c0_cands %in% cn]
  list(
    c1 = if (length(c1hit)) c1hit[[1]] else NA_character_,
    c0 = if (length(c0hit)) c0hit[[1]] else NA_character_
  )
}

# Strip display asterisks from compositional cluster labels.
.clean_composition_cluster_label <- function(x) {
  trimws(gsub("\\*\\s*", "", as.character(x)))
}

# Harmonize propeller or chisq compres to fixed columns for flat output.
.harmonize_compres_table <- function(compres, c1, c0) {
  df <- as.data.frame(compres, stringsAsFactors = FALSE)
  if ("BaselineProp.clusters" %in% colnames(df)) {
    display <- as.character(df$BaselineProp.clusters)
    pair <- .propeller_propmean_pair_cols(colnames(df), c1, c0)
    data.frame(
      cluster = .clean_composition_cluster_label(display),
      significant = df$P.Value < 0.05,
      PropMean_c1 = if (!is.na(pair$c1)) df[[pair$c1]] else NA_real_,
      PropMean_c0 = if (!is.na(pair$c0)) df[[pair$c0]] else NA_real_,
      PropRatio = df$PropRatio,
      Tstatistic = df$Tstatistic,
      P.Value = df$P.Value,
      FDR = df$FDR,
      stringsAsFactors = FALSE
    )
  } else if (all(c("cluster", "p") %in% colnames(df))) {
    display <- as.character(df$cluster)
    data.frame(
      cluster = .clean_composition_cluster_label(display),
      significant = df$p < 0.05,
      PropMean_c1 = df$c1prop,
      PropMean_c0 = df$c0prop,
      PropRatio = df$asin_ratio,
      Tstatistic = NA_real_,
      P.Value = df$p,
      FDR = df$FDR,
      stringsAsFactors = FALSE
    )
  } else {
    NULL
  }
}

.composition_results_empty <- function() {
  data.frame(
    label = character(), formula = character(), c0 = character(), c1 = character(),
    cluster = character(), significant = logical(),
    PropMean_c1 = numeric(), PropMean_c0 = numeric(), PropRatio = numeric(),
    Tstatistic = numeric(), P.Value = numeric(), FDR = numeric(),
    stringsAsFactors = FALSE
  )
}

.flatten_composition_results <- function(comps, composition_comps_list) {
  pieces <- list()
  for (i in seq_len(nrow(comps))) {
    lab <- comps$label[i]
    entry <- composition_comps_list[[lab]]
    if (is.null(entry) || is.null(entry$compres)) next
    harmonized <- .harmonize_compres_table(entry$compres, comps$c1[i], comps$c0[i])
    if (is.null(harmonized) || !nrow(harmonized)) next
    meta <- data.frame(
      label = lab,
      formula = comps$formula[i],
      c0 = comps$c0[i],
      c1 = comps$c1[i],
      stringsAsFactors = FALSE
    )
    meta_rep <- meta[rep(1L, nrow(harmonized)), , drop = FALSE]
    pieces[[length(pieces) + 1L]] <- cbind(meta_rep, harmonized, stringsAsFactors = FALSE)
  }
  if (!length(pieces)) return(.composition_results_empty())
  out <- dplyr::bind_rows(pieces)
  meta_cols <- c("label", "formula", "c0", "c1", "cluster", "significant",
                 "PropMean_c1", "PropMean_c0", "PropRatio", "Tstatistic", "P.Value", "FDR")
  out[, meta_cols, drop = FALSE]
}

.composition_results_by_comparison_table <- function(composition_results, label) {
  sub <- composition_results[composition_results$label == label, , drop = FALSE]
  if (!nrow(sub)) return(sub)
  stat_cols <- setdiff(colnames(sub), c("label", "formula", "c0", "c1"))
  sub[, stat_cols, drop = FALSE]
}

# Collapse fgsea leadingEdge (list or character vector per row) to one string.
.collapse_leading_edge <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NA_character_)
  genes <- if (is.list(x)) unlist(x, use.names = FALSE) else x
  genes <- as.character(genes)
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (!length(genes)) return(NA_character_)
  paste(genes, collapse = "/")
}

.pathway_results_empty <- function() {
  data.frame(
    label = character(), c0 = character(), c1 = character(),
    cluster = character(), pathway_category = character(),
    pathway = character(), pval = numeric(), padj = numeric(), log2err = numeric(),
    ES = numeric(), NES = numeric(), size = numeric(), leadingEdge = character(),
    stringsAsFactors = FALSE
  )
}

.flatten_pathway_gsea_results <- function(comps, pathway_analysis_mainlist_comps) {
  pieces <- list()
  for (i in seq_len(nrow(comps))) {
    lab <- comps$label[i]
    comp_list <- pathway_analysis_mainlist_comps[[lab]]
    if (is.null(comp_list) || !length(comp_list)) next
    for (pwaycat in names(comp_list)) {
      clust_list <- comp_list[[pwaycat]]
      if (is.null(clust_list) || !length(clust_list)) next
      for (clust in names(clust_list)) {
        entry <- clust_list[[clust]]
        if (is.null(entry)) next
        gseares <- if (is.list(entry) && !is.data.frame(entry)) entry$gseares else entry
        if (is.null(gseares) || !is.data.frame(gseares) || !nrow(gseares)) next
        gseares <- as.data.frame(gseares, stringsAsFactors = FALSE)
        if ("leadingEdge" %in% colnames(gseares)) {
          gseares$leadingEdge <- vapply(gseares$leadingEdge, .collapse_leading_edge, character(1))
        }
        meta <- data.frame(
          label = lab, c0 = comps$c0[i], c1 = comps$c1[i],
          cluster = clust, pathway_category = pwaycat,
          stringsAsFactors = FALSE
        )
        meta_rep <- meta[rep(1L, nrow(gseares)), , drop = FALSE]
        pieces[[length(pieces) + 1L]] <- cbind(meta_rep, gseares, stringsAsFactors = FALSE)
      }
    }
  }
  if (!length(pieces)) return(.pathway_results_empty())
  out <- dplyr::bind_rows(pieces)
  want <- c("label", "c0", "c1", "cluster", "pathway_category",
            "pathway", "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge")
  keep <- intersect(want, colnames(out))
  out[, keep, drop = FALSE]
}

.ora_results_empty <- function() {
  data.frame(
    label = character(), c0 = character(), c1 = character(),
    cluster = character(), pathway_category = character(),
    Direction = character(),
    ID = character(), Description = character(),
    GeneRatio = character(), BgRatio = character(),
    pvalue = numeric(), p.adjust = numeric(), qvalue = numeric(),
    geneID = character(), Count = integer(),
    stringsAsFactors = FALSE
  )
}

.flatten_ora_results <- function(comps, ora_mainlist_comps) {
  pieces <- list()
  for (i in seq_len(nrow(comps))) {
    lab <- comps$label[i]
    comp_list <- ora_mainlist_comps[[lab]]
    if (is.null(comp_list) || !length(comp_list)) next
    for (pwaycat in names(comp_list)) {
      clust_list <- comp_list[[pwaycat]]
      if (is.null(clust_list) || !length(clust_list)) next
      for (clust in names(clust_list)) {
        ora_res <- clust_list[[clust]]
        if (is.null(ora_res) || !is.data.frame(ora_res) || !nrow(ora_res)) next
        ora_res <- as.data.frame(ora_res, stringsAsFactors = FALSE)
        meta <- data.frame(
          label = lab, c0 = comps$c0[i], c1 = comps$c1[i],
          cluster = clust, pathway_category = pwaycat,
          stringsAsFactors = FALSE
        )
        meta_rep <- meta[rep(1L, nrow(ora_res)), , drop = FALSE]
        pieces[[length(pieces) + 1L]] <- cbind(meta_rep, ora_res, stringsAsFactors = FALSE)
      }
    }
  }
  if (!length(pieces)) return(.ora_results_empty())
  dplyr::bind_rows(pieces)
}

#' Filter flat compositional results to one comparison
#'
#' @param composition_results output of `compositional_analysis_module()`
#' @param label comparison label from `comps$label`
#' @return data.frame of compositional statistics for that comparison
#' @export
composition_results_by_comparison <- function(composition_results, label) {
  .composition_results_by_comparison_table(composition_results, label)
}

#' Filter flat pathway GSEA results
#'
#' @param pathway_results output of `pathwayanalysis_crosscondition_module()`
#' @param label comparison label
#' @param cluster optional cluster filter
#' @param pathway_category optional MSigDB category filter
#' @return filtered data.frame
#' @export
pathway_results_by_cluster <- function(pathway_results, label, cluster = NULL, pathway_category = NULL) {
  sub <- pathway_results[pathway_results$label == label, , drop = FALSE]
  if (!is.null(cluster)) sub <- sub[sub$cluster == cluster, , drop = FALSE]
  if (!is.null(pathway_category)) sub <- sub[sub$pathway_category == pathway_category, , drop = FALSE]
  sub
}

#' Filter flat ORA results
#'
#' @param ora_results output of `ORA_crosscondition_module()`
#' @param label comparison label
#' @param cluster optional cluster filter
#' @param pathway_category optional MSigDB category filter
#' @return filtered data.frame
#' @export
ora_results_by_cluster <- function(ora_results, label, cluster = NULL, pathway_category = NULL) {
  sub <- ora_results[ora_results$label == label, , drop = FALSE]
  if (!is.null(cluster)) sub <- sub[sub$cluster == cluster, , drop = FALSE]
  if (!is.null(pathway_category)) sub <- sub[sub$pathway_category == pathway_category, , drop = FALSE]
  sub
}

.ora_marker_msigdb_celltype_results_empty <- function() {
  data.frame(
    context = character(),
    cluster = character(),
    pathway_category = character(),
    ID = character(), Description = character(),
    GeneRatio = character(), BgRatio = character(),
    pvalue = numeric(), p.adjust = numeric(), qvalue = numeric(),
    geneID = character(), Count = integer(),
    stringsAsFactors = FALSE
  )
}

.select_cluster_marker_genes <- function(marker_results,
                                         cluster,
                                         marker_padj_thres,
                                         top_markers_per_cluster,
                                         min_genes) {
  req <- c("cluster", "gene", "p_val_adj", "score")
  miss <- setdiff(req, colnames(marker_results))
  if (length(miss)) {
    stop(
      "marker_results missing columns: ", paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  sub <- marker_results[marker_results$cluster == cluster, , drop = FALSE]
  sub <- sub[!is.na(sub$p_val_adj) & sub$p_val_adj < marker_padj_thres, , drop = FALSE]
  if (!nrow(sub)) return(character())
  sub <- sub[order(sub$score, decreasing = TRUE), , drop = FALSE]
  sub <- utils::head(sub, top_markers_per_cluster)
  genes <- unique(sub$gene)
  if (length(genes) < min_genes) return(character())
  genes
}

.top_pathways_per_cluster_msigdb_celltype <- function(ora_results, dotplot_top_n = 5L) {
  empty <- data.frame(
    context = character(),
    cluster = character(),
    pathway_category = character(),
    ID = character(),
    Description = character(),
    GeneRatio = character(),
    BgRatio = character(),
    pvalue = numeric(),
    p.adjust = numeric(),
    qvalue = numeric(),
    geneID = character(),
    Count = integer(),
    stringsAsFactors = FALSE
  )
  if (is.null(ora_results) || !is.data.frame(ora_results) || !nrow(ora_results)) {
    return(empty)
  }
  plot_df <- ora_results
  plot_df$pathway <- plot_df$Description
  if ("ID" %in% colnames(plot_df)) {
    empty_desc <- !nzchar(plot_df$pathway)
    plot_df$pathway[empty_desc] <- plot_df$ID[empty_desc]
  }
  plot_df$neglog10padj <- -log10(pmax(plot_df$p.adjust, .Machine$double.xmin))
  plot_df %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(.data$neglog10padj, n = dotplot_top_n, with_ties = FALSE) %>%
    dplyr::slice_max(.data$Count, n = dotplot_top_n, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    as.data.frame()
}

.save_celltype_marker_prediction_outputs <- function(ora_results,
                                                     context_label,
                                                     outdir,
                                                     dotplot_top_n = 5L,
                                                     cp.font.size = 5) {
  out_subdir <- file.path(outdir, "celltype_marker_prediction", context_label)
  dir.create(out_subdir, recursive = TRUE, showWarnings = FALSE)

  write.csv(
    ora_results,
    file.path(out_subdir, "ora_results.csv"),
    quote = FALSE,
    row.names = FALSE
  )

  top_pathways <- .top_pathways_per_cluster_msigdb_celltype(
    ora_results = ora_results,
    dotplot_top_n = dotplot_top_n
  )
  write.csv(
    top_pathways,
    file.path(out_subdir, "top_pathways_per_cluster.csv"),
    quote = FALSE,
    row.names = FALSE
  )

  dotplot <- .build_msigdb_celltype_ora_dotplot(
    ora_results = ora_results,
    context_label = context_label,
    dotplot_top_n = dotplot_top_n,
    cp.font.size = cp.font.size
  )
  if (!is.null(dotplot)) {
    pdf(
      file.path(out_subdir, "celltype_marker_prediction_dotplot.pdf"),
      width = 7,
      height = 7
    )
    print(dotplot)
    while (!is.null(grDevices::dev.list())) grDevices::dev.off()
  }

  list(ora_results = ora_results, dotplot = dotplot, context_label = context_label)
}

.build_msigdb_celltype_ora_dotplot <- function(ora_results,
                                  context_label,
                                  dotplot_top_n = 5,
                                  cp.font.size = 5) {
  if (is.null(ora_results) || !is.data.frame(ora_results) || !nrow(ora_results)) {
    return(NULL)
  }
  plot_df <- .top_pathways_per_cluster_msigdb_celltype(
    ora_results = ora_results,
    dotplot_top_n = dotplot_top_n
  )
  if (!nrow(plot_df)) return(NULL)

  top_pathway_labels <- plot_df$Description
  if ("ID" %in% colnames(plot_df)) {
    empty_desc <- !nzchar(top_pathway_labels)
    top_pathway_labels[empty_desc] <- plot_df$ID[empty_desc]
  }
  top_pathway_labels <- unique(top_pathway_labels)

  plot_df <- ora_results
  plot_df$pathway <- plot_df$Description
  if ("ID" %in% colnames(plot_df)) {
    empty_desc <- !nzchar(plot_df$pathway)
    plot_df$pathway[empty_desc] <- plot_df$ID[empty_desc]
  }
  plot_df <- plot_df[plot_df$pathway %in% top_pathway_labels, , drop = FALSE]
  plot_df$neglog10padj <- -log10(pmax(plot_df$p.adjust, .Machine$double.xmin))

  if (!nrow(plot_df)) return(NULL)

  plot_df$cluster <- factor(plot_df$cluster, levels = unique(as.character(plot_df$cluster)))
  plot_df$pathway <- factor(plot_df$pathway, levels = rev(unique(as.character(plot_df$pathway))))

  ggplot2::ggplot(plot_df, ggplot2::aes(x = cluster, y = pathway, size = Count, col = neglog10padj)) +
    ggplot2::geom_point() +
    ggplot2::theme_linedraw() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = cp.font.size),
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)
    ) +
    ggplot2::scale_color_gradient(low = "grey80", high = "#B2182B", name = "-log10(padj)") +
    ggplot2::scale_size(range = c(2, 6), name = "Gene count") +
    ggplot2::xlab("Cluster") +
    ggplot2::ylab("") +
    ggplot2::ggtitle("MSigDB cell-type signatures", subtitle = context_label)
}

#' ORA of cluster markers against MSigDB cell-type signatures
#'
#' Runs `clusterProfiler::enricher()` on cluster marker genes (top N by score among
#' markers passing `marker_padj_thres`) against MSigDB cell-type gene sets from a prepared
#' msigdbr table. Saves CSV tables and a summary dotplot with inclusive
#' top-pathway filtering per cluster (same idea as GSEA summary dotplots).
#'
#' @param marker_results data.frame from `Seurat::FindAllMarkers()` with `score` column.
#' @param pathways output of `preppathways_pathwayanalysis_crosscondition_module()`.
#' @param context_label string label for outputs (e.g. sample Code or `"integrated"`).
#' @param outdir directory; writes to `{outdir}/celltype_marker_prediction/{context_label}/`.
#' @param marker_padj_thres adjusted p-value cutoff for marker genes (default 0.05).
#' @param top_markers_per_cluster max markers per cluster after padj filter (default 100).
#' @param min_genes minimum genes required to run ORA per cluster (default 7).
#' @param pathway_padj_thres q-value cutoff passed to `enricher()` (default 0.1).
#' @param dotplot_top_n top pathways per cluster for dotplot when many are significant (default 5).
#' @param workernum number of parallel workers (default 1).
#' @param cp.font.size axis text size for dotplot (default 5).
#'
#' @return List with `ora_results` (data.frame), `dotplot` (ggplot or NULL), and `context_label`.
#' @export
ORA_cluster_markers_msigdb_celltype_module <- function(marker_results,
                                          pathways,
                                          context_label,
                                          outdir,
                                          marker_padj_thres = 0.05,
                                          top_markers_per_cluster = 100L,
                                          min_genes = 7L,
                                          pathway_padj_thres = 0.1,
                                          dotplot_top_n = 5L,
                                          workernum = 1L,
                                          cp.font.size = 5) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Install clusterProfiler to run MSigDB cell-type ORA.", call. = FALSE)
  }
  if (!is.data.frame(marker_results)) {
    stop("marker_results must be a data.frame.", call. = FALSE)
  }
  if (!is.data.frame(pathways)) {
    stop("pathways must be a data.frame.", call. = FALSE)
  }

  celltype_subcat <- .msigdb_celltype_gs_subcat()
  term2gene <- pathways[pathways$gs_subcat == celltype_subcat, c("gs_name", "gene_symbol"), drop = FALSE]
  if (!nrow(term2gene)) {
    warning("No MSigDB cell-type gene sets found in pathways table.", call. = FALSE)
    empty <- .ora_marker_msigdb_celltype_results_empty()
    return(.save_celltype_marker_prediction_outputs(
      ora_results = empty,
      context_label = context_label,
      outdir = outdir,
      dotplot_top_n = dotplot_top_n,
      cp.font.size = cp.font.size
    ))
  }

  if ("gene" %in% colnames(marker_results)) {
    marker_genes <- unique(as.character(marker_results$gene))
    marker_genes <- marker_genes[!is.na(marker_genes) & nzchar(marker_genes)]
    pathway_genes <- unique(as.character(term2gene$gene_symbol))
    overlap_n <- length(intersect(marker_genes, pathway_genes))
    if (overlap_n < min_genes) {
      warning(
        "Low overlap (", overlap_n, " genes) between cluster marker genes and MSigDB ",
        celltype_subcat, " gene symbols for context ", context_label,
        ". Check that marker genes use the same identifier type as msigdbr (e.g. symbols, not Ensembl IDs).",
        call. = FALSE
      )
    }
  }

  clusters <- unique(as.character(marker_results$cluster))
  clusters <- clusters[!is.na(clusters)]

  if (workernum > 1L) {
    cl <- parallel::makeCluster(workernum, rscript_args = c("--no-init-file", "--no-site-file", "--no-environ"))
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    clust_res <- foreach::foreach(
      clust = clusters,
      .packages = c("clusterProfiler"),
      .export = c(
        ".select_cluster_marker_genes", "marker_results", "term2gene",
        "marker_padj_thres", "top_markers_per_cluster", "min_genes",
        "pathway_padj_thres", "context_label", "celltype_subcat"
      ),
      .noexport = c("pathways"),
      .verbose = FALSE
    ) %dopar% {
      genenames <- .select_cluster_marker_genes(
        marker_results = marker_results,
        cluster = clust,
        marker_padj_thres = marker_padj_thres,
        top_markers_per_cluster = top_markers_per_cluster,
        min_genes = min_genes
      )
      if (!length(genenames)) return(NULL)

      ora_res <- tryCatch(
        clusterProfiler::enricher(
          genenames,
          TERM2GENE = term2gene,
          qvalueCutoff = pathway_padj_thres,
          pvalueCutoff = 1
        ),
        error = function(e) NULL
      )
      if (is.null(ora_res) || !length(ora_res)) return(NULL)
      ora_res <- as.data.frame(ora_res, stringsAsFactors = FALSE)
      if (!nrow(ora_res)) return(NULL)
      ora_res$cluster <- clust
      ora_res$context <- context_label
      ora_res$pathway_category <- celltype_subcat
      ora_res
    }
    names(clust_res) <- clusters
  } else {
    clust_res <- stats::setNames(lapply(clusters, function(clust) {
      genenames <- .select_cluster_marker_genes(
        marker_results = marker_results,
        cluster = clust,
        marker_padj_thres = marker_padj_thres,
        top_markers_per_cluster = top_markers_per_cluster,
        min_genes = min_genes
      )
      if (!length(genenames)) return(NULL)

      ora_res <- tryCatch(
        clusterProfiler::enricher(
          genenames,
          TERM2GENE = term2gene,
          qvalueCutoff = pathway_padj_thres,
          pvalueCutoff = 1
        ),
        error = function(e) NULL
      )
      if (is.null(ora_res) || !length(ora_res)) return(NULL)
      ora_res <- as.data.frame(ora_res, stringsAsFactors = FALSE)
      if (!nrow(ora_res)) return(NULL)
      ora_res$cluster <- clust
      ora_res$context <- context_label
      ora_res$pathway_category <- celltype_subcat
      ora_res
    }), clusters)
  }

  clust_res <- clust_res[lengths(clust_res) > 0L]
  if (!length(clust_res)) {
    empty <- .ora_marker_msigdb_celltype_results_empty()
    return(.save_celltype_marker_prediction_outputs(
      ora_results = empty,
      context_label = context_label,
      outdir = outdir,
      dotplot_top_n = dotplot_top_n,
      cp.font.size = cp.font.size
    ))
  }

  ora_results <- dplyr::bind_rows(clust_res)
  want <- c(
    "context", "cluster", "pathway_category", "ID", "Description",
    "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"
  )
  keep <- intersect(want, colnames(ora_results))
  ora_results <- ora_results[, keep, drop = FALSE]

  .save_celltype_marker_prediction_outputs(
    ora_results = ora_results,
    context_label = context_label,
    outdir = outdir,
    dotplot_top_n = dotplot_top_n,
    cp.font.size = cp.font.size
  )
}

#' Batch MSigDB cell-type ORA on named per-sample cluster marker tables
#'
#' @param marker_results_list named list of FindAllMarkers data.frames.
#' @param pathways output of `preppathways_pathwayanalysis_crosscondition_module()`.
#' @param outdir base output directory (typically `outdir_indi`).
#' @param ... passed to `ORA_cluster_markers_msigdb_celltype_module()`.
#'
#' @return Named list of results from `ORA_cluster_markers_msigdb_celltype_module()` per sample.
#' @export
ORA_cluster_markers_msigdb_celltype_batch_module <- function(marker_results_list,
                                                pathways,
                                                outdir,
                                                ...) {
  if (!is.list(marker_results_list) || is.null(names(marker_results_list))) {
    stop("marker_results_list must be a named list.", call. = FALSE)
  }
  res_list <- lapply(names(marker_results_list), function(code) {
    ORA_cluster_markers_msigdb_celltype_module(
      marker_results = marker_results_list[[code]],
      pathways = pathways,
      context_label = code,
      outdir = outdir,
      ...
    )
  })
  stats::setNames(res_list, names(marker_results_list))
}







#' Prep MSIGDB pathways for pathway analysis
#'
#' This is a modular component of the scRNAseq analysis pipeline. Prep pathways from MSIGDB via the msigdbr package. We include the Hallmarks category; Gene Ontology BP, MF and CC; Reactome; KEGG; transcription factor CHIP-seq targets in the Gene Transcription Regulation Database (TFT_GTRD); inferred transcription factor targets via motif analysis from Xie et al Nature 2005 (TFT_Legacy); and MSigDB cell-type signatures (cached for cluster-marker ORA, not used in cross-condition GSEA/ORA). Non-cell-type gene sets with < 500 genes are included; cell-type sets may contain up to 1000 genes.
#'
#' Prepared pathways are cached on disk (see [scDAPP::resolve_msigdbr_cache_dir()]) and are not written to `outdir_int` by default.
#'
#' @param species string, species such as "Homo sapiens" or "Mus musculus"
#' @param outdir_int string, pipeline integration output directory; used as fallback cache location when user/XDG/R cache dirs are not writable.
#' @param msigdbr_cache_dir optional string, directory for cached msigdbr tables (typically .../scDAPP). When NULL, uses `XDG_CACHE_HOME/scDAPP` if set, else [tools::R_user_dir()] cache.
#' @param pwaycats optional character vector of normalized `gs_subcat` values to retain. Defaults to Hallmark, GO, Reactome, KEGG, and TFT categories used by the pipeline.
#' @param refresh_msigdbr_cache logical; if TRUE, ignore existing cache and re-download from msigdbr.
#'
#' @return a data.frame similar to the output of `msigdbr::msigdbr`, but filtering for some specific categories / subcategories.
#' @export
#'
#' @examples
#' \dontrun{
#' pathways <- preppathways_pathwayanalysis_crosscondition_module(
#' species = 'Mus musculus',
#' outdir_int = 'path/to/directory')
#' }
preppathways_pathwayanalysis_crosscondition_module <- function(species,
                                                               outdir_int,
                                                               msigdbr_cache_dir = NULL,
                                                               pwaycats = NULL,
                                                               refresh_msigdbr_cache = FALSE)
{
  if (is.null(pwaycats)) {
    pwaycats <- .default_msigdbr_pwaycats()
  } else {
    pwaycats <- gsub(":", "_", pwaycats)
    names(pwaycats) <- pwaycats
  }

  if (!refresh_msigdbr_cache) {
    pathways <- .load_msigdbr_cache(
      species = species,
      pwaycats = pwaycats,
      cache_dir = msigdbr_cache_dir,
      fallback_dir = outdir_int
    )
    if (!is.null(pathways)) {
      invisible(gc(full = TRUE, reset = FALSE, verbose = FALSE))
      return(pathways)
    }
  }

  message("Accessing MSIGDBR database")
  pathways <- msigdbr::msigdbr(species = species)
  .validate_msigdbr_raw(pathways)
  pathways <- .normalize_msigdbr_pathways(pathways, pwaycats = pwaycats)
  .validate_prepared_pathways(pathways, pwaycats = pwaycats, species = species)

  .save_msigdbr_cache(
    pathways = pathways,
    species = species,
    pwaycats = pwaycats,
    cache_dir = msigdbr_cache_dir,
    fallback_dir = outdir_int
  )

  invisible(gc(full = TRUE, reset = FALSE, verbose = FALSE))
  pathways
}





#' GSEA analysis for cross condition DE analysis
#'
#' This function is a modular component of the scRNAseq pipeline. Perform GSEA analysis via the FGSEA package on the results of differential expression (DE) analysis for cross-condition comparison. Multiple conditions are supported. Plots and tables are saved.
#'
#' @param de_results harmonized data.frame from `scDAPP::de_across_conditions_module()`.
#' @param pathways data.frame, the output of `scDAPP::preppathways_pathwayanalysis_crosscondition_module()`
#' @param sample_metadata data.frame with sample names and conditions, same as in `scDAPP::de_across_conditions_module()`, see that function's documentation for description.
#' @param comps data.frame with conditions to test in GSEA, same as in `scDAPP::de_across_conditions_module()`, see that function for description
#' @param pathway_padj_thres numeric, threshold for significance of pathway enrichment after multiple test correction
#' @param pwaycats UNTESTED CURRENTLY. character vector of msigdb pathways data.frame in gs_subcat to run.
#' @param workernum integer. number of CPUs. default = 1.
#' @param outdir_int string, directory to save pathways to. Will create a sub-directory called "pathwayanalysis_crosscondition" and save inside of there.
#' @param cp.font.size numeric. size of pathway name test in summary plots. larger than 5 will likely result in name overlap and unreadable plots, currently not easy to solve
#'
#' @return A list with `pathway_results` (flat data.frame), `pathwaysummplots_comps`
#'   (summary plots per comparison), and `pathway_cluster_plots` (nested gseares/dotplot
#'   objects per cluster for HTML reporting; not a tabular API).
#' @export
#'
#' @examples
#' \dontrun{
#'
#'
#' # FIRST: run `scDAPP::de_across_conditions_module()`. the output object of that is used as the main input for this pathway analysis function.
#'
#' # SECOND: prep pathways before running
#' pathways <- scDAPP::preppathways_pathwayanalysis_crosscondition_module(species = species,
#' outdir_int = outdir_int)
#'
#' # THIRD: run the pathway analysis.
#' # see `scDAPP::de_across_conditions_module()` for a description of the sample_metadata and comps files.
#' pways_output_list <- pathwayanalysis_crosscondition_module(
#' de_results = de_results,
#' pathways = pathways,
#' sample_metadata = sample_metadata,
#' comps = comps,
#' workernum = 6,
#' outdir_int = outdir_int
#' )
#'
#' # Access the output files (also saved under outdir_int):
#' pathway_results <- pways_output_list$pathway_results
#' pathwaysummplots_comps <- pways_output_list$pathwaysummplots_comps
#'
#' }
pathwayanalysis_crosscondition_module <- function(de_results,
                                                  pathways,
                                                  sample_metadata,
                                                  comps,
                                                  pathway_padj_thres,
                                                  pwaycats,
                                                  workernum,
                                                  outdir_int,
                                                  cp.font.size
){
  
  
  require(fgsea)
  require(foreach)
  require(doParallel)
  require(parallel)
  
  if (!is.data.frame(de_results)) {
    stop("de_results must be a data.frame from de_across_conditions_module().", call. = FALSE)
  }
  comps <- .normalize_comps(comps)
  
  
  #   UPDATE DECEMBER 7 2023 deg.weight has been deprecated, we will stick with -log10pval * sign FC, but weight value can be modified for res before this
  #   if( missing(deg.weight) ){deg.weight <- 'pval'}
  #   if( !(deg.weight %in% c('auto', 'pval') ) ){
  #     stop('deg.weight must be either "auto" or "pval" ')
  #   }
  
  if( missing(cp.font.size) ) {
    
    #set font size
    # this is the best size, any bigger there will be overlap...
    cp.font.size <- 5
  }
  
  
  
  
  #### picking default categories
  # because hallmark is a "category" and rest are "subcategories", it is hard to make this automated
  # guess it may be possible if we set missing subcat as cat...
  # table( pathways[pathways$gs_subcat=='',"gs_cat"] )
  # for now hardcode these
  if(missing(pwaycats)){
    pwaycats <- c("HALLMARK", "GO_BP", "GO_MF", "GO_CC", "CP_REACTOME", "CP_KEGG", "TFT_GTRD", "TFT_TFT_Legacy")
  }
  
  if(missing(pathway_padj_thres)){
    pathway_padj_thres <- 0.1
  }
  if(missing(workernum)){
    workernum <- 1
  }
  
  
  ### pathway analysis for each condition comparison
  
  # for each condition comparison,
  # for each cluster
  # do pathway analysis and make plot
  
  
  
  
  
  compslen <- seq_len(nrow(comps))
  
  pathway_analysis_mainlist_comps <- lapply(compslen, function(compidx){
    
    c1 <- comps$c1[compidx]
    c0 <- comps$c0[compidx]
    lab <- comps$label[compidx]
    
    message(lab)
    
    de_sub <- de_results[de_results$label == lab, , drop = FALSE]
    clusters <- unique(de_sub$cluster)
    
    pwayoutdir <- paste0(outdir_int, '/pathwayanalysis_crosscondition/', lab, '/')
    if( !dir.exists(pwayoutdir) ){ dir.create(pwayoutdir, recursive = T) }
    
    
    ### loop thru pathway categories
    names(pwaycats) <- pwaycats
    
    #set gene universe
    pwaycat <- pwaycats[1] #for testing
    
    pathway_analysis_mainlist <- lapply(pwaycats, function(pwaycat){
      
      message('\n\n', pwaycat, '\n\n')
      
      
      #get pways and genes in this category
      term2gene <- pathways[pathways$gs_subcat == pwaycat,c('gs_name', 'gene_symbol')]
      
      
      #pways as list for gsea
      pwayl = split(term2gene$gene_symbol, term2gene$gs_name)
      rm(term2gene)
      
      
      #get list of pathways upreg in each cluster
      
      cl <- parallel::makeCluster(workernum, rscript_args = c("--no-init-file", "--no-site-file", "--no-environ"))
      doParallel::registerDoParallel(cl)
      
      
      
      #pwayres_DE_across_conditions_per_cluster <- lapply(clusters, function(clust){
      pwayres_DE_across_conditions_per_cluster <- foreach(clust = clusters,
                                                          .packages = c('fgsea', 'ggplot2'),
                                                          .export = c('de_sub', 'pathway_padj_thres', 'cp.font.size', '.de_results_cluster_table'),
                                                          .noexport = c('pathways'),
                                                          .verbose = T) %dopar%
        {
          
          
          
          invisible(gc(full = T, reset = F, verbose = F))
          
          
          res <- .de_results_cluster_table(de_sub, de_sub$label[1], clust)
          if (is.null(res) || !nrow(res)) return(NULL)
          
          
          
          
          # UPDATE DEC 7 2023 WE DEPRECATED DEG.WEIGHT
          # AND MADE DEFAULT WEIGHT FROM DE MODULE AS -LOG10(PVAL ) * SIGN(L2FC)
          # if(deg.weight == 'pval'){
          #
          #   ## prep weighted list ##
          #
          #   #we have to deal with underflow...
          #   # get -log10 pvalues, sort
          #   scores <- -log10(res$PValue)
          #   scores <- scores * sign(res$logFC)
          #   names(scores) <- res$gene_symbol
          #
          #   #sort by log pval with names
          #   scores <- sort(scores,decreasing = T)
          #
          #
          #   #also get logFC vector; for INf, we will sort them by LFC...
          #   logFC_vec <- res$logFC; names(logFC_vec) <- res$gene_symbol
          #
          #
          #   # fix the underflow...
          #   scores <- fix_underflow(scores, logFC_vec)
          #
          #
          #   #put in res
          #   res$weight <- scores
          #   rm(scores, logFC_vec)
          #
          #
          # }
          
          #input for gsea is weight named by gene
          res <- res[order(res$weight, decreasing = T),]
          gl <- res$weight; names(gl) <- res$gene_symbol
          
          #clean env
          rm(res)
          
          ## run GSEA ##
          # first run multilevel, then try npermsimple = 1000
          gseares <- fgsea::fgsea(pathways=pwayl, stats=gl, nproc = 1)
          
          invisible(gc(full = T, reset = F, verbose = F))
          
          
          #sometimes there are NAs due to "severely unbalanced pathways", try to fix
          if( any(is.na(gseares$NES)) ){
            rm(gseares)
            
            gseares <- fgsea::fgsea(pathways=pwayl, stats=gl, nPermSimple=10000, nproc = 1)
            
            invisible(gc(full = T, reset = F, verbose = F))
            
          }
          
          #clean env
          rm(gl)
          
          
          #format as data.frame instead of data.table
          gseares <- as.data.frame(gseares)
          
          #ensure no NAs are kept
          gseares <- gseares[complete.cases(gseares[,1:7]),,drop=F]
          
          #order by NES
          gseares <- gseares[order(gseares$NES, decreasing = T),]
          
          #apply cutoff of pathway_padj_thres
          gseares <- gseares[gseares$padj < pathway_padj_thres, ,drop=F]
          
          #select pathways with more than just 1 gene in the list
          gseares <- gseares[gseares$size > 2,,drop=F]
          
          #skip if no significant results
          if( nrow(gseares)==0){ return() }
          
          
          #prep for plot, leave out leading edge
          gseares_plot <- gseares[,-8]
          
          #if more than 20 ,select just 20
          gseares_plot <- rbind( head( gseares_plot[gseares_plot$NES>0,,drop=F], 10) ,
                                 tail( gseares_plot[gseares_plot$NES<0,,drop=F], 10) )
          
          #make pathway names more readable by using spaces instead of underscores
          gseares_plot$pathway <- gsub(gseares_plot$pathway, pattern = '_', replacement = ' ')
          
          #make pathway names more readable by splitting long ones to multiple lines
          gseares_plot$pathway <- stringr::str_wrap(gseares_plot$pathway, width = 35)
          
          #make sure order is by -log(padj) * NES
          gseares_plot$weight <- -log(gseares_plot$padj) * sign(gseares_plot$NES)
          gseares_plot <- gseares_plot[order(gseares_plot$weight, decreasing = T),]
          gseares_plot$pathway <- factor(gseares_plot$pathway, levels = rev(gseares_plot$pathway)  )
          
          
          #plot it
          
          #fix color issue when just 1 obs
          if(nrow(gseares_plot) == 1){
            dp_single_col <- sign(gseares_plot$NES)
            dp_single_col <- ifelse(dp_single_col==1, yes = 'red', no = 'steelblue')
            
            
            dp <- ggplot(gseares_plot, aes(-log10(padj), pathway, col=NES, size = size))+
              geom_point()+
              theme_linedraw()+
              theme(axis.text=element_text(size=cp.font.size) )+
              scale_color_gradientn(colors = dp_single_col) +
              scale_size(range=c(2,6))
            
          } else{
            
            
            dp <- ggplot(gseares_plot, aes(-log10(padj), pathway, col=NES, size = size))+
              geom_point()+
              theme_linedraw()+
              theme(axis.text=element_text(size=cp.font.size) )+
              scale_color_gradient2(low = 'steelblue', high = 'red', mid = 'white', midpoint = 0, name = 'Normalized\nEnrichment\nScore')+
              scale_size(range=c(2,6))
            
          }
          
          #return the result table and the plot
          return(list(gseares = gseares, dp = dp))
          
          
          
          
        } #per-cluster loop for this pathway category loop end
      
      parallel::stopCluster(cl)
      
      
      names(pwayres_DE_across_conditions_per_cluster) <- clusters
      
      
      return(pwayres_DE_across_conditions_per_cluster)
      
    }) #per-pathway category loop end
    
    
    
    
    ### remove all NULLS (clusters with no pathways)
    
    # remove null categories, ie entire category had no significant pathways
    
    #recursively set all missing to 0
    pathway_analysis_mainlist = lapply(pathway_analysis_mainlist, function(pwayres_cats){
      
      pwayres_cats <- lapply(pwayres_cats, function(pwayres_clusts){
        pwayres_clusts[lengths(pwayres_clusts) > 0]
      })
      
      pwayres_cats[lengths(pwayres_cats) > 0]
      
      
    })
    
    #remove any missing categories
    pathway_analysis_mainlist <- pathway_analysis_mainlist[lengths(pathway_analysis_mainlist)>0]
    
    
    invisible(gc(full = T, reset = F, verbose = F))
    
    
    #save pway analysis for this comparison
    
    # loop over this new pwayruns, since some categories theoretically don't have any enriched though unlikely
    pwayruns <- names(pathway_analysis_mainlist)
    
    
    #for each category, get clster res in that cateogry,
    # for each cluster, save the up/down csv and plots
    invisible(
      finalpwayouts <- lapply(pwayruns, function(pwaycat){
        
        
        
        # message(pwaycat)
        
        
        subcatout <- paste0(pwayoutdir, '/', pwaycat, '/')
        
        # dir.create(subcatout) --> do this with recursive later, maybe prevent even making it if all don't work
        
        clustres <- pathway_analysis_mainlist[[pwaycat]]
        
        clusters <- names(clustres)
        
        
        #for each cluster, get up/down csv, up/dwon plot, and save
        numpways <- lapply(clusters, function(clust){
          
          
          
          
          pwayres_DE_across_conditions_per_cluster <- clustres[[clust]]
          
          
          
          if(is.null(pwayres_DE_across_conditions_per_cluster)){return()}
          
          gseares <- pwayres_DE_across_conditions_per_cluster$gseares
          dp <- pwayres_DE_across_conditions_per_cluster$dp
          
          #save PDFs and CSVs
          
          subcatout_clustdir <- paste0(subcatout, '/', clust, '/')
          
          suppressWarnings(dir.create(subcatout_clustdir, recursive = T))
          
          
          #upcsv
          subcatout_clustdir_gseares <- paste0(subcatout_clustdir, '/pathwaytable.csv')
          
          #gseares, leading edge needs to be adjusted...
          gseares$leadingEdge <- sapply(gseares$leadingEdge, function(x){ paste(x, collapse = '/') })
          
          write.csv(gseares, subcatout_clustdir_gseares, quote = F, row.names = F)
          
          
          subcatout_clustdir_dp <- paste0(subcatout_clustdir, '/dotplot_toppathways.pdf')
          
          pdf(subcatout_clustdir_dp)
          print( dp )
          
          while (!is.null(dev.list()))  dev.off()
          
          
          
          nrow(gseares)
          
          
          
          
          
          
          
        } ) # close clusters lapply
        
        
      }) # close saving loop for all categories
      
    ) # close invisible wrap around lapply
    
    
    
    
    
    
    
    
    
    
    ### close any open devices
    while (!is.null(dev.list()))  dev.off()
    
    
    
    
    return(pathway_analysis_mainlist)
    
    
  }) # close cross condition lapply
  
  
  names(pathway_analysis_mainlist_comps) <- comps$labels
  
  
  
  
  #remove big objects
  rm(pathways)
  invisible(gc(full = T, reset = F, verbose = F))
  
  
  
  
  ### prep summary plots for each category
  
  compslen <- 1:nrow(comps)
  pathwaysummplots_comps <- lapply(compslen, function(compidx){
    
    #get pway analysis
    pathway_analysis_mainlist <- pathway_analysis_mainlist_comps[[compidx]]
    
    c1 <- comps$c1[compidx]
    c0 <- comps$c0[compidx]
    lab <- comps$label[compidx]
    
    
    ### extract the table from all categories
    
    cat_cpres_list <- lapply(pathway_analysis_mainlist, function(pwaycatlist){
      
      #in each cluster:
      # get the tables from each up/dn
      
      # use cluster index, we need the cluster name
      
      clust_cpres <- lapply(seq_along(pwaycatlist), function(clustidx){
        
        clustname <- names(pwaycatlist)[clustidx]
        pwayres_DE_across_conditions_per_cluster <- pwaycatlist[[clustidx]]
        
        #use the dotplot data for the table
        gseares_plot <- pwayres_DE_across_conditions_per_cluster$dp$data
        
        gseares_plot$cluster = clustname
        gseares_plot$condition = c1
        gseares_plot[sign(gseares_plot$NES) == -1, "condition"] = c0
        
        return(gseares_plot)
        
        
      })
      
      dplyr::bind_rows(clust_cpres)
      
    })
    
    
    #loop thru each category's result data.frame, splitting c1 and c0, and plotting
    
    # Use seq_along() so empty categories yield an empty result instead of 1:0 indexing.
    summplots_cats <- lapply(seq_along(cat_cpres_list), function(catdex){
      
      cpres_cat <- cat_cpres_list[[catdex]]
      catname <- names(cat_cpres_list)[catdex]
      
      
      #make plots for c1 and c0 direction
      # some categories have no pathways significant for condition, just return null
      
      summplots_conds <- lapply( c(c1, c0) , function(cond){
        
        #get result tbale for this condition
        cpres_cat_cond <- cpres_cat[cpres_cat$condition == cond,,drop=F]
        
        #if no conditions, it will haev nrow=0 so just return null
        if(nrow(cpres_cat_cond) == 0){ return() }
        
        # subselect categories if more than 30 total, use just top 5 per pathway
        if(nrow(cpres_cat_cond) > 30){
          #select the ones to pick
          cpres_cat_cond_sub <- cpres_cat_cond %>%
            group_by(cluster) %>%
            top_n(n=5, wt = -log10(padj)) %>%
            top_n(n=5, wt = abs(NES)) %>%
            as.data.frame()
          
          #get them, doing it this way allows viewing shared pathways
          cpres_cat_cond <- cpres_cat_cond[cpres_cat_cond$pathway %in% cpres_cat_cond_sub$pathway,]
          
        }
        
        
        #make sure orders are proper
        # for clusters:
        cpres_cat_cond$cluster <- factor(cpres_cat_cond$cluster, levels = unique(cpres_cat_cond$cluster))
        
        
        #for pathways
        cpres_cat_cond$pathway <- factor(cpres_cat_cond$pathway, levels = rev(unique(cpres_cat_cond$pathway)))
        
        
        ggplot(cpres_cat_cond, aes(x=cluster, y=pathway ,size = -log10(padj), col = NES))+
          geom_point()+
          theme_linedraw()+
          theme(axis.text=element_text(size=cp.font.size),
                axis.text.x = element_text(angle = 45, vjust = 1, hjust=1) )+
          scale_color_gradient2(low = 'steelblue', high = 'red', mid = 'white', midpoint = 0, name = 'Normalized\nEnrichment\nScore')+
          scale_size(range=c(2,6), name = '-log10(padj)')+
          xlab('Cluster')+ylab('')+
          ggtitle(catname, subtitle = cond)
        
        
      }) # close cross-condition loop for summary plots
      
      names(summplots_conds) <- c(c1, c0)
      
      summplots_conds
      
      
    })
    
    
    names(summplots_cats) <- names(cat_cpres_list)
    
    
    
    #print them to pdfs...
    
    pwayoutdir <- paste0(outdir_int, '/pathwayanalysis_crosscondition/', lab, '/')
    
    
    summarypdf <- paste0(pwayoutdir, '/SummaryDotPlots.pdf')
    pdf(summarypdf, width = 7, height = 7)
    
    print(summplots_cats)
    
    while (!is.null(dev.list()))  dev.off()
    
    return(summplots_cats)
    
    
    
  }) # close summary plot across conditons loop
  
  
  names(pathwaysummplots_comps) <- comps$labels
  
  
  
  
  ### for easily reproducing plots and etc, save them as R objects...
  pwayoutdir <- paste0(outdir_int, '/pathwayanalysis_crosscondition/')
  pathway_results <- .flatten_pathway_gsea_results(comps, pathway_analysis_mainlist_comps)

  DE_pathways_plot_objects_list <- list(
    comps = comps,
    de_results = de_results,
    pathway_results = pathway_results,
    pathwaysummplots_comps = pathwaysummplots_comps,
    pathway_cluster_plots = pathway_analysis_mainlist_comps
  )
  
  
  #save object sizes...
  # pwayobjsizedf <- data.frame(obj = names(DE_pathways_plot_objects_list))
  # pwayobjsizedf$size_bytes <- sapply(DE_pathways_plot_objects_list, object.size, simplify = T)
  #
  #
  # pwayoutdir <- paste0(outdir_int, '/pathwayanalysis_crosscondition/')
  #
  # objsizefile <- paste0(pwayoutdir, '/OBJSIZES_DE_pathways_plot_objects_list.csv')
  #
  # write.csv(pwayobjsizedf, objsizefile, quote = F, row.names = F)
  
  
  
  DE_pathways_plot_objects_list_file <- paste0(pwayoutdir, '/DE_pathways_plot_objects_list.rds')
  
  saveRDS(DE_pathways_plot_objects_list, DE_pathways_plot_objects_list_file)
  
  
  
  pways_output_list <- list(
    pathway_results = pathway_results,
    pathwaysummplots_comps = pathwaysummplots_comps,
    pathway_cluster_plots = pathway_analysis_mainlist_comps
  )

  return(pways_output_list)
  
  
}












#' Overrepresentation analysis for cross condition pathway analysis
#'
#' This function is a modular component of the scRNAseq pipeline. Perform OverRepresentation Analysis (ORA) via the ClusterProfiler package on the results of differential expression (DE) analysis for cross-condition comparison. Multiple conditions are supported. ClusterProfiler objects and tables are saved.
#'
#' @param de_results harmonized data.frame from `scDAPP::de_across_conditions_module()`.
#' @param pathways data.frame, the output of `scDAPP::preppathways_pathwayanalysis_crosscondition_module()`
#' @param sample_metadata data.frame with sample names and conditions, same as in `scDAPP::de_across_conditions_module()`, see that function's documentation for description.
#' @param comps data.frame with conditions to test in GSEA, same as in `scDAPP::de_across_conditions_module()`, see that function for description
#' @param pathway_padj_thres numeric, threshold for significance of pathway enrichment after multiple test correction, passed to qvalueCutoff in `clusterProfiler::enricher`
#' @param pwaycats UNTESTED CURRENTLY. character vector of msigdb pathways data.frame in gs_subcat to run.
#' @param workernum integer. number of CPUs. default = 1.
#' @param outdir_int string, directory to save pathways to. Will create a sub-directory called "overrepresentation_pathway_analysis" and save inside of there.
#'
#' @return Flat `ora_results` data.frame (pathway x cluster x category x comparison x direction).
#'   Per-cluster CSVs are still written under `outdir_int`.
#' @export
#'
#' @examples
#' \dontrun{
#'
#'
#' # FIRST: run `scDAPP::de_across_conditions_module()`. the output object of that is used as the main input for this pathway analysis function.
#'
#' # SECOND: prep pathways before running
#' pathways <- scDAPP::preppathways_pathwayanalysis_crosscondition_module(species = species,
#' outdir_int = outdir_int)
#'
#' # THIRD: run the pathway analysis.
#' # see `scDAPP::de_across_conditions_module()` for a description of the sample_metadata and comps files.
#' pways_output_list <- scDAPP::ORA_crosscondition_module(
#' de_results = de_results,
#' pathways = pathways,
#' sample_metadata = sample_metadata,
#' comps = comps,
#' workernum = 6,
#' outdir_int = outdir_int
#' )
#'
#' # Access the output files
#' # note that all results will be saved to outdir_int.
#' ora_results <- ORA_crosscondition_module(...)
#' subset(ora_results, label == "KO1_vs_Control")
#'
#' }
ORA_crosscondition_module <- function(de_results,
                                      pathways,
                                      sample_metadata,
                                      comps,
                                      crossconditionDE_padj_thres,
                                      crossconditionDE_lfc_thres,
                                      crossconditionDE_min.pct,
                                      pathway_padj_thres,
                                      pwaycats,
                                      workernum,
                                      outdir_int
){
  
  
  require(clusterProfiler)
  require(foreach)
  require(doParallel)
  require(parallel)
  
  if (!is.data.frame(de_results)) {
    stop("de_results must be a data.frame from de_across_conditions_module().", call. = FALSE)
  }
  comps <- .normalize_comps(comps)
  
  
  
  
  
  #### picking default categories
  # for now hardcode these
  if(missing(pwaycats)){
    # pwaycats <- c("HALLMARK", "GO_BP", "GO_MF", "GO_CC", "CP_REACTOME", "CP_KEGG", "TFT_GTRD", "TFT_TFT_Legacy")
    pwaycats <- c("HALLMARK", "GO_BP", "GO_MF", "GO_CC", "CP_REACTOME", "CP_KEGG", "TFT_GTRD", "TFT_TFT_Legacy")
  }
  
  if(missing(pathway_padj_thres)){
    pathway_padj_thres <- 0.1
  }
  if(missing(workernum)){
    workernum <- 1
  }
  
  
  ### pathway analysis for each condition comparison
  
  # for each condition comparison,
  # for each cluster
  # do pathway analysis and make plot
  
  
  
  
  
  compslen <- seq_len(nrow(comps))
  
  pathway_analysis_mainlist_comps <- lapply(compslen, function(compidx){
    
    c1 <- comps$c1[compidx]
    c0 <- comps$c0[compidx]
    lab <- comps$label[compidx]
    
    message(lab)
    
    de_sub <- de_results[de_results$label == lab, , drop = FALSE]
    clusters <- unique(de_sub$cluster)
    
    pwayoutdir <- paste0(outdir_int, '/overrepresentation_pathway_analysis/', lab, '/')
    if( !dir.exists(pwayoutdir) ){ dir.create(pwayoutdir, recursive = T) }
    
    
    ### loop thru pathway categories
    names(pwaycats) <- pwaycats
    
    #set gene universe
    pwaycat <- pwaycats[1] #for testing
    
    pathway_analysis_mainlist <- lapply(pwaycats, function(pwaycat){
      
      message('\n\n', pwaycat, '\n\n')
      
      
      #get pways and genes in this category
      term2gene <- pathways[pathways$gs_subcat == pwaycat,c('gs_name', 'gene_symbol')]
      
      
      #get list of pathways upreg in each cluster
      
      cl <- parallel::makeCluster(workernum, rscript_args = c("--no-init-file", "--no-site-file", "--no-environ"))
      doParallel::registerDoParallel(cl)
      
      
      clust = clusters[1] #for test
      
      #pwayres_DE_across_conditions_per_cluster <- lapply(clusters, function(clust){
      pwayres_DE_across_conditions_per_cluster <- foreach(clust = clusters,
                                                          .packages = c('clusterProfiler'),
                                                          .export = c('de_sub', 'pathway_padj_thres', 'crossconditionDE_padj_thres', 'crossconditionDE_lfc_thres', 'crossconditionDE_min.pct', 'c1', 'c0', '.de_results_cluster_table', 'filter_significant_de_results'),
                                                          .noexport = c('pathways'),
                                                          .verbose = T) %dopar%
        {
          
          
          
          invisible(gc(full = T, reset = F, verbose = F))
          
          
          res <- .de_results_cluster_table(de_sub, de_sub$label[1], clust)
          if (is.null(res) || !nrow(res)) return(NULL)
          
          
          
          ## split by lfc sign and run ##
          sign_nums <- c(1, -1)
          ts <- 1 #test
          signres_l <- lapply(sign_nums, function(ts){
            
            # Same DEG rules as de_across_conditions_module / count_crosscondition_degs
            sig <- filter_significant_de_results(
              res,
              crossconditionDE_padj_thres,
              crossconditionDE_lfc_thres,
              crossconditionDE_min.pct
            )
            subres <- sig[sign(sig$logFC) == ts, , drop = FALSE]
            
            
            
            #important, select min num genes; let's say 7 genes min for good luck
            if(nrow(subres) < 7){
              return()
            }
            
            
            ## if not, we can proceed with clusterProfiler
            genenames <- subres$gene_symbol
            
            
            
            
            ## run cluster profiler; any error, let's just return a null, maybe dangerous
            tryCatch(
              
              expr = {
                
                ora_res <- enricher(genenames,
                                    TERM2GENE = term2gene,
                                    qvalueCutoff = pathway_padj_thres,
                                    pvalueCutoff = 1
                                    
                )
                
              },
              
              error = function(e){
                
                ora_res <- data.frame()
                
              }
              
            )
            
            
            
            #clusterProfiler has its own trycatch, which returns NULL... i guess we'll use ours though
            if(length(ora_res) == 0){ora_res <- data.frame()}
            
            
            #return null if none 
            if(nrow(ora_res) == 0){return()}
            
            #don't use their useless class, just a less useful data.frame
            ora_res <- as.data.frame(ora_res)
            
            
            # add a direction column
            ora_res$Direction <- ifelse(ts == 1, c1, c0)
            
            
            
            return(ora_res)
            
          })
          
          
          
          #remove empty res
          signres_l <- signres_l[lengths(signres_l) > 0]
          
          if(length(signres_l) == 0){return()} # if there are just no pathways, return null
          
          signres <- dplyr::bind_rows(signres_l)
          
          
          #return the result table and the plot
          return(signres)
          
          
          
          
        } #per-cluster loop for this pathway category loop end
      
      parallel::stopCluster(cl)
      
      
      names(pwayres_DE_across_conditions_per_cluster) <- clusters
      
      
      return(pwayres_DE_across_conditions_per_cluster)
      
    }) #per-pathway category loop end
    
    
    
    
    ### remove all NULLS (clusters with no pathways)
    
    # remove null categories, ie entire category had no significant pathways
    
    #recursively set all missing to 0
    pathway_analysis_mainlist = lapply(pathway_analysis_mainlist, function(pwayres_cats){
      
      pwayres_cats <- lapply(pwayres_cats, function(pwayres_clusts){
        pwayres_clusts[lengths(pwayres_clusts) > 0]
      })
      
      pwayres_cats[lengths(pwayres_cats) > 0]
      
      
    })
    
    #remove any missing categories
    pathway_analysis_mainlist <- pathway_analysis_mainlist[lengths(pathway_analysis_mainlist)>0]
    
    
    invisible(gc(full = T, reset = F, verbose = F))
    
    
    #save pway analysis for this comparison
    
    # loop over this new pwayruns, since some categories theoretically don't have any enriched though unlikely
    pwayruns <- names(pathway_analysis_mainlist)
    
    
    #for each category, get clster res in that cateogry,
    # for each cluster, save the up/down csv and plots
    invisible(
      finalpwayouts <- lapply(pwayruns, function(pwaycat){
        
        
        
        # message(pwaycat)
        
        
        subcatout <- paste0(pwayoutdir, '/', pwaycat, '/')
        
        # dir.create(subcatout) --> do this with recursive later, maybe prevent even making it if all don't work
        
        clustres <- pathway_analysis_mainlist[[pwaycat]]
        
        clusters <- names(clustres)
        
        
        #for each cluster, get csv and save
        numpways <- lapply(clusters, function(clust){
          
          
          
          
          pwayres_DE_across_conditions_per_cluster <- clustres[[clust]]
          
          
          
          if(is.null(pwayres_DE_across_conditions_per_cluster)){return()}
          
          #this is the table
          signres <- pwayres_DE_across_conditions_per_cluster
          
          #save csv
          
          # prep subfolder
          suppressWarnings(dir.create(subcatout, recursive = T))
          
          
          #csv file
          subcatout_clustdir_gseares <- paste0(subcatout, '/significantpathways_', clust, '.csv')
          
          
          
          
          write.csv(signres, subcatout_clustdir_gseares, quote = F, row.names = F)
          
          
          
          
          nrow(signres)
          
          
          
          
          
          
          
        } ) # close clusters lapply
        
        
      }) # close saving loop for all categories
      
    ) # close invisible wrap around lapply
    
    
    
    
    
    
    
    
    
    
    ### close any open devices
    while (!is.null(dev.list()))  dev.off()
    
    
    
    
    return(pathway_analysis_mainlist)
    
    
  }) # close cross condition lapply
  
  
  names(pathway_analysis_mainlist_comps) <- comps$labels
  
  
  
  
  #remove big objects
  rm(pathways)
  invisible(gc(full = T, reset = F, verbose = F))
  
  
  
  
  ora_results <- .flatten_ora_results(comps, pathway_analysis_mainlist_comps)

  pwayoutdir <- paste0(outdir_int, '/overrepresentation_pathway_analysis/')
  DE_pathways_plot_objects_list_file <- paste0(pwayoutdir, '/DE_ORA_list_object.rds')

  saveRDS(
    list(
      comps = comps,
      de_results = de_results,
      ora_results = ora_results
    ),
    DE_pathways_plot_objects_list_file
  )

  return(ora_results)
  
  
}

















#' Compositional analysis comparing proportional abundnace across conditions for integrated Seurat objects
#'
#' This is a modular component of the scDAPP scRNAseq pipeline. Perform compositional analysis across conditions to compare the proportion of cell types. Supports multiple condtions (A vs B vs C). Supports "pseudobulk" replicate-aware analysis as implemented in the propeller test with arcsin transformation, as recommended by Simmons 2022 (https://doi.org/10.1101/2022.02.04.479123). When \code{comps$formula} includes an intercept random effect such as \code{(1|Patient)}, paired testing uses limma \code{duplicateCorrelation} with blocked \code{lmFit}. Alternatively supports old-school, non-replicate aware chisq test as implemented by the `prop.test()` function.
#'
#' @param sobjint integrated Seurat object. Metadata should have two columns: "Condition" corresponding to the A vs B conditions to compare across, and "Code" corresponding to sample / replicate names. A third column for clusters or celltypes should also be in the metadata and the name of that column will be passed to the `grouping_variable` parameter.
#' @param comps data.frame with c0 (reference), c1 (test), and optional formula, contrast, label columns. See `normalize_comps()`. Formulas with \code{(1|var)} enable paired / blocked propeller.
#' @param sample_metadata data.frame with Sample, Condition, Code, and any covariates (including random-effect block columns) in `comps$formula`.
#' @param outdir_int string, path to save results to. Will create a sub-directory called "compositional_proportion_analysis" and save inside of there.
#' @param grouping_variable string, column name of identity in Seurat object meta.data to stratify DE by. For example, clusters or celltype. Will perform A vs B DE in each of these groupings. Default is "seurat_clusters"
#' @param compositional_test string, which test to use, either "propeller" or "chisq"; will use chisq if not set and issue a warning
#' @param fill_barplots T/F, whether to "fill" the annotation barplots on the side of the heatmap, normalizing to 1; default = T
#' @param DE_test string, same as pipeline DE_test; selects contrast grammar for propeller
#'   (\code{edger_expr} for EdgeR*/Dream, \code{deseq2} for DESeq2*). Default \code{"EdgeR-LRT"}.
#'
#' @return A list with `composition_results` (flat data.frame: one row per cluster x comparison),
#'   `composition_plots` (named list of per-comparison ComplexHeatmap objects), and
#'   `globalcomposition` (cell counts and proportion tables).
#' @export
#'
#' @examples
#' \dontrun{
#'
#'
#' # `sample_metadata` looks like this:
#' Sample,Condition,Code
#' SampleXYZ1,Control,Control1
#' SampleXYZ2,Control,Control2
#' SampleABC1,KO1,KO1_1
#' SampleABC2,KO1,KO1_2
#' SampleJKL1,KO2,KO2_1
#' SampleJKL2,KO2,KO2_1
#'
#'
#' # `comps` looks like this:
#' c1,c2
#' KO1,Control
#' KO2,Control
#' KO1,K2
#'
#'
#' ## Run the analysis ##
#' comp_result <- compositional_analysis_module(sobjint,
#' comps,
#' sample_metadata,
#' outdir_int,
#' grouping_variable,
#' compositional_test = 'propeller',
#' DE_test = 'EdgeR-LRT')
#'
#' # Flat table for one comparison:
#' subset(comp_result$composition_results, label == "KO1_vs_Control")
#' # Heatmap for that comparison:
#' comp_result$composition_plots$KO1_vs_Control
#'
#' #Get some global information including cell numbers, proportions
#' # table of cell numbers
#' comp_result$globalcomposition$cellstab
#'
#' #table of cell proportions, ie the table above with each column divided by column sum
#' comp_result$globalcomposition$proptab
#'
#' #heatmap of cell proportions, without any statistical testing
#' comp_result$globalcomposition$hmprop
#' }
#'
compositional_analysis_module <- function(sobjint,
                                          comps,
                                          sample_metadata,
                                          outdir_int,
                                          grouping_variable,
                                          compositional_test,
                                          fill_barplots,
                                          DE_test = "EdgeR-LRT"
){

  require(ComplexHeatmap)
  require(circlize)
  require(speckle) #package with propeller test

  if(missing(grouping_variable)){grouping_variable = "seurat_clusters"}
  if(missing(compositional_test)){warning("No compositional test selected, will use chisq test"); compositional_test = 'chisq'}
  if(missing(fill_barplots)){fill_barplots = T}
  if (missing(DE_test) || is.null(DE_test)) DE_test <- "EdgeR-LRT"

  comps <- .apply_contrast_defaults(.normalize_comps(comps), DE_test)
  if (identical(compositional_test, "propeller")) {
    .preflight_propeller_comps(sample_metadata, comps, DE_test)
  } else {
    .preflight_comps_design(
      sample_metadata, comps, DE_test,
      Pseudobulk_mode = FALSE,
      check_de_formula_rules = FALSE
    )
  }
  contrast_style <- .contrast_style(DE_test)

  outdir_comp <- paste0(outdir_int, '/compositional_proportion_analysis/')
  dir.create(outdir_comp, recursive = T)

  ### get the composition table

  #num cells per cluster table
  md <- sobjint@meta.data
  cellstab <- table(md[,grouping_variable], md$Code)



  #prop table, divide num cells by total cells for each samp
  # ensure no div by zero
  cellstab2 <- cellstab[,!Matrix::colSums(cellstab) == 0]
  sample_metadata <- sample_metadata[sample_metadata$Code %in% colnames(cellstab2),]
  proptab <- t(t(cellstab2)/Matrix::colSums(cellstab2))

  #make a heatmap of the prop table
  col_fun = circlize::colorRamp2(c(0, max(proptab)), c( "white", "red"))

  hmprop <- ComplexHeatmap::Heatmap(proptab, name = "Proportion", col = col_fun,
                                    rect_gp = gpar(col = "black", lwd = 0.1),
                                    border_gp = gpar(col = "black", lwd = 1),
                                    column_split = sample_metadata$Condition,
                                    cluster_rows = T, cluster_columns = F,
                                    cell_fun = function(j, i, x, y, width, height, fill) {
                                      grid.text(sprintf("%.2f", proptab[i, j]), x, y, gp = gpar(fontsize = 10))
                                    })



  ## save global tables/heatmap
  cellnumfile <- paste0(outdir_comp, '/NumberCells.csv')
  write.csv(cellstab, cellnumfile, quote = F)

  cellpropfile <- paste0(outdir_comp, '/ProportionCells.csv')
  write.csv(proptab, cellpropfile, quote = F)

  hmprop_file <- paste0(outdir_comp, '/HeatmapProportions.pdf')
  pdf(hmprop_file, height = 7, width = 7)
  print(hmprop)
  dev.off()


  # proplong <- reshape2::melt(proptab)
  # proplong$Var1 <- factor(proplong$Var1, levels = str_sort(unique(proplong$Var1), numeric = T))

  # ggplot(proplong, aes(Var1, value, fill=Var2))+
  #   geom_col(position = 'fill')





  #for each comparison, do the compositional analysis

  compslen <- 1:nrow(comps)
  compidx = 1 #for testing





  composition_comps <- lapply(compslen, function(compidx){

    c0 <- comps$c0[compidx]
    c1 <- comps$c1[compidx]
    lab <- comps$label[compidx]
    design_formula <- .parse_comp_formula(comps$formula[compidx])
    has_re <- .formula_has_random_effect(design_formula)
    re_vars <- .formula_random_effect_vars(design_formula)
    if (has_re && length(re_vars) > 1L) {
      stop(
        "Propeller paired mode supports one (1|var) block; got: ",
        paste(re_vars, collapse = ", "),
        ". Use a single random-effect term in comps$formula for comparison '", lab, "'.",
        call. = FALSE
      )
    }
    block_var <- if (has_re) re_vars[1] else NULL
    fixed_formula <- if (has_re) .fixed_effects_formula(design_formula) else design_formula

    message(lab)

    subpmd <- sample_metadata[sample_metadata$Condition %in% c(c0, c1),]

    code_comps_order <- subpmd[subpmd$Condition == c1, "Code"]
    code_comps_order <- c(code_comps_order, subpmd[subpmd$Condition == c0, "Code"])

    condition_vector_ordering <- factor(
      subpmd[match(code_comps_order, subpmd$Code), "Condition"],
      levels = c(c1, c0)
    )

    ### use propeller if multiple samples ###

    if(compositional_test == 'propeller'){
      use_simple <- !has_re && .propeller_use_simple_wrapper(fixed_formula, sample_metadata)

      if (use_simple) {
        bigmd <- sobjint@meta.data
        md <- bigmd[bigmd$Condition %in% c(c0, c1), ]
        subpmd <- sample_metadata[sample_metadata$Condition %in% c(c0, c1), ]

        md[, grouping_variable] <- factor(
          md[, grouping_variable],
          levels = stringr::str_sort(unique(md[, grouping_variable]), numeric = TRUE)
        )
        md$Code <- factor(md$Code, levels = unique(subpmd$Code))
        md$Condition <- factor(md$Condition, levels = c(c1, c0))

        pres <- speckle::propeller(
          clusters = md[, grouping_variable],
          sample = md$Code,
          group = md$Condition,
          transform = 'asin'
        )
        pres <- .propeller_format_pres(pres, c1, c0)
      } else {
        require(limma)
        coldata <- .prepare_pseudobulk_coldata(sample_metadata, design_formula)
        bigmd <- sobjint@meta.data
        md <- bigmd[bigmd$Code %in% coldata$Code, ]
        md[, grouping_variable] <- factor(
          md[, grouping_variable],
          levels = stringr::str_sort(unique(md[, grouping_variable]), numeric = TRUE)
        )
        md$Code <- factor(md$Code, levels = coldata$Code)
        prop.list <- speckle::getTransformedProps(
          clusters = md[, grouping_variable],
          sample = md$Code,
          transform = "asin"
        )
        pd <- .propeller_design(sample_metadata, design_formula)
        contr <- tryCatch(
          .propeller_contrast_matrix(
            pd$design, comps$contrast[compidx], c1, c0, style = contrast_style
          ),
          error = function(e) {
            warning(
              "Skipping propeller for comparison '", lab, "': ",
              conditionMessage(e),
              call. = FALSE
            )
            NULL
          }
        )
        if (is.null(contr)) return(NULL)
        if (has_re) {
          block <- pd$coldata[[block_var]][match(rownames(pd$design), pd$coldata$Code)]
          if (anyNA(block)) {
            stop(
              "Missing values in block variable '", block_var,
              "' for propeller paired comparison '", lab, "'.",
              call. = FALSE
            )
          }
          # Each Code must map to a single block level
          code_block <- unique(pd$coldata[, c("Code", block_var), drop = FALSE])
          if (any(duplicated(code_block$Code))) {
            stop(
              "Block variable '", block_var,
              "' is not unique per Code for comparison '", lab, "'.",
              call. = FALSE
            )
          }
          pres <- .propeller_ttest_blocked(
            prop.list,
            design = pd$design,
            contrasts = contr,
            block = block,
            robust = TRUE,
            trend = FALSE,
            sort = TRUE
          )
        } else {
          pres <- .propeller_ttest(
            prop.list,
            design = pd$design,
            contrasts = contr,
            robust = TRUE,
            trend = FALSE,
            sort = TRUE
          )
        }
        pres <- .propeller_format_pres(pres, c1, c0)
      }

      #subset table of proportions
      comp_proptab <- proptab[,subpmd$Code]


      #make sure columns are in right order
      comp_proptab <- comp_proptab[,code_comps_order]

      #match table of proporitons order with result
      comp_proptab <- comp_proptab[match(pres$BaselineProp.clusters, rownames(comp_proptab)),]

      #order by pval * sign of diff...
      # pres <- pres[ order( -log10(pres$P.Value) * log(pres$PropRatio) , decreasing = T) ,]
      # comp_proptab <- comp_proptab[match(pres$BaselineProp.clusters, rownames(comp_proptab)),]


      #add signficance marks
      pres$BaselineProp.clusters <- as.character(pres$BaselineProp.clusters)
      pres[pres$P.Value<0.05, "BaselineProp.clusters"] <- paste0(
        '* ',
        pres[pres$P.Value<0.05, "BaselineProp.clusters"],
        ' *'
      )

      rownames(comp_proptab) <- pres$BaselineProp.clusters


      # row annotation: pvalue (eprecated)
      #deal with underflow
      # neglog <- -log10(pres$P.Value)
      # if(any(neglog == Inf)){neglog[neglog==Inf] <- 210}
      # row_ha = rowAnnotation( "-log10P" = anno_barplot( neglog ) )

      #negative p value * propratio
      # directional_pval <- -log10(pres$P.Value) * sign(log(pres$PropRatio))
      # row_ha = rowAnnotation( "-log10P * PropRatio" = anno_barplot( directional_pval ) )

      #row annotation: PropRatio
      # row_ha = rowAnnotation( "LogPropRatio" = anno_lines( log(pres$PropRatio), smooth = T,
      #                                                      axis_param = list(direction = "reverse"))
      #                          )

      #row annotation: barplots of observed mean proportions (c1, c0)
      propmat_annot <- .composition_sample_prop_means(comp_proptab, subpmd, c1, c0)
      if(fill_barplots==T){
        propmat_annot[sign(propmat_annot)==-1] <- propmat_annot[sign(propmat_annot)==-1] * -1
        propmat_annot_sum <- rowSums(propmat_annot)
        propmat_annot[,1] <- propmat_annot[,1] / propmat_annot_sum
        propmat_annot[,2] <- propmat_annot[,2] / propmat_annot_sum
      }
      propmat_annot <- as.matrix(propmat_annot[,1:2])
      row_ha = rowAnnotation( "Proportions" = anno_barplot( propmat_annot,
                                                            # axis_param = list(direction = "reverse"),
                                                            gp = gpar(fill = c('firebrick', 'steelblue'), col = c('firebrick', 'steelblue') ))
      )


      hmprop_comp <- Heatmap(comp_proptab, name = "Proportion", col = col_fun,
                             rect_gp = gpar(col = "black", lwd = 0.1),
                             border_gp = gpar(col = "black", lwd = 1),
                             column_split = condition_vector_ordering,
                             cluster_rows = F, cluster_columns = F,
                             right_annotation = row_ha,
                             # width = ncol(comp_proptab)*unit(20, "mm"),
                             cell_fun = function(j, i, x, y, width, height, fill) {
                               grid.text(sprintf("%.2f", comp_proptab[i, j]), x, y, gp = gpar(fontsize = 10))
                             })


      compres <- pres

    } else{



      # if no replicates, use prop.test (formula ignored)

      md <- sobjint@meta.data
      md <- md[md$Condition %in% c(c0, c1), ]
      cells_condtab <- table(md[, grouping_variable], md$Condition)

      tots <- table(md$Condition)
      c1tot <- tots[c1]
      c0tot <- tots[c0]

      cells_condtab <- cells_condtab[, c(c1, c0), drop = FALSE]

      #remove all zero clusters
      cells_condtab = cells_condtab[Matrix::rowSums(cells_condtab)>0,]


      #in each cluster, calculate proportion of c1 and c2 from overall c1 and c2, then compare with prop.test
      clust_proptest_res <- lapply(1:nrow(cells_condtab), function(i){

        clust <- rownames(cells_condtab)[i]

        #get this clusters' props
        vec <- cells_condtab[i,,drop=F]

        x = c(vec[1, 1], vec[1, 2])
        n = c(c1tot, c0tot)

        pt <- prop.test(x = x, n = n)

        ptdf <- data.frame(cluster = clust,
                           c1prop = pt$estimate[1],
                           c0prop = pt$estimate[2],
                           asin_ratio = plogis(pt$estimate[1]) / plogis(pt$estimate[2]),
                           difference = pt$estimate[1] - pt$estimate[2],
                           p = pt$p.value,
                           row.names = NULL)


      })




      #make results to data.frame
      clust_proptest_resdf <- dplyr::bind_rows(clust_proptest_res)

      # asin ratio = Inf, means c0 was zero..
      clust_proptest_resdf$asin_ratio[clust_proptest_resdf$asin_ratio == Inf] <- 1

      #add FDR
      clust_proptest_resdf$FDR <- p.adjust(clust_proptest_resdf$p)

      #sort by directional pvalue?
      # dirpval <- -log1p(clust_proptest_resdf$p) * sign(clust_proptest_resdf$difference)
      # dirpval <- order(dirpval, decreasing = T)
      # clust_proptest_resdf <- clust_proptest_resdf[order(clust_proptest_resdf$difference, decreasing = T),]

      #sort by difference, deprecated
      clust_proptest_resdf <- clust_proptest_resdf[order(clust_proptest_resdf$difference, decreasing = T),]

      #sort by "asin ratio"
      # clust_proptest_resdf <- clust_proptest_resdf[order(clust_proptest_resdf$asin_ratio, decreasing = T),]

      #subset table of proportions
      comp_proptab <- proptab[,subpmd$Code]

      #make sure columns are in right order
      comp_proptab <- comp_proptab[,code_comps_order]

      #match table of proporitons order with result
      comp_proptab <- comp_proptab[match(clust_proptest_resdf$cluster, rownames(comp_proptab)),]

      #put asterisk if significant
      clust_proptest_resdf[clust_proptest_resdf$p < 0.05, "cluster"] <- paste0( '* ', clust_proptest_resdf[clust_proptest_resdf$p < 0.05, "cluster"], ' *')

      rownames(comp_proptab) <- clust_proptest_resdf$cluster


      #annotate with pvalue, deprecated
      # #deal with underflow
      # neglog <- -log10(clust_proptest_resdf$p)
      # if(any(neglog == Inf)){neglog[neglog==Inf] <- 210}
      # row_ha = rowAnnotation( "-log10P" = anno_barplot( neglog ) )

      #annotate with pval * sign of difference
      # neglog <- -log10(clust_proptest_resdf$p)
      # if(any(neglog == Inf)){neglog[neglog==Inf] <- 210}
      # directional_pval <- neglog * sign(clust_proptest_resdf$difference)
      # row_ha = rowAnnotation( "-log10P * PropRatio" = anno_barplot( directional_pval ) )

      #row annotation: barplots of proportions
      propmat_annot <- clust_proptest_resdf[,2:3]
      if(fill_barplots==T){
        propmat_annot$sum <- rowSums(propmat_annot)
        propmat_annot[,1] <- propmat_annot[,1] / propmat_annot[,3]
        propmat_annot[,2] <- propmat_annot[,2] / propmat_annot[,3]
      }
      propmat_annot <- as.matrix(propmat_annot[,1:2])
      row_ha = rowAnnotation( "Proportions" = anno_barplot( propmat_annot,
                                                            # axis_param = list(direction = "reverse"),
                                                            gp = gpar(fill = c('firebrick', 'steelblue'), col = c('firebrick', 'steelblue') ))
      )


      hmprop_comp <- Heatmap(comp_proptab, name = "Proportion", col = col_fun,
                             rect_gp = gpar(col = "black", lwd = 0.1),
                             border_gp = gpar(col = "black", lwd = 1),
                             column_split = condition_vector_ordering,
                             cluster_rows = F, cluster_columns = F,
                             right_annotation = row_ha,
                             # width = ncol(comp_proptab)*unit(20, "mm"),
                             cell_fun = function(j, i, x, y, width, height, fill) {
                               grid.text(sprintf("%.2f", comp_proptab[i, j]), x, y, gp = gpar(fontsize = 10))
                             })



      compres <- clust_proptest_resdf

    }




    ### save it ###
    if(compositional_test == 'propeller'){
      outdir_comp_thiscomparison <- paste0(outdir_comp,
                                           '/Propellertest/',
                                           lab, '/'
      )
    } else{
      outdir_comp_thiscomparison <- paste0(outdir_comp,
                                           '/2PropZtest/',
                                           lab, '/'
      )
    }

    dir.create(outdir_comp_thiscomparison, recursive = T)

    compresfile <- paste0(outdir_comp_thiscomparison, '/CompositionAnalysis.csv')
    write.csv(compres, file = compresfile, quote = F)

    hmprop_comp_file <- paste0(outdir_comp_thiscomparison, '/HeatmapProportions.pdf')
    pdf(hmprop_comp_file, height = 7, width = 7)
    print(hmprop_comp)
    dev.off()

    return(list(compres = compres,
                hmprop_comp = hmprop_comp))

  })

  names(composition_comps) <- comps$labels

  composition_results <- .flatten_composition_results(comps, composition_comps)
  composition_plots <- lapply(composition_comps, function(x) x$hmprop_comp)
  names(composition_plots) <- names(composition_comps)

  globalcomposition <- list(cellstab = cellstab,
                            proptab = proptab,
                            hmprop = hmprop)

  comp_out <- list(
    composition_results = composition_results,
    composition_plots = composition_plots,
    globalcomposition = globalcomposition
  )

  return(comp_out)

}








# # ### testing ###
# #
# library(Seurat)
# library(tidyverse)
# library(edgeR)
#
# # rm(sobj)
# sobjint <- readRDS('~/Dropbox/Result_from_Alex/deyoudata/hto_prostate/outs/scRNAseqpipeline/multisample_integration/data_objects/Seurat-object_integrated.rds')
# sample_metadata <- read.csv('~/Dropbox/Result_from_Alex/deyoudata/hto_prostate/data/metadata/sample_metadata.csv')
# comps <- read.csv('~/Dropbox/Result_from_Alex/deyoudata/hto_prostate/data/metadata/comps.csv')
# grouping_variable <- 'seurat_clusters'
# Pseudobulk_mode <- T
# outdir_int <- 'tests/de_module_test' ; dir.create(outdir_int, recursive = T)
# cluster_prefix = T
# # min.pct <- 0
# assay = 'RISC'
# slot = 'data'
# crossconditionDE_padj_thres = 0.1
# crossconditionDE_lfc_thres = 0
# crossconditionDE_min.pct = 0.1
#
#
# set.seed(2022)
#
#
# m_bycluster_crosscondition_de_comps <- de_across_conditions_module(
#   sobjint = sobjint,
#   sample_metadata = sample_metadata,
#   comps = comps,
#   outdir_int = outdir_int,
#   grouping_variable = 'seurat_clusters',
#   Pseudobulk_mode = T
#
# )
#
#
# m_bycluster_crosscondition_de_comps_WILCOX <- de_across_conditions_module(
#   sobjint = sobjint,
#   sample_metadata = sample_metadata,
#   comps = comps,
#   grouping_variable = 'seurat_clusters',
#   Pseudobulk_mode = T
#
# )
#
#
#
# library(Seurat)
# library(tidyverse)
# library(fgsea)
# library(foreach)
# library(doParallel)
# library(parallel)
#
# rm(sobj)
# rm(sobjint)
#
# sample_metadata <- read.csv('~/Dropbox/Result_from_Alex/deyoudata/hto_prostate/data/metadata/sample_metadata.csv')
# comps <- read.csv('~/Dropbox/Result_from_Alex/deyoudata/hto_prostate/data/metadata/comps.csv')
#
#
# outdir_int <- 'tests/de_module_test' ; dir.create(outdir_int, recursive = T)
#
# deg.weight = 'auto'
# pathway_padj_thres <- 0.1
# workernum <- 4
# cp.font.size <- 5
# pwaycats <- c("HALLMARK", "GO_BP", "GO_MF", "GO_CC", "CP_REACTOME", "CP_KEGG", "TFT_GTRD", "TFT_TFT_Legacy")
#
#
#
#
#
#
# pathways <- preppathways_pathwayanalysis_crosscondition_module(species = species,
#                                                                outdir_int = outdir_int)
#
#
#
#
# pways_output_list <- pathwayanalysis_crosscondition_module(
#   m_bycluster_crosscondition_de_comps = m_bycluster_crosscondition_de_comps,
#   pathways = pathways,
#   sample_metadata = sample_metadata,
#   comps = comps,
#   workernum = 6,
#   outdir_int = outdir_int
# )
#
# beepr::beep()




# pways_output_list <- pathwayanalysis_crosscondition_module(
#   m_bycluster_crosscondition_de_comps = m_bycluster_crosscondition_de_comps,
#   pathways = pathways,
#   sample_metadata = sample_metadata,
#   deg.weight = "pval",
#   comps = comps,
#   workernum = workernum,
#   outdir_int = outdir_int
# )
