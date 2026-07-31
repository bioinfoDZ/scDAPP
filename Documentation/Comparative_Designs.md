# Comparative designs: `sample_metadata`, `comps`, and `scRNAseq_pipeline_runner()`

This cookbook shows how to set up common cross-condition comparisons in scDAPP. For general pipeline inputs and parameter lists, see [Usage.md](Usage.md).

## Quick rules

- **`c0`** = reference condition; **`c1`** = test condition. Positive log2FC (and higher proportions) favor **c1**.
- Legacy `c1,c2` files still work (`c2` is renamed to `c0`).
- **`DE_test` is pipeline-wide** — one value per `scRNAseq_pipeline_runner()` call, not per comps row.
- In **single-cell** mode, DE is pairwise `FindMarkers` on `c0` vs `c1`; `formula` is ignored for DE.
- **`formula` / covariates / random effects** apply to **pseudobulk DE** and **propeller** compositional analysis.
- **Compositional analysis:** `Pseudobulk_mode = TRUE` → propeller (respects `formula`, including blocked pairing); `FALSE` → two-proportion `prop.test()` per cluster.

### Required vs optional `comps` columns

**Required for every analysis type** (Wilcox, EdgeR/DESeq2, Dream, multi-row, confounders, interactions): only **`c0`** and **`c1`**.

**`contrast` and `label` are not required for any comparative analysis type.** Neither is `formula` unless you need a non-default design. Defaults if omitted or blank:

| Column | Default if omitted/blank |
|--------|--------------------------|
| `formula` | `~ Condition` |
| `contrast` | Condition **c1 vs c0** (same as `Condition;c1;c0`) |
| `label` | `{c1}_vs_{c0}` |

### What `contrast` does

Given your `formula` and columns in `sample_metadata`, **`contrast` selects which factor’s levels to report results for** in pseudobulk DE (and in propeller when the design is not a simple `~ Condition`).

- Format: semicolon triple `Variable;numerator;denominator` (DESeq2-style), e.g. `Condition;KO;Control` or `Batch;BatchB;BatchA`.
- The default Condition contrast is enough for simple A vs B and most Condition-focused runs.
- Non-default contrasts matter most for **multivariable** and **interaction** models: the same fitted formula can yield different reported comparisons depending on `contrast` (e.g. Condition effect adjusted for Batch vs Batch effect adjusted for Condition). See [§4](#4-multivariable-tests-with-confounders) and [§5](#5-interactions).
- Even when `contrast` targets another variable (e.g. Batch), keep **`c0`/`c1` as Condition levels** — they are still used for replicate checks, descriptive expression percentages, and default labeling.

**Examples below always include `contrast` and `label` (and often `formula`) to demonstrate the full comps table shape.** You can omit those columns in real files and get the defaults above.

---

## 1. Single-cell 1v1 (Wilcox / FindMarkers)

**When to use:** Few or no biological replicates per condition; exploratory cell-level DE; or when you intentionally want Seurat’s single-cell tests instead of pseudobulk.

**Compositional analysis:** two-proportion Z-test via `prop.test()` (not propeller).

### `sample_metadata.csv`

| Sample | Condition | Code |
|--------|-----------|------|
| Control_1 | Control | Control_1 |
| KO_1 | KO | KO_1 |

(`Code` is optional.)

### `comps.csv`

| c0 | c1 | contrast | label |
|----|----|----------|-------|
| Control | KO | Condition;KO;Control | KO_vs_Control |

### Runner

```r
scDAPP::scRNAseq_pipeline_runner(
  datadir = "path/to/datadir",
  outdir = "path/to/outdir",
  sample_metadata = "path/to/sample_metadata.csv",
  comps = "path/to/comps.csv",
  Pseudobulk_mode = FALSE,
  DE_test = "wilcox"   # default when Pseudobulk_mode = FALSE
)
```

### Other single-cell `DE_test` options

With `Pseudobulk_mode = FALSE`, `DE_test` may be any of:

`wilcox`, `wilcox_limma`, `bimod`, `t`, `negbinom`, `poisson`, `LR`, `MAST`, `DESeq2`

These map to Seurat `FindMarkers(..., test.use = ...)`. See `?Seurat::FindMarkers`. The Seurat `roc` test is **not** supported. Some options (e.g. `MAST`, `DESeq2`) require extra packages.

Example using MAST:

```r
scDAPP::scRNAseq_pipeline_runner(
  datadir = "path/to/datadir",
  outdir = "path/to/outdir",
  sample_metadata = "path/to/sample_metadata.csv",
  comps = "path/to/comps.csv",
  Pseudobulk_mode = FALSE,
  DE_test = "MAST"
)
```

---

## 2. Simple A vs B pseudobulk

**When to use:** At least ~2 biological replicates per condition; standard replicate-aware DE and cell-proportion testing.

**Compositional analysis:** propeller (default simple `~ Condition` design).

### `sample_metadata.csv`

| Sample | Condition | Code |
|--------|-----------|------|
| Control_1 | Control | Control_1 |
| Control_2 | Control | Control_2 |
| KO_1 | KO | KO_1 |
| KO_2 | KO | KO_2 |

### `comps.csv`

| c0 | c1 | contrast | label |
|----|----|----------|-------|
| Control | KO | Condition;KO;Control | KO_vs_Control |

### Runner — EdgeR LRT (default)

```r
scDAPP::scRNAseq_pipeline_runner(
  datadir = "path/to/datadir",
  outdir = "path/to/outdir",
  sample_metadata = "path/to/sample_metadata.csv",
  comps = "path/to/comps.csv",
  Pseudobulk_mode = TRUE,
  DE_test = "EdgeR-LRT"   # default when Pseudobulk_mode = TRUE
)
```

### Same files, other pseudobulk engines

Only `DE_test` changes; metadata and comps stay the same.

| `DE_test` | Notes |
|-----------|--------|
| `"EdgeR-LRT"` | Default; GLM + likelihood ratio test |
| `"EdgeR"` | Exact test for simple 2-group designs; LRT for more complex designs |
| `"DESeq2"` | Wald test; requires Bioconductor `DESeq2` |
| `"DESeq2-LRT"` | Likelihood ratio test; requires `DESeq2` |

```r
# DESeq2 Wald
scDAPP::scRNAseq_pipeline_runner(
  datadir = "path/to/datadir",
  outdir = "path/to/outdir",
  sample_metadata = "path/to/sample_metadata.csv",
  comps = "path/to/comps.csv",
  Pseudobulk_mode = TRUE,
  DE_test = "DESeq2"
)

# DESeq2 LRT
scDAPP::scRNAseq_pipeline_runner(
  datadir = "path/to/datadir",
  outdir = "path/to/outdir",
  sample_metadata = "path/to/sample_metadata.csv",
  comps = "path/to/comps.csv",
  Pseudobulk_mode = TRUE,
  DE_test = "DESeq2-LRT"
)
```

Do **not** put `(1|var)` random effects in `formula` with EdgeR/DESeq2 — use [Dream](#6-random-effects--paired-dream) instead.

---

## 3. A vs B vs C (multiple pairwise comparisons)

**When to use:** Three or more conditions; you want several pairwise contrasts in one pipeline run.

**Compositional analysis:** propeller (or `prop.test` if single-cell) once per comps row.

### `sample_metadata.csv`

| Sample | Condition | Code |
|--------|-----------|------|
| Control_1 | Control | Control_1 |
| Control_2 | Control | Control_2 |
| KO1_1 | KO1 | KO1_1 |
| KO1_2 | KO1 | KO1_2 |
| KO2_1 | KO2 | KO2_1 |
| KO2_2 | KO2 | KO2_2 |

### `comps.csv`

| c0 | c1 | contrast | label |
|----|----|----------|-------|
| Control | KO1 | Condition;KO1;Control | KO1_vs_Control |
| Control | KO2 | Condition;KO2;Control | KO2_vs_Control |
| KO1 | KO2 | Condition;KO2;KO1 | KO2_vs_KO1 |

All rows share the same pipeline-wide `DE_test`.

### Runner

```r
scDAPP::scRNAseq_pipeline_runner(
  datadir = "path/to/datadir",
  outdir = "path/to/outdir",
  sample_metadata = "path/to/sample_metadata.csv",
  comps = "path/to/comps.csv",
  Pseudobulk_mode = TRUE,
  DE_test = "EdgeR-LRT"
)
```

For single-cell mode across the same three conditions, keep the same comps file and set `Pseudobulk_mode = FALSE` (and optionally `DE_test = "wilcox"`).

---

## 4. Multivariable tests with confounders

**When to use:** Adjust for known technical or biological covariates (batch, sex, site, etc.) while testing Condition — or, with a non-default `contrast`, report the covariate effect adjusted for Condition.

**Compositional analysis:** propeller uses the same fixed-effects `formula` (via `propeller.ttest` when covariates are present) and respects `contrast`.

### `sample_metadata.csv`

Every variable named in `formula` must be a column here:

| Sample | Condition | Code | Batch |
|--------|-----------|------|-------|
| Control_1 | Control | Control_1 | BatchA |
| Control_2 | Control | Control_2 | BatchB |
| KO_1 | KO | KO_1 | BatchA |
| KO_2 | KO | KO_2 | BatchB |

### `comps.csv`

Same formula, two contrasts — usual Condition effect vs non-default Batch effect:

| c0 | c1 | formula | contrast | label |
|----|----|---------|----------|-------|
| Control | KO | `~ Condition + Batch` | Condition;KO;Control | KO_vs_Control_batch_adj |
| Control | KO | `~ Condition + Batch` | Batch;BatchB;BatchA | BatchB_vs_BatchA_cond_adj |

Keep **`c0`/`c1` as Condition levels** on both rows (replicate checks and descriptive pct still use Condition). Row 1 reports Condition **KO vs Control** adjusted for Batch. Row 2 reports Batch **BatchB vs BatchA** adjusted for Condition.

### Runner

```r
scDAPP::scRNAseq_pipeline_runner(
  datadir = "path/to/datadir",
  outdir = "path/to/outdir",
  sample_metadata = "path/to/sample_metadata.csv",
  comps = "path/to/comps.csv",
  Pseudobulk_mode = TRUE,
  DE_test = "EdgeR-LRT"   # or "DESeq2-LRT", "DESeq2", "EdgeR"
)
```

---

## 5. Interactions

**When to use:** The Condition effect may differ by another factor (e.g. Batch, Genotype × Treatment). Fit an interaction model with EdgeR/DESeq2 (fixed effects only).

**Compositional analysis:** propeller uses the fixed-effects design derived from `formula` and respects `contrast`.

### `sample_metadata.csv`

| Sample | Condition | Code | Batch |
|--------|-----------|------|-------|
| Mild_1 | COVID_MILD | Mild_1 | Ja001E |
| Mild_2 | COVID_MILD | Mild_2 | Ja005E |
| Sev_1 | COVID_SEV | Sev_1 | Ja001E |
| Sev_2 | COVID_SEV | Sev_2 | Ja005E |

### `comps.csv`

Single interaction-model row reporting the **main Condition effect** (not the interaction coefficient):

| c0 | c1 | formula | contrast | label |
|----|----|---------|----------|-------|
| COVID_MILD | COVID_SEV | `~ Condition * Batch` | Condition;COVID_SEV;COVID_MILD | SEV_vs_MILD_main_effect_in_interaction_model |

You can put simple, additive, and interaction formulas in one comps file (same run, same `DE_test`), and vary `contrast` under the interaction formula:

| c0 | c1 | formula | contrast | label |
|----|----|---------|----------|-------|
| COVID_MILD | COVID_SEV | `~ Condition` | Condition;COVID_SEV;COVID_MILD | SEV_vs_MILD_simple |
| COVID_MILD | COVID_SEV | `~ Condition + Batch` | Condition;COVID_SEV;COVID_MILD | SEV_vs_MILD_batch_adj |
| COVID_MILD | COVID_SEV | `~ Condition * Batch` | Condition;COVID_SEV;COVID_MILD | SEV_vs_MILD_main_effect_in_interaction_model |
| COVID_MILD | COVID_SEV | `~ Condition * Batch` | not implemented yet | SEV_x_Batch_interaction_term |

- Rows 1–2: Condition effect under simpler designs.
- Row 3: still the **main Condition effect** under `~ Condition * Batch` (reference-level Batch coding), **not** the Condition×Batch interaction coefficient.
- Row 4: intended target is the **Condition×Batch interaction coefficient** (LFC/p); `contrast` for that is **not implemented yet**.

### Runner

```r
scDAPP::scRNAseq_pipeline_runner(
  datadir = "path/to/datadir",
  outdir = "path/to/outdir",
  sample_metadata = "path/to/sample_metadata.csv",
  comps = "path/to/comps.csv",
  Pseudobulk_mode = TRUE,
  DE_test = "EdgeR-LRT"
)
```

**Contrast caveat:** The semicolon `contrast` form compares two levels of **one** factor in the design. It does **not** currently extract the Condition×Batch **interaction coefficient** (whether the Condition effect differs by Batch). That is listed under [Future work](#future-work--todos).

---

## 6. Random effects / paired (Dream)

**When to use:** Matched donors/patients across conditions (paired or blocked designs). Pseudobulk mixed models via Bioconductor **variancePartition** (`dream`).

**Compositional analysis:** when `formula` contains `(1|var)`, propeller automatically uses blocked limma (`duplicateCorrelation` + `lmFit`) with that blocking factor.

**Requirements:**

- `DE_test = "Dream"`
- At least one comps `formula` with an intercept random effect, e.g. `(1|Patient)`
- `sample_metadata` must include the block column (`Patient` below)
- Install: `BiocManager::install("variancePartition")`
- EdgeR/DESeq2 **reject** formulas with `(1|var)` — use Dream instead
- Only **intercept** random effects `(1|var)` are supported today; slopes like `(Condition|Patient)` are not

### `sample_metadata.csv`

Each patient ideally appears in both conditions:

| Sample | Condition | Code | Patient |
|--------|-----------|------|---------|
| Ctrl_P1 | Control | Ctrl_P1 | P1 |
| Treat_P1 | Treatment | Treat_P1 | P1 |
| Ctrl_P2 | Control | Ctrl_P2 | P2 |
| Treat_P2 | Treatment | Treat_P2 | P2 |
| Ctrl_P3 | Control | Ctrl_P3 | P3 |
| Treat_P3 | Treatment | Treat_P3 | P3 |

### `comps.csv`

| c0 | c1 | formula | contrast | label |
|----|----|---------|----------|-------|
| Control | Treatment | `~ Condition + (1\|Patient)` | Condition;Treatment;Control | Treatment_vs_Control_paired |

Paired + fixed confounder:

| c0 | c1 | formula | contrast | label |
|----|----|---------|----------|-------|
| Control | Treatment | `~ Condition + Batch + (1\|Patient)` | Condition;Treatment;Control | Treatment_vs_Control_paired_batch |

(Include `Batch` in `sample_metadata` if used.)

### Runner

```r
scDAPP::scRNAseq_pipeline_runner(
  datadir = "path/to/datadir",
  outdir = "path/to/outdir",
  sample_metadata = "path/to/sample_metadata.csv",
  comps = "path/to/comps.csv",
  Pseudobulk_mode = TRUE,
  DE_test = "Dream",
  workernum = 4   # BiocParallel workers for Dream
)
```

Interaction + pairing in one formula is allowed for Dream DE, e.g. `~ Condition * Batch + (1|Patient)`, with the same contrast caveats as in [§5](#5-interactions). Propeller still uses a single blocking factor from `(1|var)`.

---

## Summary

| Design | `Pseudobulk_mode` | Typical `DE_test` | Typical `formula` | Typical `contrast` |
|--------|-------------------|-------------------|-------------------|--------------------|
| Single-cell 1v1 | `FALSE` | `wilcox` (or other FindMarkers tests) | ignored for DE (default `~ Condition`) | `Condition;c1;c0` (informational; DE uses c0/c1) |
| Simple A vs B | `TRUE` | `EdgeR-LRT`, `EdgeR`, `DESeq2`, `DESeq2-LRT` | `~ Condition` | `Condition;c1;c0` |
| A vs B vs C | `TRUE` or `FALSE` | same as above / wilcox | `~ Condition`; multiple comps rows | `Condition;c1;c0` per row |
| Confounders (Condition effect) | `TRUE` | EdgeR / DESeq2 family | `~ Condition + Batch` | `Condition;c1;c0` |
| Confounders (Batch effect) | `TRUE` | EdgeR / DESeq2 family | `~ Condition + Batch` | `Batch;level1;level0` |
| Interactions (Condition main effect) | `TRUE` | EdgeR / DESeq2 family | `~ Condition * Batch` | `Condition;c1;c0` |
| Interactions (Interaction term) | `TRUE` | EdgeR / DESeq2 family | `~ Condition * Batch` | not implemented yet |
| Paired / random intercept | `TRUE` | `Dream` | `~ Condition + (1\|Patient)` | `Condition;c1;c0` |

---

## Future work / TODOs

- Document running comparative modules (`de_across_conditions_module`, pathway/GSEA, compositional analysis) outside the full pipeline runner.
- Support and document contrasts for interaction **coefficients** (Condition×covariate), not only factor-level triples (`Variable;numerator;denominator`).
