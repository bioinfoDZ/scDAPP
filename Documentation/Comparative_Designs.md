# Comparative designs: `sample_metadata`, `comps`, and `scRNAseq_pipeline_runner()`

This cookbook shows how to set up common cross-condition comparisons in scDAPP. For general pipeline inputs and parameter lists, see [Usage.md](Usage.md).

## Quick rules

- **`c0`** = reference condition; **`c1`** = test condition. Positive log2FC (and higher proportions) favor **c1**. Examples below use **`Control`** / **`Treatment`**.
- Legacy `c1,c2` files still work (`c2` is renamed to `c0`).
- **`DE_test` is pipeline-wide** — one value per `scRNAseq_pipeline_runner()` call, not per comps row.
- In **single-cell** mode, DE is pairwise `FindMarkers` on `c0` vs `c1`; **`formula` and `contrast` are ignored for DE**.
- **`formula` / covariates / random effects** apply to **pseudobulk DE** and **propeller** compositional analysis.
- **Compositional analysis:** `Pseudobulk_mode = TRUE` → propeller (same contrast grammar as the run’s `DE_test`); `FALSE` → two-proportion `prop.test()` per cluster.
- **Preflight** validates formulas/contrasts early (Rmd + DE/compositional modules). Sparse clusters that cannot estimate a contrast are **soft-skipped** (warning + empty result for that cluster), not a full abort.
- **Reference levels:** `Condition` (and other design factors) use level order = **first appearance in `sample_metadata`**. Under treatment coding, the first level is the design reference (intercept). Put Control (or the reference genotype) rows first when you want them as the reference. Pairwise comps still use `c0`/`c1` (blank `contrast` → Condition c1 vs c0).

### Required vs optional `comps` columns

**Required:** only **`c0`** and **`c1`**.

| Column | Default if omitted/blank |
|--------|--------------------------|
| `formula` | `~ Condition` |
| `contrast` | Style depends on `DE_test` (see below) |
| `label` | `{c1}_vs_{c0}` |

### Two contrast grammars (by `DE_test`)

| `DE_test` | Contrast style |
|-----------|----------------|
| `EdgeR`, `EdgeR-LRT`, `EdgeR-QLF`, `Dream` | **EdgeR-style expression** over design coefficient names |
| `DESeq2`, `DESeq2-LRT` | **DESeq2-style** semicolon triple |
| Wilcox / other single-cell | **Ignored** |

**Keep `c0`/`c1` on every row** for plots, replicate checks, descriptive pct, and labels — even when `contrast` targets a Batch main effect or an interaction coefficient.

Blank `contrast` auto-fills as Condition **c1 vs c0** in the active style (`Condition{c1} - Condition{c0}` for EdgeR*/Dream; `Condition;{c1};{c0}` for DESeq2*).

#### EdgeR-style expressions

- Prefer **`A - B` contrasts** for factor level comparisons (e.g. `ConditionTreatment - ConditionControl`, `BatchBatchB - BatchBatchA`).
- Under **treatment coding**, the reference level is absorbed in the intercept, so `ConditionTreatment - ConditionControl` resolves to the `ConditionTreatment` coefficient when `ConditionControl` is not a column.
- A **single coefficient name** is appropriate when testing that coef vs 0 — especially the **interaction term** (e.g. `ConditionTreatment:GenotypeKO`).
- Names must match `colnames(model.matrix(...))` (or Dream’s cell-means Condition coding). Interactions use `:` as in R’s design matrix.

#### DESeq2-style contrasts

- `Factor;numerator;denominator` is the comps string form of DESeq2’s length-3 contrast vector: `contrast = c(Factor, numerator, denominator)` passed to `DESeq2::results(...)`.

**`EdgeR` (exactTest)** only allows simple `~ Condition` with two groups — use **`EdgeR-LRT`** or **`EdgeR-QLF`** for covariates, >2 levels, or interactions. **`Dream`** requires `(1|var)` in `formula`; EdgeR/DESeq2 reject random effects. **`DESeq2-LRT`:** nested LRT per comps row (p-values are LRT; LFC follows the contrast).

**EdgeR-QLF vs EdgeR-LRT (same formulas and contrasts):**

| | EdgeR-QLF | EdgeR-LRT |
|--|-----------|-----------|
| Pros | Better type-I error control with overdispersion; modern edgeR default recommendation for many designs | Familiar LRT; slightly simpler / historically common in sc pipelines |
| Cons | Often fewer “significant” genes (more conservative) | Can be anti-conservative when quasi-dispersion is large |

Default remains `EdgeR-LRT`.

---

## 1. Single-cell 1v1 (Wilcox / FindMarkers)

**When to use:** Few or no biological replicates; exploratory cell-level DE.

**Compositional analysis:** `prop.test()` (not propeller).

**Contrast:** can be left blank (or omitted); it is **ignored** — DE uses `c0`/`c1` only.

### `sample_metadata.csv`

| Sample | Condition | Code |
|--------|-----------|------|
| Control_1 | Control | Control_1 |
| Treatment_1 | Treatment | Treatment_1 |

### `comps.csv`

| c0 | c1 | label |
|----|----|-------|
| Control | Treatment | Treatment_vs_Control |

### Runner

```r
scDAPP::scRNAseq_pipeline_runner(
  datadir = "path/to/datadir",
  outdir = "path/to/outdir",
  sample_metadata = "path/to/sample_metadata.csv",
  comps = "path/to/comps.csv",
  Pseudobulk_mode = FALSE,
  DE_test = "wilcox"
)
```

With `Pseudobulk_mode = FALSE`, `DE_test` may also be: `wilcox_limma`, `bimod`, `t`, `negbinom`, `poisson`, `LR`, `MAST`, `DESeq2` (Seurat `FindMarkers` tests; `roc` not supported).

---

## 2. Simple A vs B pseudobulk

**Contrast:** can be left blank; defaults to Condition **c1 vs c0** in the style of `DE_test`.

### `sample_metadata.csv`

| Sample | Condition | Code |
|--------|-----------|------|
| Control_1 | Control | Control_1 |
| Control_2 | Control | Control_2 |
| Treatment_1 | Treatment | Treatment_1 |
| Treatment_2 | Treatment | Treatment_2 |

### `comps.csv`

| c0 | c1 | label |
|----|----|-------|
| Control | Treatment | Treatment_vs_Control |

### Runner

```r
scDAPP::scRNAseq_pipeline_runner(
  datadir = "path/to/datadir",
  outdir = "path/to/outdir",
  sample_metadata = "path/to/sample_metadata.csv",
  comps = "path/to/comps.csv",
  Pseudobulk_mode = TRUE,
  DE_test = "EdgeR-LRT"   # or "EdgeR-QLF", "EdgeR", "DESeq2", "DESeq2-LRT"
)
```

| `DE_test` | Notes |
|-----------|--------|
| `"EdgeR-LRT"` | Default; GLM + per-contrast likelihood-ratio test |
| `"EdgeR-QLF"` | GLM + per-contrast quasi-likelihood F-test (same comps/contrasts as LRT) |
| `"EdgeR"` | exactTest for simple 2-group `~ Condition` only |
| `"DESeq2"` | Wald; needs Bioconductor `DESeq2` |
| `"DESeq2-LRT"` | Nested LRT per comps contrast; LRT p-values + contrast LFC |

Do **not** put `(1|var)` in `formula` with EdgeR/DESeq2 — use [Dream](#6-random-effects--paired-dream).

---

## 3. A vs B vs C (multiple pairwise comparisons)

**Contrast:** can be left blank on each row (defaults to that row’s Condition c1 vs c0).

### `sample_metadata.csv`

| Sample | Condition | Code |
|--------|-----------|------|
| Control_1 | Control | Control_1 |
| Control_2 | Control | Control_2 |
| Treatment_1 | Treatment | Treatment_1 |
| Treatment_2 | Treatment | Treatment_2 |
| Other_1 | Other | Other_1 |
| Other_2 | Other | Other_2 |

### `comps.csv`

| c0 | c1 | label |
|----|----|-------|
| Control | Treatment | Treatment_vs_Control |
| Control | Other | Other_vs_Control |
| Treatment | Other | Other_vs_Treatment |

Shared `~ Condition` fit; each row extracts its contrast (blank → Condition c1 vs c0 for that row).

---

## 4. Multivariable tests with confounders

**Contrast:** can be left blank for the default Condition **c1 vs c0** (adjusted for covariates in `formula`). Set `contrast` explicitly to extract **other terms** from the same formula (e.g. a Batch effect adjusted for Condition).

### `sample_metadata.csv`

| Sample | Condition | Code | Batch |
|--------|-----------|------|-------|
| Control_1 | Control | Control_1 | BatchA |
| Control_2 | Control | Control_2 | BatchB |
| Treatment_1 | Treatment | Treatment_1 | BatchA |
| Treatment_2 | Treatment | Treatment_2 | BatchB |

### `comps.csv` — EdgeR-LRT / EdgeR-QLF

| c0 | c1 | formula | contrast | label |
|----|----|---------|----------|-------|
| Control | Treatment | `~ Condition + Batch` | ConditionTreatment - ConditionControl | Treatment_vs_Control_batch_adj |
| Control | Treatment | `~ Condition + Batch` | BatchBatchB - BatchBatchA | BatchB_vs_BatchA_cond_adj |

### `comps.csv` — DESeq2 / DESeq2-LRT

| c0 | c1 | formula | contrast | label |
|----|----|---------|----------|-------|
| Control | Treatment | `~ Condition + Batch` | Condition;Treatment;Control | Treatment_vs_Control_batch_adj |
| Control | Treatment | `~ Condition + Batch` | Batch;BatchB;BatchA | BatchB_vs_BatchA_cond_adj |

Use **`EdgeR-LRT`**, **`EdgeR-QLF`**, or DESeq2* here — not `EdgeR` exactTest.

---

## 5. Interactions

Fit `~ Condition * Genotype` (equivalent to `Condition + Genotype + Condition:Genotype`). Report either the **main Condition effect** or the **interaction term** (often the scientific focus).

**Contrast:** can be left blank for the default Condition **c1 vs c0** main-effect contrast under the interaction formula. Set `contrast` explicitly to the **interaction coefficient** (EdgeR `FactorA:FactorB`). For DESeq2 interaction coefficients, use `name=ResultsName` (DESeq2’s single-coefficient form; names often use `.` instead of `:`).

### `sample_metadata.csv`

| Sample | Condition | Code | Genotype |
|--------|-----------|------|----------|
| Control_WT_1 | Control | Control_WT_1 | WT |
| Control_KO_1 | Control | Control_KO_1 | KO |
| Treatment_WT_1 | Treatment | Treatment_WT_1 | WT |
| Treatment_KO_1 | Treatment | Treatment_KO_1 | KO |
| Control_WT_2 | Control | Control_WT_2 | WT |
| Control_KO_2 | Control | Control_KO_2 | KO |
| Treatment_WT_2 | Treatment | Treatment_WT_2 | WT |
| Treatment_KO_2 | Treatment | Treatment_KO_2 | KO |

(Prefer ≥2 samples per Condition×Genotype cell so the interaction is estimable.)

### `comps.csv` — EdgeR-LRT / EdgeR-QLF / Dream (expressions)

| c0 | c1 | formula | contrast | label |
|----|----|---------|----------|-------|
| Control | Treatment | `~ Condition * Genotype` | ConditionTreatment - ConditionControl | Treatment_vs_Control_main_in_interaction_model |
| Control | Treatment | `~ Condition * Genotype` | ConditionTreatment:GenotypeKO | Treatment_x_Genotype_interaction_term |

### `comps.csv` — DESeq2 / DESeq2-LRT

| c0 | c1 | formula | contrast | label |
|----|----|---------|----------|-------|
| Control | Treatment | `~ Condition * Genotype` | Condition;Treatment;Control | Treatment_vs_Control_main_in_interaction_model |
| Control | Treatment | `~ Condition * Genotype` | name=ConditionTreatment.GenotypeKO | Treatment_x_Genotype_interaction_term |

Exact interaction coefficient names depend on factor levels and coding — check `colnames(model.matrix(...))` or `DESeq2::resultsNames(dds)` after a dry-run fit if unsure. Propeller uses the same contrast grammar (including single-coef interaction terms).

---

## 6. Random effects / paired (Dream)

**Requirements:** `DE_test = "Dream"`; at least one `(1|var)` in `formula`; Bioconductor `variancePartition`.

**Contrast:** can be left blank; defaults to EdgeR-style Condition **c1 vs c0** under Dream’s cell-means Condition coding.

### `sample_metadata.csv`

| Sample | Condition | Code | Patient |
|--------|-----------|------|---------|
| Control_P1 | Control | Control_P1 | P1 |
| Treatment_P1 | Treatment | Treatment_P1 | P1 |
| Control_P2 | Control | Control_P2 | P2 |
| Treatment_P2 | Treatment | Treatment_P2 | P2 |

### `comps.csv`

| c0 | c1 | formula | label |
|----|----|---------|-------|
| Control | Treatment | `~ Condition + (1\|Patient)` | Treatment_vs_Control_paired |

```r
scDAPP::scRNAseq_pipeline_runner(
  datadir = "path/to/datadir",
  outdir = "path/to/outdir",
  sample_metadata = "path/to/sample_metadata.csv",
  comps = "path/to/comps.csv",
  Pseudobulk_mode = TRUE,
  DE_test = "Dream",
  workernum = 4
)
```

Compositional propeller with `(1|var)` uses blocked limma independently of Dream DE. Interaction + pairing (e.g. `~ Condition * Genotype + (1|Patient)`) is allowed for Dream DE with EdgeR-style contrasts as in §5.

---

## Summary

| Design | `Pseudobulk_mode` | Typical `DE_test` | Typical `contrast` |
|--------|-------------------|-------------------|--------------------|
| Single-cell 1v1 | `FALSE` | `wilcox` | ignored (`c0`/`c1` only) |
| Simple A vs B | `TRUE` | `EdgeR-LRT` / `EdgeR-QLF` | blank or `Condition{c1} - Condition{c0}` |
| Simple A vs B | `TRUE` | `DESeq2` / `DESeq2-LRT` | blank or `Condition;{c1};{c0}` |
| A vs B vs C | `TRUE` | EdgeR* / DESeq2* | blank per row (Condition c1 vs c0) |
| Confounders | `TRUE` | `EdgeR-LRT` / `EdgeR-QLF` / DESeq2* | Condition or Batch `A - B` / triple |
| Interaction term | `TRUE` | `EdgeR-LRT` / `EdgeR-QLF` | `ConditionL1:GenotypeL2` |
| Interaction term | `TRUE` | DESeq2* | `name=ConditionL1.GenotypeL2` |
| Paired | `TRUE` | `Dream` | blank or EdgeR-style Condition c1 vs c0 |

---

## Future work / TODOs

- Document running comparative modules outside the full pipeline runner.
- Richer EdgeR expression parsing beyond single-coef and `A - B` if needed.
