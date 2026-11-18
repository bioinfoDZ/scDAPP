# Comparative designs: `sample_metadata`, `comps`, and `scRNAseq_pipeline_runner()`

This cookbook shows how to set up cross-condition comparisons. Parameter lists are in [Usage.md](Usage.md).

Comparative analyses in scDAPP include: 
- Differential expression (DE) testing via pseudobulk methods with replicates or Seurat FindMarkers wthout
- Compositional analysis (comparing proportions) via Propeller or 2-proportion Z-tests
- Optionally (for replicate-aware analyses only), differential neighborhood composition testing via Milo.

## The two metadata files

`sample_metadata.csv` has one row per sample. `Sample` matches the Cell Ranger folder or the Seurat file name. `Condition` is the group. `Code` is an optional short name; if you omit it, the pipeline uses `Condition_Sample`. Add columns such as `Batch` or `Patient` only when a comparison uses them.

`comps.csv` has one row per comparison.

- `c0` is the reference. `c1` is the test. Positive log2 fold change means higher in **c1**.
- `formula`, `contrast`, and `label` are optional. Leave them out for a plain A versus B test. A blank formula is `~ Condition`. A blank contrast is Condition c1 versus c0.
- Legacy files with `c2` still run: `c2` is renamed to `c0`.

`formula` and `contrast` are used by pseudobulk DE and by propeller (`Pseudobulk_mode = TRUE`). They are ignored by Wilcoxon and by `prop.test()`.

`Condition` levels follow first appearance in `sample_metadata`. Put Control (or whichever group should be the intercept) in the first rows. Other covariate columns read from the CSV become factors in alphabetical order. Pairwise `c0` and `c1` still name the two levels you want, whatever the intercept is.

<br />

## Preflight tests

Before a long run, you can check the comparison with `preflight_comparisons()`. Pass `sample_metadata` and `comps` as data frames or CSV paths, plus `Pseudobulk_mode` and `DE_test`. It runs the same checks the pipeline runs before integration: the allowed test, installed packages, labels, formulas, contrasts, and the propeller design when pseudobulk mode is on. The first problem stops the function with the pipeline's message. A passing check returns the comps table with labels and default contrasts filled in. A cluster that is too sparse to estimate a contrast can still be skipped without errors later, during differential expression testing.

```r
comps <- scDAPP::preflight_comparisons(
  sample_metadata = "path/to/sample_metadata.csv",
  comps = "path/to/comps.csv",
  Pseudobulk_mode = TRUE,
  DE_test = "EdgeR-LRT"
)
```


<br />

## 1. Single-cell 1v1 (Wilcox / FindMarkers)

**When to use:** Few or no biological replicates.

**Compositional analysis:** `prop.test()`, not propeller. The contrast is ignored. DE uses `c0` and `c1` only.

### `sample_metadata.csv`

| Sample      | Condition | Code        |
| ----------- | --------- | ----------- |
| Control_1   | Control   | Control_1   |
| Treatment_1 | Treatment | Treatment_1 |

### `comps.csv`

| c0      | c1        |
| ------- | --------- |
| Control | Treatment |

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

`DE_test` may also be `wilcox_limma`, `bimod`, `t`, `negbinom`, `poisson`, `LR`, `MAST`, or `DESeq2` (Seurat `FindMarkers` tests). `roc` is not supported.

<br />


## 2. Simple A vs B pseudobulk

**When to use:** Biological replicates in each condition.

**Compositional analysis:** Propeller, Condition c1 versus c0. Leave `contrast` blank.

### `sample_metadata.csv`

| Sample      | Condition | Code        |
| ----------- | --------- | ----------- |
| Control_1   | Control   | Control_1   |
| Control_2   | Control   | Control_2   |
| Treatment_1 | Treatment | Treatment_1 |
| Treatment_2 | Treatment | Treatment_2 |

### `comps.csv`

| c0      | c1        |
| ------- | --------- |
| Control | Treatment |

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

The comparison is named `Treatment_vs_Control`. A formula other than `~ Condition`, or a contrast you set, is added to that name.

<br />


## 3. A vs B vs C

**When to use:** Three or more conditions with replicates. Each row is one pairwise test. Leave `contrast` blank on each row.

Each row of comps denotes a specific comparison. As of v2.0, all samples are normalized together, and gene dispersions are estimated jointly.

### `sample_metadata.csv`

| Sample      | Condition | Code        |
| ----------- | --------- | ----------- |
| Control_1   | Control   | Control_1   |
| Control_2   | Control   | Control_2   |
| Treatment_1 | Treatment | Treatment_1 |
| Treatment_2 | Treatment | Treatment_2 |
| Other_1     | Other     | Other_1     |
| Other_2     | Other     | Other_2     |

### `comps.csv`

| c0        | c1        |
| --------- | --------- |
| Control   | Treatment |
| Control   | Other     |
| Treatment | Other     |

Each row has a different `c0`/`c1` pair, so the automatic names (`Treatment_vs_Control`, `Other_vs_Control`, `Other_vs_Treatment`) do not collide. A formula other than `~ Condition`, or a contrast you set, is also part of the automatic name, so those rows do not collide either.

<br />


## 4. Confounders

**When to use:** The condition effect should be adjusted for a covariate such as batch, or you want a covariate contrast from the same model.

Add the covariate to `sample_metadata`. Put it in `formula`. Leave `contrast` blank for Condition c1 versus c0 adjusted for that covariate. Set `contrast` to pull out another term, such as batch adjusted for condition.



### `sample_metadata.csv`

| Sample      | Condition | Code        | Batch  |
| ----------- | --------- | ----------- | ------ |
| Control_1   | Control   | Control_1   | BatchA |
| Control_2   | Control   | Control_2   | BatchB |
| Treatment_1 | Treatment | Treatment_1 | BatchA |
| Treatment_2 | Treatment | Treatment_2 | BatchB |

Below the rows show contrasts passed according to EdgeR-like or DESeq2-like notation. The first row in each gets you the condition comparison, adjusted by batch (most common); the second row is how you could get the Batch A vs Batch B comparison, adjusted by condition (less common). The `label` values in the tables are optional.


### `comps.csv` — EdgeR-LRT / EdgeR-QLF

| c0      | c1        | formula               | contrast                              | label                          |
| ------- | --------- | --------------------- | ------------------------------------- | ------------------------------ |
| Control | Treatment | `~ Condition + Batch` | ConditionTreatment - ConditionControl | Treatment_vs_Control_batch_adj |
| Control | Treatment | `~ Condition + Batch` | BatchBatchB - BatchBatchA             | BatchB_vs_BatchA_cond_adj      |

### `comps.csv` — DESeq2 / DESeq2-LRT

| c0      | c1        | formula               | contrast                    | label                          |
| ------- | --------- | --------------------- | --------------------------- | ------------------------------ |
| Control | Treatment | `~ Condition + Batch` | Condition;Treatment;Control | Treatment_vs_Control_batch_adj |
| Control | Treatment | `~ Condition + Batch` | Batch;BatchB;BatchA         | BatchB_vs_BatchA_cond_adj      |

Use `EdgeR-LRT`, `EdgeR-QLF`, or a DESeq2 test here. `EdgeR` exactTest only fits a plain `~ Condition` with two groups.

<br />


## 5. Interactions

**When to use:** the Condition effect differs by a second factor. Prefer at least two samples in each Condition by factor cell.

Fit `~ Condition * Genotype` (that is `Condition + Genotype + Condition:Genotype`). A blank contrast is the Condition c1 versus c0 main effect. Set `contrast` to the interaction coefficient when that term is the result you want.

A main-effect contrast still reports `PropMean_c1` and `PropMean_c0` in the Propeller compositonal analysis. An interaction-coefficient contrast reports the coefficient test and leaves those mean columns empty.

### `sample_metadata.csv`

| Sample         | Condition | Code           | Genotype |
| -------------- | --------- | -------------- | -------- |
| Control_WT_1   | Control   | Control_WT_1   | WT       |
| Control_KO_1   | Control   | Control_KO_1   | KO       |
| Treatment_WT_1 | Treatment | Treatment_WT_1 | WT       |
| Treatment_KO_1 | Treatment | Treatment_KO_1 | KO       |
| Control_WT_2   | Control   | Control_WT_2   | WT       |
| Control_KO_2   | Control   | Control_KO_2   | KO       |
| Treatment_WT_2 | Treatment | Treatment_WT_2 | WT       |
| Treatment_KO_2 | Treatment | Treatment_KO_2 | KO       |

### `comps.csv` — EdgeR-LRT / EdgeR-QLF / Dream

Below the rows show contrasts passed according to EdgeR-like or DESeq2-like notation. 
- The first row in each gets you the condition "main effect"
- The second row gets you the genotype "main effect"
- The third row is how you could get the interaction coeffients per-gene.


The `label` values in the tables are optional.

| c0      | c1        | formula                  | contrast                              | label                                 |
| ------- | --------- | ------------------------ | ------------------------------------- | ------------------------------------- |
| Control | Treatment | `~ Condition * Genotype` | ConditionTreatment - ConditionControl | Treatment_vs_Control_main             |
| Control | Treatment | `~ Condition * Genotype` | GenotypeKO - GenotypeWT               | KO_vs_WT_main                         |
| Control | Treatment | `~ Condition * Genotype` | ConditionTreatment:GenotypeKO         | Treatment_x_Genotype                  |

### `comps.csv` — DESeq2 / DESeq2-LRT

| c0      | c1        | formula                  | contrast                           | label                     |
| ------- | --------- | ------------------------ | ---------------------------------- | ------------------------- |
| Control | Treatment | `~ Condition * Genotype` | Condition;Treatment;Control        | Treatment_vs_Control_main |
| Control | Treatment | `~ Condition * Genotype` | Genotype;KO;WT                     | KO_vs_WT_main             |
| Control | Treatment | `~ Condition * Genotype` | name=ConditionTreatment.GenotypeKO | Treatment_x_Genotype      |

Coefficient names depend on factor levels. Check `colnames(model.matrix(...))` or `DESeq2::resultsNames(dds)` if a contrast does not match.

<br />


## 6. Paired samples (Dream)

**When to use:** The same patient or animal is in both conditions; repeated sampling, time-course studies

`sample_metadata` needs the block column (`Patient` here). `formula` needs `(1|Patient)`. `DE_test` must be `"Dream"` (Bioconductor `variancePartition`). Propeller uses the same `(1|var)` term for a blocked test


### `sample_metadata.csv`

| Sample       | Condition | Code         | Patient |
| ------------ | --------- | ------------ | ------- |
| Control_P1   | Control   | Control_P1   | P1      |
| Treatment_P1 | Treatment | Treatment_P1 | P1      |
| Control_P2   | Control   | Control_P2   | P2      |
| Treatment_P2 | Treatment | Treatment_P2 | P2      |

### `comps.csv`

| c0      | c1        | formula                        |
| ------- | --------- | ------------------------------ |
| Control | Treatment | `~ Condition + (1\|Patient)`   |

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

`~ Condition * Genotype + (1|Patient)` is also allowed. Use EdgeR-style contrasts as in section 5. Below the comps rows show contrasts passed according to EdgeR-like notation, but with pairing, which can only be done with Dream. 
- The first row in each gets you the condition "main effect"
- The second row gets you the genotype "main effect"
- The third row is how you could get the interaction coeffients per-gene.

### `sample_metadata.csv`

| Sample          | Condition | Code            | Genotype | Patient |
| --------------- | --------- | --------------- | -------- | ------- |
| Control_WT_P1   | Control   | Control_WT_P1   | WT       | P1      |
| Treatment_WT_P1 | Treatment | Treatment_WT_P1 | WT       | P1      |
| Control_WT_P2   | Control   | Control_WT_P2   | WT       | P2      |
| Treatment_WT_P2 | Treatment | Treatment_WT_P2 | WT       | P2      |
| Control_KO_P3   | Control   | Control_KO_P3   | KO       | P3      |
| Treatment_KO_P3 | Treatment | Treatment_KO_P3 | KO       | P3      |
| Control_KO_P4   | Control   | Control_KO_P4   | KO       | P4      |
| Treatment_KO_P4 | Treatment | Treatment_KO_P4 | KO       | P4      |

### `comps.csv`

| c0      | c1        | formula                                 | contrast                              | label                     |
| ------- | --------- | --------------------------------------- | ------------------------------------- | ------------------------- |
| Control | Treatment | `~ Condition * Genotype + (1\|Patient)` | ConditionTreatment - ConditionControl | Treatment_vs_Control_main |
| Control | Treatment | `~ Condition * Genotype + (1\|Patient)` | GenotypeKO - GenotypeWT               | KO_vs_WT_main             |
| Control | Treatment | `~ Condition * Genotype + (1\|Patient)` | ConditionTreatment:GenotypeKO         | Treatment_x_Genotype      |

<br />
<br />
<br />


## `DE_test` and contrast grammar

All samples are fit together, by cluster. Each comps row extracts one contrast from that fit.

| `DE_test` | Notes |
| --- | --- |
| `"EdgeR-LRT"` | Default when `Pseudobulk_mode = TRUE`. GLM plus a per-contrast likelihood-ratio test. |
| `"EdgeR-QLF"` | Same formulas and contrasts as LRT, with a quasi-likelihood F-test. Often more conservative. |
| `"EdgeR"` | exactTest. Simple two-group `~ Condition` only. |
| `"DESeq2"` | Wald test. Needs Bioconductor `DESeq2`. |
| `"DESeq2-LRT"` | Nested likelihood-ratio test per contrast. |
| `"Dream"` | Paired model. Requires `(1|var)` in `formula`. EdgeR-style contrasts. |
| `wilcox` and other Seurat `FindMarkers` tests | Used when `Pseudobulk_mode = FALSE`. Contrast is ignored. `roc` is not supported. |

| `DE_test` | How to write `contrast` |
| --- | --- |
| `EdgeR`, `EdgeR-LRT`, `EdgeR-QLF`, `Dream` | A coefficient expression: `ConditionTreatment - ConditionControl`, or a single name such as `ConditionTreatment:GenotypeKO`. |
| `DESeq2`, `DESeq2-LRT` | `Factor;numerator;denominator`, or `name=ResultsName` for one coefficient. Interaction names usually use `.` instead of `:`. |
| Wilcoxon and other single-cell tests | Ignored. |

`Factor;numerator;denominator` is the string form of `contrast = c(Factor, numerator, denominator)` passed to `DESeq2::results()`.

Preflight checks formulas, contrasts, and that labels are unique before the long run. A cluster that cannot estimate a DE contrast is skipped with a warning. A contrast that cannot be mapped onto the propeller design stops the run.

<br />
<br />

## Summary

| Design | `Pseudobulk_mode` | Typical `DE_test` | Typical `contrast` |
| --- | --- | --- | --- |
| Single-cell 1v1 | `FALSE` | `wilcox` | ignored (`c0` / `c1` only) |
| Simple A vs B | `TRUE` | `EdgeR-LRT` / `EdgeR-QLF` | blank, or `Condition{c1} - Condition{c0}` |
| Simple A vs B | `TRUE` | `DESeq2` / `DESeq2-LRT` | blank, or `Condition;{c1};{c0}` |
| A vs B vs C | `TRUE` | EdgeR* / DESeq2* | blank per row |
| Confounders | `TRUE` | `EdgeR-LRT` / `EdgeR-QLF` / DESeq2* | Condition or Batch `A - B`, or a DESeq2 triple. Label optional; duplicates stop the run. |
| Interaction term | `TRUE` | `EdgeR-LRT` / `EdgeR-QLF` | `ConditionL1:GenotypeL2`. Label optional; duplicates stop the run. |
| Interaction term | `TRUE` | DESeq2* | `name=ConditionL1.GenotypeL2`. Label optional; duplicates stop the run. |
| Paired | `TRUE` | `Dream` | blank, or an EdgeR-style Condition contrast |
| Milo | `TRUE` | same test as DE | same contrast; `run_MILO = TRUE` |
