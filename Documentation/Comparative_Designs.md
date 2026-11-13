# Comparative designs: `sample_metadata`, `comps`, and `scRNAseq_pipeline_runner()`

This cookbook shows how to set up cross-condition comparisons. Parameter lists are in [Usage.md](Usage.md).

## The two files

`sample_metadata.csv` has one row per sample. `Sample` matches the Cell Ranger folder or the Seurat file name. `Condition` is the group. `Code` is an optional short name; if you omit it, the pipeline uses `Condition_Sample`. Add columns such as `Batch` or `Patient` only when a comparison uses them.

`comps.csv` has one row per comparison.

- `c0` is the reference. `c1` is the test. Positive log2 fold change means higher in **c1**.
- `formula`, `contrast`, and `label` are optional. Leave them out for a plain A versus B test. A blank formula is `~ Condition`. A blank contrast is Condition c1 versus c0.
- `label` names the comparison in the report and in output folders. If you leave it blank, the name is `{c1}_vs_{c0}` (for example `Treatment_vs_Control`). Set `label` when two rows share the same `c0` and `c1`. Those rows would otherwise write the same folder under pathways, composition, and Milo, and the second row would overwrite the first. Section 4 is the example.
- Legacy files with `c2` still run: `c2` is renamed to `c0`.

`formula` and `contrast` are used by pseudobulk DE and by propeller (`Pseudobulk_mode = TRUE`). They are ignored by Wilcoxon and by `prop.test()`.

`Condition` levels follow first appearance in `sample_metadata`. Put Control (or whichever group should be the intercept) in the first rows. Other covariate columns read from the CSV become factors in alphabetical order. Pairwise `c0` and `c1` still name the two levels you want, whatever the intercept is.

## Characters

Letters, digits, and `_` are safe in every column. Hyphens in `Sample` and `Code` are fine. Those values are file and folder names, and Cell Ranger names often look like `S00033-Ja005E`.

These characters cause trouble:

- `/` or `\` in `Sample`, `Code`, or `label` breaks the path. Pathway and composition folders use `label` as written. Milo replaces any character outside `A-Z`, `a-z`, `0-9`, `.`, `_`, and `-` with `_`.
- `-` in a condition or covariate **level** breaks EdgeR contrasts. The contrast parser splits on `-`, so a level `COVID-SEV` is read as subtraction. DESeq2 rewrites `-` to `_` in levels, so the contrast has to use the underscore form. Write `COVID_SEV`.
- `:` in a level collides with interaction names such as `ConditionTreatment:GenotypeKO`.
- Spaces, `+`, `*`, and parentheses in a level or a column name are formula syntax.
- A column in `(1|var)` must match `[A-Za-z.][A-Za-z0-9_.]*`. Otherwise the random effect is not recognized.
- An unquoted comma or a line break inside a CSV cell splits the row. A file that does not end with a newline (typical of Excel) is still read, with a warning.

---

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

---

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

The comparison is named `Treatment_vs_Control`.

---

## 3. A vs B vs C

**When to use:** Three or more conditions with replicates. Each row is one pairwise test. Leave `contrast` blank on each row.

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

Each row has a different `c0`/`c1` pair, so the automatic names (`Treatment_vs_Control`, `Other_vs_Control`, `Other_vs_Treatment`) do not collide.

---

## 4. Confounders

**When to use:** The condition effect should be adjusted for a covariate such as batch, or you want a covariate contrast from the same model.

Add the covariate to `sample_metadata`. Put it in `formula`. Leave `contrast` blank for Condition c1 versus c0 adjusted for that covariate. Set `contrast` to pull out another term, such as batch adjusted for condition.

Both rows below are Control versus Treatment, so each row needs its own `label`.

### `sample_metadata.csv`

| Sample      | Condition | Code        | Batch  |
| ----------- | --------- | ----------- | ------ |
| Control_1   | Control   | Control_1   | BatchA |
| Control_2   | Control   | Control_2   | BatchB |
| Treatment_1 | Treatment | Treatment_1 | BatchA |
| Treatment_2 | Treatment | Treatment_2 | BatchB |

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

---

## 5. Interactions

**When to use:** You want the condition effect to differ by a second factor. Prefer at least two samples in each Condition by factor cell.

Fit `~ Condition * Genotype` (that is `Condition + Genotype + Condition:Genotype`). A blank contrast is the Condition c1 versus c0 main effect. Set `contrast` to the interaction coefficient when that term is the result you want. The same `label` rule as section 4 applies: both rows are Control versus Treatment.

A main-effect contrast still reports `PropMean_c1` and `PropMean_c0`. An interaction-coefficient contrast reports the coefficient test and leaves those mean columns empty.

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

| c0      | c1        | formula                  | contrast                              | label                                 |
| ------- | --------- | ------------------------ | ------------------------------------- | ------------------------------------- |
| Control | Treatment | `~ Condition * Genotype` | ConditionTreatment - ConditionControl | Treatment_vs_Control_main             |
| Control | Treatment | `~ Condition * Genotype` | ConditionTreatment:GenotypeKO         | Treatment_x_Genotype                  |

### `comps.csv` — DESeq2 / DESeq2-LRT

| c0      | c1        | formula                  | contrast                           | label                     |
| ------- | --------- | ------------------------ | ---------------------------------- | ------------------------- |
| Control | Treatment | `~ Condition * Genotype` | Condition;Treatment;Control        | Treatment_vs_Control_main |
| Control | Treatment | `~ Condition * Genotype` | name=ConditionTreatment.GenotypeKO | Treatment_x_Genotype      |

Coefficient names depend on factor levels. Check `colnames(model.matrix(...))` or `DESeq2::resultsNames(dds)` if a contrast does not match.

---

## 6. Paired samples (Dream)

**When to use:** The same patient or animal is in both conditions.

`sample_metadata` needs the block column (`Patient` here). `formula` needs `(1|Patient)`. `DE_test` must be `"Dream"` (Bioconductor `variancePartition`). EdgeR and DESeq2 reject random effects. Propeller uses the same `(1|var)` term for a blocked test, whether or not DE is Dream.

Leave `contrast` blank for Condition c1 versus c0. One row does not need `label`; the name is `Treatment_vs_Control`.

### `sample_metadata.csv`

| Sample       | Condition | Code         | Patient |
| ------------ | --------- | ------------ | ------- |
| Control_P1   | Control   | Control_P1   | P1      |
| Treatment_P1 | Treatment | Treatment_P1 | P1      |
| Control_P2   | Control   | Control_P2   | P2      |
| Treatment_P2 | Treatment | Treatment_P2 | P2      |

### `comps.csv`

```csv
c0,c1,formula
Control,Treatment,~ Condition + (1|Patient)
```

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

`~ Condition * Genotype + (1|Patient)` is allowed. Use EdgeR-style contrasts as in section 5, and set `label` if more than one row shares `c0` and `c1`.

---

## 7. Neighbourhood abundance (Milo)

**When to use:** `run_MILO = TRUE` and `Pseudobulk_mode = TRUE`. Milo tests overlapping neighborhoods on the integrated embedding. It uses the same `comps` rows as pseudobulk DE. A blank contrast is still Condition c1 versus c0. Set `contrast` to the interaction coefficient when that term is the result you want.

**Install:** Bioconductor `miloR`.

**k:** `milo_k = "auto"` uses `min(5 * n_samples, milo_k_max, n_cells - 1)` with `milo_k_max = 100`. An integer `milo_k` replaces that formula.

**Significance:** SpatialFDR below `milo_spatial_fdr` (default 0.1). The report prints the top 20 significant neighborhoods, plus an alluvial plot and a two-way heatmap of those 20 against the integrated clusters when that overlap table is non-empty. Tables are under `multisample_integration/differential_abundance_Milo/`.

**Paired formulas:** a `(1|var)` term is fit with Milo's mixed model and `force = TRUE`. That model only returns the last fixed coefficient. A contrast that is not that coefficient is fit without the random effect, and the report says so. If the mixed model errors, the fit is repeated with `(1|var)` removed.

```r
scDAPP::scRNAseq_pipeline_runner(
  datadir = "path/to/datadir",
  outdir = "path/to/outdir",
  sample_metadata = "path/to/sample_metadata.csv",
  comps = "path/to/comps.csv",
  Pseudobulk_mode = TRUE,
  DE_test = "EdgeR-LRT",
  run_MILO = TRUE
)
```

---

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

Preflight checks formulas and contrasts before the long run. A cluster that cannot estimate a DE contrast is skipped with a warning. A contrast that cannot be mapped onto the propeller design stops the run.

---

## Summary

| Design | `Pseudobulk_mode` | Typical `DE_test` | Typical `contrast` |
| --- | --- | --- | --- |
| Single-cell 1v1 | `FALSE` | `wilcox` | ignored (`c0` / `c1` only) |
| Simple A vs B | `TRUE` | `EdgeR-LRT` / `EdgeR-QLF` | blank, or `Condition{c1} - Condition{c0}` |
| Simple A vs B | `TRUE` | `DESeq2` / `DESeq2-LRT` | blank, or `Condition;{c1};{c0}` |
| A vs B vs C | `TRUE` | EdgeR* / DESeq2* | blank per row |
| Confounders | `TRUE` | `EdgeR-LRT` / `EdgeR-QLF` / DESeq2* | Condition or Batch `A - B`, or a DESeq2 triple. Set `label`. |
| Interaction term | `TRUE` | `EdgeR-LRT` / `EdgeR-QLF` | `ConditionL1:GenotypeL2`. Set `label`. |
| Interaction term | `TRUE` | DESeq2* | `name=ConditionL1.GenotypeL2`. Set `label`. |
| Paired | `TRUE` | `Dream` | blank, or an EdgeR-style Condition contrast |
| Milo | `TRUE` | same test as DE | same contrast; `run_MILO = TRUE` |
