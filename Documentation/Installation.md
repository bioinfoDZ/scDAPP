# Installation

If you run into any issues, you can check the [common bugs and fixes FAQ](CommonBugs.md). If the error is not reported there, please save the error message and open a Github Issue in this repository.

**Optional packages (v2.0):** [`conda_env/2026_scdappv2.yml`](../conda_env/2026_scdappv2.yml) includes the v2 feature packages that are published on conda. `RISC` and `DoubletFinder` are not on conda-forge or bioconda; `devtools::install_github()` still installs them. See [Optional packages (v2.0)](#optional-packages-v20) below and [Usage.md — Integration method dependencies](Usage.md#integration-method-dependencies).



### This guide

The general steps for installation described by this guide are:
 1. install non-R and some R dependency packages within a Conda virtual environment
 2. install R dependencies and the pipeline itself in R



If using a Linux server / HPC, please use the steps below.
If installing on Mac, the conda steps may be difficult - please see the bottom of this page.


<br />
<br />


### 1. Install Conda: 

https://docs.conda.io/en/latest/miniconda.html

You will likely need to restart the terminal for the installation to finish.

#### Using Conda on HPC:

After installation, you may need to "source" Conda each time you want to use it, for example in submission scripts. Adding a source command to something like your bash_rc file will load it in your login node, but may not work in submitted jobs.

Example:
```
#source the path to your conda.sh file
source <path/to/your/>conda.sh
```



#### Note for Einstein HPC users:

On Einstein HPC, I source my locally installed Conda like this. It will change for you depending on which folder you install Conda in:
```source /gs/gsfs0/home/aferrena/packages/miniconda3/miniconda3/etc/profile.d/conda.sh```


I do not recommend using the pre-installed conda on HPC, but rather starting from fresh. To use pre-installed conda, simply do not source the newly installed conda.
 




<br />



### 2. Create a conda environment using yaml file:
Navigate to yaml file: [conda_env/2026_scdappv2.yml](../conda_env/2026_scdappv2.yml)

Then download or copy and paste it into a file on the server / HPC. You can click the file and click on the "overlapping boxes" to copy raw file contents then paste it as a file in your location.

The file asks for R >= 4.5 and Seurat >= 5. Without those lower bounds, conda can select R 4.1 and Seurat 4, which scDAPP v2 cannot use.

[`conda_env/2025scdapp.yml`](../conda_env/2025scdapp.yml) remains for existing environments. It does not include the v2 feature packages listed below.


Create the environment using the .yml file:
```
conda env create -f 2026_scdappv2.yml 
```


It may ask your permission, just say yes.

It will take a while to install after you give the permissions.

If it still does not work, please write down the dependencies that cause trouble and the error messages, and open an issue in this repo.


<br />

### 3. Using R, install the package and dependencies in the conda environment:



Activate the conda environment, then within load R:
```
#load conda
conda activate 2026_scdappv2

#activate R
R

#to quit R, use: 
q('no') #remove hashtag in front to use this
```


####  In R, install scDAPP:


You can then install this package and its dependencies from Github with:

```
#install.packages("devtools")
devtools::install_github("bioinfoDZ/scDAPP", dependencies = TRUE)
```

It will ask your permission, just say yes.

It will take a while to install all the dependencies and this package.

Some packages (such as "XML") may already be present due to Conda, but cannot be updated in R. It is okay, you can ignore it.

It may not work with some packages causing errors. R packages can be installed manually from within R or from conda. If you used the conda virtual environment steps, many R packages are also available via conda (just google package name + "conda") which may be a lot easier than manual installation.

If you run into any issues, you can check the [common bugs and fixes FAQ](CommonBugs.md). If the error is not reported there, please save the error message and open a Github Issue in this repository.


####  Test installation:

Finally, if you did not get any errors during the package installation steps, you can test the installation by running in R:

```
scDAPP::r_package_test()
```

`r_package_test()` attaches all pipeline libraries (including RISC) via `attach_scDAPP_pipeline_libraries()` and returns package versions.

For scripts or interactive use without the version table:

```r
# Default pipeline (RISC integration)
scDAPP::attach_scDAPP_pipeline_libraries(load_risc = TRUE)

# Seurat-only integration smoke tests (no RISC attach)
scDAPP::attach_scDAPP_pipeline_libraries(load_risc = FALSE)
```

By default, `attach_scDAPP_pipeline_libraries()` also calls `set_parallel_blas_threads()` (`configure_parallel = TRUE`) to limit BLAS threading during multi-core work via `RhpcBLASctl::blas_set_num_threads(1)`.

If you created `2026_scdappv2`, **clusterProfiler** is already installed. Verify it with:

```r
requireNamespace("clusterProfiler", quietly = TRUE)
```

Bioconda currently ships the Bioconductor 3.22 build (clusterProfiler 4.18). Do not also run `BiocManager::install("clusterProfiler")` inside this env. R cannot upgrade a package that conda installed, the same way it cannot update XML (see above). The yaml does not pin an exact clusterProfiler version, so a later `conda env create` picks up a newer bioconda build when one is published.

If you are not using `2026_scdappv2`, install the feature packages in R after the conda env and `devtools::install_github()`; see [Optional packages (v2.0)](#optional-packages-v20).

You should see messages about packages being activated followed by a data.frame showing the key dependency packages and their versions such as below (version numbers do not need to match the example below, just make sure there are no errors when you run the command):
```
              pkg   vers
1       tidyverse  2.0.0
2          Seurat  5.0.2
3       patchwork  1.2.0
4        ggdendro  0.2.0
5         foreach  1.5.2
6         msigdbr  7.5.1
7      ggalluvial 0.12.5
8       ggfittext 0.10.2
9         ggrepel  0.9.5
10          hdf5r  1.3.9
11          edgeR 4.0.16
12      glmGamPoi 1.14.3
13          fgsea 1.28.0
14 ComplexHeatmap 2.18.0
15  DoubletFinder  2.0.4
16           RISC  1.7.0
17         scDAPP  2.0.0
```

You may need to exit R (via `q('no')`) and reload and rerun this if it does not work immediately following installation.


<br />

## Optional packages (v2.0)

`conda_env/2026_scdappv2.yml` installs these from conda-forge or bioconda:

| Package | When needed | Conda package |
|---------|-------------|----------------|
| `harmony` | `integration_method = "HarmonyIntegration"` | `r-harmony` |
| `mclust` | `pcs_int` or `res_int = "auto"` (cluster-stability ARI) | `r-mclust` |
| `clusterProfiler` | `run_msigdb_celltype_ora = TRUE` (default) or `run_ORA = TRUE` | `bioconductor-clusterprofiler` |
| `variancePartition` | `DE_test = "Dream"` | `bioconductor-variancepartition` |
| `miloR` | `run_MILO = TRUE` and `Pseudobulk_mode = TRUE` | `bioconductor-milor` |
| `SingleCellExperiment` | required by `miloR` | `bioconductor-singlecellexperiment` |
| `DESeq2` | `DE_test` is `"DESeq2"` or `"DESeq2-LRT"` | `bioconductor-deseq2` |
| `presto` | faster Wilcoxon tests when Seurat finds it | `r-presto` |

`bioconductor-milor` is published for `linux-64` and `osx-64`. There is no Apple Silicon (`osx-arm64`) build. This conda file is the Linux HPC path.

`RISC` and `DoubletFinder` are not on conda-forge or bioconda. `devtools::install_github("bioinfoDZ/scDAPP", dependencies = TRUE)` installs them from GitHub.

If you skip conda, or you are still using `conda_env/2025scdapp.yml`, install the feature packages in R when you need them:

| Package | Install |
|---------|---------|
| `harmony` | `install.packages("harmony")` |
| `mclust` | `install.packages("mclust")` |
| `clusterProfiler` | `if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("clusterProfiler")` |
| `variancePartition` | `BiocManager::install("variancePartition")` |
| `miloR` | `BiocManager::install("miloR")` |
| `DESeq2` | `BiocManager::install("DESeq2")` |

`presto` is a GitHub Suggest (`immunogenomics/presto`). `dependencies = TRUE` on the scDAPP install pulls it when that remote is available.

See [Usage.md](Usage.md) for parameter details and [Upgrading from v1.3](Usage.md#upgrading-from-v13) if migrating from scDAPP v1.3.x.


<br />
<br />



# If using Mac, start from here

If using Mac, the conda virtual environment steps will likely not work. Instead, just skip the conda virtual environment steps and go directly to R. As long as you can get R >= v4.0, Seurat >= v5.0, and RISC >= v1.7 working, you should not have problems running the pipeline.

If on Mac then very likely, what you will need to do is install Apple Xcode from the app store, then open a terminal and run the following:
- `xcode-select --install`
- `sudo xcodebuild -license accept`

You may also need to download the GNU Fortan compiler, accessible from this website: 
- https://mac.r-project.org/tools/

These steps will install key compiler tools that are not easy to install any other way. Then, you should be able to install R >= 4.0 for your system [from Cran](https://cran.r-project.org/), and finally within R install [Seurat using its instructions](https://satijalab.org/seurat/articles/install.html) and [RISC from github](https://github.com/bioinfoDZ/RISC).


Then, in R, you can then install scDAPP with:

```
devtools::install_github("bioinfoDZ/scDAPP", dependencies = TRUE)
```

<br />
<br />


# Development branch

The [`dev`](https://github.com/bioinfoDZ/scDAPP/tree/dev) branch currently carries **scDAPP v2.0.0** (integration backends, cluster-stability auto-tuning, RISC autoV2 reference selection, MSigDB cache, cross-condition designs, Milo, and `subset_module()`). See the [v2.0 Changelog](Changelog.md) and [Usage — Upgrading from v1.3](Usage.md#upgrading-from-v13).

Install from GitHub:

```
devtools::install_github("bioinfoDZ/scDAPP@dev", dependencies = TRUE)
```

The dev branch may not be 100% stable and can be subject to frequent updates. If you need a known-stable v1.3.x build, install a tagged release instead (see [Releases](https://github.com/bioinfoDZ/scDAPP/releases)).




