# scDAPP single-cell RNA-seq analysis pipeline usage

This pipeline will perform individual QC and clustering, label transfer from a reference scRNAseq dataset to guess cell types (optional), integration with RISC, and cross-condition compositional (cell proportion) analysis and differetential experession (DE) analysis. Multi-condition comparison (condition A vs B vs C, ie WT vs KO vs Drug) are supported. 

If multiple replicates are present (ie WT 1 and WT2 vs KO1 and KO2), this pipeline can make use of pseudobulk methods for compositional analysis (Propeller) and DE (EdgeR-LRT). If no replicates are present, will use chi-square test of proportions (R proportion.test) for compositional analysis and Wilcoxon test for DE.

<br />

# Usage


If you have any issues, please email BOTH alexander.ferrena@einsteinmed.edu and alexanderferrena@gmail.com

Or open an issue in this github repo.


Minimally, this pipeline needs three inputs: the raw UMI counts data in .h5 files or Seurat objects, a file called `sample_metadata.csv` that contains info about the samples, and a file called `comps.csv` that tells the pipeline which cross-condition comparison to perform.


<br />

# Prep the required files

### 1. `datadir`: raw UMI count matrices directly from Cellranger outputs, or a folder of Seurat objects

#### Option 1: Run directly on Cellranger output (ie folders containing .h5 files)

Run Cellranger and keep the output folders from all samples together in a single folder.
The parameter `datadir` is the path to a folder containing the Cellranger outputs.

Cellranger produces many output files for each sample. Minimally, the folders in `datadir` must contain one item, the file called `filtered_feature_bc_matrix.h5`. For example, if you have four samples, you will need (or if you have run Cellranger, already have) a folder with four sub-folders with the sample names. `datadir` should be the path pointing to folder holding everything. It will search the subfolders for the `filtered_feature_bc_matrix.h5` file. The sample names (sample sub-folder names) should match the "Sample" column of the `sample_metadata.csv` file as explained below.



#### Option 2: Inputting Seurat Objects: `input_seurat_obj` = T

Alternatively the pipeline accepts Seurat objects if `input_seurat_obj` is set to TRUE. This can be useful for hashed / multiplexed samples, pre-filtered samples, or published data for which .h5 files are not easily available.

Save each Seurat object as individual .rds files using the `saveRDS()` R command. Each sample should be named "SampleXYZ1.rds" and so on. The sample names should match the "Sample" column of the `sample_metadata.csv` file as explained below. The pipeline will use the "RNA" assay and the 

For [hashed](https://cite-seq.com/cell-hashing/) or [multiplexed](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-3p-multi) inputs, this option can also be used. We recommended to split all of the samples apart into separate Seurat objects even if they are from the same hash / CMO pool.


<br />

### 2. `sample_metadata`: sample information

The parameter `sample_metadata` is a path to a .csv file that looks like this:



```
Sample,Condition,Code
SampleXYZ1,Control,Control1
SampleXYZ2,Control,Control2
SampleABC1,KO1,KO1_1
SampleABC2,KO1,KO1_2
SampleJKL1,KO2,KO2_1
SampleJKL2,KO2,KO2_1
```

The Sample column here must match the folder names in `datadir` if `input_seurat_obj` is set to FALSE (Cellranger input). For example, `datadir` is a path to a folder with some sub-folders called SampleXYZ1, SampleXYZ2, etc. Or, they must match the names of the Seurat objects if `input_seurat_obj`. If so, `datadir` will instead be a path to a folder containing Seurat objects saved like SampleXYZ1.rds, SampleXYZ2.rds, etc.

The Code column is optional and can be used to give nicknames to samples if the sample names are ugly. If not provided, "Code" will be set to "Sample_Condition", ie SampleXYZ1_Control in the example above.

The order of the samples in this file will determine the plotting order in the report.

You can use the bash/zsh command `nano` to quickly create and save this file as well as the `comps.csv` file if working in a unix shell context.


<br />


### 3. `comps`: tell the pipeline which cross-condition comparisons to check

The parameter `comps` leads to a .csv file that looks like this:

```
c1,c2
KO1,Control
KO2,Control
KO1,KO2
```

This is used to tell the pipeline which conditions to compare. It also allows multiple cross-condition comparisons to be tested. Each row sets up an "A vs B" comparison, with A on the left (c1 column) and B on the right (c2 column). 
In the differential expression results, genes overexpressed in the A condition will have a positive log fold change, while genes in the B direction will have a negative log fold change and are said to be "underexpressed in A" by convention. Additionally, in the gene set enrichment analysis results, pathways enriched in A will have a positive normalized enrichment score and vice versa.

If comparing, for example, Disease vs Healthy, or KO vs WT, or Treated vs Untreated, it is convention to set the experiment group as A (Disease, or KO, or Treated) and the control as B (Healthy, WT, Untreated). In other words, the group of interest is denoted with positive fold changes, while the control group is denoted with negative fold changes.








<br />
<br />



# Running the pipeline

We will invoke the pipeline from the unix terminal as below.

### 1. first create a file called pipeline_runner.R containing the following:

```
#test packages
scDAPP::r_package_test()


#run pipeline with options

#run pipeline with options
scDAPP::scRNAseq_pipeline_runner(
               datadir = 'path/to/cellranger/outputs',
               outdir = 'path/to/output/folder',
               sample_metadata = 'path/to/sample_metadata.csv'
               comps = 'path/to/comps.csv',
               Pseudobulk_mode = T, #set to F if not replicates

               use_labeltransfer = F,
               refdatapath = 'path/to/reference/Seuratobject_SCTnormalized.rds',
               m_reference = 'path/to/reference/reference_FindAllMarkers.rds',

               species = 'Mus musculus',

               workernum = 1,
               input_seurat_obj = F
               )

```

### 2. Next, execute the R file from unix terminal (bash, zsh, etc) with the following commands:

```
#run the pipeline runner script
R CMD BATCH --no-save --no-restore pipeline_runner.R
```

It will produce a log file called `pipeline_runner.Rout`, which you can monitor, for example via the command `tail -f pipeline_runner.Rout`


<br />

### Running on HPC systems

On HPC, it works by submitting a job, which runs the R script “pipeline_runner.R”, which in turn runs the .Rmd file.


Here is an example HPC submission script:

```
#!/bin/bash
#SBATCH -p normal
#SBATCH --job-name=SDAP
#SBATCH -N 1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=150gb
#SBATCH -t 48:00:00
#SBATCH -o /path/to/job/report/%x-%A_%a.out



#my conda
# you may need to install your own local conda on HPC
# https://docs.conda.io/en/latest/miniconda.html
source /gs/gsfs0/home/aferrena/packages/miniconda3/miniconda3/etc/profile.d/conda.sh

#activate r_env
# see detailled install instructions for this conda env
conda activate r_env

#run R file 
# assumes there is the file "pipeline_runner.sh" in the current working directory
rfile=pipeline_runner.R

echo Submitting $rfile

R CMD BATCH --no-save --no-restore $rfile

conda deactivate


printf '\n\n\nAll Done!!!\n\n\n'


echo $SLURM_JOB_NAME
echo $SLURM_JOB_ID
```


### Running on powerful servers

If your computer has high memory and can keep running for hours, you can execute from Unix command line (bash, zsh) like so:

```
nohup R CMD BATCH --no-save --no-restore pipeline_runner.R &
```


<br />
<br />

## Key options: Cell type prediction with Label Transfer, and manual RISC reference selection

### 1. Cell type prediction with Label Transfer

Label transfer from a single-cell RNA-seq dataset to guess the cell types detected in your data is supported in the pipeline as an option. To use it, you need two files:

- `refdatapath` string, path to a Seurat object .rds file, pre-processed with `Seurat::SCTransform()`, with a column called "Celltype" in its meta.data. Ignored if `use_labeltransfer = F`
- `m_reference` string, path to .rds file containing output of `Seurat::FindAllMarkers` run on the reference object specified above. Ignored if `use_labeltransfer = F`

Once you have these, set `use_labeltransfer = T` and provide the paths with the parameters above.

To select a label transfer reference, it is recommended to use a tissue as similar as possible, from the same species. There are "atlas" databases such as Tabula Muris or Human Cell Atlas that may have the tissue you need. Alternatively, you can search for published datasets from papers, but you need to make sure these provide information about celltype in their metadata, or else you need to redefine it yourself using the paper's markers. **Make sure to carefully interrogate all label transfer results, you should view them as data-driven suggestions rather than definitive cell type annotations.**


### 2. RISC reference manual selection

For sample integration we use the ["Robust Integration of Single Cell RNA-seq" (RISC)](https://www.nature.com/articles/s41587-021-00859-x) workflow, as implemented in the [RISC R package](https://github.com/bioinfoDZ/RISC).

RISC requires selection of a reference sample. We provisionally developed an automated method for selection of the reference to ease application of this pipeline. However, the recommended approach by RISC is to inspect the "InPlot" figure.

Optionally, users can provide the reference sample they prefer rather than use auto-selection in the pipeline, for example if the auto selected sample strongly deviates from the reference you would have selected by inspection of "InPlot".





All other parameters will be explained below.


<br />
<br />




# All pipeline parameters


## Key input/output parameters

These parameters are required to run the pipeline and tell it where the data is and where to save the report

- `datadir` string, path to folder containing Cellranger output folders for each sample
- `outdir` string, path to output folder, will be created if doesn't already exist
- `sample_metadata` string, path to a .csv file containing at least two columns: "Sample", matching exactly the sample names in `datadir`, and "Condition", giving the experiment status of that sample, such as WT or KO, Case vs Control, etc. Optionally, can provide a third column "Code" giving a nickname for each sample; this is set to "Sample_Condition" for each sample if not.
- `comps` string, path to a .csv file containing two columns called "c1" and "c2". Each row will be used to compare conditions from the "sample_metadata" csv; multiple comparisons are supported.
- `input_seurat_obj` T/F. If true, will read in Seurat objects from `datadir` with names matching the sample column of `sample_metadata`. Ie, if datadir contains objects called "Sample1.rds", "Sample2.rds", and psuedobulk_metadata has "Sample1" in the Sample column, only Sample1.rds will be read in. Useful for data with some preprocessing or hashed data input.


## Key analysis parameters

These parameters allow for control of sophisticated analysis methods.

- `use_labeltransfer` T/F, whether to use labeltransfer for cell type prediction, default = F
- `refdatapath` string, path to a Seurat object .rds file, pre-processed with `Seurat::SCTransform()`, with a column called "Celltype" in its meta.data. Ignored if `use_labeltransfer` = F.
- `m_reference` string, path to .rds file containing output of `Seurat::FindAllMarkers` run on the reference object specified above. Ignored if `use_labeltransfer` = F

- `risc_reference` - string, name of sample to use as RISC reference sample, if not provided will automate the choice

- `Pseudobulk_mode` - T/F. Sets the cross-conditional analysis mode. TRUE uses pseudobulk EdgeR for DE testing and propeller for compositional analysis. FALSE uses single-cell wilcox test within Seurat for DE testing and 2-prop Z test within the `prop.test()` function for compositional analysis


- `species` string, for example 'Homo sapiens' or 'Mus musculus', default = 'Homo sapiens'; this is for pathway analysis, see `msigdbr::msigdbr_species()`



## QC filtering parameters

We apply an automated filtering approach to remove low-quality / dead cells.
By default, this will perform some baseline filtering removing cells with lower than 500 UMIs, 200 unique genes over 25% mito content, or over 25% hemoglobin content. Additionally, there is some automated filtering for the percent of mito, number of UMIs, and the "cell complexity", or the number of unique genes expected given the number of UMIs.
Finally, we use [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) to remove doublets.


**If using single-nuclei RNA-seq, turn off `autofilter_mito`**. Usually the mito disribution will be so low that it may cause the filtering to select an extremely low cutoff and throw out nuclei unecessarily.

- `min_num_UMI` - numeric, default is 500, if no filter is desired set to -Inf
- `min_num_Feature` - numeric, default is 200, if no filter is desired set to -Inf
- `max_perc_mito` - numeric, default is 25, if no filter is desired set to Inf
- `max_perc_hemoglobin` - numeric, default is 25, if no filter is desired set to Inf
- `autofilter_complexity` - T/F, default T, whether to filter cells with lower than expected number of genes given number of UMIs
- `autofilter_mito` - T/F, default T, whether to filter cells with higher than normal mito content
- `autofilter_nUMI` - T/F, default T, whether to filter cells with lower than normal UMI content
- `autofilter_medianabsolutedev_threshold` - numeric, default is 3, threshold for median abs deviation thresholding, ie cutoffs set to ⁠median +/- mad * threshold⁠
- `autofilter_loess_negative_residual_threshold` - numeric, cutoff for loess residuals applied in complexity filtering, default is -5, if you set it high (ie any higher than -2) you will probably remove many good cells.
- `doubletFinder` - T/F, default is T, whether to filter doublets with DoubletFinder





## Tunable analysis parameters

These parameters affect 


- `pcs_indi` integer, default = 30; number of PCs to use in individual sample processing / clustering
- `res_indi` numeric, default = 0.5; Louvain resolution for individual sample clustering
- `pcs_int`  integer, default = 30; number of PCs to use in integrated data processing / clustering
- `res_int` numeric, default = 0.5 ; Louvain resolution for louvain clustring of RISC integrated dataset; see `SDAP::scCluster_louvain_res()`
- `RISC_louvain_neighbors` integer, default = 10; number of nearest neighbors to consider during clustering; see `RISC::scCluster()` or `SDAP::scCluster_louvain_res()` where implementation of this is unchanged


- `crossconditionDE_padj_thres` numeric, numeric; adjusted p value threshold for significant DE genes in cross condition DE; if `Pseudobulk_mode` is set to T default is 0.1; if `Pseudobulk_mode` is F default is 0.05
- `crossconditionDE_lfc_thres` numeric, absolute value of LFC threshold for significant DE genes in cross condition DE; if `Pseudobulk_mode` is set to F default is 0 (no minimum lFC); if `Pseudobulk_mode` is F default is 0.25



## Computing Resource Allocation


### CPUs: parallelization to increase speed

Running with multiple CPU threads ("workers") can speed up the analysis, *especially if DoubletFinder is used*, but can cause the pipeline to fail due to overuse of memory.

- `workernum` integer, number of CPU threads, default = 1

If on HPC, make sure to also request the appropaite number of CPUs, for example by adding the following two lines to the SBATCH header for 11 cpus:
```
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=11
```

Equivalent Grid Engine / qsub command:
```
#$ -pe smp 11
```

Paralellization is generally implemented across samples, so do not set `workernum` higher than the number of samples.



### Memory Allocation

Memory usage can be high. One run with 11 non-multiplexed samples ran successfully with `workernum` = 11 and SBATCH memory set to 150gb.

Ask for this on SLURM-based HPC schedulers:
```
#SBATCH --mem=150gb
```

Equivalent qsub command, need to divide total mem in GB 150 by number of CPUS.
For 150 over 11 CPUs, ask for 13.64 GB for each CPU.
```
#$ -l h_vmem=13.645g
```
