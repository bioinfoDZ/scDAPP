# Installation


If you have any issues, please email BOTH alexander.ferrena@einsteinmed.edu and alexanderferrena@gmail.com

Or open an issue in this github repo.

<br />

### This guide

The general steps for installing described by this guide are:
 1. install R and non-R dependencies within a conda virtual environment (steps 1-3, for linux systems)
 2. install R dependencies and the pipeline itself in R (step 4)



If using a linux server / HPC, please use the steps below.
If installing on Mac, please see the bottom.


<br />
<br />


### 1. Install Conda: 

https://docs.conda.io/en/latest/miniconda.html

You will likely need to restart the terminal for the install to finish.


#### Note for Einstein HPC users:

On Einstein HPC, I source my local installed Conda this like this. It will change for you depending on which folder you install Conda in:
```source /gs/gsfs0/home/aferrena/packages/miniconda3/miniconda3/etc/profile.d/conda.sh```


I do not recommend using the pre-installed conda on HPC, but rather starting from fresh. To use pre-installed conda, simply do not source the newly installed conda.
 

#### General SLURM HPC usage:

You will need to source Conda each time you want to use it, for example in submission scripts. Adding to bash_rc file will load it in your login node, but not in submitted jobs.



<br />



### 2. Highly recommended: install mamba: 

This is not necessary but *highly* recommended to speed up the installation.

```
conda install -c conda-forge mamba
```

With mamba installed, you can then use "mamba" in place of "conda" for any command for faster execution.
Without mamba, just use "conda" instead of "mamba" below.

<br />

### 3. Create conda environment using yaml file:
Navigate to yaml file: https://github.com/FerrenaAlexander/scDAPP/blob/main/conda_env/r_env.yml

Then download or copy and paste it into a file on the server / HPC. You can click file and click on the "overlapping boxes" to copy raw file contents then paste it as a file in your location.


Create the environment using the .yml file:
```
mamba env create -f r_env.yml 
```


It may ask your permission, just say yes.

It will take a while to install after you give the permissions.

If it still does not work, please write down the packages that give issues and email Alex and/or open an issue in this repo.


<br />

### 4. Using R, install the package and dependencies in the conda environment:



Activate the conda environment, then within load R:
```
# load conda
conda activate r_env

#switch to R
R

#to quit R, use: 
# q('no') #remove hashtag in front to use this
```




Then, in R, you can then install the package with:

```
devtools::install_github("FerrenaAlexander/scDAPP")
```

It will ask your permissions, just say yes.

It will take a while to install all the dependencies and this package.

Some packages (such as "XML") may be present due to Conda, but cannot be updated in R. It is okay, you can ignore it.

It may not work with some packages causing errors. R packages can be installed manually with R. If you used the conda virtual environment steps, many R packages are also available via conda/mamba (just google package name + conda).

If you run into any issues, please open an issue in the github issues page for this repo or email Alex.



Finally, if you did not get any errors during the package installation steps, you can test the install by running:

```
scDAPP::r_package_test()
```

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
16           RISC  1.6.0
17         scDAPP  1.0.0
```

You may need to exit R (via `q('no')`) and reload and rerun this if it does not work.


<br />
<br />



# If using Mac, start from here

If using Mac, the conda virtual environment steps will likely not work. Instead, just skip the conda virtual environment steps and go directly to R. As long as you can get R > v4.0, Seurat > v4.0 or v5.0, and RISC > v1.6 working, you should not have problems running the pipeline.

If on Mac then very likely, what you will need to do is install Apple Xcode from the app store, then open a terminal and run the following:
- `xcode-select --install `
- `sudo xcodebuild -license accept`

You may also need to download the GNU Fortan compiler, accessible from this website: 
- https://mac.r-project.org/tools/

These steps will install key compiler tools that are not easy to install any other way. Then, you should be able to install R > 4.0 for your system [from Cran](https://cran.r-project.org/), and finally within R install [Seurat using its instructions](https://satijalab.org/seurat/articles/install.html) and [RISC from github](https://github.com/bioinfoDZ/RISC).


Then, in R, you can then install scDAPP with:

```
devtools::install_github("FerrenaAlexander/scDAPP")
```







