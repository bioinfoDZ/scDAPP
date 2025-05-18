# Installation

If you run into any issues, you can check the [common bugs and fixes FAQ](https://github.com/bioinfoDZ/scDAPP/blob/main/Documentation/CommonBugs.md). If the error is not reported there, please save the error message and open a Github Issue in this repository.



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
Navigate to yaml file: https://github.com/bioinfoDZ/scDAPP/blob/main/conda_env/2025scdapp.yml

Then download or copy and paste it into a file on the server / HPC. You can click the file and click on the "overlapping boxes" to copy raw file contents then paste it as a file in your location.


Create the environment using the .yml file:
```
conda env create -f 2025scdapp.yml 
```


It may ask your permission, just say yes.

It will take a while to install after you give the permissions.

If it still does not work, please write down the dependencies that cause trouble and the error messages, and open an issue in this repo.


<br />

### 3. Using R, install the package and dependencies in the conda environment:



Activate the conda environment, then within load R:
```
#load conda
conda activate 2025scdapp

#activate R
R

#to quit R, use: 
# q('no') #remove hashtag in front to use this
```


####  In R, install scDAPP:


You can then install this package and its dependencies from Github with:

```
# install.packages("devtools")
devtools::install_github("bioinfoDZ/scDAPP")
```

It will ask your permission, just say yes.

It will take a while to install all the dependencies and this package.

Some packages (such as "XML") may already be present due to Conda, but cannot be updated in R. It is okay, you can ignore it.

It may not work with some packages causing errors. R packages can be installed manually from within R or from conda. If you used the conda virtual environment steps, many R packages are also available via conda (just google package name + "conda") which may be a lot easier than manual installation.

If you run into any issues, you can check the [common bugs and fixes FAQ](https://github.com/bioinfoDZ/scDAPP/blob/main/Documentation/CommonBugs.md). If the error is not reported there, please save the error message and open a Github Issue in this repository.


####  Test installation:

Finally, if you did not get any errors during the package installation steps, you can test the installation by running in R:

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
16           RISC  1.7.0
17         scDAPP  1.0.0
```

You may need to exit R (via `q('no')`) and reload and rerun this if it does not work immediately following installation.


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
devtools::install_github("bioinfoDZ/scDAPP")
```

<br />
<br />


# Development branch

There is also a development branch of this package where new features / big fixes are first pushed and tested before release. 
In case of any breaking errors from upstream dependencies, if you are in a time crunch, you can try checking out the [dev branch changelog](https://github.com/bioinfoDZ/scDAPP/blob/dev/Documentation/Changelog.md) and installing like below:

```
devtools::install_github("bioinfoDZ/scDAPP@dev")
```

Note the dev branch may not be 100% stable and can be subject to frequent updates.




