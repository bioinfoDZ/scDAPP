# Installation

If you run into any issues, you can check the [common bugs and fixes FAQ](CommonBugs.md). If the error is not reported there, please save the error message and open a Github Issue in this repository.

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

On Einstein HPC, conda can be sourced like this. It will change for you depending on which folder you install Conda in:

```
source /gs/gsfs0/home/aferrena/packages/miniconda3/miniconda3/etc/profile.d/conda.sh
```


I do not recommend using the pre-installed conda on HPC, but rather starting from fresh. To use pre-installed conda, simply do not source the newly installed conda.
 




<br />



### 2. Create a conda environment using yaml file:
Navigate to yaml file: [conda_env/2026_scdappv2.yml](../conda_env/2026_scdappv2.yml)

Then download or copy and paste it into a file on the server / HPC. You can click the file and click on the "overlapping boxes" to copy raw file contents then paste it as a file in your location.


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

You should see messages about packages being activated followed by a data.frame showing the key dependency packages and their versions such as below. If any errors occur, it means something has gone wrong with one or more dependency packages.

You may need to exit R (via `q('no')`) and reload and rerun this if it does not work immediately following installation.

<br />



# Mac tips

If using Mac, the conda virtual environment steps may not work out of the box. As long as you can get R >= v4.0, Seurat >= v5.0, and RISC >= v1.7 working, you should not have problems running the pipeline.

If on Mac then very likely what you will need to do is install Apple Xcode from the app store, then open a terminal and run the following:
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

The [`dev`](https://github.com/bioinfoDZ/scDAPP/tree/dev) branch carries new features and bug fixes. Install from GitHub:

```
devtools::install_github("bioinfoDZ/scDAPP@dev", dependencies = TRUE)
```

The dev branch may not be 100% stable and can be subject to frequent updates. If you need a known-stable v1.3.x build, install a tagged release instead (see [Releases](https://github.com/bioinfoDZ/scDAPP/releases)).




