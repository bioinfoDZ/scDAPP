# scDAPP downstream: Make an app with ShinyCell for exploratory data analysis


It can often be useful to create an app for exploratory data analysis, for example to quickly check gene expression patterns. ShinyCell [(Ouyang et al Bioinformatics 2021)](https://academic.oup.com/bioinformatics/article/37/19/3374/6198103) makes this quite easy. We provide below a short vignette for using scDAPP output files with ShinyCell and uploading to the web and hosting it on [shinyapps.io](https://www.shinyapps.io/).


Steps:
1. Inputs from pipeline
2. optionally, rename, subset, or reorder the metadata
3. create ShinyCell configuration object
4. optionally, control colors of ShinyCell by modifying configuration object
5. Prepare shiny app, will save the app and scripts to a directory of your choice
6. Optional, recommended: modify shiny ui.R script to show random order instead of max 1st
7. Optional, recommended: "deploy" (upload) the app to the internet



## 1. Inputs from pipeline: 

From the pipeline, you will need the file `Seurat-object_integrated.rds`, which can be found below:

```
scDAPP_outdir/
├── individualsample_analysis
├── multisample_integration
│   └── data_objects
│       └── ** Seurat-object_integrated.rds **
├── pipeline_parameters.csv
└── scRNAseq_clustering_integration.html
```




Read in the seurat object:

```
seu <- readRDS("multisample_integration/data_objects/Seurat-object_integrated.rds")
```



## 2. optionally, rename, subset, or reorder the metadata

Shinycell uses the the information in `seu@meta.data`. You should inspect it:

```
head(seu@meta.data,2)
```


Often there can be extraneous metadata variables that are not important to share. Alternatively, you may want to rename or reorder the metadata columns, to control what is plotted first. `seu@meta.data` is essentially a data.frame object that can be saved, or overwritten.

As an example, below we access the metadata as a data.frame called `md`, remove some columns, rename some columns, and reorder some columns. Then, we overwrite `seu@meta.data` with the updated `md` data.frame.


```
library(tidyverse)
library(Seurat)
library(ShinyCell)


#access metadata
md <- seu@meta.data

#remove a column called "sizefactor"
md <- md[,!(colnames(md) %in% 'sizefactor')]

#rename column 3 to 'nUMI'
colnames(md) #check colnames of md
colnames(md)[3] <- 'nUMI' #update colname number 3

#reorder columns to put RISC_Louvain_npc30_res0.5 first
cn <- colanmes(md)
cn <- cn[!(cn %in% c('RISC_Louvain_npc30_res0.5'))]
cn <- c('RISC_Louvain_npc30_res0.5', cn)
md <- md[,cn]


#finally, overwrite seu@meta.data using md. Note the number of rows and row orders should be the same.
table(rownames(md) == rownames(seu@meta.data))
seu@meta.data <- md
```




## 3. create ShinyCell configuration object


Then, create the ShinyCell configuration object. This essentially stores information about the information in `seu@meta.data`

```
#create metadata config file
scConf = createConfig(seu)
```



## 4. optionally, control default metadata variables and colors of ShinyCell by modifying configuration object.

ShinyCell will give random default colors to all metadata object levels (such as clusters or conditions). However, if you want to change the default colors, that can be done also. 

For example, if we have three levels of the Condition column such as Healthy, Mild Disease, or Severe Disease, and we want them to be colored blue, orange, and red, we can do that like so:

```
#set default metadata columns to show. Below will show Clusters first.
scConf <- modDefault(scConf, 'seurat_clusters', 'Condition')


#check configuration object to find out which row corresponds to Condition:
scConf[,1]

#if "Condition" is row 20, we can check the factor level ordering:
scConf[20,'fID']

# Healthy|Covid_Mild|Covid_Critical

#if we want healthy to be blue, mild to be orange, and critical to be red, we can do something like below:
library(gplots)

cols <- c('steelblue', 'orange', 'firebrick')
colscodes <- gplots::col2hex(cols)
scConf[20,c('fCL')] <- paste(colscodes, collapse = '|')
```



## 5. Prepare shiny app, will save the app and scripts to a directory of your choice:


```
makeShinyApp(seu, scConf, gene.mapping = TRUE,
             shiny.title = "Example_Shiny_App_Title",
             shiny.dir = 'apps_ShinyCellExample/',
             gex.assay = 'RISC',
             default.gene1 = 'PTPRC', default.gene2 = 'CD3E',
             default.multigene = c('PTPRC', 'CD3E', 'MS4A1', 'CD68')
             )
```

Then you can use RStudio, open one of the files created in the folder indicated in the `shiny.dir` parameter, such as "ui.R". If you open this file in RStudio, you should see a button with a green arrow that says "Run App". You can click this button to run the app locally (on your own computer). You can also share this others by sharing the app folder, as long as they have RStudio they should be able to open and run it locally.


## 6. Optional, recommended: modify shiny ui.R script to show random order instead of max 1st

By default, gene expression is shown on ShinyCell FeaturePlots using the "Max-1st" option, plotting cells with the highest expression on top. However, this can cause misleading plots, where all cells look positive for all genes. To reduce this, it is recommended to use random order plotting. This can be set to the default by modifying the "ui.R" file.

You can replace all instances of `selected = "Max-1st"` with `selected = "Random"`.
The RStudio word search and replace tool makes this easy (control F or command F to use the tool).


## 7. Optional, recommended: "deploy" (upload) the app to the internet

You can "deploy" (upload) the app to [shinyapps.io](https://www.shinyapps.io/) with a free or premium account like below:

```
library(rsconnect)

#modify below code to conect your account
rsconnect::setAccountInfo(name='InputYourUsername',
                          token='InputYourToken',
                          secret='InputYourPassword')

#modify path to match where your app is saved
rsconnect::deployApp("apps_ShinyCellExample")
```


Other options are described on the ShinyCell website: [ShinyCell deploy guide](https://htmlpreview.github.io/?https://github.com/SGDDNB/ShinyCell/blob/master/docs/4cloud.html)




For more information, see the [ShinyCell website](https://github.com/SGDDNB/ShinyCell).
