# SeuratViewer
R Shiny website for viewing single cell RNA-Seq data analysed using [Seurat.](https://satijalab.org/seurat/) 
Seurat is also hosted on GitHub. You can view the repository at

- https://github.com/satijalab/seurat

## Introduction
SeuratViewer reads in the expression data, sample data, feature annotation, dimensionality reduction/ clustering, and marker gene information as an RData object and enables users to view and interact with their single cell RNAseq data

## Requirements
- R (version > 3.4)
- RStudio Server
- Shiny Server (if you need to host it online)

If you need help installing the above or getting started, refer to [this](https://deanattali.com/2015/05/09/setup-rstudio-shiny-server-digital-ocean/#install-r)

## Installation
For Linux, run the following commands in terminal 
```
sudo apt-get install libcurl4-openssl-dev lib-ssl dev
sudo apt-get install xorg
sudo apt-get install xorg libx11-dev mesa-common-dev libglu1-mesa-dev
sudo apt-get install libxml2-dev
sudo apt-get install libftgl2 freetype2-demos libfreetype6-dev
sudo apt-get install libhdf5-dev
sudo apt-get install r-cran-rcppeigen
```
Run the following commands in R to install all required packages
```
install.packages(c("devtools","shiny","shinydashboard","shinyjs","shinyBS","shinyBS","RColorBrewer","reshape2","ggplot2",
                   "dplyr","tidyr","openssl","httr","plotly","htmlwidgets","DT","shinyRGL","rgl","rglwidget","Seurat","cowplot",
                    "data.table","NMF","tibble","network","igraph","visNetwork"))

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite(c("biomaRt","Biobase"))
```
For linux users, other R dependencies include
- RcppEigen
- lme4
- flexmix

## Input Data format
### Creating your dataset 

Analyse your single cell data using the [Seurat](https://satijalab.org/seurat/) package. Run the following code to run FindMarkers function on all clusters. Please note that the object should always be saved as **scrna**.
```
scrna@misc=NA
scrna@misc <-  vector(mode="list", length=length(levels(scrna@ident)))
names(scrna@misc)=levels(scrna@ident)
for(c in levels(scrna@ident)){
  scrna@misc[[c]] <- FindMarkers(scrna,ident.1 = c) %>% tibble::rownames_to_column('gene_name')
  rownames(scrna@misc[[c]])=scrna@misc[[c]]$gene_name
} 
```
Save Seurat object  
```
save(scrna,file="scrna_data1.RData")
```

### Adding your dataset

Add your data to the param.csv file and move it to the data directory. You can find an example dataset [here.](http://165.123.69.6/SeuratViewer/scrna_data1.RData) Please note that the data directory must be in the same location as your server.R, ui.R and function.R files (rename the Example data folder into data). The param.csv file should also be saved in the data directory as the RData files.

### NOTE
Please note that this script requires a username and a password. Before running it, either comment out the Authentication section in server.R or add the username and password in authentication.csv file in the data folder. The username has to be entered in the param.csv file as well so that the user can view only specific datasets.
