# SeuratViewer
R Shiny website for viewing single cell RNA-Seq data analysed using [Seurat.](https://satijalab.org/seurat/) (version 3)
Seurat is also hosted on GitHub. You can view the repository at

- https://github.com/satijalab/seurat

## Introduction
SeuratViewer reads in the expression data, sample data, feature annotation, dimensionality reduction/ clustering, and marker gene information as an RData object and enables users to view and interact with their single cell RNAseq data

## Requirements
- R (version > 3.5)
- RStudio Server
- Shiny Server (if you need to host it online)

If you need help installing the above or getting started, refer to [this](https://deanattali.com/2015/05/09/setup-rstudio-shiny-server-digital-ocean/#install-r)

Please note that certain sections might use functions that require the installation of the following R packages. The installation instructions have been provided below. The packages are 
 - scExtras - provides additional functions for single cell data processing like running dimension reduction methods like tsne, umap and diffusion maps and integrating seurat with [monocle](http://cole-trapnell-lab.github.io/monocle-release/) and [slingshot](https://bioconductor.org/packages/release/bioc/vignettes/slingshot/inst/doc/vignette.html)
 - ligrec - function to compute ligand receptor pairs
 
## Installation
For Linux, run the following commands in terminal 
```
sudo apt-get install libcurl4-openssl-dev libssl-dev
sudo apt-get install xorg libx11-dev mesa-common-dev libglu1-mesa-dev
sudo apt-get install libxml2-dev
sudo apt-get install libftgl2 freetype2-demos libfreetype6-dev
sudo apt-get install libhdf5-dev
sudo apt-get install r-cran-rcppeigen
sudo apt install libgeos-dev
```
Run the following commands in R to install all required packages
```
install.packages(c("devtools","shiny","shinydashboard","shinyjs","shinyBS","shinyBS","RColorBrewer","reshape2","ggplot2",
                   "dplyr","tidyr","openssl","httr","plotly","htmlwidgets","DT","shinyRGL","rgl","rglwidget","Seurat","cowplot",
                    "data.table","NMF","tibble","network","igraph","visNetwork"))

#Install packages from bioconductor
install.packages("BiocManager")
BiocManager::install(c("biomaRt","Biobase","slingshot","ComplexHeatmap","destiny"))


##This package contains helper functions 
require(devtools)
install_github("Morriseylab/ligrec")
install_github("chris-mcginnis-ucsf/DoubletFinder")
install_github("Morriseylab/scExtras")

```
For linux users, other R dependencies include
- RcppEigen
- lme4
- flexmix

An alternative is to use the Dockerfile using [Shinyproxy.](https://github.com/openanalytics/shinyproxy)Docker container can be found in [dockerhub](https://hub.docker.com/repository/docker/apoorvababu/seurat_viewer)

## Creating Input data

### Set parameters
```
outdir <-'~/Seurat' 
projectname<-'project' # specify project name,this will also be the Rdata file name
input10x <- c('LAM_rep1/filtered_feature_bc_matrix/','LAM_rep2/filtered_feature_bc_matrix') # dir(s) of the 10x output files, genes.tsv,barcodes.tsv
org<-'human' 

mouseorthologfile <- 'Example data/mouse_human.csv'
npcs<-50 #How many inital PC dimensions to compute. 
k=30 #This for nearest neighbors, 30 is default
```
### Preprocess the data and create a seurat object
```
dir.create(outdir,recursive = T)
scrna <- processExper(dir=outdir,org=org,name=projectname,files=input10x ,ccscale = T,filter=T)
scrna <- ClusterDR(scrna,npcs=npcs,maxdim='auto',k=k)

scrna= ligrec(object=scrna,org=org)
saveRDS(scrna,file=paste0(projectname,'.RDS'))
```

Or you can analyse your single cell data using the [Seurat](https://satijalab.org/seurat/) package. Run the following functions to find the markers genes in all clusters and find the ligand receptor pairs. Please note that the object should always be saved as **scrna**. 
```
org = "mouse" #use mouse or human based on your dataset
scrna@misc[["findallmarkers"]] <- FindAllMarkers(object = scrna, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scrna= scExtras::ligrec(object=scrna,org=org)
save(scrna,file="scrna_v3_data.RData")
```

Note : You have to specify filetype as RDS or RData in the param file

### Adding your dataset

Add your data to the param.csv file and move it to the data directory. You can find an example dataset [here.](http://165.123.69.6/SeuratViewer/scrna_v3_data.RData) Please note that the data directory must be in the same location as your server.R, ui.R and function.R files (rename the Example data folder into data). The param.csv file should also be saved in the data directory as the RData files.

### NOTE
Please note that this script requires a username and a password. Before running it, either comment out the Authentication section in server.R or add the username and password in authentication.csv file in the data folder. The username has to be entered in the param.csv file as well so that the user can view only specific datasets.
