# SeuratViewer
R Shiny website for viewing single cell RNA-Seq data analysed using [Seurat](https://satijalab.org/seurat/) 
Seurat is also hosted on GitHub. You can view the repository at

- https://github.com/satijalab/seurat

## Introduction
SeuratViewer reads in the expression data, sample data, feature annotation, dimensionality reduction/ clustering, and marker gene information as an RData object and enables users to view and interact with their single cell RNAseq data

## Requirements
- R
- RStudio Server
- Shiny Server

If you need help installing the above or getting started, refer to [this](https://deanattali.com/2015/05/09/setup-rstudio-shiny-server-digital-ocean/#install-r)

## Installation
Run the following command to install all required packages
```
install.packages(c("shiny","shinydashboard","shinyjs","shinyBS","shinyBS","RColorBrewer","reshape2","ggplot2",
                   "dplyr","tidyr","plotly","htmlwidgets","DT","shinyRGL","rgl","rglwidget","Seurat","cowplot",
                    "data.table","NMF","tibble","network","igraph"))
                    
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite(c("biomaRt","Biobase"))
```

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
save(scrna,file="projectname.RData")
```

### Adding your dataset

Add your data to the param.csv file and move it to the data directory. Please note that the data directory must be in the same location as your server.R, ui.R and function.R files. The param.csv file should also be saved in the data directory as the RData files.

Please note that this script requires a username and a password. Before running it, either comment out the Authentication section or add the username password to lines 34-35 of server.R. Also, the username has to be entered in the param.csv file as well.
