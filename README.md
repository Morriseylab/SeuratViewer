# SeuratViewer
R Shiny website for viewing single cell RNA-Seq data analysed using [Seurat](https://satijalab.org/seurat/) 
Seurat is also hosted on GitHub. You can view the repository at

- https://github.com/satijalab/seurat

# Introduction
SeuratViewer reads in the expression data, sample data, feature annotation, dimensionality reduction/ clustering, and marker gene information as an RData object and enables users to view and interact with their single cell RNAseq data

# Requirements
- R
- RStudio Server
- Shiny Server

If you need help installing the above or getting started, refer to [this](https://deanattali.com/2015/05/09/setup-rstudio-shiny-server-digital-ocean/#install-r)

# Installation
Run the config file to install all required packages

# Input Data format
- Creating your dataset
Analyse your single cell data using the [Seurat](https://satijalab.org/seurat/) package. 

- Adding your dataset
Add your data to the param.csv file and move it to the data directory. Please note that the data directory must be in the same location as your server.R, ui.R and function.R files.

# Example
