library(shinydashboard)
library(shiny)
library(shinyBS)
library(plotly)
library(d3heatmap)
library(shinyjs)
library(rglwidget)
library(reshape2)
options(shiny.sanitize.errors = FALSE)
ui <- dashboardPage(
  dashboardHeader(title = "sEuRaT",titleWidth = 350),
  dashboardSidebar(width = 350,
                   div(style="overflow-y: scroll"),
                   # tags$style(type="text/css",
                   #            ".shiny-output-error { visibility: hidden; }",
                   #            ".shiny-output-error:before { visibility: hidden; }"
                   # ),
                   tags$head(tags$style(HTML(".sidebar { height: 250vh; overflow-y: auto; }" ))),
                   sidebarMenu(
                     menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
                     uiOutput("projects"),
                     #actionButton('load',"Load Project", icon = NULL),
                     menuItem('Project Summary', tabName = 'summ', icon = icon('hand-o-right')),
                     menuItem('Variable Genes', tabName = 'vargenes', icon = icon('hand-o-right'),badgeLabel = "new", badgeColor = "green"),
                     menuItem('Principle component Analysis', tabName = 'pca', icon = icon('hand-o-right'),badgeLabel = "new", badgeColor = "green"),
                     menuItem('tSNE Plots', tabName = 'tplot', icon = icon('hand-o-right'),
                              menuSubItem("Compare tSNE Plots", tabName = "tsneplot"),
                              menuSubItem("Interactive tSNE/uMap Plot", tabName = "intertsne")
                              ),
                     menuItem('Biplot', tabName = 'biplot', icon = icon('hand-o-right')),
                     menuItem('Differential Expression', tabName = 'deg', icon = icon('hand-o-right')),
                     menuItem('Gene Expression Plots', tabName = 'geplot', icon = icon('hand-o-right'),
                              fluidRow(
                                column(1,h4("")),
                                column(6,menuSubItem("Gene Expression Plots", tabName = "geplots")),
                                column(2,bsButton("q1", label = "", icon = icon("question"), style = "info", size = "extra-small"))),
                              bsTooltip(id = "q1", title = "Enter Gene name and view tSNE, violin and Ridge plots",placement = "right",trigger = "hover", options = NULL),
                              fluidRow(
                                column(1,h4("")),
                                column(6,menuSubItem('Cluster-wise Gene Expression', tabName = 'clusplot')),
                                column(2,bsButton("q2", label = "", icon = icon("question"), style = "info", size = "extra-small"))),
                              bsTooltip(id = "q2", title = "View tSNE/umap of genes expressed in a Cellgroup",placement = "right",trigger = "hover", options = NULL),
                              fluidRow(
                                column(1,h4("")),
                              column(6,menuSubItem('Gene-Cellgroup Dotplot', tabName = 'dotplot')),
                              column(2,bsButton("q3", label = "", icon = icon("question"), style = "info", size = "extra-small"))),
                              bsTooltip(id = "q3", title = "Upload genelist to view the gene expression across cellgroups as a dotplot",placement = "right",trigger = "hover", options = NULL)
                              ),
                     
                     menuItem('Heatmap', tabName = 'heatmap', icon = icon('hand-o-right')),
                     menuItem('Ligand Receptor Pairs', tabName = 'ligrec', icon = icon('hand-o-right')),
                     menuItem('Ligand Receptor Network', tabName = 'network', icon = icon('hand-o-right'))
                   )#end of sidebar menu
  ),#end dashboardSidebar
  
  dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
    ),
    useShinyjs(),
    tabItems(
      tabItem(tabName = "dashboard",
              box(
                width = 12, status = "primary",solidHeader = TRUE,
                title = "scRNA Data Sets",
                tableOutput("datasetTable")
              )
      ),
      ######################################################################################################################################
      
      tabItem(tabName = "summ",
              box(
                width = 12, status = "primary",solidHeader = TRUE,
                title = "Project Summary",
                htmlOutput("prjsumm")
              )
      ),
      
    ######################################################################################################################################
    tabItem(tabName = "vargenes",
            box(title = "Variable Genes",solidHeader = TRUE,width=12,status='primary',
                DT::dataTableOutput('vargenes')
            )#End box
    ),#end of vargenes
    ######################################################################################################################################
    tabItem(tabName = "pca",
            box(title = "Gene Plots",solidHeader = TRUE,width=8,status='primary',
                plotOutput("vizplot", height = 800)
            ),#End box
            box(title = "Controls",solidHeader = TRUE,width=4,status='primary',
                uiOutput("ndim"),
                sliderInput("ngenes", "Number of genes:",min = 10, max = 50, value = 10,step=5))
    ),#end of pca
    ######################################################################################################################################
    tabItem(tabName = "tsneplot",
            box(title = "Compare tSNE plots",solidHeader = TRUE,width=12,status='primary',
                fluidRow(
                  column(6,uiOutput("umapa")), #Dimensionality reduction method of left plot
                  column(6,uiOutput("umapb")) #Dimensionality reduction method of left plot
                ),
                fluidRow(
                  column(6,selectInput("categorya2", "Select one",c('Categories' = "var",'Cluster' = "clust", 'Gene Expression' = "geneexp"),selected = "clust")),
                  column(6,selectInput("categoryb2", "Select one",c('Categories' = "var",'Cluster' = "clust", 'Gene Expression' = "geneexp"),selected = "clust"))
                ),
                sliderInput("pointa2", "Point Size:",min = 0, max = 5, value = 1,step=.25),
                fluidRow(
                  column(6,conditionalPanel(
                    condition = "input.categorya2 == 'var'",
                    uiOutput("tsnea2") # Generate list of variables for left plot
                  ),
                  conditionalPanel(
                    condition = "input.categorya2 == 'geneexp'",textInput("gene1a", label = "Gene Name",value = "Axin2")
                  )
                  ),
                  column(6,conditionalPanel(
                    condition = "input.categoryb2 == 'var'",
                    uiOutput("tsneb2")  # Generate list of variables for right plot
                  ),
                  conditionalPanel(
                    condition = "input.categoryb2 == 'geneexp'",textInput("gene2a", label = "Gene Name",value = "Axin2")
                  )
                  )),
                fluidRow(
                  column(6,checkboxInput("subsa", label = "Check to subselect cells", value = FALSE)),
                  column(6,checkboxInput("subsb", label = "Check to subselect cells", value = FALSE))
                ),
                
                fluidRow(
                  column(6,conditionalPanel(
                    condition = "input.subsa ==true",uiOutput("subsaui") #generate ident list for left plot
                  )),
                  column(6,conditionalPanel(
                    condition = "input.subsb ==true",uiOutput("subsbui") #generate ident list for right plot
                  )
                  )),
                plotOutput("comptsne2", height = 600),
                downloadButton('downloadtsneplot', 'Download tSNE plot')
            )
    ),#end of degtab
    
    ###################################################################################################################################### 
    tabItem(tabName = "intertsne",
            fluidRow(
              box(
                title = "Controls",solidHeader = TRUE,width=12,status='primary',
                fluidRow(
                  column(6,uiOutput("umapint")),
                  column(6,selectInput("intercat", "Select one",c('Categories' = "var",'Gene Expression' = "geneexp"),selected = "geneexp"))),
                fluidRow(
                  column(6,uiOutput("setcategory")),
                  column(6,conditionalPanel(
                    condition = "input.intercat == 'var'",
                    uiOutput("intervar")
                  ),
                  conditionalPanel(
                    condition = "input.intercat == 'geneexp'",textInput("geneinter", label = "Gene Name",value = "Axin2")
                  )
                  )),
                sliderInput("umap_pointsize", "Point Size:",min = 0, max = 5, value = 1,step=.25)
              ),
              box(plotlyOutput("intertsne", height = 600),width=6, status='primary',title = "Interactive tSNE/uMap Plot",solidHeader = TRUE),
              box(plotlyOutput("intergene", height = 600),width=6, status='primary',title = "Interactive Gene Plot",solidHeader = TRUE)
              
            )
    ),#endbigeneplotTab
    ###################################################################################################################################### 
    
    tabItem(tabName = "biplot",
            fluidRow(
              box(plotOutput("bigeneplot", height = 600),width=8, status='primary',title = "Bigene Plot",solidHeader = TRUE),

              box(
                title = "Controls",solidHeader = TRUE,width=4,status='primary',
                textInput("bigene_genea", label = "Gene A",value = "Sox2"),
                #sliderInput("bigene_rangea", "Expression Limits Gene A(log2(UMI))",min = 0, max = 10, value = 0.5,step=.25),
                uiOutput("bigene_rangea"),
                textInput("bigene_geneb", label = "Gene B",value = "Sox9"),
                #sliderInput("bigene_rangeb", "Expression Limits Gene B(log2(UMI))",min = 0, max = 10, value = 1.5,step=.25),
                uiOutput("bigene_rangeb"),
                sliderInput("bigene_pointsize", "Point Size:",min = 0, max = 5, value = 1,step=.25),
                downloadButton('downloadbiplot', 'Download Bigene plot')
              )
            )
    ),#endbigeneplotTab
    ######################################################################################################################################
    tabItem(tabName = "deg",
            box(title = "Differentially expressed Markers",solidHeader = TRUE,width=9,status='primary',
                plotOutput("comptsne", height = 1000)
            ),

            fluidRow(
              box(title = "Controls",solidHeader = TRUE,width=3,status='primary',
                  uiOutput("umapdeg"),
                uiOutput("tsnea"),
               
                sliderInput("pointa", "Point Size:",min = 0, max = 5, value = 1,step=.25),
                checkboxInput("checkviolin", label = "Check to remove points from violin plot", value = TRUE),
                hr(),
                uiOutput("identdef"),
                checkboxInput("setident", label = "Check to choose a different category to compare", value = FALSE),
                conditionalPanel(
                  condition = "input.setident ==true",uiOutput("setidentlist"),
                  uiOutput("identa"),
                  uiOutput("identb"),
                  actionButton("goButton", "Go!"),
                sliderInput("lfc", "Log FC threshold:",min = 0.25, max = 6, value = 0.25,step=.25),
                selectInput("test", "Select test to use",c('Wilcox' = "wilcox",'T-test' = "t", 'Poisson' = "poisson",'Negative Binomial'="negbinom"),selected = "wilcox"),
                sliderInput("minpct", "Minimum Percent of cells:",min = 0.1, max = 10, value = 0.25)),
                downloadButton('downloaddeg', 'Download table'),
                downloadButton('downloadplot', 'Download Plot')
              ),
              box(DT::dataTableOutput('markergenes'),width=12, status='primary',solidHeader = TRUE,title="Marker genes")
            )#End FluidRow
    ),#end of degtab
    ######################################################################################################################################
    tabItem(tabName = "geplots",
            box(title = "Gene Expression Plots",solidHeader = TRUE,width=9,status='primary',
                plotOutput("geplots", height = 1000)
            ),
            
            fluidRow(
              box(title = "Controls",solidHeader = TRUE,width=3,status='primary',
                  uiOutput("umapge"),
                  textInput("geneid", label = "Enter Gene",value = ""),
                  sliderInput("genenid_pointsize", "Point Size:",min = 0, max = 5, value = 1,step=.25),
                  downloadButton('downloadplotge', 'Download Plot'))
            )#End FluidRow
    ),#end of geplot
    ######################################################################################################################################
    tabItem(tabName = "clusplot",
            box(title = "Gene Expression Plots",solidHeader = TRUE,width=9,status='primary',
                plotOutput("clustplots", height = 700)
            ),
            fluidRow(
              box(title = "Controls",solidHeader = TRUE,width=3,status='primary',
                  uiOutput("umapclust"),
                  uiOutput("setvar"),
                  uiOutput("selectcluster"),
                  sliderInput("pctslider", "Percent Expressed:",min =0, max = 1, value =c(0.1,1)),
                  uiOutput("avgexpslider"),
                  sliderInput("pointclust", "Point Size:",min = 0, max = 5, value = 1,step=.25),
                  downloadButton('downloadclustplot', 'Download Plot'),
                  downloadButton('downloadclustertab', 'Download Table'))
            ),
            box(title = "Table",solidHeader = TRUE,width=9,status='primary',
                DT::dataTableOutput('clustable')
            )#End box
    ),#end of clusplot
    ######################################################################################################################################
    tabItem(tabName = "dotplot",
            box(title = "Controls",solidHeader = TRUE,width=12,status='primary',
                fluidRow(
                  column(6,fileInput('genelistfile', 'Upload Text File',accept=c('text/csv','text/comma-separated-values,text/plain','.txt'))),
                  column(6,uiOutput("setdotvar"))
                )
            ),
            box(title = "Dot Plot",solidHeader = TRUE,width=12,status='primary',
                plotOutput("dotplot", height = 500)
            )
    ),#end of dotplot
        ######################################################################################################################################
    tabItem(tabName = "heatmap",
            box(
              title = "Controls",solidHeader = TRUE,width=12,status='primary',
              uiOutput("heatmapgenes"),
              uiOutput("hmpgrp"),
              selectInput("hmpcol", "Select one",c('PurpleYellow' = "PuYl",'BlueGreen' = "BuGn", 'RedYellow' = "RdYl", 'RedBlue'="RdBu"),selected = "geneexp"),
              downloadButton('downloadheatmap', 'Download Heatmap')
             ),
            box(plotOutput("heatmap", height = 900),width=12, status='primary',solidHeader = TRUE,title="Single cell heatmap of gene expression")
    ),#end of tab
    ######################################################################################################################################
    tabItem(tabName = "ligrec",
            box(width=9, status='primary',title = "Bigene Plot",solidHeader = TRUE,
              plotOutput("bigeneplot2", height = 800)
              ),
            box(width = 3, status = "primary",solidHeader = TRUE,title = "Controls",
                uiOutput("pairby"),
                radioButtons("clust","Select Cluster", c("All clusters"="all","Select Cluster"="clust"),selected = "all"),
                radioButtons("gene","Select Genes", c("All genes"="allgene","Enter Genelist"="genelist"),selected = "allgene"),
                
                conditionalPanel(
                  condition = "input.clust == 'all' && input.gene == 'genelist'" ,
                  fileInput('genelist1', 'Upload Receptor Genelist',accept=c('text/csv','text/comma-separated-values,text/plain','.txt')),
                  fileInput('genelist2', 'Upload Ligand Genelist',accept=c('text/csv','text/comma-separated-values,text/plain','.txt'))
                ),
                conditionalPanel(
                  condition = "input.clust == 'clust' && input.gene == 'allgene'" ,
                  uiOutput("clust1"),
                  uiOutput("clust2")
                ),
                conditionalPanel(
                  condition = "input.clust == 'clust' && input.gene == 'genelist'" ,
                  fileInput('genelist1.1', 'Upload Receptor Genelist',accept=c('text/csv','text/comma-separated-values,text/plain','.txt')),
                  fileInput('genelist2.1', 'Upload Ligand Genelist',accept=c('text/csv','text/comma-separated-values,text/plain','.txt')),
                  uiOutput("clust1.1"),
                  uiOutput("clust2.1")
                ),
                fluidRow(
                  column(6,checkboxInput("checksource", label = "Check to select by source", value = FALSE)),
                  column(6,checkboxInput("checkevi", label = "Check to select by evidence", value = FALSE)),
                  conditionalPanel(
                    condition = "input.checksource ==true",
                    column(6,uiOutput('source'))
                  ),
                  conditionalPanel(
                    condition = "input.checkevi ==true",
                    column(6,uiOutput('evidence'))
                  )
                ),
                uiOutput("bigene_rangea2"),
                uiOutput("bigene_rangeb2"),
                sliderInput("bigene_pointsize2", "Point Size:",min = 0, max = 5, value = 1,step=.25)
            ),
            box(
              width = 12, status = "primary",solidHeader = TRUE,
              title = "Ligand Receptor pairs",
              DT::dataTableOutput('pairs_res')
            )#end of box
            
           
    ),#end of tabitem
    ######################################################################################################################################
    tabItem(tabName = "network",
            box(title = "Network",solidHeader = TRUE,width=12,status='primary',
                plotOutput("lrnetwork", height = 800)
            ),#End box
            box(title = "Ligand-Receptor pairs",solidHeader = TRUE,width=12,status='primary',
                DT::dataTableOutput('pairs_res2')
            ))#end of network
    ######################################################################################################################################
    
    )#end of tabitems
  )#end of dashboard body
)#end of dashboard page

