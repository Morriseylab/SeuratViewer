library(shiny)
library(shinyBS)
library(RColorBrewer)
library(biomaRt)
library(Biobase)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)
library(shinyjs)
library(htmlwidgets)
library(DT)
#options(rgl.useNULL=TRUE)
#library(shinyRGL)
#library(rgl)
#library(rglwidget)
library(Seurat)
#library(Seurat, lib.loc="~/R_lib")
library(cowplot)
#library(cowplot,lib.loc = "~/R_lib")
library(data.table)
library(NMF)
library(tibble)
library(network)
library(igraph)
#library(igraph,lib.loc = "~/R_lib")
library(shinyBS)
#library(slingshot)
source("functions.R")

#Specify color palette for the tSNE and UMAP plsots
cpallette=c("#64B2CE", "#DA5724", "#74D944", "#CE50CA", "#C0717C", "#CBD588", "#5F7FC7",
            "#8B4484", "#D3D93E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
            "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D",
            "#8A7C64", "#599861","#000099","#FFCC66","#99CC33","#CC99CC","#666666", "#695F74")

#Specify user-ids and passwords
auth=read.csv("data/authentication.csv")
my_username <- auth$user
my_password <- auth$pwd

server <- function(input, output,session) {
  
  ###################################################
  ###################################################
  ################### AUTHENTICATION  ##############
  ###################################################
  ###################################################
   values <- reactiveValues(authenticated = FALSE)

   # Return the UI for a modal dialog with data selection input. If 'failed'
   # is TRUE, then display a message that the previous value was invalid.
   dataModal <- function(failed = FALSE) {
     modalDialog(
       textInput("username", "Username:"),
       passwordInput("password", "Password:"),
       footer = tagList(
         # modalButton("Cancel"),
         actionButton("ok", "OK")
       )
     )
   }

   # Show modal when button is clicked.
   # This `observe` is suspended only whith right user credential

   obs1 <- observe({
     showModal(dataModal())
   })

  # When OK button is pressed, attempt to authenticate. If successful,
  # remove the modal.
  obs2 <- observe({
    req(input$ok)
    isolate({
      Username <- input$username
      Password <- input$password
    })
    Id.username <- which(my_username == Username)
    Id.password <- which(my_password == Password)
    if (length(Id.username) > 0 & length(Id.password) > 0) {
      if (Id.username == Id.password) {
        Logged <<- TRUE
        values$authenticated <- TRUE
        obs1$suspend()
        removeModal()

      } else {
        values$authenticated <- FALSE
      }
    }
  })

  ###################################################
  ###################################################
  ####### Display username in notification bar  ####
  ###################################################
  ###################################################
  output$userloggedin = renderMenu({
    msg= paste("Logged in as ",input$username,sep="")
    prj= paste("Data: ",projectname(),sep="")
    dropdownMenu(type = "notifications", badgeStatus = "info",
                 notificationItem(icon = icon("user"), status = "info",msg),
                 notificationItem(icon = icon("book"), status = "info",prj))
  })
  ###################################################
  ###################################################
  ####### Display project list and load data  #######
  ###################################################
  ###################################################
  #Read the parameter file
  readexcel = reactive({
     user=input$username
     file = read.csv("data/param.csv")
    if(user=="allusers" | user=="admin" ){
      file=file
    }else{
      file=file[file$user==user,]
    }
  })
  
  #Get Project list and populate drop-down
  output$projectlist = renderUI({
    excel=readexcel()
    prj=as.character(excel$projects)
    selectInput("projects","Select a project",as.list(sort(as.character(prj))))
  })
  
  #display project list in Dashboard
  output$datasetTable<- renderTable({
    user=input$username
    file=read.csv('data/param.csv',stringsAsFactors = F)
    colnames(file)=c("Project Name","Project Description","Organism","Username")
    file=file[order(file$`Project Name`),]
    if(user=="allusers" | user=="admin"){
      file=file
    }else{
      file=file[file$Username==user,] %>% dplyr::select(-Username)
      colnames(file)=c("Project Name","Project Description","Organism")
      file=file[order(file$`Project Name`),]
    }
  }, digits = 1)
  
  #scrna <- reactiveValues(scrna = NULL)
  
  # observeEvent(input$load, {
  #   inFile = paste('data/',as.character(input$projects),'.RData',sep = '')
  #   withProgress(session = session, message = 'Loading Data...',detail = 'Please Wait...',{
  #   load(inFile)
  #   })
  #   scrna <<- scrna
  # })  
  
  #Load Rdata
  fileload <- reactive({
    if(input$filetype == 'list'){
    inFile = paste('data/',as.character(input$projects),'.RData',sep = '')
    load(inFile)}
    else{
      file=input$rdatafileupload
      scrna=readRDS(file$datapath)
    }
    return(scrna)
  })
  
  #Create variable for project name
  projectname = reactive({
    if(input$filetype == 'list'){
      project= input$projects
    }else if(input$filetype == 'upload'){
      project= input$rdatafileupload
      project=project$name
      project=strsplit(project,".RDS")
      project=sapply(project,"[",1)
    }
      return(project)
  })
  ###################################################
  ###################################################
  ################## Project Summary  ###############
  ###################################################
  ###################################################
  
  #Get all information from the scrna object (input file) and generate some basic project summary for the summary
  prjsumm <- reactive({
    user=input$username
    #user="allusers"
   prj= read.csv("data/param.csv")
   if(user=="allusers" | user=="admin"){
     prj=prj
   }else{
     prj=prj[prj$user==user,] 
   }
   prj=prj[prj$projects==projectname(),]
   pname=projectname()
   if(input$filetype == 'list'){
   pdesc=prj$desc
   porg=prj$organism}
   else{
     pdesc=""
     porg="prj$organism"
   }
   scrna=fileload()
   numcells.nf=dim(GetAssayData(object = scrna, slot = "counts"))[2]
   tcells=dim(GetAssayData(object=scrna))[2]
   tgenes=dim(GetAssayData(object=scrna))[1]
   if(is.null(scrna[["pca"]])){
   maxdim=dim(scrna@dr$cca.aligned@cell.embeddings)[2]
   }else{maxdim=length(scrna[["pca"]]@stdev)}
   c.cnt=as.data.frame(table(Idents(object=scrna)))
   df=as.data.frame(c(as.character(pname),as.character(pdesc),as.character(porg),numcells.nf,tcells,tgenes,maxdim,"","",c.cnt$Freq))
   rownames(df)=c("Project name","Project Description","Organism","Total number of cells before filtration","Total nummber of cells after filtering","Total number of genes","Dimension","","Cluster-wise number of genes",as.character(c.cnt$Var1))
   colnames(df)<- NULL
   return(df)
  })
  
  output$prjsumm <- renderPrint({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
    prjsumm()
    #return(df)
    })
  })
  
  ######################################################################################################
  ######################################################################################################
  ###################### CALC PARAMETERS TAB ############################################################
  ######################################################################################################
  ######################################################################################################
  #Create subtabs for the QC plots
  output$plotsubtab <- renderUI({
    tabsetPanel(id = "subTabPanel1",
                tabPanel("Create Seurat Object",tableOutput("seurobj")),
                tabPanel("Filter Cells and Scale Data",tableOutput('filtobj')),
                tabPanel("Dimension Reduction",tableOutput('dimred')),
                tabPanel("Find Clusters",DT::dataTableOutput('vargenes'))
    )
    
  })
  
  #Seuarat Object
  seurobj= reactive({
    scrna=fileload()
    ulist=scrna@calc.params$CreateSeuratObject
    ulist=ulist[1:11]
    dd  <-  as.data.frame(matrix(unlist(ulist), nrow=length(unlist(ulist[1]))))
    t=unlist(ulist)
    colnames(dd)=names(t)
    return(dd)
  })
  
  output$seurobj <- renderTable({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      seurobj()
    })
  })
  
  #Filter criteria
  filtobj= reactive({
    scrna=fileload()
    ulist=scrna@calc.params$FilterCells
    t=unlist(ulist)
    vars=scrna@calc.params$ScaleData$vars.to.regress
    vars=paste(vars,sep="",collapse=",")
    vars=as.data.frame(vars)
    rownames(vars)="variables regressed"
    colnames(vars)="t"
    t=as.data.frame(t)
    t=rbind(t,vars)
    return(t)
  })
  
  output$filtobj = renderTable({
    withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
      filtobj()    })
  })
  
  #Dimension reduction
  dimred= reactive({
    scrna=fileload()
    pca=scrna@calc.params$RunPCA$pcs.compute
    tsne=length(scrna@calc.params$RunTSNE$dims.use)
    umap=length(scrna@calc.params$RunUMAP$dims.use)
    umap.nei=scrna@calc.params$RunUMAP$n_neighbors
    df=as.data.frame(c(pca,tsne,umap,umap.nei))
    rownames(df)=c("No.PCs.computed","tSNE.dims.used","UMAP.dims.used","UMAP.n_neighbors")
    colnames(df)="values"
    return(df)
  })
  
  output$dimred = renderTable({
    withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
      dimred()    })
  })
  ######################################################################################################
  ######################################################################################################
  ###################### VARIABLE GENES TAB ############################################################
  ######################################################################################################
  ######################################################################################################
  
  #Load the scrna input RData file and extract the variable genes and its stats
  vargenes= reactive({
    scrna=fileload()
    var=as.data.frame(VariableFeatures(object = scrna))
    #var=as.data.frame(scrna@var.genes)
    colnames(var)="Gene"
    stat=HVFInfo(object = scrna)
    stat$Gene=rownames(stat)
    var=left_join(var,stat,by="Gene")
  })
  
  output$vargenes = DT::renderDataTable({
    withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
      DT::datatable(vargenes(),
                    extensions = c('Buttons','Scroller'),
                    options = list(dom = 'Bfrtip',
                                   searchHighlight = TRUE,
                                   pageLength = 10,
                                   lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                   scrollX = TRUE,
                                   buttons = c('copy', 'print')
                    ),rownames=FALSE,caption= "Variable genes",selection = list(mode = 'single', selected =1),escape = F)
    })
  })
  
  ######################################################################################################
  ######################################################################################################
  ################################# PCA TAB ############################################################
  ######################################################################################################
  ######################################################################################################
  
  #Generate drop down menu to populate the max number of dimensions used in the scRNA analysis
  output$ndim = renderUI({
    scrna=fileload()
    maxdim="NA"
    maxdim=length(scrna[['pca']]@stdev)
    validate(need(is.na(maxdim)==F,"PCA dimensional Reduction has not been computed"))
    var=1:maxdim
    selectInput("ndim","Choose number of dimensions",var,selected = 1)
  })
  
  #Plot the PCA/Viz plot for the number of dimensions and number of genes chosen
  vizplot= reactive({
    scrna=fileload()
    dim=input$ndim
    validate(need(dim,"PCA dimensional Reduction has not been computed"))
    par(mar=c(4,5,3,3))
    g1=VizDimLoadings(object = scrna, dims = dim:dim,nCol=1,nfeatures = input$ngenes)
    return(g1) 
  })
  
  #Render the vizplot
  output$vizplot = renderPlot({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      vizplot()
    })
  })
  
  
  ###################################################
  ###################################################
  ####### Compare Tsne plot with controls  ##########
  ###################################################
  ###################################################
  
  #generate variable list for left plot
  output$tsnea2 = renderUI({
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data) 
    #metadata=metadata %>% select(starts_with("var"))
    var=colnames(metadata)
    selectInput("tsnea2","Select a Variable",var,"pick one")
  })
  
  #generate variable list for right plot
  output$tsneb2 = renderUI({
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data) 
    #metadata=metadata %>% select(starts_with("var"))
    var=colnames(metadata)
    selectInput("tsneb2","Select a Variable",var,"pick one")
  })
  
  #Conditional panel. When subselect cells is chosen, if category is 'var', generate a slider with min and max values for numerical categories
  # and drop down menu with factors of that category if its non-numeric. If category is 'cluster', generate drop-down with cluster names and if
  # category is 'geneexp', display error that says cannot subselect
  output$subsaui = renderUI({
    scrna=fileload()
    clusts=levels(Idents(object=scrna))
    if(input$categorya2=="clust"){
      selectInput("selclust","Select a Cluster",clusts)
    }else if(input$categorya2=="var"){
      metadata=as.data.frame(scrna@meta.data)
      met= sapply(metadata,is.numeric)
      feature=names(met[met==TRUE])
      tsne=names(met[met==FALSE])
      t=paste("scrna@meta.data$",input$tsnea2,sep="")
      if(input$tsnea2 %in% tsne){
        opt1=unique(eval(parse(text=t)))
        selectInput("selclust2","Select one of the options",opt1)
      }else if(input$tsnea2 %in% feature){
        min=min(eval(parse(text=t)))
        max=max(eval(parse(text=t)))
        sliderInput("tsnea2lim", label = h5("Select Range"), min = min,max =max, value =c(min,max)) 
      }
    }else if(input$categorya2=="geneexp"){
      validate(need(input$categorya2!="geneexp","Cannot subselect gene expression values"))
    }
  })
  
  #Conditional panel. When subselect cells is chosen, if category is 'var', generate a slider with min and max values for numerical categories
  # and drop down menu with factors of that category if its non-numeric. If category is 'cluster', generate drop-down with cluster names and if
  # category is 'geneexp', display error that says cannot subselect
  output$subsbui = renderUI({
    scrna=fileload()
    clusts=levels(Idents(object=scrna))
    if(input$categoryb2=="clust"){
      selectInput("selclustb","Select a Cluster",clusts)
    }else if(input$categoryb2=="var"){
      metadata=as.data.frame(scrna@meta.data)
      met= sapply(metadata,is.numeric)
      feature=names(met[met==TRUE])
      tsne=names(met[met==FALSE])
      t=paste("scrna@meta.data$",input$tsneb2,sep="")
      if(input$tsneb2 %in% tsne){
        opt2=unique(eval(parse(text=t)))
        selectInput("selclustb2","Select one of the options",opt2)
      }else if(input$tsneb2 %in% feature){
        min=min(eval(parse(text=t)))
        max=max(eval(parse(text=t)))
        sliderInput("tsneb2lim", label = h5("Select Range"), min = min,max =max, value =c(min,max)) 
      }
    }else if(input$categoryb2=="geneexp"){
      validate(need(input$categoryb2!="geneexp","Cannot subselect gene expression values"))
    }
  })
  
  #Dimensionality reduction options for left plot
  output$umapa = renderUI({
    scrna=fileload()
    dimr=names(scrna@reductions)
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
    selectInput("umapa","Dimensionality Reduction",dimr,selected = "umap")})
  })
  
  #Dimensionality reduction options for right plot
  output$umapb = renderUI({
    scrna=fileload()
    dimr=names(scrna@reductions)
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
    selectInput("umapb","Dimensionality Reduction",dimr,selected = "umap")})
  })
  
  #Genelist for tsne plot A
  output$gene1aui = renderUI({
    scrna=fileload()
    options=sort(rownames(GetAssayData(object = scrna)))
    withProgress(session = session, message = 'Generating gene list...',detail = 'Please Wait...',{
      selectInput('gene1a', label='Gene Name',options,multiple=FALSE, selectize=TRUE,selected=options[1])})
  })
  
  #Genelist for tsne plot B
  output$gene2aui = renderUI({
    scrna=fileload()
    options=sort(rownames(GetAssayData(scrna)))
    withProgress(session = session, message = 'Generating gene list...',detail = 'Please Wait...',{
      selectInput('gene2a', label='Gene Name',options,multiple=FALSE, selectize=TRUE,selected=options[1])})
  }) 
  
  #Based on all use input, generate plots using the right category and dimensionality reduction methods
  comptsne2 = reactive({
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data)
    met= sapply(metadata,is.numeric)
    tsnea=input$tsnea2
    tsneb=input$tsneb2
    feature=names(met[met==TRUE])
    tsne=names(met[met==FALSE])
    
    if(input$categorya2 =="clust" & input$subsa==F){
      plot1=DimPlot(object = scrna,reduction=input$umapa,group.by = "ident",no.legend = FALSE,label = input$checklabel1,vector.friendly = T, do.return=T, pt.size = input$pointa2,label.size = 7, cols=cpallette)
    }else if(input$categorya2 =="clust" & input$subsa==TRUE){
      cells=names(Idents(object=scrna)[Idents(object=scrna)==input$selclust])
      plot1=DimPlot(object = scrna,reduction=input$umapa,cells.highlight=cells,group.by = "ident",vector.friendly = T,no.legend = FALSE,label = F, do.return=T, pt.size = input$pointa2, cols=cpallette)
    }else if(input$categorya2=="geneexp"){
      validate(need(input$gene1a %in% rownames(GetAssayData(object=scrna)),"Incorrect Gene name.Gene names are case-sensitive.Please check for typos."))
      plot1=FeaturePlot(object = scrna,reduction=input$umapa, features = input$gene1a, cols = c(input$genecolor1, input$genecolor2),pt.size = input$pointa2)
      #plot1=eval(parse(text=paste("plot1$`",input$gene1a,"`",sep="")))
    }else if(input$categorya2 =="var" & input$tsnea2 %in% tsne & input$subsa==FALSE){
      plot1=DimPlot(object = scrna,reduction=input$umapa,group.by = tsnea,no.legend = FALSE,label = input$checklabel1,vector.friendly = T, do.return=T,pt.size = input$pointa2,label.size = 7, cols=cpallette)
    }else if(input$categorya2 =="var" & input$tsnea2 %in% tsne & input$subsa==TRUE){
      t=paste("rownames(scrna@meta.data[scrna@meta.data$",input$tsnea2,"==\"",input$selclust2,"\",])",sep="")
      cells=eval(parse(text=t))
      plot1=DimPlot(object = scrna,reduction=input$umapa,group.by = tsnea,cells.highlight=cells,vector.friendly = T,no.legend = FALSE,label =F, do.return=T,pt.size = input$pointa2, cols=cpallette)
    }else if(input$categorya2 =="var" & input$tsnea2 %in% feature & input$subsa==FALSE){
      plot1=FeaturePlot(object = scrna,reduction=input$umapa, features = tsnea, cols = c(input$genecolor1, input$genecolor2),pt.size = input$pointa2)
      #plot1=eval(parse(text=paste("plot1$`",tsnea,"`",sep="")))
    }else if(input$categorya2 =="var" & input$tsnea2 %in% feature & input$subsa==TRUE){
      t=paste('rownames(scrna@meta.data[scrna@meta.data$',input$tsnea2, '>',input$tsnea2lim[1], ' & metadata$',input$tsnea2, '<', input$tsnea2lim[2],',])',sep="")
      cells=eval(parse(text=t))
      plot1=FeaturePlot(object = scrna,reduction=input$umapa, features = tsnea,cells.use = cells, cols = c(input$genecolor1, input$genecolor2),pt.size = input$pointa2)
      #plot1=eval(parse(text=paste("plot1$`",tsnea,"`",sep="")))
    }
    
    if(input$categoryb2 =="clust" & input$subsb==F){
      plot2=DimPlot(object = scrna,reduction=input$umapb,group.by = "ident",no.legend = FALSE,vector.friendly = T,label = input$checklabel2, do.return=T,pt.size = input$pointa2,label.size = 7, cols=cpallette)
    }else if(input$categoryb2 =="clust" & input$subsb==TRUE){
      cells=names(Idents(object=scrna)[Idents(object=scrna)==input$selclustb])
      plot2=DimPlot(object = scrna,reduction=input$umapb,cells.highlight=cells,group.by = "ident",vector.friendly = T,no.legend = FALSE,label = F, do.return=T, pt.size = input$pointa2, cols=cpallette)
    }else if(input$categoryb2=="geneexp"){
      validate(need(input$gene2a %in% rownames(GetAssayData(object=scrna)),"Incorrect Gene name.Gene names are case-sensitive.Please check for typos."))
      plot2=FeaturePlot(object = scrna,reduction=input$umapb, features = input$gene2a, cols = c(input$genecolor1, input$genecolor2),pt.size = input$pointa2)
      #plot2=eval(parse(text=paste("plot2$`",input$gene2a,"`",sep="")))
    }else if(input$categoryb2 =="var" & input$tsneb2 %in% tsne & input$subsb==F){
      plot2=DimPlot(object = scrna,reduction=input$umapb,group.by = tsneb,no.legend = FALSE,vector.friendly = T,label = input$checklabel2, do.return=T,pt.size = input$pointa2,label.size = 7, colse=cpallette)
    }else if(input$categoryb2 =="var" & input$tsneb2 %in% tsne & input$subsb==TRUE){
      t=paste("rownames(scrna@meta.data[scrna@meta.data$",input$tsneb2,"==\"",input$selclustb2,"\",])",sep="")
      cells=eval(parse(text=t))
      plot2=DimPlot(object = scrna,reduction=input$umapb,group.by = tsneb,cells.highlight=cells,vector.friendly = T,no.legend = FALSE,label = F, do.return=T,pt.size = input$pointa2, cols=cpallette)
    }else if(input$categoryb2 =="var" & input$tsneb2 %in% feature & input$subsb==F){
      plot2=FeaturePlot(object = scrna,reduction=input$umapb, features = tsneb, cols = c(input$genecolor1, input$genecolor2),pt.size = input$pointa2)
      #plot2=eval(parse(text=paste("plot2$`",tsneb,"`",sep="")))
    }else if(input$categoryb2 =="var" & input$tsneb2 %in% feature & input$subsb==TRUE){
      t=paste('rownames(scrna@meta.data[scrna@meta.data$',input$tsneb2, '>',input$tsneb2lim[1], ' & metadata$',input$tsneb2, '<', input$tsneb2lim[2],',])',sep="")
      cells=eval(parse(text=t))
      plot2=FeaturePlot(object = scrna,reduction=input$umapb, features = tsneb,cells.use = cells, cols = c(input$genecolor1, input$genecolor2),pt.size = input$pointa2)
      #plot2=eval(parse(text=paste("plot2$`",tsneb,"`",sep="")))
    }
    
    p=plot_grid(plot1,plot2)
    #p2 <- add_sub(p, paste(input$projects,"_CompareTsne",sep=""), x = 0.87,vpadding = grid::unit(1, "lines"),size=11)
    p2 <- add_sub(p, paste(projectname(),"_CompareTsne",sep=""), x = 0.87,vpadding = grid::unit(1, "lines"),size=11)
    ggdraw(p2)
    ggdraw(p2)
  })
  
  #render final plot
  output$comptsne2 = renderPlot({
    input$load
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      comptsne2()
    })
  })
  
  #Handler to download plot
  output$downloadtsneplot <- downloadHandler(
    filename = function() {
      paste0(projectname(),"_CompareTsne.pdf",sep="")
    },
    content = function(file){
      pdf(file,width=14,height = 8,useDingbats=FALSE)
      plot(comptsne2())
      dev.off()
    })
  
  ###################################################
  ###################################################
  ########### Interactive Tsne plots  ###############
  ###################################################
  ###################################################
  #Generate drop down menu for variables from meta-data for right plot
  output$intervar = renderUI({
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data) 
    #metadata=metadata %>% select(starts_with("var"))
    var=colnames(metadata)
    selectInput("intervar","Select a Variable",var,"pick one")
  })
  
  #Generate drop down menu for categoried starting with Var from meta-data for left plot
  output$setcategory = renderUI({
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data)
    metadata=metadata %>% dplyr::select(starts_with("var_"))
    var=c(colnames(metadata))
    selectInput("setcategory","Choose category",var,"pick one")
  })
  
  #Generate drop down menu for Dimensionality reduction options for both plots
  output$umapint = renderUI({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
    scrna=fileload()
    dimr=names(scrna@reductions)
    selectInput("umapint","Dimensionality Reduction",dimr,selected = "umap")})
  })
  
  #Genelist for tsne plot 
  output$geneinterui = renderUI({
    scrna=fileload()
    options=sort(rownames(GetAssayData(object=scrna)))
    withProgress(session = session, message = 'Generating gene list...',detail = 'Please Wait...',{
      selectInput('geneinter', label='Gene Name',options,multiple=FALSE, selectize=TRUE,selected=options[1])})
  })
  
 #Render the tsne plot using plotly
  output$intertsne = renderPlotly({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      pdf(NULL)
      scrna=fileload()
      plot1=DimPlot(object = scrna,reduction=input$umapint,group.by = input$setcategory,no.legend = FALSE,label = TRUE, do.return=T,pt.size = input$umap_pointsize,label.size = 5,vector.friendly = T, cols=cpallette)
      plot=ggplotly(plot1)
      dev.off()
      return(plot)
    })
  })

  #Render the gene plot using plotly
  output$intergene = renderPlotly({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      pdf(NULL)
      scrna=fileload()
      metadata=as.data.frame(scrna@meta.data)
      met= sapply(metadata,is.numeric)
      #metadata=metadata %>% select(starts_with("var"))
      tsnea=input$intervar
      feature=names(met[met==TRUE])
      #feature=c("nGene","nUMI","percent.mito","S.Score","G2M.Score","var.ratio.pca")
      tsne=names(met[met==FALSE])
      
      if(input$intercat=="geneexp"){
        validate(need(input$geneinter %in% rownames(GetAssayData(object=scrna)),"Incorrect Gene name.Gene names are case-sensitive.Please check for typos."))
        plot1=FeaturePlot(object = scrna,reduction=input$umapint, features = input$geneinter, cols = c("grey", "blue"),pt.size = input$umap_pointsize)
        #plot1=eval(parse(text=paste("plot1$`",input$geneinter,"`",sep="")))
      }else if(input$intercat =="var" & tsnea %in% tsne){
        plot1=DimPlot(object = scrna,reduction=input$umapint,group.by = tsnea,no.legend = FALSE,label = TRUE, do.return=T,pt.size = input$umap_pointsize,label.size = 7, cols=cpallette,vector.friendly = T)
      }else if(input$intercat =="var" & tsnea %in% feature){
        plot1=FeaturePlot(object = scrna,reduction=input$umapint, features = tsnea, cols = c("grey", "blue"),pt.size = input$umap_pointsize)
        #plot1=eval(parse(text=paste("plot1$`",tsnea,"`",sep="")))
      }
      plot=ggplotly(plot1)
      dev.off()
      return(plot)
    })
  })
  
  ######################################################################################################
  ######################################################################################################
  ####### Display Biplot plot with controls ############################################################
  ######################################################################################################
  ######################################################################################################
  
  #generate Expression limit for gene A
  output$bigene_rangea <- renderUI({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      #textInput("bigene_genea", label = "Gene A",value = bigene_genea)
      r<-getGeneRange(fileload(),input$bigene_genea)
      sliderInput("bigene_rangea", "Expression Limit Gene A (log2 UMI)",
                  min = 0, max = r[2], value = c(r[1],r[2]),step=.25)
    })
  })
  
  #generate Expression limit for gene B
  output$bigene_rangeb <- renderUI({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      #textInput("bigene_geneb", label = "Gene B",value = bigene_geneb)
      r<-getGeneRange(fileload(),input$bigene_geneb)
      sliderInput("bigene_rangeb", "Expression Limit Gene B (log2 UMI)",
                  min = 0, max = r[2], value = c(r[1],r[2]),step=.25)
    })
  })
  
  #Generate dropdown to pick dimensionality reduction method for bigeneplot
  output$bigenedim = renderUI({
    scrna=fileload()
    dimr=names(scrna@reductions)
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      selectInput("bigenedim","Dimensionality Reduction",dimr,selected = "umap")})
  })
  
  #Genelist A for bigene plot 
  output$bigene_geneaui = renderUI({
    scrna=fileload()
    options=sort(rownames(GetAssayData(object = scrna)))
    withProgress(session = session, message = 'Generating gene list...',detail = 'Please Wait...',{
      selectInput('bigene_genea', label='Gene Name',options,multiple=FALSE, selectize=TRUE,selected=options[1])})
  })
  
  #Genelist B for bigene plot 
  output$bigene_genebui = renderUI({
    scrna=fileload()
    options=sort(rownames(GetAssayData(object = scrna)))
    withProgress(session = session, message = 'Generating gene list...',detail = 'Please Wait...',{
      selectInput('bigene_geneb', label='Gene Name',options,multiple=FALSE, selectize=TRUE,selected=options[2])})
  })
  
#plot the bi-gene plot
  output$bigeneplot <- renderPlot({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      p=bigene_plot(fileload(),
                    c(input$bigene_genea,input$bigene_geneb),
                    limita=input$bigene_rangea,
                    limitb=input$bigene_rangeb,
                    marker_size = input$bigene_pointsize,
                    type=input$bigenedim)
      p2 <- add_sub(p, paste(projectname(),"_Bigeneplot",sep=""), x = 0.87,vpadding = grid::unit(1, "lines"),size=11)
      ggdraw(p2)
    })
  })
  
  # Generate cluster-wise counts for the two genes
  genecounts <- reactive({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
    scrna=fileload()
    validate(need(input$bigene_genea %in% rownames(GetAssayData(object = scrna)),"Incorrect Gene name.Gene names are case-sensitive.Please check for typos."))
    validate(need(input$bigene_geneb %in% rownames(GetAssayData(object = scrna)),"Incorrect Gene name.Gene names are case-sensitive.Please check for typos."))
    genes=c(input$bigene_genea,input$bigene_geneb)
    my.data=FetchData(scrna,c("ident",genes))
    colnames(my.data)=c("ident","gene1","gene2")
    my.data= my.data %>% group_by(ident) %>% summarize(Gene1_ct=sum(gene1>0), Gene2_ct=sum(gene2 > 0),Both_ct=sum(gene2 > 0 & gene1>0)) 
    colnames(my.data)=c("Cell Group",genes[1],genes[2],"both")
    return(my.data)
    })
  })
  
  #Generate table for cluster wise gene counts
  output$bigene_genecount = DT::renderDataTable({
    input$bigene_genea
    input$bigene_geneb
    DT::datatable(genecounts(),
                  extensions = 'Buttons', options = list(
                    dom = 'Bfrtip',
                    pageLength = 10,
                    lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                    buttons = list()),
                  rownames=FALSE,selection = list(mode = 'single', selected =1),escape=FALSE)
    
  })
  
  
  #Download bi-gene plot
  output$downloadbiplot <- downloadHandler(
    filename = function() {
      paste0(projectname(),"_",input$bigene_genea,"_",input$bigene_geneb,"_Bigene.pdf",sep="")
    },
    content = function(file){
      pdf(file,width=9,height = 9,useDingbats=FALSE)
      plot(bigene_plot(fileload(),
                       c(input$bigene_genea,input$bigene_geneb),
                       limita=input$bigene_rangea,
                       limitb=input$bigene_rangeb,
                       marker_size = input$bigene_pointsize,type=input$bigenedim))
      dev.off()
    })
  
  ######################################################################################################
  ######################################################################################################
  ######################### Display 3D plot ############################################################
  ######################################################################################################
  ######################################################################################################
  #generate dropdown for variable
  output$var3d <- renderUI({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data) 
    metadata=metadata %>% select(starts_with("var"))
    var=colnames(metadata)
    selectInput("var3d","Select a Variable",var,"pick one")
    })
  })
  
  #Generate list of reductions
  output$dimr3d <- renderUI({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
    scrna=fileload()
    opt=names(scrna@reductions)
    selectInput("dimr3d", "Select Reduction",opt,selected = "umap")
      })
  })
  
  #Create 3D plot
  output$plot3d = renderPlotly({
   withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
    scrna.sub=fileload()
    reduction=input$dimr3d
    groupby=input$var3d
    maxdim <- getMaxDim(scrna.sub)
    if(reduction=='umap'){
      scrna.sub <- RunUMAP(scrna.sub,dims=1:maxdim,n.components = 3)
    }else if(reduction=='tsne'){
      scrna.sub <- RunTSNE(scrna.sub,dims=1:maxdim,dim.embed= 3)
    }
    
    dims=1:3
    dims <- paste0(Key(object = scrna.sub[[reduction]]), dims)
    data <- FetchData(object = scrna.sub, vars = c(dims,groupby))
    colnames(data)[1:4] = c("DM_1","DM_2","DM_3","var_cluster")
    a3=aggregate(data$DM_3, by=list(data$var_cluster), FUN=mean)
    a2=aggregate(data$DM_2, by=list(data$var_cluster), FUN=mean)
    a1=aggregate(data$DM_1, by=list(data$var_cluster), FUN=mean)
    centers=inner_join(a1,a2,by="Group.1")
    centers=inner_join(centers,a3,by="Group.1")
    colnames(centers)=c("var_cluster","x","y","z")
    
    a <- list()
    for (i in 1:nrow(centers)) {
      a[[i]] <- list(x= centers$x[i],y= centers$y[i],z= centers$z[i],text= centers$var_cluster[i],showarrow= T,arrowhead=4,arrowsize=0.5)
    }
    if(input$check3d == T){
      validate(
        need(is.na(scrna.sub@misc$sds)==F,"Lineage curve information not found. Please run slingshot on the dataset and upload to website again")
      )
      curved <- bind_rows(lapply(names(scrna.sub@misc$sds$data@curves), function(x){c <- slingCurves(scrna.sub@misc$sds$data)[[x]]
      d <- as.data.frame(c$s[c$ord,seq_len(2)])
      d$curve<-x
      return(d)}))
      colnames(data)[1:3] = c("DM_1","DM_2","DM_3")
      plot=plot_ly(side=I(3)) %>%
        add_trace(x = data$DM_1,y = data$DM_2,z = data$DM_3,colors=cpallette,color=data$var_cluster,type = "scatter3d") %>% 
        add_paths(x = curved$DM_1,y = curved$DM_2,z = curved$DM_3, mode="lines",color=I("black"),size=I(7)) %>% 
        layout(scene = list(
          aspectratio = list(x = 1,y = 1,z = 1),
          dragmode = "turntable",
          xaxis = list(title = dims[1]),yaxis = list(title = dims[2]),zaxis = list(title = dims[3]),annotations = a))
    }else{
    plot=plot_ly(side=I(3)) %>%
      add_trace(x = data$DM_1,y = data$DM_2,z = data$DM_3,colors=cpallette,color=data$var_cluster,type = "scatter3d") %>% 
      #add_paths(x = curved$DM_1,y = curved$DM_2,z = curved$DM_3, mode="lines",color=I("black"),size=I(7)) %>% 
      layout(scene = list(
          aspectratio = list(x = 1,y = 1,z = 1),
          dragmode = "turntable",
          xaxis = list(title = dims[1]),yaxis = list(title = dims[2]),zaxis = list(title = dims[3]),annotations = a))
    }
    plot
    })
  })
  ####################################################
  ###################################################
  ########## Setup Control Panel for DEG ############
  ###################################################
  ###################################################
  #Generate drop down for dimensionality reduction in DEG tab
  output$umapdeg = renderUI({
    scrna=fileload()
    dimr=names(scrna@reductions)
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      selectInput("umapdeg","Dimensionality Reduction",dimr,selected = "umap")})
  })
  
  #Generate drop down menu for Cell group/ other variables to be displayed in the DEG tab
  output$tsnea = renderUI({
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data) 
    #metadata=metadata %>% select(starts_with("var"))
    var=c(colnames(metadata),'Cell.group')
    selectInput("tsnea","Select Group to display",var,selected = "Cell.group")
  })
  
  #Generate drop down menu for the default ident/cluster/variable of comparison
  output$identdef = renderUI({
    scrna=fileload()
    options=sort(unique(scrna@misc$findallmarkers$cluster))
    selectInput("identdef", "First cluster/variable of comparison (Cell Group 1)",options)
  })
  
  #Generate drop down menu for the categories starting with "var" that can be set as the new ident/variable of comparison
  output$setidentlist = renderUI({
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data)
    metadata=metadata %>% dplyr::select(starts_with("var_"))
    var=c(colnames(metadata))
    selectInput("setidentlist","Choose category to compare",var,"pick one")
    
  })
  
  #Generate drop down menu for the variables in the new ident against which the rest will be compared to find markers
  output$identa = renderUI({
    scrna=fileload()
    if(input$setident==T){
      Idents(object = scrna) = input$setidentlist
      options=sort(unique(Idents(object=scrna)))
    }else{
      options=sort(levels(Idents(object=scrna)))
    }
    selectInput("identa", "First Cell group to compare",options)
  })
  
  #Generate checkboxes for the variables in the new ident. Choose all those that you want to compare with the first 
  output$identb = renderUI({
    scrna=fileload()
    if(input$setident==T){
      Idents(object = scrna) = input$setidentlist
      options=sort(unique(Idents(object=scrna)))
    }else{
      options=sort(levels(Idents(object=scrna)))
    }
    checkboxGroupInput("identb", label="Second Cell group to compare",choices=options)
  })
  
  
  ###################################################
  ###################################################
  ####### Display DEG plot with controls  ###########
  ###################################################
  ###################################################
  
  #Based on the input selected from the control panel, run Seurat's findMarkers to find marker genes that distinguish one
  #chosen variable of a chosen category from other(s)
  markergenes = reactive({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      scrna=fileload()
      if(input$setident==T){ #set new ident if user choses to change category
        validate(need(input$goButton != 0,"Make your selections and click the GO button"))
        
        if(input$goButton == 0)
          return()
        isolate({
          Idents(object = scrna)  = input$setidentlist #set ident
          validate(need(input$identb,"Select at least one option from Second cell group to compare to first cell group. If you want to compare to all, uncheck the 'Check to choose a different category to compare' option"))
          validate(need(input$identb!=input$identa,"First and second cell groups can't be the same"))
            identb=input$identb
            p=unlist(strsplit(identb,","))
            markers=FindMarkers(object = scrna, ident.1 = input$identa, ident.2 = p, min.pct = input$minpct,logfc.threshold=input$lfc,test.use=input$test)
            geneid=rownames(markers)
            url= paste("http://www.genecards.org/cgi-bin/carddisp.pl?gene=",geneid,sep = "")
            markers$Link=paste0("<a href='",url,"'target='_blank'>",rownames(markers),"</a>")
        })
      }
      if(input$setident==F){
        markers=scrna@misc$findallmarkers
        markers=markers[markers$cluster==input$identdef,]
        geneid=rownames(markers)
        url= paste("http://www.genecards.org/cgi-bin/carddisp.pl?gene=",geneid,sep = "")
        markers$Link=paste0("<a href='",url,"'target='_blank'>",rownames(markers),"</a>")
      }
    })
    colnames(markers)[4]="Percentage expressed in Cell Group 1"
    colnames(markers)[5]="Percentage expressed in Cell Group 2"
    return(markers)
  })
  
  #Display results in a table
  output$markergenes = DT::renderDataTable({
    #input$identa
    #input$identb
    input$goButton
    DT::datatable(markergenes(),
                  extensions = 'Buttons', options = list(
                    dom = 'Bfrtip',
                    pageLength = 10,
                    lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                    buttons = list()),
                  rownames=TRUE,selection = list(mode = 'single', selected =1),escape=FALSE)
    
  })
  
  #Download function to download the table of marker genes
  output$downloaddeg <- downloadHandler(
    filename = function() { paste(projectname(), '.csv', sep='') },
    content = function(file) {
      write.csv(markergenes(), file)
    })
  
  ###################################################
  ###################################################
  ### View Differentially expressed marker genes  ###
  ###################################################
  ###################################################
  
  #Generate the plots in deg tab
  #plot1 is the tsne/umap/pca/phate/diffusion map (based on chosen dimensionality reduction)
  #plot 2 is feature plot of gene selected from table
  #plot 2 is violin of gene selected from table
  #plot 2 is ridge plot of gene selected from table
  comptsne = reactive({
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data)
    met= sapply(metadata,is.numeric)
    validate(need(is.null(scrna@meta.data$var_cluster)==F,"var_cluster not found in meta data"))
    scrna@meta.data$var_cluster=as.numeric(as.character(scrna@meta.data$var_cluster))
    tsnea=input$tsnea
    feature=names(met[met==TRUE])
    tsne=names(met[met==FALSE])
    
    if(input$tsnea =="Cell.group"){
      plot1=DimPlot(object = scrna,reduction=input$umapdeg,no.legend = FALSE,label = input$checklabel3, do.return=T, pt.size = input$pointa,label.size = 7,cols=cpallette,vector.friendly=TRUE) + theme(legend.position="bottom")
    }else if(input$tsnea %in% tsne){
      plot1=DimPlot(object = scrna,reduction=input$umapdeg,group.by = tsnea,no.legend = FALSE,label = input$checklabel3,vector.friendly=TRUE, do.return=T,pt.size = input$pointa,label.size = 7,cols=cpallette) + theme(legend.position="bottom")
    }else if(input$tsnea %in% feature){
      plot1=FeaturePlot(object = scrna, features = tsnea, cols = c("grey", "blue"),reduction = input$umapdeg,pt.size = input$pointa,combine = T)
      #plot1=eval(parse(text=paste("plot1$`",tsnea,"`",sep="")))
    }
   
    markers=markergenes()
      s=input$markergenes_rows_selected # get  index of selected row from table
      markers=markers[s, ,drop=FALSE]
      plot2=FeaturePlot(object = scrna, features = rownames(markers), cols = c("grey","blue"),reduction = input$umapdeg,pt.size = input$pointa)
      #plot2=eval(parse(text=paste("plot2$`",rownames(markers),"`",sep="")))
  if(input$setident==T){
        setident=input$setidentlist  
      if(input$checkviolin ==T){
      plot3=VlnPlot(object = scrna, features = rownames(markers),group.by = setident,pt.size=0,cols=cpallette)
      }else{plot3=VlnPlot(object = scrna, features = rownames(markers),group.by = setident,cols=cpallette)}
      plot4=RidgePlot(object = scrna, features = rownames(markers),group.by = setident,cols=cpallette)
  }else{
    if(input$checkviolin ==T){
      plot3=VlnPlot(object = scrna, features = rownames(markers),pt.size=0,cols=cpallette)
    }else{plot3=VlnPlot(object = scrna, features = rownames(markers),cols=cpallette)}
    plot4=RidgePlot(object = scrna, features = rownames(markers),cols=cpallette)
  }
    
      row1=plot_grid(plot1,plot2,align = 'h', rel_heights = c(1, 1),axis="lr", nrow=1)
      row2=plot_grid(plot3,plot4,align = 'h', rel_heights = c(1, 1),axis="lr", nrow=1)
    p=plot_grid(row1,row2,align = 'v', rel_heights = c(1.7, 1),axis="tb",ncol=1)
    p2 <- add_sub(p, paste(projectname(),"_Differential_Exp",sep=""), x = 0.87,vpadding = grid::unit(1, "lines"),size=11)
    ggdraw(p2)
  })
  
  #Render the plot to display
  output$comptsne = renderPlot({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      comptsne()
    })
  })
  
  #Download plot function
  output$downloadplot <- downloadHandler(
    filename = function() {
      markers=markergenes()
      s=input$markergenes_rows_selected # get  index of selected row from table
      markers=markers[s, ,drop=FALSE]
      paste0(projectname(),"_",rownames(markers),"_DEG.pdf",sep="")
    },
    content = function(file){
      pdf(file, width = 12, height = 11,useDingbats=FALSE)
      plot(comptsne())
      dev.off()
    })
  
  ###################################################
  ###################################################
  ####### Display Heatmap plot with controls  #######
  ###################################################
  ###################################################
  #get the list of the differentially expressed markers from the differential expression tab and compute min and max
  #number of genes to show in the dotpot
  output$heatmapgenes = renderUI({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      markers=markergenes()
      validate(
        need(nrow(markers)>0, "No Marker genes found")
      )
      if(nrow(markers)<10){
        min=1
        max=nrow(markers)
      }else{
        min=10
        max=nrow(markers)
      }
      sliderInput("heatmapgenes", "Number of top genes to plot:",min = min, max = max,value = min)
    })
  })
  
  #Generate drop-down for list of variables to group cells by
  output$hmpgrp = renderUI({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      scrna=fileload()
      metadata=as.data.frame(scrna@meta.data)
      metadata=metadata %>% dplyr::select(starts_with("var_"))
      var=c("ident",colnames(metadata))
      selectInput("hmpgrp","Select a Variable",var,"pick one")
    })
  })
  
  #generate the heatmap
  heatmap <- reactive({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      scrna=fileload()
      if(input$shmptype =="deggene"){
      markers=markergenes()
      markergenes=rownames(markers)[1:input$heatmapgenes]
      }else if(input$shmptype =="topgene"){
        markers <- FindAllMarkers(object = scrna, only.pos = TRUE, min.pct = 0.25,thresh.use = 0.25)
        markers %>% group_by(cluster) %>% top_n(input$topn, avg_logFC)
        markergenes=markers$gene
      }
#       if(input$hmpcol=="PuYl"){
#         lowcol="darkmagenta"
#         midcol="black"
#         highcol="yellow"
#       }else if(input$hmpcol=="BuGn"){
#         lowcol="yellow"
#         midcol="green"
#         highcol="blue"
#       }else if(input$hmpcol=="RdYl"){
#         lowcol="yellow"
#         midcol="red"
#         highcol="black"
#       }else if(input$hmpcol=="RdBu"){
#         lowcol="red"
#         midcol="white"
#         highcol="blue"}
      p=DoHeatmap(object = scrna, features = markergenes,group.by = input$hmpgrp, group.bar= T,label=TRUE)
      p2 <- add_sub(p, paste(projectname(),"_Heatmap",sep=""), x = 0.87,vpadding = grid::unit(1, "lines"),size=11)
      ggdraw(p2)
    })
  })
  
  #Render the heatmap
  output$heatmap <- renderPlot({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      heatmap()
    })
  })
  
  #Download function for the heatmap
  output$downloadheatmap <- downloadHandler(
    filename = function() {
      paste0(projectname(),"_Heatmap",sep="")
    },
    content = function(file){
      pdf(file, width = 13, height = 8,useDingbats=FALSE)
      plot(heatmap())
      dev.off()
    })
  
  ###################################################
  ###################################################
  ########### Plot gene expression  ################
  ###################################################
  ###################################################
  
  #Generate drop down for dimensionality reduction in GeneExpression Plots tab
  output$umapge = renderUI({
    scrna=fileload()
    dimr=names(scrna@reductions)
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      selectInput("umapge","Dimensionality Reduction",dimr,selected = "umap")})
  })
  
  #Genelist for Gene Expression Plot
  output$geneidui = renderUI({
    scrna=fileload()
    options=sort(rownames(GetAssayData(object=scrna)))
    withProgress(session = session, message = 'Generating gene list...',detail = 'Please Wait...',{
      selectInput('geneid', label='Gene Name',options,multiple=FALSE, selectize=TRUE,selected=options[1])})
  })
  
  #For any gene entered, generate a feature plot, violin plot and a ridge plot
  geplots = reactive({
    scrna=fileload()
    validate(need(input$geneid,"Enter the gene symbol"))
    validate(need(input$geneid %in% rownames(GetAssayData(object=scrna)),"Incorrect Gene name.Gene names are case-sensitive.Please check for typos."))
    plot2=FeaturePlot(object = scrna, features = input$geneid, cols = c("grey","blue"),reduction = input$umapge,pt.size = input$genenid_pointsize)
    #plot2=eval(parse(text=paste("plot2$`",input$geneid,"`",sep="")))
    if(input$checkviolin2 ==T){
      plot3=VlnPlot(object = scrna, features = input$geneid,pt.size=0,cols=cpallette)
    }else{plot3=VlnPlot(object = scrna, features = input$geneid,cols=cpallette)}
    plot4=RidgePlot(object = scrna, features = input$geneid,cols=cpallette)
    row1=plot_grid(plot2,align = 'h', rel_heights = c(1, 1),axis="lr", nrow=1)
    row2=plot_grid(plot3,plot4,align = 'h', rel_heights = c(1, 1),axis="lr", nrow=1)
    plot_grid(row1,row2,align = 'v', rel_heights = c(1.7, 1),axis="tb",ncol=1)
  })
  
  #Render the plot
  output$geplots = renderPlot({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      geplots()
    })
  })
  
  #Download the plots from gene expression plots tab
  output$downloadplotge <- downloadHandler(
    filename = function() {
      paste0(input$geneid,"_Geneexp_plot.pdf",sep="")
    },
    content = function(file){
      pdf(file,width=8,height = 13,useDingbats=FALSE)
      plot(geplots())
      dev.off()
    })
  
  ###################################################
  ###################################################
  ########### Plot gene expression in clusters  #####
  ###################################################
  ###################################################
  
  #Generate drop down for dimensionality reduction in Cluster-wise gene expression
  output$umapclust = renderUI({
    scrna=fileload()
    dimr=names(scrna@reductions)
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      selectInput("umapclust","Dimensionality Reduction",dimr,selected = "umap")})
  })
  
  #Generate drop down for categories starting with "var" in Cluster-wise gene expression
  output$setvar = renderUI({
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data)
    metadata=metadata %>% dplyr::select(starts_with("var_"))
    var=c(colnames(metadata))
    selectInput("setvar","Choose category",var,"pick one")
  })
  
  #Generate drop down for unique variable within the chose category in Cluster-wise gene expression
  output$selectcluster = renderUI({
    scrna=fileload()
    options <- paste0("scrna@meta.data$",input$setvar,sep="")
    options=unique(eval(parse(text=options)))
    options=options[order(options)]
    selectInput("selectcluster", "Select Cell group",options)
  })
  
  
  #Set ident to the chosen category and calculate average expression of genes per category and the percentage of cells
  #they are expressed in
  clusts= reactive({
    scrna=fileload()
    Idents(object = scrna) = input$setvar
    avgexp=AverageExpression(object = scrna)
    avgexp= avgexp %>% dplyr::select(input$selectcluster)
    genes.use=rownames(avgexp)
    data.use <- GetAssayData(object = scrna,slot = "data")
    cells <- WhichCells(object = scrna, ident = input$selectcluster)
    thresh.min=0
    data.temp <- round(x = apply(X = data.use[genes.use, cells, drop = F],
                                 MARGIN = 1,
                                 FUN = function(x) {
                                   return(sum(x > thresh.min) / length(x = cells))
                                 }),digits = 3)
    names(data.temp)=genes.use
    data.temp=as.data.frame(data.temp)
    data.temp$id=rownames(data.temp)
    avgexp$id=rownames(avgexp)
    df=inner_join(avgexp,data.temp,by="id") 
    rownames(df)=df$id
    df= df %>% dplyr::select(input$selectcluster,data.temp) 
    df=df[order(-df[,1],-df[,2]),]
    colnames(df)= c("Average Expression","Percentage of cells expressed in")
    df$max_avg=max(df$`Average Expression`)
    df$min_avg=min(df$`Average Expression`)
    return(df)
  })
  
  #use above table to compute min and maximum average expression and put it in a slider
  output$avgexpslider <- renderUI({
    withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
      df=clusts()
      min=unique(df$min_avg)
      max=unique(df$max_avg)
      sliderInput("avgexpslider", "Average Expression:",min = min, max = max, value = c(min,max))
    })
  })
  
  #filter the cluster-wise gene expression table based on user-defined filters of avg expression and percentage expressed
  clustable= reactive({
    df=clusts()
    df= df %>% dplyr::select(-max_avg:-min_avg) 
    df=df[df$`Percentage of cells expressed in` >input$pctslider[1] & df$`Percentage of cells expressed in` <input$pctslider[2],]
    df=df[df$`Average Expression` >input$avgexpslider[1] & df$`Average Expression` <input$avgexpslider[2],]
    return(df)
  })
  
  #Display cluster-wise gene expression table
  output$clustable = DT::renderDataTable({
    input$setvar
    input$selectcluster
    input$avgexpslider
    input$pctslider
    input$umapclust
    withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
      DT::datatable(clustable(),
                    extensions = c('Buttons','Scroller'),
                    options = list(dom = 'Bfrtip',
                                   searchHighlight = TRUE,
                                   pageLength = 10,
                                   lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                   scrollX = TRUE,
                                   buttons = c('copy', 'print')
                    ),rownames=TRUE,caption= "Cluster-wise Gene expression",selection = list(mode = 'single', selected =1),escape = F)
    })
  })
  
  #Generate plots of the chosen gene
  clustplots= reactive({
    scrna=fileload()
    tab=clustable()
    s=input$clustable_rows_selected
    tab=tab[s, ,drop=FALSE]
    gene=rownames(tab)
    #cells <- WhichCells(object = scrna, ident = input$selectcluster)
    plot1=DimPlot(object = scrna,reduction=input$umapclust,group.by = input$setvar,no.legend = FALSE,label = input$checklabel4,vector.friendly = T, do.return=T, pt.size = input$pointclust,label.size = 7, cols=cpallette)
    plot2=FeaturePlot(object = scrna,reduction=input$umapclust, features = gene, cols = c("grey", "blue"),pt.size = input$pointclust)
    #plot2=eval(parse(text=paste("plot2$`",gene,"`",sep="")))
    plot_grid(plot1,plot2)
  })
  
  #Render above plot
  output$clustplots = renderPlot({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      clustplots()
    })
  })
  
  #Download plot
  output$downloadclustplot <- downloadHandler(
    filename = function() {
      paste0("Clustexp_plot.pdf",sep="")
    },
    content = function(file){
      pdf(file,width=13,height = 9,useDingbats=FALSE)
      plot(clustplots())
      dev.off()
    })
  
  #Download cluster-wise gene expression table
  output$downloadclustertab <- downloadHandler(
    filename = function() { paste(input$selectcluster, '_cluster-wiseGeneExp.csv', sep='') },
    content = function(file) {
      write.csv(clustable(), file)
    })
  
  ###################################################
  ###################################################
  ################### Plot dot plot ################
  ###################################################
  ###################################################
  
  #Generate drop-down for list of variables to group cells by
  output$setdotvar = renderUI({
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data)
    metadata=metadata %>% dplyr::select(starts_with("var_"))
    var=c(colnames(metadata))
    selectInput("setdotvar","Choose category",var,"pick one")
  })
  
  #Genelist for Dot plot
  output$enterchoice = renderUI({
    scrna=fileload()
    options=sort(rownames(GetAssayData(object=scrna)))
    withProgress(session = session, message = 'Generating gene list...',detail = 'Please Wait...',{
      selectInput('genelistfile2', label='Gene Name',options,multiple=TRUE, selectize=TRUE,selected=options[1])})
  })
  
  #read in user-defined list of genes and generate the dot plot
  dotplot= reactive({
    scrna=fileload()
    
    if(input$radiofileup=="upload"){
      validate(
        need(input$genelistfile, "Please Upload Genelist")
      )
    file=input$genelistfile
    df=fread(file$datapath,header = FALSE) #get complete gene list as string
    genes=as.vector(df$V1)
    }else if(input$radiofileup=="enter"){
      validate(
        need(input$genelistfile2, "Please Select Genes")
      )
      genes=input$genelistfile2
    }
    g1=DotPlot(object = scrna, features                                                          = genes, plot.legend = TRUE,group.by=input$setdotvar,do.return=TRUE) 
    return(g1) 
  })
  
  #render the dot plot
  output$dotplot = renderPlot({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      dotplot()
    })
  })
  
  #Download plot
  output$downloaddotplot <- downloadHandler(
    filename = function() {
      paste0(input$project,"_Dotplot.pdf",sep="")
    },
    content = function(file){
      pdf(file,width=10,height = 9,useDingbats=FALSE)
      plot(dotplot())
      dev.off()
    })
   ###################################################
   ###################################################
   ##### CONTROL PANEL FOR LIGAND-RECEPTOR PAIRS #####
   ###################################################
   ###################################################
   #Generate drop down for dimensionality reduction for bi-gene plot
   output$bigenedimr = renderUI({
     scrna=fileload()
     dimr=names(scrna@reductions)
     withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
       selectInput("bigenedimr","Dimensionality Reduction",dimr,selected = "umap")})
   })
   
   #Generate drop-down to generate variables based on which you want to find pairs 
   #This is important for choosing the 'select cluster' option
   output$pairby <- renderUI({
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
       scrna=fileload()
       metadata=as.data.frame(scrna@meta.data)
       metadata=metadata %>% dplyr::select(starts_with("var_"))
       options=colnames(metadata)
       selectInput("pairby","Select cell group ",options,selected=options[1])
     })
   })
   
  #If user chooses to select cluster and all genes, generate drop down menus to pick receptor cluster
   output$clust1 <- renderUI({
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
     scrna=fileload()
     t=paste("scrna@meta.data$",input$pairby,sep="")
     options=unique(eval(parse(text=t)))
     if(input$pairby=="ident"){options=levels(Idents(object=scrna))}
     selectInput("clust1","Pick cellgroup for Receptor",options,selected=options[1])
     })
   })
   
   #If user chooses to select cluster and all genes, generate drop down menus to pick ligand cluster
   output$clust2 <- renderUI({
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
     scrna=fileload()
     t=paste("scrna@meta.data$",input$pairby,sep="")
     options=unique(eval(parse(text=t)))
     if(input$pairby=="ident"){options=levels(Idents(object=scrna))}
     selectInput("clust2","Pick cellgroup for Ligand",options, selected=options[2])
     })
   })
   
   #If user chooses to select both cluster and genes, generate drop down menus to pick receptor cluster
   output$clust1.1 <- renderUI({
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
     scrna=fileload()
     t=paste("scrna@meta.data$",input$pairby2,sep="")
     options=unique(eval(parse(text=t)))
     selectInput("clust1.1","Pick cellgroup for Receptor",options,selected=options[1])
     })
   })
   
   #If user chooses to select both cluster and genes, generate drop down menus to pick ligand cluster
   output$clust2.1 <- renderUI({
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
     scrna=fileload()
     t=paste("scrna@meta.data$",input$pairby2,sep="")
     options=unique(eval(parse(text=t)))
     selectInput("clust2.1","Pick cellgroup for Ligand",options,selected=options[2])
     })
   })
   
   #generate Expression limit for Ligand gene
   output$bigene_rangea2 <- renderUI({
     withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
       table=finalres()
       s=input$pairs_res_rows_selected
       table=table[s, ,drop=FALSE]
       validate(need(nrow(table) > 0,"No Ligand-receptor pairs found"))
       bigene_genea=table$ligand
       #textInput("bigene_genea", label = "Gene A",value = bigene_genea)
       r<-getGeneRange(fileload(),bigene_genea)
       sliderInput("bigene_rangea", "Expression Limit of Ligand Gene (log2 UMI)",
                   min = 0, max = r[2], value = c(r[1],r[2]),step=.25)
     })
   })
   
   #generate Expression limit for receptor gene
   output$bigene_rangeb2 <- renderUI({
     withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
       table=finalres()
       validate(need(nrow(table) > 0,"No Ligand-receptor pairs found"))
       s=input$pairs_res_rows_selected
       table=table[s, ,drop=FALSE]
       bigene_geneb=table$receptor
       #textInput("bigene_geneb", label = "Gene B",value = bigene_geneb)
       r<-getGeneRange(fileload(),bigene_geneb)
       sliderInput("bigene_rangeb", "Expression Limit of Receptor Gene (log2 UMI)",
                   min = 0, max = r[2], value = c(r[1],r[2]),step=.25)
     })
   })
   
   #Check if the data is from mouse/human and use the approprite file to list the options for source
   output$source <- renderUI({
     file = read.csv("data/param.csv")
     if (input$filetype=="list"){
     org=as.character(file$organism[file$projects==input$projects])}
     else if(input$filetype == "upload"){
       scrna=fileload()
       validate(need(scrna@meta.data$org,"Organism not found in meta data. Cannot proceed"))
       org = unique(scrna@meta.data$org)
     }
     if(org=="mouse"){rl=read.csv("data/Mm_PairsLigRec.csv")}else if(org=="human"){rl=read.csv("data/Hs_PairsLigRec.csv")}
     options=as.character(unique(rl$Pair.Source))
     checkboxGroupInput('source', label='Select source(s)',choices=options,selected=options[2])
   })
   
   #Check if the data is from mouse/human and use the approprite file to list the options for evidence
   output$evidence <- renderUI({
     # file = read.csv("data/param.csv")
     # org=as.character(file$organism[file$projects==input$projects])
     # if(org=="mouse"){rl=read.csv("data/Mm_PairsLigRec.csv")}else if(org=="human"){rl=read.csv("data/Hs_PairsLigRec.csv")}
     # options=as.character(unique(rl$Pair.Evidence))
     options=c("EXCLUDED not ligand","literature supported","putative","EXCLUDED not receptor","EXCLUDED")
     checkboxGroupInput('evidence',label='Select Evidence(s)',choices=options,selected=options[2])
   })
   
   
   ###################################################
   ###################################################
   #### Load lig-receptor list and display results ###
   ###################################################
   ###################################################
   
   #For selected project and grouping variable, generate all possible ligand receptor pairs
   datasetInput = reactive({
#     results=ligrec(fileload(),pair=input$pairby,prj=projectname(),input$perc_cells,filetype=input$filetype)
     scrna=fileload()
     results =  scrna@misc$ligrecres
   })

      #Subselect lig-rec pairs based on user input
   finalres= reactive({
     validate(need(input$lrpgo != 0,"Make your selections and click the Run button"))
     if(input$lrpgo == 0)
       return()
     isolate({
     result=datasetInput()
     if(input$clust=="clust" & input$gene=="allgene"){
       #clusters=c(input$clust1,input$clust2)
       result=result[(result$Receptor_cluster %in% input$clust1) & (result$Lig_cluster%in% input$clust2),]
     }else if(input$clust=="clust" & input$gene=="genelist"){
       #clusters.1=c(input$clust1.1,input$clust2.1)
       result=result[(result$Receptor_cluster %in% input$clust1.1) & (result$Lig_cluster%in% input$clust2.1),]
     }else{result=result
     }
     
     if(input$gene=="genelist"){
       if(input$clust=="all"){
         g1=input$genelist1
         g2=input$genelist2
       }else if(input$clust=="clust"){
         g1=input$genelist1.1
         g2=input$genelist2.1
       }
       genes=read.table(g1$datapath,stringsAsFactors = F)#get complete gene list as string
       #to avoid errors due to case, change user-defined genelist to lowecase
       g1=as.vector(genes$V1)
       g1=tolower(g1)
       firstup <- function(x) {
         substr(x, 1, 1) <- toupper(substr(x, 1, 1))
         x
       }
       g1=firstup(g1)
       genes2=read.table(g2$datapath,stringsAsFactors = F)
       g2=as.vector(genes2$V1)
       g2=tolower(g2)
       g2=firstup(g2)
       result=result[(result$receptor %in% g1) & (result$ligand %in% g2),]
     }else{
       result=result
     }
     if(input$checksource==T){result=result[result$Pair.Source %in% input$source,]}
     if(input$checkevi==T){result=result[result$Pair.Evidence %in% input$evidence,]}
     result=result %>% dplyr::select(pairname,receptor,ligand,Pair.Source:Lig_cluster)
     })
     return(result)
   })
   
   #print table with lig-rec pairs
   output$pairs_res = DT::renderDataTable({
     input$pairby
     input$clust1
     input$clust2
     input$genelist1
     input$genelist2
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
       DT::datatable(finalres(),
                     extensions = c('Buttons','Scroller'),
                     options = list(dom = 'Bfrtip',
                                    searchHighlight = TRUE,
                                    pageLength = 10,
                                    lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                    scrollX = TRUE,
                                    buttons = c('copy', 'print')
                     ),rownames=FALSE,caption= "Result",selection = list(mode = 'single', selected =1),escape = F)
     })
   })
   
   #plot the bi-gene plot based on the row selected from lig-rec table
   output$bigeneplot2 <- renderPlot({
     withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
       table=finalres()
       s=input$pairs_res_rows_selected
       table=table[s, ,drop=FALSE]
       validate(need(nrow(table) > 0,"No Ligand-receptor pairs found"))
       bigene_genea=as.character(table$ligand)
       bigene_geneb=as.character(table$receptor)
       p=bigene_plot(fileload(),
                     c(bigene_genea,bigene_geneb),
                     limita=input$bigene_rangea,
                     limitb=input$bigene_rangeb,
                     marker_size = input$bigene_pointsize2,type=input$bigenedimr)
       p2 <- add_sub(p, paste(projectname(),"_Bigeneplot",sep=""), x = 0.87,vpadding = grid::unit(1, "lines"),size=11)
       ggdraw(p2)
     })
   })
   
   ######################################################################################################
   ######################################################################################################
   ################################# Network ############################################################
   ######################################################################################################
   ######################################################################################################
   #Generate drop-down to generate variables based on which you want to find pairs 
   output$pairbynet <- renderUI({
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
       scrna=fileload()
       metadata=as.data.frame(scrna@meta.data)
       metadata=metadata %>% dplyr::select(starts_with("var_"))
       options=colnames(metadata)
       selectInput("pairbynet","Select cell group ",options,selected=options[1])
     })
   })
   
   #Generate slider to filter ligand receptor pairs by frequency of occurence
   output$filternet <- renderUI({
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
       validate(need(input$lrnset != 0,"Make your selections and click the Set filters button"))
       
       if(input$lrnset == 0)
         return()
       isolate({
         scrna=fileload()
         result=scrna@misc$ligrecres
       #result=ligrec(fileload(),pair=input$pairbynet,prj=input$projects,input$perc_cells2,filetype=input$filetype)
       if(input$checksource2==T){result=result[result$Pair.Source %in% input$source2,]}
       if(input$checkevi2==T){result=result[result$Pair.Evidence %in% input$evidence2,]}
       if(input$checkgrp==T){result=result[result$Lig_cluster %in% input$checkgrp2,]
       result=result[result$Receptor_cluster %in% input$checkgrp2,]}
       validate(need(is.na(result)==F,"No pairs found. Change Filtering options"))
       edges=result %>% dplyr::select(Receptor_cluster,Lig_cluster)
       colnames(edges)=c("from","to")
       e2=as.data.frame(table(edges[,1:2]))
       min=min(e2$Freq)
       n=ifelse(min<4,4,min)
       max=max(e2$Freq)
       sliderInput("filternet", "Filter by number of interactions",
                   min = min, max = max, value = c(n,max),step=1)
     })
   })})
   
   #Check if the data is from mouse/human and use the approprite file to list the options for source
   output$source2 <- renderUI({
     file = read.csv("data/param.csv")
     if (input$filetype=="list"){
       org=as.character(file$organism[file$projects==input$projects])}
     else if(input$filetype == "upload"){
       scrna=fileload()
       validate(need(scrna@meta.data$org,"Organism not found in meta data. Cannot proceed"))
       org = unique(scrna@meta.data$org)
     }
     if(org=="mouse"){rl=read.csv("data/Mm_PairsLigRec.csv")}else if(org=="human"){rl=read.csv("data/Hs_PairsLigRec.csv")}
     options=as.character(unique(rl$Pair.Source))
     #checkboxGroupInput('source2', label='Select source(s)',choices=options,selected=options[1])
     selectInput('source2', label='Select source(s)',options,multiple=TRUE, selectize=TRUE,selected=options[2])
   })
   
   #Check if the data is from mouse/human and use the approprite file to list the options for evidence
   output$evidence2 <- renderUI({
     options=c("EXCLUDED not ligand","literature supported","putative","EXCLUDED not receptor","EXCLUDED")
     #checkboxGroupInput('evidence2',label='Select Evidence(s)',choices=options,selected=options[2])
     selectInput('evidence2',label='Select Evidence(s)',options,multiple=TRUE, selectize=TRUE,selected=options[2])
   })

   #Select groups to view in lig-rec pairs
   output$checkgrp <- renderUI({
     scrna=fileload()
     groups=unique(eval(parse(text=paste("scrna@meta.data$",input$pairbynet,sep=""))))
     selectInput('checkgrp2', 'Group options',groups, multiple=TRUE, selectize=TRUE)
   })
   
   #For selected project and grouping variable, generate all possible ligand receptor pairs and filter based on user input
   datasetInputnet = reactive({
     #result=NA
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
     validate(need(input$lrngo != 0,"Make your selections and click the Run button"))
     
     if(input$lrngo == 0)
       return()
     isolate({
       scrna=fileload()
       result=scrna@misc$ligrecres
     #result=ligrec(fileload(),pair=input$pairbynet,prj=projectname(),input$perc_cells2,filetype=input$filetype)
     #validate(need(is.na(result)==F,"Invalid Cell group. Pick a different option"))
     if(input$checksource2==T){result=result[result$Pair.Source %in% input$source2,]}
     if(input$checkevi2==T){result=result[result$Pair.Evidence %in% input$evidence2,]}
     if(input$checkgrp==T){result=result[result$Lig_cluster %in% input$checkgrp2,]
     result=result[result$Receptor_cluster %in% input$checkgrp2,]}
     edges=result %>% dplyr::select(Receptor_cluster,Lig_cluster)
     e2=as.data.frame(table(edges[,1:2]))
     e2=e2[e2$Freq>= input$filternet[1] & e2$Freq<= input$filternet[2],]
     e2$pair=paste(e2$Receptor_cluster,"_",e2$Lig_cluster,sep="")
     result$pair=paste(result$Receptor_cluster,"_",result$Lig_cluster,sep="")
     result=result[result$pair %in% e2$pair,]
     result=result %>% dplyr::select(pairname,receptor,ligand,Pair.Source:Lig_cluster)
     })
     return(result)
     })
   })
   
   #Render the same lig-rec pairs data table again to create network
   output$pairs_res2 = DT::renderDataTable({
     input$pairbynet
     input$filternet
     input$source2
     input$evidence2
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
       DT::datatable(datasetInputnet(),
                     extensions = c('Buttons','Scroller'),
                     options = list(dom = 'Bfrtip',
                                    searchHighlight = TRUE,
                                    pageLength = 10,
                                    lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                    scrollX = TRUE,
                                    buttons = c('copy', 'print')
                     ),rownames=FALSE,caption= "Ligand Receptor Pairs Result",selection = list(mode = 'single', selected =1),escape = F)
     })
   })
   
   
   #create lig-receptor network plot
   lrnetwork = reactive({
     withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
       result=datasetInputnet()
       rec <- result %>% distinct(Receptor_cluster) %>% rename(label = Receptor_cluster)
       lig <- result %>% distinct(Lig_cluster) %>% rename(label = Lig_cluster)
       nodes <- full_join(rec,lig, by = "label")
       nodes <- nodes %>% rowid_to_column("id")
       col=cpallette[1:nrow(nodes)]
       nodes$color=col
       perpair <- result %>% group_by(Receptor_cluster, Lig_cluster) %>% summarise(freq = n()) %>% ungroup()
       edges <- perpair %>%  left_join(nodes, by = c("Receptor_cluster" = "label")) %>% rename(from = id)
       
       edges <- edges %>% left_join(nodes, by = c("Lig_cluster" = "label")) %>% rename(to = id)
       edges <- dplyr::select(edges, from, to, freq)
       edges=left_join(edges,nodes,by=c("from"="id")) %>% dplyr::select(-label)
       edge.col=edges$color
       edge.lab=as.character(edges$freq)
       OldRange = (max(edges$freq) - min(edges$freq))  
       NewRange = 8-2  
       edges$width = (((edges$freq - min(edges$freq)) * NewRange) / OldRange) + 1.5
       width=edges$width
       #network <- network(edges, vertex.attr = nodes, matrix.type = "edgelist", ignore.eval = FALSE)
       
       routes_igraph <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)
       curves <-autocurve.edges2(routes_igraph)
       plot(routes_igraph, edge.arrow.size = 0.2,vertex.label.color="black",edge.label.color="black",vertex.color=col,edge.color=edge.col,edge.width=width,edge.arrow.width=9.5,edge.arrow.size=9.5,edge.curved=curves)
     })
   })
   
   
   #Render the ligand-receptor network plot
   output$lrnetwork = renderPlot({
        lrnetwork()
   })
   
   #Download network
   output$dwldnet <- downloadHandler(
       filename = function() {
         paste0("Ligand-Receptorpair_network.pdf")
       },
       content = function(file){
         pdf(file, width = 13, height = 8,useDingbats=FALSE)
         lrnetwork()
         dev.off()
       })
     
   ######################################################################################################
   ######################################################################################################
   ####################################### GO ANALYSIS ##################################################
   ######################################################################################################
   ######################################################################################################
   #Create dropdown for cluster names
   output$grp1 <- renderUI({
     result=datasetInputnet()
     rec <- result %>% distinct(Receptor_cluster)
     lig <- result %>% distinct(Lig_cluster)
     fluidRow(
     column(6,selectInput('grp1', label='Select Receptor Group for GO Analysis',rec)),
     column(6,selectInput('grp2', label='Select Ligand Group for GO Analysis',lig))
     )
   })
   
   #Get gene list and use gene list to get GO terms
   go_genelist <- reactive({
     result=datasetInputnet()
     genes= as.character(result$ligand[result$Receptor_cluster==input$grp1 & result$Lig_cluster==input$grp2])
     file = read.csv("data/param.csv")
     if (input$filetype=="list"){
       org=as.character(file$organism[file$projects==input$projects])}
     else if(input$filetype == "upload"){
       scrna=fileload()
       validate(need(scrna@meta.data$org,"Organism not found in meta data. Cannot proceed"))
       org = unique(scrna@meta.data$org)
     }
     if(org=="human"){
       dataset="hsapiens_gene_ensembl"
     }
     else if(org=="Rat"){
       dataset="rnorvegicus_gene_ensembl"
     }
     else{
       dataset="mmusculus_gene_ensembl"
     }
     ensembl = useEnsembl(biomart="ensembl", dataset=dataset)
     GO <- getBM(attributes=c('go_id','name_1006','definition_1006'), filters ='external_gene_name', values =genes, mart = ensembl)
     return(GO)
   })
   
   #Render the GO term table for ligand-receptor pairs 
   output$gotable = DT::renderDataTable({
     input$grp1
     input$grp2
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
       DT::datatable(go_genelist(),
                     extensions = c('Buttons','Scroller'),
                     options = list(dom = 'Bfrtip',
                                    searchHighlight = TRUE,
                                    pageLength = 10,
                                    lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                    scrollX = TRUE,
                                    buttons = c('copy', 'print')
                     ),rownames=FALSE,caption= "GO Terms",selection = list(mode = 'single', selected =1),escape = F)
     })
   })
   ######################################################################################################
   ######################################################################################################
   ################################# Ligand Receptor Heatmap ############################################
   ######################################################################################################
   ######################################################################################################
   #Generate drop-down to generate variables based on which you want to find pairs 
   output$pairbyheatnet <- renderUI({
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
       scrna=fileload()
       metadata=as.data.frame(scrna@meta.data)
       metadata=metadata %>% dplyr::select(starts_with("var_"))
       options=colnames(metadata)
       selectInput("pairbyheatnet","Select cell group ",options,selected=options[1])
     })
   })
   
   #Check if the data is from mouse/human and use the approprite file to list the options for source
   output$source3 <- renderUI({
     file = read.csv("data/param.csv")
     if (input$filetype=="list"){
       org=as.character(file$organism[file$projects==input$projects])}
     else if(input$filetype == "upload"){
       scrna=fileload()
       validate(need(scrna@meta.data$org,"Organism not found in meta data. Cannot proceed"))
       org = unique(scrna@meta.data$org)
     }
     if(org=="mouse"){rl=read.csv("data/Mm_PairsLigRec.csv")}else if(org=="human"){rl=read.csv("data/Hs_PairsLigRec.csv")}
     options=as.character(unique(rl$Pair.Source))
     checkboxGroupInput('source3', label='Select source(s)',choices=options,selected=options[2])
   })
   
   #Check if the data is from mouse/human and use the approprite file to list the options for evidence
   output$evidence3 <- renderUI({
     # file = read.csv("data/param.csv")
     # org=as.character(file$organism[file$projects==input$projects])
     # if(org=="mouse"){rl=read.csv("data/Mm_PairsLigRec.csv")}else if(org=="human"){rl=read.csv("data/Hs_PairsLigRec.csv")}
     # options=as.character(unique(rl$Pair.Evidence))
     options=c("EXCLUDED not ligand","literature supported","putative","EXCLUDED not receptor","EXCLUDED")
     checkboxGroupInput('evidence3',label='Select Evidence(s)',choices=options,selected=options[2])
   })
   
   
   #Generate lig-receptor pairs table
   ligrecheat = reactive({
     validate(need(input$lrhgo != 0,"Make your selections and click the Run button"))
     
     if(input$lrhgo == 0)
       return()
     isolate({
       scrna=fileload()
       result=scrna@misc$ligrecres
     #result=ligrec(fileload(),pair=input$pairbyheatnet,prj=projectname(),input$perc_cells3,filetype=input$filetype)
     result=result %>% dplyr::select(pairname,receptor,ligand,Pair.Source:Lig_cluster)
     if(input$checksourceheat==T){result=result[result$Pair.Source %in% input$source3,]}
     if(input$checkeviheat==T){result=result[result$Pair.Evidence %in% input$evidence3,]}
     })
     return(result)
   })
   
   #Render the same lig-rec pairs data table again to create network
   output$pairs_res3 = DT::renderDataTable({
     input$pairbyheatnet
     input$source3
     input$evidence3
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
       DT::datatable(ligrecheat(),
                     extensions = c('Buttons','Scroller'),
                     options = list(dom = 'Bfrtip',
                                    searchHighlight = TRUE,
                                    pageLength = 10,
                                    lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                    scrollX = TRUE,
                                    buttons = c('copy', 'print')
                     ),rownames=FALSE,caption= "Ligand Receptor Pairs Result",selection = list(mode = 'single', selected =1),escape = F)
     })
   })
   
  #Get ligand-receptor pairs and compute frequency and use it to generate a heatmap
   netheatmap = reactive({
     dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
     result=ligrecheat()
     tab=table(result[,6:7])
     tab=as.data.frame(tab)
     tab = tab %>% spread(Receptor_cluster,Freq)
     rownames(tab)=tab$Lig_cluster
     tab= tab %>% dplyr::select(-Lig_cluster)
     tab=Filter(var, tab)
     if(input$clusterby=="both"){
       row=TRUE
       column=TRUE
     }else if(input$clusterby=="row"){
       row=TRUE
       column=FALSE
     }else if(input$clusterby=="column"){
       row=FALSE
       column=TRUE
     }else if(input$clusterby=="none"){
       row=NA
       column=NA
     }
     if(input$checkbox==TRUE){
       aheatmap(as.matrix(tab),distfun=dist2,Rowv=row,Colv=column,col = colorRampPalette(brewer.pal(n = 9,input$hmpcolnet))(30),main= "Receptor genes (x) vs Ligand genes (y)")}
     else{aheatmap(as.matrix(tab),distfun=dist2,Rowv=row,Colv=column,col = colorRampPalette(rev(brewer.pal(n = 9,input$hmpcolnet)))(30),main= "Receptor genes (x) vs Ligand genes (y)")}
   })
   #Render heatmap
   output$netheatmap = renderPlot({
     input$hmpcolnet
     input$clusterby
     input$checkbox
     netheatmap()
   })
   
   #Download function for the heatmap
   output$downloadlrheatmap <- downloadHandler(
     filename = function() {
       paste0("Ligand-Receptorpair_Heatmap.pdf")
     },
     content = function(file){
       pdf(file, width = 13, height = 8,useDingbats=FALSE)
       netheatmap()
       dev.off()
     })
   ######################################################################################################
   ######################################################################################################
   ################################# TROUBLESHOOT #######################################################
   ######################################################################################################
   ######################################################################################################
   #list devices in use
   output$device <- renderText({ 
     dev.list()
   })
   
   #Reset graphics device
   observeEvent(input$devoff, {
     graphics.off()
   })

}#end of server