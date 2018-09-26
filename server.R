library(shiny)
library(shinyBS)
library(RColorBrewer)
library(Biobase)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)
library(shinyjs)
library(htmlwidgets)
library(DT)
library(shinyRGL)
library(rgl)
library(rglwidget)
library(Seurat)
library(cowplot)
library(data.table)
#library(ggnetwork)
#library(igraph)
library(visNetwork)
source("functions.R")

#Specify color palette for the tSNE and UMAP plots
cpallette=c("#64B2CE", "#DA5724", "#74D944", "#CE50CA", "#C0717C", "#CBD588", "#5F7FC7",
            "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
            "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D",
            "#8A7C64", "#599861","#000099","#FFCC66","#99CC33","#CC99CC","#666666")

#Specify user-ids and passwords
my_username <- c("Sealelab","Morriseylab","Jainlab","allusers")
my_password <- c("pseale@999#","emorrisey$123","rjain@2018","allusers@1")

server <- function(input, output,session) {
  
  ###################################################
  ###################################################
  ################### AUTHENTICATION  ##############
  ###################################################
  ###################################################
  #  values <- reactiveValues(authenticated = FALSE)
  # 
  #  # Return the UI for a modal dialog with data selection input. If 'failed'
  #  # is TRUE, then display a message that the previous value was invalid.
  #  dataModal <- function(failed = FALSE) {
  #    modalDialog(
  #      textInput("username", "Username:"),
  #      passwordInput("password", "Password:"),
  #      footer = tagList(
  #        # modalButton("Cancel"),
  #        actionButton("ok", "OK")
  #      )
  #    )
  #  }
  # 
  #  # Show modal when button is clicked.
  #  # This `observe` is suspended only whith right user credential
  # 
  #  obs1 <- observe({
  #    showModal(dataModal())
  #  })
  # 
  # # When OK button is pressed, attempt to authenticate. If successful,
  # # remove the modal.
  # obs2 <- observe({
  #   req(input$ok)
  #   isolate({
  #     Username <- input$username
  #     Password <- input$password
  #   })
  #   Id.username <- which(my_username == Username)
  #   Id.password <- which(my_password == Password)
  #   if (length(Id.username) > 0 & length(Id.password) > 0) {
  #     if (Id.username == Id.password) {
  #       Logged <<- TRUE
  #       values$authenticated <- TRUE
  #       obs1$suspend()
  #       removeModal()
  # 
  #     } else {
  #       values$authenticated <- FALSE
  #     }
  #   }
  # })

  ###################################################
  ###################################################
  ####### Display username in notification bar  ####
  ###################################################
  ###################################################
  # output$userloggedin = renderMenu({
  #   msg= paste("Logged in as ",input$username,sep="")
  #   prj= paste("Data: ",input$projects,sep="")
  #   dropdownMenu(type = "notifications", badgeStatus = "info",
  #                notificationItem(icon = icon("user"), status = "info",msg),
  #                notificationItem(icon = icon("book"), status = "info",prj))
  # })
  ###################################################
  ###################################################
  ####### Display project list and load data  #######
  ###################################################
  ###################################################
  #Read the parameter file
  readexcel = reactive({
     #user=input$username
     file = read.csv("data/param.csv")
    # if(user=="allusers"){
    #   file=file
    # }else{
    #   file=file[file$user==user,]
    # }
  })
  
  #Get Project list and populate drop-down
  output$projects = renderUI({
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
    if(user=="allusers"){
      file=file
    }else{
      file=file[file$user==user,] %>% dplyr::select(-user)
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
    # validate(need(input$load != 0,"Load dataset"))
    # 
    # if(input$load == 0)
    #   return()
    # isolate({
    inFile = paste('data/',as.character(input$projects),'.RData',sep = '')
    load(inFile)
    #loaddata=scrna
    return(scrna)
    # scrna <<- scrna
    # })
  })

  ###################################################
  ###################################################
  ################## Project Summary  ###############
  ###################################################
  ###################################################
  
  #Get all information from the scrna object (input file) and generate some basic project summary for the summary
  prjsumm <- reactive({
    #user=input$username
    user="allusers"
   prj= read.csv("data/param.csv")
   if(user=="allusers"){
     prj=prj
   }else{
     prj=prj[prj$user==user,] 
   }
   prj=prj[prj$projects==input$projects,]
   pname=prj$projects
   pdesc=prj$desc
   porg=prj$organism
   scrna=fileload()
   tcells=dim(scrna@data)[2]
   tgenes=dim(scrna@data)[1]
   if(is.null(scrna@dr$pca)){
   maxdim=dim(scrna@dr$cca.aligned@cell.embeddings)[2]
   }else{maxdim=length(scrna@dr$pca@sdev)}
   c.cnt=as.data.frame(table(scrna@ident))
   df=as.data.frame(c(as.character(pname),as.character(pdesc),as.character(porg),tcells,tgenes,maxdim,"","",c.cnt$Freq))
   rownames(df)=c("Project name","Project Description","Organism","Total nummber of cells","Total number of genes","Dimension","","Cluster-wise number of genes",as.character(c.cnt$Var1))
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
  ###################### VARIABLE GENES TAB ############################################################
  ######################################################################################################
  ######################################################################################################
  
  #Load the scrna input RData file and extract the variable genes and its stats
  vargenes= reactive({
    scrna=fileload()
    var=as.data.frame(scrna@var.genes)
    colnames(var)="Gene"
    stat=scrna@hvg.info
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
    maxdim=length(scrna@dr$pca@sdev)
    var=1:maxdim
    selectInput("ndim","Choose number of dimensions",var,selected = 1)
  })
  
  #Plot the PCA/Viz plot for the number of dimensions and number of genes chosen
  vizplot= reactive({
    scrna=fileload()
    dim=input$ndim
    par(mar=c(4,5,3,3))
    g1=VizPCA(object = scrna, pcs.use = dim:dim,nCol=1,font.size = 1,num.genes = input$ngenes)
    return(g1) 
  })
  
  #Render the vizplot
  output$vizplot = renderPlot({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      vizplot()
    })
    dev.off()
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
    clusts=levels(scrna@ident)
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
    clusts=levels(scrna@ident)
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
    dimr=names(scrna@dr)
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
    selectInput("umapa","Dimensionality Reduction",dimr,selected = "tsne")})
  })
  
  #Dimensionality reduction options for right plot
  output$umapb = renderUI({
    scrna=fileload()
    dimr=names(scrna@dr)
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
    selectInput("umapb","Dimensionality Reduction",dimr,selected = "tsne")})
  })
  
  #Based on all use input, generate plots using the right category and dimensionality reduction methods
  comptsne2 = reactive({
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data)
    met= sapply(metadata,is.numeric)
    #metadata=metadata %>% select(starts_with("var"))
    tsnea=input$tsnea2
    tsneb=input$tsneb2
    feature=names(met[met==TRUE])
    #feature=c("nGene","nUMI","percent.mito","S.Score","G2M.Score","var.ratio.pca")
    tsne=names(met[met==FALSE])
    
    if(input$categorya2 =="clust" & input$subsa==F){
      plot1=DimPlot(object = scrna,reduction.use=input$umapa,group.by = "ident",no.legend = FALSE,do.label = TRUE,vector.friendly = T, do.return=T, pt.size = input$pointa2,label.size = 7, cols.use=cpallette)
    }else if(input$categorya2 =="clust" & input$subsa==TRUE){
      cells=names(scrna@ident[scrna@ident==input$selclust])
      plot1=DimPlot(object = scrna,reduction.use=input$umapa,cells.highlight=cells,group.by = "ident",vector.friendly = T,no.legend = FALSE,do.label = F, do.return=T, pt.size = input$pointa2, cols.use=cpallette)
    }else if(input$categorya2=="geneexp"){
      plot1=FeaturePlot(object = scrna,reduction.use=input$umapa,vector.friendly = T, features.plot = input$gene1a, cols.use = c("grey", "blue"),do.return=T,pt.size = input$pointa2,no.legend = FALSE)
      plot1=eval(parse(text=paste("plot1$",input$gene1a,sep="")))
    }else if(input$categorya2 =="var" & input$tsnea2 %in% tsne & input$subsa==FALSE){
      plot1=DimPlot(object = scrna,reduction.use=input$umapa,group.by = tsnea,no.legend = FALSE,do.label = TRUE,vector.friendly = T, do.return=T,pt.size = input$pointa2,label.size = 7, cols.use=cpallette)
    }else if(input$categorya2 =="var" & input$tsnea2 %in% tsne & input$subsa==TRUE){
      t=paste("rownames(scrna@meta.data[scrna@meta.data$",input$tsnea2,"==\"",input$selclust2,"\",])",sep="")
      cells=eval(parse(text=t))
      plot1=DimPlot(object = scrna,reduction.use=input$umapa,group.by = tsnea,cells.highlight=cells,vector.friendly = T,no.legend = FALSE,do.label =F, do.return=T,pt.size = input$pointa2, cols.use=cpallette)
    }else if(input$categorya2 =="var" & input$tsnea2 %in% feature & input$subsa==FALSE){
      plot1=FeaturePlot(object = scrna,reduction.use=input$umapa,vector.friendly = T, features.plot = tsnea, cols.use = c("grey", "blue"),do.return=T,pt.size = input$pointa2,no.legend = FALSE)
      plot1=eval(parse(text=paste("plot1$",tsnea,sep="")))
    }else if(input$categorya2 =="var" & input$tsnea2 %in% feature & input$subsa==TRUE){
      t=paste('rownames(scrna@meta.data[scrna@meta.data$',input$tsnea2, '>',input$tsnea2lim[1], ' & metadata$',input$tsnea2, '<', input$tsnea2lim[2],',])',sep="")
      cells=eval(parse(text=t))
      plot1=FeaturePlot(object = scrna,reduction.use=input$umapa, features.plot = tsnea,cells.use = cells,vector.friendly = T, cols.use = c("grey", "blue"),do.return=T,pt.size = input$pointa2,no.legend = FALSE)
      plot1=eval(parse(text=paste("plot1$",tsnea,sep="")))
    }
    
    if(input$categoryb2 =="clust" & input$subsb==F){
      plot2=DimPlot(object = scrna,reduction.use=input$umapb,group.by = "ident",no.legend = FALSE,vector.friendly = T,do.label = TRUE, do.return=T,pt.size = input$pointa2,label.size = 7, cols.use=cpallette)
    }else if(input$categoryb2 =="clust" & input$subsb==TRUE){
      cells=names(scrna@ident[scrna@ident==input$selclustb])
      plot2=DimPlot(object = scrna,reduction.use=input$umapb,cells.highlight=cells,group.by = "ident",vector.friendly = T,no.legend = FALSE,do.label = F, do.return=T, pt.size = input$pointa2, cols.use=cpallette)
    }else if(input$categoryb2=="geneexp"){
      plot2=FeaturePlot(object = scrna,reduction.use=input$umapb,vector.friendly = T, features.plot = input$gene2a, cols.use = c("grey", "blue"),do.return=T,pt.size = input$pointa2,no.legend = FALSE)
      plot2=eval(parse(text=paste("plot2$",input$gene2a,sep="")))
    }else if(input$categoryb2 =="var" & input$tsneb2 %in% tsne & input$subsb==F){
      plot2=DimPlot(object = scrna,reduction.use=input$umapb,group.by = tsneb,no.legend = FALSE,vector.friendly = T,do.label = TRUE, do.return=T,pt.size = input$pointa2,label.size = 7, cols.use=cpallette)
    }else if(input$categoryb2 =="var" & input$tsneb2 %in% tsne & input$subsb==TRUE){
      t=paste("rownames(scrna@meta.data[scrna@meta.data$",input$tsneb2,"==\"",input$selclustb2,"\",])",sep="")
      cells=eval(parse(text=t))
      plot2=DimPlot(object = scrna,reduction.use=input$umapb,group.by = tsneb,cells.highlight=cells,vector.friendly = T,no.legend = FALSE,do.label = F, do.return=T,pt.size = input$pointa2, cols.use=cpallette)
    }else if(input$categoryb2 =="var" & input$tsneb2 %in% feature & input$subsb==F){
      plot2=FeaturePlot(object = scrna,reduction.use=input$umapb, features.plot = tsneb,vector.friendly = T, cols.use = c("grey", "blue"),do.return=T,pt.size = input$pointa2,no.legend = FALSE)
      plot2=eval(parse(text=paste("plot2$",tsneb,sep="")))
    }else if(input$categoryb2 =="var" & input$tsneb2 %in% feature & input$subsb==TRUE){
      t=paste('rownames(scrna@meta.data[scrna@meta.data$',input$tsneb2, '>',input$tsneb2lim[1], ' & metadata$',input$tsneb2, '<', input$tsneb2lim[2],',])',sep="")
      cells=eval(parse(text=t))
      plot2=FeaturePlot(object = scrna,reduction.use=input$umapb, features.plot = tsneb,vector.friendly = T,cells.use = cells, cols.use = c("grey", "blue"),do.return=T,pt.size = input$pointa2,no.legend = FALSE)
      plot2=eval(parse(text=paste("plot2$",tsneb,sep="")))
    }
    
    p=plot_grid(plot1,plot2)
    p2 <- add_sub(p, paste(input$projects,"_CompareTsne",sep=""), x = 0.87,vpadding = grid::unit(1, "lines"),size=11)
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
      paste0(input$projects,"_CompareTsne.pdf",sep="")
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
    metadata=metadata %>% select(starts_with("var_"))
    var=c(colnames(metadata))
    selectInput("setcategory","Choose category",var,"pick one")
  })
  
  #Generate drop down menu for Dimensionality reduction options for both plots
  output$umapint = renderUI({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
    scrna=fileload()
    dimr=names(scrna@dr)
    selectInput("umapint","Dimensionality Reduction",dimr,selected = "tsne")})
  })
  
  #Based on input options, generate left interative plot 
  # intertsne = reactive({
  #   scrna=fileload()
  #   plot1=DimPlot(object = scrna,reduction.use=input$umapint,group.by = input$setcategory,no.legend = FALSE,do.label = TRUE, do.return=T,pt.size = input$umap_pointsize,label.size = 5,vector.friendly = T, cols.use=cpallette)
  #   plot=ggplotly(plot1)
  #   return(plot)
  # })
  
  
  output$intertsne = renderPlotly({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      pdf(NULL)
      scrna=fileload()
      plot1=DimPlot(object = scrna,reduction.use=input$umapint,group.by = input$setcategory,no.legend = FALSE,do.label = TRUE, do.return=T,pt.size = input$umap_pointsize,label.size = 5,vector.friendly = T, cols.use=cpallette)
      plot=ggplotly(plot1)
      dev.off()
      return(plot)
    })
  })
  
  #Based on input options, generate right interative plot 
  # intergene = reactive({
  #   pdf(NULL)
  #   scrna=fileload()
  #   metadata=as.data.frame(scrna@meta.data)
  #   met= sapply(metadata,is.numeric)
  #   #metadata=metadata %>% select(starts_with("var"))
  #   tsnea=input$intervar
  #   feature=names(met[met==TRUE])
  #   #feature=c("nGene","nUMI","percent.mito","S.Score","G2M.Score","var.ratio.pca")
  #   tsne=names(met[met==FALSE])
  #   
  #   if(input$intercat=="geneexp"){
  #     plot1=FeaturePlot(object = scrna,reduction.use=input$umapint, features.plot = input$geneinter,vector.friendly = T, cols.use = c("grey", "blue"),do.return=T,pt.size = input$umap_pointsize,no.legend = FALSE)
  #     plot1=eval(parse(text=paste("plot1$",input$geneinter,sep="")))
  #   }else if(input$intercat =="var" & tsnea %in% tsne){
  #     plot1=DimPlot(object = scrna,reduction.use=input$umapint,group.by = tsnea,no.legend = FALSE,do.label = TRUE, do.return=T,pt.size = input$umap_pointsize,label.size = 7, cols.use=cpallette,vector.friendly = T)
  #   }else if(input$intercat =="var" & tsnea %in% feature){
  #     plot1=FeaturePlot(object = scrna,reduction.use=input$umapint, features.plot = tsnea,vector.friendly = T, cols.use = c("grey", "blue"),do.return=T,pt.size = input$umap_pointsize,no.legend = FALSE)
  #     plot1=eval(parse(text=paste("plot1$",tsnea,sep="")))
  #   }
  #   plot=ggplotly(plot1)
  #   dev.off()
  #   return(plot)
  # })
  
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
        plot1=FeaturePlot(object = scrna,reduction.use=input$umapint, features.plot = input$geneinter,vector.friendly = T, cols.use = c("grey", "blue"),do.return=T,pt.size = input$umap_pointsize,no.legend = FALSE)
        plot1=eval(parse(text=paste("plot1$",input$geneinter,sep="")))
      }else if(input$intercat =="var" & tsnea %in% tsne){
        plot1=DimPlot(object = scrna,reduction.use=input$umapint,group.by = tsnea,no.legend = FALSE,do.label = TRUE, do.return=T,pt.size = input$umap_pointsize,label.size = 7, cols.use=cpallette,vector.friendly = T)
      }else if(input$intercat =="var" & tsnea %in% feature){
        plot1=FeaturePlot(object = scrna,reduction.use=input$umapint, features.plot = tsnea,vector.friendly = T, cols.use = c("grey", "blue"),do.return=T,pt.size = input$umap_pointsize,no.legend = FALSE)
        plot1=eval(parse(text=paste("plot1$",tsnea,sep="")))
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
    dimr=names(scrna@dr)
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      selectInput("bigenedim","Dimensionality Reduction",dimr,selected = "tsne")})
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
      p2 <- add_sub(p, paste(input$projects,"_Bigeneplot",sep=""), x = 0.87,vpadding = grid::unit(1, "lines"),size=11)
      ggdraw(p2)
    })
  })
  
  
  #Download bi-gene plot
  output$downloadbiplot <- downloadHandler(
    filename = function() {
      paste0(input$projects,"_",input$bigene_genea,"_",input$bigene_geneb,"_Bigene.pdf",sep="")
    },
    content = function(file){
      pdf(file,width=9,height = 9,useDingbats=FALSE)
      plot(bigene_plot(fileload(),
                       c(input$bigene_genea,input$bigene_geneb),
                       limita=input$bigene_rangea,
                       limitb=input$bigene_rangeb,
                       marker_size = input$bigene_pointsize,type=input$bigenedimr))
      dev.off()
    })
  
  ####################################################
  ###################################################
  ########## Setup Control Panel for DEG ############
  ###################################################
  ###################################################
  #Generate drop down for dimensionality reduction in DEG tab
  output$umapdeg = renderUI({
    scrna=fileload()
    dimr=names(scrna@dr)
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      selectInput("umapdeg","Dimensionality Reduction",dimr,selected = "tsne")})
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
    options=names(scrna@misc)
    selectInput("identdef", "First cluster/variable of comparison",options)
  })
  
  #Generate drop down menu for the categories starting with "var" that can be set as the new ident/variable of comparison
  output$setidentlist = renderUI({
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data)
    metadata=metadata %>% select(starts_with("var_"))
    var=c(colnames(metadata))
    selectInput("setidentlist","Choose category to compare",var,"pick one")
    
  })
  
  #Generate drop down menu for the variables in the new ident against which the rest will be compared to find markers
  output$identa = renderUI({
    scrna=fileload()
    if(input$setident==T){
      scrna <- SetAllIdent(object = scrna, id = input$setidentlist)
      options=unique(scrna@ident)
    }else{
      options=levels(scrna@ident)
    }
    selectInput("identa", "First Cell group to compare",options)
  })
  
  #Generate checkboxes for the variables in the new ident. Choose all those that you want to compare with the first 
  output$identb = renderUI({
    scrna=fileload()
    if(input$setident==T){
      scrna <- SetAllIdent(object = scrna, id = input$setidentlist)
      options=unique(scrna@ident)
    }else{
      options=levels(scrna@ident)
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
          scrna <- SetAllIdent(object = scrna, id = input$setidentlist) #set ident
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
        markers=eval(parse(text=paste("scrna@misc$`",input$identdef,"`",sep="")))
        geneid=rownames(markers)
        url= paste("http://www.genecards.org/cgi-bin/carddisp.pl?gene=",geneid,sep = "")
        markers$Link=paste0("<a href='",url,"'target='_blank'>",rownames(markers),"</a>")
      }
    })
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
    filename = function() { paste(input$projects, '.csv', sep='') },
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
    scrna@meta.data$var_cluster=as.numeric(as.character(scrna@meta.data$var_cluster))
    tsnea=input$tsnea
    tsneb=input$tsneb
    feature=names(met[met==TRUE])
    tsne=names(met[met==FALSE])
    
    if(input$tsnea =="Cell.group"){
      plot1=DimPlot(object = scrna,reduction.use=input$umapdeg,group.by = "ident",no.legend = FALSE,do.label = TRUE, do.return=T, pt.size = input$pointa,label.size = 7,cols.use=cpallette,vector.friendly=TRUE) + theme(legend.position="bottom")
    }else if(input$tsnea %in% tsne){
      plot1=DimPlot(object = scrna,reduction.use=input$umapdeg,group.by = tsnea,no.legend = FALSE,do.label = TRUE,vector.friendly=TRUE, do.return=T,pt.size = input$pointa,label.size = 7,cols.use=cpallette) + theme(legend.position="bottom")
    }else if(input$tsnea %in% feature){
      plot1=FeaturePlot(object = scrna, features.plot = tsnea, cols.use = c("grey", "blue"),vector.friendly = T,reduction.use = input$umapdeg,do.return=T,pt.size = input$pointa)
      plot1=eval(parse(text=paste("plot1$`",tsnea,"`",sep="")))
    }
   
    markers=markergenes()
      s=input$markergenes_rows_selected # get  index of selected row from table
      markers=markers[s, ,drop=FALSE]
      plot2=FeaturePlot(object = scrna, features.plot = rownames(markers), cols.use = c("grey","blue"),reduction.use = input$umapdeg,vector.friendly = T,
                        no.legend = FALSE,pt.size = input$pointa,do.return = T)
      plot2=eval(parse(text=paste("plot2$`",rownames(markers),"`",sep="")))
      if(input$checkviolin ==T){
      plot3=VlnPlot(object = scrna, features.plot = rownames(markers),group.by = input$setidentlist,do.return = T,x.lab.rot=TRUE,point.size.use=0,cols.use=cpallette)
      }else{plot3=VlnPlot(object = scrna, features.plot = rownames(markers),group.by = input$setidentlist,do.return = T,x.lab.rot=TRUE,cols.use=cpallette)}
      plot4=RidgePlot(object = scrna, features.plot = rownames(markers),group.by = input$setidentlist,do.return = T,x.lab.rot=TRUE,cols.use=cpallette)
      
    
      row1=plot_grid(plot1,plot2,align = 'h', rel_heights = c(1, 1),axis="lr", nrow=1)
      row2=plot_grid(plot3,plot4,align = 'h', rel_heights = c(1, 1),axis="lr", nrow=1)
    p=plot_grid(row1,row2,align = 'v', rel_heights = c(1.7, 1),axis="tb",ncol=1)
    p2 <- add_sub(p, paste(input$projects,"_Differential_Exp",sep=""), x = 0.87,vpadding = grid::unit(1, "lines"),size=11)
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
      paste0("Plot.pdf")
    },
    content = function(file){
      pdf(file, width = 12, height = 11,useDingbats=FALSE)
      plot(comptsne())
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
    dimr=names(scrna@dr)
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      selectInput("umapge","Dimensionality Reduction",dimr,selected = "tsne")})
  })
  
  #For any gene entered, generate a feature plot, violin plot and a ridge plot
  geplots = reactive({
    scrna=fileload()
    validate(need(input$geneid,"Enter the gene symbol"))
    plot2=FeaturePlot(object = scrna, features.plot = input$geneid, cols.use = c("grey","blue"),reduction.use = input$umapge,vector.friendly = T,
                      no.legend = FALSE,pt.size = input$genenid_pointsize,do.return = T)
    plot2=eval(parse(text=paste("plot2$`",input$geneid,"`",sep="")))
    plot3=VlnPlot(object = scrna, features.plot = input$geneid,group.by = "ident",do.return = T,x.lab.rot=TRUE,cols.use=cpallette)
    plot4=RidgePlot(object = scrna, features.plot = input$geneid,group.by = "ident",do.return = T,x.lab.rot=TRUE,cols.use=cpallette)
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
    dimr=names(scrna@dr)
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      selectInput("umapclust","Dimensionality Reduction",dimr,selected = "tsne")})
  })
  
  #Generate drop down for categories starting with "var" in Cluster-wise gene expression
  output$setvar = renderUI({
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data)
    metadata=metadata %>% select(starts_with("var_"))
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
    scrna <- SetAllIdent(object = scrna, id = input$setvar)
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
    plot1=DimPlot(object = scrna,reduction.use=input$umapclust,group.by = input$setvar,no.legend = FALSE,do.label = TRUE,vector.friendly = T, do.return=T, pt.size = input$pointclust,label.size = 7, cols.use=cpallette)
    plot2=FeaturePlot(object = scrna,reduction.use=input$umapclust, features.plot = gene, cols.use = c("grey", "blue"),vector.friendly = T,do.return=T,pt.size = input$pointclust,no.legend = FALSE)
    plot2=eval(parse(text=paste("plot2$`",gene,"`",sep="")))
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
    metadata=metadata %>% select(starts_with("var_"))
    var=c(colnames(metadata))
    selectInput("setdotvar","Choose category",var,"pick one")
  })
  
  #read in user-defined list of genes and generate the dot plot
  dotplot= reactive({
    scrna=fileload()
    validate(
      need(input$genelistfile, "Please Upload Genelist")
    )
    file=input$genelistfile
    df=fread(file$datapath,header = FALSE) #get complete gene list as string
    genes=as.vector(df$V1)
    g1=DotPlot(object = scrna, genes.plot = genes, plot.legend = TRUE,group.by=input$setdotvar,do.return=TRUE) 
    return(g1) 
  })
  
  #render the dot plot
  output$dotplot = renderPlot({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      dotplot()
    })
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
       metadata=metadata %>% select(starts_with("var_"))
       var=c("ident",colnames(metadata))
       selectInput("hmpgrp","Select a Variable",var,"pick one")
     })
   })
   
   #generate the heatmap
   heatmap <- reactive({
     withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
       scrna=fileload()
       markers=markergenes()
       if(input$hmpcol=="PuYl"){
         lowcol="darkmagenta"
         midcol="black"
         highcol="yellow"
       }else if(input$hmpcol=="BuGn"){
         lowcol="yellow"
           midcol="green"
           highcol="blue"
       }else if(input$hmpcol=="RdYl"){
         lowcol="yellow"
           midcol="red"
           highcol="black"
       }else if(input$hmpcol=="RdBu"){
         lowcol="red"
         midcol="white"
         highcol="blue"}
       p=DoHeatmap(object = scrna, genes.use = rownames(markers)[1:input$heatmapgenes],group.by = input$hmpgrp, draw.line= T,
                 group.label.rot= T, col.low=lowcol, col.mid =midcol ,col.high = highcol,slim.col.label=TRUE)
       p2 <- add_sub(p, paste(input$projects,"_Heatmap",sep=""), x = 0.87,vpadding = grid::unit(1, "lines"),size=11)
       ggdraw(p2)
     })
   })
   
   #Render the heatmap
   output$heatmap <- renderPlot({
     withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
       heatmap()
     })
     dev.off()
   })
   
   #Download function for the heatmap
   output$downloadheatmap <- downloadHandler(
     filename = function() {
       paste0("Heatmap.pdf")
     },
     content = function(file){
       pdf(file, width = 13, height = 8,useDingbats=FALSE)
       plot(heatmap())
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
     dimr=names(scrna@dr)
     withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
       selectInput("bigenedimr","Dimensionality Reduction",dimr,selected = "tsne")})
   })
   
   #Generate drop-down to generate variables based on which you want to find pairs 
   #This is important for choosing the 'select cluster' option
   output$pairby <- renderUI({
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
       scrna=fileload()
       metadata=as.data.frame(scrna@meta.data)
       metadata=metadata %>% select(starts_with("var_"))
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
     if(input$pairby=="ident"){options=levels(scrna@ident)}
     selectInput("clust1","Pick cellgroup for Receptor",options,selected=options[1])
     })
   })
   
   #If user chooses to select cluster and all genes, generate drop down menus to pick ligand cluster
   output$clust2 <- renderUI({
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
     scrna=fileload()
     t=paste("scrna@meta.data$",input$pairby,sep="")
     options=unique(eval(parse(text=t)))
     if(input$pairby=="ident"){options=levels(scrna@ident)}
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
   
   #If user chooses to select select both cluster and genes, generate drop down menus to pick ligand cluster
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
     org=as.character(file$organism[file$projects==input$projects])
     if(org=="mouse"){rl=read.csv("data/Mm_PairsLigRec.csv")}else if(org=="human"){rl=read.csv("data/Hs_PairsLigRec.csv")}
     options=as.character(unique(rl$Pair.Source))
     checkboxGroupInput('source', label='Select source(s)',choices=options,selected=options[1])
   })
   
   #Check if the data is from mouse/human and use the approprite file to list the options for evidence
   output$evidence <- renderUI({
     file = read.csv("data/param.csv")
     org=as.character(file$organism[file$projects==input$projects])
     if(org=="mouse"){rl=read.csv("data/Mm_PairsLigRec.csv")}else if(org=="human"){rl=read.csv("data/Hs_PairsLigRec.csv")}
     options=as.character(unique(rl$Pair.Evidence))
     checkboxGroupInput('evidence',label='Select Evidence(s)',choices=options,selected=options[1])
   })
   
   
   ###################################################
   ###################################################
   #### Load lig-receptor list and display results ###
   ###################################################
   ###################################################
   
   #For selected project and grouping variable, generate all possible ligand receptor pairs
   datasetInput = reactive({
     results=ligrec(fileload(),pair=input$pairby,prj=input$projects)
   })
   
   #Subselect lig-rec pairs based on user input
   finalres= reactive({
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
       p2 <- add_sub(p, paste(input$projects,"_Bigeneplot",sep=""), x = 0.87,vpadding = grid::unit(1, "lines"),size=11)
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
       metadata=metadata %>% select(starts_with("var_"))
       options=colnames(metadata)
       selectInput("pairbynet","Select cell group ",options,selected=options[1])
     })
   })
   
   #Generate slider to filter ligand receptor pairs by frequency of occurence
   output$filternet <- renderUI({
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
       result=ligrec(fileload(),pair=input$pairbynet,prj=input$projects)
       edges=result %>% dplyr::select(Receptor_cluster,Lig_cluster)
       colnames(edges)=c("from","to")
       e2=as.data.frame(table(edges[,1:2]))
       min=min(e2$Freq)
       max=max(e2$Freq)
       sliderInput("filternet", "Frequency of occurence of ligand-receptor pairs",
                   min = min, max = max, value = c(min,max),step=2)
     })
   })
   
   #Check if the data is from mouse/human and use the approprite file to list the options for source
   output$source2 <- renderUI({
     file = read.csv("data/param.csv")
     org=as.character(file$organism[file$projects==input$projects])
     if(org=="mouse"){rl=read.csv("data/Mm_PairsLigRec.csv")}else if(org=="human"){rl=read.csv("data/Hs_PairsLigRec.csv")}
     options=as.character(unique(rl$Pair.Source))
     checkboxGroupInput('source2', label='Select source(s)',choices=options,selected=options[1])
   })
   
   #Check if the data is from mouse/human and use the approprite file to list the options for evidence
   output$evidence2 <- renderUI({
     file = read.csv("data/param.csv")
     org=as.character(file$organism[file$projects==input$projects])
     if(org=="mouse"){rl=read.csv("data/Mm_PairsLigRec.csv")}else if(org=="human"){rl=read.csv("data/Hs_PairsLigRec.csv")}
     options=as.character(unique(rl$Pair.Evidence))
     checkboxGroupInput('evidence2',label='Select Evidence(s)',choices=options,selected=options[1])
   })
   
   #For selected project and grouping variable, generate all possible ligand receptor pairs and filter based on user input
   datasetInputnet = reactive({
     result=ligrec(fileload(),pair=input$pairbynet,prj=input$projects)
     edges=result %>% dplyr::select(Receptor_cluster,Lig_cluster)
     e2=as.data.frame(table(edges[,1:2]))
     # validate(need(input$appfil != 0,"Make your selections and click the Apply filter button"))
     # 
     # if(input$appfil == 0)
     #   return()
     # isolate({
     e2=e2[e2$Freq>= input$filternet[1] & e2$Freq<= input$filternet[2],]
     e2$pair=paste(e2$Receptor_cluster,"_",e2$Lig_cluster,sep="")
     result$pair=paste(result$Receptor_cluster,"_",result$Lig_cluster,sep="")
     result=result[result$pair %in% e2$pair,]
     result=result %>% dplyr::select(pairname,receptor,ligand,Pair.Source:Lig_cluster)
     if(input$checksource2==T){result=result[result$Pair.Source %in% input$source2,]}
     if(input$checkevi2==T){result=result[result$Pair.Evidence %in% input$evidence2,]}
     #})
     return(result)
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
   
   #Get nodes 
   nodes = reactive({
     withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
       result=datasetInputnet()
       nodes=as.data.frame(unique(c(result$Receptor_cluster,result$Lig_cluster)))
       colnames(nodes)="id"
       col=cpallette[1:nrow(nodes)]
       nodes$color=col
       nodes$groups=nodes$id
       return(nodes)
     })})
   
   #get edges
   edges = reactive({
     withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
       result=datasetInputnet()
       edges=result %>% dplyr::select(Receptor_cluster,Lig_cluster)
       colnames(edges)=c("from","to")
       e2=as.data.frame(table(edges[,1:2]))
       e2$title=paste(e2$from,"_to_",e2$to,"_freq_",e2$Freq,sep="")
       return(e2)
     })})
   
   #create lig-receptor network plot
   lrnetwork = reactive({
     withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
       e2=edges()
       nodes=nodes()
       # if(input$freeze == 0){
       #   visNetwork(nodes, e2) %>% visEdges(arrows ="to") %>% 
       #     visLegend(useGroups = FALSE, addNodes = data.frame(label = "Nodes", shape = "circle"), 
       #               addEdges = data.frame(label = "Edges", color = "black")) %>% 
       #     visOptions(highlightNearest =list(enabled=TRUE,hover=TRUE,degree=list(from=0.5,to=0.5)), 
       #                nodesIdSelection = list(enabled = TRUE, 
       #                                        style = 'width: 200px; height: 26px;background: #f8f8f8;color: darkblue;border:none;outline:none;')) %>%
       #     visPhysics(stabilization=TRUE,timestep=0.2, adaptiveTimestep = T,barnesHut = list(avoidOverlap=0))
       # }else if(input$freeze != 0){
         visNetwork(nodes, e2) %>% visEdges(arrows ="to") %>% 
           visLegend(useGroups = FALSE, addNodes = data.frame(label = "Nodes", shape = "circle"), 
                     addEdges = data.frame(label = "Edges", color = "black")) %>% 
           visOptions(highlightNearest =list(enabled=TRUE,hover=TRUE,degree=list(from=0.5,to=0.5)), 
                      nodesIdSelection = list(enabled = TRUE, 
                                              style = 'width: 200px; height: 26px;background: #f8f8f8;color: darkblue;border:none;outline:none;')) %>%
           visPhysics(stabilization=TRUE,timestep=0.2, adaptiveTimestep = T,barnesHut = list(avoidOverlap=0,gravitationalConstant=-25000)) %>% 
           visInteraction(dragNodes = FALSE, dragView = FALSE, zoomView = FALSE)
       #}
     })
   })
   
   
   #Render the ligand-receptor network plot
   output$lrnetwork = renderVisNetwork({
     lrnetwork()
   })
   
   observe({
     visNetworkProxy("lrnetwork") %>%
       visStabilize(iterations = 1000)
   })
   

   # get position info
   observeEvent(input$store_position, {
     visNetworkProxy("lrnetwork") %>% visGetPositions()
   })
   
   # format positions
   nodes_positions <- reactive({
     positions <- input$lrnetwork_positions
     if(!is.null(positions)){
       nodes_positions <- do.call("rbind", lapply(positions, function(x){ data.frame(x = x$x, y = x$y)}))
       nodes_positions$id <- names(positions)
       nodes_positions
     } else {
       NULL
     }
   })
     
   #Download plot
   output$dwldnet <- downloadHandler(
     filename = function() {
       paste0(input$projects,"_lig-rec_Network.html",sep="")
     },
     content = function(con) {
       nodes=nodes()
         edges=edges()
       # nodes_positions <- nodes_positions()
       # if(!is.null(nodes_positions)){
       #   nodes_save <- merge(nodes, nodes_positions, by = "id", all = T)
       # } else  {
       #   nodes_save <- nodes
       # }
       visNetwork(nodes, edges) %>% visEdges(arrows ="to") %>% 
         visLegend(useGroups = FALSE, addNodes = data.frame(label = "Nodes", shape = "circle"), 
                   addEdges = data.frame(label = "Edges", color = "black")) %>% 
         visOptions(highlightNearest =list(enabled=TRUE,hover=TRUE,degree=list(from=0.5,to=0.5)), 
                    nodesIdSelection = list(enabled = TRUE, 
                    style = 'width: 200px; height: 26px;background: #f8f8f8;color: darkblue;border:none;outline:none;')) %>%
         visExport(type="pdf",name=paste(input$projects,"_Lig-Rec_network.pdf",sep="")) %>% 
         visPhysics(stabilization=TRUE,timestep=0.2, adaptiveTimestep = T,barnesHut = list(avoidOverlap=0,gravitationalConstant=-25000)) %>%
         visInteraction(dragNodes = FALSE, dragView = FALSE, zoomView = FALSE) %>% visEdges(smooth = FALSE) %>% visSave(con)

     }
   )
   
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
       metadata=metadata %>% select(starts_with("var_"))
       options=colnames(metadata)
       selectInput("pairbyheatnet","Select cell group ",options,selected=options[1])
     })
   })
   
   #Check if the data is from mouse/human and use the approprite file to list the options for source
   output$source3 <- renderUI({
     file = read.csv("data/param.csv")
     org=as.character(file$organism[file$projects==input$projects])
     if(org=="mouse"){rl=read.csv("data/Mm_PairsLigRec.csv")}else if(org=="human"){rl=read.csv("data/Hs_PairsLigRec.csv")}
     options=as.character(unique(rl$Pair.Source))
     checkboxGroupInput('source3', label='Select source(s)',choices=options,selected=options[1])
   })
   
   #Check if the data is from mouse/human and use the approprite file to list the options for evidence
   output$evidence3 <- renderUI({
     file = read.csv("data/param.csv")
     org=as.character(file$organism[file$projects==input$projects])
     if(org=="mouse"){rl=read.csv("data/Mm_PairsLigRec.csv")}else if(org=="human"){rl=read.csv("data/Hs_PairsLigRec.csv")}
     options=as.character(unique(rl$Pair.Evidence))
     checkboxGroupInput('evidence3',label='Select Evidence(s)',choices=options,selected=options[1])
   })
   
   #Generate lig-receptor pairs table
   ligrecheat = reactive({
     result=ligrec(fileload(),pair=input$pairbyheatnet,prj=input$projects)
     result=result %>% dplyr::select(pairname,receptor,ligand,Pair.Source:Lig_cluster)
     if(input$checksource2==T){result=result[result$Pair.Source %in% input$source3,]}
     if(input$checkevi2==T){result=result[result$Pair.Evidence %in% input$evidence3,]}
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
     result=ligrecheat()
     tab=table(result[,6:7])
     tab=as.data.frame(tab)
     tab = tab %>% spread(Receptor_cluster,Freq)
     rownames(tab)=tab$Lig_cluster
     tab= tab %>% dplyr::select(-Lig_cluster)
     #aheatmap(as.matrix(tab),col = colorRampPalette(brewer.pal(n = 9,input$hmpcolnet))(30))
     heatmap.2(as.matrix(tab),col = colorRampPalette(brewer.pal(n = 9,input$hmpcolnet))(30),xlab = "Receptor Cluster",ylab= "Ligand Cluster") 
   })
   #Render heatmap
   output$netheatmap = renderPlot({
     input$hmpcolnet
     netheatmap()
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