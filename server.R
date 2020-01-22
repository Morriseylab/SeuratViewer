library(shiny)
library(shinyBS)
library("AnnotationDbi")
library("org.Mm.eg.db")
library(gage)
library(gageData)
library(RColorBrewer)
library(NMF)
#library(KEGGgraph)
library(Biobase)
library(reshape2)
library(ggplot2)
library(biomaRt)
library(KEGGREST)
library(png)
library(GO.db)
library(d3heatmap)
library(dplyr)
library(tidyr)
library(plotly)
library(shinyjs)
library(htmlwidgets)
library(DT)
library(FactoMineR)
library(factoextra)
library(shinyRGL)
library(rgl)
library(rglwidget)
library(SPIA)
library(ReactomePA)
library(limma)
library(ggrepel)
library(readxl)
library(biomaRt)
library(data.table)
source("functions.R")
#Specify user-ids and passwords
auth=read.csv("data/authentication.csv")
my_username <- auth$user
my_password <- auth$pwd

#Create a theme for all plots.
plotTheme <-theme_bw() + theme(axis.title.x = element_text(face="bold", size=12),
                               axis.text.x  = element_text(angle=35, vjust=0.5, size=12),
                               axis.title.y = element_text(face="bold", size=12),
                               axis.text.y  = element_text(angle=0, vjust=0.5, size=12))

server <- function(input, output, session) {
  
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
  
  ################################################################
  ################################################################
  ####### LOAD EXCEL AND POPULATE DROP DOWN FOR PROJECTS #########
  ################################################################
  ################################################################
  
  #Read the parameter file
  readexcel = reactive({
    user=input$username
    file = read.csv(paste("data/param.csv",sep=""))
    if(user=="allusers"){
      file=file
    }else{
      file=file[file$user==user,]
    }
  })
  
  #Get Project list and populate drop-down
  output$projects = renderUI({
    excel=readexcel()
    prj=excel$projects
    selectInput("projects","Select a project",as.list(sort(as.character(prj))))
  })
  
  ################################################################
  ################# DISPLAY FILE LIST IN DASHBOARD ###############
  ################################################################
  
  #Display file in dashboard
  output$dashdata<- renderTable({
    user=input$username
    file=read.csv('data/param.csv',stringsAsFactors = F)
    if(user=="allusers"){
      file=file
      colnames(file)=c("Project Name","Project Description","Old?(Y/N)","Type(RNA/Microarray)","Username")
      file=file[order(file$`Project Name`),]
    }else{
      file=file[file$user==user,] %>% dplyr::select(-old:-user)
      colnames(file)=c("Project Name","Project Description")
      file=file[order(file$`Project Name`),]
    }
  }, digits = 1)
  
  ###################################################
  ###################################################
  ####### LOAD RDATA FILE AND GET CONTRASTS##########
  ###################################################
  ###################################################
  
  #Load Rdata
  fileload <- reactive({
    inFile = paste('data/',as.character(input$projects),'.RData',sep = '')
    load(inFile)
    loaddata=results
    return(loaddata)
  })
  
  #Get contrast list and populate drop-down
  output$contrasts = renderUI({
    results=fileload()
    lim=results$limma
    contrasts=as.list(as.character(unlist(lapply((names(lim)),factor))))
    selectInput("contrast","Select a comparison",contrasts,"pick one")
  })
  
  
  #######################################################################################################################################################
  #######################################################################################################################################################
  ####################################################################### PCA PLOT ######################################################################
  #######################################################################################################################################################
  #######################################################################################################################################################
  #Populate drop-down for PC to plot on x-axis
  output$pcaxoptions <- renderUI({
    selectInput("pcaxaxes","Select Principle Component to plot on the X-axis ",c(1:10))
  })
  
  #Populate drop-down for PC to plot on y-axis
  output$pcayoptions <- renderUI({
    selectInput("pcayaxes","Select Principle Component to plot on the Y-axis",c(1:10),selected=2)
  })
  
  #PRint the PC's chosen to be plotted
  output$biplottitle <- renderText({
    text=as.character(paste("Dim",input$pcaxaxes," vs Dim",input$pcayaxes,sep=""))
    return(text)
  })
  
  #Textbox to enter number of genes to use to plot
  output$pcipslide <- renderUI({
    textInput(inputId = 'pcipslide', label = "Enter top number of input genes that show maximum variance", value = '500')
  })
  
  #Textbox to enter number of genes to view in the biplot
  output$pcslide <- renderUI({
    textInput(inputId = 'pcslide', label = "Enter number of genes to view in the biplot", value = '0')
  })
  
  #Drop down menu for pc-plot colorby option
  output$pcacolorby = renderUI({ 
    results=fileload()
    eset=results$eset
    pd=pData(eset) #get pheno-data
    pd=pd %>% select(starts_with("var")) #get columns from phenodata that start with "var"
    kt=as.data.frame(t(na.omit(t(pd)))) #omit columns that have only NA's
    bpcols=c("maineffect",colnames(kt))
    selectInput("pcacolorby","Color By",bpcols) #populate drop down menu with the phenodata columns
  })
  
  #Checkbox to view ellipses in the PCA plot
  output$ellipse <- renderUI({
    checkboxInput("ellipse", label = "Check to view ellipses", value = FALSE)
  })
  
  
  #Function for PCA plot
  plotbiplot = reactive({
    res.pca = res_pca()
    x=as.numeric(input$pcaxaxes)
    y=as.numeric(input$pcayaxes)
    results=fileload()
    v = results$eset
    pData<-pData(v)
    colorby=input$pcacolorby
    hab=eval(parse(text = paste0("pData$",colorby,sep="")))
    validate(
      need(input$pcslide, "Enter number of genes to view in biplot")
    )
    if(input$pcslide==0 & input$ellipse==F){
      fviz_pca_ind(res.pca, repel=T,geom='point',label='var',addEllipses=FALSE, habillage = as.factor(hab),pointsize = 3.35,axes=c(x,y))+scale_shape_manual(values = c(rep(19,length(unique(hab)))))+theme(axis.title.x = element_text(face="bold", size=14),
                                                                                                                                                                                                           axis.title.y = element_text(face="bold", size=14),
                                                                                                                                                                                                           legend.text  = element_text(angle=0, vjust=0.5, size=14),
                                                                                                                                                                                                           legend.title  = element_text(angle=0, vjust=0.5, size=14),
                                                                                                                                                                                                           plot.title  = element_text(angle=0, vjust=0.5, size=16))
    }
    else if(input$pcslide==0 & input$ellipse==T){
      fviz_pca_ind(res.pca, repel=T,geom='point',label='var',addEllipses=T,ellipse.type="confidence",ellipse.alpha=0.2, habillage = as.factor(hab),pointsize = 3.35,axes=c(x,y))+scale_shape_manual(values = c(rep(19,length(unique(hab)))))+theme(axis.title.x = element_text(face="bold", size=14),
                                                                                                                                                                                                                                                   axis.title.y = element_text(face="bold", size=14),
                                                                                                                                                                                                                                                   legend.text  = element_text(angle=0, vjust=0.5, size=14),
                                                                                                                                                                                                                                                   legend.title  = element_text(angle=0, vjust=0.5, size=14),
                                                                                                                                                                                                                                                   plot.title  = element_text(angle=0, vjust=0.5, size=16))
      
    }
    
    #fviz_pca_ind(res.pca, geom = c("point", "text"))}
    else if(input$pcslide!=0 & input$ellipse==F){fviz_pca_biplot(res.pca,repel=T, label=c("var","ind"),habillage = as.factor(hab),pointsize = 3.35,axes=c(x,y),select.var = list(contrib = as.numeric(input$pcslide)))+scale_shape_manual(values = c(rep(19,length(unique(hab)))))+theme(axis.title.x = element_text(face="bold", size=14),
                                                                                                                                                                                                                                                                                         axis.title.y = element_text(face="bold", size=14),
                                                                                                                                                                                                                                                                                         legend.text  = element_text(angle=0, vjust=0.5, size=14),
                                                                                                                                                                                                                                                                                         legend.title  = element_text(angle=0, vjust=0.5, size=14),
                                                                                                                                                                                                                                                                                         plot.title  = element_text(angle=0, vjust=0.5, size=16))
    }
    
    else{fviz_pca_biplot(res.pca,repel=T, label=c("var","ind"),addEllipses=T,ellipse.type="confidence",ellipse.alpha=0.1,habillage = as.factor(hab),pointsize = 3.35,axes=c(x,y),select.var = list(contrib = as.numeric(input$pcslide)))+scale_shape_manual(values = c(rep(19,length(unique(hab)))))+theme(axis.title.x = element_text(face="bold", size=14),
                                                                                                                                                                                                                                                                                                           axis.title.y = element_text(face="bold", size=14),
                                                                                                                                                                                                                                                                                                           legend.text  = element_text(angle=0, vjust=0.5, size=14),
                                                                                                                                                                                                                                                                                                           legend.title  = element_text(angle=0, vjust=0.5, size=14),
                                                                                                                                                                                                                                                                                                           plot.title  = element_text(angle=0, vjust=0.5, size=16))
    }
  })
  
  #plotting function for pca plot
  output$biplot = renderPlot({
    plotbiplot()
  })
  
  #Button for dwnloading PCA plot
  output$dwldbiplot = renderUI({
    downloadButton('downloadbiplot', 'Download Biplot')
  }) 
  
  #Download function for pca plot
  output$downloadbiplot <- downloadHandler(
    filename = function() {
      paste0("biplot.pdf")
    },
    content = function(file){
      pdf(file,width=14,height = 9,useDingbats=FALSE)
      plot(plotbiplot())
      dev.off()
    })
  
  ###################################################
  ###################################################
  ########### VARIANCES OF PCA PLOT #################
  ###################################################
  ###################################################
  #Text explaining PCA variances
  output$pcatitle <- renderText({
    text="The proportion of variances retained by the principal components can be viewed in the scree plot. The scree plot is a graph of the eigenvalues/variances associated with components"
    return(text)
  })
  
  #PLot scree plot of all PC's
  output$pcaplot_ip = renderPlot({
    res.pca = res_pca()
    fviz_screeplot(res.pca, ncp=10)
  })
  
  #get expression data and perform PCA
  res_pca = reactive({
    n=as.numeric(input$pcipslide)
    validate(
      need(as.numeric(input$pcipslide) > 199, "Minimum value of input genes that show maximum variance should at least be 200")
    )
    results=fileload()
    v = results$eset
    keepGenes <- v@featureData@data
    pData<-phenoData(v)
    v.filter = v[rownames(v@assayData$exprs) %in% rownames(keepGenes),]
    Pvars <- apply(v.filter@assayData$exprs,1,var)
    select <- order(Pvars, decreasing = TRUE)[seq_len(min(n,length(Pvars)))]
    v.var <-v.filter[select,]
    m<-v.var@assayData$exprs
    rownames(m) <- v.var@featureData@data$SYMBOL
    m=as.data.frame(m)
    m=unique(m)
    res.pca = PCA(t(m), graph = FALSE)
  })
  
  #Extract PCA information like eigan values, variance of each PC 
  pcaplo_tab = reactive({
    res.pca =res_pca()
    eigenvalues = res.pca$eig
    return(eigenvalues)
  })
  
  #Display above PC information in a table
  output$pcaplot_tab = DT::renderDataTable({
    DT::datatable(pcaplo_tab(),
                  extensions = c('Scroller'),
                  options = list(
                    searchHighlight = TRUE,
                    scrollX = TRUE
                  ))
  })
  
  ###################################################
  ###################################################
  ##################3D PCA PLOT #####################
  ###################################################
  ###################################################
  #PLot 3D PCA plot
  output$pcaplot3d = renderRglwidget({
    graphics.off()
    pdf(NULL)
    v=datasetInput3()
    results=fileload()
    pData=pData(results$eset)
    v=t(v)
    v= v[,apply(v, 2, var, na.rm=TRUE) != 0]
    pca <- res_pca()
    vars <- apply(pca$var$coord, 2, var)
    props <- round((vars / sum(vars))*100,1)
    groups=factor(gsub('-','_',pData$maineffect))
    try(rgl.close())
    open3d()
    # resize window
    par3d(windowRect = c(100, 100, 612, 612))
    palette(c('blue','red','green','orange','cyan','black','brown','pink'))
    plot3d(pca$ind$coord[,1:3], col =as.numeric(groups), type='s',alpha=1.75,axes=F,
           xlab=paste('PC1 (',props[1],'%)',sep=''),
           ylab=paste('PC2 (',props[2],'%)',sep=''),
           zlab=paste('PC3 (',props[3],'%)',sep='')
    )
    axes3d(edges=c("x--", "y--", "z"), lwd=2, expand=10, labels=FALSE,box=T)
    grid3d("x")
    grid3d("y")
    grid3d("z")
    l=length(levels(groups))
    ll=1:l
    y=1+(ll*15)
    legend3d("topright", legend = levels(groups), pch = 16, col=palette(),cex=1, inset=c(0.02))
    rglwidget()
  })
  
  ###################################################
  ###################################################
  ######## GET PROJECT DESC AND DISPLAY   ###########
  ###################################################
  ###################################################
  #Read parameter file and get project description for the project selected
  prjdesc = reactive({
    file = readexcel()
    prj=input$projects
    desc=file$desc[file$projects %in% prj]
    desc=as.character(desc)
  })
  
  #Display text in main project description panel
  output$pdesc <- renderText({
    desc=prjdesc()
  })
  
  #######################################################################################################################################################
  #######################################################################################################################################################
  ####################################################################### DOT PLOT ######################################################################
  #######################################################################################################################################################
  #######################################################################################################################################################
  
  #Drop down menu for dot-plot x-axis grouping
  output$boxplotcol = renderUI({ 
    results=fileload()
    eset=results$eset
    pData=pData(eset) #get pheno-data
    kc=pData[ , grepl( "var_" , colnames(pData) ) ] #get columns from phenodata that start with "var"
    kt=as.data.frame(t(na.omit(t(kc)))) #omit columns that have only NA's
    kc=data.frame(maineffect=pData$maineffect,sample_name=pData$sample_name,kt) #create a dataframe with maineffect and sample_name and non-null columns starting with var
    bpcols=as.list(as.character(unlist(lapply((colnames(kc)),factor))))
    selectInput("color","Select an Attribute for the X-axis",bpcols) #populate drop down menu with the phenodata columns
  })
  
  #Drop down menu for dot-plot color
  output$boxplotcol2 = renderUI({ 
    results=fileload()
    eset=results$eset
    pData=pData(eset) #get pheno-data
    kc=pData[ , grepl( "var_" , colnames(pData) ) ] #get columns from phenodata that start with "var"
    kt=as.data.frame(t(na.omit(t(kc)))) #omit columns that have only NA's
    kc=data.frame(maineffect=pData$maineffect,sample_name=pData$sample_name,kt) #create a dataframe with maineffect and sample_name and non-null columns starting with var
    bpcols=as.list(as.character(unlist(lapply((colnames(kc)),factor))))
    selectInput("color2","Color By",bpcols) #populate drop down menu with the phenodata columns
  })
  
  #Checkbox for whether or not to display the minimum expression line 
  output$minexprline = renderUI({ 
    tagList(
      checkboxInput("minexprline", label = "Show expression threshold line", value = FALSE),
      bsTooltip("minexprline","Please note that not all projects have this option currently", placement = "bottom", trigger = "hover",options = NULL)
    )
  })
  
  #Extract expression data to create dot-plot
  dotplot_out = reactive({
    s = input$table_rows_selected #select rows from table
    dt = datasetInput() #load limma data
    dt$id=rownames(dt)
    dt=data.frame(dt$id,dt[,-ncol(dt)])
    validate(
      need((is.data.frame(dt) && nrow(dt))!=0, "No data in table")
    )
    dt1 = dt[s, , drop=FALSE]#get limma data corresponding to selected row in table
    id = as.character(dt[s,1]) 
    results=fileload()
    eset <- results$eset
    pData=pData(eset) #get pheno-data
    if(is.factor(pData$sample_name)==T){lev=levels(pData$sample_name)}
    minexpr=pData$minexpr[1]
    signal=as.data.frame(eset@assayData$exprs[id,])
    colnames(signal)="signal"
    signal$id=rownames(signal)
    e=left_join(pData,signal,by=c('sample_name'='id'))
    if(is.factor(pData$sample_name)==T){e$sample_name= factor(e$sample_name, levels = levels(pData$sample_name))}
    if(is.na(dt1$SYMBOL)) #if gene symbol does not exist,use ENSEMBL id
    {genesymbol=dt1$ENSEMBL}
    else{
      genesymbol=dt1$SYMBOL} #get the gene symbol of the row selected
    if(input$minexprline==T){
      gg=ggplot(e,aes_string(x=input$color,y="signal",col=input$color2))+plotTheme+guides(color=guide_legend(title=as.character(input$color2)))+
        labs(title=genesymbol, x="Condition", y="Expression Value") + geom_point(size=5,position=position_jitter(w = 0.1))+ geom_smooth(method=lm,se=FALSE) +
        stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax= "mean", size= 0.3, geom = "crossbar",width=.2) + geom_hline(yintercept=minexpr, linetype="dashed", color = "red")}
    else{
      gg=ggplot(e,aes_string(x=input$color,y="signal",col=input$color2))+plotTheme+guides(color=guide_legend(title=as.character(input$color2)))+
        labs(title=genesymbol, x="Condition", y="Expression Value") + geom_point(size=5,position=position_jitter(w = 0.1))+ geom_smooth(method=lm,se=FALSE) +
        stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax= "mean", size= 0.3, geom = "crossbar",width=.2)
    }
    gg
  })
  
  # plot dotplot
  output$dotplot = renderPlot({
    dotplot_out()
  })
  
  #function to download dot plot
  output$downloaddotplot <- downloadHandler(
    filename = function() {
      paste0(input$projects, '_dotplot.jpg', sep='') 
      #paste0("dotplot.jpg")
    },
    content = function(file){
      jpeg(file, quality = 100, width = 800, height = 800)
      plot(dotplot_out())
      dev.off()
    })
  
  ###################################################
  ###################################################
  ###########LOAD LIMMA FILE AND DISPLAY#############
  ###################################################
  ###################################################
  #Read limma data from eset
  datasetInput0.5 = reactive({
    contrast=input$contrast
    results=fileload()
    k=paste('results$limma$',contrast,sep='')
    limmadata=eval(parse(text = k))
  })
  
  #Update limma results based on gene selection (upregulated, downregulated, both or none)
  datasetInput = reactive({
    contrast=input$contrast #select contrast
    limmadata=datasetInput0.5()
    limmadata=subset(limmadata, select=-c(logFC))
    lfc=as.numeric(input$lfc) #get logFC
    apval=as.numeric(input$apval)#get adjusted P.Vals
    if(is.null(input$radio))
    {
      d = limmadata
    }
    else if(input$radio=='none')
    {
      d=limmadata
    }
    else if(input$radio=='down')
    {
      d=limmadata
      d = d[which(d$fc < (-1*(lfc)) & d$adj.P.Val < apval),]
    }
    else if(input$radio=='up')
    {
      d=limmadata
      d = d[which(d$fc>lfc & d$adj.P.Val < apval),]
    }
    else if(input$radio=='both')
    {
      d=limmadata
      d = d[which(abs(d$fc) > lfc & d$adj.P.Val < apval),]
    }
    #limma = as.data.frame(d) #load limma data
    geneid=d$SYMBOL
    url= paste("http://www.genecards.org/cgi-bin/carddisp.pl?gene=",geneid,sep = "")
    if(url=="http://www.genecards.org/cgi-bin/carddisp.pl?gene="){
      d$link<-NULL
    }else{
      d$link=paste0("<a href='",url,"'target='_blank'>","Link to GeneCard","</a>")}
    d=as.data.frame(d) 
    #d=as.data.frame(d) %>% select(-t)
    return(d)
  })
  
  #print limma results in data table
  output$table = DT::renderDataTable({
    input$lfc
    input$apval
    input$project
    input$contrast
    DT::datatable(datasetInput(),
                  extensions = 'Buttons', options = list(
                    dom = 'Bfrtip',
                    buttons = list()),
                  rownames=FALSE,selection = list(mode = 'single', selected =1),escape=FALSE)
  })
  
  #Display text (contrast name) above limma table
  output$contrdesc <- renderText({
    contrastname=input$contrast
    text=paste('CONTRAST:  ',contrastname,sep="   ")
    return(text)
  })
  
  #download limma results data as excel sheet
  output$dwld <- downloadHandler(
    filename = function() { paste(input$projects, '.csv', sep='') },
    content = function(file) {
      write.csv(datasetInput(), file)
    })
  
  ###################################################
  ###################################################
  ############# DISPLAY VOLCANO PLOT  ###############
  ###################################################
  ###################################################
  #Get limma data
  datasetInputvol = reactive({
    limmadata=datasetInput()
    return(limmadata)
  })
  
  #Drop down to choose what genes to display on volcano plot
  output$volcdrop <- renderUI({
    selectInput("volcdrop", "Select input type",c('Significant genes' = "signi",'GO genes' = "go"))
    
  })
  
  #Slider to choose number of genes to display on volcano plot
  output$volcslider <- renderUI({
    conditionalPanel(
      condition = "input.volcdrop == 'signi'",
      fluidRow(
        column(6,sliderInput("volcslider", label = h4("Select top number of genes"), min = 0,max = 25, value = 5))
      ))
    
  })
  
  #Function to assign values to volcano plot points
  vpt = reactive({
    diff_df=datasetInput0.5()
    
    FDR=input$apval
    lfc=input$lfc
    
    if(input$volcdrop=="signi"){
      diff_df$group <- "NotSignificant"
      # change the grouping for the entries with significance but not a large enough Fold change
      diff_df[which(diff_df['adj.P.Val'] < FDR & abs(diff_df['logFC']) < lfc ),"group"] <- "Filtered by FDR"
      
      # change the grouping for the entries a large enough Fold change but not a low enough p value
      diff_df[which(diff_df['adj.P.Val'] > FDR & abs(diff_df['logFC']) > lfc ),"group"] <- "Filtered by FC"
      
      # change the grouping for the entries with both significance and large enough fold change
      diff_df[which(diff_df['adj.P.Val'] < FDR & abs(diff_df['logFC']) > lfc ),"group"] <- "Significant (Filtered by both FDR and FC)"
    }
    else if(input$volcdrop=="go"){
      top_peaks2=GOHeatup()
      diff_df$group <- "All genes"
      diff_df[which(diff_df$SYMBOL %in% top_peaks2$SYMBOL ),"group"] <- "Selected_genes"
    }
    return(diff_df)
  })
  
  #Function to draw the volcano plot
  volcanoplot_out = reactive({
    diff_df=vpt()
    
    if(input$volcdrop=="signi"){
      # Find and label the top peaks..
      n=input$volcslider
      if(n>0){
        top_peaks <- diff_df[with(diff_df, order(adj.P.Val,logFC)),][1:n,]
        top_peaks <- rbind(top_peaks, diff_df[with(diff_df, order(adj.P.Val,-logFC)),][1:n,])
        
        a <- list()
        for (i in seq_len(nrow(top_peaks))) {
          m <- top_peaks[i, ]
          a[[i]] <- list(x = m[["logFC"]],y = -log10(m[["adj.P.Val"]]),text = m[["SYMBOL"]],xref = "x",yref = "y",showarrow = FALSE,arrowhead = 0.5,ax = 20,ay = -40)
        }
        p <- plot_ly(data = diff_df, x = diff_df$logFC, y = -log10(diff_df$adj.P.Val),text = diff_df$SYMBOL, mode = "markers", color = diff_df$group) %>% layout(title ="Volcano Plot",xaxis=list(title="Log Fold Change"),yaxis=list(title="-log10(FDR)")) %>%
          layout(annotations = a)
      }
      else{
        p <- plot_ly(data = diff_df, x = diff_df$logFC, y = -log10(diff_df$adj.P.Val),text = diff_df$SYMBOL, mode = "markers", color = diff_df$group) %>% layout(title ="Volcano Plot",xaxis=list(title="Log Fold Change"),yaxis=list(title="-log10(FDR)"))
      }
    }
    else if(input$volcdrop=="go"){
      # Find and label the top peaks..
      top_peaks <- diff_df[diff_df$SYMBOL %in% top_peaks2$SYMBOL,]
      a <- list()
      for (i in seq_len(nrow(top_peaks))) {
        m <- top_peaks[i, ]
        a[[i]] <- list(x = m[["logFC"]],y = -log10(m[["adj.P.Val"]]),text = m[["SYMBOL"]],xref = "x",yref = "y",showarrow = FALSE,arrowhead = 0.5,ax = 20,ay = -40)
      }
      p <- plot_ly(data = diff_df, x = diff_df$logFC, y = -log10(diff_df$adj.P.Val),text = diff_df$SYMBOL, mode = "markers", color = diff_df$group) %>% layout(title ="Volcano Plot",xaxis=list(title="Log Fold Change"),yaxis=list(title="-log10(FDR)")) 
    }
    p
  })
  
  #Make non-interactive plot for volcano plot download
  volcanoplot_dout = reactive({
    diff_df=vpt()
    if(input$volcdrop=="signi"){
      n=input$volcslider
      if(n>0){
        top_peaks <- diff_df[with(diff_df, order(adj.P.Val,logFC)),][1:n,]
        top_peaks <- rbind(top_peaks, diff_df[with(diff_df, order(adj.P.Val,-logFC)),][1:n,])
        p <- ggplot(data = diff_df, aes(x = diff_df$logFC, y = -log10(diff_df$adj.P.Val))) + geom_point_rast(aes(color=diff_df$group)) +ggtitle("Volcano Plot") + xlab("Log Fold Change") + ylab("-log10(FDR)") +labs(color="")+ geom_label_repel(data=top_peaks,aes(x = top_peaks$logFC, y = -log10(top_peaks$adj.P.Val),label=top_peaks$SYMBOL)) + theme_bw()
      }
      else{
        p <- ggplot(data = diff_df, aes(x = diff_df$logFC, y = -log10(diff_df$adj.P.Val))) + geom_point_rast(aes(color=diff_df$group)) +ggtitle("Volcano Plot") + xlab("Log Fold Change") + ylab("-log10(FDR)") +labs(color="") + theme_bw()
      }
    }
    else if(input$volcdrop=="go"){
      top_peaks <- diff_df[diff_df$SYMBOL %in% top_peaks2$SYMBOL,]
      p <- ggplot(data = diff_df, aes(x = diff_df$logFC, y = -log10(diff_df$adj.P.Val))) + geom_point_rast(aes(color=diff_df$group)) +ggtitle("Volcano Plot") + xlab("Log Fold Change") + ylab("-log10(FDR)") +labs(color="") + theme_bw()
    }
    p
  })
  
  #Render and display interactive volcano plot
  output$volcanoplot = renderPlotly({
    input$radio
    input$lfc
    input$apval
    input$volcslider
    input$volcdrop
    volcanoplot_out()
  })
  
  #Display limma results 
  output$table_volc = DT::renderDataTable({
    DT::datatable(datasetInput(),
                  extensions = c('Buttons','Scroller'),
                  options = list(dom = 'Bfrtip',
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                 scrollX = TRUE,
                                 buttons = c('copy', 'print')
                  ),rownames=TRUE,selection = list(mode = 'single', selected =1),escape=FALSE)
  })
  
  
  #Download non-interactive volcano plot
  output$dwldvolcanoplot <- downloadHandler(
    filename = function() {
      paste0("volcano.pdf")
    },
    content = function(file){
      pdf(file,width=14,height = 9,useDingbats=FALSE)
      plot(volcanoplot_dout())
      dev.off()
    })
  
  ###################################################
  ###################################################
  #######CONDITIONAL PANEL FOR Limma ################
  ###################################################
  ###################################################
  
  #Create checkboxes with contrasts corresponding to the project (displayed only when multiple contrast checkbox is selected)
  output$contrastslimma <- renderUI({
    results=fileload()
    lim=results$limma
    contrasts=as.list(as.character(unlist(lapply((names(lim)),factor))))
    checkboxGroupInput("multicontrast",label="Pick Contrasts",choices=contrasts)
  })
  
  #create table with p.value and FC value for the contrasts selected
  multilimma = reactive({
    validate(
      need(input$multicontrast, "Please Select at least one comparison ")
    )
    contr=input$multicontrast
    results=fileload()
    full_limma = data.frame(id=as.character())
    for(i in 1:length(contr)){
      k=paste('results$limma$',contr[i],sep='')
      limmadata=eval(parse(text = k))
      limmadata2=data.frame(id=rownames(limmadata),logFC=limmadata$logFC,adj.P.Val=limmadata$adj.P.Val)
      colnames(limmadata2)[-1]=paste(colnames(limmadata2[,c(-1)]),contr[i], sep = "_")
      full_limma=full_join(full_limma,limmadata2,by='id')
    }
    k=data.frame(id=rownames(limmadata),SYMBOL=limmadata$SYMBOL)
    m=full_join(k,full_limma,by='id')
    return(m)
  })
  
  
  #update table with the dataframe
  output$table_TRUE = DT::renderDataTable({
    input$project
    input$contrast
    DT::datatable(multilimma(),
                  extensions = c('Buttons','Scroller'),
                  options = list(dom = 'Bfrtip',
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                 scrollX = TRUE,
                                 buttons = list('copy')
                  ),rownames=TRUE,selection = list(mode = 'single', selected =1),escape=FALSE)
  })
  
  #action button to download the table
  output$dwldmultitab = renderUI({
    downloadButton('multidwld','Download Table')
  }) 
  
  #fucntion to download multi-contrast limma table
  output$multidwld <- downloadHandler(
    filename = function() { paste(input$projects, '_multiple_contrasts.csv', sep='') },
    content = function(file) {
      write.csv(multilimma(), file,row.names=FALSE)
    })
  
  #######################################################################################################################################################
  #######################################################################################################################################################
  ####################################################### DISPLAY RAW EXPRESSION (VOOM) DATA ############################################################
  #######################################################################################################################################################
  #######################################################################################################################################################
  
  #load voom data from eset
  datasetInput3 = reactive({
    results=fileload()
    exprsdata=results$eset@assayData$exprs
  })
  
  #annotate voom data using featuresdata 
  datasetInput33 = reactive({
    results=fileload()
    exprsdata=as.data.frame(results$eset@assayData$exprs)
    features=as.data.frame(pData(featureData(results$eset)))
    features$id=rownames(features)
    exprsdata$id=rownames(exprsdata)
    genes <- inner_join(features,exprsdata,by=c('id'='id'))
    return(genes)
  })
  
  #print voom or expression data file
  output$table3 = DT::renderDataTable({
    DT::datatable(datasetInput33(),
                  extensions = c('Buttons','Scroller'),
                  options = list(dom = 'Bfrtip',
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                 scrollX = TRUE,
                                 buttons = c('copy', 'print')
                  ),rownames=FALSE,caption= "Voom data")
  })
  
  #action button to download the raw expression matrix
  output$dwldrawtab = renderUI({
    downloadButton('rawdwld','Download Raw Data')
  }) 
  
  #fucntion to download voom expression data table
  output$rawdwld <- downloadHandler(
    filename = function() { paste(input$projects, '_rawdata.csv', sep='') },
    content = function(file) {
      write.csv(datasetInput33(), file,row.names=FALSE)
    })
  
  #######################################################################################################################################################
  #######################################################################################################################################################
  ################################################################ DISPLAY PHENO DATA ###################################################################
  #######################################################################################################################################################
  #######################################################################################################################################################
  
  #load pheno from eset
  phenofile = reactive({
    results=fileload()
    pd=pData(results$eset)
    if("minexpr" %in% colnames(pData)){
      pd=pd %>% dplyr::select(-minexpr)
    }
    else{pd=pd}
  })
  
  #print pheno data file 
  output$phenofile = DT::renderDataTable({
    DT::datatable(phenofile(),
                  extensions = c('Buttons','Scroller'),
                  options = list(dom = 'Bfrtip',
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                 scrollX = TRUE,
                                 buttons = c('copy', 'print')
                  ),rownames=FALSE,caption= "Sample data")
  })
  
  #######################################################################################################################################################
  #######################################################################################################################################################
  ############################################################### CAMERA OUTPUT DISPLAY #################################################################
  #######################################################################################################################################################
  #######################################################################################################################################################
  
  #populate camera dropdown menu in the sidebar with the genesets based on the project RData
  output$cameradd = renderUI({
    results=fileload()
    contrast=input$contrast
    cam=paste("results$camera$",contrast,sep="")
    cam=eval(parse(text=cam))
    cameradd=as.list(as.character(unlist(lapply((names(cam)),factor))))
    selectInput("cameradd","Select a Gene Set",cameradd)
  })
  
  #Get camera data from Rdata file for the chosen contrast
  geneid = reactive({
    results=fileload()
    cameradd=input$cameradd
    contrast=input$contrast #get user input for contrast/comparison
    c=paste('results$camera$',contrast,'$',cameradd,'$camera_result',sep='') #get camera data corresponding to the contrast chosen
    cam=eval(parse(text = c)) #convert string to variable
    cam=data.frame(name=rownames(cam),cam)
    name=cam$name
    if (cameradd == "GO")
    {
      url= paste("http://amigo.geneontology.org/amigo/term/",name,sep = "") #create link to Gene Ontology Consortium
      cam$link=paste0("<a href='",url,"'target='_blank'>","Link to Gene Ontology Consortium","</a>")
      cam=as.data.frame(cam)
    }else{
      url= paste("http://software.broadinstitute.org/gsea/msigdb/cards/",name,".html",sep = "")
      cam$link=paste0("<a href='",url,"'target='_blank'>","Link to Molecular Dignature Database","</a>")
      cam=as.data.frame(cam)}
    return(cam) # return datatable with camera results
  })
  
  # print out camera results in a table
  output$tablecam = DT::renderDataTable({
    input$camera
    input$cameradd
    input$contrast
    isolate({
      DT::datatable(geneid(),
                    extensions = c('Buttons','Scroller'),
                    options = list(dom = 'Bfrtip',
                                   searchHighlight = TRUE,
                                   pageLength = 10,
                                   lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                   scrollX = TRUE,
                                   buttons = c('copy', 'print')
                    ),rownames= FALSE,selection = list(mode = 'single', selected =1),escape=FALSE,caption = "Camera Results")
    })
  })
  
  #Generate text title for the gene list table
  output$camdesc <- renderText({
    s = input$tablecam_rows_selected
    dt = geneid() 
    dt = as.character(dt[s, , drop=FALSE]) 
    camname=dt[1]
    text=paste('Gene list for Camera term :',camname,sep="")
    return(text)
  })
  
  #get the gene-list for every row in camera results table
  campick2 = reactive({
    results=fileload()
    cameradd=input$cameradd
    contrast=input$contrast #get user input for contrast/comparison
    c=paste('results$camera$',contrast,'$',cameradd,'$indices',sep='') #get camera indices corresponding to the contrast chosen
    cameraind=eval(parse(text = c))
    cam=geneid() #get datatable with camera data from reactive
    s=input$tablecam_rows_selected # get  index of selected row from table
    cam=cam[s, ,drop=FALSE]
    res=datasetInput0.5()
    res2=datasetInput33()
    if("ENTREZID" %in% colnames(res2)){
      res2=res2
    }
    else{res2=res}
    
    #get gene list from indices
    if (cameradd == "GO")
    {
      k=paste('res2$ENTREZID[cameraind$`',cam$name,'`]',sep='')}
    else{
      k=paste('res2$ENTREZID[cameraind$',cam$name,']',sep='')
    }
    genes=eval(parse(text = k)) #get entrez id's corresponding to indices
    genesid=res[res$ENTREZID %in% genes,] #get limma data corresponding to entrez id's
    return(data.frame(genesid)) #return the genelist
  })
  
  #print data table with gene list corresponding to each row in camera datatable
  output$campick3 = DT::renderDataTable({
    input$cameradd
    input$contrast
    input$projects
    DT::datatable(campick2(),
                  extensions = c('Buttons','Scroller'),
                  options = list(dom = 'Bfrtip',
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                 scrollX = TRUE,
                                 buttons = c('copy', 'print')
                  ),rownames=FALSE,escape=FALSE,caption="GENE LIST")
  })
  
  
  #download camera datatable
  output$downloadcam <- downloadHandler(
    filename = function() { paste('Camera_',input$projects,'_',input$contrast,'.csv', sep='') },
    content = function(file) {
      write.csv(geneid(), file)
    })
  
  ###################################################
  ###################################################
  ###### CREATE ENRICHMENT PLOT FROM CAMERA #########
  ###################################################
  ###################################################
  #Create enrichment plot for the camera term
  eplotcamera = reactive({
    results=fileload()
    cameradd=input$cameradd
    contrast=input$contrast #get user input for contrast/comparison
    s = input$camres_rows_selected
    dt = geneid() 
    dt = as.character(dt[s, , drop=FALSE]) 
    cat= dt$name
    c=paste('results$camera$',contrast,'$',cameradd,'$indices$',category,sep='') #get camera indices corresponding to the contrast chosen
    cameraind=eval(parse(text = c))
    exprsdata=as.data.frame(results$eset@assayData$exprs)
    features=as.data.frame(pData(featureData(results$eset)))
    features$id=rownames(features)
    exprsdata$id=rownames(exprsdata)
    res2<- inner_join(features,exprsdata,by=c('id'='id'))
    k=res2$ENTREZID[cameraind]
    limma_all=datasetInput0.5()
    #limma=limma[limma$ENTREZID %in% k,]
    
  })
  
  #Render enrichment plot
  output$eplotcamera = renderPlot({
    eplotcamera()
  })
  
  # print out camera results in a table
  output$camres = DT::renderDataTable({
    input$camera
    input$cameradd
    input$contrast
    isolate({
      DT::datatable(geneid(),
                    extensions = c('Buttons','Scroller'),
                    options = list(dom = 'Bfrtip',
                                   searchHighlight = TRUE,
                                   pageLength = 10,
                                   lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                   scrollX = TRUE,
                                   buttons = c('copy', 'print')
                    ),rownames= FALSE,selection = list(mode = 'single', selected =1),escape=FALSE,caption = "Camera Results")
    })
  })
  
  ###################################################
  ###################################################
  ######### CREATE HEATMAP FROM CAMERA ##############
  ###################################################
  ###################################################
  #extract voom expression data of all genes corresponding to selected row in camera datatable
  heatmapcam <- reactive({
    genesid=campick2()  #gene list from camera
    voom=as.data.frame(datasetInput3())#voom data
    genes_cam<-voom[rownames(voom) %in% rownames(genesid),]
    
  })
  
  #Set limit for number of genes that can be viewed in the heatmap
  output$hmplimcam <- renderUI({
    pval=campick2() 
    top_expr=datasetInput3()
    top_expr=top_expr[rownames(top_expr) %in% rownames(pval),]
    mx=nrow(top_expr)
    sliderInput("hmplimcam", label = h5("Select number of genes to view in the heatmap"), min = 2,max =mx, value = mx)
  })
  
  #Create scale for heatmap
  output$hmpscale_out2 = renderPlot({
    hmpscaletest(hmpcol=input$hmpcol2,voom=datasetInput3(),checkbox=input$checkbox2)
  })
  
  #create heatmap for heatmap
  camheatmap = reactive({
    dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    top_expr=heatmapfun(results=fileload(),expr=heatmapcam(),pval=campick2(),file = readexcel(),prj=input$projects,hmplim=input$hmplimcam,hmpsamp=input$hmpsamp2,
                        contrast=input$contrast)
    sym=rownames(top_expr)
    #Remove rows that have variance 0 (This will avoid the Na/Nan/Inf error in heatmap)
    ind = apply(top_expr, 1, var) == 0
    top_expr <- top_expr[!ind,]
    if(input$checkbox2==TRUE){
      d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterby2,xaxis_font_size = 10,colors = colorRampPalette(brewer.pal(n = 9, input$hmpcol2))(30),labRow = sym)}
    else{d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterby2,xaxis_font_size = 10,colors = colorRampPalette(rev(brewer.pal(n = 9, input$hmpcol2)))(30),labRow = sym)}
  })
  
  # Render heatmap for camera genes
  output$camheatmap <- renderD3heatmap({
    input$hmpcol #user input-color palette
    input$clusterby #user input-cluster by
    input$checkbox #user input-reverse colors
    input$gene #user input-slider input for number of genes
    input$genelist
    input$makeheat
    input$gage
    input$go_dd
    input$table4_rows_selected
    input$tablecam_rows_selected
    input$projects
    input$contrast
    input$cameradd
    input$hmpsamp2
    input$hmplimcam
    camheatmap()
  })
  
  #Create non-interactive heatmap for download
  camheatmapalt = reactive({
    dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    top_expr=heatmapfun(results=fileload(),expr=heatmapcam(),pval=campick2(),file = readexcel(),prj=input$projects,hmplim=input$hmplimcam,hmpsamp=input$hmpsamp2,
                        contrast=input$contrast)
    sym=rownames(top_expr)
    #Remove rows that have variance 0 (This will avoid the Na/Nan/Inf error in heatmap)
    ind = apply(top_expr, 1, var) == 0
    top_expr <- top_expr[!ind,]
    if(input$checkbox2==TRUE){
      aheatmap(as.matrix(top_expr),distfun=dist2,scale="row",Rowv=TRUE,Colv=TRUE,fontsize = 10,color = colorRampPalette(brewer.pal(n = 9, input$hmpcol2))(30),labRow = sym)}
    else{aheatmap(as.matrix(top_expr),distfun=dist2,scale="row",Rowv=TRUE,Colv=TRUE,fontsize = 10,color = colorRampPalette(rev(brewer.pal(n = 9, input$hmpcol2)))(30),labRow = sym)}
  })
  
  #Download camera heatmap
  output$downloadcamheatmap <- downloadHandler(
    filename = function(){
      paste0('camera_heatmap','.pdf',sep='')
    },
    content = function(file){
      pdf(file,width=9,height = 14,useDingbats=FALSE, onefile = F)
      camheatmapalt()
      dev.off()
    })
  
  ########################################################################################################################################################
  ########################################################################################################################################################
  ################################################################## SPIA PATHWAY ANALYSIS################################################################
  ########################################################################################################################################################
  ########################################################################################################################################################
  #For the chosen contrast, get SPIA results from the RData 
  spia_op <- reactive({
    results=fileload()
    contrast=input$contrast #get user input for contrast/comparison
    c=paste('results$spia$',contrast,sep='') #get SPIA data corresponding to the contrast chosen
    sp=eval(parse(text = c)) #convert string to variable
    spia_result=data.frame(sp)
    validate(
      need(nrow(spia_result) > 1, "No Results")
    )
    spia_result$KEGGLINK <- paste0("<a href='",spia_result$KEGGLINK,"' target='_blank'>","Link to KEGG","</a>")
    return(spia_result) 
  })
  
  #Display SPIA results in a table
  output$spiaop <- DT::renderDataTable({
    input$runspia
    input$contrast
    input$projects
    isolate({
      DT::datatable(spia_op(),escape = FALSE,selection = list(mode = 'single', selected =1),
                    extensions = c('Buttons','Scroller'),
                    options = list(
                      dom = 'RMDCT<"clear">lfrtip',
                      searchHighlight = TRUE,
                      pageLength = 10,
                      lengthMenu = list(c(5, 10, 15, 20, 25, -1), c('5', '10', '15', '20', '25', 'All')),
                      scrollX = TRUE,
                      buttons = c('copy', 'print')
                    ),rownames=FALSE)
    })
  })
  
  #Display the SPIA term selected from table above the genelist
  output$spiadesc <- renderText({
    s = input$spiaop_rows_selected
    dt = spia_op() 
    dt = dt[s, , drop=FALSE]
    camname=dt$Name
    text=paste('Gene list for SPIA term :',camname,'-',dt[2],sep="")
    return(text)
  })
  
  #Get genelist for SPIA term selected from the table of SPIA results
  spiagenes = reactive({
    spiaid=spia_op() 
    final_res=datasetInput()
    s=input$spiaop_rows_selected 
    row=spiaid[s, ,drop=FALSE]
    results=fileload()
    pd=pData(results$eset)
    org=unique(pd$organism)
    if(org %in% c("Mus musculus", "Mouse", "Mm","Mus_musculus", "mouse")){
      id=paste("mmu",row$ID,sep="")
      allgenelist=keggLink("mmu",id) #for each kegg id, get gene list
    }else{
      id=paste("hsa",row$ID,sep="")
      allgenelist=keggLink("hsa",id) #for each kegg id, get gene list
    }
    p=strsplit(allgenelist,":")
    genes_entrez=sapply(p,"[",2)
    genelist=final_res[final_res$ENTREZID %in% genes_entrez,]
    return(genelist) #return the genelist
  })
  
  #Render table to display the genelist per SPAI term
  output$spiagenes = DT::renderDataTable({
    DT::datatable(spiagenes(),
                  extensions = c('Buttons','Scroller'),
                  options = list(dom = 'Bfrtip',
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                 scrollX = TRUE,
                                 buttons = c('copy', 'print')
                  ),rownames=FALSE,escape=FALSE,selection = list(mode = 'single', selected =1,caption="Genelist"))
  })
  
  #Download function to download SPIA results as a csv file
  output$dwldspia <- downloadHandler(
    filename = function() { paste(input$projects,'_',input$contrast, '_spia.csv', sep='') },
    content = function(file) {
      write.csv(spia_op(), file)
    })
  
  #######################################################
  #######################################################
  ######### CREATE HEATMAP FROM SPIA ####################
  #######################################################
  #######################################################
  
  #extract voom expression data of all genes corresponding to selected row in spia datatable
  heatmapspia <- reactive({
    genesid=spiagenes()  #gene list from camera
    voom=as.data.frame(datasetInput3())#voom data
    genes_spia<-voom[rownames(voom) %in% rownames(genesid),]
    
  })
  
  #get max and min genes per SPIA term to show on slider
  output$hmplimspia <- renderUI({
    pval=spiagenes() 
    top_expr=datasetInput3()
    top_expr=top_expr[rownames(top_expr) %in% rownames(pval),]
    mx=nrow(top_expr)
    sliderInput("hmplimspia", label = h5("Select number of genes to view in the heatmap"), min = 2,max =mx, value = mx)
  })
  
  #Generate a heatmap color scale
  output$hmpscale_out2spia = renderPlot({
    hmpscaletest(hmpcol=input$hmpcolspia,voom=datasetInput3(),checkbox=input$checkboxspia)
  })
  
  #Function to generate d3 camera heatmap
  camheatmap = reactive({
    dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    top_expr=heatmapfun(results=fileload(),expr=heatmapcam(),pval=campick2(),file = readexcel(),prj=input$projects,hmplim=input$hmplimcam,hmpsamp=input$hmpsamp2,
                        contrast=input$contrast)
    sym=rownames(top_expr)
    #Remove rows that have variance 0 (This will avoid the Na/Nan/Inf error in heatmap)
    ind = apply(top_expr, 1, var) == 0
    top_expr <- top_expr[!ind,]
    if(input$checkbox2==TRUE){
      d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterby2,xaxis_font_size = 10,colors = colorRampPalette(brewer.pal(n = 9, input$hmpcol2))(30),labRow = sym)}
    else{d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterby2,xaxis_font_size = 10,colors = colorRampPalette(rev(brewer.pal(n = 9, input$hmpcol2)))(30),labRow = sym)}
  })
  
  # Render SPIA heatmap 
  output$spiaheatmap <- renderD3heatmap({
    input$hmpcolspia #user input-color palette
    input$clusterbyspia #user input-cluster by
    input$checkboxspia #user input-reverse colors
    input$gene #user input-slider input for number of genes
    input$genelist
    input$spiaop_rows_selected
    input$projects
    input$contrast
    input$hmpsamp2spia
    input$hmplimspia
    spiaheatmap()
  })
  
  #create SPIA heatmap function
  spiaheatmap <- reactive({
    dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    top_expr=heatmapfun(results=fileload(),expr=heatmapspia(),pval=spiagenes(),file = readexcel(),prj=input$projects,hmplim=input$hmplimspia,hmpsamp=input$hmpsamp2spia,
                        contrast=input$contrast)
    sym=rownames(top_expr)
    #Remove rows that have variance 0 (This will avoid the Na/Nan/Inf error in heatmap)
    ind = apply(top_expr, 1, var) == 0
    top_expr <- top_expr[!ind,]
    if(input$checkboxspia==TRUE){
      d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterbyspia,xaxis_font_size = 10,colors = colorRampPalette(brewer.pal(n = 9, input$hmpcolspia))(30),labRow = sym)}
    else{d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterbyspia,xaxis_font_size = 10,colors = colorRampPalette(rev(brewer.pal(n = 9, input$hmpcolspia)))(30),labRow = sym)}
  })
  
  #Create non-interactive SPIA heatmap function for download
  spiaheatmapalt <- reactive({
    dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    top_expr=heatmapfun(results=fileload(),expr=heatmapspia(),pval=spiagenes(),file = readexcel(),prj=input$projects,hmplim=input$hmplimspia,hmpsamp=input$hmpsamp2spia,
                        contrast=input$contrast)
    sym=rownames(top_expr)
    #Remove rows that have variance 0 (This will avoid the Na/Nan/Inf error in heatmap)
    ind = apply(top_expr, 1, var) == 0
    top_expr <- top_expr[!ind,]
    if(input$checkboxspia==TRUE){
      aheatmap(as.matrix(top_expr),distfun=dist2,scale="row",Rowv=TRUE,Colv=TRUE,fontsize = 10,color = colorRampPalette(brewer.pal(n = 9, input$hmpcolspia))(30),labRow = sym)}
    else{aheatmap(as.matrix(top_expr),distfun=dist2,scale="row",Rowv=TRUE,Colv=TRUE,fontsize = 10,color = colorRampPalette(rev(brewer.pal(n = 9, input$hmpcolspia)))(30),labRow = sym)}
  })
  
  #Download SPIA heatmap
  output$downloadspiaheatmap <- downloadHandler(
    filename = function(){
      paste0('SPIA_heatmap','.pdf',sep='')
    },
    content = function(file){
      pdf(file,width=9,height = 14,useDingbats=FALSE, onefile = F)
      spiaheatmapalt()
      dev.off()
    })
  ########################################################################################################################################################
  ########################################################################################################################################################
  ################################################################## REACTOME PA ANALYSIS################################################################
  ########################################################################################################################################################
  ########################################################################################################################################################
  #Get list of enriched pathways
  enrichpath = reactive({
    results=fileload()
    pd=pData(results$eset)
    org=unique(pd$organism)
    if(org %in% c("Mus musculus", "Mouse", "Mm","Mus_musculus","mouse")){
      org="mouse"
    }else{
      org="human"
    }
    deg= datasetInput0.5()
    deg=deg[abs(deg$fc) >2,]
    res <- enrichPathway(gene=deg$ENTREZID,pvalueCutoff=0.05, readable=T,organism=org)
  })
  
  #create different table to display
  enrichpath2 = reactive({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      res= enrichpath()
      res=as.data.frame(res)
      validate(
        need(nrow(res) > 0, "No results")
      )
      res = res %>% dplyr::select(-geneID)
    })
  })
  
  #get list of enriched pathways and display in table
  output$enrichpath = DT::renderDataTable({
    input$project
    input$contrast
    DT::datatable(enrichpath2(),
                  extensions = 'Buttons', options = list(
                    dom = 'Bfrtip',
                    buttons = list()),
                  rownames=FALSE,selection = list(mode = 'single', selected =1),escape=FALSE)
  })
  
  #Display list of genes in each enrichment pathway
  enrichgenes = reactive({
    res=enrichpath()
    validate(
      need(nrow(as.data.frame(res))>0,"No Enriched Pathways")
    )
    res=as.data.frame(res) 
    s = input$enrichpath_rows_selected
    genes = res[s, , drop=FALSE]
    genes = genes$geneID
    genes=gsub("/",", ",genes)
    return(genes)
  })
  
  #print genelist
  output$enrichgenes = renderPrint({
    enrichgenes()
  })
  
  #Create plot for visualizing enrichment results
  enrichplot = reactive({
    res= enrichpath()
    shiny::validate(
      need(nrow(as.data.frame(res))>0,"No Enriched Pathways")
    )
    if(input$enrichradio=='barplot'){
      barplot(res, showCategory = input$ncat)
    }else if(input$enrichradio=='dotplot'){
      dotplot(res,showCategory= input$ncat)
    }else if(input$enrichradio=='enrich'){
      emapplot(res)
    }
  })
  
  #Render the plot
  output$enrichplot <- renderPlot({
    enrichplot()
  })
  
  #Render the plot
  output$cnetplot <- renderPlot({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      res= enrichpath()
      validate(
        need(nrow(as.data.frame(res))>0,"No Enriched Pathways")
      )
      limmares= datasetInput0.5()
      genelist= limmares$fc
      names(genelist)=limmares$ENTREZID
      cnetplot(res, categorySize="pvalue", foldChange=genelist)
    })
  })
  
  ########################################################################################################################################################
  ########################################################################################################################################################
  ################################################################## REACTOME PA GSEA ################################################################
  ########################################################################################################################################################
  ########################################################################################################################################################
  #Get list of enriched pathways from GSEA
  gseapath = reactive({
    results=fileload()
    pd=pData(results$eset)
    org=unique(pd$organism)
    if(org %in% c("Mus musculus", "Mouse", "Mm","Mus_musculus","mouse")){
      org="mouse"
    }else{
      org="human"
    }
    limmares= datasetInput0.5()
    genelist= limmares$fc
    names(genelist)=limmares$ENTREZID
    genelist = sort(genelist, decreasing = TRUE)
    y <- gsePathway(genelist, nPerm=10000,pvalueCutoff=0.2,pAdjustMethod="BH", verbose=FALSE,organism=org)
  })
  
  #create different table to display
  gseapath2 = reactive({
    res= gseapath()
    res=as.data.frame(res) 
  })
  
  #Create Results table
  output$gseares = DT::renderDataTable({
    input$project
    input$contrast
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      DT::datatable(gseapath2(),
                    extensions = 'Buttons', options = list(
                      dom = 'Bfrtip',
                      buttons = list()),
                    rownames=FALSE,selection = list(mode = 'single', selected =1),escape=FALSE)
    })
  })
  
  #Render the plot emap
  output$plotemap <- renderPlot({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      res= gseapath()
      emapplot(res, color="pvalue")
    })
  })
  
  #Render the plot gsea
  output$plotgsea <- renderPlot({
    res= gseapath()
    gseares=gseapath2()
    s = input$gseares_rows_selected
    gseares = gseares[s, , drop=FALSE]
    id = gseares$ID
    gseaplot(res, geneSetID = id)
  })
  
  #Render the gsea pathway
  output$plotpath <- renderPlot({
    gseares=gseapath2()
    s = input$gseares_rows_selected
    gseares = gseares[s, , drop=FALSE]
    id = gseares$Description
    limmares= datasetInput0.5()
    limmares=limmares[is.na(limmares$ENTREZID)==F,]
    limmares=limmares[!duplicated(limmares$ENTREZID),]
    genelist= limmares$fc
    names(genelist)=limmares$ENTREZID
    genelist = sort(genelist, decreasing = TRUE)
    results=fileload()
    pd=pData(results$eset)
    org=unique(pd$organism)
    if(org %in% c("Mus musculus", "Mouse", "Mm","Mus_musculus","mouse")){
      org="mouse"
    }else{
      org="human"
    }
    viewPathway(id, readable=TRUE, foldChange=genelist, organism = org)
  })
  
  ######################################################################################################################################################
  ######################################################################################################################################################
  ################################################################ GAGE GENE ONTOLOGY ##################################################################
  ######################################################################################################################################################
  ######################################################################################################################################################
  
  #Run gage and get results
  datasetInput7 = reactive({
    final_res=datasetInput0.5() #get limma data
    logfc=final_res$fc #get FC values from limma data
    names(logfc)=final_res$ENTREZID # get entrez ids for each row
    results=fileload()
    pd=pData(results$eset)
    organism=pd$organism
    prjs=c("DS_FalcorFoxA2","YT_mir302","RJ_ESC_Laminin","RJ_CardiacHdac7_updated","DS_FalcorKO")
    prj2=c("DK_IPSC_lungepi","ZA_Boa_PKM2")
    if(!input$projects %in% prjs){
      if(!input$projects %in% prj2){
        validate(
          need(length(unique(organism))==1,"Please check pData file for errors in organism column. Does it have more than one organism or is it empty?")
        )
        organism=unique(pd$organism)[1]
      }}
    if(input$projects %in% prjs){
      organism="mouse"
    }
    else if(input$projects %in% prj2){
      organism="human"
    }
    if(organism=="human")
    {
      data(go.sets.hs) #load GO data from gage
      data(go.subs.hs)
      if(input$gage=='BP')
      {
        gobpsets = go.sets.hs[go.subs.hs$BP]
        go_res = gage(logfc, gsets=gobpsets)
      }
      else if(input$gage=='cc')
      {
        goccsets = go.sets.hs[go.subs.hs$CC]
        go_res = gage(logfc, gsets=goccsets, same.dir=TRUE)
      }
      else if(input$gage=='MF')
      {
        gomfsets = go.sets.hs[go.subs.hs$MF]
        go_res = gage(logfc, gsets=gomfsets, same.dir=TRUE)
      }}
    else if(organism=="Rat")
    {
      data(go.sets.rn) #load GO data from gage
      data(go.subs.rn)
      if(input$gage=='BP')
      {
        gobpsets = go.sets.rn[go.subs.rn$BP]
        go_res = gage(logfc, gsets=gobpsets)
      }
      else if(input$gage=='cc')
      {
        goccsets = go.sets.rn[go.subs.rn$CC]
        go_res = gage(logfc, gsets=goccsets, same.dir=TRUE)
      }
      else if(input$gage=='MF')
      {
        gomfsets = go.sets.rn[go.subs.rn$MF]
        go_res = gage(logfc, gsets=gomfsets, same.dir=TRUE)
      }
    }
    else 
    {
      data(go.sets.mm) #load GO data from gage
      data(go.subs.mm)
      
      if(input$gage=='BP')
      {
        gobpsets = go.sets.mm[go.subs.mm$BP]
        go_res = gage(logfc, gsets=gobpsets)
      }
      else if(input$gage=='cc')
      {
        goccsets = go.sets.mm[go.subs.mm$CC]
        go_res = gage(logfc, gsets=goccsets, same.dir=TRUE)
      }
      else if(input$gage=='MF')
      {
        gomfsets = go.sets.mm[go.subs.mm$MF]
        go_res = gage(logfc, gsets=gomfsets, same.dir=TRUE)
      }
    }
    return(go_res)
  })
  
  #Get all GO terms based on user-selection (upregulated/downregulated)
  datasetInput8 = reactive({
    go_res=datasetInput7()
    go_dd=input$go_dd
    if(go_dd=="upreg"){
      res=data.frame(go_res$greater)} #load limma data
    else if(go_dd=="downreg"){
      res=data.frame(go_res$less)
    }
    res = data.frame(GOterm=rownames(res),res)
    
    #Get GO id from GO terms
    row=data.frame(lapply(res,as.character),stringsAsFactors = FALSE)
    p=strsplit(row[,1], " ")
    m=sapply(p,"[",1)
    go_up=data.frame(GO_id=m,res)
    go_term=go_up$GO_id
    url= paste("http://amigo.geneontology.org/amigo/term/",go_term,sep = "") #create link to Gene Ontology Consortium
    go_up$link=paste0("<a href='",url,"'target='_blank'>","Link to Gene Ontology Consortium","</a>")
    go_up=as.data.frame(go_up)
    return(go_up)
  })
  
  #Print GO results in datatable
  output$table4 = DT::renderDataTable({
    input$go_dd
    input$gage
    input$radio
    input$project
    input$contrast
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      DT::datatable(datasetInput8(),
                    extensions = c('Buttons','Scroller'),
                    options = list(dom = 'Bfrtip',
                                   searchHighlight = TRUE,
                                   pageLength = 10,
                                   lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                   scrollX = TRUE,
                                   buttons = c('copy','print')
                    ),rownames=FALSE,escape=FALSE,selection = list(mode = 'single', selected =1))
    })
  })
  
  # Download function to get GO results in csv file
  output$downloadgo <- downloadHandler(
    filename = function() { paste('GO_',input$projects,'_',input$contrast,'_',input$gage,'_',input$go_dd,'.csv', sep='') },
    content = function(file) {
      write.csv(datasetInput8(), file)
    })
  
  ###################################################
  ###################################################
  ############## GET GENES FROM  GO #################
  ###################################################
  ###################################################
  #Text title for gene list table
  output$godesc <- renderText({
    s = input$table4_rows_selected
    dt = datasetInput8() #load GO data
    dt = dt[s, , drop=FALSE] #get GO data corresponding to selected row in table
    goid=dt$GO_id
    text=paste('Gene list for GO term :',goid,sep="")
    return(text)
  })
  
  # get GO associated genes
  GOHeatup = reactive({
    s = input$table4_rows_selected
    dt = datasetInput8() #load GO data
    dt = dt[s, , drop=FALSE] #get GO data corresponding to selected row in table
    results=fileload()
    pd=pData(results$eset)
    organism=pd$organism[1]
    prjs=c("DS_FalcorFoxA2","YT_mir302","RJ_ESC_Laminin","RJ_CardiacHdac7_updated","DS_FalcorKO")
    prj2=c("DK_IPSC_lungepi","ZA_Boa_PKM2")
    if(input$projects %in% prjs){
      organism="mouse"
    }
    else if(input$projects %in% prj2){
      organism="human"
    }
    goid=dt$GO_id
    if(organism=="human"){
      enterezid=paste("go.sets.hs$`",goid,"`",sep="")
    }
    else if(organism=="Rat"){
      enterezid=paste("go.sets.rn$`",goid,"`",sep="")
    }
    else{
      enterezid=paste("go.sets.mm$`",goid,"`",sep="")
    }
    entrezid=eval(parse(text=enterezid))
    limma=datasetInput0.5()
    lim_vals=limma[limma$ENTREZID %in% entrezid,]
  })
  
  #Print datatable with gene list
  output$x4 = DT::renderDataTable({
    input$gage
    input$go_dd
    input$radio
    input$project
    input$contrast
    goheatup=GOHeatup()
  },caption="Gene List",escape=FALSE)
  
  #Download function to get GO gene list as csv file
  output$downloadgogene <- downloadHandler(
    filename = function() { paste('GO_',input$projects,'_',input$contrast,'_',input$gage,'_',input$go_dd,'.csv', sep='') },
    content = function(file) {
      write.csv(GOHeatup(), file)
    })
  
  ###################################################
  ###################################################
  ########## MAKE HEATMAP WITH GO ###################
  ###################################################
  ###################################################
  #Set limit for number of genes that can be viewed in the heatmap
  output$hmplimgo <- renderUI({
    pval=GOHeatup()
    top_expr=datasetInput3()
    top_expr=top_expr[rownames(top_expr) %in% rownames(pval),]
    mx=nrow(top_expr)
    sliderInput("hmplimgo", label = h5("Select number of genes to view in the heatmap"), min = 2,max =mx, value = mx)
  })
  
  #Generate a heatmap color scale
  output$hmpscale_out3 = renderPlot({
    hmpscaletest(hmpcol=input$hmpcol3,voom=datasetInput3(),checkbox=input$checkbox3)
  })
  
  #plot heatmap
  goheatmapup <- reactive({
    dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    top_expr=datasetInput3() 
    pval=GOHeatup()
    top_expr=top_expr[rownames(top_expr) %in% rownames(pval),]#voom expression data of all genes corresponding to selected row in GO datatable
    top_expr=heatmapfun(results=fileload(),expr=as.data.frame(top_expr),pval=GOHeatup(),file = readexcel(),prj=input$projects,hmplim=input$hmplimgo,hmpsamp=input$hmpsamp3,
                        contrast=input$contrast)
    #Remove rows that have variance 0 (This will avoid the Na/Nan/Inf error in heatmap)
    ind = apply(top_expr, 1, var) == 0
    top_expr <- top_expr[!ind,]
    sym=rownames(top_expr)
    if(input$checkbox3==TRUE){
      d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterby3,xaxis_font_size = 10,colors = colorRampPalette(brewer.pal(n = 9, input$hmpcol3))(30),labRow = rownames(top_expr))}
    else{d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterby3,xaxis_font_size = 10,colors = colorRampPalette(rev(brewer.pal(n = 9, input$hmpcol3)))(30),labRow =rownames(top_expr))}
  })
  
  # render D3heatmap for GO genes
  output$goheatmap <- renderD3heatmap({
    input$hmpcol #user input-color palette
    input$clusterby #user input-cluster by
    input$checkbox #user input-reverse colors
    input$gene #user input-slider input for number of genes
    input$genelist
    input$makeheat
    input$gage
    input$go_dd
    input$table4_rows_selected
    input$tablecam_rows_selected
    input$projects
    input$contrast
    input$cameradd
    input$hmpsamp3
    input$hmplimgo
    goheatmapup()
  })
  
  #function for non-interactive heatmap for download
  goheatmapupalt <- reactive({
    dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    top_expr=datasetInput3() 
    top_expr=top_expr[rownames(top_expr) %in% rownames(pval),]#voom expression data of all genes corresponding to selected row in GO datatable
    top_expr=heatmapfun(results=fileload(),expr=top_expr,pval=GOHeatup(),file = readexcel(),prj=input$projects,hmplim=input$hmplimgo,hmpsamp=input$hmpsamp3,
                        contrast=input$contrast)
    
    #Remove rows that have variance 0 (This will avoid the Na/Nan/Inf error in heatmap)
    ind = apply(top_expr, 1, var) == 0
    top_expr <- top_expr[!ind,]
    if(input$checkbox3==TRUE){
      aheatmap(as.matrix(top_expr),distfun=dist2,scale="row",Rowv=TRUE,Colv =TRUE,fontsize = 10,color = colorRampPalette(brewer.pal(n = 9, input$hmpcol3))(30),labRow = rownames(top_expr))}
    else{aheatmap(as.matrix(top_expr),distfun=dist2,scale="row",Rowv=TRUE,Colv = TRUE,fontsize = 10,color = colorRampPalette(rev(brewer.pal(n = 9, input$hmpcol3)))(30),labRow = rownames(top_expr))}
  })
  
  #Download GO heatmap
  output$downloadgoheatmap <- downloadHandler(
    filename = function(){
      paste0('GO_heatmap','.pdf',sep='')
    },
    content = function(file){
      pdf(file,width=9,height = 14,useDingbats=FALSE, onefile = F)
      goheatmapupalt()
      dev.off()
    })
  
  #########################################################################################################################################################
  #########################################################################################################################################################
  ########################################################## CREATE HEATMAP FOR LIMMA DATA#################################################################
  #########################################################################################################################################################
  #########################################################################################################################################################
  #Text title for type of heatmap being displayed in the heatmap tab
  output$htitle <- renderText({
    hmip=input$hmip
    if(input$hmip=="genenum"){text="Heatmap of Top Genes "}
    else if(input$hmip=="geneli"){text="Heatmap of Genelist "}
    else if(input$hmip=="vargenes"){text="Heatmap of top n variable genes "}
  })
  
  #manually create scale (colorkey) for heatmap
  output$hmpscale_out = renderPlot({
    hmpscaletest(hmpcol=input$hmpcol,voom=datasetInput3(),checkbox=input$checkbox)
  })
  
  ###################################################
  ###################################################
  #################### TOP GENES ####################
  ###################################################
  ###################################################
  output$dropdown <- renderUI({
    radio=input$radio
    if(radio=="none"){
      selectInput("sortby", "Sort By",c('FDR'="sortnone",'Absolute Fold Change' = "sortab",'Positive Fold Change' = "sortpos",'Negative Fold Change' = "sortneg"))
    }
    else if(radio=="up"){
      selectInput("sortby", "Sort By",c('FDR'="sortnone",'Fold Change' = "sortab"))
    }
    else if(radio=="down"){
      selectInput("sortby", "Sort By",c('FDR'="sortnone",'Fold Change' = "sortab"))
    }
    else if(radio=="both"){
      selectInput("sortby", "Sort By",c('FDR'="sortnone",'Absolute Fold Change' = "sortab",'Positive Fold Change' = "sortpos",'Negative Fold Change' = "sortneg"))
    }
    
  })
  
  #create heatmap function for top number of genes as chosen from the slider
  datasetInput4 <- reactive({
    validate(
      need(input$gene, "Please Enter number of genes to plot heatmap ")
    )
    #sort by pval
    n<-input$gene #number of genes selected by user (input from slider)
    d<-datasetInput()
    sortby=input$sortby
    if(sortby=='sortnone'){
      res<-d[order(d$adj.P.Val),]
    }else if(sortby=='sortab'){
      res<-d[order(-abs(d$fc)),]
    }else if(sortby=='sortpos'){
      res<-d[order(-d$fc),]
    }else if(sortby=='sortneg'){
      res<-d[order(d$fc),]
    }
    if(n>nrow(d)){
      reqd_res=res[1:nrow(d),]} #get top n number of genes
    else{
      reqd_res=res[1:n,]
    }
    return(reqd_res)
  })
  
  
  #create heatmap function for top n genes
  heatmap <- reactive({
    dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    pval=datasetInput4()
    top_expr= createheatmap(results=fileload(),expr=datasetInput3(),pval=pval,hmpsamp=input$hmpsamp,contrast=input$contrast)
    top_expr=as.data.frame(top_expr)
    col=colnames(top_expr)
    top_expr$ENSEMBL=rownames(top_expr)
    top_expr=inner_join(top_expr,pval,by="ENSEMBL")
    rownames(top_expr)=top_expr$SYMBOL
    top_expr=top_expr %>% dplyr::select(col)
    validate(
      need(nrow(top_expr) > 1, "No results")
    )
    if(input$checkbox==TRUE){
      d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterby,xaxis_font_size = 10,colors = colorRampPalette(brewer.pal(n = 9, input$hmpcol))(30))}
    else{d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterby,xaxis_font_size = 10,colors = colorRampPalette(rev(brewer.pal(n = 9, input$hmpcol)))(30))}
  })
  
  #alternate hearmap function for download
  heatmapalt <- reactive({
    dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    pval=datasetInput4()
    top_expr= createheatmap(results=fileload(),expr=datasetInput3(),pval=pval,hmpsamp=input$hmpsamp,contrast=input$contrast)
    sym=pval$SYMBOL
    validate(
      need(nrow(top_expr) > 1, "No results")
    )
    if(input$checkbox==TRUE){
      aheatmap(as.matrix(top_expr),distfun=dist2,scale="row",Rowv = TRUE,Colv = TRUE,fontsize = 10,color = colorRampPalette(brewer.pal(n = 9, input$hmpcol))(30),labRow = sym)}
    else{aheatmap(as.matrix(top_expr),distfun=dist2,scale="row",Rowv = TRUE,Colv = TRUE,fontsize = 10,color = colorRampPalette(rev(brewer.pal(n = 9, input$hmpcol)))(30),labRow = sym)}
  })
  
  ###################################################
  ###################################################
  ####### ENTER GENELIST ############################
  ###################################################
  ###################################################
  # Get gene list from user, annotate to ENSEMBL id and get their expression values
  datasetInput41 = reactive({
    file=input$genelistfile
    genes=read.table(file=file$datapath, stringsAsFactors = F) #get complete gene list as string
    df=as.vector(genes$V1)
    df=tolower(df)
    firstup <- function(x) {
      substr(x, 1, 1) <- toupper(substr(x, 1, 1))
      x
    }
    genelist=firstup(df)
    results=fileload()
    #     pd=pData(results$eset)
    #     org=unique(pd$organism)
    #     
    #     if(org=="human"){
    #       dataset="hsapiens_gene_ensembl"
    #     }
    #     else if(org=="Rat"){
    #       dataset="rnorvegicus_gene_ensembl"
    #     }
    #     else{
    #       dataset="mmusculus_gene_ensembl"
    #     }
    #     ensembl = useEnsembl(biomart="ensembl", dataset=dataset)
    #load limma and voom data
    limma=datasetInput()
    voom=datasetInput3()
    #get expression values of the genes in the gene list
    # user-defined identifier for the gene list
    if(input$selectidentifier=='ensembl')
    {
      sym=limma[limma$ENSEMBL %in% genelist,] 
      sym= sym %>% dplyr::select(ENSEMBL,SYMBOL)
      #       genes <- getBM(attributes=c('ensembl_gene_id','external_gene_name'), filters ='ensembl_gene_id', values =df, mart = ensembl)
      #       genelist=genes$ensembl_gene_id
    }
    else if(input$selectidentifier=='entrez')
    {
      sym=limma[limma$ENTREZID %in% genelist,] 
      sym= sym %>% dplyr::select(ENSEMBL,SYMBOL)
      #       genes <- getBM(attributes=c('ensembl_gene_id','entrezgene'), filters ='entrezgene', values =df, mart = ensembl)
      #       genelist=genes$ensembl_gene_id
    }
    else if(input$selectidentifier=='genesym')
    {
      sym=limma[limma$SYMBOL %in% genelist,] 
      sym= sym %>% dplyr::select(ENSEMBL,SYMBOL)
      #       genes <- getBM(attributes=c('ensembl_gene_id','external_gene_name'), filters ='external_gene_name', values =df, mart = ensembl)
      #       genelist=genes$ensembl_gene_id
    }
    expr_vals=merge(voom,sym,by="row.names")
    rownames(expr_vals)=expr_vals$SYMBOL
    expr_vals = expr_vals %>% dplyr::select(-Row.names,-SYMBOL,-ENSEMBL)
    #expr_vals=data.frame(expr_vals[,-c(1,(ncol(expr_vals)-1))])
    validate(
      need(nrow(expr_vals) > 1, "Please Check Identifier chosen or Select genelist from Raw Expression Data tab")
    )
    return(expr_vals)
  })
  
  #create heatmap function for gene-list given by user
  heatmap2 = function(){
    dist2 = function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    limma=datasetInput()
    expr = datasetInput41()
    #expr=data.frame(expr[,-ncol(expr)])
    #     genelist= rownames(expr)
    #     sym=limma[limma$ENSEMBL %in% genelist,] %>% dplyr::select(SYMBOL)
    expr2= createheatmap(results=fileload(),expr=expr,hmpsamp=input$hmpsamp,contrast=input$contrast)
    validate(
      need(nrow(expr2)>1, "No results")
    )
    if(input$checkbox==TRUE){
      d3heatmap(as.matrix(expr2),distfun=dist2,scale="row",dendrogram=input$clusterby,xaxis_font_size = 10,colors = colorRampPalette(brewer.pal(n = 9, input$hmpcol))(30))}
    else{d3heatmap(as.matrix(expr2),distfun=dist2,scale="row",dendrogram=input$clusterby,xaxis_font_size = 10,colors = colorRampPalette(rev(brewer.pal(n = 9, input$hmpcol)))(30))}
  }
  
  heatmap2alt = function(){
    dist2 = function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    expr = datasetInput41()
    #expr2=data.frame(expr[,-ncol(expr)])
    top_expr= createheatmap(results=fileload(),expr=expr2,hmpsamp=input$hmpsamp,contrast=input$contrast)
    if(input$checkbox==TRUE){
      aheatmap(as.matrix(expr2),distfun=dist2,scale="row",Rowv=TRUE,Colv=TRUE,fontsize = 10,color = colorRampPalette(brewer.pal(n = 9, input$hmpcol))(30))}
    else{aheatmap(as.matrix(expr2),distfun=dist2,scale="row",Rowv=TRUE,Colv=TRUE,fontsize = 10,color = colorRampPalette(rev(brewer.pal(n = 9, input$hmpcol)))(30))}
  }
  
  ###################################################
  ###################################################
  ####### TOP VARIABLE GENES  #######################
  ###################################################
  ###################################################
  #Extract top n (user-selected) variable genes
  var.genes = reactive({
    n=as.numeric(input$vgene)
    results=fileload()
    v = results$eset
    keepGenes <- v@featureData@data
    #keepGenes <- v@featureData@data %>% filter(!(seq_name %in% c('X','Y')) & !(is.na(SYMBOL)))
    pData<-phenoData(v)
    v.filter = v[rownames(v@assayData$exprs) %in% rownames(keepGenes),]
    Pvars <- apply(v.filter@assayData$exprs,1,var)
    select <- order(Pvars, decreasing = TRUE)[seq_len(min(n,length(Pvars)))]
    v.var <-v.filter[select,]
    m<-v.var@assayData$exprs
    rownames(m) <- v.var@featureData@data$SYMBOL
    m=as.data.frame(m)
    m=unique(m)
    return(m)
  })
  
  #D3 heatmap for top n variable genes
  varheatmap <- reactive({
    dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    top_expr= createheatmap(results=fileload(),expr=var.genes(),hmpsamp=input$hmpsamp,contrast=input$contrast)
    validate(
      need(nrow(top_expr) > 1, "No results")
    )
    if(input$checkbox==TRUE){
      d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterby,xaxis_font_size = 10,colors = colorRampPalette(brewer.pal(n = 9, input$hmpcol))(30))}
    else{d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterby,xaxis_font_size = 10,colors = colorRampPalette(rev(brewer.pal(n = 9, input$hmpcol)))(30))}
  })
  
  # Alternate function to download non-interactive heatmap of top n variable genes
  varheatmapalt <- reactive({
    dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    top_expr= createheatmap(results=fileload(),expr=var_genes(),hmpsamp=input$hmpsamp,contrast=input$contrast)
    validate(
      need(nrow(top_expr) > 1, "No results")
    )
    if(input$checkbox==TRUE){
      aheatmap(as.matrix(top_expr),distfun=dist2,scale="row",Rowv = TRUE,Colv = TRUE,fontsize = 10,color = colorRampPalette(brewer.pal(n = 9, input$hmpcol))(30))}
    else{aheatmap(as.matrix(top_expr),distfun=dist2,scale="row",Rowv = TRUE,Colv = TRUE,fontsize = 10,color = colorRampPalette(rev(brewer.pal(n = 9, input$hmpcol)))(30))}
  })
  
  
  # Render d3 heatmap function 
  output$heatmap <- renderD3heatmap({
    input$hmpcol #user input-color palette
    input$clusterby #user input-cluster by
    input$checkbox #user input-reverse colors
    input$gene #user input-slider input for number of genes
    input$genelist
    input$hmip
    input$makeheat
    input$gage
    input$go_dd
    input$ga
    input$table4_rows_selected
    input$tablecam_rows_selected
    input$radio
    input$projects
    input$contrast
    input$cameradd
    input$hmpsamp
    input$hmplim
    input$lfc
    input$apval
    input$sortby
    input$vgene
    #if user selected enter n num of genes, call heatmap() and if user entered genelist, call heatmap2()
    isolate({
      if(input$hmip == 'genenum'){heatmap()}
      else if(input$hmip == 'geneli'){heatmap2()}
      else if(input$hmip == 'vargenes' ){varheatmap()}
    })
  })
  
  #Download function for heatmaps 
  output$downloadheatmap <- downloadHandler(
    filename = function(){
      paste0('heatmap','.pdf',sep='')
    },
    content = function(file){
      pdf(file,width=9,height =14,useDingbats=FALSE, onefile = F)
      if(input$hmip == 'genenum'){heatmapalt()}
      else if(input$hmip == 'geneli'){heatmap2alt()}
      else if(input$hmip == 'vargenes' ){varheatmapalt()}
      dev.off()
    })
  
}#end of server
