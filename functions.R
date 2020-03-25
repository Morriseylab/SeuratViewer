require(dplyr)
library(slingshot)

getGeneRange <- function(scrna,gene_probes){
  gene_values=as.data.frame(FetchData(scrna,gene_probes[1]))
  minr<- round(min(gene_values),2) 
  maxr<- round(max(gene_values),2)
  
  return(c(ifelse(minr==0,.1,minr-.1),maxr))
  #return(c(1,12))
  
}

bigene_getValues <- function(scrna,gene_probes,limita,limitb){
  gene_values=FetchData(scrna,c(gene_probes[1],gene_probes[2]))
  colnames(gene_values) <- c('genea','geneb')
  # gene_values=as.data.frame(gene_values) %>% 
  #   mutate(value =ifelse(genea>=limita[1] & geneb <limitb[1], gene_probes[1], ifelse(genea<limita[1] & geneb >=limitb[1],
  #                   gene_probes[2],ifelse(genea>=limita[1] & geneb >=limitb[1],"DoublePos","NULL"))))
  as.data.frame(gene_values) %>% 
    mutate(value = ifelse(genea>=limita[1] & geneb>=limitb[1],
                          'both',
                          ifelse(genea>=limita[1] & geneb<limitb[1],
                                 gene_probes[1],
                                 ifelse(genea<=limita[1] & geneb>=limitb[1],
                                        gene_probes[2],
                                        'none')))
    ) #%>% select(value)
}

monocle_theme_opts <- function()
{
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}

bigene_plot <- function (scrna, gene_probes, x=1,y=2, limita=c(1,100), limitb=c(1,100), marker_size = 0.1,
                         title = NULL,type="tsne")
{
  
  gene_values <- bigene_getValues(scrna,gene_probes,limita,limitb)
  projection=as.data.frame(eval(parse(text=paste("Embeddings(scrna, reduction =\"",type,"\")",sep=""))))
  colnames(projection) <- c("Component.1", "Component.2")
  proj_gene <- data.frame(cbind(projection, gene_values))
  #proj_gene$value = factor(proj_gene$value,levels=unique(proj_gene$value))
  proj_gene$value = factor(proj_gene$value,levels=c('both',gene_probes[1],gene_probes[2],'none'))
  proj_gene <- arrange(proj_gene, desc(value))
  
  p <- ggplot(proj_gene, aes(Component.1, Component.2)) +
    geom_point(aes(colour = value), size = marker_size) +
    scale_color_manual(values=c("#E41A1C","#377EB8","#4DAF4A", 'grey75'),drop=F) +
    theme(legend.key.size = unit(10,"point")) + xlab(paste("Component", x)) +
    ylab(paste("Component", y))
  
  
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  p <- p + monocle_theme_opts() + theme(plot.title = element_text(hjust = 0.5),
                                        legend.position="bottom",
                                        legend.title=element_blank(),
                                        legend.text=element_text(size=14),
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank())
  
  return(p)
}


#given scrna data, variable/ ident and project name, find all lig rec pairs
ligrec <- function(scrna,pair,prj,perc,filetype){
  #get grouping variable
  var=as.character(pair)
  tt=rownames(GetAssayData(object = scrna, slot = "counts"))
  #Read ligrec file based on organism
  file = read.csv("data/param.csv")
  if (filetype=="list"){
    org=as.character(file$organism[file$projects==prj])}
  else if(filetype == "upload"){
    org = unique(scrna@meta.data$org)
  }
  #genes=fread("data/ligrecgenes.txt",header = TRUE)
  if(org=="mouse"){rl=read.csv("data/Mm_PairsLigRec.csv")}else if(org=="human"){rl=read.csv("data/Hs_PairsLigRec.csv")}
  genes=unique(c(as.character(rl$ligand),as.character(rl$receptor)))
  # if(org=="human"){
  #   genes$genes=toupper(genes$genes)
  # }
  genes2=tt[tt %in% genes]
  #For all unique genes in the ligrec list, get their expression value for all cells and the groups the cells belong to
  my.data=FetchData(scrna,c(var,genes2))
  colnames(my.data)[1]= "clust"
  #my.data$clust=factor(my.data$clust,levels=unique(my.data$clust))
  perc=perc/100
  #if(org=="mouse"){rl=read.csv("data/Mm_PairsLigRec.csv")}else if(org=="human"){rl=read.csv("data/Hs_PairsLigRec.csv")}
  result=data.frame()
  res=data.frame()
  #loop over each cluster to find pairs
  for(i in 1:(length(levels(my.data$clust)))){
    for(j in 1:(length(levels(my.data$clust)))){
      # for(i in levels(my.data$clust)){
      #   for(j in levels(my.data$clust)){
      #if(i!=j){
        #from the large martix, subselect receptor and lig subgoups (if i=1 and j=2, keep cells in grps 1 and 2)
        test=my.data[my.data$clust==levels(my.data$clust)[i] | my.data$clust==levels(my.data$clust)[j],]
        #Subselect genes in receptor list in cells in rec subgroup (say 1)
        R_c1=test[test$clust==levels(my.data$clust)[i] ,(colnames(test) %in% rl$receptor)]
        #Subselect genes in ligand list in cells in lig subgroup (say 2)
        L_c2=test[test$clust==levels(my.data$clust)[j] , (colnames(test) %in% rl$ligand)]
        if(nrow(R_c1)!=0 &nrow(L_c2)!=0){
          #keep genes that are expressed in more than user-input percent of the cells
          keep1 = colSums(R_c1>0)>=perc*dim(R_c1)[1]
          keep2 = colSums(L_c2>0)>=perc*dim(L_c2)[1]
          R_c1=R_c1[,keep1]
          L_c2=L_c2[,keep2]
          #get list of lig-rec pairs
          res=rl[(rl$ligand %in% colnames(L_c2)) & (rl$receptor %in% colnames(R_c1)),]
        }else{}
        
      # }
      # else{}
      if(nrow(res)!=0){
        res$Receptor_cluster=levels(my.data$clust)[i]
        res$Lig_cluster=levels(my.data$clust)[j]
        result=rbind(result,res)
      }else{result=result}
    }
  }
  # get final list of all lig-rec pairs
  #result=result[result$Receptor_cluster!=result$Lig_cluster,]
  return(result)
}




#function to curve graph edges in igraph
autocurve.edges2 <-function (graph, start = 0.5)
{
  cm <- count.multiple(graph)
  mut <-is.mutual(graph)  #are connections mutual?
  el <- apply(get.edgelist(graph, names = FALSE), 1, paste,
              collapse = ":")
  ord <- order(el)
  res <- numeric(length(ord))
  p <- 1
  while (p <= length(res)) {
    m <- cm[ord[p]]
    mut.obs <-mut[ord[p]] #are the connections mutual for this point?
    idx <- p:(p + m - 1)
    if (m == 1 & mut.obs==FALSE) { #no mutual conn = no curve
      r <- 0
    }
    else {
      r <- seq(-start, start, length = m)
    }
    res[ord[idx]] <- r
    p <- p + m
  }
  res
}

#Get maximum dimensions run
getMaxDim <- function(object){
  
  object@reductions$pca@jackstraw$overall.p.values %>% 
    as.data.frame(.) %>% 
    mutate(adj = p.adjust(Score,method='bonferroni')) %>% 
    filter(adj <0.05) %>% 
    summarise(max=max(PC)) %>% 
    pull(max)
  
  
}

#Seurat Extras functions
CurvePlot = function(object,
                     sds=NULL,
                     group.by = NULL,
                     reduction = 'umap',
                     dims = 1:2,
                     cols=NULL,
                     label=T
) {
  object[['ident']] <- Idents(object = object)
  group.by <- group.by %||% 'ident'
  dims <- paste0(Key(object = object[[reduction]]), dims)
  
  curved <-
    bind_rows(lapply(names(slingCurves(sds)), function(x) {
      c <- slingCurves(sds)[[x]]
      d <- as.data.frame(c$s[c$ord, dims])
      d$curve <- x
      return(d)
    }))
  
  
  DimPlot(object,cols=cols,label = label,group.by = group.by,reduction = reduction) +
    geom_path(aes_string(dims[1], dims[2], linetype = "curve"), curved, size =1)
}

`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}