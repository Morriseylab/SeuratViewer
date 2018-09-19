require(dplyr)


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
                         title = NULL)
{
  
  gene_values <- bigene_getValues(scrna,gene_probes,limita,limitb)
  projection=as.data.frame(scrna@dr$tsne@cell.embeddings)
  colnames(projection) <- c("Component.1", "Component.2")
  proj_gene <- data.frame(cbind(projection, gene_values))
  #proj_gene$value = factor(proj_gene$value,levels=unique(proj_gene$value))
  proj_gene$value = factor(proj_gene$value,levels=c('both',gene_probes[1],gene_probes[2],'none'))
  proj_gene <- arrange(proj_gene, desc(value))
  
  p <- ggplot(proj_gene, aes(Component.1, Component.2)) +
    geom_point(aes(colour = value), size = marker_size) +
    scale_color_manual(values=c("#E41A1C","#377EB8","#4DAF4A", 'grey90'),drop=F) +
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