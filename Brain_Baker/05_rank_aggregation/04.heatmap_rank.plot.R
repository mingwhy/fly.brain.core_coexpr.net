
#https://stackoverflow.com/questions/23938889/how-to-create-discrete-legend-in-pheatmap
changeLegend<-function(breaks, color) {
  tree <- grid.ls(viewport=T, print=F)
  #find legend
  legendvp <- tail( grep("GRID.VP", tree$name), 1)
  #get rid of grobs in this viewport
  drop <- tree$name[grepl(tree$vpPath[legendvp],tree$vpPath) & 
                      grepl("grob",tree$type)]
  sapply(drop, grid.remove)
  
  #calculate size/position of labels
  legend_pos = seq(0,to=1,length.out=length(breaks))
  brat = seq(0,to=1,length.out=length(breaks))
  h = 1/(length(breaks)-1)
  
  #render legend
  seekViewport(tree$name[legendvp])    
  grid.rect(x = 0, y = brat[-length(brat)], 
            width = unit(10, "bigpts"), height = h, hjust = 0, 
            vjust = 0, gp = gpar(fill = color, col = "#FFFFFF00"))
  grid.text(breaks, x = unit(12, "bigpts"), 
            y = legend_pos, hjust = 0,)
}

#####################################################
## heatmap
options(stringsAsFactors = F)
library(igraph)
library(Matrix);
library(ggplot2);library(gridExtra)
library(ggpubr)
library(tidyverse)
library(RColorBrewer)
library(ComplexHeatmap)

file='cor.rank.list.rds'
cor.rank.list=readRDS(file); #read in list is fast
length(cor.rank.list) #sc 67
sapply(cor.rank.list,length); #ranked edges in each cluster


# edge presence/absence across cell clusters
edge.list=readRDS('top_0.05_edge_in.cell.types.rds');
length(edge.list); #numebr of cell clusters  
sapply(edge.list,dim) #all have the same row


tmp=edge.list[[1]]
edges <- do.call(paste, c(tmp[ , 1:2 ], sep=";")) #much faster
#edges <- apply( tmp[ , 1:2 ] , 1 , paste , collapse = ";" )

values=lapply(edge.list, '[',,3) #extract the 1st column in each list matrix member
val.mat=Reduce(`cbind`,values)
dim(val.mat)
rownames(val.mat)=edges
colnames(val.mat)=names(edge.list)

commonality=apply(val.mat,1,sum)
max(commonality) #sc 42
table(commonality) 
length(commonality) # 1509453 number of edges per network

#edge.i=43;sample.nedge.for.this.group=100;
#edge.i=35;sample.nedge.for.this.group=100; #sn
edge.i=10;sample.nedge.for.this.group=100; #cocaine brain
  #edge.group=val.mat[commonality==edge.i,,drop=F]
  edge.group=val.mat[commonality>=edge.i,,drop=F]
  dim(edge.group) #323  47; 163 edges, 0/1 binary matrix, n.edge X n.cell.cluster
  
  # sample 100 edges for sc
  if(T){
    set.seed(220107)
    edge.group=edge.group[sample(1:nrow(edge.group),sample.nedge.for.this.group,replace=F),,drop=F]
    dim(edge.group)
  }
  
  # replace 0/1 with rank or quantile number
  rownames(edge.group) #edge names: gene1;gene2
  head(cor.rank.list[[1]]) #the same name convention
  out=lapply(1:length(cor.rank.list),function(i){
    x=cor.rank.list[[i]]
    y=((1:length(x))-1)/(length(x)-1);
    names(y)=x;
    y[rownames(edge.group)]
  })
  df.out=as.data.frame(Reduce(`cbind`,out))
  dim(df.out) # 323 edges x 47 cell clusters
  colnames(df.out)=colnames(edge.group)
  rownames(df.out)=rownames(edge.group)
  
  #my.breaks <- seq(0, 1, by=0.2)
  #my.breaks <- c(0,0.005,0.01,0.02,0.05,0.1,1)
  my.breaks <- c(0,0.001,0.002,0.005,0.01,0.02,0.05,0.1,1)
  #my.colors <- colorRampPalette(colors = c("blue4", "white"))(length(my.breaks)-1)
  my.colors <-viridis::magma(length(my.breaks)-1)
  #my.colors<-c('red','blue','yellow')
  #my.colors[length(my.breaks)-1]='white'
  my.colors[length(my.breaks)]='white'
  
  #pdf(paste0('./cocaine_top0.05_edge_',edge.i,'.pdf'),width=30,height = 16)
  pdf(paste0('./sn_top0.05_edge_',edge.i,'.pdf'),width=14,height = 8)
  pheatmap(t(df.out),
           color = my.colors,
           breaks = my.breaks,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean", 
           clustering_method='ward.D2',
           legend=T,
           fontsize_row = 10, fontsize_col=6);
           #treeheight_row = 0, treeheight_col = 0 )
  changeLegend(my.breaks, my.colors)
  
  dev.off()
  
  
  #https://stackoverflow.com/questions/23938889/how-to-create-discrete-legend-in-pheatmap
  changeLegend<-function(breaks, color) {
    tree <- grid.ls(viewport=T, print=F)
    #find legend
    legendvp <- tail( grep("GRID.VP", tree$name), 1)
    #get rid of grobs in this viewport
    drop <- tree$name[grepl(tree$vpPath[legendvp],tree$vpPath) & 
                        grepl("grob",tree$type)]
    sapply(drop, grid.remove)
    
    #calculate size/position of labels
    legend_pos = seq(0,to=1,length.out=length(breaks))
    brat = seq(0,to=1,length.out=length(breaks))
    h = 1/(length(breaks)-1)
    
    #render legend
    seekViewport(tree$name[legendvp])    
    grid.rect(x = 0, y = brat[-length(brat)], 
              width = unit(10, "bigpts"), height = h, hjust = 0, 
              vjust = 0, gp = gpar(fill = color, col = "#FFFFFF00"))
    grid.text(breaks, x = unit(12, "bigpts"), 
              y = legend_pos, hjust = 0,)
  }
  
  
if(F){
  gplots::heatmap.2(t(df.out),
                    #scale='column',
                    scale='none',
                    #cellnote = my.df,  # same data set for cell labels
                    notecol="black",      # change font color of cell labels to black
                    main = paste(nrow(edge.group),'edges at edge',i), # heat map title
                    density.info="none",  # turns off density plot inside color legend
                    trace="none",         # turns off trace lines inside the heat map
                    margins =c(20,20),  # ("margin.Y", "margin.X"), give space for names
                    #col=my_palette,       # use on color palette defined earlier
                    col=c('white','blue'),
                    #dendrogram="row",     # only draw a row dendrogram
                    dendrogram="both",
                    Rowv=T,Colv=T)  
}


