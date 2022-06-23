
# compute cell type-specific network summary statistics
options(stringsAsFactors = F)
library(igraph)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(RColorBrewer)

## generate network summary statistics
#topK=0.05; #sc
#(files=Sys.glob('./brain_scRNA-seq_n15c0.005_bigScale2/top_*_networks_cutoff_cell.type.rds'))
#path='./brain_scRNA-seq_n15c0.005_bigScale2/net_summary/';

topK=0.1; #sn
(files=Sys.glob('./brain_snRNA-seq_n15c0.005_bigScale2/top_*_networks_cutoff_cell.type.rds'))
path='./brain_snRNA-seq_n15c0.005_bigScale2/net_summary/';

if(!dir.exists(path)){dir.create(path)}

(files=files[grep(topK,files)])

for(file in files){
  (cutoff=strsplit(basename(file),'\\_')[[1]][[2]])
  (outfile=paste0(path,'/top_',cutoff,'_per.gene_net.summary.rds'));
  (outfile2=paste0(path,'/top_',cutoff,'_net.summary.txt'));
  
  if(!file.exists(outfile)){
    keep.list=readRDS(file);
    names(keep.list);
    network.stat=list();
    
    for(cell.type in names(keep.list)){
      mat=keep.list[[cell.type]];
      diag(mat)<-0;
      i=Matrix::rowSums(mat)
      rm.i=which(i==0)
      if(length(rm.i)!=0){mat=mat[-rm.i,-rm.i]}
      mat[lower.tri(mat)]=0;
      
      G <- graph_from_adjacency_matrix(mat,mode='undirected')
      #is.weighted(G)
      #is.directed(G)
      
      ############ begin compute summary stats ##################
      (ngene<-length(V(G)))
      (nedge<-gsize(G))
      
      ## mean_distance(G)
      (ave.path.length<-average.path.length(G,directed=F))
      
      ## the clustering coefficient
      (clustering.coeff=transitivity(G,type='global')) #Transitivity measures the probability that the adjacent vertices of a vertex are connected. This is sometimes also called the clustering coefficient.
      # per gene, clustering
      Cluster.coeff.local=transitivity(G,type='local')
      Cluster.coeff.barrat=transitivity(G, type="barrat")
     
      ## per gene, betweeness 
      Betweenness=betweenness(graph=G,directed=F,weights=NULL,normalized = TRUE)
      
      ## per gene, degree and pagerank
      Degree=igraph::degree(graph = G)
      PAGErank=igraph::page_rank(graph=G,directed = F)$vector
      
      #Connected components of a graph
      clu <- components(G)
      (largest.connected.comp<-max(clu$csize))
      
      ## modularity
      wtc <- cluster_walktrap(G)
      (module<-modularity(G, membership(wtc)))
      
      ## density, The density of a graph is the ratio of the number of edges and the number of possible edges.
      (density<-edge_density(G, loops=F) )
    
      out<-data.frame(cell.type,gene=names(Degree),
                      ngene, nedge, clustering.coeff, ave.path.length,
                      Degree,Cluster.coeff.local,Cluster.coeff.barrat,Betweenness,PAGErank,
                      largest.connected.comp, module,density);
     
      network.stat[[cell.type]]<-out;
    }
    saveRDS(network.stat,outfile);
    
    #network.stat=readRDS(outfile);
    names(network.stat); #33 cell types
    colnames(network.stat[[1]])
    net.summary=t(sapply(network.stat,function(x){
      x1=x[1,c(1,3:6,12:14)]
      x2=apply(x[,7:11],2,function(i) mean(i,na.rm=T)) #these stats were computed on a per-gene basis
      c(x1,x2)
    }))
    dim(net.summary) #33 cell types X 13 stats
    #write.table(net.summary,file=outfile2,quote=F)
    data.table::fwrite(net.summary,file=outfile2,quote=F)
  }else{
    
    network.stat=readRDS(outfile);
    net.summary=read.table(outfile2);
  } 
    
}

##############################################################3
## use top_0.05 as main result in the manuscript and make plot
(outfile=paste0(path,'/top_0.05_per.gene_net.summary.rds'));
(outfile2=paste0(path,'/top_0.05_net.summary.txt'));
df=data.table::fread(outfile2)
head(df)
p1=ggplot(df,aes(x=ngene,y=clustering.coeff))+geom_point()+theme_bw()+
  theme(axis.title = element_text(size=20))
p2=ggplot(df,aes(x=module,y=clustering.coeff))+geom_point()+theme_bw()+
  theme(axis.title = element_text(size=20))+xlab('modularity')
grid.arrange(p1,p2,ncol=2)

summary(df$largest.connected.comp/df$ngene)

df=df[order(df$ngene,decreasing = T),]
rownames(df)<-NULL
x=df[,c(1,2,3,4,5,6)]
colnames(x)[c(4,5,6)]=c('clustering.coefficient','average.path.length','largest.connected.component')
write.table(x,quote=F,paste0(path,'/top_0.05_net.summary2.txt'))

if(F){
pdf(paste0(path,'top_0.05_net.summary.pdf'),height = 16,width = 16)
#grid.table(df[,c(1,2,3,4,5,6)])
grid.table(x,rows=NULL)
dev.off()
}
