
# compaute rank aggregated score for each edge 
# using all cell types

library(Matrix);library(tidyverse)
library(ggplot2);library(gridExtra)
library(RColorBrewer); #library(RankAggreg) #too slow
library(RobustRankAggreg)

nets=readRDS('~/Documents/fly.brain_coexpr.net_2022.06/Brain_Davie/brain_scRNA-seq_n15c0.005_bigScale2/top_0.05_networks_cutoff_cell.type.rds')
length(nets) #67
names(nets)

## cell type tree
if(F){
  dist.mat=matrix(0,ncol=length(nets),nrow=length(nets))
  for(i in 1:(nrow(dist.mat)-1)){
    for(j in ((i+1):nrow(dist.mat))){
      net1=nets[[i]]
      net2=nets[[j]]
      union.net= net1+net2
      #isSymmetric(union.net)
      overlap.edge=sum(union.net==2)-length(diag(union.net))
      similarity=overlap.edge/sum(net1)
      distance=1-similarity;
      dist.mat[i,j]=distance
      dist.mat[j,i]=distance
    }
  }
  rownames(dist.mat)=colnames(dist.mat)=names(nets)
  saveRDS(dist.mat,'src_cell.type.tree_based.on.network.similarity.rds')
}

dist.mat=readRDS('src_cell.type.tree_based.on.network.similarity.rds')
rownames(dist.mat)
my_dist<-as.dist(dist.mat)
#Distance-based phylogenetic trees in R, https://rpubs.com/WalshJake75/674724
my_upgma <- phangorn::upgma(my_dist)

#https://stackoverflow.com/questions/47874486/pheatmap-default-distance-metric-r
#https://stats.stackexchange.com/questions/10347/making-a-heatmap-with-a-precomputed-distance-matrix-and-data-matrix-in-r
my_similarity=1-dist.mat
diag(my_similarity)<-NA
pheatmap::pheatmap(1-my_similarity)
pheatmap::pheatmap(my_similarity)
#https://support.bioconductor.org/p/69971/

pdf('src_cell.type.tree_based.on.network.similarity.pdf',useDingbats = T,width = 12,height = 12)
#plot(my_upgma,cex=0.5)
print(pheatmap::pheatmap(my_similarity),
      cluster_rows = TRUE,
      clustering_method = "complete",
      clustering_distance_rows = "correlation",
      #clustering_distance = "euclidean",
      fontsize_row = 0.1,fontsize_column=0.1)
dev.off()
