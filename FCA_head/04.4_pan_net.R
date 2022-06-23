
# combine cell type-specific gene co-expression networks to 
# generate a pan-network at a given topN cutoff

options(stringsAsFactors = F)
library(igraph)
library(Matrix);library(tidyverse)
library(ggplot2);library(gridExtra)
library(RColorBrewer);library(ggpubr)
path='./brain_snRNA-seq_n15c0.005_bigScale2/';

# read in networks at different cutoffs to edge commonality distribution
(files=Sys.glob(paste0(path,'/top_*_networks_cutoff_cell.type.rds')))

for(file in files){
  (cutoff=strsplit(basename(file),'\\_')[[1]][[2]])
  
  (outfile=paste0(path,'/top_',cutoff,'_pan.rds'));
  if(!file.exists(outfile)){
    keep.list=readRDS(file);
    #isSymmetric(keep.list[[1]]) #TRUE
    
    # as different cell types have differnet dimension of cor.mat
    sapply(keep.list,dim)
    all.genes=unique(unlist(lapply(keep.list,function(x) rownames(x))))
    length(all.genes)
    
    pan<-Reduce(`+`,keep.list)
    dim(pan);     
    saveRDS(pan,outfile)

  }else{
    pan=readRDS(outfile);
    #isSymmetric(pan); #TRUE
    # does every gene in pan have at least one edge
    g=graph_from_adjacency_matrix(pan,mode='undirected')
    cat('top',cutoff,sum(degree(g)==0),'genes in the pan have 0 degree\n')
  }
}

