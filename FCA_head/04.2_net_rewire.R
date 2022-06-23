
# network randomization for each cell type separately at different cutoffs
library(ggplot2);library(igraph)
library(gridExtra)
library(grid);library(tidyverse)
library(RColorBrewer)

n.simu=20;

path.out=paste0('net_rewire_rep',n.simu,'/');

if(!dir.exists(path.out)){dir.create(path.out)}

(files=Sys.glob('./brain_snRNA-seq_n15c0.005_bigScale2/top_*_networks_cutoff_cell.type.rds'))
(files=files[grep('top_0.05',files)]) #sc,sn


#############################################################
## read in networks

for(file in files){
  (cutoff=strsplit(basename(file),'\\_')[[1]][[2]])
  
  keep.list=readRDS(file);
  names(keep.list);
  plots=list();
  for(i in 1:n.simu){
    outfile=paste0(path.out,'/top_',cutoff,'_shuffle_rep',i,'.rds');
    if(!file.exists(outfile)){
      simu.modular=list();
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
          
        n.edge=length(E(G))
        
        if(n.edge>=999999){ #for computation time consideration
          net0<- rewire(G,keeping_degseq(loops=F,niter=2*n.edge))
        }else{
          net0<- rewire(G,keeping_degseq(loops=F,niter=10*n.edge))
        }
        simu.modular[[cell.type]]=net0;
      }
      saveRDS(simu.modular, outfile)
    } 
  }
}

