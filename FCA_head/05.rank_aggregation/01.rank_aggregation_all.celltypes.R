
# compaute rank aggregated score for each edge 
# using all cell types

library(Matrix);library(tidyverse)
library(ggplot2);library(gridExtra)
library(RColorBrewer); #library(RankAggreg) #too slow
library(RobustRankAggreg)

## read in kept cell types
x=unlist(read.table('../02_choose_top.percentage_threshold/sn_keep.cell.types.txt',sep='\n'))
length(x) #39
#df.keep=data.frame(Reduce(`rbind`,sapply(x,function(i) strsplit(i,'\\_'))))
#colnames(df.keep)=c('sex','cell.type')
df.keep=x

topK='top_0.05';
input.dir='../brain_snRNA-seq_n15c0.005/bigScale2_coexpr/';
out.dir='./';

#topK='top_0.1';
#input.dir='../02.coexpr_net/brain_snRNA-seq_n15c0.005/bigScale2_coexpr/';
#out.dir=paste0('brain_snRNA-seq_n15c0.005_bigScale2_',topK);

dir.create(out.dir);
out.file=paste0(out.dir,'/cor.rank.list.rds');

# collect bigscale output nets  
(files=Sys.glob(paste0(input.dir,'/*.rds')))
tcs=sapply(files,function(file){
  gsub('ncell_\\d+_|_pearson.rds','',basename(file))
})
sum(tcs %in% df.keep) #39
files=files[tcs %in% df.keep]
length(files) #39

if(!file.exists(out.file)){
  cor.rank.list=list()
  for(i in 1:length(files)){
    filename=files[i]
    cluster.name=paste0(strsplit(basename(filename),'\\_')[[1]][c(3,4)],collapse='_')
    cat('cell.cluster',cluster.name,'\n')
    
    x=readRDS(filename)
    Dp=x$Dp.mat
    ## check for NA
    x=Dp[upper.tri(Dp)]
    cat('cell.cluster',cluster.name,'NA.percentage =',sum(is.na(x))/length(x),'\n');
    ## if there is NA, replace them with 0
    if(sum(is.na(x))!=0){
      Dp[is.na(Dp)]=0
      x=Dp[upper.tri(Dp)]
      cat('cell.cluster',cluster.name,'change NA into 0, na.prop=',sum(is.na(x))/length(x),'\n');
    }
    
    # change Dp matrix into data.frame  
    ind <- which(upper.tri(Dp, diag = F), arr.ind = TRUE)
    nn <- dimnames(Dp)
    tmp=data.frame(row = nn[[1]][ind[, 1]],
                   col = nn[[2]][ind[, 2]],
                   val = Dp[ind])
    edges <- do.call(paste, c(tmp[ , 1:2 ], sep=";")) #much faster
    #edges <- apply( tmp[ , 1:2 ] , 1 , paste , collapse = ";" )
    abs.cor=abs(tmp[,3])
    edges=edges[order(abs.cor,decreasing = T)]
    
    cor.rank.list[[cluster.name]]<-edges;
    cat('cell.cluster',cluster.name,'is done\n')
  }
  sapply(cor.rank.list,length) 
  saveRDS(cor.rank.list,file=out.file);
  
  }else{
    cor.rank.list=readRDS('cor.rank.list.rds');
}
sapply(cor.rank.list,length) 

## use RankAggregate to aggregate ranks
if(F){
  cor.rank.df=Reduce(`rbind`,cor.rank.list)
  dim(cor.rank.df)
  rownames(cor.rank.df)=names(cor.rank.list)
  
  rm('cor.rank.list') #too big file
  
  cat(ncol(cor.rank.df),'edges in total\n')
  #pick=c(100000,1000000,ncol(cor.rank.df))
  pick=ncol(cor.rank.df);
  for(i in pick){
    dat=cor.rank.df[,1:i]
    dat.list=lapply(1:nrow(dat),function(i){dat[i,]})
    #edge.rank=RankAggreg(dat,k=ncol(dat),seed=100,method='CE',distance='Kendall',verbose=F);
    edge.rank.RRA=aggregateRanks(glist=dat.list,N=ncol(dat),method='RRA',full=T)
    #edge.rank.stuart=aggregateRanks(glist=dat.list,N=ncol(dat),method='stuart')
    
    cat('top',i,'is done\n')
    saveRDS(edge.rank.RRA,paste0(out.dir,'/edge.rank_all',i,'.rds'))
  }
}
