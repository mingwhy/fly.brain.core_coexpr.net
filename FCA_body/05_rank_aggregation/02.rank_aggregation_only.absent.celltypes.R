# on server
# calculate RankAggreg for edges in absent cell clusters
options(stringsAsFactors = F)
library(igraph)
library(Matrix);
library(RobustRankAggreg)
library(parallel)

sample.nedge=100;
edge.groups=seq(2,100,2) #sn_FCA_body

n.core=4;
#n.core=length(edge.groups);

# read in networks at different cutoffs to edge commonality distribution
topK='top_0.05';
#(files=Sys.glob('../brain_snRNA-seq_n15c0.005_bigScale2/top_*_networks_cutoff_cell.type.rds'))
file='../brain_snRNA-seq_n15c0.005_bigScale2/top_0.05_networks_cutoff_cell.type.rds'
#file='top_0.05_networks_cutoff_cell.type.rds';

#out.dir0=paste0('brain_scRNA-seq_n15c0.005_bigScale2_',topK,'/');
out.dir0=paste0('RAA_output');
dir.create(out.dir0);

output.dir=paste0(out.dir0,'/result_RAA_sub.mat');
dir.create(output.dir)

# read in edge abs PCC value for each cell type, this may take some time
cor.rank.list=readRDS('cor.rank.list.rds'); #output of `01.rank_aggregation_all.33.celltypes.R`

cor.rank.df=Reduce(`rbind`,cor.rank.list)
dim(cor.rank.df) #31 377146 # #cluster X #edge
rownames(cor.rank.df)=names(cor.rank.list)

rm('cor.rank.list')


(cutoff=strsplit(basename(file),'\\_')[[1]][[2]])

# get edge 0/1 binary pattern across cell types
#(outfile=paste0(out.dir0,'/top_',cutoff,'_edge_in.cell.types.rds'));  

(outfile=paste0('top_',cutoff,'_edge_in.cell.types.rds'));  
if(!file.exists(outfile)){
  
  keep.list=readRDS(file);
  #isSymmetric(keep.list[[1]]) #TRUE
  
  edge.list=lapply(keep.list,function(Dp){
    ind <- which(upper.tri(Dp, diag = F), arr.ind = TRUE)
    nn <- dimnames(Dp)
    tmp=data.frame(row = nn[[1]][ind[, 1]],
                   col = nn[[2]][ind[, 2]],
                   val = Dp[ind])
    tmp
  })
  
  sapply(edge.list,dim)
  saveRDS(edge.list,outfile)
  
}else{
  edge.list=readRDS(outfile);
}


cat('begin RAA\n');
# for each edge, select cell types it is absent, get their abs.PCC values, aggregate
names(edge.list) #cell.type names
#sum(edge.list[[1]][,3]) # number.of.edges with the topK cutoff in this cell cluster
tmp=edge.list[[1]]
edges <- do.call(paste, c(tmp[ , 1:2 ], sep=";")) #much faster
#edges <- apply( tmp[ , 1:2 ] , 1 , paste , collapse = ";" )

values=lapply(edge.list, '[',,3) #extract the 1st column in each list matrix member
val.mat=Reduce(`cbind`,values)
dim(val.mat) #377146     31
rownames(val.mat)=edges #edge names
colnames(val.mat)=names(edge.list) #cell type names

commonality=apply(val.mat,1,sum)
cat('max.commonality.score',max(commonality),'\n') #24
length(commonality)
unique(commonality)

## begin rank aggregation for edges in absent cell types 
for(i in edge.groups){
  
  if(sum(commonality==i)==0){next} #some commonality group doesn't have any edges

  edge.group=val.mat[commonality==i,,drop=F] #edge by cell type matrix
  sample.nedge.for.this.group=sample.nedge
  if(nrow(edge.group)<sample.nedge){sample.nedge.for.this.group=nrow(edge.group)}
  
  set.seed(220107)
  edge.group=edge.group[sample(1:nrow(edge.group),sample.nedge.for.this.group,replace=F),,drop=F]

  (filename=paste0(out.dir0,'/RAA.per.edge_at.edge.',i,'_sample',sample.nedge.for.this.group,'edges.rds'))  
  dir=paste0(output.dir,'/RAA_sub.mat_',i);
  if(!dir.exists(dir)){dir.create(dir)}

  df=as.numeric()
  r<-mclapply(1:nrow(edge.group), function(edge.i){
    one.edge=rownames(edge.group)[edge.i]
    pick.cells=names(which(edge.group[one.edge,]==0)) #cell types that lack this edge
    track.name=paste(which(edge.group[one.edge,]==0),collapse = '_')
    #(out.file=paste0(dir,'/RRA_',track.name,'.rds'))
    (out.file=paste0(dir,'/RRA_',edge.i,'.rds'))
    cat('edge.i ',edge.i,', ',one.edge,' ',track.name,'\n');
    
    cat('save file',out.file,'\n')
    
    dat=cor.rank.df[pick.cells,]
    dat.list=lapply(1:nrow(dat),function(i){dat[i,]})
    #edge.rank=RankAggreg(dat,k=ncol(dat),seed=100,method='CE',distance='Kendall',verbose=F);
    edge.rank.RRA=aggregateRanks(glist=dat.list,N=ncol(dat),method='RRA',full=T)
    #edge.rank.RRA=aggregateRanks(glist=dat.list,N=ncol(dat),method='stuart',full=T)
    saveRDS(edge.rank.RRA,out.file)
    edge.rank.RRA[one.edge,]
  },mc.cores=n.core)

  df= Reduce(`rbind`,r); 
  cat('edge at',i,'is done\n');
  saveRDS(df,filename)
  
}



