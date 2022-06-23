# construct pairwise gene Pearson's correlation matrix for each cell type

# the bigScale2 algorithm was published in: Iacono, Giovanni, Ramon Massoni-Badosa, and Holger Heyn. "Single-cell transcriptomics unveils gene regulatory network plasticity." Genome biology 20, no. 1 (2019): 1-20.
# the source code can be found: https://github.com/iaconogi/bigSCale2
# I modify the original 'Functions.R' script to output and save Pearson's correlation matrix.
source('../modify_bigScale2/ming-Functions.R')
library(Seurat)
library(tidyverse)
library(Matrix);
library(ggplot2);library(gridExtra)
library(RColorBrewer)

genes=readRDS('../01.common.gene.extraction/brain_scRNA-seq_n15c0.005/all_common_genes.rds')
#genes=readRDS('../01.common.gene.extraction/brain_snRNA-seq_n15c0.005/all_common_genes.rds')
common.genes=genes$common.genes
out.dir1='brain_scRNA-seq_n15c0.005/';
#out.dir1='brain_snRNA-seq_n15c0.005/';
if(!dir.exists(out.dir1)) dir.create(out.dir1)

out.dir=paste0(out.dir1,'/bigScale2_coexpr/');
if(!dir.exists(out.dir)) dir.create(out.dir)

## read in processed wholebrain data
if(length(grep('sn',out.dir1))==1){
  file="../../single.cell_datasets/FCA_head/whole_head_filtered_valid.rds";
}else{
  file="../../single.cell_datasets/fly.brain.atlas/wholebrain_filtered_valid.rds";
}

dat=readRDS(file);
df.expr=dat@assays$RNA@data 


## separate female and male, select cell clusters that contain>=500cells in both sexes
min.cell=200;max.cell=Inf; #separate sex
i=apply(table(dat$sex,dat$annotation),2,function(i)sum(i>=min.cell & i<max.cell))
#pick.cell.clusters=names(which(i==2)) 
pick.cell.clusters=names(which(i>=1))
pick.cell.clusters # cell_100: sn,54; sc,60.

dat=subset(dat,annotation %in% pick.cell.clusters)
df.expr=dat@assays$RNA@counts;
dim(df.expr); #sn:12602 50588

sparsity=sum(df.expr==0)/nrow(df.expr)/ncol(df.expr)
cat(out.dir,sparsity,'\n') #sn=0.9537,sc=0.9055


## network construction begin
## effective.expressed.gene >=1umi in max(10,ncell*10%) cells

for(i.sex in c('male','female')){
  for(i.cluster in pick.cell.clusters){
    mat=df.expr[,dat$sex==i.sex & dat$annotation==i.cluster]
    if(ncol(mat)<min.cell){next}
    gene.filter=rownames(mat) %in% common.genes;
    mat=mat[gene.filter,]
    ncell=ncol(mat)
    
    i.cluster2=gsub('\\/','\\.',i.cluster)#there is _ and -, use . to replace /
    out.file=paste0(out.dir,'/ncell_',ncell,'_',i.sex,'_',i.cluster2,'_pearson.rds')
  
    if(!file.exists(out.file)){
      # two things enter next stage
      gene.names=rownames(mat) #vector
      expr.data=mat #dgCMatrix geneXcell
      
      net.out=compute.network(expr.data = expr.data,gene.names = gene.names,
                  clustering='recursive',speed.preset = "fast")
      #Dp=net.out$Dp.mat;
      #Ds=net.out$Ds.mat;
      saveRDS(net.out,out.file);
    }
    
    cat('cell.cluster',i.cluster,'done\n')
  }
}

