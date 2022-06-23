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

genes=readRDS('./brain_scRNA-seq_n15c0.005/all_common_genes.rds')
#genes=readRDS('../01.common.gene.extraction/brain_snRNA-seq_n15c0.005/all_common_genes.rds')
common.genes=genes$common.genes #1738 genes
out.dir1='brain_scRNA-seq_n15c0.005/';
#out.dir1='brain_snRNA-seq_n15c0.005/';
if(!dir.exists(out.dir1)) dir.create(out.dir1)

out.dir=paste0(out.dir1,'/bigScale2_coexpr/');
if(!dir.exists(out.dir)) dir.create(out.dir)

## read in processed wholebrain data
file='../../single.cell_datasets/fly.brain_Cocaine/Brain.control_filtered.rds'


dat=readRDS(file);
df.expr=dat@assays$RNA@data 

dim(df.expr); #10386 41520
min.cell=200;
## network construction begin
## effective.expressed.gene >=1umi in max(10,ncell*10%) cells

for(i.sex in c('male','female')){
  for(i.cluster in unique(dat$annotation)){
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

