
library(Seurat)
input.folder='../01.common.gene.extraction/brain_scRNA-seq_n15c0.005/';
#input.folder='../01.common.gene.extraction/brain_snRNA-seq_n15c0.005/';

## read in genes
genes=readRDS(paste0(input.folder,'all_common_genes.rds'))
c.genes=unlist(genes$genes)
all.genes=genes$all.genes
common.genes=genes$common.genes
n.cell.cluster=length(genes$genes)

gene.names=names(which(table(c.genes)==n.cell.cluster))
gene_count=length(gene.names) #2088 for sc, 842 for sn
cat(input.folder,'#gene',gene_count,'express in 100% #cell.cluster',n.cell.cluster,'\n');

saveRDS(gene.names,paste0('common.genes_',basename(input.folder),'.rds'))

## select a cell cluster
## read in processed wholebrain data
if(length(grep('sc',input.folder))>0){
  file="../../single.cell_datasets/fly.brain.atlas/wholebrain_filtered_valid.rds";
}else{
  file="../../single.cell_datasets/FCA_head/whole_head_filtered_valid.rds";
}
dat=readRDS(file);
df.expr=dat@assays$RNA@counts 

## separate female and male, select cell clusters that contain>=200cells in both sexes
min.cell=200;max.cell=Inf;
i=apply(table(dat$sex,dat$annotation),2,function(i)sum(i>=min.cell & i<max.cell))
#pick.cell.clusters=names(which(i==2)) #79
pick.cell.clusters=names(which(i>=1)) #79
pick.cell.clusters #min200: sn,40. sc,38. min100: sn,54; sc,60.


out.dir=paste0('common.genes_',basename(input.folder),'/');
dir.create(out.dir)

for(sex in c('female','male')){
  for(cell.type in pick.cell.clusters){
    
    mat=df.expr[, dat$annotation==cell.type & dat$sex==sex] #sc

    mat=mat[gene.names,]
    dim(mat)  #2876 2023; 
    mat.count=mat
    
    #set.seed(13579)
    #i=sample(1:ncol(mat),replace = F)
    #mat=mat[,i]
    
    ## if you want to use scran to normalize data
    if(T){
      #https://bioconductor.org/packages/devel/bioc/vignettes/scran/inst/doc/scran.html
      library(scran)
      sce <- SingleCellExperiment(list(counts=mat),
                                  colData=DataFrame(cell.type=rep(cell.type,ncol(mat)),sex=sex),
                                  rowData=DataFrame(gene=rownames(mat)) )
      sce
      clusters <- quickCluster(sce)
      sce <- computeSumFactors(sce, clusters=clusters)
      summary(sizeFactors(sce))
      sce <- logNormCounts(sce)
      log.sce=sce@assays@data$logcounts
      dim(log.sce)
      sum(mat==0) #the same number of 0
      sum(log.sce==0) #the same number of 0
      mat.scran=log.sce
    }
    class(mat.scran);
    class(mat.count);
    ncell=ncol(mat.scran);
    cell.type2=gsub('\\/','\\.',cell.type)#there is _ and -, use . to replace /
    out.file=paste0(out.dir,'/ncell_',ncell,'_',sex,'_',cell.type2,'.rds')
    
    saveRDS(list(mat.count=mat.count,mat.scran=mat.scran),
            file=out.file)
  }
}

