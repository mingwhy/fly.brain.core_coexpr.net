
library(Seurat)
input.folder='../brain_scRNA-seq_n15c0.005/'

## read in genes
genes=readRDS(paste0(input.folder,'all_common_genes.rds'))
common.genes=genes$common.genes
length(common.genes) #1738
gene.names=common.genes;

## read in processed wholebrain data
file='../../../single.cell_datasets/fly.brain_Cocaine/Brain.control_filtered.rds'

dat=readRDS(file);
df.expr=dat@assays$RNA@counts 
min.cell=200;

out.dir=paste0('common.genes_',basename(input.folder),'/');
dir.create(out.dir)

for(sex in c('female','male')){
  for(cell.type in unique(dat$annotation)){
    
    mat=df.expr[, dat$annotation==cell.type & dat$sex==sex] #sc
    if(ncol(mat)<min.cell){next}
    
    mat=mat[gene.names,]
    dim(mat)  #2876 2023; 
    mat.count=mat
    
    #set.seed(13579)
    #i=sample(1:ncol(mat),replace = F)
    #mat=mat[,i]
    
    ## if you want to use scran to normalize data
    if(F){
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
    #class(mat.scran);
    class(mat.count);
    ncell=ncol(mat.count);
    cell.type2=gsub('\\/','\\.',cell.type)#there is _ and -, use . to replace /
    out.file=paste0(out.dir,'/ncell_',ncell,'_',sex,'_',cell.type2,'.rds')
    
    #saveRDS(list(mat.count=mat.count,mat.scran=mat.scran),file=out.file)
    saveRDS(list(mat.count=mat.count),file=out.file)
  }
}

