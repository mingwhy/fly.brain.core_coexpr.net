
library(Seurat)

## read in raw data
file="../../../single.cell_datasets/fly.brain.atlas/wholebrain_filtered_valid.rds";
dat=readRDS(file);
df.expr=dat@assays$RNA@data 
unique(dat$annotation)

## read in common expr genes
genes=readRDS('../01.common.gene.extraction/brain_scRNA-seq_n15c0.005/all_common_genes.rds')
common.genes=genes$common.genes
length(common.genes) #2088
set.seed(20220120)
pick=matrix(sample(common.genes,8,replace = F),nrow=4)

## read in cell type co.expr mat and make plots
out.dir='brain_scRNA-seq_n15c0.005/bigScale2_coexpr';
(files=Sys.glob(paste0(out.dir,'/*.rds')))
files=files[grep('female_G-KC|female_Mip|_male_T2',files)]
length(files) # 67 


pdf(paste0('example_4gene.pairs.pdf'),width=14,height = 10, useDingbats = TRUE)
par(mfrow=c(4,6),mar=c(5,5,3,1)) #33 cell types, 7x5

for(j in 1:nrow(pick)){
  file=files[1]
 
  # title
  x=strsplit(basename(file),'_')[[1]]
  ncell=x[2];sex=x[3]
  cell.cluster=gsub('male_|_pearson','',stringr::str_extract(basename(file),'male_.*_pearson'))
  
  net=readRDS(file)
  tot.score=net$tot.scores
  Dp=net$Dp.mat
  cat('cell.cluster',cell.cluster,'ncluster',max(net$clusters),'\n');
  
  gene1=pick[j,1];gene2=pick[j,2]
  pos=match(pick[j,],rownames(Dp))
  
  for(i in 1:length(files)){
    filename=files[i]
    
    x=strsplit(basename(filename),'_')[[1]]
    ncell=x[2];sex=x[3]
    cell.cluster=gsub('male_|_pearson','',stringr::str_extract(basename(filename),'male_.*_pearson'))
    cat('cell.cluster',cell.cluster,'\n')
    
    ## raw count
    mat=df.expr[,dat$sex==sex & dat$annotation==cell.cluster]
    pair.raw=mat[c(gene1,gene2),]
    cor.value=round(cor(pair.raw[1,],pair.raw[2,]),3)
    plot(jitter(pair.raw[1,]),jitter(pair.raw[2,]),
         pch=16,cex.lab=1.5,col=rgb(0,0,0,alpha = 0.4),
         xlab=paste(gene1,'(raw count)'),
         ylab=paste(gene2,'(raw count)'),main=paste0(sex,', ',cell.cluster,'\nPCC=',cor.value))
    
    ## read in Z score from bigScale2
    net=readRDS(filename)
    cat('cell.cluster',cell.cluster,'ncluster',max(net$clusters),'\n');
    
    tot.score=net$tot.scores
    Dp=net$Dp.mat
    ## check for NA
    x=Dp[upper.tri(Dp)]
    cat('cell.cluster',cell.cluster,'NA.percentage =',sum(is.na(x))/length(x),'\n');
    ## if there is NA, replace them with 0
    if(sum(is.na(x))!=0){
      Dp[is.na(Dp)]=0
      x=Dp[upper.tri(Dp)]
      cat('cell.cluster',cell.cluster,'change NA into 0, na.prop=',sum(is.na(x))/length(x),'\n');
    }
    #cor(tot.score[,pos])
    (cor.value=round(Dp[gene1,gene2],3))
    plot(jitter(tot.score[,pos[1]]),jitter(tot.score[,pos[2]]),pch=16,cex.lab=1.5,col='steelblue',
        xlab=paste(gene1,'(Z-score)'),
        ylab=paste(gene2,'(Z-score)'),main=paste0(sex,', ',cell.cluster,'\nPCC=',cor.value))
  }
  #plot.new();plot.new()
}
dev.off()



