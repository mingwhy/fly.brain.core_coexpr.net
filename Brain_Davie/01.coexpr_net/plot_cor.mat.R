
library(Seurat)
library(ggplot2);library(gridExtra);
options(stringsAsFactors = F)

top=0.05;
genes=readRDS('../01.common.gene.extraction/brain_scRNA-seq_n15c0.005/all_common_genes.rds')
out.dir1='brain_scRNA-seq_n15c0.005/';

#top=0.05;
#genes=readRDS('../01.common.gene.extraction/brain_snRNA-seq_n15c0.005/all_common_genes.rds')
#out.dir1='brain_snRNA-seq_n15c0.005/';

common.genes=genes$common.genes
if(!dir.exists(out.dir1)) dir.create(out.dir1)

out.dir=paste0(out.dir1,'/bigScale2_coexpr/');
if(!dir.exists(out.dir)) dir.create(out.dir)


## plot PCC distribution
(files=Sys.glob(paste0(out.dir,'/*.rds')))

pdf(paste0(out.dir1,'/check_pearson.cor.mat_hist.pdf'),
    useDingbats = T,width = 14)
par(mfrow=c(4,4),oma=c(3,3,3,3),mai=c(1,1,1,0.5),mar=c(2.5,1,1.5,0.5))
for(file in files){
  tmp=readRDS(file)
  cor.mat=tmp$Dp.mat
  cor.mat[is.na(cor.mat)]=0
  x=cor.mat[upper.tri(cor.mat)]
  cut.val=as.numeric(quantile(abs(x),1-top))
  my_hist=hist(x,plot=F,breaks=40)
  my.breaks=sort(c(my_hist$breaks,-1*cut.val,cut.val))
  my_hist=hist(x,breaks=my.breaks,plot=F)
  
  # title
  x=strsplit(basename(file),'_')[[1]]
  ncell=x[2];sex=x[3]
  cell.cluster=gsub('male_|_pearson','',stringr::str_extract(basename(file),'male_.*_pearson'))
  
  # Color vector
  my_color=ifelse(my_hist$breaks<(-1*cut.val),'purple',#rgb(0.2,0.8,0.5,0.5) , 
                   ifelse (my_hist$breaks >cut.val, "purple", rgb(0.2,0.2,0.2,0.2) ))
  my_hist$density=my_hist$counts/sum(my_hist$counts)
  # Final plot
  plot(my_hist, col=my_color , border=F , xlim=c(-1,1),
       #main=paste0(sex,', ',cell.cluster,', ',ncell,'cells'),
       main=gsub('_pearson.rds','',basename(file)),
       ylab='Frequency',
       xlab="")
  abline(v=cut.val,col='darkred',lwd=1)
  abline(v=-1*cut.val,col='darkred',lwd=1)
}
dev.off()

