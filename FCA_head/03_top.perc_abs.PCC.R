
# plot pearson.corr.coeff distribution for each cell type
# use different top percentile cutoffs on absolute values of pearson.cor.matrix 
# to get cell type-specific gene co-expression networks

library(Matrix);library(tidyverse)
library(ggplot2);library(gridExtra)
library(RColorBrewer)

## read in kept cell types
x=unlist(read.table('02_choose_top.percentage_threshold/sn_keep.cell.types.txt',sep='\n'))
length(x) #39
#df.keep=data.frame(Reduce(`rbind`,sapply(x,function(i) strsplit(i,'\\_'))))
#colnames(df.keep)=c('sex','cell.type')
df.keep=x

topn.list=c(0.005,0.01,0.05,0.1,0.3,0.5,0.7) 
path.in='./brain_snRNA-seq_n15c0.005/bigScale2_coexpr/';
path='./brain_snRNA-seq_n15c0.005_bigScale2/';
dir.create(path)

# collect bigscale output nets 
(files=Sys.glob(paste0(path.in,'/*.rds')))
length(files) #76 for sn

cor.list=list()
for(i in 1:length(files)){
  filename=files[i]
  x=strsplit(basename(filename),'\\_|.rds')[[1]]
  sex=x[3]
  cell.cluster=x[4]
  cat('cell.cluster',sex,cell.cluster,'\n')
  name=paste0(sex,'_',cell.cluster)
  
  if(sum(name %in% df.keep)==0){next}
  
  x=readRDS(files[i])
  Dp=x$Dp.mat; #there is NA
  cor.list[[name]]<-Dp;
  cat('cell.type',name,'is done\n')
}

sapply(cor.list,dim) 


# get the the corresponding abs(pearson.cor.coeff) value for each top percentile cutoff
# save result in a data.frame: cell.type, topN, cutoff.value, ngene, nedge
# one data.frame for one cell.type
names(cor.list)

for(cutoff in topn.list){
  (out.file=paste0(path,'/top_',cutoff,'_networks_cutoff_cell.type.rds'))
  topn.out=c();
  if(!file.exists(out.file)){
    keep.list=list();
    for(cell.type in names(cor.list)){
      ## remember it's absolute pearson.cor.coeff values
      #abs.cor=abs(cor2.list[[cell.type]])
      abs.cor=abs(cor.list[[cell.type]])  
      abs.cor[is.na(abs.cor)]=0;
      all.cor=abs.cor[upper.tri(abs.cor)]
      
      all.ngene=nrow(abs.cor);
      all.nedge=length(all.cor);
      (cutoff.value=quantile(all.cor,1-cutoff)) 
      #(cutoff.value=sort(all.cor,decreasing = T)[cutoff])
      cat('cell.type',cell.type,'top',cutoff,'value',cutoff.value,'\n')
      
      tmp=abs.cor
      tmp[abs.cor>=cutoff.value]=1;
      tmp[abs.cor<cutoff.value]=0;
      keep.list[[cell.type]]=tmp;
      
      diag(tmp)=0
      gene.degree=apply(tmp,1,sum)
      ngene=length(gene.degree)-sum(gene.degree==0)
      nedge=sum(tmp)/2 #symmetric matrix
      #cell.type, topN, cutoff.value, ngene, nedge
      
      x=c(cell.type,all.ngene,all.nedge,cutoff,cutoff.value,ngene,nedge);
      topn.out=rbind(topn.out,x);
    }
    
    colnames(topn.out)=c('cell.type','all.ngene','all.nedge','topN.cutoff','cutoff.value','ngene','nedge')
    write.table(topn.out,paste0(path,'/top_',cutoff,"_cor.value.txt"),quote=F,row.names = F)
    saveRDS(keep.list,out.file);
    #sink();
  }
}  


