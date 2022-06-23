
library(ggplot2)
files=Sys.glob('common.genes_brain_snRNA-seq_n15c0.005/*.rds')

out.file='sn_sparsity.rds';

if(!file.exists(out.file)){
  out.sparse=lapply(files,function(file){
    mat.list=readRDS(file);    
    mat.count=as.matrix(mat.list$mat.count);
    obs.sparse=sum(mat.count!=0)/nrow(mat.count)/ncol(mat.count)#0.6528626
    #cat(basename(file),'obs.sparse',obs.sparse,'\n')
    x=strsplit(basename(file),'\\_')[[1]]
    ncell=x[2];sex=x[3];cluster=gsub('ncell_\\d*\\_|.rds','',basename(file));
    c(basename(file),ncell,sex,cluster,obs.sparse)
  })
  
  df.out.sparse=as.data.frame(Reduce(`rbind`,out.sparse))
  colnames(df.out.sparse)=c('cell.cluster','ncell','sex','cluster','sparsity')
  df.out.sparse$ncell=as.numeric(df.out.sparse$ncell)
  df.out.sparse$sparsity=as.numeric(df.out.sparse$sparsity)
  #hist(df.out.sparse$sparsity)
  df.out.sparse=df.out.sparse[df.out.sparse$ncell>=200,]
  summary(df.out.sparse$sparsity)
  df.out.sparse=df.out.sparse[order(df.out.sparse$sparsity,decreasing = T),]
  
  saveRDS(df.out.sparse,out.file)
}

df.out.sparse=readRDS(out.file)
df.out.sparse$cluster=factor(df.out.sparse$cluster,df.out.sparse$cluster)
max(df.out.sparse$sparsity) #0.5263604

pdf(paste0(out.file,'.pdf'),height = 15,width = 12)
print( ggplot(df.out.sparse,aes(x=cluster,y=sparsity))+
         geom_bar(stat='identity',fill=NA,col='darkblue')+
         scale_y_continuous(breaks=c(0,0.2,0.4,0.5,0.55),labels=c(0,0.2,0.4,0.5,0.55))+
         theme_classic()+
         theme(
           axis.title=element_text(size=16),
           axis.text.y=element_text(size=15),
           axis.text.x=element_text(size=17,angle=0,vjust=1,hjust=1))+
         coord_flip()+
         xlab('')
)
print( ggplot(df.out.sparse,aes(x=sparsity,y=ncell))+
        geom_point()+theme_classic()+
         scale_x_log10()+scale_y_log10()+
         ylab('ncell')
)
dev.off()

summary(df.out.sparse$sparsity)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.2350  0.2661  0.3058  0.3334  0.3790  0.5360 

