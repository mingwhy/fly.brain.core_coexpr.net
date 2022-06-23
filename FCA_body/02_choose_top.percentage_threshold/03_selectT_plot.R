
library(ggplot2);library(gridExtra)

#input.file='sc_sparsity.rds';
input.file='sn_sparsity.rds';
##################################################
## read in ncell per cell type, run server
all.ncell=list();
#files=Sys.glob('common.genes_brain_scRNA-seq_n15c0.005/*.rds');
files=Sys.glob('common.genes_brain_snRNA-seq_n15c0.005/*.rds');
length(files) #35

for(file in files){
  x=basename(file)
  ncell=as.numeric(strsplit(x,'\\_')[[1]][[2]])
  cell.type=gsub('ncell_\\d*_|.rds','',x)
  if(ncell>=200){
    all.ncell[[cell.type]]=ncell
  }
}
length(all.ncell) #35

## begin plot
df.sparse=readRDS(input.file)
max(df.sparse$sparsity)
#df.sparse=df.sparse[df.sparse$sparsity<0.5,]
dim(df.sparse) #35

library(tidyverse)
library(COGENT)
#out.file='sc_bigscale_rep10.rds'
out.file='sn_bigscale_rep10.rds'

if(!file.exists(out.file)){
  #inp.dir='sc_bigscale_rep10/';
  inp.dir='sn_bigscale_rep10/';
  (files=Sys.glob(paste0(inp.dir,'/*.rds')))
  files=files[-grep('splitData',files)]
  
  thresholds<- c(0.7,0.9,0.95,0.99,0.995)
  result=c();
  for(file in files){
    rep=stringr::str_extract(file,'rep\\d+_ncell')
    rep=gsub('_ncell','',rep)
    
    cell.type=gsub('.*_rep\\d+_ncell_\\d+_|.rds','',basename(file));
    cat('start',cell.type, rep,'\n')
    
    out=readRDS(file)
    mat1=out$net.out1$Dp.mat;
    mat2=out$net.out2$Dp.mat;
    mat.out=list(mat1,mat2)
    abs.mat.out=lapply(mat.out,abs)
    sapply(abs.mat.out,dim)
    
    for(th in thresholds){
      As=lapply(abs.mat.out,function(A){
        threshold <- quantile(A[upper.tri(A)], th, na.rm=TRUE)
        A <- 1*(A>=threshold); diag(A) <- 0; 
        A[is.na(A)]=0
        A
      })
      result=rbind(result,c(cell.type,rep,th,
                            getEdgeSimilarityCorrected(As,type="expected")$correctedSimilarity))
    }
  }
  saveRDS(result,out.file)
}else{
  df=as.data.frame(readRDS(out.file))
  colnames(df)=c('cell.type','rep','th','value') #correctedSimilarity
  length(unique(df$cell.type))
  df=df[df$cell.type %in% names(all.ncell),]
  #df$th=as.numeric(df$th);
  dim(df)
  df$value=as.numeric(df$value);
  df$threshold=factor(1- as.numeric(df$th))
  df.all.sum=df %>% group_by(cell.type,th) %>% 
    summarise(median.ratio=median(value))
  df.all.sum$threshold=factor(1-as.numeric(df.all.sum$th))
  

  plots=list();i.plot=0;
  for(cell.type in unique(df.all.sum$cell.type)){
    x=df[df$cell.type==cell.type,]
    x$threshold=as.numeric(as.character(x$threshold))
    i.plot=i.plot+1;
    plots[[i.plot]]<- ggplot(x,aes(x=threshold,y=value))+
      geom_point()+
      geom_smooth()+
      #ggtitle("Threshold choice for gene co-expression networks") +
      #scale_y_continuous("Density-adjusted consistency") +
      ylab("Signal to noise score") +
      theme_classic()+ylim(0,max(x$value))+
      scale_x_log10("Threshold cutoff", breaks=x$threshold)+
      ggtitle(cell.type)
    
  }
  
  p.median<- ggplot(df.all.sum,aes(x=threshold,y=median.ratio,group=cell.type))+
           geom_point(size=0.5)+geom_line()+theme_classic(base_size = 20)+
           ggtitle("Threshold choice for gene co-expression networks") +
           #scale_y_continuous("Density-adjusted consistency") +
           scale_y_continuous("Signal to noise score")+ #scale_y_log10()+
           xlab("Threshold cutoff")
  
 
}

df.all.sum[df.all.sum$threshold==0.05,]->tmp
tmp=tmp[order(tmp$median.ratio),]
dim(tmp)#35 cell cluster
sum(tmp$median.ratio<0.02) #4
rm.cell.type=tmp[tmp$median.ratio<0.02,]

## p.median
print(p.median);
df.all.sum$type=1;
#rm.cell.types=unique(df.all.sum[df.all.sum$median.ratio<0.02,]$cell.type)
#df.all.sum[df.all.sum$cell.type %in% rm.cell.types,]$type=0
df.all.sum[df.all.sum$cell.type %in% rm.cell.type$cell.type,]$type=0;

pdf('sn_median.score_over.thresholds.pdf',useDingbats = T,width = 12)
print( ggplot(df.all.sum,aes(x=threshold,y=median.ratio,group=cell.type))+
         geom_point(size=0.5)+geom_line(aes(col=factor(type)))+theme_classic(base_size = 20)+
         scale_color_manual(values=c("#999999", "#56B4E9"))+
         ggtitle("Threshold choice for gene co-expression networks") +
         #scale_y_continuous("Density-adjusted consistency") +
         scale_y_continuous("Signal to noise score")+ #scale_y_log10()+
         xlab("Threshold cutoff")+
         theme(legend.position = 'none') )
dev.off()

keep.cell.types=unique(df.all.sum[df.all.sum$type==1,]$cell.type) #31 remain
write.table(keep.cell.types,'sn_keep.cell.types.txt',quote=F,row.names = F,col.names = F)


####
layout_matrix = matrix(1:16,4,4,byrow = T)
length(plots)
#pdf('sc_bigscale_rep10.pdf',useDingbats = T,width = 12,height = 9)
pdf('sn_bigscale_rep10.pdf',useDingbats = T)
#print(p.median);
grid.arrange(arrangeGrob(grobs= plots[1:16],layout_matrix = layout_matrix))
grid.arrange(arrangeGrob(grobs= plots[17:32],layout_matrix = layout_matrix))
grid.arrange(arrangeGrob(grobs= plots[33:35],layout_matrix = layout_matrix))
#grid.arrange(arrangeGrob(grobs= plots[49:54],layout_matrix = layout_matrix))
#grid.arrange(arrangeGrob(grobs= plots[65:67],layout_matrix = layout_matrix))
dev.off()

min(subset(df.all.sum,th==0.95)$median.ratio)
#  0.01181517

##################################################
## plot score VS sparsity, run local
df.out.sparse=readRDS(input.file)
head(df.out.sparse);
dim(df.out.sparse) #35

head(df.all.sum)
unique(df.all.sum$cell.type) #35 sn
sum(df.all.sum$cell.type %in% df.out.sparse$cluster) #all exist

dfc=merge(df.out.sparse,df.all.sum,by.x='cluster',by.y='cell.type')
head(dfc);
summary(dfc$sparsity);

T=0.05; #sc,which threshold to use in ploting S2N score and ncell, sparsity
#T=0.1; #sn
#dfc=subset(dfc,threshold==T & sparsity<0.5)
dfc=subset(dfc,threshold==T)
dim(dfc) #35tc

library(viridis)
# ncell ~ score
p= ggplot(dfc,aes(x=ncell,y=median.ratio))+
  geom_point()+theme_classic(base_size = 20)+
  scale_x_log10()+
  #scale_color_viridis() +
  ylab('Signal-to-Noise Score')


lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  f <- summary(m)$fstatistic
  pval <- as.numeric(pf(f[1],f[2],f[3],lower.tail=F))
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2~","~~italic(p)~"="~pval, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        pval = format(pval,digits=3)))
  as.character(as.expression(eq));
}
df=data.frame(x=dfc$ncell,y=dfc$median.ratio)
p1=p + geom_smooth(method='lm')+
  geom_text(aes(x=520,y=0.5),size=6,label = lm_eqn(df), parse = TRUE,stat = "unique")
  #geom_text(aes(x=400,y=0.3),label = lm_eqn(df), parse = TRUE) #for sn

# sparsity ~ score
p.sparse= ggplot(dfc,aes(x=sparsity,y=median.ratio))+
  geom_point()+theme_classic(base_size = 20)+
  scale_x_log10()+
  #scale_color_viridis() +
  ylab('Signal-to-Noise Score')
df=data.frame(x=dfc$sparsity,y=dfc$median.ratio)
p.sparse1=p.sparse + geom_smooth(method='lm')+
  geom_text(aes(x=0.3,y=0.5),size=6,label = lm_eqn(df), parse = TRUE,stat = "unique")

df2=subset(dfc,sparsity>0.3)
df=data.frame(x=df2$sparsity,y=df2$median.ratio)
p.sparse2= ggplot(df2,aes(x=sparsity,y=median.ratio))+
  geom_point()+theme_classic(base_size = 20)+
  scale_x_log10()+
  #scale_color_viridis() +
  ylab('Signal-to-Noise Score')+
  geom_smooth(method='lm')+
  geom_text(aes(x=0.38,y=0.3),size=6,label = lm_eqn(df), parse = TRUE,stat = "unique")


pdf(paste0(input.file,'_T',T,'_S2N.ncell.pdf'),useDingbats = T,
    width = 8,height = 12)
#grid.arrange(p,p1,ncol=2)
#grid.arrange(p.sparse,p.sparse1,p.sparse2,ncol=3)
grid.arrange(p1,p.sparse1,ncol=1)
dev.off()

if(F){
  print( ggplot(dfc,aes(x=sparsity,y=median.ratio,col=ncell))+
    facet_wrap(.~threshold,ncol=5)+geom_point()+theme_classic()+
      ylab('Signal-to-Noise Score')+
    scale_color_viridis()  )
  print( ggplot(subset(dfc,ncell>=500),aes(x=sparsity,y=median.ratio,col=ncell))+
           facet_wrap(.~th,ncol=5)+geom_point()+theme_classic()+
           ylab('Signal-to-Noise Score')+
           scale_color_viridis()  )
  
  tmp=df.all.sum[df.all.sum$cell.type %in% 
                   unique(df.all.sum$cell.type)[c(61,60,29)] ,]
  print( ggplot(tmp,aes(x=threshold,y=median.ratio,group=cell.type))+
           geom_point(size=0.5)+geom_line()+theme_classic(base_size = 20)+
           #ggtitle("Threshold choice for gene co-expression networks") +
           #scale_y_continuous("Density-adjusted consistency") +
           scale_y_continuous("Signal to noise score")+ #scale_y_log10()+
           xlab("Threshold cutoff")
  ) 

}


