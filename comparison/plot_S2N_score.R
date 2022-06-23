
library(tidyverse)

path2=c('~/Documents/fly.brain_coexpr.net_2022.06/Brain_Davie/02_choose_top.percentage_threshold/',
        '~/Documents/fly.brain_coexpr.net_2022.06/Brain_Baker/02_choose_top.percentage_threshold/',
        '~/Documents/fly.brain_coexpr.net_2022.06/FCA_head/02_choose_top.percentage_threshold/',
        '~/Documents/fly.brain_coexpr.net_2022.06/FCA_body/02_choose_top.percentage_threshold/');

datasets=c('brain, Davie et al. 2018','brain, Baker et al. 2021',
           'head, Li et al. 2022', 'body, Li et al. 2022')

plots<-lapply(path2,function(path.i){
  
  all.ncell=list();
  files=Sys.glob(paste0(path.i,'/common.genes_brain_s*RNA-seq_n15c0.005/*.rds'));
  length(files)
  
  for(file in files){
    x=basename(file)
    ncell=as.numeric(strsplit(x,'\\_')[[1]][[2]])
    cell.type=gsub('ncell_\\d*_|.rds','',x)
    if(ncell>=200){
      all.ncell[[cell.type]]=ncell
    }
  }
  length(all.ncell) 
  
  out.file=Sys.glob(paste0(path.i,'s*_bigscale_rep10.rds'))
  print(out.file)

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
  

  df.all.sum[df.all.sum$threshold==0.05,]->tmp
  tmp=tmp[order(tmp$median.ratio),]
  dim(tmp)#35 cell cluster
  sum(tmp$median.ratio<0.02) #0
  rm.cell.type=tmp[tmp$median.ratio<0.02,]

  df.all.sum$type=1;
  df.all.sum[df.all.sum$cell.type %in% rm.cell.type$cell.type,]$type=0;
  df.all.sum$type=factor(df.all.sum$type,levels=c('1','0'))
  print( ggplot(df.all.sum,aes(x=threshold,y=median.ratio,group=cell.type))+
           geom_point(size=0.5)+geom_line(aes(col=factor(type)))+theme_classic(base_size = 12)+
           scale_color_manual(values=c("#56B4E9","#999999"))+
           #scale_color_manual(values=c("#56B4E9"))+
           #ggtitle("Threshold choice for gene co-expression networks") +
           #scale_y_continuous("Density-adjusted consistency") +
           scale_y_continuous("Signal to noise score")+ #scale_y_log10()+
           xlab("Threshold cutoff")+
           theme(legend.position = 'none') )
})
plots[[1]]


pdf('plot_S2N.pdf',useDingbats = T,height = 5,width = 19)
plot_grid(plots[[1]],plots[[2]],plots[[3]],plots[[4]],labels = datasets, ncol = 4)
dev.off()

