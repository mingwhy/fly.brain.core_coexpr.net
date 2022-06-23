
library(ggplot2);library(gridExtra);
options(stringsAsFactors = F)

datasets=c('brain, Davie et al. 2018','brain, Baker et al. 2021',
           'head, Li et al. 2022', 'body, Li et al. 2022')

## common expr gene plot
path=c('~/Documents/fly.brain_coexpr.net_2022.06/Brain_Davie/',
       '~/Documents/fly.brain_coexpr.net_2022.06/Brain_Baker/',
       '~/Documents/fly.brain_coexpr.net_2022.06/FCA_head/',
       '~/Documents/fly.brain_coexpr.net_2022.06/FCA_body/');

plots1<-lapply(path,function(i){
  x=readRDS(paste0(i,'/all_common_genes.rds'));
  df.ngene=data.frame(cell.cluster=names(x$genes),ngene=sapply(x$genes,length))
  head(df.ngene)
  df.ngene=df.ngene[order(df.ngene$ngene),]
  df.ngene$cell.cluster=factor(df.ngene$cell.cluster,df.ngene$cell.cluster)
  p1=ggplot(df.ngene,aes(x=cell.cluster,y=ngene))+
    geom_bar(stat='identity',fill="#56B4E9")+
    geom_text(aes(x=cell.cluster,y=ngene+350,label=ngene))+
    coord_flip()+theme_classic()+
    ylab('Number of expressed genes per cell cluster')+xlab('Cell cluster')
  
  
  df.ncell=as.data.frame(table((table(unlist(x$genes)))))
  colnames(df.ncell)=c('ncell','ngene')
  p2=ggplot(df.ncell,aes(x=ncell,y=ngene))+
    geom_bar(stat='identity',fill="#E69F00")+
    geom_text(aes(x=ncell,y=ngene+100,label=ngene))+
    theme_classic()+coord_flip()+
    xlab('Number of expressed cell clusters per gene')+ylab('Number of genes')
  
  file=Sys.glob(paste0(i,'*_ncell_per.cluster.txt'))
  print(file)
  df=data.table::fread(file)
  df=df[df$female>=200 | df$male>=200,]
  p3 <- tableGrob(df, rows = NULL, theme = ttheme_default())
  list(p1,p2,p3)
})

## sparsity plot
path2=c('~/Documents/fly.brain_coexpr.net_2022.06/Brain_Davie/02_choose_top.percentage_threshold/',
       '~/Documents/fly.brain_coexpr.net_2022.06/Brain_Baker/02_choose_top.percentage_threshold/',
       '~/Documents/fly.brain_coexpr.net_2022.06/FCA_head/02_choose_top.percentage_threshold/',
       '~/Documents/fly.brain_coexpr.net_2022.06/FCA_body/02_choose_top.percentage_threshold/');
plots2<-lapply(path2,function(path.i){
  file=Sys.glob(paste0(path.i,'./s*_sparsity.rds'))
  print(file)
  df.out.sparse=readRDS(file)
  df.out.sparse$cluster=factor(df.out.sparse$cluster,df.out.sparse$cluster)
  max(df.out.sparse$sparsity) #0.5263604
  
   ggplot(df.out.sparse,aes(x=cluster,y=sparsity))+
           geom_bar(stat='identity',fill=NA,col='darkblue')+
           scale_y_continuous(breaks=c(0,0.2,0.4,0.5,0.55),labels=c(0,0.2,0.4,0.5,0.55))+
           theme_classic()+
           #theme(
             #axis.title=element_text(size=16),
             #axis.text.y=element_text(size=15),
             #axis.text.x=element_text(size=17,angle=0,vjust=1,hjust=1))+
           coord_flip()+
           xlab('')
})

labels=c('A','B','C','D')
library(cowplot)

pdf('plot_common.expr.gene_sparsity.pdf',useDingbats = T,height = 12,width = 26)
for(i in 1:4){
  print( plot_grid(plots1[[i]][[1]],plots1[[i]][[2]],plots1[[i]][[3]],plots2[[i]], 
                #labels = labels[i], 
                ncol = 4,rel_widths = c(0.25,0.25,0.18,0.36))
   )
}
dev.off()

pdf('plot_sparsity.pdf',useDingbats = T,height = 12,width = 26)
plot_grid(plots2[[1]],plots2[[2]],plots2[[3]],plots2[[4]],labels = '', ncol = 4)
dev.off()


