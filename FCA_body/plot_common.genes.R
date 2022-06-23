
library(Seurat)
library(ggplot2);library(gridExtra);
options(stringsAsFactors = F)

#######################################################################
library(ggplot2)

x=readRDS(paste0('~/Documents/fly.brain_coexpr.net_2022.06/FCA_body/all_common_genes.rds'));
df.ngene=data.frame(cell.cluster=names(x$genes),ngene=sapply(x$genes,length))
head(df.ngene)
df.ngene=df.ngene[order(df.ngene$ngene),]
df.ngene$cell.cluster=factor(df.ngene$cell.cluster,df.ngene$cell.cluster)
p1=ggplot(df.ngene,aes(x=cell.cluster,y=ngene))+
  geom_bar(stat='identity',fill="#56B4E9")+
  geom_text(aes(x=cell.cluster,y=ngene+350,label=ngene))+
  coord_flip()+theme_classic()+
  ylab('Number of expressed genes per cell cluster')+xlab('Cell cluster')
p1

df.ncell=as.data.frame(table((table(unlist(x$genes)))))
colnames(df.ncell)=c('ncell','ngene')
p2=ggplot(df.ncell,aes(x=ncell,y=ngene))+
  geom_bar(stat='identity',fill="#E69F00")+
  geom_text(aes(x=ncell,y=ngene+100,label=ngene))+
  theme_classic()+coord_flip()+
  xlab('Number of expressed cell clusters per gene')+ylab('Number of genes')

df=data.table::fread('FCA.body_ncell_per.cluster.txt')
head(df)
df=df[df$female>=200 | df$male>=200,]
grid.table(df)

gt <- tableGrob(df, rows = NULL, theme = ttheme_default())

pdf("plot_common.genes.pdf",useDingbats = T,width=20,height=12)
grid.arrange(p1,p2,gt,ncol=3)
#grid.arrange(p1,p2,ncol=2)
#plot.new()
#grid.table(df)
dev.off()

