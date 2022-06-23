
## core net
x1=read.table('./brain_scRNA-seq_n15c0.005/cut43_net.txt',header=T);
x2=read.table('./brain_snRNA-seq_n15c0.005/cut35_net.txt',header=T);
head(x1)
head(x2)

gene1=unique(c(x1[,1],x1[,2]))
gene2=unique(c(x2[,1],x2[,2]))

k1=length(gene1) #205
k2=length(gene2) #20

obs=sum(gene2 %in% gene1) #7

## common genes 
x=readRDS('../01.common.gene.extraction/brain_scRNA-seq_n15c0.005/all_common_genes.rds')
common.gene1=x$common.genes
x=readRDS('../01.common.gene.extraction/brain_snRNA-seq_n15c0.005/all_common_genes.rds')
common.gene2=x$common.genes

n1=length(common.gene1); #2088
n2=length(common.gene2); #842

## simu
n.simu=1000;
simu.out=sapply(1:n.simu,function(i.simu){
  p1=sample(common.gene1,k1,replace = F)
  p2=sample(common.gene2,k2,replace = F)
  sum(p1 %in% p2)
})

summary(simu.out)
pvalue=(sum(simu.out>obs)+1)/(n.simu+1)
cat('pvalue=',pvalue)
#pvalue= 0.000999001

df=as.data.frame(table(simu.out))
df$prop=df$Freq/sum(df$Freq)
library(ggplot2)
df$simu.out=as.numeric(as.character(df$simu.out))

pdf('overlap_core.genes.pdf',useDingbats = T)
ggplot(df,aes(x=simu.out,y=prop))+
  geom_bar(stat='identity',fill=NA,col='steelblue',width =1)+
  scale_x_continuous(breaks=0:8,labels=0:8)+
  ylab('Frequency')+xlab('#overlapped core genes')+
  theme_classic(base_size = 20)+
  ggtitle('#overlapped core genes: obs VS simu')+
  geom_vline(xintercept = obs,col='darkred',lwd=2)
dev.off()

#hist(simu.out,xlim=c(0,max(obs,simu.out)),main='#overlapped core genes: obs VS simu',xlab='#overlapped core genes')
#abline(v=obs,col='darkred',lwd=2)

