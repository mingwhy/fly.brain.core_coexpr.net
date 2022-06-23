
## core net
x1=read.table('~/Documents/fly.brain_coexpr.net_2022.02/2022_02_fly.brain.core_coexpr.net/07.detect_core/brain_scRNA-seq_n15c0.005/cut43_net.txt',header=T);
x2=read.table('./brain_snRNA-seq_n15c0.005/cut12_net.txt',header=T);
head(x1)
head(x2)
dim(x1) #2140
dim(x2) #1894

# #overlap edge
edges1=paste(x1$from,x1$to) #brain
edges2=paste(x2$from,x2$to) #body
sum(edges1 %in% edges2) #572
overlap.edges=edges1[edges1 %in% edges2] #all RpS, "DnaJ-1 Hsp83"

gene1=unique(c(x1[,1],x1[,2]))
gene2=unique(c(x2[,1],x2[,2]))

(k1=length(gene1)) #205, brain
(k2=length(gene2)) #206, body
overlap.genes=gene2[gene2 %in% gene1]

write.table(overlap.edges,'overlap.edges.txt',quote=F,row.names = F)
write.table(overlap.genes,'overlap.genes.txt',quote=F,row.names = F)

(obs=sum(gene2 %in% gene1)) #74

########################################################################
## common genes 
x=readRDS('~/Documents/fly.brain_coexpr.net_2022.02/2022_02_fly.brain.core_coexpr.net/01.common.gene.extraction/brain_scRNA-seq_n15c0.005/all_common_genes.rds')
common.gene1=x$common.genes
x=readRDS('./brain_snRNA-seq_n15c0.005/all_common_genes.rds')
common.gene2=x$common.genes

(n1=length(common.gene1)); #2088
(n2=length(common.gene2)); #869

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

pdf('overlap_core.genes.pdf',useDingbats = T,height = 4)
ggplot(df,aes(x=simu.out,y=prop))+
  geom_bar(stat='identity',fill=NA,col='steelblue',width =1)+
  scale_x_continuous(breaks=seq(0,76,4),labels=seq(0,76,4))+
  ylab('Frequency')+xlab('#overlapped core genes')+
  theme_classic(base_size = 12)+
  ggtitle('#overlapped core genes: obs VS simu')+
  geom_vline(xintercept = obs,col='darkred',lwd=2)
dev.off()

#hist(simu.out,xlim=c(0,max(obs,simu.out)),main='#overlapped core genes: obs VS simu',xlab='#overlapped core genes')
#abline(v=obs,col='darkred',lwd=2)

