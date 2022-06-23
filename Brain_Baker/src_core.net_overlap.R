# due to symbol: gamma and \{gamma}
# use FBgn to match data

## common genes 
x=readRDS('../01.common.gene.extraction/brain_scRNA-seq_n15c0.005/all_common_genes.rds')
common.gene1=x$common.genes
tmp=data.table::fread('../../single.cell_datasets/fly.brain.atlas/gene.meta_brain.txt')
head(tmp)
sum(common.gene1 %in% tmp$SYMBOL) #2088
common.gene1=tmp[tmp$SYMBOL %in% common.gene1,]$validated_id

x=readRDS('./brain_scRNA-seq_n15c0.005_all.clusters.005/all_common_genes.rds')
common.gene2=gsub('Dmel-','',x$common.genes)
length(common.gene2) #1738
#write.table(common.gene2,'common.gene_CGgenes.txt',quote=F,row.names = F,col.names = F)

# validate on flybase: common.gene_CGgenes_FlyBase_IDs.txt
x=read.table('common.gene_CGgenes_FlyBase_IDs.txt',header=F)
length(unique(x[,3])) #1748
hit1=x[x[,1] %in% names(which(table(x[,1])==1)),] #1725
hit2=x[!(x[,1] %in% hit1[,1]),]
hit2=hit2[c(2,4,5,8,10,11,13,15,17,19,22,23,25),]

common.gene=rbind(hit1,hit2)
length(unique(common.gene$V3)) #1738
#write.table(common.gene,'common.gene_CGgenes2symbols.txt',quote=F,row.names = F,col.names = F)

common.gene2=common.gene[,2]
n1=length(common.gene1); #2088
n2=length(common.gene2); #1738


## core net
x1=read.table('../07.detect_core/brain_scRNA-seq_n15c0.005/cut43_net.txt',header=T);
#x2=read.table('./brain_scRNA-seq_n15c0.005/cut37_net.txt',header=T);
x2=read.table('./brain_scRNA-seq_n15c0.005/cut33_net.txt',header=T);
head(x1)
head(x2)

gene1=unique(c(x1[,1],x1[,2]))
gene2=unique(c(x2[,1],x2[,2]))
gene1=tmp[tmp$SYMBOL %in% gene1,]$validated_id #use FBgn
gene2=common.gene[common.gene[,3] %in% gene2,][,2] #use FBgn

k1=length(gene1) #205
k2=length(gene2) #88; 65

obs=sum(gene2 %in% gene1) #62; 46


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
  #scale_x_continuous(breaks=seq(1,47,5),labels=seq(1,47,5))+
  scale_x_continuous(breaks=(1:13)*5,labels=(1:13)*5)+
  ylab('Frequency')+xlab('#overlapped core genes')+
  theme_classic(base_size = 20)+
  ggtitle('#overlapped core genes: obs VS simu')+
  geom_vline(xintercept = obs,col='darkred',lwd=2)
dev.off()

#hist(simu.out,xlim=c(0,max(obs,simu.out)),main='#overlapped core genes: obs VS simu',xlab='#overlapped core genes')
#abline(v=obs,col='darkred',lwd=2)

