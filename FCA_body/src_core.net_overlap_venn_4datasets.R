
library(VennDiagram)
## core net
x1=read.table('~/Documents/fly.brain_coexpr.net_2022.02/2022_02_fly.brain.core_coexpr.net/07.detect_core/brain_scRNA-seq_n15c0.005/cut43_net.txt',header=T);
x2=read.table('./brain_snRNA-seq_n15c0.005/cut12_net.txt',header=T);
x3=read.table('~/Documents/fly.brain_coexpr.net_2022.02/2022_02_fly.brain.core_coexpr.net/cocaine_fly.brain/brain_scRNA-seq_n15c0.005/cut33_net.txt',header=T);
x4=read.table('~/Documents/fly.brain_coexpr.net_2022.02/2022_02_fly.brain.core_coexpr.net/07.detect_core/brain_snRNA-seq_n15c0.005/cut35_net.txt',header = T)
head(x1)
dim(x1) #2140
dim(x2) #1894
dim(x3) #323
dim(x4) #22

# #overlap edge
edges1=paste(x1$from,x1$to) #brain
edges2=paste(x2$from,x2$to) #body, FCA
edges3=paste(x3$from,x3$to) #cocanie brain
edges4=paste(x4$from,x4$to) #head, FCA
sum(edges1 %in% edges2) #572
overlap.edges=edges1[edges1 %in% edges2] #all RpS, "DnaJ-1 Hsp83"

gene1=unique(c(x1[,1],x1[,2]))
gene2=unique(c(x2[,1],x2[,2]))
gene3=unique(c(x3[,1],x3[,2]))
gene4=unique(c(x4[,1],x4[,2]))

(k1=length(gene1)) #205, brain
(k2=length(gene2)) #206, body
sum(gene1 %in% gene2) #74
overlap.genes=gene2[gene2 %in% gene1]

library(gplots)
pdf('core.net_overlap_venn_4datasets.pdf',useDingbats = T,height = 11)
par(mfrow=c(2,1),mar=c(0,0,0,0))
my.gene.sets<-list(fly.brain=gene1,FCA.body=gene2,cocaine.brain=gene3,FCA.head=gene4)
v.genes<-venn(my.gene.sets)
text(200,400,'Venn diagram of core genes',cex=2)
#title='Venn diagram of core genes'

my.edge.sets<-list(fly.brain=edges1,FCA.body=edges2,cocaine.brain=edges3,FCA.head=edges4)
v.edges<-venn(my.edge.sets)
text(200,400,'Venn diagram of core edges',cex=2)
dev.off()

dim(v.genes) #16 combos
v.genes[1:2,]
x=attr(v.genes,"intersections")
names(x)
x[grep(".+:.+:.+",names(x))] #three

y=x=attr(v.edges,"intersections")
y[grep(".+:.+:.+",names(y))] #three

shared.3.edges=y[grep(".+:.+:.+",names(y))][[1]]
df=as.data.frame(Reduce(`rbind`,sapply(shared.3.edges,strsplit,' ')))
head(df)       
colnames(df)=c('from','to')
write.table(df,file='overlap.genes.edges.4datasets.txt',quote=F,row.names = F,sep='\t')

unique(c(df$from,df$to)) #13 genes
write.table(unique(c(df$from,df$to)),'overlap.genes.4datasets.txt',quote=F,row.names = F,col.names = F,sep='\t')

gene.meta=data.table::fread('~/Documents/Data_fly_FCA/gene.meta.txt',header=T,sep='\t')
head(gene.meta)
tmp=gene.meta[gene.meta$current_symbol %in% unique(c(df$from,df$to)),]
colnames(tmp)
tmp[,c('SYMBOL','NAME')]
library(ggplot2);library(gridExtra);library(grid)

pdf('overlap.genes.4datasets.pdf',useDingbats = T)
plot.new()
grid.table(tmp[,c('SYMBOL','NAME')])
dev.off()

 