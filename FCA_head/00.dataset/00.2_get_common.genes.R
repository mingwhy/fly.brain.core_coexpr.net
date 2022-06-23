
library(Seurat)
library(ggplot2);library(gridExtra);
options(stringsAsFactors = F)

#out.dir1='brain_scRNA-seq_n15c0.005/';
out.dir1='brain_snRNA-seq_n15c0.005/';
if(!dir.exists(out.dir1)) dir.create(out.dir1)

## read in processed wholebrain data
file='~/Documents/Data_fly_FCA/FCA_head/whole_head_filtered_valid.rds';
dat=readRDS(file);

unique(dat@meta.data$sex)
table(dat@meta.data$annotation)

## ncell per sex per cell cluster
library(tidyverse)
tc=unique(dat$annotation)
tc.number<-tc[grep('^\\d+$',tc)] #41
tc.name<-tc[!tc %in% tc.number] #74

cell.meta=dat@meta.data
x=cell.meta %>% group_by(sex,annotation) %>% summarise(n=n())
x1=x %>% spread(sex,n)
sum(x1$female>=200,na.rm=T)+sum(x1$male>=200,na.rm=T) #68
sum(x1$female>=200,na.rm=T) #37
sum(x1$male>=200,na.rm=T) #31
sum(x1[x1$female>=200,]$annotation %in% tc.name)#17 cell types
sum(x1[x1$male>=200,]$annotation %in% tc.name)  #11 cell types
data.table::fwrite(as.data.frame(x1),quote=F,'brain_ncell_per.cluster.txt')

## separate female and male, select cell clusters that contain>=500cells in both sexes
min.cell=200;max.cell=Inf;
i=apply(table(dat$sex,dat$annotation),2,function(i)sum(i>=min.cell & i<max.cell))
#pick.cell.clusters=names(which(i==2))
pick.cell.clusters=names(which(i>=1)) 
pick.cell.clusters #min200: sn,40

dat=subset(dat,annotation %in% pick.cell.clusters)
df.expr=dat@assays$RNA@counts;
dim(df.expr); #sn:12602 46619

sparsity=sum(df.expr==0)/nrow(df.expr)/ncol(df.expr)
cat(sparsity,'\n') 
#min200: sc=0.9051, sn=0.9544

genes=list();
for(i.sex in c('male','female')){
  for(i.cluster in pick.cell.clusters){
    mat=df.expr[,dat$sex==i.sex & dat$annotation==i.cluster]
    if(ncol(mat)<min.cell){next}
    gene.filter=Matrix::rowSums(mat>0) >= max(15,ncol(mat)*0.005)
    
    mat=mat[gene.filter,]
    genes[[paste0(i.sex,'_',i.cluster)]]=rownames(mat)
  }
}

c.genes=unlist(genes)
common.genes=names(which(table(c.genes)==length(genes)))
all.genes=unique(c.genes)
length(all.genes) #8867
length(common.genes) #842

saveRDS(list(genes=genes,common.genes=common.genes,all.genes=all.genes),
        file=paste0('all_common_genes.rds'))

length(grep('female',names(genes)))
length(grep('^male',names(genes)))


library(ggplot2)
x=readRDS(paste0(out.dir1,'all_common_genes.rds'));
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

#pdf("sc_ngene_per.cluster_ncluster_per.gene.pdf",useDingbats = T,width=12,height=12)
pdf("sn_ngene_per.cluster_ncluster_per.gene.pdf",useDingbats = T,width=14,height=12)
grid.arrange(p1,p2,ncol=2)
dev.off()
