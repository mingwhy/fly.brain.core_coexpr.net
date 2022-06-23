library(tidyverse)
library(Seurat)
library(ggplot2);library(gridExtra);
options(stringsAsFactors = F)

##################################
#out.dir1='brain_scRNA-seq_n15c0.005/';
out.dir1='brain_snRNA-seq_n15c0.005/';
if(!dir.exists(out.dir1)) dir.create(out.dir1)

## read in processed wholebrain data
file='~/Documents/Data_fly_FCA/xFCA_body/whole_body_filtered_valid.rds';
dat=readRDS(file);

unique(dat@meta.data$sex)
table(dat@meta.data$annotation)
unique(dat@meta.data$annotation) 

rm.tc=c("germline cell", "female reproductive system",
        'male accessory gland','spermatocyte');

cell.meta=dat@meta.data
x=cell.meta[!cell.meta$annotation %in% rm.tc,] %>% group_by(sex,annotation) %>% summarise(n=n())
x1=x %>% spread(sex,n)
sum(x1$female>=200,na.rm=T)+sum(x1$male>=200,na.rm=T) #35
sum(x1$female>=200,na.rm=T) #18
sum(x1$male>=200,na.rm=T) #17
data.table::fwrite(as.data.frame(x1),quote=F,'FCA.body_ncell_per.cluster.txt')

dat=dat[,!dat$annotation %in% rm.tc]
## separate female and male, select cell clusters that contain>=500cells in both sexes
min.cell=200;max.cell=Inf;
i=apply(table(dat$sex,dat$annotation),2,function(i)sum(i>=min.cell & i<max.cell))
#pick.cell.clusters=names(which(i==2)) #16
pick.cell.clusters=names(which(i>=1)) #19
pick.cell.clusters #19

dat=subset(dat,annotation %in% pick.cell.clusters)
df.expr=dat@assays$RNA@counts;
dim(df.expr); #sn:14996 80845

sparsity=sum(df.expr==0)/nrow(df.expr)/ncol(df.expr)
cat(out.dir1,sparsity,'\n') #0.9599326  

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
length(all.genes) #11141
length(common.genes) #869

saveRDS(list(genes=genes,common.genes=common.genes,all.genes=all.genes),
        file=paste0('all_common_genes.rds'))

length(grep('^female',names(genes)))
length(grep('^male',names(genes)))

##################################
library(ggplot2)
x=readRDS(paste0('all_common_genes.rds'));
df.ngene=data.frame(cell.cluster=names(x$genes),ngene=sapply(x$genes,length))
head(df.ngene)
df.ngene=df.ngene[order(df.ngene$ngene),]
df.ngene$cell.cluster=factor(df.ngene$cell.cluster,df.ngene$cell.cluster)
df.ngene=df.ngene[!df.ngene$cell.cluster %in% rm.tc,]
unique(df.ngene$cell.cluster) #35

p1=ggplot(df.ngene,aes(x=cell.cluster,y=ngene))+
  geom_bar(stat='identity',fill="#56B4E9")+
  geom_text(aes(x=cell.cluster,y=ngene+400,label=ngene))+
  coord_flip()+theme_classic()+
  ylab('Number of expressed genes per cell cluster')+xlab('Cell cluster')
p1

df.ncell=as.data.frame(table((table(unlist(x$genes)))))
colnames(df.ncell)=c('ncell','ngene')
p2=ggplot(df.ncell,aes(x=ncell,y=ngene))+
  geom_bar(stat='identity',fill="#E69F00")+
  geom_text(aes(x=ncell,y=ngene+20,label=ngene))+
  theme_classic()+coord_flip()+
  xlab('Number of expressed cell clusters per gene')+ylab('Number of genes')
p2

pdf("FCA.body_ngene_per.cluster_ncluster_per.gene.pdf",useDingbats = T,width=14,height=12)
grid.arrange(p1,p2,ncol=2)
dev.off()
