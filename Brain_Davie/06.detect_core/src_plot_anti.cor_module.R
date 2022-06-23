library(grid)
library(gridExtra)
library(pheatmap)
library(ggplot2)

gene.module=data.table::fread('~/Documents/fly.brain_coexpr.net_2022.06/comparison/cytospace/cut12_net_module_gene.txt')
gene.module=gene.module[!duplicated(gene.module),]
head(gene.module)
dim(gene.module) #1292 genes with a module membership

edge.info=data.table::fread('~/Documents/fly.brain_coexpr.net_2022.06/comparison/cytospace/cut12_net_module.txt')
dim(edge.info) #69463
sum(gene.module$gene %in% unique(c(edge.info$from,edge.info$to))) #1292

df=merge(edge.info,gene.module,by.x='from',by.y='gene')
df=df[,c('from','to','module')]
table(df$module)
sum(df$from=='Ald1');sum(df$to=='Ald1')
#df[df$from=='Ald1',]$from='Ald';
#df[df$to=='Ald1',]$to='Ald';

genes=unique(c(df$from,df$to))
length(genes) #1292

# collect bigscale output nets 
(files=Sys.glob(paste0('../brain_scRNA-seq_n15c0.005/bigScale2_coexpr/*.rds')))
length(files) # 67

#################################################################
## look all core genes
core.net=list()
for(i.file in 1:length(files)){
  #name=basename(dirname(files[i.file]));
  tmp=stringr::str_extract(basename(files[i.file]),'_.+male_.*_');
  tmp=gsub('^_\\d+_','',tmp)
  name=substr(tmp,1,nchar(tmp)-1)
  net=readRDS(files[i.file]);
  Dp=net$Dp.mat;
  #sum(gene.module$gene %in% rownames(Dp)) #125 
  
  x=Dp[gene.module$gene,gene.module$gene]
  core.net[[name]]=x
}
length(core.net)

#https://stackoverflow.com/questions/54140152/pheatmap-change-text-color
#image(core.net[[1]])
e=pheatmap::pheatmap(core.net[[1]],treeheight_row = F,treeheight_col = F,
                   cluster_rows=FALSE, cluster_cols=FALSE,
                   fontsize_row=4,fontsize_col = 4,
                   #show_colnames = F,show_rownames = F,
                   main=names(core.net)[1])
e$gtable$grobs[[3]]$gp=gpar(col=gene.module$color,fontsize=4)
e$gtable$grobs[[4]]$gp=gpar(col=gene.module$color,fontsize=4)
e


plots=list();
for(i in 1:length(core.net)){
#for(i in 1:2){
  e=pheatmap::pheatmap(core.net[[i]],treeheight_row = F,treeheight_col = F,
                       cluster_rows=FALSE, cluster_cols=FALSE,
                       fontsize_row=4,fontsize_col = 4,
                       #show_colnames = F,show_rownames = F,
                       main=names(core.net)[i])
  e$gtable$grobs[[3]]$gp=gpar(col=gene.module$color,fontsize=1)
  e$gtable$grobs[[4]]$gp=gpar(col=gene.module$color,fontsize=1)
  #print(e)
  plots[[i]]<-e[[4]]
}
  
layout_matrix = matrix(1:16,4,4,byrow = T)


#png('test.png',res=90,units='in',width = 22,height = 20)
#grid.arrange(grobs= plots)
#dev.off()

if(F){
  # folder: module.heatmap_cell.cluster
  pdf('anti.cor_module.pdf',useDingbats = T,width = 22,height = 20)
  grid.arrange(arrangeGrob(grobs= plots[1:16],layout_matrix = layout_matrix))
  grid.arrange(arrangeGrob(grobs= plots[17:32],layout_matrix = layout_matrix))
  grid.arrange(arrangeGrob(grobs= plots[33:48],layout_matrix = layout_matrix))
  grid.arrange(arrangeGrob(grobs= plots[49:64],layout_matrix = layout_matrix))
  grid.arrange(arrangeGrob(grobs= plots[65:67],layout_matrix = layout_matrix))
  dev.off()
}

png('anti.cor_module%02d.png',res=90,units='in',width = 22,height = 20)
grid.arrange(arrangeGrob(grobs= plots[1:16],layout_matrix = layout_matrix))
grid.arrange(arrangeGrob(grobs= plots[17:32],layout_matrix = layout_matrix))
grid.arrange(arrangeGrob(grobs= plots[33:48],layout_matrix = layout_matrix))
grid.arrange(arrangeGrob(grobs= plots[49:64],layout_matrix = layout_matrix))
grid.arrange(arrangeGrob(grobs= plots[65:67],layout_matrix = layout_matrix))
dev.off()

#################################################################
## only look at modules 
# collect core.edge.cor.values in each gene co-expression network
if(F){
module.net=list()
for(i.file in 1:length(files)){
  name=basename(dirname(files[i.file]));
  net=readRDS(files[i.file]);
  Dp=net$Dp.mat;
  sum(genes %in% rownames(Dp)) #125 
  
  cor.values=sapply(1:nrow(df),function(j){
    pos=match(df[j,c(1,2)],rownames(Dp))
    Dp[pos[1],pos[2]]
  })
  module.net[[name]]=cbind(df,cor.values)
}
for(i in 1:length(module.net)){
  dfc=module.net[[i]]
  dfc$module=factor(dfc$module);
  ggplot(dfc,aes(y=cor.values,fill=module,col=module))+
    geom_scatter()+
    #geom_boxplot(outlier.shape = NA)+
    theme_classic(base_size = 20)+
    ggtitle(names(module.net)[i])
}
}

