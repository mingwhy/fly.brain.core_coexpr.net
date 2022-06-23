
# locate the core, find modules and perform GO enrichment
options(stringsAsFactors = F)
library(igraph)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(RColorBrewer)
library(readxl)
library(org.Dm.eg.db,verbose=F,quietly=T)
library(gplots)
library(UpSetR)
library(ComplexHeatmap)

(mycol=c("#C4961A",'#E7298A',"#52854C", "#4E84C4")) #"#D16103" red
barplot(1:4,col=mycol)
datasets=c('brain, Davie et al. 2018','brain, Baker et al. 2021',
           'head, Li et al. 2022','body, Li et al. 2022')

# validate on flybase: common.gene_CGgenes_FlyBase_IDs.txt
gene.meta=read.table('~/Documents/fly.brain_coexpr.net_2022.06/Brain_Baker/common.gene_CGgenes_FlyBase_IDs.txt',header=F)
colnames(gene.meta)<-c('CG','FBgn','symbol')

## common expr gene plot
path=c('~/Documents/fly.brain_coexpr.net_2022.06/Brain_Davie/',
       '~/Documents/fly.brain_coexpr.net_2022.06/Brain_Baker/',
       '~/Documents/fly.brain_coexpr.net_2022.06/FCA_head/',
       '~/Documents/fly.brain_coexpr.net_2022.06/FCA_body/');

genes<-lapply(path,function(i){
  x=readRDS(paste0(i,'/all_common_genes.rds'));
  common.genes=x$common.genes
  if(length(grep('Dmel-',common.genes[1]))>0){
    CG.names=gsub('Dmel-','',common.genes)
    tmp=match(CG.names,gene.meta$CG)
    gene.meta2=gene.meta[tmp,]
    if(sum(gene.meta2$CG==CG.names)!=length(common.genes)){cat('gene.name match not right\n');break}
    gene.meta2$CG
  }else{
    common.genes
  }
})

v.genes<-venn(genes)

names(genes)=datasets
#UpSet plot: https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html
m=make_comb_mat(genes)
m
set_name(m)
comb_name(m)
set_size(m)  
comb_size(m)

p1<-UpSet(m, comb_col = "#0000FF", bg_col = c("#F0F0FF", "#FFF0F0"), bg_pt_col = "#CCCCFF",
          top_annotation = upset_top_annotation(m, add_numbers = TRUE,
                                                gp = gpar(col = comb_degree(m),fontsize = 20),height = unit(5, "cm")),
          right_annotation = upset_right_annotation(m, add_numbers = TRUE,
                                                    gp = gpar(fill = "white"),
                                                    annotation_name_side = "top",
                                                    axis_param = list(side = "top"),
                                                    width = unit(5, "cm"))
)
p1

pdf('intersect_genes.pdf',useDingbats = T,width = 9,height = 5)
print(p1)
dev.off()
