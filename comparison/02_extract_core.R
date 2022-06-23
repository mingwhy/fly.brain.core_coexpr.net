
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

if(T){
  # validate on flybase: common.gene_CGgenes_FlyBase_IDs.txt
  gene.meta=read.table('~/Documents/fly.brain_coexpr.net_2022.06/Brain_Baker/common.gene_CGgenes_FlyBase_IDs.txt',header=F)
  colnames(gene.meta)<-c('CG','FBgn','symbol')
  
  cutoffs=c(12,10,6,11)
  files=c('~/Documents/fly.brain_coexpr.net_2022.06/Brain_Davie/brain_scRNA-seq_n15c0.005_bigScale2/top_0.05_pan.rds',
          '~/Documents/fly.brain_coexpr.net_2022.06/Brain_Baker/brain_scRNA-seq_n15c0.005_bigScale2/top_0.05_pan.rds',
          '~/Documents/fly.brain_coexpr.net_2022.06/FCA_head/brain_snRNA-seq_n15c0.005_bigScale2/top_0.05_pan.rds',
          '~/Documents/fly.brain_coexpr.net_2022.06/FCA_body/brain_snRNA-seq_n15c0.005_bigScale2/top_0.05_pan.rds')
  
  cores<-lapply(1:4,function(ii){
    cut=cutoffs[ii]
    pan=readRDS(files[ii])
    diag(pan)=0;pan[is.na(pan)]=0;dim(pan)
    i=Matrix::rowSums(pan,na.rm=T)
    sum(i==0)
    if(sum(i==0)>0){pan=pan[i!=0,i!=0]}
    pan.genes=names(which(i>0))
    length(pan.genes)
    if(length(grep('Dmel-',pan.genes[1]))>0){
      CG.names=gsub('Dmel-','',pan.genes)
      tmp=match(CG.names,gene.meta$CG)
      gene.meta=gene.meta[tmp,]
      if(sum(gene.meta$CG==CG.names)!=length(pan.genes)){cat('gene.name match not right\n');break}
      rownames(pan)=colnames(pan)=gene.meta$symbol
    }
    net=pan>=cut; #True or False binary matrix
    net0=net; #unweight
    diag(net0)=FALSE;
    i=apply(net0,1,sum)
    sum(i==0);length(i) #869 common expr genes in total
    net1=net0[i!=0,i!=0]
    cat('nlink.cut',cut,',ngene',nrow(net1),',nedge',sum(net1)/2,'\n')
    net1
  })
  #nlink.cut 12 ,ngene 1428 ,nedge 69664 
  #nlink.cut 10 ,ngene 1171 ,nedge 39076 
  #nlink.cut 6 ,ngene 630 ,nedge 29413 
  #nlink.cut 11 ,ngene 244 ,nedge 2357 
  edges.df<-lapply(cores,function(net1){
    # use gene symbols in graph
    g<-igraph::graph_from_adjacency_matrix(net1,mode='undirected',weighted = NULL,diag=F) 
    df.g<-igraph::as_data_frame(g)
    df.g$edge=1
    df.g
  })
  edges<-lapply(edges.df,function(x){
    apply(x[,c(1,2)],1,function(j) paste(sort(as.character(j)),collapse = ';')) 
  })
  
  names(edges)=datasets
  sapply(edges,length)
  saveRDS(list(cores=cores,edges=edges,edges.df=edges.df),'extract_core.rds')
}
x=readRDS('extract_core.rds')
cores=x$cores
edges=x$edges
edges.df=x$edges.df
v.genes<-venn(edges)
text(200,400,'Venn diagram of core edges',cex=2)

for(i in 1:4){
  tmp=edges.df[[i]]
  data.table::fwrite(tmp,paste0(datasets[[i]],'.txt'))
}

#UpSet plot: https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html
m=make_comb_mat(edges)
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

pdf('intersect_core.pdf',useDingbats = T,width = 9,height = 5)
print(p1)
dev.off()

core.edges=extract_comb(m, "1111")
df=as.data.frame(Reduce(`rbind`,sapply(core.edges,strsplit,';')))
dim(df) #773
head(df)
colnames(df)=c('from','to')
data.table::fwrite(df,'shared_edges.txt')

x=data.frame(shared.genes=unique(c(df$from,df$to)),'shared'=1)
data.table::fwrite(x,'shared_genes.txt')

##############################################################
## GO
if(F){
source('src_fly_go_enrichment.R')
library(org.Dm.eg.db,verbose=F,quietly=T)
shared.genes=data.table::fread('shared_genes.txt')
shared.genes=shared.genes$shared.genes
length(shared.genes) #64 genes
if(!file.exists('shared_genes_GO.rds')){
  out=GOenrich(shared.genes,category='BP',cutoff=0.6) #larger cutoff, smaller returned GO terms
  saveRDS(out,'shared_genes_GO.rds')
}
out=readRDS('shared_genes_GO.rds')
x1=out;
x1$GeneRatio1=x1$GeneRatio
x1$GeneRatio=sapply(x1$GeneRatio,function(x){
  p=as.numeric(unlist(strsplit(x,'/')))
  p[1]/p[2]
})
summary(x1$GeneRatio)
sum(is.na(x1$Description))
x1=x1[!is.na(x1$Description),]

x1=x1[order(x1$p.adjust),]
#x1=x1[x1$Count>=5,]
x1$desp=factor(x1$Description,levels=rev(x1$Description))

pick.row.cols=c( "#D16103","#FFDB6D","#4E84C4","#52854C")
x1$log.p.adjust=-1*log(x1$p.adjust,base=10)
x2=x1;
x2$desp=as.character(x2$desp)
x2$x.axis=x2$desp
if(nrow(x2)>=5){
  x2=x2[1:5,];
  x2$x.axis=factor(x2$x.axis,levels=rev(x2$x.axis))
}else{
  x2$x.axis=factor(x2$x.axis,levels=rev(x2$x.axis))
}
p2=ggplot(x2,aes(x=x.axis,y=log.p.adjust))+
  #geom_bar(fill=adjustcolor('grey10',alpha.f =0.2),stat='identity',width=0.9)+
  geom_bar(fill=pick.row.cols[1],stat='identity',width=0.05)+
  geom_point(col=pick.row.cols[1],size=6)+
  coord_flip()+
  geom_text(aes(label=desp),y=max(x2$log.p.adjust*0.01),
            vjust=-0.4,hjust=0,size=9)+
  theme_bw(base_size=20)+ylab('-log10(p.adjust)')+xlab('')+
  ggtitle(paste0('shared core: ',length(shared.genes),' genes in ',nrow(df),' edges'))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)),limits = c(0, NA))+
  theme(legend.position = 'none',
        axis.title = element_text(size=22),
        axis.text.y=element_blank(),
        axis.text.x=element_text(size=20),
        plot.title = element_text(size = 20, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
p2

pdf("shared_genes_GO.pdf",useDingbats = T,width = 10,height = 6)
print(p2)
dev.off()

}


