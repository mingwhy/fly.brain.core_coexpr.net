
# plot the gene and edge commonality distribution 
# at different topN cutoffs

library(Seurat)
library(Matrix);library(tidyverse)
library(ggplot2);library(gridExtra)
library(RColorBrewer)
input.folder='./brain_snRNA-seq_n15c0.005_bigScale2/';
output.folder='brain_snRNA-seq_n15c0.005/';
dir.create(output.folder)

# read in networks at different cutoffs to get edge commonality distribution
(files=Sys.glob(paste0(input.folder,'/top_*_networks_cutoff_cell.type.rds')))

common.edges=list();common.genes=list()
for(file in files){
  (cutoff=strsplit(basename(file),'\\_')[[1]][[2]])
  
  (outfile=paste0(output.folder,'/top_',cutoff,'_commonality.distribution.rds'));
  cat('processing file',file,'\n');
  
  if(!file.exists(outfile)){
    keep.list=readRDS(file); #raed in binary networks
    
    pan<-Reduce(`+`,keep.list)
  
    edge.common=pan[upper.tri(pan)]
    edge.common.distr=table(edge.common)
    
    x=sapply(keep.list,function(x){Matrix::rowSums(x)})
    gene.occur=apply(x,1,function(x){sum(x>1)}) #keep.list[[1]][1:2,1:2];#diag=1
    gene.common.distr=table(gene.occur)
    
    saveRDS(list(gene.common.distr=gene.common.distr,
                 edge.common.distr=edge.common.distr),outfile)
  }
}


common.edges=list();common.genes=list();
(files=Sys.glob(paste0(output.folder,'/top*commonality.distribution.rds')))
for(file in files){
  cutoff=strsplit(basename(file),'\\_')[[1]][2]
  result=readRDS(file)
  common.edges[[as.character(cutoff)]]=result$edge.common.distr
  common.genes[[as.character(cutoff)]]=result$gene.common.distr
}


# prepare data for plotting gene and edge commonality distribution
common.genes;
common.edges
common.df=list();
common.df=lapply(list(common.genes,common.edges),function(common.set){
  cutoff.list=names(common.set);
  cutoff=rep(cutoff.list,sapply(common.set,length))
  commonality=as.numeric(as.character(unlist(sapply(common.set,names))))
  sapply(common.set,sum) #all the same, as 0 also shows up in the group
  df=data.frame(cutoff=cutoff,
                commonality=commonality,
                count=unlist(common.set))
  df$commonality
  df1=df[df$commonality!=0,]
  df2 = df1 %>% group_by(cutoff) %>% 
    #summarise(total.link=sum(count),prop=count/sum(count))  
    mutate(total.link=sum(count),prop=count/sum(count))   
  sum(df2$prop)
  df2;
})


# begin plot
names(common.genes) #cutoff values
pick.cutoff='0.05'
pick.gene=subset(common.df[[1]],cutoff==pick.cutoff)
subset(pick.gene,commonality==35)$count #82
pick.edge=subset(common.df[[2]],cutoff==pick.cutoff & commonality!=0)
sum(subset(pick.edge,commonality<=4)$prop)#0.9050484
sum(subset(pick.edge,commonality>=15)$prop)#0.001664504
max(pick.edge$commonality) #27
sum(subset(pick.edge,commonality==27)$count) #1

pdf(paste0(output.folder,'/pan_n_gene_edge_commonality.pdf'),width=18,height = 12,useDingbats = T)
for(pick.cutoff in names(common.genes)){
  
  pick.gene=subset(common.df[[1]],cutoff==pick.cutoff)
  pick.edge=subset(common.df[[2]],cutoff==pick.cutoff & commonality!=0)
   
  my.level=sort(unique(c(pick.gene$commonality,pick.edge$commonality)))
  
  nedge=sum(pick.edge$count)
  if(abs(sum(pick.edge$prop)-1)>1e5){cat(pick.cutoff,"edge.prop didn't sum up to 1\n")}
  
  pick.edge[pick.edge$commonality==1,] #prop.of.edges with commonality=1
  sum(pick.edge[pick.edge$commonality>10,]$prop) # 0.005342755
  
  tmp=data.frame(commonality=my.level)
  x=merge(tmp,pick.edge,all=T)
  x$commonality=factor(x$commonality,levels=my.level)
  p1<-ggplot(x,aes(x=commonality,y=prop,fill=factor(cutoff)))+
    xlab('Edge commonality')+ylab('Frequency of edges')+
    scale_x_discrete(expand = expansion(mult = c(0.02, 0.01)),
                     breaks=seq(1,35,2),
                     labels=seq(1,35,2))+
    scale_y_continuous(trans = scales::pseudo_log_trans(1e-05, 10),
                       breaks=c(0, 1e-05,1e-04,1e-03,0.01, 0.1, 1))+
    geom_text(label=x$count,hjust=0.4,vjust=-0.2,size=4)+
    #geom_text(label=paste0(tmp$count,'\n',round(tmp$prop,4)),
    #          hjust=0.4,vjust=-0.2)+
    scale_fill_manual(name='',values=c("#E69F00"))+
    geom_bar(stat='identity',position='dodge')+
    #scale_fill_manual(values=mycol)+
    theme_bw(base_size = 14)+
    #ggtitle(paste0(nedge,' edges in the union of all cell type-specific networks\nPCC top',pick.cutoff))+
    theme(axis.text = element_text(size=14),
          axis.text.x = element_text(size=14,angle=0,hjust=0.5),
          axis.title = element_text(size=20,face="bold"),
          legend.position = 'none',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    )
    
  ngene=sum(pick.gene$count)
  if(sum(pick.gene$prop)!=1){cat("gene.prop didn't sum up to 1\n")}
  pick.gene$prop=round(pick.gene$prop,3)
  tmp=data.frame(commonality=my.level)
  x=merge(tmp,pick.gene,all=T)
  x$commonality=factor(x$commonality,levels=my.level)
  #p2<-ggplot(tmp,aes(x=factor(commonality),y=count,fill=factor(cutoff)))+

  p2<-ggplot(x,aes(x=commonality,y=prop,fill=factor(cutoff)))+
    xlab('Gene commonality')+ylab('Frequency of genes')+
    #ylim(0,7.2e05)+
    #ylim(0,max(tmp$count)*1.2)+
    scale_x_discrete(expand = expansion(mult = c(0.02, 0.01)),
                     breaks=seq(1,35,2),
                     labels=seq(1,35,2))+
    scale_y_continuous(trans = scales::pseudo_log_trans(1e-05, 10),
                       breaks=c(0, 1e-05,1e-04,1e-03,0.01, 0.1, 1))+
    #geom_text(label=paste0(tmp$count,'\n',round(tmp$prop,6)),
    #          hjust=0.4,vjust=-0.2)+
    geom_text(label=x$count,hjust=0.4,vjust=-0.2,size=4)+
    scale_fill_manual(name='',values=c("#56B4E9"))+
    geom_bar(stat='identity',position='dodge')+
    #scale_fill_manual(values=mycol)+
    theme_bw(base_size = 14)+
    #ggtitle(paste0(ngene,' genes in the union of all cell type-specific networks\nPCC top',pick.cutoff))+
    theme(axis.text = element_text(size=14),
          axis.text.x = element_text(size=14,angle = 0,hjust=0.5),
          axis.title = element_text(size=20,face="bold"),
          legend.position = 'none',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    )
  
    #grid.arrange(plots[[1]],plots[[2]],ncol=2)
    grid.arrange(p1,p2,ncol=1)
    plot.new()
    grid.table(head(pick.gene,20))
    plot.new()
    grid.table(head(pick.edge,20))
}
dev.off();

