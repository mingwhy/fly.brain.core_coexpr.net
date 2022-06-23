
# locate the core, find modules and perform GO enrichment
options(stringsAsFactors = F)
library(igraph)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(RColorBrewer)
library(readxl)
library(org.Dm.eg.db,verbose=F,quietly=T)

source('./src_fly_go_enrichment.R')
path='../03.top.perc_net/brain_scRNA-seq_n15c0.005_bigScale2/';
out.dir='brain_scRNA-seq_n15c0.005/';
#path='../03.top.perc_net/brain_snRNA-seq_n15c0.005_bigScale2/';
#out.dir='brain_snRNA-seq_n15c0.005/';
dir.create(out.dir)
topK='top_0.05';
cut=43;

#####################################################################################
## decompose core modules, GO enrichment, get module gene members age distributions
library(MCL)  
library(simplifyEnrichment)
library(ComplexHeatmap)
library(circlize) #for colorRamp2
#source('./src_make_transparent_colors_in_r.R')

(files=Sys.glob(paste0(path,'/top_*_pan.rds')))
(file=files[grep(topK,files)])

(cutoff=strsplit(basename(file),'\\_')[[1]][[2]])

pan=readRDS(file)
diag(pan)=0;pan[is.na(pan)]=0;dim(pan)
i=Matrix::rowSums(pan,na.rm=T)
sum(i==0)
if(sum(i==0)>0){pan=pan[i!=0,i!=0]}
pan.genes=names(which(i>0))
length(pan.genes)

net=pan>=cut; #True or False binary matrix
net0=net; #unweight
diag(net0)=FALSE;
i=apply(net0,1,sum)
sum(i==0);length(i)
net1=net0[i!=0,i!=0]
cat('cutoff',cutoff,'nlink.cut',cut,',ngene',nrow(net1),',nedge',sum(net1)/2,'\n')
#sc, cutoff 0.05 nlink.cut 43 ,ngene 205 ,nedge 2140 
#sn, cutoff 0.05 nlink.cut 35 ,ngene 20 ,nedge 21 

# begin module detection
g<-igraph::graph_from_adjacency_matrix(net1,mode='undirected',weighted = NULL,diag=F) 
df.g<-igraph::as_data_frame(g)
df.g$edge=1
#data.table::fwrite(df.g,file=paste0(out.dir,'cut',cut,'_net.txt'),quote=F,sep='\t',row.names = F)

#pdf(paste0(out.dir,'cut',cut,'_net.pdf'),useDingbats = T)
#plot(g)
#dev.off();

## https://rdrr.io/cran/MCL/man/mcl.html
wc=mcl(x = net1, addLoops = TRUE, allow1 = F,inflation = 2,max.iter = 5000); #sc
#wc=mcl(x = net1, addLoops = TRUE, allow1 = F,inflation = 3,max.iter = 5000); #sn
table(wc$Cluster)
sum(table(wc$Cluster)[names(which(table(wc$Cluster)>=5))])
pick.module.id=names(which(table(wc$Cluster)>=5))
pick.module.id

keep.genes=sapply(pick.module.id,function(i){
  rownames(net1)[wc$Cluster==i]
})
genes.order=unlist(keep.genes)
gene.module=rep(1:length(keep.genes),sapply(keep.genes,length))
mat=net1[genes.order,genes.order] #smaller than net1
gene.module.info=data.frame(gene=genes.order,module=gene.module)
#data.table::fwrite(gene.module.info,file=paste0(out.dir,'cut',cut,'_gene_module.txt'),quote=F,sep='\t',row.names = F)

dim(mat) #178 by 178 genes, genes with module membership
mat[mat==T]=1
mat[mat==F]=0

g<-igraph::graph_from_adjacency_matrix(mat,mode='undirected',weighted = NULL,diag=F) 
df.g<-igraph::as_data_frame(g)
df.g$edge=1
#data.table::fwrite(df.g,file=paste0(out.dir,'cut',cut,'_net_module.txt'),quote=F,sep='\t',row.names = F)

# simpliest plot
#barplot(1:8,col=brewer.pal(8,"Paired"))
#barplot(1:4,col=brewer.pal(8,"Paired")[c(6,8,4,2)])
colors = structure(c('white','#52854C'),names = c('0', '1'))
col_func=colorRamp2(sort(unique(gene.module)), 
                    col=c( "#D16103","#FFDB6D","#4E84C4","#52854C"))
pick.row.cols=col_func(1:max(gene.module)) #module color for 1,2,...
row.cols=pick.row.cols[gene.module]

ht_opt$TITLE_PADDING = unit(c(4, 4), "points")
row_ha = rowAnnotation(module=gene.module,col=list(module=col_func),show_legend =FALSE);

Heatmap(mat,col = colors,
        row_order=rownames(mat),column_order=colnames(mat),
        #top_annotation = HeatmapAnnotation(module=gene.module,col=list(module=col_func),show_legend =FALSE),
        left_annotation = row_ha,
        #row_title_gp = gpar(fill = unique(row.cols), font = 1:max(wc$membership)),
        column_names_gp = gpar(col=rep(NA,ncol(mat))),
        #row_names_gp = gpar(col=rep(NA,ncol(mat1))),
        #column_names_gp = gpar(col = row.cols,fontsize = 5),
        row_names_gp = gpar(col = row.cols,fontsize =4),
        heatmap_width = unit(12, "cm"), heatmap_height = unit(16, "cm"),
        show_heatmap_legend = FALSE,
        row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), 
        #rect_gp = gpar(col = "grey", lwd = 2),
        #cluster_columns=FALSE, cluster_column_slices = FALSE,
        #border=T)->ht
        border=T,row_split =paste('module', gene.module),
        column_split =paste('module', gene.module),name = "foo")->ht0
#column_split =factor(as.character(gene.module)),name = "foo")->ht

ht<-Heatmap(mat,col = colors,
            row_order=rownames(mat),column_order=colnames(mat),
            column_names_gp = gpar(col=rep(NA,ncol(mat))),
            row_names_gp = gpar(col = row.cols,fontsize =3),
            heatmap_width = unit(20, "cm"), heatmap_height = unit(24, "cm"),
            show_heatmap_legend = FALSE,
            row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), 
            border=T,
            column_split =paste('module', gene.module),column_title='NA',
            row_split =paste('module', gene.module),row_title_rot = 0,
            left_annotation = row_ha,
            name = "foo")
# top bar, from left to right, module 1,2,...
table(gene.module)
#draw(ht, padding = unit(c(4, 2, 10, 0), "cm")) # add space for titles

## GO annotation of the modules
go.outfile=paste0(out.dir,'cut',cut,'_MCL.module_GO.result.rds')
if(!file.exists(go.outfile)){
  GO.result=list();
  for(i in unique(gene.module.info$module)){
   
    genes=gene.module.info[gene.module.info$module==i,]$gene
    out=GOenrich(genes,category='BP',cutoff=0.6) #larger cutoff, smaller returned GO terms
    
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
    GO.result[[i]]=x1;
  }
  saveRDS(GO.result,go.outfile)
}

## plot GO result
plots=list();plots2=list();plots3=list();
library(grDevices) #for using transparent colors
GO.result=readRDS(go.outfile)
for(i in unique(gene.module.info$module)){
  genes=gene.module.info[gene.module.info$module==i,]$gene
  
  x1=GO.result[[i]];
  x1=x1[order(x1$p.adjust),]
  #x1=x1[x1$Count>=5,]
  x1$desp=factor(x1$Description,levels=rev(x1$Description))
  
  count.cut=3; #minimal GO term gene count
  #if(i==4){count.cut=0}else{count.cut=3}
  #plots[[i]]<-ggplot(subset(x1,Count>=3),aes(x=GeneRatio,y=desp,size=Count,col=p.adjust))+
  plots[[i]]<-ggplot(subset(x1,Count>=count.cut),aes(x=GeneRatio,y=desp,size=Count,col=p.adjust))+
    geom_point()+theme_bw(base_size=14)+
    scale_color_gradient(low="blue", high="red")+
    scale_size(range = c(2,8))+
    #scale_y_discrete(labels=y.lab.text)+
    #scale_size(breaks = seq(10,120,30))+
    #ylab('desp, Molecular Function')+
    #ylab(basename(file))+
    #ylab('Biological Process GO enrichment analysis of tissue specific genes')+
    ylab('')+
    ggtitle(paste0('module ',i,', #gene ',length(genes)))+
    theme(
      plot.title =element_text(size=20, face='bold'),
      panel.grid = element_blank(),
      axis.text=element_text(size=14),
      axis.title=element_text(size=20),
      axis.text.x=element_text(size=20,angle=45, hjust=1),
      axis.text.y=element_text(size=24,angle=0, hjust=1),
      axis.ticks.y = element_blank())
  #write.table(x1,'module2.go.txt',quote=F,row.names = F)
  
  x1$log.p.adjust=-1*log(x1$p.adjust,base=10)
  x2=x1;
  x2$desp=as.character(x2$desp)
  x2$x.axis=x2$desp
  if(nrow(x2)>=6){
    x2=x2[1:6,];
    x2$x.axis=factor(x2$x.axis,levels=rev(x2$x.axis))
  }else{
    x2$x.axis=factor(x2$x.axis,levels=rev(x2$x.axis))
  }
  #plots2[[i]]<- ggplot(subset(x2,Count>=count.cut),aes(x=desp,y=log.p.adjust))+
  plots2[[i]]<- ggplot(x2,aes(x=x.axis,y=log.p.adjust))+
    geom_bar(fill=adjustcolor('grey10',alpha.f =0.2),stat='identity',width=0.9)+
    coord_flip()+
    #geom_text(aes(label=desp),vjust=0,hjust = "inward",size=6)+
    #geom_text(aes(label=desp),y=0,vjust=0,hjust = 0,size=6)+
    geom_text(aes(label=desp),y=max(x2$log.p.adjust*0.01),vjust=0,hjust = 0,size=10)+
    theme_bw(base_size=24)+ylab('-log10(p.adjust)')+xlab('')+
    #scale_x_discrete(expand = expansion(mult = c(0, 0))) + 
    ggtitle(paste('module',i))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)),limits = c(0, NA))+
    theme(legend.position = 'none',
          axis.title = element_text(size=22),
          axis.text.y=element_blank(),
          axis.text.x=element_text(size=20),
          plot.title = element_text(size = 20, face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  plots3[[i]]<- ggplot(x2,aes(x=x.axis,y=log.p.adjust))+
    #geom_bar(fill=adjustcolor('grey10',alpha.f =0.2),stat='identity',width=0.9)+
    geom_bar(fill=pick.row.cols[i],stat='identity',width=0.05)+
    geom_point(col=pick.row.cols[i],size=6)+
    coord_flip()+
    geom_text(aes(label=desp),y=max(x2$log.p.adjust*0.01),
              vjust=-0.4,hjust=0,size=10)+
    theme_bw(base_size=20)+ylab('-log10(p.adjust)')+xlab('')+
    ggtitle(paste('module',i))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)),limits = c(0, NA))+
    theme(legend.position = 'none',
          axis.title = element_text(size=22),
          axis.text.y=element_blank(),
          axis.text.x=element_text(size=20),
          plot.title = element_text(size = 20, face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}


pdf(paste0(out.dir,'cut',cut,'_module_GO.pdf'),
    width = 24,height = 12,useDingbats = T)
grid.arrange(grobs = plots3, ncol=2)
dev.off()


pdf(paste0(out.dir,'cut',cut,'_MCL_modules_GO.pdf'),height = 22,width = 18)
draw(ht0, padding = unit(c(-10, 1, 5, 1), "mm")) 
draw(ht, padding = unit(c(-10, 1, 5, 1), "mm")) # add space for titles
#draw(ht) # add space for titles

lay <- rbind(c(1,1,2,4),
             c(1,1,3,5),
             c(6,6,6,7),
             c(6,6,6,7));
grid.arrange(grobs = plots, ncol=2)
grid.arrange(grobs = plots2, ncol=2)

dev.off();



## output all enriched GO terms
GO.result=readRDS(go.outfile)
df.go=Reduce(`rbind`,GO.result)
df.go=as.data.frame(df.go)
df.go$module.id=rep(1:length(GO.result),sapply(GO.result,nrow))
dim(df.go)
write.table(df.go,paste0(out.dir,'cut',cut,'_all.enriched.go.terms.txt'),quote=F,row.names = F,sep='\t')


  