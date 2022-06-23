

options(stringsAsFactors = F)
library(igraph)
library(readxl)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(clusterProfiler)
library(AnnotationDbi)
library(GO.db);library(gprofiler2)
library(RColorBrewer)
library(org.Dm.eg.db,verbose=F,quietly=T)
source('./src_fly_go_enrichment.R')

topK=0.05;
edge.cut=43;
(files=Sys.glob(paste0('../03.top.perc_net/brain_scRNA-seq_n15c0.005_bigScale2/top_*_pan.rds')))
(file=files[grep(topK,files)])

out.dir='./brain_scRNA-seq_n15c0.005/';
if(!dir.exists(out.dir)){dir.create(out.dir)}

col_rnorm = circlize::colorRamp2(c(-3, 0, 3), c("green", "white", "red"))
col_runif = circlize::colorRamp2(c(0, 3), c("white", "orange"))
col_letters = c("a" = "pink", "b" = "purple", "c" = "blue")

####################################################################
## read in fly gene age data from 'No Evidence for Phylostratigraphic Bias Impacting Inferences on Patterns of Gene Emergence and Evolution'
# more information can be found here: https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/Phylostratigraphy_fly.gene_age 
div.times=c(25,235,325,521,542,560,573,615,653,1547,2100,3900);
length(div.times)
df2017<-read_excel('../external_data/phylostratigraphy_fly.gene.age/TableS1.xlsx',sheet='Dmel_updated_phylostratigraphy',skip=3)
dim(df2017) #13794     4
table(df2017$phylostrata)

# this data contain protein ids
# check the overlap genes between fly  prot_id in 2017.dataset and fly.prot from org.Dm.eg.db.
columns(org.Dm.eg.db)
fly.prot<-keys(org.Dm.eg.db, keytype='FLYBASEPROT')
length(fly.prot) #30563 fly proteins in total of fly genome in this database
sum(df2017$prot_id %in% fly.prot) #13624 overlap genes
# get flybase ID for each fly prot_id
df.fly.prot<-AnnotationDbi::select(org.Dm.eg.db,
                                   keys=df2017$prot_id,
                                   keytype="FLYBASEPROT",c("FLYBASECG","FLYBASE",'SYMBOL',"GENENAME"))
dim(df.fly.prot) #13794     5
sum(is.na(df.fly.prot$FLYBASE)) #170, equals 13794-13624
sum(is.na(df.fly.prot$SYMBOL)) #170, equals 13794-13624
gene.age=merge(df.fly.prot,df2017,by.x='FLYBASEPROT',by.y='prot_id')
dim(gene.age) #13794
table(gene.age$phylostrata)
sum(is.na(gene.age$FLYBASE)) #170 genes without FLYBASE id.
sum(table(gene.age$phylostrata)) #13794 genes in df2017 publication
gene.age2=gene.age[!is.na(gene.age$FLYBASE),]
dim(gene.age2) #13624 genes
head(gene.age2)
length(unique(gene.age2$SYMBOL)) #13624 fly symbols

gene.age=gene.age2

#############################################################################
## choose edge commonality cutoff to characterize core

      plots=list();iplot=0;
      
      pan=readRDS(file);
      dim(pan)
      net= pan >= edge.cut; #True or False binary matrix
          
      net0=net;
      diag(net0)=FALSE;
      i=apply(net0,1,sum)
      sum(i==0);length(i)
      net1=net[i!=0,i!=0]
      dim(net1)
      cat('cutoff',topK,'nlink.cut',edge.cut,',ngene',nrow(net1),'\n') #cutoff 0.01, nlink.cut 14, ngene 179
      #cutoff 0.05 nlink.cut 43 ,ngene 205 
      
      pan.genes=rownames(pan)
      core.genes=rownames(net1)
      length(pan.genes); #2088 genes
      length(core.genes)    # 205 genes   
      
      sum(pan.genes %in% gene.age$SYMBOL) # 2222 common.expr
      df=data.frame(genes=pan.genes,status='not.common')
      df[df$genes %in% core.genes,]$status='common';
      df1=merge(gene.age,df,by.x='SYMBOL',by.y='genes')
      dim(df1) #  2002 common.expr
      head(df1)
      table(df1$phylostrata)
      #1   2   3   4   5   6   7   8   9  10  11  12 
      #911 641  87 122  32  73  14  16  13  43  39  11 
      
      iplot=iplot+1; #2222 common.expr VS 179 core.genes
      plots[[iplot]]<-ggplot(df1,aes(x=factor(phylostrata),fill=status))+
        geom_bar(position = 'dodge')+theme_bw()+xlab('Gene age (old<->young)')+
        #scale_x_discrete(expand = c(0, 0))+
        scale_y_continuous(expand = c(0, 100), limits = c(0, NA))+
        scale_fill_manual(labels=c('core','non.core'),values=c("#D95F02","#56B4E9"))
      
      x1=df1 %>% group_by(phylostrata,phylostrata_name) %>% summarise(ngene.in.age.class=n())
      x1$total.gene=sum(x1$ngene.in.age.class)
      x1$prop=x1$ngene.in.age.class/x1$total.gene;
      x1$status='background';
      
      x2=subset(df1,status=='common') %>% group_by(phylostrata,phylostrata_name) %>% summarise(ngene.in.age.class=n())
      x2$total.gene=sum(x2$ngene.in.age.class)
      x2$prop=x2$ngene.in.age.class/x2$total.gene;
      x2$status='common';
      dim(x1);dim(x2)
      x=rbind(x1,x2)
      
      x[x$status=='background',]$total.gene
      sum(x[x$status=='background',]$ngene.in.age.class)
      sum(x[x$status=='background',]$ngene.in.age.class)
      iplot=iplot+1; #2222 common.expr VS 179 core.genes
      plots[[iplot]]<-ggplot(x,aes(x=factor(phylostrata),y=prop,fill=status))+
        geom_bar(stat='identity',position = 'dodge')+theme_bw()+xlab('Gene age (old<->young)')+
        #scale_x_discrete(expand = c(0, 0))+
        scale_y_continuous(expand = c(0, 0.01), limits = c(0, NA))+
        scale_fill_manual(labels=c('background','core'),values=c("#999999","#D95F02"))
      
      ## add p value as marks above the bars
      df.com12=x;
      pvalues=lapply(unique(df.com12$phylostrata),function(group){
        x=df.com12[df.com12$phylostrata==group,]
        if(nrow(x)==1){return(NA)}
        tmp1=x[x$status=='common',]
        tmp2=x[x$status!='common',]
        x1=tmp1$ngene.in.age.class; #TF gene and common
        x2=tmp2$ngene.in.age.class-x1 #TF gene and not-common
        y1=tmp1$total.gene-x1 #common and not-TF
        y2=tmp2$total.gene-tmp2$ngene.in.age.class-y1
        #sum(x1,x2,y1,y2)==tmp2$total.gene
        mat=matrix(c(x1,y1,x2,y2),ncol=2)
        tmp=fisher.test(mat,alternative='greater') #one-sided
        #tmp=fisher.test(mat) #two-sided
        tmp$p.value})
      names(pvalues)=unique(df.com12$phylostrata)
      pvalues.adjust=p.adjust(pvalues,method='BH')
      pvalues.adjust
      
      dfp=data.frame(phylostrata =names(pvalues.adjust),adjust.pvals=pvalues.adjust)
      df.com12.p=merge(dfp,df.com12,by='phylostrata')  
      #df.com12.p$lab = ifelse(df.com12.p$adjust.pvals <= .001, '**', ifelse(df.com12.p$adjust.pvals <= .05, '*', ''))
      df.com12.p$lab = ifelse(df.com12.p$adjust.pvals <= .001, '***', ifelse(df.com12.p$adjust.pvals <= .01, '**', ifelse(df.com12.p$adjust.pvals<0.05,"*",'')))
      
      df.com12.p[df.com12.p$status!='common',]$lab=''
      head(df.com12.p)
      table(df.com12.p$lab)
      df.com12.p$phylostrata=factor(df.com12.p$phylostrata,
                                    levels=sort(unique(as.numeric(df.com12.p$phylostrata))))
      # posistion='dodge'
      #pos = position_dodge(.9)
      pos = position_dodge(0)
      iplot=iplot+1;
      y=df.com12.p$prop*1.2
      
      df.com12.p[df.com12.p$phylostrata_name=='cellular_organisms',]$phylostrata_name='CellLife'
      tmp=df.com12.p[,c('phylostrata','phylostrata_name')]
      my.x.lab=unique(tmp[order(tmp$phylostrata),]$phylostrata_name)
      
      plots[[iplot]]<-ggplot(df.com12.p,aes(x=factor(phylostrata),y=prop,fill=status))+
        geom_bar(stat='identity',position = 'dodge')+theme_bw(base_size = 16)+xlab('Gene age (old<->young)')+
        scale_x_discrete(labels=my.x.lab)+
        #scale_x_discrete(expand = c(0, 1.2))+
        scale_y_continuous(limits = c(0, NA),expand= expansion(mult = c(0, 0.4)))+
        scale_fill_manual(name='',labels=c('commonly.expressed.genes','core genes'),values=c("#999999","#D95F02"))+
        #scale_fill_manual(labels=c('background','core'),values=c("#999999","#D95F02"))+
        geom_text(aes(y=y,label=df.com12.p$lab), 
                  hjust='center', vjust=-0.05, size=10, angle=0, position=pos)+
        #ylab('Percentage of genes')+
        ylab('Frequency of genes')+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              axis.text.x=element_text(size=18,angle=45,hjust=1),
              axis.text.y=element_text(size=12),
              axis.title=element_text(size=18,face="bold"),
              legend.key.size = unit(1, 'cm'),
              legend.key.height = unit(1, 'cm'),
              legend.key.width = unit(1, 'cm'), 
              legend.position = 'top',
              legend.title=element_text(size=16,face="bold"),
              legend.text=element_text(size=14),
              plot.margin=unit(c(1.5, 2, 2, 5.5),"cm"))
      
      
      pdf(paste0(out.dir,'/top_',topK,'_edge.cut',edge.cut,'_gene.age_coreVSnot.core.pdf'),
          width=12,height = 6,useDingbats = T)
      for(i in plots){print(i)}
      dev.off()
      
      
