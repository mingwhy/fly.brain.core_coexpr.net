
library(Seurat)
library(tidyverse)
library(Matrix);
library(ggplot2);library(gridExtra)
library(RColorBrewer)
## raw files
if(!file.exists('datasets_overview.txt')){
  files=c('~/Documents/Data_fly_FCA/fly.brain.atlas/raw_stats.rds',
          '~/Documents/Data_fly_FCA/fly.brain_Cocaine/raw_stats.rds',
          '~/Documents/Data_fly_FCA/FCA_head/raw_stats.rds',
          '~/Documents/Data_fly_FCA/xFCA_body/raw_stats.rds')
  out=list();
  for(file in files){
    x=readRDS(file)
    ncell=nrow(x$cell.annotation.df)
    ntype=length(unique(x$cell.annotation.df$annotation))
    ngene=length(x$genes.df)
    print(table(x$cell.annotation.df$sex)) #some 'mixed sex' involved
    cat(file,ncell,ngene,ntype,'\n') 
    out[[file]]=c(ncell,ngene,ntype,sum(x$cell.annotation.df=='female'),
                  sum(x$cell.annotation.df=='male'),sum(x$cell.annotation.df=='mix'));
  }
  #~/Documents/Data_fly_FCA/fly.brain.atlas/raw_stats.rds 56902 17473 116 
  #~/Documents/Data_fly_FCA/fly.brain_Cocaine/raw_stats.rds 43824 10949 39 
  #~/Documents/Data_fly_FCA/FCA_head/raw_stats.rds 100527 13056 82 
  #~/Documents/Data_fly_FCA/xFCA_body/raw_stats.rds 96926 15267 34 
  df.out=as.data.frame(Reduce(`rbind`,out))
  colnames(df.out)=c('#cell','#gene','#cell.cluster','#sex.label=female','#sex.label=male','#sex.label=mix')
  names(out)
  df.out$dataset=c('brain, Davie et al. 2018','brain, Baker et al. 2021',
              'head, Li et al. 2022', 'body, Li et al. 2022');
  colnames(df.out);rownames(df.out)<-NULL
  data.table::fwrite(df.out[,c(7,1:6)],'datasets_overview.txt',sep='\t')
}
df.out=data.table::fread('datasets_overview.txt',sep='\t')
plot.new()
grid.table(df.out)

## 
min.cell=200;
if(!file.exists('filtered_cell.clusters.rds')){
  datasets=c('brain, Davie et al. 2018','brain, Baker et al. 2021',
             'head, Li et al. 2022', 'body, Li et al. 2022')
  out=list();
  for(dat.name in datasets){
    if(dat.name=='body, Li et al. 2022'){
      ## FCA.body
      genes=readRDS('~/Documents/fly.brain_coexpr.net_2022.06/FCA_body/all_common_genes.rds')
      common.genes=genes$common.genes
      file='~/Documents/Data_fly_FCA/xFCA_body/whole_body_filtered_valid.rds';
    }else if(dat.name=='head, Li et al. 2022'){
      ## FCA.head
      genes=readRDS('~/Documents/fly.brain_coexpr.net_2022.06/FCA_head/all_common_genes.rds');
      common.genes=genes$common.genes
      file='~/Documents/Data_fly_FCA/FCA_head/whole_head_filtered_valid.rds';
    }else if(dat.name=='brain, Davie et al. 2018'){
      genes=readRDS('~/Documents/fly.brain_coexpr.net_2022.06/Brain_Davie/all_common_genes.rds');
      common.genes=genes$common.genes
      file='~/Documents/Data_fly_FCA/fly.brain.atlas/wholebrain_filtered_valid.rds';
    }else if(dat.name=='brain, Baker et al. 2021'){
      genes=readRDS('~/Documents/fly.brain_coexpr.net_2022.06/Brain_Baker/all_common_genes.rds');
      common.genes=genes$common.genes
      file='~/Documents/Data_fly_FCA/fly.brain_Cocaine/Brain.control_filtered.rds';
    }
    
    ##
    min.cell=200;
    dat=readRDS(file);
    
    tc_sex=paste(dat$sex,dat$annotation);
    #i=names(which(table(tc_sex)>min.cell))
    dat$tc_sex=tc_sex;
    
    df=dat@meta.data %>% group_by(tc_sex) %>% summarise(n=n())
    cat(dat.name,grep('artefact|unanno',df$tc_sex,ignore.case = T),'\n')
    df=df[order(df$n),]
    
    tmp=dat[,dat$tc_sex %in% df[df$n>min.cell,]$tc_sex]
    cat(dat.name,'ngene',nrow(tmp),'ncell',ncol(tmp),'\n')
    out[[dat.name]]=df
  }
  saveRDS(out,'filtered_cell.clusters.rds')
}
#brain, Davie et al. 2018 ngene 12094 ncell 43659 
#brain, Baker et al. 2021 ngene 10445 ncell 41520 
#head, Li et al. 2022 ngene 12602 ncell 45945 
#body, Li et al. 2022 ngene 14996 ncell 80362 

out=readRDS('filtered_cell.clusters.rds')
names(out)
sapply(out,nrow)

sapply(out,function(x){ length(unique(gsub('^male|^female','',x$tc_sex))) })
#115, 39, 80, 28

sapply(out,function(x){ length(grep('^male',x$tc_sex)) })
sapply(out,function(x){ length(grep('^female',x$tc_sex)) })

sapply(out,function(x){  x=x[x$n>=min.cell,]; 
        length(unique(gsub('^male|^female','',x$tc_sex))) })
#38, 28, 40, 19
sapply(out,function(x){ x=x[x$n>=min.cell,]; length(grep('^male',x$tc_sex)) })
#brain, Davie et al. 2018 brain, Baker et al. 2021     head, Li et al. 2022     body, Li et al. 2022 
#31                       28                       36                       17 
sapply(out,function(x){  x=x[x$n>=min.cell,]; length(grep('^female',x$tc_sex)) })
#brain, Davie et al. 2018 brain, Baker et al. 2021     head, Li et al. 2022     body, Li et al. 2022 
#37                       26                       40                       18 

#for(dat.name in names(out)){
plots<-lapply(names(out),function(dat.name){
  df=out[[dat.name]]
  df$col='grey';
  df[df$n>=200,]$col='yellow'
  df$tc_sex=factor(df$tc_sex,levels=df$tc_sex)
  
    ggplot(df,aes(x=n,y=tc_sex,fill=col))+geom_bar(stat='identity')+
      scale_fill_manual(values=c("#999999","#E69F00"))+theme_bw()+
      ylab('')+xlab('')+
      #scale_x_log10()+
      ggtitle(dat.name)+theme(axis.text.x = element_text(size=12),
                              axis.text.y = element_blank(),
                              legend.position = 'none')
  
})

pdf('filtered_cell.clusters.pdf',useDingbats = T,width = 12,height = 8)
plot.new()
grid.table(df.out)
do.call('grid.arrange',c(plots,ncol=4))
dev.off()



