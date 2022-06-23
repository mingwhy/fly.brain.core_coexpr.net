

library(tidyverse)
library(ggplot2)
library(gridExtra)

(mycol=c("#C4961A",'#E7298A',"#52854C", "#4E84C4")) #"#D16103" red
barplot(1:4,col=mycol)
datasets=c('brain, Davie et al. 2018','brain, Baker et al. 2021',
           'head, Li et al. 2022','body, Li et al. 2022')
all.out=list();

####
for(dat.name in datasets){
  if(dat.name=='body, Li et al. 2022'){
    ## FCA.body
    out.dir='~/Documents/fly.brain_coexpr.net_2022.06/FCA_body/05_rank_aggregation/RAA_output/';
    (files=Sys.glob(paste0(out.dir,'/*sample*edges.rds')));
    length(files)
    ncell.cluster=31
  }else if(dat.name=='head, Li et al. 2022'){
    ## FCA.head
    out.dir='~/Documents/fly.brain_coexpr.net_2022.06/FCA_head/05.rank_aggregation/RAA_output/';
    (files=Sys.glob(paste0(out.dir,'/*sample*edges.rds')));
    length(files)
    ncell.cluster=39
  }else if(dat.name=='brain, Davie et al. 2018'){
    out.dir='~/Documents/fly.brain_coexpr.net_2022.06/Brain_Davie/05.rank_aggregation/RAA_output/';
    (files=Sys.glob(paste0(out.dir,'/*sample*edges.rds')));
    length(files)
    ncell.cluster=67
  }else if(dat.name=='brain, Baker et al. 2021'){
    out.dir='~/Documents/fly.brain_coexpr.net_2022.06/Brain_Baker/05_rank_aggregation/RAA_output/';
    (files=Sys.glob(paste0(out.dir,'/*sample*edges.rds')));
    length(files)
    ncell.cluster=47
  }
  
  
  ###################################################################
  ## cutoff on x-axis, -log10(P) value on y axis
  plots=list();
  out=lapply(files, function(file){
    #x=readRDS("./result_RAA_sub.mat/RAA.per.edge_at.edge.10_sample100edges.rds")
    x=readRDS(file)
    if(length(x)==0){next} #edge.commonality at 29
    (commonality=as.numeric(strsplit(file,'at.edge.|_sample')[[1]][[2]]))
    max(x$Score)
    x$derive.p.value=x$Score
    #x$derive.p.value=x$Score*(ncell.cluster-commonality);
    x$adjust.p= p.adjust(x$derive.p.value,method='BH')  
    max(x$adjust.p)
    x$logP=-1*log10(x$adjust.p)
    x$cutoff=commonality
    x
  })
  df=as.data.frame(Reduce(`rbind`,out))
  df$cutoff=factor(df$cutoff)  
  x=df %>% group_by(cutoff) %>% summarise(p=1-sum(adjust.p<0.05)/length(adjust.p))
  dfx=merge(df,x,by='cutoff')
  
  dfx$dataset=dat.name
  dfx;
  all.out[[dat.name]]=dfx;
}
names(all.out)
saveRDS(all.out,'plot_rank.aggre_all.four.rds')

####
all.out=readRDS('plot_rank.aggre_all.four.rds')
dfx=as.data.frame(Reduce(`rbind`,all.out))
head(dfx)
table(dfx$dataset)

dfx$cutoff=factor(dfx$cutoff,unique(sort(as.numeric(as.character(dfx$cutoff)))))
#dfx=dfx[as.numeric(dfx$cutoff) %% 2!=0,]
dfx$p=round(dfx$p,2)
Time=factor(dfx$p, levels=unique(dfx$p))
dfx$dataset=factor(dfx$dataset,levels=datasets)

#################################################
## fisher method combine p value
names(all.out)
head(all.out[[1]])

fdr.out=list()
for(dat.name in names(all.out)){
  x=all.out[[dat.name]]
  p=sapply(sort(unique(x$cutoff)),function(i){
    pvals= x[x$cutoff==i,]$adjust.p;
    #if(length(pvals)==1){return(NA)}
    if(length(pvals)<5){return(NA)}
    df = 2*length(pvals)
    pchisq( -2*sum(log(pvals)), df, lower.tail=FALSE)
  })
  #fdr=p.adjust(p,method='BH')
  fdr=p.adjust(p,method='bonferroni')
  df=data.frame(edge.commonality.group=sort(unique(x$cutoff)),fdr=fdr,dataset=dat.name)
  fdr.out[[dat.name]]=df
}
df=as.data.frame(Reduce(`rbind`,fdr.out))
head(df)
df$dataset=factor(df$dataset,levels=datasets)
df$edge.commonality.group=factor(df$edge.commonality.group,
                                 levels=unique(sort(as.numeric(as.character(df$edge.commonality.group)))))
df$log_fdr= -1*log10(df$fdr)

## combine two plots
datasets
scale.factors=c(20,10,2.5,2.5);
plots<-lapply(1:length(datasets),function(i){
  tmp1=dfx[dfx$dataset==datasets[[i]],]
  tmp2=df[df$dataset==datasets[[i]],]
  #tmp1=dfx[dfx$dataset=='head, Li et al. 2022',]
  #tmp2=df[df$dataset=='head, Li et al. 2022',]
  scaleFactor=scale.factors[i]
  p=ggplot()+
    geom_violin(aes(x=tmp1$cutoff,y=tmp1$logP),col=mycol[i],scale = "width")+
    geom_jitter(aes(x=tmp1$cutoff,y=tmp1$logP),col=mycol[i],width = 0.2,size=0.4)+ 
    geom_line(aes(x=tmp2$edge.commonality.group,
              y=tmp2$log_fdr/scaleFactor,group=tmp2$dataset),col='black',size=0.6)+
    #scale_color_manual(values=mycol[i])+
    stat_summary(aes(x=tmp1$cutoff,y=tmp1$logP),fun.y=median, geom="point", size=1, color="black")+
    xlab('Edge commonality')+
    scale_y_continuous(
      name = "-log10(Adjusted P value)", # Features of the first axis
      sec.axis = sec_axis(trans=~.*scaleFactor,
            name="-log10(Bonferroni corrected P)"))+  # Add a second axis and specify its features
    theme_classic(base_size = 12)+
    theme(panel.grid = element_blank(),
          axis.text = element_text(size=8),
          axis.title = element_text(size=10),
          legend.position = 'none')
  if(i==1){
    p=p+scale_x_discrete(breaks=seq(1,63,2))
  }
  p
} )

pdf('plot_rank.aggre_all.four_fdr.pdf',useDingbats = T,width = 6,height = 12)
#print( grid.arrange(grobs=plots,ncol=2) )
print( grid.arrange(grobs=plots,ncol=1) )
dev.off()

pdf('RAA_for_brain_davie.pdf',width=7,height = 3,useDingbats = T)
print( plots[[1]]+ theme_bw(base_size = 8)+
         xlab('Edge commonality')+
  theme(axis.text = element_text(size=10),
        axis.text.x = element_text(size=10,angle=0,hjust=0.5),
        axis.title = element_text(size=12,face="bold"),
        legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
)
dev.off()

cutoffs=lapply(datasets,function(i){
  tmp=df[df$dataset==i,];
  as.numeric(tmp[which(tmp$fdr<0.01),]$edge.commonality.group[1])
})
names(cutoffs)=datasets
cutoffs

data.table::fwrite(as.data.frame(cutoffs),'core_cutoffs.txt')


