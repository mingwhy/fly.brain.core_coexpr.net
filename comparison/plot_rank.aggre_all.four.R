

library(tidyverse)
library(ggplot2)
library(gridExtra)

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

library(RColorBrewer)
(mycol=brewer.pal(4,"Paired"))
# "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C"
(mycol=c("#C4961A","#D16103","#52854C", "#4E84C4"))

#datasets[c(3,4,1,2)]
#[1] "brain, Davie et al. 2018" "brain, Baker et al. 2021" "body, Li et al. 2022"    
#[4] "head, Li et al. 2022"  

dfx$dataset=factor(dfx$dataset,levels=datasets)
p1=ggplot(dfx,aes(x=cutoff,y=logP,col=dataset))+
  #ggplot(tmp,aes(x=cutoff,y=logP,col=dataset))+
  facet_wrap(.~dataset,scale='free')+
  geom_violin(scale = "width")+
  #geom_text(aes(y=10,label=paste(Time)),
  geom_jitter(width = 0.2,size=0.4)+ #sc
  scale_color_manual(values=mycol)+
  stat_summary(fun.y=median, geom="point", size=1, color="black")+
  ylab('-log10(Adjusted P value)')+xlab('Edge commonality group')+
  #geom_hline(yintercept=-1*log10(0.01),linetype="dashed", color = "red")+
  #geom_hline(yintercept=-1*log10(0.05),linetype="dashed", color = "red")+
  theme_classic(base_size = 12)+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=8),
        axis.title = element_text(size=12),
        legend.position = 'none')

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
    if(length(pvals)<10){return(NA)}
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

p2<-ggplot(df,aes(x=edge.commonality.group,y=log_fdr,col=dataset,group=dataset))+geom_point()+
  geom_line()+theme_bw()+scale_color_manual(values=mycol)+
  #ylab('-log10(FDR)')+
  xlab('Edge commonality group')+
  ylab('-log10(Bonferroni correction P value)')
 
min(df$fdr,na.rm=T)
pdf('plot_rank.aggre_all.four.pdf',useDingbats = T,width = 12)
print(p1);print(p2);
#print(p2+ scale_y_log10()+geom_abline(slope=0,intercept = 1,linetype="dashed", color = "black"))

dev.off()

View(df[grep('head',df$dataset),])
View(df[grep('body',df$dataset),])

