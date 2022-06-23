

library(ggplot2)
library(gridExtra)
out.dir='RAA_output/';
(files=Sys.glob(paste0(out.dir,'/*sample*edges.rds')));
length(files)


plots=list();
for(file in files){
  #x=readRDS("./result_RAA_sub.mat/RAA.per.edge_at.edge.10_sample100edges.rds")
  x=readRDS(file)
  if(length(x)==0){next} #edge.commonality at 29
  (commonality=as.numeric(strsplit(file,'at.edge.|_sample')[[1]][[2]]))
  
  x$derive.p.value=x$Score
  x$adjust.p= p.adjust(x$derive.p.value,method='BH')
  max(x$adjust.p)
  breaks=c(0,1e-06,1e-05,1e-04,1e-03,0.01,0.05,1)
  
  groups=cut(x$adjust.p,breaks = breaks,include.lowest = T)
  x$groups=groups
  tmp=as.data.frame(table(x$groups))
  sum(tmp$Freq);nrow(x)
  
  plots[[as.character(commonality)]]<-ggplot(tmp,aes(x=Var1,y=Freq))+geom_bar(stat='identity',fill=NA,col='black')+
    geom_text(aes(label=Freq),vjust=-0.4)+theme_bw(base_size = 20)+
    xlab('Adjusted P value')+ylab('Frequency')+
    scale_y_continuous(expand = expansion(mult = c(0, 0.2)), limits = c(0, NA))+
    #ggtitle(paste0('aggregated scores of sampled ',nrow(x),
    #        ' edges at commonality ',commonality,' in (33-',commonality,')=',33-commonality,' cell types'))+
    #ggtitle(paste0('edges at commonality ',commonality,' in (33-',commonality,')=',33-commonality,' cell types'))+
    ggtitle(paste0('edges at commonality ',commonality))+
    theme(panel.grid = element_blank(),
          axis.text = element_text(size=20),
          axis.text.x = element_text(size=20,angle=45,hjust=1))
  
}

plots1=plots[order(as.numeric(names(plots)))]
length(plots1)

pdf(paste0(out.dir,'rank_aggre_score_all.pdf'),width=20,height=16,useDingbats =T)
do.call('grid.arrange',c(plots1,ncol=4))
dev.off()

###################################################################
## cutoff on x-axis, -log10(P) value on y axis
plots=list();
out=lapply(files, function(file){
  #x=readRDS("./result_RAA_sub.mat/RAA.per.edge_at.edge.10_sample100edges.rds")
  x=readRDS(file)
  if(length(x)==0){next} #edge.commonality at 29
  (commonality=as.numeric(strsplit(file,'at.edge.|_sample')[[1]][[2]]))  
  x$derive.p.value=x$Score
  x$adjust.p= p.adjust(x$derive.p.value,method='BH') 
  x$logP=-1*log10(x$adjust.p)
  x$cutoff=commonality
  x
})
df=as.data.frame(Reduce(`rbind`,out))

df$cutoff=factor(df$cutoff)  

library(tidyverse)
x=df %>% group_by(cutoff) %>% 
  summarise(p=1-sum(adjust.p<0.05)/length(adjust.p))

dfx=merge(df,x,by='cutoff')

#pdf(paste0(out.dir,'rank_aggreg.pdf'),width=24,height=12,useDingbats =T)
pdf(paste0(out.dir,'rank_aggreg.pdf'),width=12,height=8,useDingbats =T)

dfx$p=round(dfx$p,2)
Time=factor(dfx$p, levels=unique(dfx$p))
print(
ggplot(dfx,aes(x=cutoff,y=logP))+
  #geom_boxplot(outlier.shape = NA)+
  #geom_text(aes(x=cutoff,y=10,label=dfx$p))+
  geom_violin(scale = "width")+
  #geom_text(aes(y=10,label=paste(Time)),
  geom_text(aes(y=7.5,label=paste(Time)),
            position = position_dodge(width = 1), 
            vjust = -0.5, size = 5, stat = "unique", parse = TRUE)+
  geom_jitter(width = 0.2,size=1.5)+ #sc
  #geom_jitter(width = 0.2,size=0.8)+
  #stat_summary(fun.y=median, geom="point", size=2, color="red")+
  ylab('-log10(Adjusted P value)')+xlab('Edge commonality group')+
  #geom_hline(yintercept=-1*log10(0.01),linetype="dashed", color = "red")+
  geom_hline(yintercept=-1*log10(0.05),linetype="dashed", color = "red")+
  theme_classic(base_size = 24)+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18))
)

dev.off();
