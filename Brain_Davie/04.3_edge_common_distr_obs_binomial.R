
options(stringsAsFactors = F)
library(igraph)
library(ggplot2)
library(gridExtra)
library(tidyverse)

# file from code: edge_common_distr_obs.R
(files=Sys.glob('./brain_scRNA-seq_n15c0.005/top_*_commonality.distribution.rds'));
(files=files[grep('0.05',files)])

pdf('./brain_scRNA-seq_n15c0.005/observed_binomial.pdf',width = 9,height = 6)

n.gene=2088;
n.cell.type=67;
for(file in files){
  (cutoff=strsplit(basename(file),'\\_')[[1]][[2]])
  obs=readRDS(file)
  obs=obs$edge.common.distr
  
  prob=as.numeric(cutoff)
  theo=dbinom(0:n.cell.type,size=n.cell.type,prob=prob)
  #plot(0:n.cell.type,theo)
  sprintf('%.04f',theo)
  choose(n.gene,2)*theo[1]
  choose(n.gene,2)*(1-theo[1])
  
  sum(obs) #number of unique edges
  sum(as.numeric(names(obs))*obs)/n.cell.type
  choose(n.gene,2)*prob #should be equal to above
  
  obs.freq=theo
  names(obs.freq)=0:n.cell.type
  obs.freq[names(obs)]=obs/sum(obs)

  #df=data.frame(score=0:n.cell.type,obs=obs.freq,theo=theo)
  df=data.frame(score=rep(0:n.cell.type,2),value=c(obs.freq,theo),
                type=c(rep('observed',1+n.cell.type),rep('binomial',1+n.cell.type)));
  print(
    ggplot(df,aes(x=score,y=value,col=type))+geom_line()+geom_point()+
    theme_bw(base_size = 16)+
    xlab('Edge commonality')+ylab('Frequency of edges')+
    scale_y_continuous(trans = scales::pseudo_log_trans(1e-05, 10),
                       breaks=c(0, 1e-05,1e-04,1e-03,0.01, 0.1, 1))+
    theme(axis.text = element_text(size=16),
          axis.text.x = element_text(size=16,angle = 0),
          axis.title = element_text(size=16,face="bold"),
          #legend.position = 'none',
          legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    )+ggtitle(paste('top',cutoff))
  )
}
dev.off()


