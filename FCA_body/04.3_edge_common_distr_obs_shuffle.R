
# compare observed edge.commonality.distribution 
# with the randomized ones

library(Matrix);
library(ggplot2);library(gridExtra)
library(RColorBrewer)
library(igraph);library(tidyverse)

path='brain_snRNA-seq_n15c0.005_bigScale2/';
path.shuffle='net_rewire_rep20/';
out.dir='brain_snRNA-seq_n15c0.005/';

## collect cutoff values
topK='top_0.05';
(files=Sys.glob(paste0(path,'/top_*_networks_cutoff_cell.type.rds')))
(files=files[grep(topK,files)])


for(file in files){
  (cutoff=strsplit(basename(file),'top\\_|\\_networks',)[[1]][2])
  
  simu.files=Sys.glob(paste0(path.shuffle,'/top_',cutoff,'_shuffle_rep*.rds'))
  simu.files
  
  (output.file=paste0(out.dir,'/top_',cutoff,'_edge.common_distr_obs_shuffle.rds'));
  
  if(file.exists(output.file)){
    #simu.result=readRDS(output.file);
  }else{
    nrep=length(simu.files)
    simu.nets=readRDS(simu.files[1])
    my.common.genes=unique(as.character(unlist(sapply(simu.nets,function(x)names(igraph::V(x))))))
    length(my.common.genes) #2368
    
    adj0=matrix(0,ncol=length(my.common.genes),nrow=length(my.common.genes));
    rownames(adj0)=my.common.genes;
    colnames(adj0)=my.common.genes;
    dim(adj0)
    
    simu.result=list();
    irep=0;
    for(simu.file in simu.files){
      irep=irep+1;
      simu=readRDS(simu.file)
      length(simu) #number of cell types
      
      simu.adj=lapply(simu,function(x){
        net=as.matrix(igraph::as_adj(x));
        #genes=rownames(net);
        #i=genes %in% my.common.genes
        #net=net[i,i]
        adj=adj0;
        adj[rownames(net),colnames(net)]<-net;
        adj})
      
      sapply(simu.adj,dim) #all the same dimension
      x<-Reduce(`+`,simu.adj); #combine matrix to get link commonality
      
      edge.common=x[upper.tri(x)]
      edge.common.distr=table(edge.common)
      simu.result[[irep]]<-edge.common.distr
      rm(list=c('simu','simu.adj'));
      cat('simu for cutoff',cutoff,'at irep',irep,'is done\n')
    }
    
    saveRDS(simu.result,paste0(output.file))
    cat('simu for cutoff',cutoff,'is done\n')
  }
}
#############################################################3
## calcualte obs VS null and plot
obs.list=list();
(files=Sys.glob(paste0(out.dir,'/top_*_commonality.distribution.rds'))); #generated from `01_edge_common_distr_obs.R`
(files=files[grep(topK,files)])

for(file in files){
  (cutoff=strsplit(basename(file),'\\_')[[1]][[2]])
  result=readRDS(file);
  names(result)
  obs=result[['edge.common.distr']];
  
  names(obs)
  
  ### compare the two
  obs# remove commonality==0 and calcualte link.commonality distribution
  x1=obs[names(obs)!=0]; #remove 0
  obs.prop=x1/sum(x1)
  obs.prop.df=data.frame(cutoff=rep(cutoff,length(obs.prop)),
                         commonality= names(obs.prop),
                         count=as.numeric(x1),
                         prop=as.numeric(obs.prop))
  obs.prop.df$commonality=factor(obs.prop.df$commonality,levels=
                 sort(as.numeric(as.character(obs.prop.df$commonality))));
  head(obs.prop.df)
  obs.list[[cutoff]]=obs.prop.df
}


# begin plot
plots=list();plots.line=list();
plots.bar=list();
dist.list=list();
(files=Sys.glob(paste0(path,'/top_*_networks_cutoff_cell.type.rds')))
(files=files[grep(topK,files)])

for(file in files){
  (cutoff=strsplit(basename(file),'\\_')[[1]][[2]])
  file2=paste0(out.dir,'/top_',cutoff,'_edge.common_distr_obs_shuffle.rds');
  
  simu.result=readRDS(file2);
  combine=data.frame()
  for(i in 1:length(simu.result)){
    x=as.data.frame(simu.result[[i]])
    x$irep=i;
    x$status='simu';
    combine=rbind(combine,x)
  }
  
  ## observsed and theoretical, frequency/prop values
  obs=obs.list[[cutoff]]
  colnames(obs);
  obs.df=data.frame(edge.common=obs$commonality,
                    Freq=obs$prop,
                    irep='obs',status='obs.prop')
  colnames(combine)
  colnames(obs.df);
  combine=rbind(combine,obs.df)
  
  unique(combine$edge.common)
  combine=combine[combine$edge.common!=0,]
  simu=tapply(combine[combine$status=='simu',]$Freq,
              combine[combine$status=='simu',]$edge.common,mean)
  simu[is.na(simu[1:length(simu)])]=0;
  simu=simu[names(simu)!='0']
  
  unique(combine$status)
  #obs.here=combine[combine$status=='obs',]$Freq
  #names(obs.here)=combine[combine$status=='obs',]$edge.common
  #https://www.r-bloggers.com/2017/01/practical-kullback-leibler-kl-divergence-discrete-case-2/
  #obs2=obs.here/sum(obs.here);
  obs2=combine[combine$status=='obs.prop',]$Freq
  names(obs2)=combine[combine$status=='obs.prop',]$edge.common
  simu2=simu/sum(simu);
  simu2=simu2[names(obs2)]
  #(dist=sprintf('%.4f',entropy::KL.empirical(simu2,obs2))); #KL divergence
  (dist=philentropy::JSD(rbind(obs2,simu2),unit='log2')) #JSD distance
  #(dist=sum(as.numeric(names(simu2))^2*(obs2-simu2)^2)) #weighted Euclidean distance,https://math.stackexchange.com/questions/917066/calculating-weighted-euclidean-distance-with-given-weights 
  #out=ks.test(obs2,simu2);dist=as.numeric(out$statistic); #Kolmogorov-Smirnov test
  (dist=sprintf('%.5f',dist))
  
  dist.list[[cutoff]]=dist;
  
  combine2<-combine %>% group_by(irep) %>% mutate(prop=Freq/sum(Freq))
  #table(combine2[combine2$status=='simu',]$edge.common)
  #max(combine2[combine2$status=='simu',]$edge.common) #36
  
  line.group=combine2$status;
  combine2$status=factor(combine2$status)
  combine2$edge.common=as.numeric(as.character(combine2$edge.common))
  
  pick.edge=obs;
  nedge=sum(pick.edge$count)
  plots.bar[[cutoff]]<-ggplot(pick.edge,aes(x=factor(commonality),y=prop,fill=factor(cutoff)))+
    xlab('Edge commonality')+ylab('Frequency of edges')+
    scale_y_continuous(trans = scales::pseudo_log_trans(1e-05, 10),
                       breaks=c(0, 1e-05,1e-04,1e-03,0.01, 0.1, 1))+
    #geom_text(label=pick.edge$count,hjust=0.4,vjust=-0.2)+
    #geom_text(label=paste0(tmp$count,'\n',round(tmp$prop,4)),
    #          hjust=0.4,vjust=-0.2)+
    scale_fill_manual(name='',values=c("#E69F00"))+
    geom_bar(stat='identity',position='dodge')+
    #scale_fill_manual(values=mycol)+
    theme_bw(base_size = 14)+
    ggtitle(paste('top',cutoff,',',nedge,' edges in total'))+
    #ggtitle(paste0(nedge,' edges in the union of all cell type-specific networks\nPCC top',cutoff))+
    theme(axis.text = element_text(size=14),
          axis.text.x = element_text(size=14,angle = 0),
          axis.title = element_text(size=14,face="bold"),
          legend.position = 'none',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    )
  
  plots.line[[cutoff]]<- ggplot(combine2 %>% group_by(status,edge.common) %>% summarize(mean.prop=mean(prop)),
                                aes(x=edge.common,y=mean.prop,group=status,col=status))+
    geom_jitter(data=combine2,aes(x=edge.common,y=prop,
                                  #shape=factor(status),
                                  size=factor(status),col=factor(status)))+
    theme_bw()+  geom_line(size=0.3)+
    #scale_y_log10()+
    scale_x_continuous(breaks=c(1,5*(1:24)),labels=c(1,5*(1:24)))+
    #scale_y_continuous(trans = scales::pseudo_log_trans(1e-07, 10),
    #                   breaks=c(0, 1e-07,1e-06,1e-05,1e-04,1e-03,0.01, 0.1, 1))+
    scale_y_continuous(trans = scales::pseudo_log_trans(1e-05, 10),
                       breaks=c(0, 1e-05,1e-04,1e-03,0.01, 0.1, 1))+
    #scale_x_log10()+
    #ggplot(combine,aes(x=edge.common,y=Freq,col=factor(status)))+
    #ggplot(subset(combine,edge.common!=0),aes(x=edge.common,y=Freq,col=factor(status)))+
    
    xlab('Edge commonality')+ylab('Frequency of edges')+
    #scale_shape_manual(values=c(16, 17))+
    #scale_x_discrete(labels=xlabels2)+
    scale_color_manual(name='',values=c('#E69F00','#999999',"#56B4E9"))+
    scale_size_manual(values=c(2,1,2),guide = FALSE)+
    ggtitle(paste('top',cutoff,dist))+
    theme(axis.text=element_text(size=14),
          axis.text.x=element_text(size=12,angle=0,hjust=0.5),
          axis.title=element_text(size=14,face="bold"),
          legend.title=element_text(size=14,face="bold"),
          legend.text=element_text(size=14));
  
  plots[[cutoff]]<-ggplot(combine2,aes(x=edge.common,y=prop,
                                       #shape=factor(status),
                  size=factor(status),col=factor(status),group=factor(status)))+
    geom_jitter()+theme_bw()+
    scale_x_continuous(breaks=c(1,5*(1:24)),labels=c(1,5*(1:24)))+
    #scale_y_continuous(trans = scales::pseudo_log_trans(1e-07, 10),
    #                   breaks=c(0, 1e-07,1e-06,1e-05,1e-04,1e-03,0.01, 0.1, 1))+
    scale_y_continuous(trans = scales::pseudo_log_trans(1e-05, 10),
                       breaks=c(0, 1e-05,1e-04,1e-03,0.01, 0.1, 1))+
    #scale_x_log10()+
    #ggplot(combine,aes(x=edge.common,y=Freq,col=factor(status)))+
    #ggplot(subset(combine,edge.common!=0),aes(x=edge.common,y=Freq,col=factor(status)))+
    
    xlab('Edge commonality')+ylab('Frequency of edges')+
    #scale_shape_manual(values=c(16, 17))+
    #scale_x_discrete(labels=xlabels2)+
    scale_color_manual(name='',values=c('#E69F00','#999999',"#56B4E9"),
                       labels=c('obs','simu'))+
    scale_size_manual(values=c(2,1,2),guide = FALSE)+
    #ggtitle(paste0('top',cutoff,', JSD=',dist))+
    theme(axis.text=element_text(size=14),
          axis.text.x=element_text(size=12,angle=0,hjust=0.5),
          axis.title=element_text(size=14,face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title=element_text(size=14,face="bold"),
          legend.text=element_text(size=14));
}      

length(plots);length(plots.line)

plots0<-plots.bar[as.character(sort(as.numeric(names(plots.bar)),decreasing = F))]
plots1<-plots.line[as.character(sort(as.numeric(names(plots.line)),decreasing = F))]
plots2<-plots[as.character(sort(as.numeric(names(plots)),decreasing = F))]
names(plots2)
plots1[[1]]
plots2[[1]]

if(T){ #only for select cutoff
pdf(paste0(out.dir,'/edge.common_obs_vs_simu_',topK,'.pdf'),useDingbats = T,width = 9,height = 6)
print(plots2[[1]])
#print(plots1[[4]]);print(plots2[[4]])
dev.off()
}

if(F){ #for all cutoffs
pdf(paste0(out.dir,'/edge.common_obs_vs_simu.pdf'),width=20,heigh=25,useDingbats = T);
#grid.arrange(grobs=plots0,ncol=2)
grid.arrange(plots0[[1]],plots2[[1]],
             plots0[[2]],plots2[[2]],
             plots0[[3]],plots2[[3]],
             plots0[[4]],plots2[[4]],
#             ncol=2,nrow=4)
#grid.arrange(
             plots0[[5]],plots2[[5]],
             plots0[[6]],plots2[[6]],
             plots0[[7]],plots2[[7]],
             #plots0[[8]],plots2[[8]],
             ncol=2,nrow=7)

dev.off()
}


