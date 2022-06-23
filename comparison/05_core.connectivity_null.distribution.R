
# for a given cutoff, compare the stat between obs and permu

options(stringsAsFactors = F)
library(igraph)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(RColorBrewer)
library(readxl)
library(org.Dm.eg.db,verbose=F,quietly=T)

out.path='connectivity_null';
dir.create(out.path)

files1=c('~/Documents/fly.brain_coexpr.net_2022.06/Brain_Davie/brain_scRNA-seq_n15c0.005_bigScale2/top_0.05_pan.rds',
        '~/Documents/fly.brain_coexpr.net_2022.06/Brain_Baker/brain_scRNA-seq_n15c0.005_bigScale2/top_0.05_pan.rds',
        '~/Documents/fly.brain_coexpr.net_2022.06/FCA_head/brain_snRNA-seq_n15c0.005_bigScale2/top_0.05_pan.rds',
        '~/Documents/fly.brain_coexpr.net_2022.06/FCA_body/brain_snRNA-seq_n15c0.005_bigScale2/top_0.05_pan.rds')

cutoffs=data.table::fread('core_cutoffs.txt')

(mycol=c("#C4961A",'#E7298A',"#52854C", "#4E84C4")) #"#D16103" red
barplot(1:4,col=mycol)
datasets=c('brain, Davie et al. 2018','brain, Baker et al. 2021',
           'head, Li et al. 2022','body, Li et al. 2022')


#files=Sys.glob('*et al.*.txt')
#(files=files[c(3,2,4,1)])

n.rep=1000;
####################################################################
## read in pan-networks at different cutoffs 
for(ii in 1:4){
  cut=cutoffs[[ii]]
  file=files1[ii];
  cat(file,cut)
  
  pan=readRDS(file)
  diag(pan)=0;
  pan[is.na(pan)]=0; #there is NA in the pan
  dim(pan)
  i=Matrix::rowSums(pan,na.rm=T)
  sum(i==0)
  if(sum(i==0)>0){pan=pan[i!=0,i!=0]}
  pan.genes=names(which(i>0))
  length(pan.genes)
  
  #obs 
  net=pan>=cut; #True or False binary matrix
  net0=net; #unweight
  #net0= pan * net ; #weight
  diag(net0)=FALSE;
  i=apply(net0,1,sum)
  sum(i==0);length(i)
  net1=net0[i!=0,i!=0]
  dim(net1) 
  dim(net1);sum(net1)/2;
  cat('nlink.cut',cut,'ngene',nrow(net1),'\n')
  #cutoff 0.05 nlink.cut 12 ngene 1428 
  
  #gene degree of net1 genes
  x=pan[rownames(net1),]
  core.gene.degree=apply(x,1,function(i){sum(i>0)})      
  all.gene.degree=apply(pan,1,function(i){sum(i>0)})
  sum(all.gene.degree %in% core.gene.degree) #381
  gene.degree.list=sapply(core.gene.degree,function(i){names(which(all.gene.degree==i))})
  sum(sapply(gene.degree.list,length)==0)
  
  g<-igraph::graph_from_adjacency_matrix(net1,mode='undirected',weighted = NULL,diag=F) 
  
  ngene=length(V(g))
  nlink=length(E(g));#sum(net[upper.tri(net)])
  gene.degree=degree(g)
  mean.gene.degree=mean(gene.degree)
  
  wtc <- cluster_walktrap(g,  steps = 4);
  (modular<-modularity(wtc))
  
  (cluster.coeff<-transitivity(g, type = c("undirected"))); #"undirected",This is the same as global.
  #0.4950044
  obs<-c(ngene, nlink, mean.gene.degree,
         modular, cluster.coeff);
  
  # begin simu
  net=pan; #True or False binary matrix
  net0=net; #unweight
  #net0= pan * net ; #weight
  diag(net0)=FALSE;
  i=apply(net0,1,sum)
  sum(i==0);length(i)
  net1=net0[i!=0,i!=0]
  dim(net1) # 2088 2088
  dim(net1);sum(net1)/2;#7299114 edges among 2088 genes 
  net1[net1>0]=1 #change into 0/1 binary matrix
  sum(net1)
  
  (simu.file=paste0(out.path,'/null_',datasets[ii],'.rds'))
  if(!file.exists(simu.file)){
    simu.out=as.numeric();
    for(i.rep in 1:n.rep){
      signal=1;
      while(signal){
        pick.genes=c();
        for(j in gene.degree.list){
          if(length(j)==1){pick.genes=c(pick.genes,j);next}
          pick.genes=c(pick.genes,sample(j,1))
        }
        net.sub=net1[pick.genes,pick.genes] #you may sample genes with no connection with other sampled genes
        if(sum(apply(net.sub,1,sum)==0)==0){signal=0}
      }
      sum(net.sub[upper.tri(net.sub)]>0) # number of edges among these genes
      #sample nlink edges with all nodes connected at least once. 
      #two step sample:1) pick one edge for each node. 2) sample(nlink-nnode) edges
      net.for.sample=net.sub;net.for.sample[net.for.sample!=0]=0
      for(i in 1:nrow(net.sub)){
        tmp=sample(which(net.sub[i,]>0),1)
        net.for.sample[i,tmp]=1;
        net.for.sample[tmp,i]=1 #symmetric
      }
      #isSymmetric(net.sub);isSymmetric(net.for.sample)
      net.remain=net.sub- net.for.sample
      #unique(as.numeric(net.remain)) #only 0 1
      edge.remain.sample=nlink-sum(net.for.sample)/2
      g<-igraph::graph_from_adjacency_matrix(net.remain,mode='undirected',weighted = NULL,diag=F) 
      df.g=igraph::as_data_frame(g)
      nrow(df.g)  
      pick.edges=sample(1:nrow(df.g),edge.remain.sample,replace = F)
      df.g.pick=df.g[pick.edges,]
      
      g0<-igraph::graph_from_adjacency_matrix(net.for.sample,mode='undirected',weighted = NULL,diag=F) 
      df.g0=igraph::as_data_frame(g0)
      nrow(df.g0)  
      
      if(nrow(df.g0)+nrow(df.g.pick)!=nlink){cat("simu edge!=obs edge number\n")}
      df.combine=rbind(df.g0,df.g.pick)
      g<-igraph::graph_from_data_frame(df.combine,directed=F) 
      cc<-transitivity(g, type = c("undirected")); #"undirected",This is the same as global.
      simu.out[[i.rep]]=cc
    }
    saveRDS(simu.out,simu.file)
  }
  
  simu.out=readRDS(simu.file)
  summary(simu.out) 
  #Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #0.1581  0.1687  0.1716  0.1718  0.1747  0.1854 
  
  df=data.frame(status=rep('simu',length(simu.out)),cc=simu.out)
  
  pdf(paste0(paste0(out.path,'/null_',datasets[ii],'.pdf')),width=9)
  print(ggplot(df,aes(x=cc))+geom_histogram(bins=400)+theme_bw()+
          geom_vline(xintercept = obs[5], #linetype="dotted", 
                     color = "red", size=1)+
          ylab('Frequency')+xlab('Clustering Coefficient')+
          theme(panel.grid = element_blank(),
                axis.text=element_text(size=12,angle = 0),
                axis.title=element_text(size=12,angle = 0,face="bold")))
  dev.off()  
}

(files=Sys.glob('connectivity_null/*.rds'))
(files=files[c(3,2,4,1)])
plots=lapply(1:4,function(ii){
  simu.out=readRDS(files[[ii]])
  summary(simu.out) 
  df=data.frame(status=rep('simu',length(simu.out)),cc=simu.out)
 
  ggplot(df,aes(x=cc))+geom_histogram(bins=400)+theme_bw()+
          geom_vline(xintercept = obs[5], #linetype="dotted", 
                     color = "red", size=1)+
          ylab('Frequency')+xlab('Clustering Coefficient')+
          ggtitle(datasets[[ii]])+
          theme(panel.grid = element_blank(),
                axis.text=element_text(size=12,angle = 0),
                axis.title=element_text(size=12,angle = 0,face="bold"))

})
pdf(paste0(out.path,'/null_obs.pdf'),useDingbats = T)
grid.arrange(grobs=plots,ncol=2)
dev.off()
