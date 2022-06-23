
# run on server
library(Seurat);
library(ggplot2);library(gridExtra);
library(parallel)
source('../../modify_bigScale2/ming-Functions.R') 

my.method='bigscale'
prefix='sn'
files=Sys.glob('common.genes_brain_snRNA-seq_n15c0.005/*.rds')

thresholds<- c(0.7,0.9,0.95,0.99,0.995)
n.rep=10;
(out.dir=paste0(prefix,'_',my.method,'_rep',n.rep,'/'))
dir.create(out.dir)

splitData <- function(df, colIdx,propShared=0.5){
  sampleCount <- ncol(df)
  #colIdx <- sample.int(sampleCount)
  sharedSampleCount <- round(propShared*sampleCount)
  sharedColsIdx <- 1:sharedSampleCount
  specificSampleCount <- floor(0.5*(sampleCount-sharedSampleCount))
  
  s1=sharedSampleCount+1;
  e1=sharedSampleCount+1+specificSampleCount;
  s2=e1+1;
  e2=min(s2+specificSampleCount,sampleCount);
  df1 <- df[,colIdx[c(sharedColsIdx, s1:e1) ]]
  df2 <- df[,colIdx[c(sharedColsIdx, s2:e2)]]
  dfList <- list(df1, df2)
  return(dfList)
}
################################################

for(file in files){
  mat.list=readRDS(file);
  #mat.scran=as.matrix(mat.list$mat.scran);
  mat.count=as.matrix(mat.list$mat.count);
  if(ncol(mat.count)==0){next}
  cat('begin file',file,'\n');
  
  set.seed(221013)
  splitData.mat=replicate(n.rep,sample.int(ncol(mat.count)))
  saveRDS(splitData.mat,paste0(out.dir,'splitData.mat_',basename(file)) )
  
  ## original data, subsample
  inp.mat=mat.count
  gene.names=rownames(inp.mat);
  mclapply(1:n.rep,function(i){
  #for(i in 1:n.rep){
    (out.file=paste0(out.dir,'Original_rep',i,'_',basename(file)))
    if(file.exists(out.file)) next 
    colIdx=splitData.mat[,i];
    dat=splitData(inp.mat, colIdx,propShared=0.5)
    net.out1=compute.network(expr.data = dat[[1]],gene.names = gene.names,speed.preset = "fast",clustering='recursive')
    net.out2=compute.network(expr.data = dat[[2]],gene.names = gene.names,speed.preset = "fast",clustering='recursive')
    out=list(net.out1=net.out1,net.out2=net.out2);
    saveRDS(out,out.file)}, mc.cores=10)
  #}
}



