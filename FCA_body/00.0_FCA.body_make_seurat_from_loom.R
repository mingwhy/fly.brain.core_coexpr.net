## data source: https://flycellatlas.org/
## fly head: 10x, Stringent, Loom

## use SCopeLoomR to read loom file
if(F){remotes::install_github("aertslab/SCopeLoomR")}
library(SingleCellExperiment)
library(SCopeLoomR)
library(Seurat)


file='~/Documents/Data_fly_FCA/xFCA_body/s_fca_biohub_body_10x.loom';


  tissue='body'
  (out.rds=paste0('whole_',tissue,'_filtered.rds'))
  cat('begin tissue',tissue,',output file to',out.rds,'\n');
  
  # data download from https://flycellatlas.org/
  #loom_path <- './s_fca_biohub_body_10x.loom';
  loom_path <- file
  loom <- open_loom(loom_path, mode="r+")
  
  cell.annotation.all=get_cell_annotation(loom)
  dim(cell.annotation.all) # 96926   414
  
  labels=colnames(cell.annotation.all)
  tmp=cell.annotation.all[,-grep('TrackRegulonsAUC|MotifRegulonsAUC',labels)]
  colnames(tmp)
  table(tmp$batch_id)
  sort(table(tmp$batch))
  sort(table(tmp$id))  #only two samples
  sort(table(tmp$batch))==sort(table(tmp$id)) #batch <=> id
  
  #write.table(tmp,file='cell.annotation.csv',quote=F,sep='\t',row.names = F) #too big
  
  cell.annotation.df=tmp
  colnames(cell.annotation.df)
  head(cell.annotation.df$annotation)
  table(cell.annotation.df$annotation)
  length(table(cell.annotation.df$annotation)) #34 (including one unannotated)
  
  # stats
  genes=get_genes(loom)
  length(genes) #13056 genes
  dim(cell.annotation.all) #100527    348
  length(unique(cell.annotation.all$annotation)) #82
  #write.table(cell.annotation.df,file='cell.annotation.csv',quote=F,sep='\t',row.names = F) #too big
  #write.table(genes,file='genes.csv',quote=F,sep='\t',row.names = F) #too big
  saveRDS(list(cell.annotation.df=cell.annotation.df,genes.df=genes), 'raw_stats.rds')
  
  raw <- get_dgem(loom)
  raw[1:5,1:5]
  dim(raw) #15267 gene by 96926 cell
  
  close_loom(loom)
  
  ## filter condition
  nGeneLowCutOff <- 200; nGeneHighCutOff <- Inf
  nUMILowCutOff <- 500; nUMIHighCutOff <- Inf
  MitoLowCutOff <- -Inf; MitoHighCutOff <- 0.30
  min.cell <- 4; # keep gene expresses at >=4 cells
  
  ## filter step: 1) calculate mito% for each cell, remove cell whose mito.perc > 0.3
  # 2) remove mito gene
  # 3) filter cell: gene (<200), umi (<400)
  # 4) filter gene:  discard genes which express <4 cells
  
  # 1) calculate mito% for each cell, remove cell whose mito.perc > 0.3
  (mito.genes <- grep(pattern = "^mt:",ignore.case = T,x = rownames(raw), value = TRUE))
  # 24 mito genes
  prop.mito <- Matrix::colSums(raw[mito.genes, ]) / Matrix::colSums(raw)
  summary(prop.mito)
  #Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
  #0.000000 0.001161 0.002257 0.003319 0.004043 0.049934 
  summary(cell.annotation.df$percent_mito) #confirm
  
  i=prop.mito>MitoHighCutOff
  cat(sum(i),'cells fail MitoHighCutOff\n');
  if(sum(i)!=0){
    raw=raw[,!i];
    cell.annotation.df=cell.annotation.df[,!i]
  }
  
  # 2) remove mito gene
  sum(rownames(raw) %in% mito.genes) #24 mito gene
  cat(sum(rownames(raw) %in% mito.genes) ,'mito genes removed\n')
  
  raw<-raw[!rownames(raw) %in% mito.genes,]
  dim(raw) #13032 genes, 100527 cells
  
  # 3) filter cell: gene(<200), umi(<400)
  x=Matrix::colSums(raw>0)
  i= x < nGeneLowCutOff
  sum(i) #0, all cells express at least 200 genes 
  cat(sum(i),'cells fail nGeneLowCutOff\n');
  
  if(sum(i) !=0 ){
    raw=raw[,!i]
    cell.annotation.df=cell.annotation.df[!i,]
  }
  
  x=Matrix::colSums(raw)
  sum(x<nUMILowCutOff) #0 cells, total UMI<400
  i=x<nUMILowCutOff
  cat(sum(i),'cells fail nUMILowCutOff\n');
  
  if(sum(i) !=0){
    raw=raw[,!i]
    cell.annotation.df=cell.annotation.df[!i,]
  }
  
  # 4) filter gene:  discard genes which express <4 cells
  x<-Matrix::rowSums(raw>0)
  sum(x==0) #0 gene
  sum(x<min.cell) #416 genes, express in less than 4 cells
  i=x<min.cell
  
  cat(sum(i),'genes fail nUMILowCutOff\n');
  
  if(sum(i)!=0){
    raw<-raw[!rownames(raw) %in% names(which(x<min.cell)),]
    dim(raw) 
    dim(cell.annotation.df) #100267     34
  }
  
  
  whole<- CreateSeuratObject(counts = raw,
                             min.cells=0, min.features = 0, project = "whole")
  whole #15010 features across 96259 samples within 1 assay 
  head(rownames(whole))
  head(colnames(whole))
  head(rownames(cell.annotation.df))
  x=sum(colnames(whole)==rownames(cell.annotation.df)) #37254
  cat('check barcode match,ncell',ncol(whole),x,'matched\n')
  whole@meta.data=cell.annotation.df
  

  saveRDS(whole, out.rds)
  cat('tissue',tissue,'is done\n')


