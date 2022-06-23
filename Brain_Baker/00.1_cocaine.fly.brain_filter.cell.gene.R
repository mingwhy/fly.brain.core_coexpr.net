
dat=readRDS('Brain.control.rds');
raw=dat@assays$RNA@counts
cell.annotation.df=dat@meta.data
table(cell.annotation.df$annotation)

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
  # 0 mito genes
  
  # 2) remove mito gene 
  
  # 3) filter cell: gene(<200), umi(<400)
  x=Matrix::colSums(raw>0)
  i= x < nGeneLowCutOff
  sum(i) #0, all cells express at least 200 genes 
  cat(sum(i),'cells fail nGeneLowCutOff\n');
  # 0 cells fail nGeneLowCutOff
 
  
  x=Matrix::colSums(raw)
  sum(x<nUMILowCutOff) #0 cells, total UMI<400
  i=x<nUMILowCutOff
  cat(sum(i),'cells fail nUMILowCutOff\n');
  #3 cells fail nUMILowCutOff

  if(sum(i) !=0){
    raw=raw[,!i]    
    cell.annotation.df=cell.annotation.df[!i,]
  }
  
  # 4) filter gene:  discard genes which express <4 cells
  x<-Matrix::rowSums(raw>0)
  sum(x==0) #0 gene
  sum(x<min.cell) #563 genes, express in less than 4 cells
  i=x<min.cell
  
  cat(sum(i),'genes fail nUMILowCutOff\n');
  #504 genes fail nUMILowCutOff

  if(sum(i)!=0){
    raw<-raw[!rownames(raw) %in% names(which(x<min.cell)),]
    dim(raw) # 10386 41520    
  }
  
  
  whole<- CreateSeuratObject(counts = raw,
                             min.cells=0, min.features = 0, project = "whole")
  whole #12616 features across 100267 samples
  head(rownames(whole))
  head(colnames(whole))
  
  x=sum(colnames(whole)==rownames(cell.annotation.df)) #37254
  cat('check barcode match,ncell',ncol(whole),x,'matched\n')
  
  #cell.annotation.df$annotation=cell.annotation.df$Celltype;
  tmp=stringr::str_extract(cell.annotation.df$celltype.gender_stim,'_.*_')
  sex=gsub('_','',tmp)
  cell.annotation.df$sex=sex;
  
  whole@meta.data=cell.annotation.df
  table(whole$sex)
  #female   male 
  #20236  23585 
  length(table(whole$annotation)) #39
  table(whole$sex,whole$annotation)
  
  saveRDS(whole,'./Brain.control_filtered.rds')
  

