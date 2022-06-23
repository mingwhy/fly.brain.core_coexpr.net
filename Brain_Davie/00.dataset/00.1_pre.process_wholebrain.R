
library(Seurat)
library(ggplot2)

raw<-readRDS("./wholebrain-raw.rds")
raw #17473 features across 56234 samples within 1 assay 

## filter condition
nGeneLowCutOff <- 200; nGeneHighCutOff <- Inf
nUMILowCutOff <- 500; nUMIHighCutOff <- Inf
MitoLowCutOff <- -Inf; MitoHighCutOff <- 0.30
min.cell <- 4; #a gene expresses at >=4 cells

summary(raw@meta.data$nFeature_RNA)
summary(raw@meta.data$nCount_RNA)
summary(raw@meta.data$prop.mito)

## filter step: 1) calculate mito% for each cell, remove cell whose mito.perc >= 0.3
# 2) remove mito gene
# 3) filter cell: gene, umi
# 4) filter gene:  discard genes which express <4 cells

# 1) calculate mito% for each cell, remove cell whose mito.perc >= 0.3
(mito.genes <- grep(pattern = "^mt:",ignore.case = T,x = rownames(raw), value = TRUE))
# 37 mito genes
prop.mito <- Matrix::colSums(raw@assays$RNA@counts[mito.genes, ]) / Matrix::colSums(raw@assays$RNA@counts)
summary(prop.mito)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.01994 0.03226 0.03807 0.05036 0.13491
summary(mydata@meta.data$percent.mito) #my calculation is similar to meta.data in the original publication
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.02001 0.03248 0.03823 0.05064 0.13559 

# 2) remove mito gene
(mito.genes <- grep(pattern = "^mt:",ignore.case = T,x = rownames(raw), value = TRUE))
# 37 mito.genes
sum(rownames(raw) %in% mito.genes)
#[1] 37
raw_rmMito<-raw[!rownames(raw) %in% mito.genes,]
raw_rmMito #17436 features across 56234 samples within 1 assay 

# 3) filter cell: gene, umi

x=Matrix::colSums(raw_rmMito@assays$RNA@counts>0)
sum(x<nGeneLowCutOff) #0, all cells express at least 200 genes 
x=Matrix::colSums(raw_rmMito@assays$RNA@counts)
sum(x<nUMILowCutOff) #42 cells, total UMI<500
# the filter below, only fliter cells, 42 cells were filtered out

raw_filter <- subset(x = raw_rmMito, 
                     subset= nFeature_RNA >=nGeneLowCutOff & nFeature_RNA <=nGeneHighCutOff &
                       nCount_RNA>=nUMILowCutOff & nCount_RNA<= nUMIHighCutOff &
                       prop.mito>=MitoLowCutOff & prop.mito<=MitoHighCutOff);
raw_filter  #17436 features across 56192 samples within 1 assay 

# 4) filter gene:  discard genes which express <4 cells
x<-Matrix::rowSums(raw_filter@assays$RNA@counts>0)
sum(x==0) #3018 gene
sum(x<min.cell) #4934 genes

mydata<-raw_filter[!rownames(raw_filter) %in% names(which(x<min.cell)),]
mydata #12502 features across 56192 samples within 1 assay 
head(mydata@meta.data)


saveRDS(mydata,'./wholebrain_filtered.rds')

