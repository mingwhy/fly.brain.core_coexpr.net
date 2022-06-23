## data source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107451

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)
options(future.globals.maxSize = 3.2 * 1024^3,stringsAsFactors = F)

# Load Data
df.expr=readMM("../00.raw_data/GSE107451_DGRP-551_w1118_WholeBrain_57k_0d_1d_3d_6d_9d_15d_30d_50d_10X_DGEM_MEX.mtx.tsv/matrix.mtx")
dim(df.expr); #17473gene x 56902cell
class(df.expr)
#df.expr[1:3,1:2]

barcode=scan("../00.raw_data/GSE107451_DGRP-551_w1118_WholeBrain_57k_0d_1d_3d_6d_9d_15d_30d_50d_10X_DGEM_MEX.mtx.tsv/barcodes.tsv",
             what="")
head(barcode)


## cell type info
df.cell=read.table("../00.raw_data/GSE107451_DGRP-551_w1118_WholeBrain_57k_Metadata.tsv",as.is=T,header=T);
table(df.cell$annotation)
sum(df.cell$annotation=='Hsp') #668 cells annotated to 'Hsp'
# remove them due to personal email communication, they are 'stressed' cells.

## gene ID mapping
df.gene=read.table("../00.raw_data/GSE107451_DGRP-551_w1118_WholeBrain_57k_0d_1d_3d_6d_9d_15d_30d_50d_10X_DGEM_MEX.mtx.tsv/genes.tsv",as.is=T);
dim(df.gene) #17473 genes
colnames(df.gene)=c("flybase",'symbol')

dim(df.cell);length(barcode); #56902cell
dim(df.expr); #17473gene x 56902cell
saveRDS(list(cell.annotation.df=df.cell,genes.df=df.gene$flybase), 
        'raw_stats.rds')


## generate Seurat object
dim(df.cell)
sum(barcode==df.cell$new_barcode)
rownames(df.expr)=df.gene$symbol
colnames(df.expr)=df.cell$new_barcode

df.expr2=df.expr[,!df.cell$new_barcode %in% df.cell[df.cell$annotation == 'Hsp',]$new_barcode] #remove Hsp cells
dim(df.expr2) #17473 56234
df.cell2=df.cell[df.cell$annotation!='Hsp',]

whole_brain<- CreateSeuratObject(counts = df.expr2,
          min.cells=0, min.features = 0, project = "wholebrain")

sum(rownames(whole_brain@meta.data)==df.cell2$new_barcode)#56234
dim(whole_brain)
whole_brain@meta.data=cbind(whole_brain@meta.data,df.cell2)#modify meta.data

table(whole_brain@meta.data$annotation)

whole_brain #117473 features across 56234 samples within 1 assay 
saveRDS(whole_brain,'./wholebrain-raw.rds')

