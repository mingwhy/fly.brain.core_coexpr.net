################################################
#Data source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152495

#GSM4616550	Female Sucrose Rep 1
#GSM4616551	Female Sucrose Rep 2
#GSM4616552	Male Sucrose Rep 1
#GSM4616553	Male Sucrose Rep 2
#GSM4616554	Female Cocaine Rep 1
#GSM4616555	Female Cocaine Rep 2
#GSM4616556	Male Cocaine Rep 1
#GSM4616557	Male Cocaine Rep 2

#Original paper code source: https://github.com/vshanka23/The-Drosophila-Brain-on-Cocaine-at-Single-Cell-Resolution/blob/master/Rcode_for_analyses.R
################################################

library(cowplot)
library(ggplot2)
library(dplyr)
library(Seurat)
library(rlang)
#library(TopKLists)
library(readxl)
library(scAlign)
##scTransform requires considerable RAM so add this statement before proceeding to intergration.
options(future.globals.maxSize= 53687091200)

##Read data in
mtx <- "GSE152495_RAW/GSM4616550_1_S1_L001_matrix.mtx.gz"
cells <- "GSE152495_RAW/GSM4616550_1_S1_L001_barcodes.tsv.gz"
features <- "GSE152495_RAW/GSM4616550_1_S1_L001_features.tsv.gz"
Female_sucrose_R1_Data  <- ReadMtx(mtx = mtx, cells = cells, features = features)

mtx <- "GSE152495_RAW/GSM4616551_2_S2_L001_matrix.mtx.gz"
cells <- "GSE152495_RAW/GSM4616551_2_S2_L001_barcodes.tsv.gz"
features <- "GSE152495_RAW/GSM4616551_2_S2_L001_features.tsv.gz"
Female_sucrose_R2_Data <- ReadMtx(mtx = mtx, cells = cells, features = features)

mtx <- "GSE152495_RAW/GSM4616552_3_S3_L001_matrix.mtx.gz"
cells <- "GSE152495_RAW/GSM4616552_3_S3_L001_barcodes.tsv.gz"
features <- "GSE152495_RAW/GSM4616552_3_S3_L001_features.tsv.gz"
Male_sucrose_R1_Data <- ReadMtx(mtx = mtx, cells = cells, features = features)

mtx <- "GSE152495_RAW/GSM4616553_4_S4_L001_matrix.mtx.gz"
cells <- "GSE152495_RAW/GSM4616553_4_S4_L001_barcodes.tsv.gz"
features <- "GSE152495_RAW/GSM4616553_4_S4_L001_features.tsv.gz"
Male_sucrose_R2_Data <- ReadMtx(mtx = mtx, cells = cells, features = features)

mtx <- "GSE152495_RAW/GSM4616554_5_S5_L001_matrix.mtx.gz"
cells <- "GSE152495_RAW/GSM4616554_5_S5_L001_barcodes.tsv.gz"
features <- "GSE152495_RAW/GSM4616554_5_S5_L001_features.tsv.gz"
Female_cocaine_R1_Data <- ReadMtx(mtx = mtx, cells = cells, features = features)

mtx <- "GSE152495_RAW/GSM4616555_6_S6_L001_matrix.mtx.gz"
cells <- "GSE152495_RAW/GSM4616555_6_S6_L001_barcodes.tsv.gz"
features <- "GSE152495_RAW/GSM4616555_6_S6_L001_features.tsv.gz"
Female_cocaine_R2_Data <- ReadMtx(mtx = mtx, cells = cells, features = features)


mtx <- "GSE152495_RAW/GSM4616556_7_S7_L001_matrix.mtx.gz"
cells <- "GSE152495_RAW/GSM4616556_7_S7_L001_barcodes.tsv.gz"
features <- "GSE152495_RAW/GSM4616556_7_S7_L001_features.tsv.gz"
Male_cocaine_R1_Data <- ReadMtx(mtx = mtx, cells = cells, features = features)

mtx <- "GSE152495_RAW/GSM4616557_8_S8_L001_matrix.mtx.gz"
cells <- "GSE152495_RAW/GSM4616557_8_S8_L001_barcodes.tsv.gz"
features <- "GSE152495_RAW/GSM4616557_8_S8_L001_features.tsv.gz"
Male_cocaine_R2_Data <- ReadMtx(mtx = mtx, cells = cells, features = features)


##Create Seurat objects
Female_sucrose_R1 <- CreateSeuratObject(counts=Female_sucrose_R1_Data,project="Female_Sucrose",min.cells=5)
Female_sucrose_R2 <- CreateSeuratObject(counts=Female_sucrose_R2_Data,project="Female_Sucrose",min.cells=5)
Male_sucrose_R1 <- CreateSeuratObject(counts=Male_sucrose_R1_Data,project="Male_Sucrose",min.cells = 5)
Male_sucrose_R2 <- CreateSeuratObject(counts=Male_sucrose_R2_Data,project="Male_Sucrose",min.cells = 5)
Female_cocaine_R2 <- CreateSeuratObject(counts=Female_cocaine_R2_Data,project="Female_Cocaine",min.cells=5)
Female_cocaine_R1 <- CreateSeuratObject(counts=Female_cocaine_R1_Data,project="Female_Cocaine",min.cells=5)
Male_cocaine_R1 <- CreateSeuratObject(counts=Male_cocaine_R1_Data,project="Male_Cocaine",min.cells = 5)
Male_cocaine_R2 <- CreateSeuratObject(counts=Male_cocaine_R2_Data,project="Male_Cocaine",min.cells = 5)

##Set identities for condition
Male_cocaine_R1$stim <- "cocaine"
Male_cocaine_R2$stim <- "cocaine"
Male_sucrose_R1$stim <- "sucrose"
Male_sucrose_R2$stim <- "sucrose"
Female_cocaine_R1$stim <- "cocaine"
Female_cocaine_R2$stim <- "cocaine"
Female_sucrose_R1$stim <- "sucrose"
Female_sucrose_R2$stim <- "sucrose"

##Set identities for gender and condition
Male_cocaine_R1$gender_stim <- "male_cocaine"
Male_cocaine_R2$gender_stim <- "male_cocaine"
Male_sucrose_R1$gender_stim <- "male_sucrose"
Male_sucrose_R2$gender_stim <- "male_sucrose"
Female_cocaine_R1$gender_stim <- "female_cocaine"
Female_cocaine_R2$gender_stim <- "female_cocaine"
Female_sucrose_R1$gender_stim <- "female_sucrose"
Female_sucrose_R2$gender_stim <- "female_sucrose"

##Set identities for sample
Male_cocaine_R1$sample_id <- "male_cocaine_R1"
Male_cocaine_R2$sample_id <- "male_cocaine_R2"
Male_sucrose_R1$sample_id <- "male_sucrose_R1"
Male_sucrose_R2$sample_id <- "male_sucrose_R2"
Female_cocaine_R1$sample_id <- "female_cocaine_R1"
Female_cocaine_R2$sample_id <- "female_cocaine_R2"
Female_sucrose_R1$sample_id <- "female_sucrose_R1"
Female_sucrose_R2$sample_id <- "female_sucrose_R2"

##Remove spurious features (low and high)
Male_cocaine_R2 <- subset(Male_cocaine_R2,subset=nFeature_RNA > 300 & nFeature_RNA < 2500)
Male_cocaine_R1 <- subset(Male_cocaine_R1,subset=nFeature_RNA > 300 & nFeature_RNA < 2500)
Male_sucrose_R2 <- subset(Male_sucrose_R2,subset=nFeature_RNA > 300 & nFeature_RNA <2500)
Male_sucrose_R1 <- subset(Male_sucrose_R1,subset=nFeature_RNA > 300 & nFeature_RNA <2500)
Female_cocaine_R2 <- subset(Female_cocaine_R2,subset=nFeature_RNA > 300 & nFeature_RNA < 2500)
Female_cocaine_R1 <- subset(Female_cocaine_R1,subset=nFeature_RNA > 300 & nFeature_RNA < 2500)
Female_sucrose_R2 <- subset(Female_sucrose_R2,subset=nFeature_RNA > 300 & nFeature_RNA < 2500)
Female_sucrose_R1 <- subset(Female_sucrose_R1,subset=nFeature_RNA > 300 & nFeature_RNA < 2500)

##Normalize
Male_cocaine_R1 <- SCTransform(Male_cocaine_R1,verbose = FALSE,return.only.var.genes = FALSE)
Male_cocaine_R2 <- SCTransform(Male_cocaine_R2,verbose = FALSE,return.only.var.genes = FALSE)
Male_sucrose_R1 <- SCTransform(Male_sucrose_R1,verbose = FALSE,return.only.var.genes = FALSE)
Male_sucrose_R2 <- SCTransform(Male_sucrose_R2,verbose = FALSE,return.only.var.genes = FALSE)
Female_cocaine_R1 <- SCTransform(Female_cocaine_R1,verbose = FALSE,return.only.var.genes = FALSE)
Female_cocaine_R2 <- SCTransform(Female_cocaine_R2,verbose = FALSE,return.only.var.genes = FALSE)
Female_sucrose_R1 <- SCTransform(Female_sucrose_R1,verbose = FALSE,return.only.var.genes = FALSE)
Female_sucrose_R2 <- SCTransform(Female_sucrose_R2,verbose = FALSE,return.only.var.genes = FALSE)

##Feature selection
Brain.features <- SelectIntegrationFeatures(object.list=c(Male_cocaine_R1,Male_cocaine_R2,Male_sucrose_R1,Male_sucrose_R2,Female_cocaine_R1,Female_cocaine_R2,Female_sucrose_R1,Female_sucrose_R2),nfeatures=1500)
Brain.list <- PrepSCTIntegration(object.list = c(Male_cocaine_R1,Male_cocaine_R2,Male_sucrose_R1,Male_sucrose_R2,Female_cocaine_R1,Female_cocaine_R2,Female_sucrose_R1,Female_sucrose_R2),anchor.features = Brain.features,verbose = FALSE)

##Identify anchors
Brain.anchors <- FindIntegrationAnchors(object.list = Brain.list, normalization.method = "SCT", anchor.features = Brain.features, verbose = FALSE)
Brain.integrated <- IntegrateData(anchorset = Brain.anchors, normalization.method = "SCT", verbose = FALSE)
saveRDS(Brain.integrated,'Brain.integrated.rds')
#An object of class Seurat 
#23155 features across 86223 samples within 3 assays 
#Active assay: integrated (1500 features, 1500 variable features)

#########################################################

Brain.integrated=readRDS('Brain.integrated.rds')
Brain.integrated <- RunPCA(object = Brain.integrated, verbose = FALSE)


## Use elbowplot to determine how many dimensions to use for capturing as much variability as possible without using all of the axes. Look for a point of desaturation (elbow).
pdf('test.pdf')
ElbowPlot(Brain.integrated)
dev.off()

##Use info from elbow plot to set number of dimensions for UMAP
Brain.integrated <- RunUMAP(Brain.integrated, reduction = "pca", dims = 1:15)
Brain.integrated <- FindNeighbors(Brain.integrated,reduction = "pca", dims = 1:15)
brain_clusters <- FindClusters(Brain.integrated,resolution = 0.8) #39 clusters
table(Idents(brain_clusters))
#   0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
#9241 8361 7044 6260 5993 5550 3274 3153 2538 2348 2286 2228 2178 2132 2011 1916 
#  16   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31 
#1861 1845 1562 1512 1500 1343 1309 1293 1198  933  668  608  583  515  496  483 
#  32   33   34   35   36   37   38 
# 466  403  403  243  208  168  110


ncol(Brain.integrated)
sum(names(Idents(brain_clusters))==colnames(Brain.integrated))

cell.annotation.df=Brain.integrated@meta.data
cell.annotation.df$annotation=Idents(brain_clusters)
genes=rownames(Brain.integrated)
saveRDS(list(cell.annotation.df=cell.annotation.df,genes.df=genes), 
        'raw_cocaine_sucrose_stats.rds')


#########################################################
##Create new identities for DE analysis downstream
Idents(Brain.integrated)<-Idents(brain_clusters)
Brain.integrated$celltype.stim <- paste(Idents(Brain.integrated),Brain.integrated$stim,sep="_")
Brain.integrated$celltype.gender_stim <- paste(Idents(Brain.integrated),Brain.integrated$gender_stim,sep="_")

table(Brain.integrated$celltype.gender_stim)
table(Brain.integrated$stim)
#cocaine sucrose 
#  42399   43824 
saveRDS(Brain.integrated,'Brain.integrated_cluster.rds')



## only keep 4 control samples
Brain.integrated=readRDS('Brain.integrated_cluster.rds')
Brain.control=subset(Brain.integrated,stim=='sucrose')
head(Brain.control$celltype.gender_stim)
all.clusters=unique(Brain.control$celltype.stim)
length(all.clusters)#39
all.clusters=unique(Brain.control$celltype.gender_stim)
length(all.clusters)#78


#pick.clusters=names(which(table(Brain.control$celltype.gender_stim)>=200)) #54
#Brain.control.pick=subset(Brain.control, celltype.gender_stim %in% pick.clusters)
##23155 features across 41523 samples within 3 assays 
Brain.control.pick=Brain.control;
names(Brain.control.pick@assays)
#[1] "RNA"        "SCT"        "integrated"
dim(Brain.control.pick@assays$RNA@counts)
#[1] 10949 43824  #[1] 10949 41523
dim(Brain.control.pick@assays$SCT@data)
#[1] 10706 43824  #[1] 10706 41523

umi.count=Brain.control.pick@assays$RNA@counts;
meta.data=Brain.control.pick@meta.data
meta.data$annotation=meta.data$celltype.stim
sex=gsub('_sucrose','',meta.data$gender_stim)
table(sex)
meta.data$sex=sex
#table(meta.data$annotation)

Brain.control.data <- CreateSeuratObject(counts=umi.count,project="Brain.control",min.cells=0)
Brain.control.data@meta.data=meta.data
table(Brain.control.data$annotation)

Brain.control.data
#10949 features across 43824 samples within 1 assay 
#10949 features across 41523 samples within 1 assay
saveRDS(list(cell.annotation.df=Brain.control.data@meta.data,genes.df=rownames(Brain.control.data)), 
        'raw_stats.rds')
saveRDS(Brain.control.data,'Brain.control.rds')

