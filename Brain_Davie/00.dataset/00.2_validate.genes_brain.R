
library(Seurat)
library(Matrix)

(files=Sys.glob('whole*_filtered.rds'))

## validate gene id
## FB2021_05, released October 15, 2021
all.genes=c()
for(file in files){
  dat=readRDS(file)
  genes=rownames(dat)
  all.genes=c(all.genes,genes)
}
gene.pool=unique(all.genes)

#gene.pool=rownames(dat)
length(gene.pool) #12502 genes
data.table::fwrite(as.data.frame(gene.pool),'gene.pool.txt',col.names = F)

# upload 'gene.pool.txt' to flybase ID validator
# download via 'validation table' button
df=data.table::fread('FlyBase_IDs.txt')
head(df)
colnames(df)=c('submitted_item','validated_id','current_symbol');
dim(df) #13443   

# effective genes with FBgnxxx ID, throw away FBalxxx ID
df1=df[grep('FBgn',df$validated_id),]
dim(df1) #13005 


keep.genes=list();

# begin select genes
x=table(df1$submitted_item)
unique.hits=names(which(x==1)) #done
length(unique.hits) #12094 genes have unique hits
keep.genes[[1]]=df1[df1$submitted_item %in% unique.hits,]

# handle genes with >1 hits
genes=names(which(x>1));
length(genes); #408 genes with >1 hit
tmp=df1[df1$submitted_item %in% genes,]
i=tmp$submitted_item==tmp$current_symbol
sum(i) #382 genes, submitted_item and current_symbol are consistent
safe.hits=tmp[tmp$submitted_item==tmp$current_symbol,]
dim(safe.hits) #382
keep.genes[[2]]=safe.hits


# handle the rest
tmp=tmp[!tmp$submitted_item %in% safe.hits$submitted_item,]
dim(tmp) #53 rows
data.table::fwrite(tmp,file='duplicated.genes.csv')

# manual check, save some genes by adding a column of 'select'
tmp2=readxl::read_excel('after.manual_duplicated.genes.xlsx')
head(tmp2)
tmp2=tmp2[tmp2$select==1,]
dim(tmp2) #keep 2 genes

keep.genes[[2]]=tmp2[,-4]

df.keep=Reduce(`rbind`,keep.genes)
dim(df.keep) #12096 
length(unique(df.keep$submitted_item))#12096
length(unique(df.keep$current_symbol))#12095

# handle duplicated current_symbols
x=df.keep[duplicated(df.keep$current_symbol),]$current_symbol
tmp=df.keep[df.keep$current_symbol %in% x,]
tmp
# priority: submitted_item==current_symbol

df.keep1=df.keep[!df.keep$current_symbol %in% x,]
dim(df.keep1) # 12094  genes
tmp1=tmp[tmp$submitted_item==tmp$current_symbol,]
tmp1
df.keep2=rbind(df.keep1,tmp1)
dim(df.keep2) #12094  

length(unique(df.keep2$submitted_item))#12094  
length(unique(df.keep2$current_symbol))#12094  

data.table::fwrite(df.keep2,'validated_gene.symbols_brain.txt',sep='\t',quote=F,row.names=F);


###############################################################################
## go to flybase 'batch download', upload the 
#3rd column (current_symbol) of 'validated_gene.symbols.txt'
# extract chr info (Genomic Location and Detailed Mapping Data)
# download as 'FlyBase_Fields_download.txt' file
df=data.table::fread('FlyBase_Fields_download.txt')
dim(df) #12094   genes

sum(df$SYMBOL==df$`#SUBMITTED ID`) #12094  

head(df.keep2)
head(df)
df3=merge(df.keep2,df,by.x='current_symbol',by.y='#SUBMITTED ID')
dim(df3) # 12094     11
data.table::fwrite(df3,'gene.meta_brain.txt',sep='\t',quote=F,row.names=F);

gene.meta=data.table::fread('gene.meta_brain.txt',header=T,sep='\t')
dim(gene.meta) #12094   
table(gene.meta$LOCATION_ARM)

###############################################################################
## update rds with validate gene names 
files
for(file in files){
  dat=readRDS(file)
  dim(dat) #10353  3797
  genes=rownames(dat)
  
  dat1=dat[genes %in% gene.meta$submitted_item,]
  dim(dat1) #9710 1297
  
  genes=rownames(dat1)
  
  i=match(genes,gene.meta$submitted_item)
  gene.meta2=gene.meta[i,]
  if(sum(genes==gene.meta2$submitted_item)!=length(genes)){
    cat('wrong with ',file,'\n')
  }
  rownames(dat1@assays$RNA@counts)<-gene.meta2$SYMBOL
  rownames(dat1@assays$RNA@data)<-gene.meta2$SYMBOL
  dat1 #12094 genes x  56192 cells
  
  file.name=gsub('.rds','_valid.rds',file)
  saveRDS(dat1,file.name)
  cat('file',file,'is done\n')
}
