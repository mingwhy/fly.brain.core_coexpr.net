
options(stringsAsFactors = F)

(files=Sys.glob('*filtered.rds'))#1 files

######################################
## validate gene id
## FB2021_05, released October 15, 2021
if(!file.exists('gene.meta.txt')){
  
  all.genes=c()
  for(file in files){
    dat=readRDS(file)
    genes=rownames(dat)
    all.genes=c(all.genes,genes)
  }
  gene.pool=unique(all.genes)
  
  #gene.pool=rownames(dat)
  length(gene.pool) #12616 genes
  #data.table::fwrite(as.data.frame(gene.pool),'gene.pool.txt',col.names = F)
  
  # upload gene.pool.txt to http://flybase.org/convert/id
  # download result as FlyBase_IDs.txt
  df=read.table('FlyBase_IDs.txt')
  colnames(df)=c('submitted_item','validated_id','current_symbol');
  dim(df) #17270 
  
  # effective genes with FBgnxxx ID, throw away FBalxxx ID
  df1=df[grep('FBgn',df$validated_id),]
  dim(df1) #16809 
  keep.genes=list();
  
  # begin select genes
  x=table(df1$submitted_item)
  unique.hits=names(which(x==1)) #done
  length(unique.hits) #15869 genes have unique hits
  keep.genes[[1]]=df1[df1$submitted_item %in% unique.hits,]
  
  # handle genes with >1 hits
  genes=names(which(x>1));
  length(genes); #421 genes with >1 hit
  tmp=df1[df1$submitted_item %in% genes,]
  i=tmp$submitted_item==tmp$current_symbol
  sum(i) #411 genes, submitted_item and current_symbol are consistent
  safe.hits=tmp[tmp$submitted_item==tmp$current_symbol,]
  dim(safe.hits) #411
  keep.genes[[2]]=safe.hits
  
  
  # handle the rest
  tmp=tmp[!tmp$submitted_item %in% safe.hits$submitted_item,]
  dim(tmp) #20 rows
  #data.table::fwrite(tmp,file='duplicated.genes.csv')
  
  # manual check or just remove these genes
  remove.genes=unique(tmp$submitted_item) 
  length(remove.genes) #10 genes to be resolved, or just removed
  
  
  df.keep=Reduce(`rbind`,keep.genes)
  dim(df.keep) #16280
  length(unique(df.keep$submitted_item))#16280
  length(unique(df.keep$current_symbol))#16276
  
  # handle duplicated current_symbols
  x=df.keep[duplicated(df.keep$current_symbol),]$current_symbol
  tmp=df.keep[df.keep$current_symbol %in% x,]
  # priority: submitted_item==current_symbol
  
  df.keep1=df.keep[!df.keep$current_symbol %in% x,]
  dim(df.keep1) #16272 genes
  tmp1=tmp[tmp$submitted_item==tmp$current_symbol,]
  df.keep2=rbind(df.keep1,tmp1)
  dim(df.keep2) #16276
  
  length(unique(df.keep2$submitted_item))#16276
  length(unique(df.keep2$current_symbol))#16276
  
  #data.table::fwrite(df.keep2,'validated_gene.symbols.txt',sep='\t',quote=F,row.names=F);
  
  ## go to flybase 'batch download', upload the 3rd column (current_symbol) of 'validated_gene.symbols.txt'
  # extract chr info (Genomic Location and Detailed Mapping Data)
  # download as 'FlyBase_Fields_download.txt' file
  df=data.table::fread('FlyBase_Fields_download.txt')
  dim(df) #16276 genes
  sum(df$SYMBOL==df$`#SUBMITTED ID`) #16276

  head(df.keep2)
  head(df)
  df3=merge(df.keep2,df,by.x='current_symbol',by.y='#SUBMITTED ID')
  dim(df3) # 16276   11  
  data.table::fwrite(df3,'gene.meta.txt',sep='\t',quote=F,row.names=F);
    
             
}

gene.meta=data.table::fread('gene.meta.txt',header=T,sep='\t')
dim(gene.meta) #16276 
length(unique(gene.meta$submitted_item))
rownames(gene.meta)=gene.meta$submitted_item

table(gene.meta$LOCATION_ARM)

###############################################################################
## update rds with validate gene names 
(files=Sys.glob('*filtered.rds'))

for(file in files){
  cat('begin file',file,'\n')
  dat=readRDS(file) #15010 features across 96259 samples within 1 assay 
  genes=rownames(dat)
  
  dat1=dat[genes %in% gene.meta$submitted_item,]
  dim(dat1)
  genes=rownames(dat1)
  
  i=match(genes,gene.meta$submitted_item)
  gene.meta2=gene.meta[i,]
  if(sum(genes==gene.meta2$submitted_item)!=length(genes)){
    cat('wrong with ',file,'\n')
  }
  
  rownames(dat1@assays$RNA@counts)=gene.meta2$SYMBOL
  rownames(dat1@assays$RNA@data)=gene.meta2$SYMBOL
  dat1=subset(dat1,sex!='mix' & annotation!='unannotated' & annotation!='artefact' )
  dat1 #14996 features across 84477 samples within 1 assay 
  
  rm.tc=c("germline cell", "female reproductive system",
        'male accessory gland','spermatocyte'); 
  dat2=dat1[,!dat1$annotation %in% rm.tc]
  dat2 #14996 features across 81848 samples within 1 assay 
  
  (file.name=gsub('.rds','_valid.rds',file))
  saveRDS(dat2,file.name)
  cat('file',file,'is done\n')
}

