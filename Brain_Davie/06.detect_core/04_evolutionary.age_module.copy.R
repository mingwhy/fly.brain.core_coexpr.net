
# plot module gene member age 
options(stringsAsFactors = F)
library(igraph)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(RColorBrewer)
library(readxl)
library(circlize)
library(org.Dm.eg.db,verbose=F,quietly=T)

out.dir='cytospace/';cut=12;
gene.module.info=data.table::fread('cytospace/cut12_gene_module.txt');
dim(gene.module.info)
# 1295 genes have a module membership

####################################################################
## read in fly gene age data from 'No Evidence for Phylostratigraphic Bias Impacting Inferences on Patterns of Gene Emergence and Evolution'
# more information can be found here: https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/Phylostratigraphy_fly.gene_age 
div.times=c(25,235,325,521,542,560,573,615,653,1547,2100,3900);
length(div.times)
df2017<-read_excel('../external_data/phylostratigraphy_fly.gene.age/TableS1.xlsx',sheet='Dmel_updated_phylostratigraphy',skip=3)
dim(df2017) #13794     4
table(df2017$phylostrata)

# this data contain protein ids
# check the overlap genes between fly  prot_id in 2017.dataset and fly.prot from org.Dm.eg.db.
columns(org.Dm.eg.db)
fly.prot<-keys(org.Dm.eg.db, keytype='FLYBASEPROT')
length(fly.prot) # 30719 fly proteins in total of fly genome in this database
sum(df2017$prot_id %in% fly.prot) #13624 overlap genes
# get flybase ID for each fly prot_id
df.fly.prot<-AnnotationDbi::select(org.Dm.eg.db,
                                   keys=df2017$prot_id,
                                   keytype="FLYBASEPROT",c("FLYBASECG","FLYBASE",'SYMBOL',"GENENAME"))
dim(df.fly.prot) #13794     5
sum(is.na(df.fly.prot$FLYBASE)) #170, equals 13794-13624
sum(is.na(df.fly.prot$SYMBOL)) #170, equals 13794-13624
gene.age=merge(df.fly.prot,df2017,by.x='FLYBASEPROT',by.y='prot_id')
dim(gene.age) #13794
table(gene.age$phylostrata)
sum(is.na(gene.age$FLYBASE)) #170 genes without FLYBASE id.
sum(table(gene.age$phylostrata)) #13794 genes in df2017 publication
gene.age2=gene.age[!is.na(gene.age$FLYBASE),]
dim(gene.age2) #13624 genes
head(gene.age2)
length(unique(gene.age2$SYMBOL)) #13624 fly symbols

gene.age=gene.age2

#####################################################
## age signature
age.group=list();
for(i in unique(gene.module.info$module)){
  genes=gene.module.info[gene.module.info$module==i,]$gene
  sum(genes %in% gene.age$SYMBOL )
  age.group[[as.character(i)]]=gene.age[gene.age$SYMBOL %in% genes,]
}

gene.module=gene.module.info$module
colors = structure(c('white','#52854C'),names = c('0', '1'))
#col_func=circlize::colorRamp2(sort(unique(gene.module)), brewer.pal(max(gene.module),"Paired"))
col_func=colorRamp2(sort(unique(gene.module)), 
                    col=c( "#E69F00", "#009E73",
                           "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
pick.row.cols=col_func(1:max(gene.module)) #module color for 1,2,...

sapply(age.group,dim)
age.group.df=Reduce(`rbind`,age.group)
age.group.df$module.id= rep(1:length(age.group),sapply(age.group,nrow))
age.group.df$module.name=age.group.df$module.id

df=age.group.df %>% group_by(phylostrata,phylostrata_name,module.id,module.name) %>% summarise(ngene=n())

df$module.id=factor(df$module.id)
df$color=pick.row.cols[as.numeric(as.character(df$module.id))] #module color
df$phylostrata_name=as.character(df$phylostrata_name)
df[df$phylostrata_name=='cellular_organisms',]$phylostrata_name='CellLife'

df[!duplicated(df$phylostrata),]$phylostrata_name
df$phylostrata_name=factor(df$phylostrata_name,levels=
                             df[!duplicated(df$phylostrata),]$phylostrata_name)
df1=df %>% group_by(module.id) %>% mutate(ngene.total=sum(ngene))
df1$prop=df1$ngene/df1$ngene.total

df1$module.id=as.numeric(as.character(df1$module.id))
tmp1=merge(as.data.frame(gene.module.info),df1[,c('module.id','color')],by.x='module',by.y='module.id')
tmp1$rgb=apply(col2rgb(tmp1$color),2,function(i){paste(i,collapse = '-')})


# get gene description for genes
tmp2=tmp1[!duplicated(tmp1$gene),];
dim(tmp2) #178 genes
write.table(tmp2,file=paste0(out.dir,'/cut',cut,'_net_module_gene.txt'),quote=F,sep='\t',row.names = F)

tmp2.out=AnnotationDbi::select(org.Dm.eg.db,keys=tmp2$gene,keytype='SYMBOL',c("FLYBASE","SYMBOL", "GENENAME"))

tmp3=merge(tmp2.out,tmp2[,c(1,2)],by.x='SYMBOL',by.y='gene',all.x=T)
tmp3=tmp3[order(tmp3$module),]
write.table(tmp3,file=paste0(out.dir,'cut',cut,'_net_module_gene_name.txt'),quote=F,sep='\t',row.names = F)

write.table(age.group.df[,c(9,4,3,5,6,8)],file=paste0(out.dir,'cut',cut,'_age.group.df.txt'),quote=F,sep='\t',row.names = F)
print( ggplot(age.group.df,aes(x=phylostrata))+
         geom_bar()+theme_bw(base_size = 16)+
         scale_x_continuous(breaks=1:max(age.group.df$phylostrata))+
         facet_wrap(.~module.name,ncol=1,scales="free_y") )



pdf(paste0(out.dir,'cut',cut,'_MCL_modules_age.pdf'),
    height = 9,width = 14)

print( ggplot(df1,aes(x=phylostrata_name,y=module.name,size=prop))+
         geom_point(col=df1$color)+theme_bw(base_size = 16)+
         xlab('Gene age (old<->young)')+ylab('Module')+
         scale_y_continuous(breaks=unique(df$module.name),labels=unique(df$module.name))+
         #scale_color_manual(values=rev(c('#68369b','#2e6fba','#4ead5b','#b89130')))+
         #scale_color_manual(values=c('#68369b','#2e6fba','#4ead5b','#b89130'))+
         #scale_x_continuous(breaks=1:max(age.group.df$phylostrata))+
         scale_size(name='Percentage',range = c(4,12))+
         theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(), 
               axis.text.x=element_text(size=26,angle=45,hjust=1),
               axis.text.y=element_text(size=26),
               axis.title=element_text(size=26,face='bold'),
               legend.key.size = unit(0.5, 'cm'),
               legend.key.height = unit(0.5, 'cm'),
               legend.key.width = unit(0.5, 'cm'), 
               #legend.position = 'none',
               legend.title=element_text(size=22,face="bold"),
               legend.text=element_text(size=20),
               plot.margin=unit(c(1.5, 2, 2, 5.5),"cm")) +
         guides(colour = FALSE ) );

dev.off();



