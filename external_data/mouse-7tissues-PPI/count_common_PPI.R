# I download the excel file from https://www.sciencedirect.com/science/article/pii/S0092867421007042
# Download : Download spreadsheet (4MB)
# Table S2. PCP-SILAM interactomes of seven mouse tissues, related to Figure 2.

library(readxl)
df=read_excel('./1-s2.0-S0092867421007042-mmc2.xlsx',sheet='Table S3')
head(df)
dim(df)
table(df$Tissue) #number of PPI per tissue
# Brain  Heart Kidney  Liver   Lung Muscle Thymus 
# 22673  35536  34847  19804  26138  28561  24901 
apply(df[1:2,c(1,2)],1,paste)
df$edge=apply(df[,c(1,2)],1,function(i){paste(i,collapse = '-')})
head(df)
x=table(df$edge)
length(x) #125696 unique PPI
table(x)
#1     2     3     4     5     6     7 
#91100 18689  7457  3930  2218  1313   989
sum(x==7)/length(x) 
#  0.00786819
sum(x==7)/table(df$Tissue) 
#Brain      Heart     Kidney      Liver       Lung     Muscle     Thymus 
#0.04362016 0.02783093 0.02838121 0.04993941 0.03783763 0.03462764 0.03971728 

# quote from original publication:
# The resulting maps provide a pro- teome-scale survey of interactome rewiring across mammalian tissues, 
# revealing more than 125,000 unique interactions at a quality comparable to the highest-quality human screens.

