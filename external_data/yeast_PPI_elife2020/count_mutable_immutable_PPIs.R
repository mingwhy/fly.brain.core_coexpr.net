options(stringsAsFactors = F)
df=read.csv('TableS5_validation_rates_PPI_in_each_environment.csv')
dim(df)
head(df)
x=table(df$Positive_environment_number)
barplot(x)
x[1] #7724
sum(x[-1]) #5257
x[1]/sum(x) #0.5950235 
#We partitioned PPIs by the number of environments in which they were identified and defined PPIs at opposite ends of this spectrum as ‘mutable’ PPIs (identified in only 1–3 environments) and ‘immutable’ PPIs (identified in 8–9 environments). Mutable PPIs far outnumbered immutable PPIs, with PPIs identified in only one environment outnumbering all other PPIs combined (Figure 3A and Figure 3—figure supplement 1). 

df=read.csv('TableS6_high_confidence_PPI_in_each_environment.csv');
dim(df) #5781
head(df)
x=table(df$Positive_environment_number)
barplot(x)
x[1] #7724
sum(x[-1]) #5257
x[1]/sum(x)
