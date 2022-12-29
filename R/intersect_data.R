exp=read.table("D:/Data_Multi-omics/GBM/GLIO_Gene_Expression.txt",header=TRUE,row.names = 1)
#write.csv(exp, file = "cover_breast.csv")


exp=read.csv("breast_gene.csv")
write.csv(exp[,1],"breast_name.csv")

library(tidyverse)
#dir.create("cover") 
files<-dir(path = "./cover",
           full.names = T,
           pattern = ".csv")

df<-map(files,read.csv)
class(df)
df1<-reduce(df,inner_join)
dim(df1)
View(df1)
write.csv(df1, file = "logFC_res2.csv")

##需要转置后再写出
###
df1=t(df1)
write.csv(df1, file = "differential_breast_data.csv")