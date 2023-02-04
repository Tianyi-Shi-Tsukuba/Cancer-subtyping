breast=read.table("D:/multi_omic_data/gbm/exp",header = T,row.names = 1)
breast=t(breast)
write.csv(breast,"GLIO_gene.csv")

# breast=read.table("D:/Data_Multi-omics/BIC/BREAST_Methy_Expression.txt",header = T,row.names = 1)
# breast=t(breast)
# write.csv(breast,"breast_Methy.csv")
# 
# breast=read.table("D:/Data_Multi-omics/BIC/BREAST_Mirna_Expression.txt",header = T,row.names = 1)
# breast=t(breast)
# write.csv(breast,"breast_Mirna.csv")

install.packages("BiocManager") 
BiocManager::install("WGCNA")
library(WGCNA)
