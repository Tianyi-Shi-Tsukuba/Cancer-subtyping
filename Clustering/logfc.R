library(edgeR)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")


#计算logFC需要原始的数据 不能有0

library(edgeR)
library("fdrtool")

#200个marker gene太少
#x <- read.csv("Cyto_draw.csv",row.names = 1)

x <- read.csv("D:/mutiview_GLIO_3/gene_sim_Expression.csv",row.names = 1)
x[x<0]=0

x=t(x)

group <-  read.csv("D:/mutiview_GLIO_3/GLIO_SC_cSNN_label_3.csv",row.names = 1)

group=t(group)

dim(x)
dim(group)


y<- DGEList(counts=x,group=group)




#过滤
keep <- rowSums(cpm(y)>1) >= 1
y <- y[keep,,keep.lib.sizes=FALSE]
##TMM 标准化
y <- calcNormFactors(y)
y$samples




###推测离散度，若样本是人，设置bcv = 0.4，模式生物设置0.1(此处是根据有相关教程进行设置，也可以你根据你的结果自行设置，最终得到你的理想值)

bcv <- 0.2
et <- exactTest(y, dispersion=bcv^2)
topTags(et)
summary(de <- decideTestsDGE(et))
###图形展示检验结果
png('0h_vs_2h_MAplot.png')
detags <- rownames(y)[as.logical(de)];
plotSmear(et, de.tags=detags)

#矫正P值
DE <- et$table
res <- DE
head(DE)
res$FDRP <- p.adjust(res$PValue,method = "fdr",n=length(res$PValue))
#res$logP <- -log10(res$FDRP)
head(res)
write.csv(res,"logfc_res.csv")
















