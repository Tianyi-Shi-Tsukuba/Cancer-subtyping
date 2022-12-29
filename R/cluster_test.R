fourier  <- read.table("D:/Data_Multi-omics/GBM/GLIO_Gene_Expression.txt",header=TRUE,row.names = 1)
profcorr <- read.table("D:/Data_Multi-omics/GBM/GLIO_Methy_Expression.txt", header=TRUE,row.names = 1)
pixels   <- read.table("D:/Data_Multi-omics/GBM/GLIO_Mirna_Expression.txt", header=TRUE,row.names = 1)

fourier=t(fourier)
profcorr=t(profcorr)
pixels=t(pixels)

# length(classes)
classes <- as.vector(sapply(0:9, function(x) rep(x,200)))

clust   <- mvsc(list(fourier,profcorr,pixels), k=2)

clust$clustering




write.table(clust$clustering, file = "GLIO_SC_cSNN_label_2.txt",row.names=F,col.names = F)

write.csv(clust$clustering, file = "GLIO_SC_cSNN_label_2.csv")


# $clustering member has the clustering assignment vector
knitr::kable(table(classes, clust$clustering))

## ---- fig.width=7, fig.height=7------------------------------------------
clust   <- mvsc(list(fourier, profcorr, pixels, morpho), k=10, neighbours=2)

# $clustering member has the clustering assignment vector
knitr::kable(table(classes, clust$clustering))
clust$clustering


write.csv(clust$clustering, file = "LUNG_mtView_label.csv")
clust$clustering

write.csv(fourier,"sim_gene_Expression.csv")











