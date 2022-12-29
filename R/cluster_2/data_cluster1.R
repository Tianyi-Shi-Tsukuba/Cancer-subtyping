Express  <- read.table("D:/Data_Multi-omics/GBM/GLIO_Gene_Expression.txt",header=TRUE,row.names = 1)
label <- read.csv("D:/mutiview_GLIO_change_simitarity/test_pic_documents/GLIO_SC_cSNN_label.csv",row.names = 1)
Express=t(Express)
lab_Ex=cbind(label,Express)

Total1<-lab_Ex[order(lab_Ex[,1],decreasing = TRUE),]
newdata=Total1[Total1[,1]==2,]
newdata=newdata[,-1]
dim(newdata)
write.csv(newdata,"GLIO_cluster2.csv")












