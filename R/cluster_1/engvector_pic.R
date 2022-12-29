cluster=read.csv("C:/Users/10524/Desktop/multiview-master/data_label/BREAST_mtView_label.csv",row.names = 1)
dataset=read.csv("breast_gene.csv",row.names = 1)
dataset=t(dataset)
dim(dataset)
cluster=t(cluster)
dim(cluster)
clu_data=rbind(cluster,dataset)
clu_data=t(clu_data)


Total1<-clu_data[order(clu_data[,1],decreasing = TRUE),]

newdata=Total1[Total1[,1]==5,]
newdata=newdata[,-1]

write.csv(newdata,"breast_cluster5.csv")



######all_cluster######

dim(TOM)
sum(TOM>0.1)

TOM[TOM >0.1] <- 1
TOM[TOM <=0.1] <- 0

B=TOM

tolsum=sum(B==1)
for (i in 1:dim(B)[1]) {
  for (j in 1:dim(B)[2]) {
    B[i,j] = B[i,j] - sum(B[i,])*sum(B[,j])/tolsum
  }
}

B=B+t(B)
ev=eigen(B)
View(ev$vectors)
max(B)
min(B)


L1=colSums(abs(ev$vectors))

max(L1)
min(L1)

#################

ev$vectors[3,]


plot(x = ev$vectors[,24],y = ev$vectors[,94],
     type = 'p',
     xlab = "eigenvector 24",
     ylab = "eigenvector 94",             
     main = "Cluster 1",
     col = "blue", font = 1,cex = 1,
     pch = 19
)



library(ggplot2)
qplot(x = ev$vectors[1,],y = L1)



plot(L1[1:100],type="l",main="L1 norms in cluster1",xlab="Eigenvector index",ylab="L1 norm",col=c("red","blue","pink"))




write.csv(ev$vectors[,94]>(0.05),"ev30.csv")



#######cluster_5############
#https://blog.csdn.net/weixin_43779957/article/details/118864955
data <-runif(10,40,90)
plot(data,type="p",main="绘点图",xlab="x轴",ylab="y轴",col=c("red","blue","pink"))
plot(data,type="l",main="绘点线",xlab="x轴",ylab="y轴",col=c("red","blue","pink"))
plot(data,type="b",main="同时绘制点和线",xlab="x轴",ylab="y轴",col=c("red","blue","pink"))
plot(data,type="c",main="仅绘制参数b所示的线",xlab="x轴",ylab="y轴",col=c("red","blue","pink"))
plot(data,type="o",main="同时绘制点和线，且线穿过点",xlab="x轴",ylab="y轴",col=c("red","blue","pink"))
plot(data,type="h",main="绘制出点到横坐标轴的垂直线",xlab="x轴",ylab="y轴",col=c("red","blue","pink"))
plot(data,type="s",main="绘制出阶梯图（先横后纵）",xlab="x轴",ylab="y轴",col=c("red","blue","pink"))
plot(data,type="S",main="绘制出阶梯图（先纵后竖）",xlab="x轴",ylab="y轴",col=c("red","blue","pink"))
plot(data,type="n",main="空图",xlab="x轴",ylab="y轴")




