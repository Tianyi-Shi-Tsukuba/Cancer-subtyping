dataExpr=read.csv("Cyto_draw.csv.edges.txt",sep="\t",header = T)
hist(dataExpr[,3])
x=dataExpr[,3]
x[x >=0.05] <- 1
x[x <0.05] <- 0
dataExpr[,3]=x

#write.csv(dataExpr,"change_differential.txt")
write.table (dataExpr, file ="change_differential.txt", sep ="\t", row.names =FALSE, col.names =TRUE, quote =TRUE)
