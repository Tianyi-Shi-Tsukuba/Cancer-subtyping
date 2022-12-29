dataExpr=read.csv("GLIO_intergene_cluster3.csv.edges.txt",sep="\t",header = T)
hist(dataExpr[,3])
x=dataExpr[,3]
x[x >=(-0.2)] <- 1
x[x <(-0.2)] <- 0
dataExpr[,3]=x

#write.csv(dataExpr,"change_differential.txt")
write.table (dataExpr, file ="change_differential_cluster3.txt", sep ="\t", row.names =FALSE, col.names =TRUE, quote =TRUE)

