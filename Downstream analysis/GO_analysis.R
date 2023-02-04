install.packages("remotes")
remotes::install_github("WCDZD/KEGG.db", force = TRUE)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")

library(org.Hs.eg.db)
library(openxlsx)
library(clusterProfiler)
library(ggplot2)
library(enrichplot)
library(GOplot)
library(DOSE)
library(stringr)

library(KEGG.db)
fr<-read.csv("logFC_res2.csv")
gene<-bitr(fr$Genes,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = 'org.Hs.eg.db')

GO<-enrichGO(
  gene$ENTREZID,
  OrgDb = 'org.Hs.eg.db',
  keyType = "ENTREZID",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  minGSSize = 10,
  maxGSSize = 500,
  readable = TRUE
) #GO富集
KEGG<-enrichKEGG(
  gene$ENTREZID,
  organism = "hsa",#我用到是数据是人的组织数据，所以这里选择‘hsa’
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 5,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
) #KEGG富集
# 富集基因与所在功能集/通路集的关联网络图：
enrichplot::cnetplot(GO,circular=FALSE,colorEdge = TRUE)#基因-通路关联网络图
enrichplot::cnetplot(KEGG,circular=FALSE,colorEdge = TRUE)#circluar为指定是否环化，基因过多时建议设置为FALSE



#也可以以热图形式展现关联关系:
enrichplot::heatplot(GO,showCategory = 50)#基因-通路关联热图
enrichplot::heatplot(KEGG,showCategory = 50)



#富集到的功能集/通路集之间的关联网络图：
GO2 <- pairwise_termsim(GO)
KEGG2 <- pairwise_termsim(KEGG)
enrichplot::emapplot(GO2,showCategory = 50, color = "p.adjust", layout = "kk")#通路间关联网络图
enrichplot::emapplot(KEGG2,showCategory =50, color = "p.adjust", layout = "kk")



barplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free") #GO聚类条形图，fig1

dotplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free") #GO聚类气泡图,fig2

barplot(KEGG,showCategory = 40,title = 'KEGG Pathway') #KEGG聚类条形图,fig3

dotplot(KEGG) #KEGG聚类气泡图,fig4

enrichplot::cnetplot(GO,circular=TRUE,colorEdge = TRUE)#GO通路-基因网络图，fig5

enrichplot::cnetplot(KEGG,circular=TRUE,colorEdge = TRUE)#KEGG通路-基因网络图,fig6

enrichplot::heatplot(GO,showCategory = 50) #GO富集瀑布图,fig7



MOXD1


enrichplot::heatplot(KEGG,showCategory = 50) #kegg富集瀑布图,fig8

GOplotIn<-GO[1:10,c(2,3,7,9)] #我们先提取GO富集结果的前10行，和提取ID,Description,p.adjust,GeneID四列。
GOplotIn<-GO[1:1000]

GOplotIn$geneID <-str_replace_all(GOplotIn$geneID,'/',',') #把GeneID列中的’/’替换成‘,’
names(GOplotIn)<-c('ID','Term','adj_pval','Genes')#修改列名，后面和弦图绘制的时候需要这样的格式，不然会报错
GOplotIn$Category = "BP"#因为我们提取的前10列为BP,所以再加一列分类信息

genedata<-data.frame(ID=fr$Genes,logFC=fr$logFC)
circ<-GOplot::circle_dat(GOplotIn,genedata) #GOplot导入数据格式整理
chord<-chord_dat(data = circ,genes = genedata) #生成含有选定基因的数据框

GOChord( #GO富集和弦图，fig9
  data = chord,
  title = 'GOchord plot',
  space = 0,#GO Term间距
  limit = c(1,1),
  gene.order = 'logFC',
  gene.space = 0.25,
  gene.size = 5,
  lfc.col = c('red','white','blue'), #上下调基因颜色
  #ribbon.col = brewer.pal(length(GOplotIn$Term)),#GO Term colors
  process.label = 10 #GO Term字体大小
)





GOCircle(circ) #和弦图+表格，fig10 
GOCluster(circ,GOplotIn$Term) #GO富集聚类图，fig11

