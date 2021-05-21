library(SummarizedExperiment)
library(TCGAbiolinks)
library(data.table)
####
#Download HTSeq-Count expression data of BRCA patients from TCGA database
query <- GDCquery(project = "TCGA-BRCA",
                  experimental.strategy ="RNA-Seq",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query)
write.csv(query[[1]][[1]],file = "query.csv")
dataBRCA<-GDCprepare(query)
dataexpt<-assay(dataBRCA)
pData<-SummarizedExperiment::colData(dataBRCA)
save(dataexpt,file = "dataexpt.RData")

library("clusterProfiler")
library("org.Hs.eg.db")
dataexpt<-as.data.frame(dataexpt)
geneId <- bitr(rownames(dataexpt), fromType= "ENSEMBL", toType= c("SYMBOL", "GENENAME"), OrgDb="org.Hs.eg.db")
table(duplicated(geneId$SYMBOL))

dataexpt.2 <- dataexpt[geneId$ENSEMBL,]
dim(dataexpt.2)#25856 1222
dataexpt.2 <- cbind(geneId$SYMBOL, dataexpt.2)
dim(dataexpt.2)#25856 1223
write.table(dataexpt.2, file = "C:/Users/jenny/Documents/TCGA/BRCAdata.csv",sep=",",row.names=F, quote = FALSE, na = "NA")

library(readr)
BRCA_TNBC <- read.csv("C:/Users/jenny/Documents/TCGA/BRCA_TNBC.csv",header = TRUE,sep = ",")
row.names(BRCA_TNBC)<-BRCA_TNBC[,1]
BRCA_TNBC<-BRCA_TNBC[,-1]
View(BRCA_TNBC)

rmDupID <- function(matrix){
  expSet <- BRCA_TNBC[,-1]
  rowMean <- apply(expSet, 1, function(x){mean(as.numeric(x), na.rm = TRUE)})
  BRCA_TNBC <- BRCA_TNBC[order(rowMean,decreasing = TRUE),]
  BRCA_TNBC<- BRCA_TNBC[!duplicated(BRCA_TNBC[,1]),]
  rownames(BRCA_TNBC) <- BRCA_TNBC[,1]
  BRCA_TNBC <- BRCA_TNBC[,-1]
  return(BRCA_TNBC)
}
dataexpt.3 <- rmDupID(BRCA_TNBC)
dim(dataexpt.3)