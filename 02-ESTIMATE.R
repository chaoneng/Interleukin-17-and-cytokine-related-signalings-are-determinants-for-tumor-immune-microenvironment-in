#Estimation of immune cell content in tumour tissue using RNA-seq data
##
library(estimate)
cancer_type = "TNBC"
input.f=paste0(cancer_type,'_estimate_input.txt')
output.f=paste0(cancer_type,'_estimate_gene.gct')
output.ds=paste0(cancer_type,'_estimate_score.gct')

write.table(dataexpt.3,file = input.f,sep = '\t',quote = F)

filterCommonGenes(input.f = input.f,
                  output.f = output.f ,
                  id = "GeneSymbol")

estimateScore(input.ds = output.f,
              output.ds=output.ds)
# "1 gene set: StromalSignature  overlap= 136"
# "2 gene set: ImmuneSignature  overlap= 140"
scores = read.table(output.ds,skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
scores<-as.data.frame(scores)
View(scores)
write.csv(scores,file = "scores.csv")

plotPurity(scores = "TNBC_estimate_score.gct")

整理資料
tmp<-substring(row.names(scores),1,12)
tmp<-as.data.frame(tmp)
bcr_patient_barcode<-tmp
tmp2<- cbind(tmp, scores)
row.names(tmp2)<-tmp2[,1]
tmp3<-merge(
  nTNBC,tmp2,by="bcr_patient_barcode")
write.csv(tmp3,file="tmp3.csv")