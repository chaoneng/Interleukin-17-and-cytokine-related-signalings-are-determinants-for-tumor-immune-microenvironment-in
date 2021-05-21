cibersoftHL<-as.data.frame(cibersoftHL)
rownames(cibersoftHL)=paste(cibersoftHL[,1]) 
cibersoftHL<-cibersoftHL[,-1]

library(survival)
library(survminer)
library(DESeq2)
coldata <- data.frame(condition = factor(c(rep("StromalhighImmunelow",10),rep("Normal",112)), levels = c('StromalhighImmunelow', 'Normal')))
dds <- DESeqDataSetFromMatrix(countData = cibersoftHL, colData = coldata, design= ~condition)
#
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
#
res <- results(dds1, contrast = c('condition', 'StromalhighImmunelow', 'Normal'))
res
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
##筛选差异表达基因
res1 <- res1[order(res1$padj, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
#log2FC≥1 & padj<0.01 标识 up，代表显著上调的基因
#log2FC≤-1 & padj<0.01 标识 down，代表显著下调的基因
#其余标识 none，代表非差异的基因
res1[which(res1$log2FoldChange >= 4 & res1$padj < 0.0001),'sig'] <- 'up'
res1[which(res1$log2FoldChange <= -4 & res1$padj < 0.0001),'sig'] <- 'down'
res1[which(abs(res1$log2FoldChange) <= 4 | res1$padj >= 0.0001),'sig'] <- 'none'
res1_select <- subset(res1, sig %in% c('up', 'down'))
write.csv(res1_select,file = "res1_select.csv")
#根据 up 和 down 分开输出
res1_up <- subset(res1, sig == 'up')
write.csv(res1_up,file="res1_up.csv")
res1_down <- subset(res1, sig == 'down')
write.csv(res1_down,file="res1_down.csv")
library(ggplot2)

p <- ggplot(data = res1, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(size = 1) +  #绘制散点图
  scale_color_manual(values = c('red', 'gray', 'blue'), limits = c('up', 'none', 'down')) +  #自定义点的颜色
  labs(x = 'log2 Fold Change', y = '-log10 adjust p-value', title = 'StromalhighImmunelow vs Normal', color = '') +  #坐标轴标题
  theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(), #背景色、网格线、图例等主题修改
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-1, 1), lty = 3, color = 'black') +  #添加阈值线
  geom_hline(yintercept = 2, lty = 3, color = 'black') +
  xlim(-6, 10) + ylim(0, 11)  #定义刻度边界
p
#
group_list=c(rep("StromalhighImmunelow",10),rep("Normal",112))
library(DESeq2)
rld<-vst(dds1,blind = FALSE)
expr_vst<-assay(rld)
expr_vst<-as.data.frame(expr_vst)
anno=data.frame(condition=group_list)
rownames(anno)=colnames(expr_vst)
anno1=data.frame(Expression=res1_select$sig)
rownames(anno1)=rownames(DEG_expr_matr)
library(pheatmap)
res1_select<-as.data.frame(res1_select)
DEG_expr_matr<-expr_vst[rownames(res1_select),]
z_score_matrix<-t(scale(t(DEG_expr_matr)))
#write.csv(z_score_matrix,file = "z_score_matrix.csv")
pheatmap(mat = z_score_matrix,cluster_rows = T,show_colnames = F,cluster_cols = T,
         annotation_col = anno,
         annotation_row = anno1,
         show_rownames = F,color = colorRampPalette(c("navy","white","firebrick3"))(100),
         legend = T,legend_breaks = c(-2,0,2),breaks = unique(seq(-2,2,length=100)),
)
##GO分析
library(clusterProfiler)#富集分析
library(topGO)#画GO图用的
library(Rgraphviz)
library(pathview)#看KEGG pathway的
library(org.Hs.eg.db)
library(ggplot2)
###########
genesymbol<-as.character(overexp$Gene_id)
entrezid<-mapIds(x=org.Hs.eg.db,
                 keys=genesymbol,
                 keytype = "SYMBOL",
                 column = "ENTREZID")
entrezid<-as.data.frame(entrezid)
#
GO<-enrichGO(
  gene=entrezid$entrezid,
  OrgDb = 'org.Hs.eg.db',
  keyType = "ENTREZID",
  ont = "ALL",
  pvalueCutoff = 0.01,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = T )
head(GO)
write.csv(GO,file = "GO.csv")
dotplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")
#
KEGG<-enrichKEGG(
  gene = entrezid$entrezid,
  organism = "hsa",#人的组织数据
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
  )
barplot(KEGG,showCategory = 10,title = 'KEGG Pathway')
head(KEGG)
write.csv(KEGG,file = "KEGG.csv")
##########
genesymbol2<-as.character(underexp$Gene_id)
entrezid2<-mapIds(x=org.Hs.eg.db,
                 keys=genesymbol2,
                 keytype = "SYMBOL",
                 column = "ENTREZID")
entrezid2<-as.data.frame(entrezid2)
#
GO2<-enrichGO(
  gene=entrezid2$entrezid2,
  OrgDb = 'org.Hs.eg.db',
  keyType = "ENTREZID",
  ont = "ALL",
  pvalueCutoff = 0.01,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = T )
head(GO2)
write.csv(GO2,file = "GO2.csv")
dotplot(GO2, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")
#
KEGG2<-enrichKEGG(
  gene = entrezid2$entrezid2,
  organism = "hsa",#人的组织数据
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
barplot(KEGG2,showCategory = 10,title = 'KEGG Pathway')
head(KEGG2)
write.csv(KEGG2,file = "KEGG2.csv")
##370
gene370<-as.character(row.names(res1_select))
entrezid370<-mapIds(x=org.Hs.eg.db,
                 keys=gene370,
                 keytype = "SYMBOL",
                 column = "ENTREZID")
entrezid370<-as.data.frame(entrezid370)
GO370<-enrichGO(
  gene=entrezid370$entrezid370,
  OrgDb = 'org.Hs.eg.db',
  keyType = "ENTREZID",
  ont = "ALL",
  pvalueCutoff = 0.01,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = T )
head(GO370)
write.csv(GO370,file = "GO370.csv")
dotplot(GO370, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")
KEGG370<-enrichKEGG(
  gene = entrezid370$entrezid370,
  organism = "hsa",#人的组织数据
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
barplot(KEGG370,showCategory = 10,title = 'KEGG Pathway')
head(KEGG370)
write.csv(KEGG370,file = "KEGG370.csv")
##187
gene187<-as.character(row.names(res1_up))
entrezid187<-mapIds(x=org.Hs.eg.db,
                    keys=gene187,
                    keytype = "SYMBOL",
                    column = "ENTREZID")
entrezid187<-as.data.frame(entrezid187)
GO187<-enrichGO(
  gene=entrezid187$entrezid187,
  OrgDb = 'org.Hs.eg.db',
  keyType = "ENTREZID",
  ont = "ALL",
  pvalueCutoff = 0.01,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = T )
head(GO187)
write.csv(GO187,file = "GO187.csv")
dotplot(GO187, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")
KEGG187<-enrichKEGG(
  gene = entrezid187$entrezid187,
  organism = "hsa",#人的组织数据
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
barplot(KEGG187,showCategory = 10,title = 'KEGG Pathway')
head(KEGG187)
write.csv(KEGG187,file = "KEGG187.csv")
##183
gene183<-as.character(row.names(res1_down))
entrezid183<-mapIds(x=org.Hs.eg.db,
                    keys=gene183,
                    keytype = "SYMBOL",
                    column = "ENTREZID")
entrezid183<-as.data.frame(entrezid183)
GO183<-enrichGO(
  gene=entrezid183$entrezid183,
  OrgDb = 'org.Hs.eg.db',
  keyType = "ENTREZID",
  ont = "ALL",
  pvalueCutoff = 0.01,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = T )
head(GO183)
write.csv(GO183,file = "GO183.csv")
dotplot(GO183, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")
KEGG183<-enrichKEGG(
  gene = entrezid183$entrezid183,
  organism = "hsa",#人的组织数据
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
barplot(KEGG183,showCategory = 10,title = 'KEGG Pathway')
head(KEGG183)
write.csv(KEGG183,file = "KEGG183.csv")
