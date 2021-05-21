cibersoftLH<-as.data.frame(cibersoftLH)
rownames(cibersoftLH)=paste(cibersoftLH[,1]) 
cibersoftLH<-cibersoftLH[,-1]
#
library(survival)
library(survminer)
library(DESeq2)
coldata <- data.frame(condition = factor(c(rep("StromallowImmunehigh",10),rep("Normal",112)), levels = c('StromallowImmunehigh', 'Normal')))
dds <- DESeqDataSetFromMatrix(countData = cibersoftLH, colData = coldata, design= ~condition)
#
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
#
res <- results(dds1, contrast = c('condition', 'StromallowImmunehigh', 'Normal'))
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
#默认情况下，横轴展示 log2FoldChange，纵轴展示 -log10 转化后的 padj
p <- ggplot(data = res1, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(size = 1) +  #绘制散点图
  scale_color_manual(values = c('red', 'gray', 'blue'), limits = c('up', 'none', 'down')) +  #自定义点的颜色
  labs(x = 'log2 Fold Change', y = '-log10 adjust p-value', title = 'StromallowImmunehigh vs Normal', color = '') +  #坐标轴标题
  theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(), #背景色、网格线、图例等主题修改
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-1, 1), lty = 3, color = 'black') +  #添加阈值线
  geom_hline(yintercept = 2, lty = 3, color = 'black') +
  xlim(-6, 10) + ylim(0, 11)  #定义刻度边界
p
#
group_list=c(rep("StromallowImmunehigh",10),rep("Normal",112))
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
##651
library(clusterProfiler)#富集分析
library(topGO)#画GO图用的
library(Rgraphviz)
library(pathview)#看KEGG pathway的
library(org.Hs.eg.db)
library(ggplot2)
gene651<-as.character(row.names(res1_select))
entrezid651<-mapIds(x=org.Hs.eg.db,
                    keys=gene651,
                    keytype = "SYMBOL",
                    column = "ENTREZID")
entrezid651<-as.data.frame(entrezid651)
GO651<-enrichGO(
  gene=entrezid651$entrezid651,
  OrgDb = 'org.Hs.eg.db',
  keyType = "ENTREZID",
  ont = "ALL",
  pvalueCutoff = 0.01,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = T )
head(GO651)
write.csv(GO651,file = "GO651.csv")
dotplot(GO651, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")
KEGG651<-enrichKEGG(
  gene = entrezid651$entrezid651,
  organism = "hsa",#人的组织数据
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
barplot(KEGG651,showCategory = 10,title = 'KEGG Pathway')
head(KEGG651)
write.csv(KEGG651,file = "KEGG651.csv")
##284
gene284<-as.character(row.names(res1_up))
entrezid284<-mapIds(x=org.Hs.eg.db,
                    keys=gene284,
                    keytype = "SYMBOL",
                    column = "ENTREZID")
entrezid284<-as.data.frame(entrezid284)
GO284<-enrichGO(
  gene=entrezid284$entrezid284,
  OrgDb = 'org.Hs.eg.db',
  keyType = "ENTREZID",
  ont = "ALL",
  pvalueCutoff = 0.01,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = T )
head(GO284)
write.csv(GO284,file = "GO284.csv")
dotplot(GO284, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")
KEGG284<-enrichKEGG(
  gene = entrezid284$entrezid284,
  organism = "hsa",#人的组织数据
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
barplot(KEGG284,showCategory = 10,title = 'KEGG Pathway')
head(KEGG284)
write.csv(KEGG284,file = "KEGG284.csv")
##367
gene367<-as.character(row.names(res1_down))
entrezid367<-mapIds(x=org.Hs.eg.db,
                    keys=gene367,
                    keytype = "SYMBOL",
                    column = "ENTREZID")
entrezid367<-as.data.frame(entrezid367)
GO367<-enrichGO(
  gene=entrezid367$entrezid367,
  OrgDb = 'org.Hs.eg.db',
  keyType = "ENTREZID",
  ont = "ALL",
  pvalueCutoff = 0.01,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = T )
head(GO367)
write.csv(GO367,file = "GO367.csv")
dotplot(GO367, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")
KEGG367<-enrichKEGG(
  gene = entrezid367$entrezid367,
  organism = "hsa",#人的组织数据
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
barplot(KEGG367,showCategory = 10,title = 'KEGG Pathway')
head(KEGG367)
write.csv(KEGG367,file = "KEGG367.csv")
