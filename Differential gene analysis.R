setwd("E:/a课题相关/xiaoheirna/csv/csv2/")
source("https://bioconductor.org/biocLite.R") 
biocLite('BiocInstaller') 
biocLite("DESeq2") 
install.packages(c("gplots", "amap", "ggplot2"))
library(DESeq2)
library("RColorBrewer") 
library("gplots")
library("amap")
library("ggplot2")
library("BiocParallel")

#读取基因counts文件
data<-read.csv("g_count_vol1.csv")
head(data)
#删除表头name
rownames(data)<-data[,1]
data<-data[,-1]
head(data)
#撇掉在多于两个样本中 count<1 的值，如果样本数多，这个数值可以适当增加 
#排除极低表达基因的干扰
data <- data[rowSums(data)>2,] 
head(data)

#构建样品信息矩阵
conditions <- factor(c(rep("HV",3),rep("LV",3)), levels = c("HV","LV"))
conditions
sample <- data.frame(row.names=colnames(data), conditions)
sample

#按照DESeq格式产生DESeq数据集
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = data, colData = sample, design= ~ conditions)
dds <- DESeq(ddsFullCountTable)

# ?counts 查看此函数功能 
# normalized=T, 返回标准化的数据 
normalized_counts <- counts(dds, normalized=TRUE) 
head(normalized_counts)

#根据基因在不同的样本中表达变化的差异程度 mad 值对数据排序，差异越大的基因排位越前。 
normalized_counts_mad <- apply(normalized_counts, 1, mad) 
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]

# 标准化后的数据输出 
write.table(normalized_counts, file="Count_matrix.xls.DESeq2.normalized.xls", quote=F, sep="\t", row.names=T, col.names=T)

# log 转换后的结果 
rld <- rlog(dds, blind=FALSE) 
rlogMat <- assay(rld) 
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]
write.table(rlogMat, file="Count_matrix.xls.DESeq2.normalized.rlog.xls", quote=F, sep="\t", row.names=T, col.names=T)

#样品层级聚类分析，判断样品的相似性和组间组内差异
# 生成颜色 
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

# 计算相关性pearson correlation  
pearson_cor <- as.matrix(cor(rlogMat, method="pearson"))

# 层级聚类 
hc <- hcluster(t(rlogMat), method="pearson")

# 热图绘制 
heatmap.2(pearson_cor, Rowv=as.dendrogram(hc), symm=T, trace="none", col=hmcol, margins=c(11,11), main="The pearson correlation of each sample") 

#样品PCA分析 
plotPCA(rld, intgroup=c("conditions"))

#差异基因分析
sampleA = "HV" 
sampleB = "LV"
contrastV <- c("conditions", sampleA, sampleB) 
res <- results(dds, contrast=contrastV) 
res

#给DESeq2的原始输出结果增加样品平均表达信息，使得结果更容易理解和解析。
# 获得第一组数据均值 
baseA <- counts(dds, normalized=TRUE)[, colData(dds)$conditions == sampleA] 
if (is.vector(baseA)){ 
  baseMeanA <- as.data.frame(baseA) 
  } else { 
    baseMeanA <- as.data.frame(rowMeans(baseA)) 
    } 
colnames(baseMeanA) <- sampleA 
head(baseMeanA)

# 获得第二组数据均值 
baseB <- counts(dds, normalized=TRUE)[, colData(dds)$conditions == sampleB] 
if (is.vector(baseB)){ 
  baseMeanB <- as.data.frame(baseB) 
  } else {
    baseMeanB <- as.data.frame(rowMeans(baseB)) 
    } 
colnames(baseMeanB) <- sampleB 
head(baseMeanB)

# 结果组合 
res <- cbind(baseMeanA, baseMeanB, as.data.frame(res)) 
head(res)

# 增加 ID 信息 
res <- cbind(ID=rownames(res), as.data.frame(res)) 
res$baseMean <- rowMeans(cbind(baseA, baseB))
# pvalue为 NA 赋值为 1 
res$pvalue[is.na(res$pvalue)] <- 1
# 按 pvalue 排序, 把差异大的基因放前面 
res <- res[order(res$pvalue),] 
head(res)

#整体分析结果输出到文件 
comp314 <- paste(sampleA, "_vs_", sampleB, sep=".")
# 生成文件名 
file_base <- paste("Count_matrix.xls.DESeq2", comp314, sep=".") 
file_base1 <- paste(file_base, "results.xls", sep=".") 
write.table(as.data.frame(res), file=file_base1, sep="\t", quote=F, row.names=F)

#提取差异表达基因
# 差异基因筛选，pvalue<0.1
res_de <- subset(res, res$pvalue<0.1, select=c('ID', sampleA, 
                                             sampleB, 'log2FoldChange', 'pvalue')) 
# foldchang > 1
res_de_up <- subset(res_de, res_de$log2FoldChange>=1) 
file <- paste("Count_matrix.xls.DESeq2",sampleA,"_higherThan_",sampleB, 'xls', sep=".") 
write.table(as.data.frame(res_de_up), file=file, sep="\t", quote=F, row.names=F) 

res_de_dw <- subset(res_de, res_de$log2FoldChange<=(-1)*1) 
file <- paste("Count_matrix.xls.DESeq2",sampleA, "_lowerThan_", sampleB, 'xls', sep=".") 
write.table(as.data.frame(res_de_dw), file=file, sep="\t", quote=F, row.names=F)

# 差异基因 ID 
res_de_up_id = data.frame(ID=res_de_up$ID, 
                          type=paste(sampleA,"_higherThan_", sampleB, sep=".")) 
res_de_dw_id = data.frame(ID=res_de_dw$ID, type=paste(sampleA,"_higherThan_", sampleB, sep=".")) 
de_id = rbind(res_de_up_id, res_de_dw_id) 
file <- "Count_matrix.xls.DESeq2.all.DE.xls" 
write.table(as.data.frame(de_id), file=file, sep="\t", quote=F, row.names=F, col.names=F)

#绘制火山图
logCounts <- log2(res$baseMean+1) 
logFC <- res$log2FoldChange 
FDR <- res$pvalue
plot(logFC, -1*log10(FDR), col=ifelse(FDR<=0.01, "red", "black"), xlab="logFC", ylab="-1*log1o(FDR)", main="Volcano plot", pch=".") 

#差异基因热图 
res_de_up_top20_id <- as.vector(head(res_de_up$ID,20)) 
res_de_dw_top20_id <- as.vector(head(res_de_dw$ID,20)) 
red_de_top20 <- c(res_de_up_top20_id, res_de_dw_top20_id) 
red_de_top20
red_de_top20_expr <- normalized_counts[rownames(normalized_counts) %in% red_de_top20,] 
install.packages("pheatmap")
library(pheatmap)
pheatmap(red_de_top20_expr, cluster_row=T, scale="row", annotation_col=sample)

#GO富集分析
source("https://bioconductor.org/biocLite.R") 
biocLite(c("clusterProfiler", "org.Bt.eg.db")) 
library("org.Bt.eg.db") 
readable=TRUE 
library(clusterProfiler)

#转换为ENTREZ ID
k=keys(org.Bt.eg.db,keytype = "ENSEMBL")
head(k,5)
#基于提取的ENSEMBL ID，提取对应的所有Gene ID(ENTREZID)，(以及Symbol），并查看一下提取的内容。
list=select(org.Bt.eg.db,keys=k,columns = c("ENTREZID","SYMBOL"), keytype="ENSEMBL")
dim(list)
head(list,5)
ID<-de_id$ID
ID_list=list[match(ID,list[,"ENSEMBL"]),]
ID_list
head(ID_list)
ENTREZ_ID<-na.omit(ID_list$ENTREZID)

ego <- enrichGO(gene = ENTREZ_ID, OrgDb = org.Bt.eg.db,ont = "ALL",pvalueCutoff = 0.5,qvalueCutoff = 1,readable = TRUE)
head(ego)
dotplot(ego,showCategory=20,title="EnrichmentGO_dot")
barplot(ego, showCategory=20,title="EnrichmentGO")

kk <- enrichKEGG(gene = ENTREZ_ID,
                 organism = 'bta', #KEGG可以用organism = '***'
                 pvalueCutoff = 1)

