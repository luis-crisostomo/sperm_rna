
#Help
help(package="DESeq2",help="html")
vignette("DESeq2")

#libraries
library("DESeq2")
library("vsn")
library("pheatmap")
library("ggplot2")
library(ggrepel)
library(tidyverse)
library("genefilter")
library("sqldf")
library(RColorBrewer)
library(data.table)
library(ggpubr)

#Datasets preparation

cts <- read.table("merge_tiRNAs.txt", header = T, sep = "", row.names = "Annotation")
cts <- cts[, c(2:25)]
cts[is.na(cts)] = 0

coldata <- read.csv("sample_annotation.csv", header = T, row.names = 1)

all(rownames(coldata) == colnames(cts))
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata,design = ~ Condition)
dds
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 10, ]
nrow(dds)

#Data vizualization
lambda <- 10^seq(from = -1, to = 2, length = 1000)
cts <- matrix(rpois(1000*100, lambda), ncol = 100)
meanSdPlot(cts, ranks = FALSE)
log.cts.one <- log2(cts + 1)
meanSdPlot(log.cts.one, ranks = FALSE)

#Data transformation
vsd <- varianceStabilizingTransformation(dds)         

#Heatmap plot 1000 * 600
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$name, vsd$Condition, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_method="ward.D2",
         clustering_distance_cols = sampleDists,
         col = colors, fontsize = 15)

#PCA plot 1500 * 900
pcaData <- plotPCA(vsd, intgroup = c( "name", "Condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x= PC1, y= PC2, color=Condition, shape=Condition))+
  scale_shape_manual(values=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11))+
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 35) +
  geom_point(size = 5)
geom_text_repel(aes(label=Condition), size=6, vjust=0.6, hjust=0.5)

#Making sure to set up which control samples to use
dds$Condition <- relevel(dds$Condition, ref = "HFD_F0")
as.data.frame(colData(dds))
dds <- DESeq(dds)
res <- results(dds)
res
resultsNames(dds)
summary(res)

#Running the analysis
res <- results( dds, contrast = c("Condition", "HFDT_F0", "HFD_F0") )
res
resOrdered <- res[order(res$padj),]
head(resOrdered)
mcols(res, use.names = TRUE)
summary(res)
sum( res$padj < 0.1, na.rm = TRUE)
table( is.na(res$pvalue))
sum( res$padj < 0.1, na.rm = TRUE)
resSig <- res[ which(res$padj < 0.1 & abs(res$log2FoldChange) > 1.5), ]
head( resSig)
tail( resSig )
write.csv( as.data.frame(res), file="results_tiRNAs_HFDT_F0vsHFD_F0.csv" )

df <- read.csv("results_tiRNAs_HFDT_F0vsHFD_F0.csv", header = T)
setnames(df, "X", "Gene_id")
write.csv( as.data.frame(df), file="results_tiRNAs_HFDT_F0vsHFD_F0.csv", row.names = FALSE )

#MA plot of the results 1000 * 600
ggmaplot(res, fdr = 0.1, fc = 2.828427, size = 1,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         legend="top", top = 15, font.label = c("bold", 10), 
         main = "HFDT_F0 vs HFD_F0 tiRNAs",
         label.rectangle = TRUE, font.legend = c("bold",12),font.main = "bold")

