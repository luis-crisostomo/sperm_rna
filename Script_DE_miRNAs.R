
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

cts <- read.table("merge_miRNAs.txt", header = T, sep = "", row.names="Annotation")
cts$row.names <- NULL
cts[is.na(cts)] = 0
coldata <- read.csv(file = "sample_annotation.csv", header = T, row.names = 1)
head(cts,2)
coldata
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
write.csv( as.data.frame(res), file="results_miRNAs_HFDT_F0vsHFD_F0.csv" )

df <- read.csv("results_miRNAs_HFDT_F0vsHFD_F0.csv", header = T)
setnames(df, "X", "Gene_id")
write.csv( as.data.frame(df), file="results_miRNAs_HFDT_F0vsHFD_F0.csv", row.names = FALSE )

#MA plot of the results 1000 * 600
ggmaplot(res, fdr = 0.1, fc = 2.828427, size = 1,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         legend="top", top = 10, font.label = c("bold", 10), 
         main = "HFDT_F0 vs HFD_F0 miRNAs",
         label.rectangle = TRUE, font.legend = c("bold",12),font.main = "bold")


resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
#resdatafilt <- filter(resdata, baseMean >= 100)
#resdatafilt <- dplyr::filter(resdatafilt, grepl("p", Gene))
resdata$Gene <- gsub("mature-tRNA-", "", resdata$Gene)
resdata$Gene <- gsub("mature-", "", resdata$Gene)

resdata$Significant = as.factor(abs(resdata$log2FoldChange) > 1.5 & resdata$padj < 0.1)
ggplot(resdata, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(color=Significant), size =6) +
  scale_color_manual(values = c("grey", "green")) +
  theme_bw(base_size = 35) + theme(legend.position = "none") +
  xlim(c(-6, 6)) + ylim(c(0, 3)) + 
  labs(title = "Volcano plot LF...Met vs LF...LF diet") +
  geom_vline(xintercept = 1.5, col = "red", linetype = "dotted", size = 2) +
  geom_vline(xintercept = -1.5, col = "red", linetype = "dotted", size = 2) +
  geom_hline(yintercept = -log10(0.1), col = "red", linetype = "dotted", size = 2) +
  geom_text_repel(
    data = subset(resdata, resdata$padj < 0.1 & abs(resdata$log2FoldChange) > 1.5),
    aes(label = Gene),
    size=6
  )


#create heatmap from sample variances
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 200)
mat  <- assay(rld)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
mat <- as.data.frame(mat, header=T, sep="")
write.csv(as.data.frame(mat), file="test.csv")
matfilt <- read.csv(file = "test.csv", header = T, sep=",")
write.csv( as.data.frame(resSig), file="results_tRNAs_HFHFvsChCh_Pvalue.csv" )
res <- read.csv("results_tRNAs_HFHFvsChCh_Pvalue.csv")
matfilt2 <- sqldf("select * from matfilt f1 inner join res f2 on (f1.X == f2.X) ")
matfilt2$X <- NULL
matfilt2$baseMean <- NULL
matfilt2$log2FoldChange <- NULL
matfilt2$lfcSE <- NULL
matfilt2$stat <- NULL
matfilt2$pvalue <- NULL
matfilt2$padj <- NULL
matfilt3 <- matfilt2[,c(19,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)]
write.csv(as.data.frame(matfilt3), file="test.csv")
matfilt <- read.csv(file = "test.csv", header = T, row.names = "X")
matfilt$X.1 <- NULL

anno <- as.data.frame(colData(rld)[, c("Condition","Type")])
anno$Type <- NULL
levels(anno$Condition)
anno$Condition <- factor(anno$Condition, levels = unique(anno$Condition))
levels(anno$Condition)
pheatmap(matfilt, annotation_col = anno, annotation_names_col = F, show_colnames = FALSE, fontsize = 24, cellwidth = 22, cellheight = 22, cluster_rows = FALSE, main = "HFHFvsChCh")

#For Luis data
dds <- dds[ , dds$time == "48h" ]
dds$time <- droplevels( dds$time )
dds$treatment <- relevel( dds$treatment, "Control" )
as.data.frame( colData(dds) )

#Count outliers
assays (dds) [["cooks"]]
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range = 0, las = 2)
plotDispEsts(dds)

#plot counts

#HFvsCh
plotCounts(dds, gene="mature-mt_tRNA-Glu-TTC_3_end", intgroup= "Condition")

#HFHFvsChCh
par(mfrow=c(6,3), cex = 1.1, cex.main=1.2, cex.lab=1.2)
plotCounts(dds, gene="mature-tRNA-Gln-CTG", intgroup= "Condition")
plotCounts(dds, gene="mature-mt_tRNA-Leu-TAG_5_end", intgroup= "Condition")
plotCounts(dds, gene="mature-mt_tRNA-Tyr-GTA_5_end", intgroup= "Condition")
plotCounts(dds, gene="mature-tRNA-Ala-AGC_3_end", intgroup= "Condition")
plotCounts(dds, gene="mature-tRNA-Gly-TCC_5_end", intgroup= "Condition")
plotCounts(dds, gene="mature-tRNA-iMet-CAT_5_end", intgroup= "Condition")
plotCounts(dds, gene="mature-tRNA-Ser-AGA_3_end", intgroup= "Condition")
plotCounts(dds, gene="mature-tRNA-Ser-CGA_3_end", intgroup= "Condition")

#HFRunvsHFHF
plotCounts(dds, gene="mature-tRNA-Glu-CTC_CCA_end", intgroup= "Condition")





df <- read.csv("results_tiRNAs_HFHFvsChCh_without_outliers.csv", header = T )
setnames(df, "X", "Gene_id")



