#read input file
library(readxl)
library(dplyr)
countData <- read_excel("C:/Users/luisd/Dropbox/UP/Doutoramento em Ciências Biomédicas/Resultados/RNAseq/Results/sncRNA_distribution_results/cts.xlsx", 
    col_types = c("text", "numeric", "numeric", 
        "numeric", "numeric", "numeric", 
         "numeric", "numeric", "numeric", 
         "numeric", "numeric", "numeric", 
         "numeric", "numeric", "numeric", 
         "numeric", "numeric", "numeric", 
         "numeric", "numeric", "numeric", 
         "numeric", "numeric", "numeric", 
         "numeric"))

u <- pull(countData [, 1])
countData <- subset.data.frame(countData,,c(2:25))
rownames(countData) <- u

#read sample info file
colData <- read_excel("C:/Users/luisd/Dropbox/UP/Doutoramento em Ciências Biomédicas/Resultados/RNAseq/Results/sncRNA_distribution_results/coldata.xlsx")
colData <- mutate_if(colData, is.character, as.factor)
u <- pull(colData [, 1])
rownames(colData) <- u
colData <- subset.data.frame(colData,,c(2:4))
rownames(colData) <- u
remove(u)
all(rownames(colData) == colnames(countData)) #verify consistency between tables


library(DESeq2)
#normalization of all samples seting CTRL_F0 as reference
#dds_total <- DESeqDataSetFromMatrix(countData = countData,colData = colData, design = ~ Gen + Group + Gen:Group) #create DESeqDataSet object
dds_total <- DESeqDataSetFromMatrix(countData = countData,colData = colData, design = ~ GenGroup)
dds$GenGroup <- relevel(dds$GenGroup, ref = "CTRL_F0") #set reference
dds_total2 <- DESeq(dds_total)

##obtain the F0 comparisons
res_total1 <- results(dds_total2, contrast=c("GenGroup","CTRL_F0","HFD_F0"))
res_total2 <- results(dds_total2, contrast=c("GenGroup","CTRL_F0","HFDt_F0"))
res_total3 <- results(dds_total2, contrast=c("GenGroup","HFD_F0","HFDt_F0"))
write.table(as.data.frame(res_total1), file="CTRL_F0 vs. HFD_F0.csv", sep=";", col.names = NA)
write.table(as.data.frame(res_total2), file="CTRL_F0 vs. HFDt_F0.csv", sep=";", col.names = NA)
write.table(as.data.frame(res_total3), file="HFD_F0 vs. HFDt_F0.csv", sep=";", col.names = NA)

##obtain the F1 comparisons
res_total1 <- results(dds_total2, contrast=c("GenGroup","CTRL_F1","HFD_F1"))
res_total2 <- results(dds_total2, contrast=c("GenGroup","CTRL_F1","HFDt_F1"))
res_total3 <- results(dds_total2, contrast=c("GenGroup","HFD_F1","HFDt_F1"))
write.table(as.data.frame(res_total1), file="CTRL_F1 vs. HFD_F1.csv", sep=";", col.names = NA)
write.table(as.data.frame(res_total2), file="CTRL_F1 vs. HFDt_F1.csv", sep=";", col.names = NA)
write.table(as.data.frame(res_total3), file="HFD_F1 vs. HFDt_F1.csv", sep=";", col.names = NA)

##obtain the F2 comparisons
res_total1 <- results(dds_total2, contrast=c("GenGroup","CTRL_F2","HFD_F2"))
res_total2 <- results(dds_total2, contrast=c("GenGroup","CTRL_F2","HFDt_F2"))
res_total3 <- results(dds_total2, contrast=c("GenGroup","HFD_F2","HFDt_F2"))
write.table(as.data.frame(res_total1), file="CTRL_F2 vs. HFD_F2.csv", sep=";", col.names = NA)
write.table(as.data.frame(res_total2), file="CTRL_F2 vs. HFDt_F2.csv", sep=";", col.names = NA)
write.table(as.data.frame(res_total3), file="HFD_F2 vs. HFDt_F2.csv", sep=";", col.names = NA)

##obtain CTRL transgenerational comparisons
res_total1 <- results(dds_total2, contrast=c("GenGroup","CTRL_F0","CTRL_F1"))
res_total2 <- results(dds_total2, contrast=c("GenGroup","CTRL_F0","CTRL_F2"))
res_total3 <- results(dds_total2, contrast=c("GenGroup","CTRL_F1","CTRL_F2"))
write.table(as.data.frame(res_total1), file="CTRL_F0 vs. CTRL_F1.csv", sep=";", col.names = NA)
write.table(as.data.frame(res_total2), file="CTRL_F0 vs. CTRL_F2.csv", sep=";", col.names = NA)
write.table(as.data.frame(res_total3), file="CTRL_F1 vs. CTRL_F2.csv", sep=";", col.names = NA)

##obtain HFD transgenerational comparisons
res_total1 <- results(dds_total2, contrast=c("GenGroup","HFD_F0","HFD_F1"))
res_total2 <- results(dds_total2, contrast=c("GenGroup","HFD_F0","HFD_F2"))
res_total3 <- results(dds_total2, contrast=c("GenGroup","HFD_F1","HFD_F2"))
write.table(as.data.frame(res_total1), file="HFD_F0 vs. HFD_F1.csv", sep=";", col.names = NA)
write.table(as.data.frame(res_total2), file="HFD_F0 vs. HFD_F2.csv", sep=";", col.names = NA)
write.table(as.data.frame(res_total3), file="HFD_F1 vs. HFD_F2.csv", sep=";", col.names = NA)

##obtain HFDt transgenerational comparisons
res_total1 <- results(dds_total2, contrast=c("GenGroup","HFDt_F0","HFDt_F1"))
res_total2 <- results(dds_total2, contrast=c("GenGroup","HFDt_F0","HFDt_F2"))
res_total3 <- results(dds_total2, contrast=c("GenGroup","HFDt_F1","HFDt_F2"))
write.table(as.data.frame(res_total1), file="HFDt_F0 vs. HFDt_F1.csv", sep=";", col.names = NA)
write.table(as.data.frame(res_total2), file="HFDt_F0 vs. HFDt_F2.csv", sep=";", col.names = NA)
write.table(as.data.frame(res_total3), file="HFDt_F1 vs. HFDt_F2.csv", sep=";", col.names = NA)

#obtain normalized counts
counts(estimateSizeFactors(dds_f0), normalized=TRUE)
counts(estimateSizeFactors(dds_f1), normalized=TRUE)
counts(estimateSizeFactors(dds_f2), normalized=TRUE)
counts(estimateSizeFactors(dds_total), normalized=TRUE)
write.table(counts(estimateSizeFactors(dds_total), normalized=TRUE), "normalized_counts.csv", sep=";", col.names = NA)

#obtain base means (normalizing factor)
rowMeans(counts(estimateSizeFactors(dds_f0), normalized=TRUE))
rowMeans(counts(estimateSizeFactors(dds_f1), normalized=TRUE))
rowMeans(counts(estimateSizeFactors(dds_f2), normalized=TRUE))

write.csv(as.data.frame(counts(estimateSizeFactors(dds_f0), normalized=TRUE)), file="Bar_data_F0.csv")
write.csv(as.data.frame(counts(estimateSizeFactors(dds_f1), normalized=TRUE)), file="Bar_data_F1.csv")
write.csv(as.data.frame(counts(estimateSizeFactors(dds_f2), normalized=TRUE)), file="Bar_data_F2.csv")
write.table(as.data.frame(counts(estimateSizeFactors(dds_total), normalized=TRUE)), "normalized_counts.csv", sep=";", col.names = NA)


#bar plot
library(ggplot2)
require(scales)
library(readxl)
bar_plot <- read_excel("bar_chart_data.xlsx", col_types = c("text", "text", "text", "skip", "skip", "skip", "skip", "skip", "skip", "numeric", "numeric"))
bar_plot$group <- factor(bar_plot$group, levels = c("CTRL", "HFD", "HFDt"))
#bar_plot$gen <- factor(bar_plot$gen, levels = c("F0 (n[group] = 3)", "F1 (n[group] = 2)", "F2 (n[group] = 3)"))
bar_plot$gen <- factor(bar_plot$gen, levels = c("F0", "F1", "F2"))
bar_plot$gen <- factor(bar_plot$gen, labels = c(bquote("F0"~(n[group]==3)), bquote("F1"~(n[group]==2)), bquote("F2"~(n[group]==3))))

#new graph
tiff("Figure 1 - sncRNA categories.tiff", width=3200, height=1000, units="px", res = 300)
ggplot(data=bar_plot, aes(y=mean, x=gene_grp, fill=group)) + 
geom_bar(stat="identity", position=position_dodge(preserve = 'single')) +
geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.5,position=position_dodge(.9)) + 
facet_grid(~bar_plot$gen, labeller = label_parsed) +
scale_y_continuous(trans = 'log10', breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))) +
scale_fill_brewer(palette="Set1") +
labs(x = "RNA biotype", y = "Normalized counts", fill = "Group") +
theme(axis.title = element_text(size = 12, face="bold")) + 
theme(legend.title = element_text(size = 12, face="bold")) + 
theme(strip.text.x = element_text(size = 14, face="bold"))
dev.off()