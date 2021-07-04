#Annotate sncRNA targets using GO terms

#will clear all objects includes hidden objects.
rm(list = ls(all.names = TRUE))

#Load packages
library(readr)
library(tidyr)
library(topGO) #annotation tool
library(Mus.musculus) #mouse annotation
#library(org.Mm.eg.db) #loaded by the previous

# custom filter function
 count.filter <- structure(function (allScore)                                      
 {                                                                  
     return(allScore > 1)                                        
 }, source = c("function(allScore) {", "  return(allScore > 1)", 
 "}")) 

#Generation F0: CTRL vs. HFDt
##load targets of DEGs
F0_CTRL_HFDt_Targets <- read_csv("F0_CTRL_vs_HFDt - sRNAtools  Targets.csv")
F0_CTRL_HFDt_Targets <- F0_CTRL_HFDt_Targets %>% separate_rows(Targets)

##data structure
test.genes_F0_CTRL_HFDt <- F0_CTRL_HFDt_Targets[4]
test.genes_F0_CTRL_HFDt$Targets <- gsub("\\..*","",test.genes_F0_CTRL_HFDt$Targets)
test.genes_F0_CTRL_HFDt <- as.data.frame(table(test.genes_F0_CTRL_HFDt$Targets))
colnames(test.genes_F0_CTRL_HFDt) <- c("Targets","counts")

genecounts_F0_CTRL_HFDt <- test.genes_F0_CTRL_HFDt[[2]]
names(genecounts_F0_CTRL_HFDt) <- test.genes_F0_CTRL_HFDt[[1]]

##prepare topGO dataobjects
topGO.BP_F0_CTRL_vs_HFDt <- new("topGOdata", ontology = "BP", allGenes = genecounts_F0_CTRL_HFDt, geneSel = count.filter,
description = "GO analysis BP: CTRL_F0 vs. HFDt_F0", nodeSize = 10, annotationFun = annFUN.org, mapping = "org.Mm.eg.db", ID = "ensembl")

topGO.MF_F0_CTRL_vs_HFDt <- new("topGOdata", ontology = "MF", allGenes = genecounts_F0_CTRL_HFDt, geneSel = count.filter,
description = "GO analysis MF: CTRL_F0 vs. HFDt_F0", nodeSize = 10, annotationFun = annFUN.org, mapping = "org.Mm.eg.db", ID = "ensembl")

topGO.CC_F0_CTRL_vs_HFDt <- new("topGOdata", ontology = "CC", allGenes = genecounts_F0_CTRL_HFDt, geneSel = count.filter,
description = "GO analysis CC: CTRL_F0 vs. HFDt_F0", nodeSize = 10, annotationFun = annFUN.org, mapping = "org.Mm.eg.db", ID = "ensembl")

##perform GO test - weight01 & Fisher
runGO.BP_F0_CTRL_vs_HFDt <- runTest(topGO.BP_F0_CTRL_vs_HFDt, algorithm = "weight01", statistic = "fisher")#, topNodes = 5)
runGO.MF_F0_CTRL_vs_HFDt <- runTest(topGO.MF_F0_CTRL_vs_HFDt, algorithm = "weight01", statistic = "fisher")
runGO.CC_F0_CTRL_vs_HFDt <- runTest(topGO.CC_F0_CTRL_vs_HFDt, algorithm = "weight01", statistic = "fisher")

##generate result tables
GO.BP_F0_CTRL_vs_HFDt <- GenTable(topGO.BP_F0_CTRL_vs_HFDt, Weight01 = runGO.BP_F0_CTRL_vs_HFDt, orderBy = "Weight01", topNodes = 30)
GO.MF_F0_CTRL_vs_HFDt <- GenTable(topGO.MF_F0_CTRL_vs_HFDt, Weight01 = runGO.MF_F0_CTRL_vs_HFDt, orderBy = "Weight01", topNodes = 30)
GO.CC_F0_CTRL_vs_HFDt <- GenTable(topGO.CC_F0_CTRL_vs_HFDt, Weight01 = runGO.CC_F0_CTRL_vs_HFDt, orderBy = "Weight01", topNodes = 30)

##Write results tables
write.csv2(GO.BP_F0_CTRL_vs_HFDt,"GO Biological Process CTRL-F0 vs HFDt-F0.csv", row.names = FALSE)
write.csv2(GO.MF_F0_CTRL_vs_HFDt,"GO Molecular Function CTRL-F0 vs HFDt-F0.csv", row.names = FALSE)
write.csv2(GO.CC_F0_CTRL_vs_HFDt,"GO Cellular Component CTRL-F0 vs HFDt-F0.csv", row.names = FALSE)


#Generation F1: CTRL vs. HFD
##load targets of DEGs
F1_CTRL_HFD_Targets <- read_csv("F1_CTRL_vs_HFD - sRNAtools  Targets.csv")
F1_CTRL_HFD_Targets <- F1_CTRL_HFD_Targets %>% separate_rows(Targets)

##data structure
test.genes_F1_CTRL_HFD <- F1_CTRL_HFD_Targets[4]
test.genes_F1_CTRL_HFD$Targets <- gsub("\\..*","",test.genes_F1_CTRL_HFD$Targets)
test.genes_F1_CTRL_HFD <- as.data.frame(table(test.genes_F1_CTRL_HFD$Targets))
colnames(test.genes_F1_CTRL_HFD) <- c("Targets","counts")

genecounts_F1_CTRL_HFD <- test.genes_F1_CTRL_HFD[[2]]
names(genecounts_F1_CTRL_HFD) <- test.genes_F1_CTRL_HFD[[1]]

##prepare topGO dataobjects
topGO.BP_F1_CTRL_vs_HFD <- new("topGOdata", ontology = "BP", allGenes = genecounts_F1_CTRL_HFD, geneSel = count.filter,
description = "GO analysis BP: CTRL_F1 vs. HFD_F1", nodeSize = 10, annotationFun = annFUN.org, mapping = "org.Mm.eg.db", ID = "ensembl")

topGO.MF_F1_CTRL_vs_HFD <- new("topGOdata", ontology = "MF", allGenes = genecounts_F1_CTRL_HFD, geneSel = count.filter,
description = "GO analysis MF: CTRL_F1 vs. HFD_F1", nodeSize = 10, annotationFun = annFUN.org, mapping = "org.Mm.eg.db", ID = "ensembl")

topGO.CC_F1_CTRL_vs_HFD <- new("topGOdata", ontology = "CC", allGenes = genecounts_F1_CTRL_HFD, geneSel = count.filter,
description = "GO analysis CC: CTRL_F1 vs. HFD_F1", nodeSize = 10, annotationFun = annFUN.org, mapping = "org.Mm.eg.db", ID = "ensembl")

##perform GO test - weight01 & Fisher
runGO.BP_F1_CTRL_vs_HFD <- runTest(topGO.BP_F1_CTRL_vs_HFD, algorithm = "weight01", statistic = "fisher")#, topNodes = 5)
runGO.MF_F1_CTRL_vs_HFD <- runTest(topGO.MF_F1_CTRL_vs_HFD, algorithm = "weight01", statistic = "fisher")
runGO.CC_F1_CTRL_vs_HFD <- runTest(topGO.CC_F1_CTRL_vs_HFD, algorithm = "weight01", statistic = "fisher")

#Generation F1: CTRL vs. HFDt
##load targets of DEGs
F1_CTRL_HFDt_Targets <- read_csv("F1_CTRL_vs_HFDt - sRNAtools  Targets.csv")
F1_CTRL_HFDt_Targets <- F1_CTRL_HFDt_Targets %>% separate_rows(Targets)

##data structure
test.genes_F1_CTRL_HFDt <- F1_CTRL_HFDt_Targets[4]
test.genes_F1_CTRL_HFDt$Targets <- gsub("\\..*","",test.genes_F1_CTRL_HFDt$Targets)
test.genes_F1_CTRL_HFDt <- as.data.frame(table(test.genes_F1_CTRL_HFDt$Targets))
colnames(test.genes_F1_CTRL_HFDt) <- c("Targets","counts")

genecounts_F1_CTRL_HFDt <- test.genes_F1_CTRL_HFDt[[2]]
names(genecounts_F1_CTRL_HFDt) <- test.genes_F1_CTRL_HFDt[[1]]

##prepare topGO dataobjects
topGO.BP_F1_CTRL_vs_HFDt <- new("topGOdata", ontology = "BP", allGenes = genecounts_F1_CTRL_HFDt, geneSel = count.filter,
description = "GO analysis BP: CTRL_F1 vs. HFDt_F1", nodeSize = 10, annotationFun = annFUN.org, mapping = "org.Mm.eg.db", ID = "ensembl")

topGO.MF_F1_CTRL_vs_HFDt <- new("topGOdata", ontology = "MF", allGenes = genecounts_F1_CTRL_HFDt, geneSel = count.filter,
description = "GO analysis MF: CTRL_F1 vs. HFDt_F1", nodeSize = 10, annotationFun = annFUN.org, mapping = "org.Mm.eg.db", ID = "ensembl")

topGO.CC_F1_CTRL_vs_HFDt <- new("topGOdata", ontology = "CC", allGenes = genecounts_F1_CTRL_HFDt, geneSel = count.filter,
description = "GO analysis CC: CTRL_F1 vs. HFDt_F1", nodeSize = 10, annotationFun = annFUN.org, mapping = "org.Mm.eg.db", ID = "ensembl")

##perform GO test - weight01 & Fisher
runGO.BP_F1_CTRL_vs_HFDt <- runTest(topGO.BP_F1_CTRL_vs_HFDt, algorithm = "weight01", statistic = "fisher")#, topNodes = 5)
runGO.MF_F1_CTRL_vs_HFDt <- runTest(topGO.MF_F1_CTRL_vs_HFDt, algorithm = "weight01", statistic = "fisher")
runGO.CC_F1_CTRL_vs_HFDt <- runTest(topGO.CC_F1_CTRL_vs_HFDt, algorithm = "weight01", statistic = "fisher")

##generate result tables
GO.BP_F1_CTRL_vs_HFDt <- GenTable(topGO.BP_F1_CTRL_vs_HFDt, Weight01 = runGO.BP_F1_CTRL_vs_HFDt, orderBy = "Weight01", topNodes = 30)
GO.MF_F1_CTRL_vs_HFDt <- GenTable(topGO.MF_F1_CTRL_vs_HFDt, Weight01 = runGO.MF_F1_CTRL_vs_HFDt, orderBy = "Weight01", topNodes = 30)
GO.CC_F1_CTRL_vs_HFDt <- GenTable(topGO.CC_F1_CTRL_vs_HFDt, Weight01 = runGO.CC_F1_CTRL_vs_HFDt, orderBy = "Weight01", topNodes = 30)

##Write results tables
write.csv2(GO.BP_F1_CTRL_vs_HFDt,"GO Biological Process CTRL-F1 vs HFDt-F1.csv", row.names = FALSE)
write.csv2(GO.MF_F1_CTRL_vs_HFDt,"GO Molecular Function CTRL-F1 vs HFDt-F1.csv", row.names = FALSE)
write.csv2(GO.CC_F1_CTRL_vs_HFDt,"GO Cellular Component CTRL-F1 vs HFDt-F1.csv", row.names = FALSE)


#Generation F1: HFD vs. HFDt
##load targets of DEGs
F1_HFD_HFDt_Targets <- read_csv("F1_HFD_vs_HFDt - sRNAtools  Targets.csv")
F1_HFD_HFDt_Targets <- F1_HFD_HFDt_Targets %>% separate_rows(Targets)

##data structure
test.genes_F1_HFD_HFDt <- F1_HFD_HFDt_Targets[4]
test.genes_F1_HFD_HFDt$Targets <- gsub("\\..*","",test.genes_F1_HFD_HFDt$Targets)
test.genes_F1_HFD_HFDt <- as.data.frame(table(test.genes_F1_HFD_HFDt$Targets))
colnames(test.genes_F1_HFD_HFDt) <- c("Targets","counts")

genecounts_F1_HFD_HFDt <- test.genes_F1_HFD_HFDt[[2]]
names(genecounts_F1_HFD_HFDt) <- test.genes_F1_HFD_HFDt[[1]]

##prepare topGO dataobjects
topGO.BP_F1_HFD_vs_HFDt <- new("topGOdata", ontology = "BP", allGenes = genecounts_F1_HFD_HFDt, geneSel = count.filter,
description = "GO analysis BP: HFD_F1 vs. HFDt_F1", nodeSize = 10, annotationFun = annFUN.org, mapping = "org.Mm.eg.db", ID = "ensembl")

topGO.MF_F1_HFD_vs_HFDt <- new("topGOdata", ontology = "MF", allGenes = genecounts_F1_HFD_HFDt, geneSel = count.filter,
description = "GO analysis MF: HFD_F1 vs. HFDt_F1", nodeSize = 10, annotationFun = annFUN.org, mapping = "org.Mm.eg.db", ID = "ensembl")

topGO.CC_F1_HFD_vs_HFDt <- new("topGOdata", ontology = "CC", allGenes = genecounts_F1_HFD_HFDt, geneSel = count.filter,
description = "GO analysis CC: HFD_F1 vs. HFDt_F1", nodeSize = 10, annotationFun = annFUN.org, mapping = "org.Mm.eg.db", ID = "ensembl")

##perform GO test - weight01 & Fisher
runGO.BP_F1_HFD_vs_HFDt <- runTest(topGO.BP_F1_HFD_vs_HFDt, algorithm = "weight01", statistic = "fisher")#, topNodes = 5)
runGO.MF_F1_HFD_vs_HFDt <- runTest(topGO.MF_F1_HFD_vs_HFDt, algorithm = "weight01", statistic = "fisher")
runGO.CC_F1_HFD_vs_HFDt <- runTest(topGO.CC_F1_HFD_vs_HFDt, algorithm = "weight01", statistic = "fisher")

##generate result tables
GO.BP_F1_HFD_vs_HFDt <- GenTable(topGO.BP_F1_HFD_vs_HFDt, Weight01 = runGO.BP_F1_HFD_vs_HFDt, orderBy = "Weight01", topNodes = 30)
GO.MF_F1_HFD_vs_HFDt <- GenTable(topGO.MF_F1_HFD_vs_HFDt, Weight01 = runGO.MF_F1_HFD_vs_HFDt, orderBy = "Weight01", topNodes = 30)
GO.CC_F1_HFD_vs_HFDt <- GenTable(topGO.CC_F1_HFD_vs_HFDt, Weight01 = runGO.CC_F1_HFD_vs_HFDt, orderBy = "Weight01", topNodes = 30)

##Write results tables
write.csv2(GO.BP_F1_HFD_vs_HFDt,"GO Biological Process HFD-F1 vs HFDt-F1.csv", row.names = FALSE)
write.csv2(GO.MF_F1_HFD_vs_HFDt,"GO Molecular Function HFD-F1 vs HFDt-F1.csv", row.names = FALSE)
write.csv2(GO.CC_F1_HFD_vs_HFDt,"GO Cellular Component HFD-F1 vs HFDt-F1.csv", row.names = FALSE)


#Generation F2: CTRL vs. HFDt
##load targets of DEGs
F2_CTRL_HFDt_Targets <- read_csv("F2_CTRL_vs_HFDt - sRNAtools  Targets.csv")
F2_CTRL_HFDt_Targets <- F2_CTRL_HFDt_Targets %>% separate_rows(Targets)

##data structure
test.genes_F2_CTRL_HFDt <- F2_CTRL_HFDt_Targets[4]
test.genes_F2_CTRL_HFDt$Targets <- gsub("\\..*","",test.genes_F2_CTRL_HFDt$Targets)
test.genes_F2_CTRL_HFDt <- as.data.frame(table(test.genes_F2_CTRL_HFDt$Targets))
colnames(test.genes_F2_CTRL_HFDt) <- c("Targets","counts")

genecounts_F2_CTRL_HFDt <- test.genes_F2_CTRL_HFDt[[2]]
names(genecounts_F2_CTRL_HFDt) <- test.genes_F2_CTRL_HFDt[[1]]

##prepare topGO dataobjects
topGO.BP_F2_CTRL_vs_HFDt <- new("topGOdata", ontology = "BP", allGenes = genecounts_F2_CTRL_HFDt, geneSel = count.filter,
description = "GO analysis BP: CTRL_F2 vs. HFDt_F2", nodeSize = 10, annotationFun = annFUN.org, mapping = "org.Mm.eg.db", ID = "ensembl")

topGO.MF_F2_CTRL_vs_HFDt <- new("topGOdata", ontology = "MF", allGenes = genecounts_F2_CTRL_HFDt, geneSel = count.filter,
description = "GO analysis MF: CTRL_F2 vs. HFDt_F2", nodeSize = 10, annotationFun = annFUN.org, mapping = "org.Mm.eg.db", ID = "ensembl")

topGO.CC_F2_CTRL_vs_HFDt <- new("topGOdata", ontology = "CC", allGenes = genecounts_F2_CTRL_HFDt, geneSel = count.filter,
description = "GO analysis CC: CTRL_F2 vs. HFDt_F2", nodeSize = 10, annotationFun = annFUN.org, mapping = "org.Mm.eg.db", ID = "ensembl")

##perform GO test - weight01 & Fisher
runGO.BP_F2_CTRL_vs_HFDt <- runTest(topGO.BP_F2_CTRL_vs_HFDt, algorithm = "weight01", statistic = "fisher")#, topNodes = 5)
runGO.MF_F2_CTRL_vs_HFDt <- runTest(topGO.MF_F2_CTRL_vs_HFDt, algorithm = "weight01", statistic = "fisher")
runGO.CC_F2_CTRL_vs_HFDt <- runTest(topGO.CC_F2_CTRL_vs_HFDt, algorithm = "weight01", statistic = "fisher")

##generate result tables
GO.BP_F2_CTRL_vs_HFDt <- GenTable(topGO.BP_F2_CTRL_vs_HFDt, Weight01 = runGO.BP_F2_CTRL_vs_HFDt, orderBy = "Weight01", topNodes = 30)
GO.MF_F2_CTRL_vs_HFDt <- GenTable(topGO.MF_F2_CTRL_vs_HFDt, Weight01 = runGO.MF_F2_CTRL_vs_HFDt, orderBy = "Weight01", topNodes = 30)
GO.CC_F2_CTRL_vs_HFDt <- GenTable(topGO.CC_F2_CTRL_vs_HFDt, Weight01 = runGO.CC_F2_CTRL_vs_HFDt, orderBy = "Weight01", topNodes = 30)

##Write results tables
write.csv2(GO.BP_F2_CTRL_vs_HFDt,"GO Biological Process CTRL-F2 vs HFDt-F2.csv", row.names = FALSE)
write.csv2(GO.MF_F2_CTRL_vs_HFDt,"GO Molecular Function CTRL-F2 vs HFDt-F2.csv", row.names = FALSE)
write.csv2(GO.CC_F2_CTRL_vs_HFDt,"GO Cellular Component CTRL-F2 vs HFDt-F2.csv", row.names = FALSE)


#Generation F2: HFD vs. HFDt
##load targets of DEGs
F2_HFD_HFDt_Targets <- read_csv("F2_HFD_vs_HFDt - sRNAtools  Targets.csv")
F2_HFD_HFDt_Targets <- F2_HFD_HFDt_Targets %>% separate_rows(Targets)

##data structure
test.genes_F2_HFD_HFDt <- F2_HFD_HFDt_Targets[4]
test.genes_F2_HFD_HFDt$Targets <- gsub("\\..*","",test.genes_F2_HFD_HFDt$Targets)
test.genes_F2_HFD_HFDt <- as.data.frame(table(test.genes_F2_HFD_HFDt$Targets))
colnames(test.genes_F2_HFD_HFDt) <- c("Targets","counts")

genecounts_F2_HFD_HFDt <- test.genes_F2_HFD_HFDt[[2]]
names(genecounts_F2_HFD_HFDt) <- test.genes_F2_HFD_HFDt[[1]]

##prepare topGO dataobjects
topGO.BP_F2_HFD_vs_HFDt <- new("topGOdata", ontology = "BP", allGenes = genecounts_F2_HFD_HFDt, geneSel = count.filter,
description = "GO analysis BP: HFD_F2 vs. HFDt_F2", nodeSize = 10, annotationFun = annFUN.org, mapping = "org.Mm.eg.db", ID = "ensembl")

topGO.MF_F2_HFD_vs_HFDt <- new("topGOdata", ontology = "MF", allGenes = genecounts_F2_HFD_HFDt, geneSel = count.filter,
description = "GO analysis MF: HFD_F2 vs. HFDt_F2", nodeSize = 10, annotationFun = annFUN.org, mapping = "org.Mm.eg.db", ID = "ensembl")

topGO.CC_F2_HFD_vs_HFDt <- new("topGOdata", ontology = "CC", allGenes = genecounts_F2_HFD_HFDt, geneSel = count.filter,
description = "GO analysis CC: HFD_F2 vs. HFDt_F2", nodeSize = 10, annotationFun = annFUN.org, mapping = "org.Mm.eg.db", ID = "ensembl")

##perform GO test - weight01 & Fisher
runGO.BP_F2_HFD_vs_HFDt <- runTest(topGO.BP_F2_HFD_vs_HFDt, algorithm = "weight01", statistic = "fisher")#, topNodes = 5)
runGO.MF_F2_HFD_vs_HFDt <- runTest(topGO.MF_F2_HFD_vs_HFDt, algorithm = "weight01", statistic = "fisher")
runGO.CC_F2_HFD_vs_HFDt <- runTest(topGO.CC_F2_HFD_vs_HFDt, algorithm = "weight01", statistic = "fisher")

##generate result tables
GO.BP_F2_HFD_vs_HFDt <- GenTable(topGO.BP_F2_HFD_vs_HFDt, Weight01 = runGO.BP_F2_HFD_vs_HFDt, orderBy = "Weight01", topNodes = 30)
GO.MF_F2_HFD_vs_HFDt <- GenTable(topGO.MF_F2_HFD_vs_HFDt, Weight01 = runGO.MF_F2_HFD_vs_HFDt, orderBy = "Weight01", topNodes = 30)
GO.CC_F2_HFD_vs_HFDt <- GenTable(topGO.CC_F2_HFD_vs_HFDt, Weight01 = runGO.CC_F2_HFD_vs_HFDt, orderBy = "Weight01", topNodes = 30)

##Write results tables
write.csv2(GO.BP_F2_HFD_vs_HFDt,"GO Biological Process HFD-F2 vs HFDt-F2.csv", row.names = FALSE)
write.csv2(GO.MF_F2_HFD_vs_HFDt,"GO Molecular Function HFD-F2 vs HFDt-F2.csv", row.names = FALSE)
write.csv2(GO.CC_F2_HFD_vs_HFDt,"GO Cellular Component HFD-F2 vs HFDt-F2.csv", row.names = FALSE)


#test database
# xx <- annFUN.org("BP", mapping = "org.Mm.eg.db", ID = "ensembl")
# head(xx)

#perform GO tests
# go_test <- runTest(test, algorithm = "weight01", #the more balance algorithm. 'classic' (unhierarchical), 'lea' (?), 'elim' (more conservative), 'parentchild' (close to weight-9
# statistic = "fisher")


#create charts
#input files were combined and edited in Excel
library(readxl)
library(ggplot2)
library(scales)

#Biological process - Transgenerational
GO_results_BP <- read_excel("GO_results_Biological process - All generations.xlsx", col_types = c("text", "text", "text", "text", "numeric", "numeric", "numeric"))
GO_results_BP$group <- factor(GO_results_BP$group, levels = c("CTRL vs. HFD", "CTRL vs. HFDt", "HFD vs. HFDt"))
GO_results_BP$Gen <- as.factor(GO_results_BP$Gen)

	
reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))
}

tiff("Figure 5 - GO Biological Process.tiff", width=800, height=600, units="px")	
ggplot(data=GO_results_BP, aes(y=Term, x=p.value, fill=group)) + 
geom_bar(stat="identity", position=position_dodge(preserve = 'single')) +
facet_grid(rows = GO_results_BP$Gen, scales="free_y", space="free") +
scale_x_continuous(trans = reverselog_trans(10), breaks = breaks_log(n = 7, base = 10)) +
scale_fill_brewer(palette="Set1") +
labs(x = "-log(p-value)", y = "GO terms", fill = "Group") +
theme(axis.title = element_text(size = 12, face="bold")) + 
theme(legend.title = element_text(size = 12, face="bold")) + 
theme(strip.text.y = element_text(size = 12, face="bold"))
dev.off()

#Molecular function - Transgenerational
GO_results_MF <- read_excel("GO_results_Molecular function - All generations.xlsx", col_types = c("text", "text", "text", "text", "numeric", "numeric", "numeric"))
GO_results_MF$group <- factor(GO_results_MF$group, levels = c("CTRL vs. HFD", "CTRL vs. HFDt", "HFD vs. HFDt"))
GO_results_MF$Gen <- as.factor(GO_results_MF$Gen)

tiff("Figure S1 - GO Molecular Function.tiff", width=800, height=550, units="px")	
ggplot(data=GO_results_MF, aes(y=Term, x=p.value, fill=group)) + 
geom_bar(stat="identity", position=position_dodge(preserve = 'single')) +
facet_grid(rows = GO_results_MF$Gen, scales="free_y", space="free") +
scale_x_continuous(trans = reverselog_trans(10), breaks = breaks_log(n = 7, base = 10)) +
scale_fill_brewer(palette="Set1") +
labs(x = "-log(p-value)", y = "GO terms", fill = "Group") +
theme(axis.title = element_text(size = 12, face="bold")) + 
theme(legend.title = element_text(size = 12, face="bold")) + 
theme(strip.text.y = element_text(size = 12, face="bold"))
dev.off()

#Cellular component - Transgenerational
GO_results_CC <- read_excel("GO_results_Cellular component - All generations.xlsx", col_types = c("text", "text", "text", "text", "numeric", "numeric", "numeric"))
GO_results_CC$group <- factor(GO_results_CC$group, levels = c("CTRL vs. HFD", "CTRL vs. HFDt", "HFD vs. HFDt"))
GO_results_CC$Gen <- as.factor(GO_results_CC$Gen)

tiff("Figure S2 - GO Cellular Component.tiff", width=800, height=500, units="px")	
ggplot(data=GO_results_CC, aes(y=Term, x=p.value, fill=group)) + 
geom_bar(stat="identity", position=position_dodge(preserve = 'single')) +
facet_grid(rows = GO_results_CC$Gen, scales="free_y", space="free") +
scale_x_continuous(trans = reverselog_trans(10), breaks = breaks_log(n = 7, base = 10)) +
scale_fill_brewer(palette="Set1") +
labs(x = "-log(p-value)", y = "GO terms", fill = "Group") +
theme(axis.title = element_text(size = 12, face="bold")) + 
theme(legend.title = element_text(size = 12, face="bold")) + 
theme(strip.text.y = element_text(size = 12, face="bold"))
dev.off()


#Euler diagram DEGs accross generations/diets
euler.ult <- c("F0 CTRL vs. HFDt" = 19, "F1 CTRL vs. HFD" = 86, "F1 CTRL vs. HFDt" = 24, "F1 HFD vs. HFDt" = 5, "F2 CTRL vs. HFDt" = 3, "F2 HFD vs. HFDt" = 4,
            "F1 CTRL vs. HFD&F1 CTRL vs. HFDt" = 10, "F1 CTRL vs. HFD&F1 HFD vs. HFDt" = 5, "F1 CTRL vs. HFDt&F1 HFD vs. HFDt" = 1,
                    "F2 CTRL vs. HFDt&F2 HFD vs. HFDt" = 2)
					

##print version
svg(filename = "Figure 6 - Euler diagram_v1.svg", width=9, height=6, pointsize = 12, family = "sans")
plot(
euler(euler.ult, shape = "ellipse"), counts = TRUE, 
quantities = list(type = c("counts", "percent"), round=2, cex=0.8), 
fills =list(fill=c("#377EB8", "#E41A1C", "#377EB8", "#4DAF4A", "#377EB8", "#4DAF4A")), alpha = 0.3, c("#377EB8", "#E41A1C", "#377EB8", "#4DAF4A", "#377EB8", "#4DAF4A"), alpha = 0.3,
labels = list(font = 2, cex=1),
)
dev.off()

##colorblind version
svg(filename = "Figure 6 - Euler diagram_v2.svg", width=9, height=6, pointsize = 12, family = "sans")
plot(
euler(euler.ult, shape = "ellipse"), counts = TRUE, 
quantities = list(type = c("counts", "percent"), round=2, cex=0.8), 
labels = list(font = 2, cex=1),
)
dev.off()


#código antigo
#Generation F0
GO_results_Gen0 <- read_excel("GO_results_Gen0.xlsx", col_types = c("text", "text", "text", "text", "numeric", "numeric"))
GO_results_Gen0$group <- factor(GO_results_Gen0$group, levels = c("CTRL vs. HFD", "CTRL vs. HFDt", "HFD vs. HFDt"))
GO_results_Gen0$GO.analysis <- factor(GO_results_Gen0$GO.analysis, levels = c("Biological Process", "Molecular Function", "Cellular Component"))

#bar graph (save as 1200 x 664)
ggplot(data=GO_results_Gen0, aes(y=Term, x=norm.annot, fill=group)) + 
geom_bar(stat="identity", position=position_dodge(preserve = 'single')) +
facet_grid(rows = GO_results_Gen0$GO.analysis, scales="free_y", space="free") +
scale_x_continuous(breaks = pretty(GO_results_Gen0$norm.annot, n = 10)) +
scale_fill_brewer(palette="Set1") +
labs(x = "GO annotations (per 1000 targets)", y = "GO terms", fill = "Group") +
theme(axis.title = element_text(size = 12, face="bold")) + 
theme(legend.title = element_text(size = 12, face="bold")) + 
theme(strip.text.y = element_text(size = 12, face="bold"))

#Generation F1
GO_results_Gen1 <- read_excel("GO_results_Gen1.xlsx", col_types = c("text", "text", "text", "text", "numeric", "numeric"))
GO_results_Gen1$group <- factor(GO_results_Gen1$group, levels = c("CTRL vs. HFD", "CTRL vs. HFDt", "HFD vs. HFDt"))
GO_results_Gen1$GO.analysis <- factor(GO_results_Gen1$GO.analysis, levels = c("Biological Process", "Molecular Function", "Cellular Component"))

#bar graph (save as 1200 x 764)
ggplot(data=GO_results_Gen1, aes(y=Term, x=norm.annot, fill=group)) + 
geom_bar(stat="identity", position=position_dodge(preserve = 'single')) +
facet_grid(rows = GO_results_Gen1$GO.analysis, scales="free_y", space="free") +
scale_x_continuous(breaks = pretty(GO_results_Gen1$norm.annot, n = 10)) +
scale_fill_brewer(palette="Set1") +
labs(x = "GO annotations (per 1000 targets)", y = "GO terms", fill = "Group") +
theme(axis.title = element_text(size = 12, face="bold")) + 
theme(legend.title = element_text(size = 12, face="bold")) + 
theme(strip.text.y = element_text(size = 12, face="bold"))

#Generation F2
GO_results_Gen2 <- read_excel("GO_results_Gen2.xlsx", col_types = c("text", "text", "text", "text", "numeric", "numeric"))
GO_results_Gen2$group <- factor(GO_results_Gen2$group, levels = c("CTRL vs. HFD", "CTRL vs. HFDt", "HFD vs. HFDt"))
GO_results_Gen2$GO.analysis <- factor(GO_results_Gen2$GO.analysis, levels = c("Biological Process", "Molecular Function", "Cellular Component"))

#bar graph (save as 1200 x 764)
ggplot(data=GO_results_Gen2, aes(y=Term, x=norm.annot, fill=group)) + 
geom_bar(stat="identity", position=position_dodge(preserve = 'single')) +
facet_grid(rows = GO_results_Gen2$GO.analysis, scales="free_y", space="free") +
scale_x_continuous(breaks = pretty(GO_results_Gen2$norm.annot, n = 10)) +
scale_fill_brewer(palette="Set1") +
labs(x = "GO annotations (per 1000 targets)", y = "GO terms", fill = "Group") +
theme(axis.title = element_text(size = 12, face="bold")) + 
theme(legend.title = element_text(size = 12, face="bold")) + 
theme(strip.text.y = element_text(size = 12, face="bold"))

#Biological process - Transgenerational
GO_results_BP <- read_excel("GO_results_Biological process - All generations.xlsx", col_types = c("text", "text", "text", "text", "numeric", "numeric"))
GO_results_BP$group <- factor(GO_results_BP$group, levels = c("CTRL vs. HFD", "CTRL vs. HFDt", "HFD vs. HFDt"))
GO_results_BP$GO.analysis <- factor(GO_results_BP$GO.analysis, levels = c("F0 Generation", "F1 Generation", "F2 Generation"))

#bar graph (save as 1200 x 764)
ggplot(data=GO_results_BP, aes(y=Term, x=norm.annot, fill=group)) + 
geom_bar(stat="identity", position=position_dodge(preserve = 'single')) +
facet_grid(rows = GO_results_BP$GO.analysis, scales="free_y", space="free") +
scale_x_continuous(breaks = pretty(GO_results_BP$norm.annot, n = 10)) +
scale_fill_brewer(palette="Set1") +
labs(x = "GO annotations (per 1000 targets)", y = "GO terms", fill = "Group") +
theme(axis.title = element_text(size = 12, face="bold")) + 
theme(legend.title = element_text(size = 12, face="bold")) + 
theme(strip.text.y = element_text(size = 12, face="bold"))
  
#other options:
#ylab("GO terms") + xlab("Number of annotations") +
#geom_text(aes(label=Annotated), vjust=1.0, color="black", position = position_dodge(0.0), size=3.5) +
#+ scale_fill_discrete(name="Type", breaks=c("TRUE", "FALSE", "NA"), labels=c("existing", "not existing", "missing values")) + scale_x_discrete(labels = c("bad", "good"))
#facet_wrap(~GO.analysis, scales="free_y") + 
#theme_minimal()
#coord_flip()

#Generation F0 vs. F1 - CTRL vs. HFD
GO_results_Gen0vs1_HFD <- read_excel("GO_results_Gen0_Gen1-HFD.xlsx", col_types = c("text", "text", "text", "text", "numeric", "numeric"))
GO_results_Gen0vs1_HFD$group <- factor(GO_results_Gen0vs1_HFD$group, levels = c("Generation F0", "Generation F0"))
GO_results_Gen0vs1_HFD$GO.analysis <- factor(GO_results_Gen0vs1_HFD$GO.analysis, levels = c("Biological Process", "Molecular Function", "Cellular Component"))

#bar graph (save as 1200 x 764)
ggplot(data=GO_results_Gen0vs1_HFD, aes(y=Term, x=norm.annot, fill=group)) + 
geom_bar(stat="identity", position=position_dodge(preserve = 'single')) +
facet_grid(rows = GO_results_Gen0vs1_HFD$GO.analysis, scales="free_y", space="free") +
scale_x_continuous(breaks = pretty(GO_results_Gen0vs1_HFD$norm.annot, n = 10)) +
scale_fill_brewer(palette="Paired") +
labs(x = "GO annotations (per 1000 targets)", y = "GO terms", fill = "Generation") +
theme(axis.title = element_text(size = 12, face="bold")) + 
theme(legend.title = element_text(size = 12, face="bold")) + 
theme(strip.text.y = element_text(size = 12, face="bold"))