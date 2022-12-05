setwd("E:/Data Mining for TME/TCGA-LIHC/TCGA-Immune/20210427/Heatmap/Whole")
library(ggplot2)
library(immunedeconv)
library(tidyverse)
library(pheatmap)
library(circlize)
library(GetoptLong)
library(dendextend)
library(cluster)
# Deconvolution algorithms of TIMER, CIBERSORT, and MCPCOUNTER
expr<-read.table(file = "HiSeqV2.txt", sep = "\t",header = T,row.names = 1, as.is=TRUE)
SAM <- read.table("SAMPLE.txt",sep = "\t",header=T)
expr<-expr[,SAM$sampleID]
res_mcp_counter  <- deconvolute(expr, method="mcp_counter") 
write.table(res_mcp_counter, 'mcp_counter.txt', row.names = T, sep = '\t', quote = FALSE)
res_cibersort  <- deconvolute(expr, method="cibersort") 
write.table(res_cibersort, 'cibersort.txt', row.names = T, sep = '\t', quote = FALSE)
res_TIMER  <- deconvolute(expr, method="timer",indications=20530) 
write.table(res_TIMER, 'TIMER.txt', row.names = T, sep = '\t', quote = FALSE)
# Illustrate the heatmap
# Load the row and column definition
TOTAL <- read.table("TOTAL.txt",row.names = 1,header = T)
high <- read.table("high.txt",row.names = 1,header = T)
low <- read.table("low.txt",row.names = 1,header = T)
colname <- read.table("col_name.txt",row.names = 1,header = T)
rowname <- read.table("rowname_NA3.txt",row.names = 1,header = T)
colname_high <- read.table("colname_high.txt",row.names = 1,header = T)
matrix(colname_high)
colname_low <- read.table("colname_low.txt",row.names = 1,header = T)
matrix(colname_low)

col = list(TNM_Stage = c("1" = "#66C2A5", "2" = "#0099B4B2", 
                         "3" = "#925E9FB2","4" = "#FDAF91B2"),
           TMErisk.Score = c("low" = "#91D1C2CC","high" = "#00A087CC"),
           type = c("TIMER" = "#66C2A5", "CIBERSORT" = "#FC8D62", 
                    "MCPCOUNTER" = "#8DA0CB") ,
           ESTIMATE_Score =colorRampPalette(colors = c("#8DD3C7", "#FFFFB3"))(10),
           Tumor_Purity  =colorRampPalette(colors = c("#BEBADA", "#FCCDE5"))(10),
           Stromal_Score  =colorRampPalette(colors = c("#CCEBC5", "#80B1D3"))(10),
           Immune_Score  =colorRampPalette(colors = c("#FEE0D2", "#FC9272"))(10))
p1 <- rowAnnotation(cluster = anno_block(gp = gpar(fill = c("#66C2A5","#FC8D62","#8DA0CB") ), 
                                         labels = c("TIMER", "CIBERSORT","MCPCOUNTER"), 
                                         labels_gp = gpar(cex = 0.8, col = "black"))) 

rownames(colname) <- colnames(TOTAL)
pdf(file="Immune Landscape5.pdf", width=14, height=12)
pheatmap(TOTAL,
         cluster_cols = F,
         cluster_rows = T,
         gaps_col = c(144),
         show_colnames = F, 
         fontsize = 16,
         fontsize_row = 12,
         color = colorRamp2(c(0,7.5,15), c("#F9FDFD", "red","#DF0000")),
         annotation_col  = colname,
         annotation_row = rowname,
         annotation_colors = col,border_color=NA)
dev.off()