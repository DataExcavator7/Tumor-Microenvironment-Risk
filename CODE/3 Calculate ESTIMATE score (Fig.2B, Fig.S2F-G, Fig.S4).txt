setwd("E:/TCGA-Immune/Estimate")
library(ggplot2)
library(ggpubr)
library(tidyr)
library (survminer)
library(survival)
library(ggsci)
library(estimate)
# Calculate ESTIMATE score of the TCGA-LIHC 
whole<-read.table("HiSeqV2.txt",header=T,sep="\t",row.names = 1)
filterCommonGenes(input.f="whole.txt", output.f="whole.gct", id="GeneSymbol")
estimateScore(input.ds="whole.gct", output.ds="whole.gct")

# Boxplot
ESTIMATE<-read.table("ESTIMATE.txt",header=T,sep="\t",row.names = 1)
ggplot(TOTAL,aes(x = as.factor(SAMPLE),y = ESTIMATEscore))+
  geom_boxplot(aes(fill = SAMPLE),outlier.shape = 19,outlier.colour = NA, width=0.7, cex=0.5, notch = F, notchwidth = 0.7,)+
  scale_fill_lancet(alpha = 0.8)+
  labs(fill="",y="ESTITATE score")+
  stat_compare_means(comparisons =my_comparisons,label.y = c(2200),hide.ns = TRUE,label = "p.signif")+
  theme_bw() + 
  theme(panel.grid=element_blank(),axis.text.x = element_text(size = 12),legend.text=element_text(size=12),legend.position="top",legend.title=element_text(size=20))

# Survival curve
cut <- surv_cutpoint(dat, time = "time", event = "status", variables = "TME")
cat <- surv_categorize(cut)#Categorize variables
fit <- survfit(Surv(time, status) ~ TME, data = cat)
ggsurvplot(fit, data=cat, size = 0.9, title = "", legend.labs = c("High", "Low"), legend = c(0.9,0.8), legend.title ="", xlab = "Follow up time(d)", pval = TRUE, pval.size=5, surv.median.line="hv", palette = c( "#ED0000B2","#00468BB2"))

# The same method to illustrate the boxplot and survival curve for 6 validation datasets mentioned follow
# Randomized half samples of TCGA-LIHC, FigS2F-G
# ICGC LICA-JP, FigS4A-B
# GSE144269, FigS4C-D
# GSE116174, FigS4E-F
# GSE76427, FigS4G-H
# GSE10186, FigS4I-J