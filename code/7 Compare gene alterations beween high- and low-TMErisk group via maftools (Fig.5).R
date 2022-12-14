library(BSgenome)
library(deconstructSigs)
library(maftools)
library(GenVisR)
library(forestplot)
library(tidyverse)
library(forestplot)
library(meta)
colnames(Data_str)
# Fig.5A: Illustrate gene mutation status heatmap
mut_high = read.csv("high.csv")
mut_low = read.csv("low.csv")
clin_high = read.csv("high_cli.csv")
clin_low = read.csv("low_cli.csv")
laml_high = read.maf(maf = mut_high, clinicalData = clin_high)
laml_low = read.maf(maf = mut_low, clinicalData = clin_low)

col1=c("#BEBADA" ,"#FCCDE5")
names(col1) = c("No", "Yes")
col2=c("#FFFFB3" ,"#FDB462")
names(col2) = c("No", "Yes")
col3=c("#FB8072", "#80B1D3")
names(col3) = c("Yes", "No")
col = list(Alcohol.consumption = col1, Hepatitis.B = col2, Hepatitis.C = col3)

oncoplot(maf = laml_high, draw_titv = T,
         genes = c("CTNNB1","TP53","AXIN1","TTN","MUC16","OBSCN","PCLO","DOCK2","RYR2","HMCN1"),
         keepGeneOrder = T,drawRowBar = F,drawColBar = TRUE,logColBar = T,
         clinicalFeatures = c("Alcohol.consumption","Hepatitis.B", "Hepatitis.C"),
         annotationColor = col,barcode_mar = 100,
         barcodeSrt = 90,gene_mar = 7,legend_height = 4,
         bgCol = "gray87",borderCol = "white",
         fontSize = 0.9,titleFontSize = 2,legendFontSize = 1.5,
         showTitle = TRUE,titleText = "",showPct = TRUE)
oncoplot(maf = laml_low, draw_titv = T,
         genes = c("CTNNB1","TP53","AXIN1","TTN","MUC16","OBSCN","PCLO","DOCK2","RYR2","HMCN1"),
         keepGeneOrder = T,drawRowBar = F,drawColBar = TRUE,logColBar = T,
         clinicalFeatures = c("Alcohol.consumption","Hepatitis.B", "Hepatitis.C"),
         annotationColor = col,barcode_mar = 100,
         barcodeSrt = 90,gene_mar = 7,legend_height = 4,
         bgCol = "gray87",borderCol = "white",
         fontSize = 0.9,titleFontSize = 2,legendFontSize = 1.5,
         showTitle = TRUE,titleText = "",showPct = TRUE)

#Fig.5B: Compare the difference in TMB 
TMB <- tmb(maf = laml,logScale = F)
TMB_compare <- ggplot(TOTAL,aes(x = as.factor(SAMPLE),y = TMB))+
  geom_boxplot(aes(fill = SAMPLE),outlier.shape = 19,outlier.colour = NA, width=0.7, cex=0.5, notch = F, notchwidth = 0.7,)+
  scale_fill_lancet(alpha = 0.8)+
  labs(fill="",y="TMB")+
  stat_compare_means(comparisons =my_comparisons,label.y = c(2200),hide.ns = TRUE,label = "p.signif")+
  theme_bw()+ 
  theme(panel.grid=element_blank(),axis.text.x = element_text(size = 12), legend.text=element_text(size=12),legend.position="top",legend.title=element_text(size=20))

#Fig.5C: Compare the hazard ratios of highly mutated genes 
Data_str <- read.csv("FOREST.csv",header = FALSE)
forestplot(labeltext = as.matrix(Data_str[,1:5]), 
           mean  = Data_str$V7,lower = Data_str$V8,upper = Data_str$V6,
           is.summary = c(rep(TRUE, 1), rep(FALSE, 10)),
           xticks = c(0,1,3),clip = c(0,2),
           fn.ci_norm = fpDrawCircleCI, 
           zero = 1, boxsize = 0.3, colgap = unit(6,'mm'),
           lwd.zero = 1.8,lty.ci =7,lwd.ci = 3,
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.6), xlab = gpar(cex = 0.8), cex = 0.8), 
           lineheight = "auto", 
           lwd.xaxis=2,xlog=F,grid = F,
           graphwidth = unit(0.35,"npc"), 
           col=fpColors(box='#FDB462',summary= "#BC80BD",lines = '#BA5F54',zero = '#D9D9D9'),
           xlab="Odds ratio with 95% CI",
           graph.pos = 1)

#Fig.5D: Correlation analyze for top 10 mutated genes in high TMErisk group
Interact <- somaticInteractions(maf = laml_low, top = 10, pvalue = c(0.05, 0.1),
                                fontSize = 0.8,showSum = FALSE, countType= 0.7,colNC = 10,colPal = "PiYG", sigSymbolsSize = 2,sigSymbolsFontSize = 1)

#Fig.5E:  Display the distribution of mutation types 
lollipopPlot2(m1 = laml_high, m2 = laml_low, 
              gene = "CTNNB1",AACol1 = "HGVSp",AACol2 = "HGVSp", 
              m1_name = "high TMErisk", m2_name = "low TMErisk",
              domainLabelSize = 1,labPosSize=5,legendTxtSize=1.1)
lollipopPlot2(m1 = laml_high, m2 = laml_low, 
              gene = "TP53",AACol1 = "HGVSp",AACol2 = "HGVSp", 
              m1_name = "high TMErisk", m2_name = "low TMErisk",
              domainLabelSize = 1,labPosSize=5,legendTxtSize=1.1)
lollipopPlot2(m1 = laml_high, m2 = laml_low, 
              gene = "AXIN1",AACol1 = "HGVSp",AACol2 = "HGVSp", 
              m1_name = "high TMErisk", m2_name = "low TMErisk",
              domainLabelSize = 1,labPosSize=5,legendTxtSize=1.1)



