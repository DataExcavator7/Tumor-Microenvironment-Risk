# Fig.3A-B left panel: Take the intersection via venn
setwd("E:/TCGA-LIHC/VENN")
library(ggvenn)
library(VennDiagram)
library(RColorBrewer)
blue<-read.table("blue.txt", header=T, sep="\t", check.names=F) 
downSig<-read.table("downSig.txt", header=T, sep="\t", check.names=F) 
x<-list("RRA"=downSig$Name,"WGCNA"=blue$MM)
pdf("Protect.pdf",width=5, height=5)
ggvenn(x,fill_color=c("LightCoral","PowDerBlue"),fill_alpha = .7,stroke_linetype = "solid",set_name_size = 8,text_size=5) 
dev.off()
modules<-read.table("3modules.txt", header=T, sep="\t", check.names=F) 
upSig<-read.table("upSig.txt", header=T, sep="\t", check.names=F) 
x<-list("RRA"=upSig$Name,"WGCNA"=modules$MM)
pdf("Dystroy.pdf",width=5, height=5)
ggvenn(x,fill_color=c("LightCoral","PowDerBlue"),fill_alpha = .7,stroke_linetype = "solid",set_name_size = 8,text_size=5) 
dev.off()

# Fig.3A-B right panel: Develop the TMErisk model via  
setwd("E:/TCGA-LIHC/Lasso+Cox")
library(dplyr)
library(car) 
library(corrplot) 
library(leaps)
library(glmnet) 
library(caret) 
#129 genes from the intersection of RRA down and MEblue
data <- read.table("237_129.txt",sep = "\t",header=T,row.names=1)
# Store the independent and dependent variables
x <- as.matrix(data[-c(1,2,3)])
y <- data[,2]
lasso <- glmnet(x, y, alpha = 1,family="binomial")
print(lasso)
par(mfrow=c(1,1))
plot(lasso, xvar = "lambda", label = TRUE)  
cvfit=cv.glmnet(x,y)
plot(cvfit)
model_lasso_min <- glmnet(x=x, y=y, alpha = 1, lambda=cvfit$lambda.min)
model_lasso_1se <- glmnet(x=x, y=y, alpha = 1, lambda=cvfit$lambda.1se)
head(model_lasso_min$beta)
cvfit$lambda.min
cvfit$lambda.1se
choose_gene_min=rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
choose_gene_1se=rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]
length(choose_gene_min)
length(choose_gene_1se)
as.matrix(model_lasso_min$beta)
as.matrix(model_lasso_1se$beta)
write.table(as.matrix(model_lasso_min$beta), file="13.txt", sep="\t", quote=F, row.names = T)

#96 genes from the intersection of RRA up and MEred, MEsalmon, and MEgrey60
data <- read.table("237_96.txt",sep = "\t",header=T,row.names=1)
# Store the independent and dependent variables
x <- as.matrix(data[-c(1,2,3)])
y <- data[,2]
lasso <- glmnet(x, y, alpha = 1,family="binomial")
print(lasso)
par(mfrow=c(1,1))
plot(lasso, xvar = "lambda", label = TRUE)  
cvfit=cv.glmnet(x,y)
plot(cvfit)
model_lasso_min <- glmnet(x=x, y=y, alpha = 1, lambda=cvfit$lambda.min)
model_lasso_1se <- glmnet(x=x, y=y, alpha = 1, lambda=cvfit$lambda.1se)
cvfit$lambda.min
cvfit$lambda.1se
head(model_lasso_min$beta)
choose_gene_min=rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
choose_gene_1se=rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]
length(choose_gene_min)
length(choose_gene_1se)
as.matrix(model_lasso_min$beta)
as.matrix(model_lasso_1se$beta)
write.table(as.matrix(model_lasso_min$beta), file="5.txt", sep="\t", quote=F, row.names = T)

# Cox regression
library(My.stepwise)
library(survival)
data <- read.table("protect_binary.txt",sep = "\t",header=T,row.names=1)
protect <- read.table("protect8.txt",sep = "\t",header=T)
vl <- protect$Gene
My.stepwise.coxph(Time = "time", 
                  Status = "status", 
                  variable.list = vl, 
                  data = data)
model<-coxph(formula = Surv(time, status) ~ ITPRIP + KLRB1 + SYNE1 + 
               CAMK4, data = data)
summary(model)$concordance
broom::glance(model)$AIC

# Fig.3C: Survival curve
dat<-read.table("clinical_information.txt",sep = "\t",header=T,row.names=1)
cut <- surv_cutpoint(dat,
                     time = "time", event = "status",
                     variables ="TMErisk score" )
summary(cut)
plot(cut, "TMErisk score", palette = "npg")
cat <- surv_categorize(cut)#Categorize variables
head(cat)
write.csv(cat,"TMErisk score.csv")
fit <- survfit(Surv(time, status) ~TMErisk score, data = cat)
pdf(file="TMErisk score.pdf", width=5, height=5,onefile=FALSE)
ggsurvplot(fit,
           data=cat,
           size = 0.9,
           xlab = "Follow up time(d)",
           pval = TRUE, 
           surv.median.line="hv", 
           palette = c("#ED0000B2","#00468BB2"),
           conf.int=FALSE,
           risk.table = FALSE)
dev.off()

# Fig.S3A-B: Examination of the risk and protective factors on prognosis 
# Here take CKAP2 as an example
library (survminer)
dat<-read.table("clinical_information.txt",sep = "\t",header=T,row.names=1)
cut <- surv_cutpoint(dat,
                     time = "time", event = "status",
                     variables ="CKAP2" )
summary(cut)
plot(cut, "CKAP2", palette = "npg")
cat <- surv_categorize(cut)
fit <- survfit(Surv(time, status) ~CKAP2, data = cat)
pdf(file="CKAP2.pdf", width=8, height=6,onefile=FALSE)
ggsurvplot(fit, 
           data=cat,
           size = 0.9,
           legend.labs = c("High", "Low"),
           legend = c(0.9,0.8), 
           xlab = "Follow up time(d)", 
           pval = TRUE, 
           surv.median.line="hv", 
           palette = c("#ED0000B2","#00468BB2"),
           conf.int=FALSE,
           #break.x.by = 100)  
           risk.table = FALSE)
dev.off()

# Fig.S3C-D: Estimate hazard ratios (HR) of acquired genes
rt <- read.table("together_binary.txt",sep = "\t",header=T,row.names=1)
rt$time=rt$time/365     
pFilter=0.05
outTab=data.frame()
sigGenes=c("time","status")
for(gene in colnames(rt[,3:ncol(rt)])){
  cox=coxph(Surv(time, status) ~ rt[,gene], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  
  if(coxP<pFilter){
    sigGenes=c(sigGenes,gene)
    outTab=rbind(outTab,
                 cbind(gene=gene,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxP) )
  }
}
write.table(outTab,file="uniCox18genes.txt",sep="\t",row.names=F,quote=F)
surSigExp=rt[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="uniSigExp.txt",sep="\t",row.names=F,quote=F)

# Fig.3D: Calculate correlation between the TMErisk score and ESTIMATE score
setwd("E:/TCGA-LIHC/Scatter Plot")
library(ggplot2)
library(ggpubr)
data<-read.table("Scatter Plot.txt",sep = "\t",header=T,row.names=1)
pdf(file="TMErisk score Box.pdf", width=5, height=6)
ggboxplot(data,
          x = "TMEriskGroup",
          y = "TMErisk",
          color = "black",
          fill = "TMEriskGroup",
          xlab = "",
          palette = c("#ED0000B2","#00468BB2"),
          ylab = "ESTIMATE score",
          main = "" +
            theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1),
                  legend.position="top",legend.title=element_text(size=20))+
            stat_compare_means(aes(group = TMEriskGroup,label = "kruskal.test"),
                               label = "p.signif", 
                               na.rm = TRUE)+ylim(-2.6,0.1))+geom_signif(comparisons = list(c("low","high")),map_signif_level = TRUE,test = t.test,y_position = c(-0.25,0),tip_length = c(0.05,0.05))
dev.off()

pdf(file="ESTIMATE score Scatter.pdf", width=8, height=6)
ggplot(data)+geom_smooth(data, mapping=aes(x=TMErisk,y=ESTIMATEScore),method='lm',se=T,size=1.25,color ="black",linetype="longdash")+geom_point(data, mapping=aes(x=TMErisk,y=ESTIMATEScore,color=TMEriskGroup), alpha=.9)+
  theme_bw() +scale_color_manual(values = c("#ED0000B2","#00468BB2"))+ylim(-3500,3500)+xlim(-2.6,0.1)+theme(axis.title.x =element_text(size=20), axis.title.y=element_text(size=20),axis.text = element_text(size = 12),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"),plot.title = element_text(hjust = 0.5,size = 20,face = "bold")) 
dev.off()
data$TMEriskGroup <- factor(data$TMEriskGroup,levels=c("low","high"),ordered = TRUE)
pdf(file="ESTIMATE score Box.pdf", width=3, height=6)
ggboxplot(data,
          x = "TMEriskGroup",
          y = "ESTIMATEScore",
          color = "black",
          fill = "TMEriskGroup",
          xlab = "",
          palette = c( "#00468BB2","#ED0000B2"),
          ylab = "ESTIMATE score",
          main = "" +
            theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1),
                  legend.position="top",legend.title=element_text(size=20))+
            stat_compare_means(aes(group = TMEriskGroup,label = "kruskal.test"),
                               label = "p.signif", 
                               na.rm = TRUE)+ylim(-3500,3500))+geom_signif(comparisons = list(c("low","high")),map_signif_level = TRUE,test = t.test,y_position = c(3000,0),tip_length = c(0.05,0.05))
dev.off()


pdf(file="Stromal score Scatter.pdf", width=8, height=6)
ggplot(data)+geom_smooth(data, mapping=aes(x=TMErisk,y=StromalScore),method='lm',se=T,size=1.25,color ="black",linetype="longdash")+geom_point(data, mapping=aes(x=TMErisk,y=StromalScore,color=TMEriskGroup), alpha=.9)+
  theme_bw() +scale_color_manual(values = c("#ED0000B2","#00468BB2"))+ylim(-2000,2000)+xlim(-2.6,0.1)+theme(axis.title.x =element_text(size=20), axis.title.y=element_text(size=20),axis.text = element_text(size = 12),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"),plot.title = element_text(hjust = 0.5,size = 20,face = "bold")) 
dev.off()
data$TMEriskGroup <- factor(data$TMEriskGroup,levels=c("low","high"),ordered = TRUE)
pdf(file="Stromal score Box.pdf", width=3, height=6)
ggboxplot(data,
          x = "TMEriskGroup",
          y = "StromalScore",
          color = "black",
          fill = "TMEriskGroup",
          xlab = "",
          palette = c( "#00468BB2","#ED0000B2"),
          ylab = "Stromal score",
          main = "" +
            theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1),
                  legend.position="top",legend.title=element_text(size=20))+
            stat_compare_means(aes(group = TMEriskGroup,label = "kruskal.test"),
                               label = "p.signif", 
                               na.rm = TRUE)+ylim(-2000,2000))+geom_signif(comparisons = list(c("low","high")),map_signif_level = TRUE,test = t.test,y_position = c(1500,0),tip_length = c(0.05,0.05))
dev.off() 

pdf(file="Immune score Scatter.pdf", width=8, height=6)
ggplot(data)+geom_smooth(data, mapping=aes(x=TMErisk,y=ImmuneScore),method='lm',se=T,size=1.25,color ="black",linetype="longdash")+geom_point(data, mapping=aes(x=TMErisk,y=ImmuneScore,color=TMEriskGroup), alpha=.9)+
  theme_bw() +scale_color_manual(values = c("#ED0000B2","#00468BB2"))+ylim(-1500,3000)+xlim(-2.6,0.1)+theme(axis.title.x =element_text(size=20), axis.title.y=element_text(size=20),axis.text = element_text(size = 12),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"),plot.title = element_text(hjust = 0.5,size = 20,face = "bold")) 
dev.off()
data$TMEriskGroup <- factor(data$TMEriskGroup,levels=c("low","high"),ordered = TRUE)
pdf(file="Immune score Box.pdf", width=3, height=6)
ggboxplot(data,
          x = "TMEriskGroup",
          y = "ImmuneScore",
          color = "black",
          fill = "TMEriskGroup",
          xlab = "",
          palette = c( "#00468BB2","#ED0000B2"),
          ylab = "Immune score",
          main = "" +
            theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1),
                  legend.position="top",legend.title=element_text(size=20))+
            stat_compare_means(aes(group = TMEriskGroup,label = "kruskal.test"),
                               label = "p.signif", 
                               na.rm = TRUE)+ylim(-1500,3000))+geom_signif(comparisons = list(c("low","high")),map_signif_level = TRUE,test = t.test,y_position = c(2500,0),tip_length = c(0.05,0.05))
dev.off() 

pdf(file="TumorPurity score Scatter.pdf", width=8, height=6)
ggplot(data)+geom_smooth(data, mapping=aes(x=TMErisk,y=TumorPurity),method='lm',se=T,size=1.25,color ="black",linetype="longdash")+geom_point(data, mapping=aes(x=TMErisk,y=TumorPurity,color=TMEriskGroup), alpha=.9)+
  theme_bw() +scale_color_manual(values = c("#ED0000B2","#00468BB2"))+ylim(0,1.2)+xlim(-2.6,0.1)+theme(axis.title.x =element_text(size=20), axis.title.y=element_text(size=20),axis.text = element_text(size = 12),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"),plot.title = element_text(hjust = 0.5,size = 20,face = "bold")) 
dev.off()
data$TMEriskGroup <- factor(data$TMEriskGroup,levels=c("low","high"),ordered = TRUE)
pdf(file="TumorPurity score Box.pdf", width=3, height=6)
ggboxplot(data,
          x = "TMEriskGroup",
          y = "TumorPurity",
          color = "black",
          fill = "TMEriskGroup",
          xlab = "",
          palette = c( "#00468BB2","#ED0000B2"),
          ylab = "TumorPurity",
          main = "" +
            theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1),
                  legend.position="top",legend.title=element_text(size=20))+
            stat_compare_means(aes(group = TMEriskGroup,label = "kruskal.test"),
                               label = "p.signif", 
                               na.rm = TRUE)+ylim(0,1.2))+geom_signif(comparisons = list(c("low","high")),map_signif_level = TRUE,test = t.test,y_position = c(1,0),tip_length = c(0.05,0.05))
dev.off() 

# Fig.3E: Calculate time dependent ROC of 6 datasets (TCGA-LIHC, ICGC-JP, GSE76427, GSE144269, GSE116174, and GSE10186)
setwd("E:/TCGA-LIHC/Lasso+Cox/COX/ROC")
library(timeROC)
library(survival)
dat<-read.table("data_TCGA.txt",sep = "\t",header=T,row.names=1)
time_roc_res <- timeROC(
  T = dat$time,
  delta = dat$status,
  marker = dat$TCGA,
  cause = 1,
  weighting="marginal",
  times = c(0.5 * 365,1 * 365,1.5 * 365,2 * 365,2.5 * 365,3 * 365,3.5 * 365, 4 * 365,4.5 * 365,5 * 365),
  ROC = TRUE,
  iid = TRUE)
time_roc_res$AUC
confint(time_roc_res, level = 0.95)$CI_AUC
plot(time_roc_res, time=2 * 365, col = "red", title = FALSE)  
plot(time_roc_res, time=4 * 365, add=TRUE, col="blue") 
plot(time_roc_res, time=6 * 365, add=TRUE, col="green") 
legend("bottomright",c("2 Years" ,"4 Years", "6 Years"),
       col=c("red", "blue", "green"), lty=1, lwd=2)
time_ROC_df <- data.frame(
  TP_2year = time_roc_res$TP[, 1],
  FP_2year = time_roc_res$FP[, 1],
  TP_4year = time_roc_res$TP[, 2],
  FP_4year = time_roc_res$FP[, 2],
  TP_6year = time_roc_res$TP[, 3],
  FP_6year = time_roc_res$FP[, 3]
)
library(ggplot2)
ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_2year, y = TP_2year), size = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_4year, y = TP_4year), size = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_6year, y = TP_6year), size = 1, color = "#E18727FF") +
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 2 years = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#BC3C29FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 4 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#0072B5FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 6 years = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#E18727FF"
  ) +
  labs(x = "False positive rate", y = "True positive rate") +
  theme(
    axis.text = element_text(face = "bold", size = 11, color = "black"),
    axis.title.x = element_text(face = "bold", size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
    axis.title.y = element_text(face = "bold", size = 14, color = "black", margin = margin(c(0, 15, 0, 0)))
  )
dat_GSE76247<-read.table("data_GSE76247.txt",sep = "\t",header=T,row.names=1)
time_roc_res_GSE76247 <- timeROC(
  T = dat_GSE76247$time,
  delta = dat_GSE76247$status,
  marker = dat_GSE76247$GSE76247,
  cause = 1,
  weighting="marginal",
  times = c(0.5 * 365,1 * 365,1.5 * 365,2 * 365,2.5 * 365,3 * 365,3.5 * 365, 4 * 365,4.5 * 365,5 * 365),
  ROC = TRUE,
  iid = TRUE)
time_roc_res_GSE76247$AUC
plotAUCcurve(time_roc_res, conf.int=F, col="red")
plotAUCcurve(time_roc_res_GSE76247, conf.int=F, col="blue", add=TRUE)
legend("bottomright",c("TCGA", "GSE76247"), col = c("red","blue"), lty=1, lwd=2)

dat_GSE116174<-read.table("data_GSE116174.txt",sep = "\t",header=T,row.names=1)
time_roc_res_GSE116174 <- timeROC(
  T = dat_GSE116174$time,
  delta = dat_GSE116174$status,
  marker = dat_GSE116174$GSE116174,
  cause = 1,
  weighting="marginal",
  times = c(0.5 * 365,1 * 365,1.5 * 365,2 * 365,2.5 * 365,3 * 365,3.5 * 365, 4 * 365,4.5 * 365,5 * 365),
  ROC = TRUE,
  iid = TRUE)
time_roc_res_GSE116174$AUC

dat_GSE144269<-read.table("data_GSE144269.txt",sep = "\t",header=T,row.names=1)
time_roc_res_GSE144269 <- timeROC(
  T = dat_GSE144269$time,
  delta = dat_GSE144269$status,
  marker = dat_GSE144269$GSE144269,
  cause = 1,
  weighting="marginal",
  times = c(0.5 * 365,1 * 365,1.5 * 365,2 * 365,2.5 * 365,3 * 365,3.5 * 365, 4 * 365),
  ROC = TRUE,
  iid = TRUE)
time_roc_res_GSE144269$AUC

dat_GSE10186<-read.table("data_GSE10186.txt",sep = "\t",header=T,row.names=1)
time_roc_res_GSE10186 <- timeROC(
  T = dat_GSE10186$time,
  delta = dat_GSE10186$status,
  marker = dat_GSE10186$GSE10186,
  cause = 1,
  weighting="marginal",
  times = c(0.75 * 365,1 * 365,1.5 * 365,2 * 365,2.5 * 365,3 * 365,3.5 * 365, 4 * 365,4.5 * 365,5 * 365),
  ROC = TRUE,
  iid = TRUE)
time_roc_res_GSE10186$AUC

dat_ICGCJP<-read.table("data_ICGCJP.txt",sep = "\t",header=T,row.names=1)
time_roc_res_ICGCJP <- timeROC(
  T = dat_ICGCJP$time,
  delta = dat_ICGCJP$status,
  marker = dat_ICGCJP$ICGC,
  cause = 1,
  weighting="marginal",
  times = c(0.5 * 365,1 * 365,1.5 * 365,2 * 365,2.5 * 365,3 * 365,3.5 * 365, 4 * 365,4.5 * 365,5 * 365),
  ROC = TRUE,
  iid = TRUE)
time_roc_res_ICGCJP$AUC

write.table(time_roc_res$AUC,file="TCGAAUC.txt",sep="\t",quote=F,row.names=T)
write.table(time_roc_res_GSE76247$AUC,file="GSE76247AUC.txt",sep="\t",quote=F,row.names=T)
write.table(time_roc_res_GSE116174$AUC,file="GSE116174AUC.txt",sep="\t",quote=F,row.names=T)
write.table(time_roc_res_GSE144269$AUC,file="GSE144269AUC.txt",sep="\t",quote=F,row.names=T)
write.table(time_roc_res_GSE10186$AUC,file="GSE10186AUC.txt",sep="\t",quote=F,row.names=T)
write.table(time_roc_res_ICGCJP$AUC,file="ICGCJPAUC.txt",sep="\t",quote=F,row.names=T)

# Fig.3F: Estimate hazard ratios (HR) of 5 datasets
rt <- read.csv("TME_ICGCJP.csv",row.names = 1)
outTab=data.frame()
sigGenes=c("time","status")
cox=coxph(Surv(time, status) ~ TMErisk score, data = rt)
coxSummary = summary(cox)
coxP=coxSummary$coefficients[,"Pr(>|z|)"]
HR=coxSummary$conf.int[,"exp(coef)"]
HR.95L=coxSummary$conf.int[,"lower .95"]
HR.95H=coxSummary$conf.int[,"upper .95"]
pvalue=coxP
outTab=rbind(cbind(HR=coxSummary$conf.int[,"exp(coef)"],
                   HR.95L=coxSummary$conf.int[,"lower .95"],
                   HR.95H=coxSummary$conf.int[,"upper .95"],
                   pvalue=coxP) )
write.table(outTab,file="ICGCJP.txt",sep="\t",row.names=F,quote=F)

rt <- read.csv("TME_TCGA.csv",row.names = 1)
outTab=data.frame()
sigGenes=c("time","status")
cox=coxph(Surv(time, status) ~ TMErisk score, data = rt)
coxSummary = summary(cox)
coxP=coxSummary$coefficients[,"Pr(>|z|)"]
HR=coxSummary$conf.int[,"exp(coef)"]
HR.95L=coxSummary$conf.int[,"lower .95"]
HR.95H=coxSummary$conf.int[,"upper .95"]
pvalue=coxP
outTab=rbind(cbind(HR=coxSummary$conf.int[,"exp(coef)"],
                   HR.95L=coxSummary$conf.int[,"lower .95"],
                   HR.95H=coxSummary$conf.int[,"upper .95"],
                   pvalue=coxP) )
write.table(outTab,file="TCGA.txt",sep="\t",row.names=F,quote=F)

rt <- read.csv("TME_GSE144269.csv",row.names = 1)
outTab=data.frame()
sigGenes=c("time","status")
cox=coxph(Surv(time, status) ~ TMErisk score, data = rt)
coxSummary = summary(cox)
coxP=coxSummary$coefficients[,"Pr(>|z|)"]
HR=coxSummary$conf.int[,"exp(coef)"]
HR.95L=coxSummary$conf.int[,"lower .95"]
HR.95H=coxSummary$conf.int[,"upper .95"]
pvalue=coxP
outTab=rbind(cbind(HR=coxSummary$conf.int[,"exp(coef)"],
                   HR.95L=coxSummary$conf.int[,"lower .95"],
                   HR.95H=coxSummary$conf.int[,"upper .95"],
                   pvalue=coxP) )
write.table(outTab,file="GSE144269.txt",sep="\t",row.names=F,quote=F)

rt <- read.csv("TME_GSE116174.csv",row.names = 1)
outTab=data.frame()
sigGenes=c("time","status")
cox=coxph(Surv(time, status) ~ TMErisk score, data = rt)
coxSummary = summary(cox)
coxP=coxSummary$coefficients[,"Pr(>|z|)"]
HR=coxSummary$conf.int[,"exp(coef)"]
HR.95L=coxSummary$conf.int[,"lower .95"]
HR.95H=coxSummary$conf.int[,"upper .95"]
pvalue=coxP
outTab=rbind(cbind(HR=coxSummary$conf.int[,"exp(coef)"],
                   HR.95L=coxSummary$conf.int[,"lower .95"],
                   HR.95H=coxSummary$conf.int[,"upper .95"],
                   pvalue=coxP) )
write.table(outTab,file="GSE116174.txt",sep="\t",row.names=F,quote=F)

rt <- read.csv("TME_GSE10186.csv",row.names = 1)
outTab=data.frame()
sigGenes=c("time","status")
cox=coxph(Surv(time, status) ~ TMErisk score, data = rt)
coxSummary = summary(cox)
coxP=coxSummary$coefficients[,"Pr(>|z|)"]
HR=coxSummary$conf.int[,"exp(coef)"]
HR.95L=coxSummary$conf.int[,"lower .95"]
HR.95H=coxSummary$conf.int[,"upper .95"]
pvalue=coxP
outTab=rbind(cbind(HR=coxSummary$conf.int[,"exp(coef)"],
                   HR.95L=coxSummary$conf.int[,"lower .95"],
                   HR.95H=coxSummary$conf.int[,"upper .95"],
                   pvalue=coxP) )
write.table(outTab,file="GSE10186.txt",sep="\t",row.names=F,quote=F)