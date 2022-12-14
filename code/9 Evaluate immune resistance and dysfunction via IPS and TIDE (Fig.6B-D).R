# Evaluate IPS of TCGA-LIHC 
gene_expression<-read.table(file = "TPM_expr.txt", sep = "\t",header = T,row.names = 1)
library(ggplot2)
library(grid)
library(gridExtra)
# Calculate Immunophenoscore
ipsmap<- function (x) {
  if (x<=0) {
    ips<-0
  } else {
    if (x>=3) {
      ips<-10
    } else {
      ips<-round(x*10/3, digits=0)
    }
  }
  return(ips)
}

# Assign colors 
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)
mapcolors<-function (x) {
  za<-NULL
  if (x>=3) {
    za=1000
  } else {
    if (x<=-3) {
      za=1
    } else {
      za=round(166.5*x+500.5,digits=0)
    }
  }
  return(my_palette[za])
}
my_palette2 <- colorRampPalette(c("black", "white"))(n = 1000)
mapbw<-function (x) {
  za2<-NULL
  if (x>=2) {
    za2=1000
  } else {
    if (x<=-2) {
      za2=1
    } else {
      za2=round(249.75*x+500.5,digits=0)
    }
  }
  return(my_palette2[za2])
}
# Read expression data from tab-delimited text file, with official human gene symbols (HGNC) in the first columns and expression values (i.e. log2(TPM+1)) for each sample in the other columns
sample_names<-colnames(gene_expression)
# Read IPS genes and corresponding weights from tab-delimited text file "IPS_genes.txt"
IPSG<-read.table("IPS_genes.txt",header=TRUE, sep="\t", dec = ".",check.names=FALSE)
unique_ips_genes<-as.vector(unique(IPSG$NAME))

IPS<-NULL
MHC<-NULL
CP<-NULL
EC<-NULL
SC<-NULL
AZ<-NULL

# Gene names in expression file
GVEC<-rownames(gene_expression)
# Genes names in IPS genes file
VEC<-IPSG$GENE
# Match IPS genes with genes in expression file
ind<-which(is.na(match(VEC,GVEC)))
# List genes missing or differently named
MISSING_GENES<-VEC[ind]
dat<-IPSG[ind,]
if (length(MISSING_GENES)>0) {
  cat("differently named or missing genes: ",MISSING_GENES,"\n")
}
for (x in 1:length(ind)) {
  print(IPSG[ind,])
}

for (i in 1:length(sample_names)) {	
  GE<-gene_expression[[i]]
  mGE<-mean(GE)
  sGE<-sd(GE)
  Z1<-(gene_expression[as.vector(IPSG$GENE),i]-mGE)/sGE
  W1<-IPSG$WEIGHT
  WEIGHT<-NULL
  MIG<-NULL
  k<-1
  for (gen in unique_ips_genes) {
    MIG[k]<- mean(Z1[which (as.vector(IPSG$NAME)==gen)],na.rm=TRUE)
    WEIGHT[k]<- mean(W1[which (as.vector(IPSG$NAME)==gen)])
    k<-k+1
  }
  WG<-MIG*WEIGHT
  MHC[i]<-mean(WG[1:10])
  CP[i]<-mean(WG[11:20])
  EC[i]<-mean(WG[21:24])
  SC[i]<-mean(WG[25:26])
  AZ[i]<-sum(MHC[i],CP[i],EC[i],SC[i])
  IPS[i]<-ipsmap(AZ[i])
  
  ## Plot Immunophenogram
  data_a <- data.frame (start = c(0,2.5,5,7.5,10,15,seq(20,39),0,10,20,30), end = c(2.5,5,7.5,10,15,seq(20,40),10,20,30,40), y1=c(rep(2.6,26),rep(0.4,4)),y2=c(rep(5.6,26),rep(2.2,4)),z=c(MIG[c(21:26,11:20,1:10)],EC[i],SC[i],CP[i],MHC[i]),vcol=c(unlist(lapply(MIG[c(21:26,11:20,1:10)],mapcolors)), unlist(lapply(c(EC[i],SC[i],CP[i],MHC[i]),mapbw))), label = c(unique_ips_genes[c(21:26,11:20,1:10)],"EC","SC","CP","MHC"))
  data_a$label <- factor(data_a$label, levels=unique(data_a$label))
  plot_a1<-ggplot() + geom_rect(data=data_a, mapping=aes(xmin=start, xmax=end, ymin=y1, ymax=y2, fill=label), size=0.5,color="black", alpha=1) +  coord_polar() + scale_y_continuous(limits = c(0, 6)) + scale_fill_manual(values =as.vector(data_a$vcol),guide=FALSE) + theme_bw() + theme(panel.margin = unit(0, 'mm'), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "white"), axis.text=element_blank(), axis.ticks= element_blank()) + geom_text(aes(x=5, y=1.3, label="EC"), size=4) + geom_text(aes(x=15, y=1.3, label="SC"), size=4) + geom_text(aes(x=25, y=1.3, label="CP"), size=4) + geom_text(aes(x=35, y=1.3, label="MHC"), size=4)
  plot_a2<-plot_a1+geom_text(aes(x=1.25, y=4.1, label="+ Act CD4"), angle=78.75, size=4)+geom_text(aes(x=3.75, y=4.1, label="+ Act CD8"),angle=56.25, size=4)+geom_text(aes(x=6.25, y=4.1, label="+ Tem CD4"), angle=33.75,size=4)+geom_text(aes(x=8.75, y=4.1, label="+ Tem CD8"), angle=11.25,size=4)+geom_text(aes(x=17.5, y=4.1, label="- MDSC"), angle=-67.5,size=4)+geom_text(aes(x=12.5, y=4.1, label="- Treg"), angle=-22.5,size=4)	
  plot_a3<-plot_a2+geom_text(aes(x=20.5, y=4.1, label="PD-1 -"), angle=85.5, size=4)+geom_text(aes(x=21.5, y=4.1, label="CTLA4 -"), angle=76.5, size=4)+geom_text(aes(x=22.5, y=4.1, label="LAG3 -"), angle=67.5, size=4)+geom_text(aes(x=23.5, y=4.1, label="TIGIT -"), angle=58.5, size=4)+geom_text(aes(x=24.5, y=4.1, label="TIM3 -"), angle=49.5, size=4)+geom_text(aes(x=25.5, y=4.1, label="PD-L1 -"), angle=40.5, size=4)+geom_text(aes(x=26.5, y=4.1, label="PD-L2 -"), angle=31.5, size=4)+geom_text(aes(x=27.5, y=4.1, label="CD27 +"), angle=22.5, size=4)+geom_text(aes(x=28.5, y=4.1, label="ICOS +"), angle=13.5, size=4)+geom_text(aes(x=29.5, y=4.1, label="IDO1 -"), angle=4.5, size=4)
  plot_a4<-plot_a3+geom_text(aes(x=30.5, y=4.1, label="B2M +"), angle=-4.5, size=4)+geom_text(aes(x=31.5, y=4.1, label="TAP1 +"), angle=-13.5, size=4)+geom_text(aes(x=32.5, y=4.1, label="TAP2 +"), angle=-22.5, size=4)+geom_text(aes(x=33.5, y=4.1, label="HLA-A +"), angle=-31.5, size=4)+geom_text(aes(x=34.5, y=4.1, label="HLA-B +"), angle=-40.5, size=4)+geom_text(aes(x=35.5, y=4.1, label="HLA-C +"), angle=-49.5, size=4)+geom_text(aes(x=36.5, y=4.1, label="HLA-DPA1 +"), angle=-58.5, size=4)+geom_text(aes(x=37.5, y=4.1, label="HLA-DPB1 +"), angle=-67.5, size=4)+geom_text(aes(x=38.5, y=4.1, label="HLA-E +"), angle=-76.5, size=4)+geom_text(aes(x=39.5, y=4.1, label="HLA-F +"), angle=-85.5, size=4)
  plot_a5<-plot_a4+geom_text(aes(x=0, y=6, label=paste("Immunophenoscore: ",IPS[i],sep="")), angle=0,size=6,vjust=-0.5)+ theme(axis.title=element_blank())
  plot_a <-plot_a5 + theme(plot.margin=unit(c(0,0,0,0),"mm")) + geom_text(vjust=1.15,hjust=0,aes(x=25.5, y=6,label="\n\n\n\n   MHC: Antigen Processing                                 EC: Effector Cells\n   CP: Checkpoints | Immunomodulators              SC: Suppressor Cells\n\n", hjust = 0), size=4)
  
  ## Legend sample-wise (averaged) z-scores
  data_b <- data.frame (start = rep(0,23), end = rep(0.7,23), y1=seq(0,22,by=1), y2=seq(1,23,by=1),z=seq(-3,3,by=6/22),vcol=c(unlist(lapply(seq(-3,3,by=6/22),mapcolors))), label = LETTERS[1:23])
  data_b_ticks <- data.frame(x = rep(1.2, 7), value = seq(-3,3, by=1), y = seq(0,6, by=1)*(22/6) +0.5)
  legendtheme <- theme(plot.margin = unit(c(2,0,2,0),"inch"), panel.margin = unit(0,"null"), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "white"), axis.text=element_blank(), axis.ticks= element_blank(), axis.title.x=element_blank())
  plot_b<-ggplot(hjust=0) + geom_rect(data=data_b, mapping=aes(xmin=start, xmax=end, ymin=y1, ymax=y2, fill=label), size=0.5,color="black", alpha=1) + scale_x_continuous(limits = c(0, 1.5),expand = c(0,0)) + scale_fill_manual(values =as.vector(data_b$vcol),guide=FALSE) + geom_text(data=data_b_ticks, aes(x=x, y=y, label=value),hjust="inward", size=4) +theme_bw() + legendtheme + ylab("Sample-wise (averaged) z-score") 
  
  ## Legend weighted z-scores
  data_c <- data.frame (start = rep(0,23), end = rep(0.7,23), y1=seq(0,22,by=1), y2=seq(1,23,by=1),z=seq(-2,2,by=4/22),vcol=c(unlist(lapply(seq(-2,2,by=4/22),mapbw))), label = LETTERS[1:23])
  data_c_ticks <- data.frame(x = rep(1.2, 5), value = seq(-2,2, by=1), y = seq(0,4, by=1)*(22/4) +0.5)
  plot_c<-ggplot() + geom_rect(data=data_c, mapping=aes(xmin=start, xmax=end, ymin=y1, ymax=y2, fill=label), size=0.5,color="black", alpha=1) + scale_x_continuous(limits = c(0, 1.5),expand = c(0,0)) + scale_fill_manual(values =as.vector(data_c$vcol),guide=FALSE) + geom_text(data=data_c_ticks, aes(x=x, y=y, label=value),hjust="inward", size=4) +theme_bw() + legendtheme + ylab("Weighted z-score") 
  
  ## Save plot to file (1 pdf file for each sample)
  file_name<-paste("IPS_",sample_names[i],".pdf",sep="")
  pdf(file_name, width=10, height=8)
  grid.arrange(plot_a,plot_b,plot_c, ncol=3, widths=c(0.8,0.1,0.1))
  dev.off()
}
DF<-data.frame(SAMPLE=sample_names,MHC=MHC,EC=EC,SC=SC,CP=CP,AZ=AZ,IPS=IPS)
write.table(DF,file="IPS.txt",row.names=FALSE, quote=FALSE,sep="\t")
barplot(DF[[7]],DF[[1]],col=c(rep("yellow",7),rep("green",4)))

# Evaluate TIDE score of TCGA-LIHC, GSE76427, GSE10186, and GSE144269
setwd("E:/TCGA-LIHC/TCGA-Immune/TIDE")
library(survminer)
library(survival)
expression<-read.table(file = "whole_287.txt", sep = "\t",header = T,row.names = 1)
CLINICAL<-read.table(file = "CLINICAL.txt", sep = "\t",header = T,row.names = 1)
SAMPLE<-read.table(file = "SAMPLE.txt", sep = "\t",header = T)
expression<-expression[,SAMPLE$Sample]
result<-read.table("exclusionsignature.txt", sep='\t', header=T) 
expression<-as.matrix(expression)
CTL_genes = c('CD8A', 'CD8B', 'GZMA', 'GZMB', 'PRF1')
writemat = function(result, output)
{
  result = result[!is.na(result[,"p"]),,drop=F]
  FDR = p.adjust(result[,"p"], method="fdr")
  result = cbind(result, FDR)
  write.table(result, file=output, sep='\t', quote=F)
}
mat = t(expression) 
pivot<-t(expression [CTL_genes,])
survival = data.frame(CLINICAL)
survival = survival[survival[,1] > 0,]
output = "LIHC"
common = Reduce(intersect, list(rownames(mat),rownames(survival),rownames(pivot)))
print(paste(length(common), "samples"))
if(length(common) < 20) stop("two few samples")
not_included = setdiff(CTL_genes, colnames(pivot))
if(length(not_included) > 0) stop(paste(c("pivot genes", not_included,"are not included"), ' '))
pivot = rowMeans(pivot[common,CTL_genes])
mat = mat[common,,drop=F]
survival = survival[common,,drop=F]
death_rate = sum(survival[,2])/dim(survival)[1]
if(death_rate < 0.1) stop(paste("death rate", death_rate, "is too low"))

surv = Surv(survival[,1], survival[,2])

if(dim(survival)[2] > 2){
  B = survival[,3:dim(survival)[2], drop=F]
}else{
  B = survival[,c(), drop=F]
}

B = cbind(B, pivot, pivot, pivot)
B = as.data.frame(B)
N_B = dim(B)[2]
coxph.pivot = summary(coxph(surv~., data=B[,c(-N_B, -(N_B-1)), drop=F]))$coef
write.csv(coxph.pivot,file=paste0(output, "_coxph.csv"))

colnames(B)[N_B-1] = "partner"
colnames(B)[N_B] = "Interaction"

features = colnames(mat)
N = length(features)
result_interaction = matrix(nrow=N, ncol=2, dimnames=list(features, c('z','p')))
result_main = matrix(nrow=N, ncol=2, dimnames=list(features, c('z','p')))
result_partial = matrix(nrow=N, ncol=2, dimnames=list(features, c('z','p')))
result_base = matrix(nrow=N, ncol=2, dimnames=list(features, c('z','p')))


for (i in 1:N){
  title = features[i]
  partner = mat[,i]
  B[,N_B-1] = partner
  B[,N_B] = partner * pivot 
  errflag = F
  coxph.fit = tryCatch(coxph(surv~., data=B),
                       error = function(e) errflag <<- T,
                       warning = function(w) errflag <<- T)
  
  if(!errflag){
    reg.summary = summary(coxph.fit)$coef
    result_interaction[i,] = reg.summary["Interaction", c("z", "Pr(>|z|)")]
    result_main[i,] = reg.summary["partner", c("z", "Pr(>|z|)")]
  }
  errflag = F
  coxph.fit = tryCatch(coxph(surv~., data=B[,-N_B]),
                       error = function(e) errflag <<- T,
                       warning = function(w) errflag <<- T)
  
  if(!errflag){
    reg.summary = summary(coxph.fit)$coef
    result_partial[i,] = reg.summary["partner", c("z", "Pr(>|z|)")]
  }
  errflag = F
  coxph.fit = tryCatch(coxph(surv~., data=B[,c(-N_B, -(N_B-2)), drop=F]),
                       error = function(e) errflag <<- T,
                       warning = function(w) errflag <<- T)
  
  if(!errflag){
    reg.summary = summary(coxph.fit)$coef
    result_base[i,] = reg.summary["partner", c("z", "Pr(>|z|)")]
  }
}

writemat(result_interaction, paste0(output, "_interaction.txt"))
writemat(result_main, paste0(output, "_main.txt"))
writemat(result_partial, paste0(output, "_partial.txt"))
writemat(result_base, paste0(output, "_base.txt"))

result_interaction<-data.frame(result_interaction)
result_diff<-result_interaction[which(result_interaction$p<0.05),]

dys_expr<-expression[row.names(result_diff),]
TIDE<-data.frame(matrix(0,237,2))
TIDE$sample<-colnames(dys_expr)
colnames(TIDE)<-c("Dysfunction","Exclusion","Sample")
for(i in 1:237){
  TIDE[i,1]<-cor(dys_expr[,i],result_diff[,1],method="pearson")
}

deg<-intersect(row.names(result),row.names(expression))
excl_ave<-result[deg,] 
excl_expr<- expression[deg,] 
for(i in 1:237){  
  TIDE[i,2]<-cor(excl_expr[,i],excl_ave[,4],method="pearson")
}

write.table(TIDE, 'TIDE.txt', row.names =T, sep = '\t', quote = FALSE)
write.table(pivot, 'pivot.txt', row.names =T, sep = '\t', quote = FALSE)

Tide_Merged<-read.table(file = "TIDE_Merged.txt", sep = "\t",header = T,row.names = 1)
Tide_Merged$scale<-scale(Tide_Merged$Tide_Merged)
write.table(Tide_Merged, 'Tide_Merged.txt', row.names =T, sep = '\t', quote = FALSE)