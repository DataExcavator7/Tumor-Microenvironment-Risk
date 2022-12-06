setwd("E:/TCGA-LIHC/TMErisk_Group")
library(limma)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggthemes)
# Fig.4B, Fig.S5A: Boxplots compare the expression of the critical genes 
raw <- read.csv("GENE_IMM.csv",row.names = 1,header = T)
# Select different range of columns to illutrate different kinds of genes
dd <- raw %>% as.data.frame() %>% rownames_to_column("sample") %>% pivot_longer(cols = 3:11,names_to = "GENE",values_to = "VALUE")
plot.info <- dd[,c(1,2,109,110)]
levels <- plot.info %>% group_by(GENE) %>% summarise(m = median(VALUE)) %>% arrange(desc(m)) %>% pull(GENE)
plot.info$GENE = factor(plot.info$GENE,levels = levels)
plot.info$TMErisk.Score <- factor(plot.info$TMErisk.Score,levels = c("high","low"))
ggboxplot(plot.info,x = "GENE",y = "VALUE",color = "black",fill = "TMErisk.Score",
          outlier.colour = NA,xlab = "", palette = c("#00468BB2" ,"#ED0000B2","#42B540B2"),
          ylab = "Log2(FPKM+1)",main = "") +
  theme_base() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1), 
        legend.position="top",legend.title=element_text(size=15), 
        legend.text = element_text(size=15))+
  stat_compare_means(aes(group = TMErisk.Score,label = "kruskal.test"),
                     label = "p.signif",na.rm = TRUE)

# Fig.4C: Search DEGs in low- and high-TMErisk group
data <- read.table("HiSeqV2.txt",sep = "\t",header=T,row.names=1)
SAMPLE<- read.table("SAMPLE.txt",sep = "\t",header=T)
expr<-data[,SAMPLE$sampleID]
write.table (expr, file ="expr237.txt", sep="\t", quote=F, row.names = T,col.names =T)
group<-c(1,1,1,1,1,2,2,2,1,1,2,2,1,1,1,1,2,1,1,1,2,1,2,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,2,2,2,1,2,1,2,2,1,1,1,1,2,2,1,1,1,1,1,1,1,2,2,1,1,2,1,2,2,2,1,1,1,1,2,1,2,1,1,2,2,2,2,2,2,1,2,1,1,2,1,1,2,2,2,2,1,1,1,2,1,2,1,2,2,1,1,1,2,1,2,1,1,1,2,1,2,2,2,1,2,2,2,1,1,1,1,1,2,1,1,1,1,2,2,1,2,2,1,2,2,1,1,2,2,1,2,1,2,1,1,1,1,1,1,2,2,2,1,2,2,1,1,2,2,2,1,2,2,2,1,1,1,2,2,2,1,2,1,1,1,1,2,1,1,1,1,2,2,2,1,2,1,2,2,1,2,1,1,2,1,1,1,1,1,1,1,1,1,2,2,1,2,2,2,1,1,1,2,1,1,1,1,1,1,1,2)
design <- model.matrix(~ -1+factor(group)) 
colnames(design)<-c("High_Risk","Low_Risk") 
contrast.matrix <- makeContrasts(High_Risk-Low_Risk,levels=design) 

fit <- lmFit(expr, design)
fit1 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit1)
#results<-decideTests(fit2, method="global", adjust.method="none", p.value=0.05, lfc=0.585)
#summary(results)
dif <- topTable(fit2, coef = 1,  n = nrow(fit2), lfc = 0) 
#dif <- topTable(fit2,coef=1,number=nrow(fit2) ,adjust.method="BH",sort.by="B",resort.by="M")
write.table (dif, file ="TMErisk.txt", sep="\t", quote=F, row.names = T,col.names =T)
dif<-read.table("TMErisk.txt",sep = "\t",header=T,row.names=1)
dif$logP<- -log10(dif$adj.P.Val)
dif$Group="Not-significant"
dif$Group[which((dif$adj.P.Val<0.05) & (dif$logFC>=1))]="Up-regulated"
dif$Group[which((dif$adj.P.Val<0.05) & (dif$logFC<=(-1)))]="Down-regulated"

dif$Label=""
dif<-dif[order(dif$adj.P.Val),]
upgenes<-head(rownames(dif)[which(dif$Group=="Up-regulated")],15)
upgenes
downgenes<-head(rownames(dif)[which(dif$Group=="Down-regulated")],15)
top20<-c(as.character(upgenes),as.character(downgenes))
dif$Label[match(top20,rownames(dif))]<- top20
write.table (dif, file ="TMErisk.txt", sep="\t", quote=F, row.names = T,col.names =T)
dif<-read.table("TMErisk.txt",sep = "\t",header=T,row.names=1)

pdf(file="volcano.pdf", width=12,height=10)
ggscatter(dif,x="logFC",y="logP",
          color = "Group",
          palette = c("#2f5688","#BBBBBB","#CC0000"),
          size=1.5,
          label=dif$Label,
          font.label=16,
          repel=T,
          max.overlaps=500,
          xlab = "logFC",ylab = "-logFDR")+
  theme_base()+
  geom_hline(yintercept=1.3,linetype="dashed")+
  geom_vline(xintercept=c(-1,1),linetype="dashed")+
  coord_cartesian(xlim = c(-2.25,2.25),ylim = c(0,20))
dev.off()

# Fig.4D: GO and KEGG enrichment analysis
setwd("E:/TCGA-LIHC/TMErisk_Group")
rvcheck="https://cran.r-project.org/src/contrib/Archive/rvcheck/rvcheck_0.1.8.tar.gz"
library(rvcheck)
library(clusterProfiler)
library(DOSE)
library(stringr)
library(org.Hs.eg.db)
library(enrichplot)
library(GOplot)
# GO enrichment analysis
pvalueFilter=0.05       
qvalueFilter=0.05       
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}
rt=read.table("Low.txt", header=T, sep="\t", check.names=F)  
colnames(rt)[1]="Gene"
genes=as.vector(rt[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]      
kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
write.table(GO, file="GO_Low.txt", sep="\t", quote=F, row.names = F)
showNum=10
if(nrow(GO)<30){
  showNum=nrow(GO)
}

pdf(file="GO_Low.pdf", width=10, height=10)
bar=dotplot(kk, showCategory=showNum, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()

#KEGG enrichment analysis
pvalueFilter=0.05       
qvalueFilter=0.05       
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}
genes=as.vector(rt[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]      
kk=enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$id[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))#此行可不要，会把geneid冲掉
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
write.table(KEGG, file="KEGG_Low.txt", sep="\t", quote=F, row.names = F)
showNum=30
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}

pdf(file="KEGG_barplotlow.pdf", width=10, height=10)
barplot(kk, drop=TRUE, showCategory=showNum, color=colorSel)
dev.off()
pdf(file="KEGG_dotplotlow.pdf", width=10, height=10)
dotplot(kk, showCategory=showNum, color=colorSel)
dev.off()