library(Seurat)
library(dplyr)
library(Matrix)
library(magrittr)
library(stringr)
library(tidyverse)
options(stringsAsFactors = F)

metadata = read.table("bcc_all_metadata.txt",sep = "\t",  header=T,check.names = F)
raw.data = read.csv("bcc_scRNA_counts.txt",sep = "\t", row.names=1, header=T,check.names = F)

seruat= CreateSeuratObject(counts = raw.data_BCC, metadata = metadata)
seruat <- AddMetaData(seruat, metadata = metadata)
seruat <- FindVariableFeatures(seruat, selection.method = "vst",verbose = FALSE)

seruat <- ScaleData(seruat,features = rownames(seruat))
seruat <- RunPCA(seruat, npcs = 20, verbose = FALSE)

seruat <- FindNeighbors(seruat, reduction = "pca", dims = 1:5)
seruat <- FindClusters(seruat, resolution = 0.4)

seruat <- RunTSNE(seruat, reduction = "pca", dims = 1:5,check_duplicates = FALSE)
pdf(file="tsnetreatment.pdf", width=10, height=10)
DimPlot(seruat,reduction = "tsne",group.by="Response",pt.size = 1)+
  theme(axis.title.x =element_text(size=26), axis.title.y=element_text(size=26),axis.text = element_text(size = 20),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"),plot.title = element_text(hjust = 0.5,size = 32,face = "bold"),plot.margin = unit(rep(2,4),'lines'))
dev.off()

scale<-seruat[["RNA"]]@counts
scale<-scale[c("SYNE1","CAMK4","KLRB1","ITPRIP"),]
write.table(as.data.frame(scale),file="4gene.txt",row.names=T, quote=FALSE,sep="\t")
TMErisk<- read.table("GSE123813_counts.txt" , sep="\t",row.names=1, header=T)
TMErisk <- as(as.matrix(TMErisk), "dgCMatrix")
seruat[["RNA"]]@counts=TMErisk

pdf(file="SYNE1.pdf", width=10, height=10)
FeaturePlot(seruat, features = c("SYNE1"),slot="counts",reduction = "tsne")+
  theme(axis.title.x =element_text(size=26), axis.title.y=element_text(size=26),axis.text = element_text(size = 20),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"),plot.title = element_text(hjust = 0.5,size = 32,face = "bold"),plot.margin = unit(rep(2,4),'lines'))
dev.off()
pdf(file="CAMK4.pdf", width=10, height=10)
FeaturePlot(seruat, features = c("CAMK4"),slot="counts",reduction = "tsne")+
  theme(axis.title.x =element_text(size=26), axis.title.y=element_text(size=26),axis.text = element_text(size = 20),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"),plot.title = element_text(hjust = 0.5,size = 32,face = "bold"),plot.margin = unit(rep(2,4),'lines'))
dev.off()
pdf(file="ITPRIP.pdf", width=10, height=10)
FeaturePlot(seruat, features = c("ITPRIP"),slot="counts",reduction = "tsne")+
  theme(axis.title.x =element_text(size=26), axis.title.y=element_text(size=26),axis.text = element_text(size = 20),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "blacazP;.k"),plot.title = element_text(hjust = 0.5,size = 32,face = "bold"),plot.margin = unit(rep(2,4),'lines'))
dev.off()
pdf(file="KLRB1.pdf", width=10, height=10)
FeaturePlot(seruat, features = c("KLRB1"),slot="counts",reduction = "tsne")+
  theme(axis.title.x =element_text(size=26), axis.title.y=element_text(size=26),axis.text = element_text(size = 20),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"),plot.title = element_text(hjust = 0.5,size = 32,face = "bold"),plot.margin = unit(rep(2,4),'lines'))
dev.off()
pdf(file="TMErisk.pdf", width=10, height=10)
FeaturePlot(seruat, features = c("TMErisk"),slot="counts",reduction = "tsne",cols=c("PowDerBlue","LightCoral"))+
  theme(axis.title.x =element_text(size=26), axis.title.y=element_text(size=26),axis.text = element_text(size = 20),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"),plot.title = element_text(hjust = 0.5,size = 32,face = "bold"),plot.margin = unit(rep(2,4),'lines'))
dev.off()
