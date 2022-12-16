library(limma)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(WGCNA)
library(reshape2)
library(stringr)
library(ggvenn)
library(VennDiagram)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
# Load samples and gene expression data 
sample<-read.table(file = "SAMPLE.txt", sep = "\t",header = T)
expr<-read.table(file = "HiSeqV2.txt", sep = "\t",header = T,row.names = 1)
expr<-expr[,sample$Sample]
trait <- read.table("trait.txt",sep = "\t",header=T,row.names = 1)
updown<-read.table(file = "up&down_5000.txt", sep = "\t",header = T)
expr_5000<-expr_5000[updown$Gene,]
# Calculate and omit the proportion of missing value data
which(is.na(expr_5000))
loss = sum(is.na(expr_5000))
have = sum(complete.cases(expr_5000))
ratio = loss/(loss+have)
expr_5000 = na.omit(expr_5000)
m.mad <- apply(expr_5000,1,mad)
expr_5000Var <- expr_5000[which(m.mad >max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
expr_5000 <- as.data.frame(t(expr_5000Var))
# Detect missing values
gsg = goodSamplesGenes(expr_5000, verbose = 3)
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(expr_5000)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(expr_5000)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  expr_5000 = expr_5000[gsg$goodSamples, gsg$goodGenes]
}
# Clarify the number of samples and genes
nGenes = ncol(expr_5000)
nSamples = nrow(expr_5000)
# Perform a systematic clustering and check whether there are outlier samples
sampleTree = hclust(dist(expr_5000), method = "average")
pdf(file="Sample clustering.pdf", width=60, height=10)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",cex=1.3,cex.main=2,cex.lab=1.5,cex.main=2)
dev.off()

powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(expr_5000, powerVector=powers,verbose=5)
pdf(file="Scale independence.pdf", width=8, height=8)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"),cex=2,cex.sub=2,cex.lab=2,cex.main=3,cex.axis=1.5)
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=2,col="red")
abline(h=0.85,col="red")
dev.off()

pdf(file="Mean connectivity.pdf", width=8, height=8)
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"),cex=2,cex.sub=2,cex.lab=2,cex.main=3,cex.axis=1.5)
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=2, col="red")
dev.off()

power = sft$powerEstimate
net = blockwiseModules(expr_5000, power =6, minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.2,
                       numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs=TRUE, loadTOMs=TRUE, saveTOMFileBase = "expr_5000", verbose = 3)
table(net$colors)
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
pdf(file="Cluster Dendrogram.pdf", width=30, height=10)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],"Module colors",dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,cex.colorLabels=1,cex.lab=2,cex.main=4,cex.axis=2)
dev.off()

sample_colors <- numbers2colors(trait, signed = FALSE)
pdf(file="Sample clustering and trait.pdf", width=30, height=20)
plotDendroAndColors(expr_5000_tree, sample_colors,
                    groupLabels = colnames(trait),
                    cex.dendroLabels = 0.6,
                    marAll = c(1, 4, 3, 1),
                    cex.rowText = 0.01,
                    main = "Sample dendrogram and trait heatmap", sub="", xlab="",cex=1.2,cex.lab=2,cex.main=4)
dev.off()

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(expr_5000, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
# Output modules and the KME value of each module
for(module in substring(colnames(MEs),3)){
  if(module == "grey") next
  ME=as.data.frame(MEs[,paste("ME",module,sep="")])
  colnames(ME)=module
  datModExpr=expr_5000[,moduleColors==module]
  datKME = signedKME(datModExpr, ME)
  datKME=cbind(datKME,rep(module,length(datKME)))
  write.table(datKME,quote = F,row.names = T,append = T,file = "All_Gene_KME.txt",col.names = F)
}
moduleTraitCor = cor(MEs, trait , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
# Display the correlation values within a heatmap plot
pdf(file="relationships.pdf", width=16, height=12)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(trait),
               cex.lab.y=0.8,
               cex.lab.x=0.85,
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = colorRampPalette(c("PowDerBlue","white","LightCoral"))(50),
               textMatrix = textMatrix,
               setStdMargins = T,
               cex.text = 0.9,
               cex.main = 2.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
# Draw the adjacency matrix between modules
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap",plotDendrograms = F,marDendro = c(4,4,2,4))
MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
pdf(file="Eigengene adjacency heatmap.pdf", width=8, height=10)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()

geneTree = net$dendrograms[[1]]
TOM<-TOMsimilarityFromExpr(expr_5000, power = 6)
dissTOM = 1 - TOM
nSelect = 400
set.seed(10)
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]
sizeGrWindow(9,9)
plotDiss = selectTOM^6
diag(plotDiss) = NA
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes",col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'))
pdf(file="Networkheatmap.pdf", width=10, height=12)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes",col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'))
dev.off()
# Export the network into the edge and node list
probes = colnames(expr_5000)
dimnames(TOM) <- list(probes, probes)
cyt = exportNetworkToCytoscape(TOM,edgeFile = paste("wgcna_cluster", ".edges.txt", sep=""),nodeFile = paste("wgcna_cluster", ".nodes.txt", sep=""),weighted = TRUE, threshold = 0,nodeNames =  probes, nodeAttr = moduleColors)

#Extract genes from the blue module
module = "blue"
column = match("MEblue", colnames(MEs_col))
moduleGenes = moduleColors==module
blue_module<-as.data.frame(dimnames(data.frame(expr_5000))[[2]][moduleGenes]) 
names(blue_module)="genename"
hub<- abs(geneModuleMembership$MM2)>0.80 & abs(geneTraitSignificance)>0.2
hubgene_blue <- dimnames(data.frame(expr_5000))[[2]][hub]
MM<-abs(geneModuleMembership[moduleGenes,column])
GS<-abs(geneTraitSignificance[moduleGenes, 1])
blue<-as.data.frame(cbind(MM,GS))
rownames(blue)=blue_module$genename
match<- blue_module$genename %in% hubgene_blue
blue$group<-match
write.table(blue, file="blue.txt", sep="\t", quote=F, row.names = T)

module = "salmon"
column = match("MEsalmon", colnames(MEs_col))
moduleGenes = moduleColors==module
salmon_module<-as.data.frame(dimnames(data.frame(expr_5000))[[2]][moduleGenes]) 
names(salmon_module)="genename"
hub<- abs(geneModuleMembership$MM2)>0.80 & abs(geneTraitSignificance)>0.2
hubgene_salmon <- dimnames(data.frame(expr_5000))[[2]][hub]
MM<-abs(geneModuleMembership[moduleGenes,column])
GS<-abs(geneTraitSignificance[moduleGenes, 1])
salmon<-as.data.frame(cbind(MM,GS))
rownames(salmon)=salmon_module$genename
match<- salmon_module$genename %in% hubgene_salmon
salmon$group<-match
write.table(salmon, file="salmon.txt", sep="\t", quote=F, row.names = T)

module = "grey60"
column = match("MEgrey60", colnames(MEs_col))
moduleGenes = moduleColors==module
grey60_module<-as.data.frame(dimnames(data.frame(expr_5000))[[2]][moduleGenes]) 
names(grey60_module)="genename"
hub<- abs(geneModuleMembership$MM2)>0.80 & abs(geneTraitSignificance)>0.2
hubgene_grey60 <- dimnames(data.frame(expr_5000))[[2]][hub]
MM<-abs(geneModuleMembership[moduleGenes,column])
GS<-abs(geneTraitSignificance[moduleGenes, 1])
grey60<-as.data.frame(cbind(MM,GS))
rownames(grey60)=grey60_module$genename
match<- grey60_module$genename %in% hubgene_grey60
grey60$group<-match
write.table(grey60, file="grey60.txt", sep="\t", quote=F, row.names = T)

module = "yellow"
column = match("MEyellow", colnames(MEs_col))
moduleGenes = moduleColors==module
yellow_module<-as.data.frame(dimnames(data.frame(expr_5000))[[2]][moduleGenes]) 
names(yellow_module)="genename"
hub<- abs(geneModuleMembership$MM2)>0.80 & abs(geneTraitSignificance)>0.2
hubgene_yellow <- dimnames(data.frame(expr_5000))[[2]][hub]
MM<-abs(geneModuleMembership[moduleGenes,column])
GS<-abs(geneTraitSignificance[moduleGenes, 1])
yellow<-as.data.frame(cbind(MM,GS))
rownames(yellow)=yellow_module$genename
match<- yellow_module$genename %in% hubgene_yellow
yellow$group<-match
write.table(yellow, file="yellow.txt", sep="\t", quote=F, row.names = T)
