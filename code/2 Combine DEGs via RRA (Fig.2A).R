library(pheatmap)
library(RobustRankAggreg)
padj=0.05
logFC=1
files=c("HCC_GSE76297.txt","HCC_GSE121248.txt","HCC_GSE136247.txt")
upList=list()
downList=list()
allFCList=list()
for(i in 1:length(files)){
  inputFile=files[i]
  rt=read.table(inputFile,header=T,sep="\t") 
  header=unlist(strsplit(inputFile,"\\."))
  downList[[header[1]]]=as.vector(rt[,1])
  upList[[header[1]]]=rev(as.vector(rt[,1]))
  fcCol=rt[,1:2]
  colnames(fcCol)=c("Gene",header[[1]])
  allFCList[[header[1]]]=fcCol
}

mergeLe=function(x,y){
  merge(x,y,by="Gene",all=T)}
newTab=Reduce(mergeLe,allFCList)
ids = unique(newTab[,1])
newTab = newTab[!duplicated(newTab[,1]),]
rownames(newTab) = ids
rownames(newTab)=newTab[,1]
newTab=newTab[,2:ncol(newTab)]
newTab[is.na(newTab)]=0
write.table(newTab,file="newTab.txt",sep="\t",quote=F,row.names=F)

upMatrix = rankMatrix(upList)
upAR = aggregateRanks(rmat=upMatrix)
colnames(upAR)=c("Name","Pvalue")
upAdj=p.adjust(upAR$Pvalue,method="BH")
upXls=cbind(upAR,adjPvalue=upAdj)
upFC=newTab[as.vector(upXls[,1]),]
upXls=cbind(upXls,logFC=rowMeans(upFC))
upSig=upXls[(upXls$adjPvalue<0.05),]
downMatrix = rankMatrix(downList)
downAR = aggregateRanks(rmat=downMatrix)
colnames(downAR)=c("Name","Pvalue")
downAdj=p.adjust(downAR$Pvalue,method="BH")
downXls=cbind(downAR,adjPvalue=downAdj)
downFC=newTab[as.vector(downXls[,1]),]
downXls=cbind(downXls,logFC=rowMeans(downFC))
downSig=downXls[(downXls$adjPvalue<0.05),]
allSig = rbind(upSig,downSig)
colnames(allSig)
allSig = allSig[,c("Name","logFC")]
write.table(allSig,file = 'allSign.txt',sep = '\t',quote = F)

hminput=newTab[c(as.vector(upSig[1:20,1]),as.vector(downSig[1:20,1])),]
pdf(file="RRA.pdf", width=6,height=9)
pheatmap(hminput,display_numbers = TRUE,
         fontsize_row=10,
         fontsize_col=12,
         color = colorRampPalette(c("PowDerBlue","white","LightCoral"))(50),
         cluster_cols = FALSE,cluster_rows = FALSE)
dev.off()