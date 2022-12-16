writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
Sys.which("make")
pkgbuild::find_rtools(debug = TRUE)

options(stringsAsFactors = F)
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)

dir='E:/TCGA-LIHC/TCGA-Immune/ICI/DataFiles/Training Data'
GDSC1_Expr = readRDS(file=file.path(dir,'GDSC1_Expr (RMA Normalized and Log Transformed).rds'))
GDSC1_Res = readRDS(file = file.path(dir,"GDSC1_Res.rds"))
GDSC1_Res <- exp(GDSC1_Res) 
testExpr<-read.table(file = "expr237.txt", sep = "\t",header = T,row.names = 1)
testExpr<-as.matrix(testExpr)
testExpr[1:4,1:4]  
colnames(testExpr)=paste0('test',colnames(testExpr))
dim(testExpr)  
calcPhenotype(trainingExprData = GDSC1_Expr,
              trainingPtype = GDSC1_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )
