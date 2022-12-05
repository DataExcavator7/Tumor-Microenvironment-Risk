#Preprocessing microarray gene expression data
library(affy) 
library(oligo)
library(annotate) 
library(hgu133plus2.db)
library(pd.hta.2.0)
library(limma)
library(dplyr)
library(stringr)
#Extracting microarray probes
setwd("E:/TCGA-LIHC/GSE121248")
celfiles_HCC <- list.files("E:/TCGA-LIHC/GSE121248",  "\\.gz$")
data.raw_HCC<- read.celfiles(filenames = file.path("E:/TCGA-LIHC/GSE121248/", celfiles_HCC))
treats_blank <- strsplit("HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent", " ")[[1]]
treats <- strsplit("HCC1 HCC2 HCC3 HCC4 HCC5 HCC6 HCC7 HCC8 HCC9 HCC10 HCC11 HCC12 HCC13 HCC14 HCC15 HCC16 HCC17 HCC18 HCC19 HCC20 HCC21 HCC22 HCC23 HCC24 HCC25 HCC26 HCC27 HCC28 HCC29 HCC30 HCC31 HCC32 HCC33 HCC34 HCC35 HCC36 HCC37 HCC38 HCC39 HCC40 HCC41 HCC42 HCC43 HCC44 HCC45 HCC46 HCC47 HCC48 HCC49 HCC50 HCC51 HCC52 HCC53 HCC54 HCC55 HCC56 HCC57 HCC58 HCC59 HCC60 HCC61 HCC62 HCC63 HCC64 HCC65 HCC66 HCC67 HCC68 HCC69 HCC70 HCC_adjacent1 HCC_adjacent2 HCC_adjacent3 HCC_adjacent4 HCC_adjacent5 HCC_adjacent6 HCC_adjacent7 HCC_adjacent8 HCC_adjacent9 HCC_adjacent10 HCC_adjacent11 HCC_adjacent12 HCC_adjacent13 HCC_adjacent14 HCC_adjacent15 HCC_adjacent16 HCC_adjacent17 HCC_adjacent18 HCC_adjacent19 HCC_adjacent20 HCC_adjacent21 HCC_adjacent22 HCC_adjacent23 HCC_adjacent24 HCC_adjacent25 HCC_adjacent26 HCC_adjacent27 HCC_adjacent28 HCC_adjacent29 HCC_adjacent30 HCC_adjacent31 HCC_adjacent32 HCC_adjacent33 HCC_adjacent34 HCC_adjacent35 HCC_adjacent36 HCC_adjacent37", " ")[[1]]
sampleNames(data.raw_HCC) <- treats
pData(data.raw_HCC)$index <- treats_blank
sampleNames(data.raw_HCC)
pData(data.raw_HCC)
fit1 <- fitProbeLevelModel(data.raw_HCC)
boxplot(fit1, names = NA, col = as.factor(treats_blank))
legend("topright", legend = unique(treats_blank), fill = as.factor(unique(treats_blank)),box.col = NA, bg = "white", inset = 0.01)
data.eset <- rma(data.raw_HCC)
data.exprs <- exprs(data.eset)
#Annotation according to the soft file
ff <- "E:/TCGA-LIHC/GSE121248/GSE121248_family.soft.gz"
nn <- grep("^[^#!^]", readLines(ff))[1] - 1
pfinfo <- read.table(ff, sep = "\t", quote = "", header = TRUE, skip = nn, fill = TRUE)
colnames(pfinfo)
pfinfo <- pfinfo[, c(1,11)]
pfinfo$Gene.Symbol <- toupper(pfinfo$Gene.Symbol)
pfinfo <- pfinfo[!duplicated(pfinfo), ]
pfinfo <- pfinfo[pfinfo$Gene.Symbol != "", ]
rownames(pfinfo) <- pfinfo$ID
nrow(pfinfo)
result <- dif[rownames(dif) %in% pfinfo$ID, ]
nrow(result)
result1$agi <- pfinfo[rownames(result1), 2]
write.table (result, file ="HCC_anno.txt", sep="\t", quote=F, row.names = T,col.names =T)
library(stringr)
split<-str_split(result$agi,"///")
agi<-sapply(split,"[",1)
Gene<-sapply(split,"[",2)
result$Gene<-Gene
#Searching for DEGs
group<-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)         
design <- model.matrix(~ -1+factor(group))
colnames(design)<-c("HCC","HCC_adjacent") 
contrast.matrix <- makeContrasts(HCC-HCC_adjacent,levels=design)
fit <- lmFit(data.exprs, design)
fit1 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit1)
dif <- topTable(fit2, coef = 1,  n = nrow(fit2), lfc = 0) 
write.table (dif, file ="HCC_GSE121248_updown.txt", sep="\t", quote=F, row.names = T,col.names =T)

#Extracting microarray probes
setwd("E:/TCGA-LIHC/GSE136247")
celfiles_HCC <- list.files("E:/TCGA-LIHC/GSE136247",  "\\.gz$")
data.raw_HCC<- read.celfiles(filenames = file.path("E:/TCGA-LIHC/GSE136247/", celfiles_HCC))
treats_blank <- strsplit("HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent HCC_adjacent", " ")[[1]]
treats <- strsplit("HCC1 HCC2 HCC3 HCC4 HCC5 HCC6 HCC7 HCC8 HCC9 HCC10 HCC11 HCC12 HCC13 HCC14 HCC15 HCC16 HCC17 HCC18 HCC19 HCC20 HCC21 HCC22 HCC23 HCC24 HCC25 HCC26 HCC27 HCC28 HCC29 HCC30 HCC31 HCC32 HCC33 HCC34 HCC35 HCC36 HCC37 HCC38 HCC39 HCC_adjacent1 HCC_adjacent2 HCC_adjacent3 HCC_adjacent4 HCC_adjacent5 HCC_adjacent6 HCC_adjacent7 HCC_adjacent8 HCC_adjacent9 HCC_adjacent10 HCC_adjacent11 HCC_adjacent12 HCC_adjacent13 HCC_adjacent14 HCC_adjacent15 HCC_adjacent16 HCC_adjacent17 HCC_adjacent18 HCC_adjacent19 HCC_adjacent20 HCC_adjacent21 HCC_adjacent22 HCC_adjacent23 HCC_adjacent24 HCC_adjacent25 HCC_adjacent26 HCC_adjacent27 HCC_adjacent28 HCC_adjacent29 HCC_adjacent30", " ")[[1]]
sampleNames(data.raw_HCC) <- treats
pData(data.raw_HCC)$index <- treats_blank
sampleNames(data.raw_HCC)
pData(data.raw_HCC)
fit1 <- fitProbeLevelModel(data.raw_HCC)
boxplot(fit1, names = NA, col = as.factor(treats_blank))
legend("topright", legend = unique(treats_blank), fill = as.factor(unique(treats_blank)),box.col = NA, bg = "white", inset = 0.01)
data.eset <- rma(data.raw_HCC)
data.exprs <- exprs(data.eset)
#Annotation according to the soft file
ff <- "E:/TCGA-LIHC/GSE136247/GSE136247_family.soft.gz"
nn <- grep("^[^#!^]", readLines(ff))[1] - 1
pfinfo <- read.table(ff, sep = "\t", quote = "", header = TRUE, skip = nn, fill = TRUE)
colnames(pfinfo)
pfinfo <- pfinfo[, c(1,11)]
pfinfo$Gene.Symbol <- toupper(pfinfo$Gene.Symbol)
pfinfo <- pfinfo[!duplicated(pfinfo), ]
pfinfo <- pfinfo[pfinfo$Gene.Symbol != "", ]
rownames(pfinfo) <- pfinfo$ID
nrow(pfinfo)
result <- dif[rownames(dif) %in% pfinfo$ID, ]
nrow(result)
result1$agi <- pfinfo[rownames(result1), 2]
write.table (result, file ="HCC_anno.txt", sep="\t", quote=F, row.names = T,col.names =T)
library(stringr)
split<-str_split(result$agi,"///")
agi<-sapply(split,"[",1)
Gene<-sapply(split,"[",2)
result$Gene<-Gene
#Searching for DEGs
group<-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)         
design <- model.matrix(~ -1+factor(group))
colnames(design)<-c("HCC","HCC_adjacent") 
contrast.matrix <- makeContrasts(HCC-HCC_adjacent,levels=design)
fit <- lmFit(data.exprs, design)
fit1 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit1)
dif <- topTable(fit2, coef = 1,  n = nrow(fit2), lfc = 0) 
write.table (dif, file ="HCC_GSE136247_updown.txt", sep="\t", quote=F, row.names = T,col.names =T)

#Extracting microarray probes
setwd("E:/TCGA-LIHC/GSE76297")
celfiles_HCC <- list.files("E:/TCGA-LIHC/GSE76297/HCC",  "\\.gz$")
data.raw_HCC<- read.celfiles(filenames = file.path("E:/TCGA-LIHC/GSE76297/HCC/", celfiles_HCC))
treats_blank <- strsplit("HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC HCC_adjacent HCC", " ")[[1]]
treats <- strsplit("HCC1 HCC_adjacent1 HCC2 HCC_adjacent2 HCC3 HCC_adjacent3 HCC4 HCC_adjacent4 HCC5 HCC_adjacent5 HCC6 HCC_adjacent6 HCC7 HCC_adjacent7 HCC8 HCC_adjacent8 HCC9 HCC_adjacent9 HCC10 HCC_adjacent10 HCC11 HCC_adjacent11 HCC12 HCC_adjacent12 HCC13 HCC_adjacent13 HCC14 HCC_adjacent14 HCC15 HCC_adjacent15 HCC16 HCC_adjacent16 HCC17 HCC_adjacent17 HCC18 HCC_adjacent18 HCC19 HCC_adjacent19 HCC20 HCC_adjacent20 HCC21 HCC22 HCC_adjacent22 HCC23 HCC_adjacent23 HCC24 HCC_adjacent24 HCC25 HCC_adjacent25 HCC26 HCC_adjacent26 HCC27 HCC_adjacent27 HCC28 HCC29 HCC_adjacent29 HCC30 HCC_adjacent30 HCC31 HCC_adjacent31 HCC32 HCC_adjacent32 HCC33 HCC_adjacent33 HCC34 HCC_adjacent34 HCC35 HCC_adjacent35 HCC36 HCC_adjacent36 HCC37 HCC_adjacent37 HCC38 HCC_adjacent38 HCC39 HCC_adjacent39 HCC40 HCC_adjacent40 HCC41 HCC_adjacent41 HCC42 HCC_adjacent42 HCC43 HCC_adjacent43 HCC44 HCC_adjacent44 HCC45 HCC_adjacent45 HCC46 HCC_adjacent46 HCC47 HCC_adjacent47 HCC48 HCC_adjacent48 HCC49 HCC_adjacent49 HCC50 HCC_adjacent50 HCC51 HCC_adjacent51 HCC52 HCC_adjacent52 HCC53 HCC_adjacent53 HCC54 HCC_adjacent54 HCC55 HCC_adjacent55 HCC56 HCC_adjacent56 HCC57 HCC_adjacent57 HCC58 HCC_adjacent58 HCC59 HCC_adjacent59 HCC60 HCC_adjacent60 HCC61", " ")[[1]]
sampleNames(data.raw_HCC) <- treats
pData(data.raw_HCC)$index <- treats_blank
sampleNames(data.raw_HCC)
pData(data.raw_HCC)
fit1 <- fitProbeLevelModel(data.raw_HCC)
boxplot(fit1, names = NA, col = as.factor(treats_blank))
legend("topright", legend = unique(treats_blank), fill = as.factor(unique(treats_blank)),box.col = NA, bg = "white", inset = 0.01)
data.eset <- rma(data.raw_HCC)
data.exprs <- exprs(data.eset)
#Annotation according to the soft file
ff <- "E:/TCGA-LIHC/GSE76297/GSE76297_family.soft"
nn <- grep("^[^#!^]", readLines(ff))[1] - 1
pfinfo <- read.table(ff, sep = "\t", quote = "", header = TRUE, skip = nn, fill = TRUE)
colnames(pfinfo)
pfinfo <- pfinfo[, c(1, 8)]
pfinfo$gene_assignment <- toupper(pfinfo$gene_assignment)
pfinfo <- pfinfo[!duplicated(pfinfo), ]
pfinfo <- pfinfo[pfinfo$gene_assignment != "", ]
rownames(pfinfo) <- pfinfo$ID
nrow(pfinfo)
result <- dif[rownames(dif) %in% pfinfo$ID, ]
nrow(result1)
result$agi <- pfinfo[rownames(result), 2]
write.table (result, file ="HCC_anno.txt", sep="\t", quote=F, row.names = T,col.names =T)
split<-str_split(result$agi,"//")
agi<-sapply(split,"[",1)
Gene<-sapply(split,"[",2)
result$Gene<-Gene
#Searching for DEGs
group<-c(1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,1,2,1,2,1,2,1,2,1,2,1,2,1,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1) 
design <- model.matrix(~ -1+factor(group))
colnames(design)<-c("HCC","HCC_adjacent") 
contrast.matrix <- makeContrasts(HCC-HCC_adjacent,levels=design)
fit <- lmFit(data.exprs, design)
fit1 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit1)
dif <- topTable(fit2, coef = 1,  n = nrow(fit2), lfc = 0) 
write.table (dif, file ="HCC_GSE76297_updown.txt", sep="\t", quote=F, row.names = T,col.names =T)