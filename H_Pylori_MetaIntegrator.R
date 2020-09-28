###############################################################################################
## Mohamed Omar
## 21/07/2019
## Goal: Discovery and validation of a small gene signature that can predict the metastasis potential in primary prostate cancer
################################################################################################

## Clean work space
rm(list = ls())

## Set the working directory
setwd("/Volumes/Macintosh/Research/Projects/H_Pylori")


## Load necessary packages
library(MetaIntegrator)
library(GEOquery)
library(pROC)
library(caret)
library(genefilter)
library(mltools)
library(illuminaHumanv3.db)
library(annotate)
library(precrec)
library(patchwork)

########################################################


# HPyloriData <- getGEOData(c("GSE27411", "GSE60427", "GSE60662", "GSE5081"))
# CellLinesData <- getGEOData(c("GSE39919", "GSE70394", "GSE74577", "GSE74492"))
# 
# # Get each dataset separately
# Dataset1 <- HPyloriData$originalData$GSE27411
# Dataset2 <- HPyloriData$originalData$GSE60427
# Dataset3 <- HPyloriData$originalData$GSE60662
# Dataset4 <- HPyloriData$originalData$GSE5081
# 
# ValDataset1 <- CellLinesData$originalData$GSE39919
# ValDataset2 <- CellLinesData$originalData$GSE70394
# ValDataset3 <- CellLinesData$originalData$GSE74577
# ValDataset4 <- CellLinesData$originalData$GSE74492
# 

#CanDataset1 <- getGEOData(c("GSE65801"))
#CanDataset1 <- CanDataset1$originalData$GSE65801
#save(CanDataset1, file = "./Data/CanDataSet1.rda")

# save(Dataset1, Dataset2, Dataset3, Dataset4, file = "./Data/HPyloriData.rda")
# save(ValDataset1, ValDataset2, ValDataset3, ValDataset4, file = "./Data/HPyloriCellLineData.rda")

## Load Training data sets (4)
load("./Data/HPyloriData.rda")

## Load the validation datasets (cell lines)
load("./Data/HPyloriCellLineData.rda")

## Load the cancer dataset 1
load("./Data/CanDataSet1.rda")

# load Cancer DataSet2 Pheno
canPheno2 <- read.delim("./Data/EMTAB1440/EMTAB1440Pheno.txt")

# load Cancer DataSet2 Expr
canExpr2 <- read.delim("./Data/EMTAB1440/ExpressionNormalized_v2.txt")

###############
#######################################################

## Getting the phenotype data for each data set
pheno1 <- Dataset1$pheno
pheno2 <- Dataset2$pheno
pheno3 <- Dataset3$pheno
pheno4 <- Dataset4$pheno

Valpheno1 <- ValDataset1$pheno
Valpheno2 <- ValDataset2$pheno
Valpheno3 <- ValDataset3$pheno
Valpheno4 <- ValDataset4$pheno

canPheno1 <- CanDataset1$pheno 
################################

## Getting the expression data for each data set
## load expr4
#load("/Users/mohamedomar/Documents/Research/Projects/Prostate/Data/expr4.rda")
#ProstateData$originalData$GSE46691$expr <- expr4
expr1 <- Dataset1$expr
expr2 <- Dataset2$expr
expr3 <- Dataset3$expr
expr4 <- Dataset4$expr

Valexpr1 <- ValDataset1$expr
Valexpr2 <- ValDataset2$expr
Valexpr3 <- ValDataset3$expr
Valexpr4 <- ValDataset4$expr

canExpr1 <- CanDataset1$expr

## Checking if the expression data are normalized and log2 transformed
# boxplot(expr1[,1:15], outline= FALSE)
# boxplot(expr2[,1:15], outline= FALSE)
# boxplot(expr3[,1:15], outline = FALSE)
# boxplot(expr4[,1:15], outline= FALSE)
# 
# boxplot(Valexpr1[,1:15], outline = FALSE)
# boxplot(Valexpr2[,1:6], outline = FALSE)
# boxplot(Valexpr3[,1:6], outline = FALSE)
# boxplot(Valexpr4[,1:6], outline = FALSE)
# 
# boxplot(canExpr1[,1:14], outline = FALSE)
#################################################################
## Create a list containing training data sets
AllDataSets <- list(Dataset1, Dataset2, Dataset3, Dataset4)
names(AllDataSets) <- c(Dataset1$formattedName, Dataset2$formattedName, Dataset3$formattedName, Dataset4$formattedName)

###################################################################

## Annotate expression

## Expr1
head(rownames(expr1))
rownames(expr1) <- Dataset1$keys
expr1 <- expr1[!is.na(rownames(expr1)), ]

#####################
## expr2
head(rownames(expr2))
rownames(expr2) <- Dataset2$keys
expr2 <- expr2[!is.na(rownames(expr2)), ]
dim(expr2)


#####################
## expr3
head(rownames(expr3))
rownames(expr3) <- Dataset3$keys
expr3 <- expr3[!is.na(rownames(expr3)), ]
dim(expr3)


# #######################
# expr4
head(rownames(expr4))
rownames(expr4) <- Dataset4$keys
expr4 <- expr4[!is.na(rownames(expr4)), ]
dim(expr4)

# #######################
# Valexpr1
head(rownames(Valexpr1))
rownames(Valexpr1) <- ValDataset1$keys
Valexpr1 <- Valexpr1[!is.na(rownames(Valexpr1)), ]
dim(Valexpr1)

# #######################
# Valexpr2
head(rownames(Valexpr2))
rownames(Valexpr2) <- ValDataset2$keys
Valexpr2 <- Valexpr2[!is.na(rownames(Valexpr2)), ]
dim(Valexpr2)

# #######################
# Valexpr3
head(rownames(Valexpr3))
rownames(Valexpr3) <- ValDataset3$keys
Valexpr3 <- Valexpr3[!is.na(rownames(Valexpr3)), ]
dim(Valexpr3)

# #######################
# Valexpr4
head(rownames(Valexpr4))
rownames(Valexpr4) <- ValDataset4$keys
Valexpr4 <- Valexpr4[!is.na(rownames(Valexpr4)), ]
dim(Valexpr4)

# #######################
# cancerExpr1
head(rownames(canExpr1))
rownames(canExpr1) <- CanDataset1$keys
canExpr1 <- canExpr1[!is.na(rownames(canExpr1)), ]
dim(canExpr1)

# #######################
# cancerExpr2
head(rownames(canExpr2))

# annotate canExpr2 (illumina)
tmp <- canExpr2$Hybridization.REF
canExpr2$GeneSymbol<- mapIds(illuminaHumanv3.db,
                              keys=tmp,
                              column="SYMBOL",
                              keytype="PROBEID",
                              multiVals="first")

canExpr2 <- canExpr2[!duplicated(canExpr2$GeneSymbol), ]
canExpr2 <- canExpr2[!is.na(canExpr2$GeneSymbol), ]

rownames(canExpr2) <- canExpr2$GeneSymbol
canExpr2$Hybridization.REF <- NULL
canExpr2$GeneSymbol <- NULL
dim(canExpr2)

# Finally modify the sample names to match those in the phenotype
colnames(canExpr2) <- gsub("X", "", colnames(canExpr2))
colnames(canExpr2) <- tolower(colnames(canExpr2))
colnames(canExpr2)

# Convert to numeric matrix
COLS <- colnames(canExpr2)
ROWS <- rownames(canExpr2)
canExpr2 <- matrix(as.numeric(unlist(canExpr2)),nrow=nrow(canExpr2))
rownames(canExpr2) <- ROWS
colnames(canExpr2) <- COLS
# ############################################################
####################################################################
#### Modify the phenotypes

# Pheno1

pheno1$DiseaseStatus <- pheno1$`disease state:ch1`
pheno1$DiseaseStatus[pheno1$DiseaseStatus == "non-infected"] <- "control"
pheno1$DiseaseStatus[pheno1$DiseaseStatus == "infected"] <- "case"
pheno1$DiseaseStatus[pheno1$DiseaseStatus == "atrophy"] <- "case"
pheno1$DiseaseStatus <- factor(pheno1$DiseaseStatus, levels = c("control","case"))
table(pheno1$DiseaseStatus)
all(rownames(pheno1) == colnames(expr1))
# 
Dataset1$pheno <- pheno1
Dataset1$expr <- expr1
Dataset1$keys <- rownames(expr1)
########## 

## Modify pheno2

pheno2$DiseaseStatus <- pheno2$`gastritis grade:ch1`
pheno2$DiseaseStatus[pheno2$DiseaseStatus == "normal"] <- "control"
pheno2$DiseaseStatus[pheno2$DiseaseStatus %in% c("mild", "IM", "severe")] <- "case"
pheno2$DiseaseStatus <- factor(pheno2$DiseaseStatus, levels = c("control", "case"))
table(pheno2$DiseaseStatus)

# 
all(rownames(pheno2) == colnames(expr2))

## Finally, replace the expression and phenotype data in the dataset with the new modified versions
Dataset2$expr <- expr2
Dataset2$pheno <- pheno2
Dataset2$keys <- rownames(expr2)

######
# Modify pheno3
pheno3$DiseaseStatus <- pheno3$title
pheno3$DiseaseStatus[1:4] <- "control"
pheno3$DiseaseStatus[5:length(pheno3$DiseaseStatus)] <- "case"

pheno3$DiseaseStatus <- factor(pheno3$DiseaseStatus, levels = c("control", "case"))
table(pheno3$DiseaseStatus)

all(rownames(pheno3) == colnames(expr3))

## Finally, replace the expression and phenotype data in the dataset with the new modified versions
Dataset3$pheno <- pheno3
Dataset3$expr <- expr3
Dataset3$keys <- rownames(expr3)

####### 
## Modify pheno4
pheno4$DiseaseStatus <- pheno4$title
pheno4$DiseaseStatus <- gsub("\\,.+", "", pheno4$DiseaseStatus)
pheno4$DiseaseStatus <- gsub("Gastric biopsy ", "", pheno4$DiseaseStatus)
pheno4$DiseaseStatus[pheno4$DiseaseStatus == "HP- ER+"] <- "case"
pheno4$DiseaseStatus[pheno4$DiseaseStatus == "HP+ ER+"] <- "case"
pheno4$DiseaseStatus[pheno4$DiseaseStatus == "HP+ ER-"] <- "case"
pheno4$DiseaseStatus[pheno4$DiseaseStatus == "HP- ER-"] <- "control"

pheno4$DiseaseStatus <- factor(pheno4$DiseaseStatus, levels = c("control", "case"))
table(pheno4$DiseaseStatus)

## Modify sample names to match sample names of pheno4
all(rownames(pheno4) == colnames(expr4))

## Finally, replace the expression and phenotype data in the dataset with the new modified versions
Dataset4$expr <- expr4
Dataset4$pheno <- pheno4
Dataset4$keys <- rownames(expr4)

################## 
## Modify valPheno1
Valpheno1$DiseaseStatus <- Valpheno1$source_name_ch1
Valpheno1 <- Valpheno1[1:8, ]
Valpheno1$DiseaseStatus[1:4] <- "case"
Valpheno1$DiseaseStatus[5:8] <- "control"

Valpheno1$DiseaseStatus <- factor(Valpheno1$DiseaseStatus, levels = c("control", "case"))
table(Valpheno1$DiseaseStatus)

## Modify sample names to match sample names of valpheno1
Valexpr1 <- Valexpr1[, colnames(Valexpr1) %in% rownames(Valpheno1)]
all(rownames(Valpheno1) == colnames(Valexpr1))

## Finally, replace the expression and phenotype data in the dataset with the new modified versions
ValDataset1$expr <- Valexpr1
ValDataset1$pheno <- Valpheno1
ValDataset1$keys <- rownames(Valexpr1)

################## 
## Modify valPheno2
Valpheno2$DiseaseStatus <- Valpheno2$`infection:ch1`
Valpheno2$DiseaseStatus[Valpheno2$DiseaseStatus == "uninfected"] <- "control"
Valpheno2$DiseaseStatus[Valpheno2$DiseaseStatus == "H. pylori strain 60190"] <- "case"

Valpheno2$DiseaseStatus <- factor(Valpheno2$DiseaseStatus, levels = c("control", "case"))
table(Valpheno2$DiseaseStatus)

## Modify sample names to match sample names of valpheno2
all(rownames(Valpheno2) == colnames(Valexpr2))

## Finally, replace the expression and phenotype data in the dataset with the new modified versions
ValDataset2$expr <- Valexpr2
ValDataset2$pheno <- Valpheno2
ValDataset2$keys <- rownames(Valexpr2)

################## 
## Modify valPheno3
Valpheno3$DiseaseStatus <- Valpheno3$title
Valpheno3$DiseaseStatus[1:3] <- "control"
Valpheno3$DiseaseStatus[4:6] <- "case"

Valpheno3$DiseaseStatus <- factor(Valpheno3$DiseaseStatus, levels = c("control", "case"))
table(Valpheno3$DiseaseStatus)

## Modify sample names to match sample names of valpheno3
all(rownames(Valpheno3) == colnames(Valexpr3))

## Finally, replace the expression and phenotype data in the dataset with the new modified versions
ValDataset3$expr <- Valexpr3
ValDataset3$pheno <- Valpheno3
ValDataset3$keys <- rownames(Valexpr3)

################## 
## Modify valPheno4
Valpheno4$DiseaseStatus <- Valpheno4$title
Valpheno4$DiseaseStatus[1:3] <- "control"
Valpheno4$DiseaseStatus[4:6] <- "case"

Valpheno4$DiseaseStatus <- factor(Valpheno4$DiseaseStatus, levels = c("control", "case"))
table(Valpheno4$DiseaseStatus)

## Modify sample names to match sample names of valpheno3
all(rownames(Valpheno4) == colnames(Valexpr4))

## Finally, replace the expression and phenotype data in the dataset with the new modified versions
ValDataset4$expr <- Valexpr4
ValDataset4$pheno <- Valpheno4
ValDataset4$keys <- rownames(Valexpr4)

################## 
## Modify cancer pheno1
canPheno1$DiseaseStatus <- canPheno1$`tissue:ch1`
canPheno1$DiseaseStatus[canPheno1$DiseaseStatus == "normal gastric tissue"] <- "control"
canPheno1$DiseaseStatus[canPheno1$DiseaseStatus == "gastric tumor"] <- "case"

canPheno1$DiseaseStatus <- factor(canPheno1$DiseaseStatus, levels = c("control", "case"))
table(canPheno1$DiseaseStatus)

## Modify sample names to match sample names of valpheno3
all(rownames(canPheno1) == colnames(canExpr1))

## Finally, replace the expression and phenotype data in the dataset with the new modified versions
CanDataset1$expr <- canExpr1
CanDataset1$pheno <- canPheno1
CanDataset1$keys <- rownames(canExpr1)

################## 
## Modify cancer pheno2

rownames(canPheno2) <- canPheno2$Source.Name


canPheno2$DiseaseStatus <- canPheno2$Factor.Value.disease.
canPheno2$DiseaseStatus[canPheno2$DiseaseStatus == "normal"] <- "control"
canPheno2$DiseaseStatus[canPheno2$DiseaseStatus == "gastric adenocarcinoma"] <- "case"

canPheno2$DiseaseStatus <- factor(canPheno2$DiseaseStatus, levels = c("control", "case"))
table(canPheno2$DiseaseStatus)

## Modify sample names to match sample names of valpheno3
all(rownames(canPheno2) == colnames(canExpr2))

## Finally, replace the expression and phenotype data in the dataset with the new modified versions
CanDataset2 <- list()
CanDataset2$expr <- matrix(as.numeric(unlist(canExpr2)),nrow=nrow(canExpr2))
CanDataset2$pheno <- canPheno2
CanDataset2$keys <- rownames(canExpr2)
CanDataset2$formattedName <- "EMTAB1440"

#########################################################################
############################################################################

## Label samples (All samples need to be assigned labels in the $class vector, 1 for ‘disease’ or 0 for ‘control’)
Dataset1 <- classFunction(Dataset1, column = "DiseaseStatus", diseaseTerms = c("case"))
Dataset2 <- classFunction(Dataset2, column = "DiseaseStatus", diseaseTerms = c("case"))
Dataset3 <- classFunction(Dataset3, column = "DiseaseStatus", diseaseTerms = c("case"))
Dataset4 <- classFunction(Dataset4, column = "DiseaseStatus", diseaseTerms = c("case"))

ValDataset1 <- classFunction(ValDataset1, column = "DiseaseStatus", diseaseTerms = c("case"))
ValDataset2 <- classFunction(ValDataset2, column = "DiseaseStatus", diseaseTerms = c("case"))
ValDataset3 <- classFunction(ValDataset3, column = "DiseaseStatus", diseaseTerms = c("case"))
ValDataset4 <- classFunction(ValDataset4, column = "DiseaseStatus", diseaseTerms = c("case"))

CanDataset1 <- classFunction(CanDataset1, column = "DiseaseStatus", diseaseTerms = c("case"))
CanDataset2 <- classFunction(CanDataset2, column = "DiseaseStatus", diseaseTerms = c("case"))

############################################################################
############################################################################
#########################################################################
## The metaanalysis

## Creating the meta object
AllDataSets <- list(Dataset1, Dataset2, Dataset3, Dataset4)
names(AllDataSets) <- c(Dataset1$formattedName, Dataset2$formattedName, Dataset3$formattedName, Dataset4$formattedName)

HPylori_Meta <- list()
HPylori_Meta$originalData <- AllDataSets

## Check the meta object before the metaanalysis
checkDataObject(HPylori_Meta, "Meta", "Pre-Analysis") ## If true, Proceed to the meta analysis

#HPylori_Meta <- geneSymbolCorrection(HPylori_Meta)

## Run the meta analysis
HPylori_metaanalysis <- runMetaAnalysis(HPylori_Meta, runLeaveOneOutAnalysis = F, maxCores = 3)

# effectSize = 0/ FDR = 0.1/ Nstudies = 3

## Filter out significant genes from the metaanalysis results (this will be the gene signature that separates Metas from No_Mets)
HPylori_metaanalysis <- filterGenes(HPylori_metaanalysis, isLeaveOneOut = F, effectSizeThresh = 0, FDRThresh = 0.05, numberStudiesThresh = 4)

## Assigning a name to the filter
filter <- HPylori_metaanalysis$filterResults[[1]]
filter

PositiveGenes <- filter$posGeneNames
NegativeGenes <- filter$negGeneNames
#save(PositiveGenes, NegativeGenes, file = "./Objs/NewSigGenes_Pre.rda")


## Summarize filter results
filter_summary <- summarizeFilterResults(metaObject = HPylori_metaanalysis, getMostRecentFilter(HPylori_metaanalysis))

## Save the filter
save(filter, file = "./Objs/filter.rda")

## Save a table of the positive genes and negative genes
#write.table(filter_summary$pos, file = "./Objs/Meta/Positive_genes_filter.csv", quote = TRUE, sep = "\t", col.names = TRUE, row.names = TRUE, dec = ".")
#write.table(filter_summary$neg, file = "./Objs/Meta/Negative_genes_filter.csv", quote = TRUE, sep = "\t", col.names = TRUE, row.names = TRUE, dec = ".")

## Modify the gene signature for more accuracy and AUC
# Using forward search 
#New_filter <- forwardSearch(metaObject = HPylori_metaanalysis, filterObject = filter)

#save(New_filter, file = "./Objs/NewFilter.rda")
load("./Objs/NewFilter.rda")

## Replace the old filter with the new smaller one
HPylori_metaanalysis$filterResults$FDR0.05_es0_nStudies4_looaFALSE_hetero0$posGeneNames <- New_filter$posGeneNames
HPylori_metaanalysis$filterResults$FDR0.05_es0_nStudies4_looaFALSE_hetero0$negGeneNames <- New_filter$negGeneNames

New_filter <- HPylori_metaanalysis$filterResults[[1]]
New_filter_summary <- summarizeFilterResults(metaObject = HPylori_metaanalysis, getMostRecentFilter(HPylori_metaanalysis))

## Save the tables of positive and negative genes
write.table(New_filter_summary$pos, file = "./Objs/NewFilter_Positive_genes.csv", quote = T, sep = "\t", col.names = T, row.names = T)
write.table(New_filter_summary$neg, file = "./Objs/NewFilter_Negative_genes.csv", quote = T, sep = "\t", col.names = T, row.names = T)

## Gene names
PositiveGenes <- New_filter$posGeneNames
PositiveGenes
NegativeGenes <- New_filter$negGeneNames
NegativeGenes

write.table(PositiveGenes, file = "./UpRegulatedGenes.txt")
write.table(NegativeGenes, file = "./DownRegulatedGenes.txt")

#New_SignatureGenes <- c(PositiveGenes, NegativeGenes)
#save(New_SignatureGenes, file = "./Objs/NewSigGenes.rda")

## Create a summary ROC curve (Training data sets)
set.seed(333)
png(filename = "./Figs/Pooled_ROC_TrainingDatasets.png", width = 2000, height = 2000, res = 300)
pooledROCPlot(metaObject = HPylori_metaanalysis, filterObject = New_filter)
dev.off()



#### ## Load the new genes from Mohamed
NewGenes <- read.delim("./Data/FromMohamed/0-Heli_compined_gene_signature.txt")
## Replace the old filter with the new smaller one
HPylori_metaanalysis$filterResults$FDR0.05_es0_nStudies4_looaFALSE_hetero0$posGeneNames <- NewGenes$Upregulated
HPylori_metaanalysis$filterResults$FDR0.05_es0_nStudies4_looaFALSE_hetero0$negGeneNames <- NewGenes$Downregulated[-c(25:31)]

New_filter2 <- HPylori_metaanalysis$filterResults[[1]]
New_filter2_summary <- summarizeFilterResults(metaObject = HPylori_metaanalysis, getMostRecentFilter(HPylori_metaanalysis))

## Create a summary ROC curve (Training data sets)
set.seed(333)
png(filename = "./Figs/Pooled_ROC_TrainingDatasets_NewFilter2.png", width = 2000, height = 2000, res = 300)
pooledROCPlot(metaObject = HPylori_metaanalysis, filterObject = New_filter2)
dev.off()

## Effect size visualization in each study
# For the up-regulated genes
#pdf(file = "./Figs/Meta/UpGenes_effect_size.pdf", title = "Effect size of up-regulated genes across the discovery data sets", width = 20, height = 10)
#par(mfrow=c(2,2))

#forestPlot(metaObject = HPylori_metaanalysis, geneName = "CAMK2N1", textColor = "black")
#forestPlot(metaObject = HPylori_metaanalysis, geneName = "ASNS", textColor = "black")
#forestPlot(metaObject = HPylori_metaanalysis, geneName = "NCOA2", textColor = "black")
#forestPlot(metaObject = HPylori_metaanalysis, geneName = "PTPN9", textColor = "black")
#forestPlot(metaObject = HPylori_metaanalysis, geneName = "GNPTAB", textColor = "black")
#forestPlot(metaObject = HPylori_metaanalysis, geneName = "IQGAP3", textColor = "black")
#forestPlot(metaObject = HPylori_metaanalysis, geneName = "RPRML", textColor = "black")
#forestPlot(metaObject = HPylori_metaanalysis, geneName = "NCOA2", textColor = "black")
#dev.off()

# For the down-regulated genes
# pdf(file = "./Figs/Meta/DownGenes_effect_size.pdf", title = "Effect size of Down-regulated genes across the discovery data sets", width = 20, height = 12)
# par(mfrow=c(3,4))
# 
# forestPlot(metaObject = HPylori_metaanalysis, geneName = "AZGP1", textColor = "black")
# forestPlot(metaObject = HPylori_metaanalysis, geneName = "CPA3", textColor = "black")
# #forestPlot(metaObject = HPylori_metaanalysis, geneName = "BRE", textColor = "black")
# forestPlot(metaObject = HPylori_metaanalysis, geneName = "UFM1", textColor = "black")
# forestPlot(metaObject = HPylori_metaanalysis, geneName = "DPT", textColor = "black")
# forestPlot(metaObject = HPylori_metaanalysis, geneName = "PART1", textColor = "black")
# forestPlot(metaObject = HPylori_metaanalysis, geneName = "CBLL1", textColor = "black")
# forestPlot(metaObject = HPylori_metaanalysis, geneName = "LUZP2", textColor = "black")
# forestPlot(metaObject = HPylori_metaanalysis, geneName = "CHRNA2", textColor = "black")
# forestPlot(metaObject = HPylori_metaanalysis, geneName = "WNT8B", textColor = "black")
# forestPlot(metaObject = HPylori_metaanalysis, geneName = "NT5E", textColor = "black")
# forestPlot(metaObject = HPylori_metaanalysis, geneName = "DCXR", textColor = "black")
# forestPlot(metaObject = HPylori_metaanalysis, geneName = "EDN3", textColor = "black")
# #forestPlot(metaObject = HPylori_metaanalysis, geneName = "NT5DC1", textColor = "black")
# #forestPlot(metaObject = HPylori_metaanalysis, geneName = "COL4A5", textColor = "black")

#dev.off()

######################################
## Heatmap of the effect sizes of the signature genes

#png("./Figs/Meta/Heatmap_EffectSizes.png", width = 4000, height = 2000, res = 300)
heatmapPlot(HPylori_metaanalysis, filterObject = New_filter2)
#dev.off()

#####################################################
# get thresholds from the training data
TrainScore1 <- calculateScore(filterObject = New_filter2, datasetObject = HPylori_Meta$originalData$GSE27411)
thr_train1 <- coords(roc(pheno1$DiseaseStatus, TrainScore1, levels = c("control", "case"), direction = "<"),"best", transpose = TRUE)["threshold"]
coords(roc(pheno1$DiseaseStatus, TrainScore1, levels = c("control", "case"), direction = "<"),"local maximas", transpose = TRUE)

#Train_Predictions1 <- ifelse(TrainScore1 >= thr_train1, "Mets", "No_Mets")
#confusionMatrix(as.factor(Train_Predictions1), pheno2$Metastasis, positive = "Mets")

#mcc(preds = as.factor(Train_Predictions1), actuals = pheno2$Metastasis)

TrainScore2 <- calculateScore(filterObject = New_filter2, datasetObject = HPylori_Meta$originalData$GSE60427)
thr_train2 <- coords(roc(pheno2$DiseaseStatus, TrainScore2, levels = c("control", "case"), direction = "<"),"best", transpose = TRUE)["threshold"]
coords(roc(pheno2$DiseaseStatus, TrainScore2, levels = c("control", "case"), direction = "<"),"local maximas", transpose = TRUE)

TrainScore3 <- calculateScore(filterObject = New_filter2, datasetObject = HPylori_Meta$originalData$GSE60662)
thr_train3 <- coords(roc(pheno3$DiseaseStatus, TrainScore3, levels = c("control", "case"), direction = "<"),"best", transpose = TRUE)["threshold"]
coords(roc(pheno3$DiseaseStatus, TrainScore3, levels = c("control", "case"), direction = "<"),"local maximas", transpose = TRUE)

TrainScore4 <- calculateScore(filterObject = New_filter2, datasetObject = HPylori_Meta$originalData$GSE5081)
thr_train4 <- coords(roc(pheno4$DiseaseStatus, TrainScore4, levels = c("control", "case"), direction = "<"),"best", transpose = TRUE)["threshold"]
coords(roc(pheno4$DiseaseStatus, TrainScore4, levels = c("control", "case"), direction = "<"),"local maximas", transpose = TRUE)

Mean_Thr <- mean(c(thr_train1, thr_train2, thr_train3, thr_train4))
Mean_Thr
###########################################################################################
#############################################################################

## The next step is validation on indepndent data set (GSE116918, and GSE16560)

##############################################

### Now lets examine the performance of our filter on the independent data set

#### Using ROC curve

# Validation dataset1
png(filename = "./Figs/ROC_Test1_Signature2.png", width = 2000, height = 2000, res = 300)
rocPlot(datasetObject = ValDataset1, filterObject = New_filter)
dev.off()

# Validation dataset2

png(filename = "./Figs/ROC_Test2_Signature2.png", width = 2000, height = 2000, res = 300)
rocPlot(datasetObject = ValDataset2, filterObject = New_filter)
dev.off()

# Validation dataset3

png(filename = "./Figs/ROC_Test3_Signature2.png", width = 2000, height = 2000, res = 300)
rocPlot(datasetObject = ValDataset3, filterObject = New_filter)
dev.off()

# Validation dataset4

png(filename = "./Figs/ROC_Test4_Signature2.png", width = 2000, height = 2000, res = 300)
rocPlot(datasetObject = ValDataset4, filterObject = New_filter)
dev.off()


## Cancer dataset1
png(filename = "./Figs/ROC_Test_Cancer1_Signature2.png", width = 2000, height = 2000, res = 300)
rocPlot(datasetObject = CanDataset1, filterObject = New_filter2)
dev.off()

## Cancer dataset2
png(filename = "./Figs/ROC_Test_Cancer2_Signature2.png", width = 2000, height = 2000, res = 300)
rocPlot(datasetObject = CanDataset2, filterObject = New_filter2)
dev.off()
################################

#### Using PRC plot

# val dataset1
png(filename = "./Figs/PRC_Test1_Signature2.png", width = 2000, height = 2000, res = 300)
prcPlot(datasetObject = ValDataset1, filterObject = New_filter)
dev.off()

# val dataset2
png(filename = "./Figs/PRC_Test2_Signature2.png", width = 2000, height = 2000, res = 300)
prcPlot(datasetObject = ValDataset2, filterObject = New_filter)
dev.off()

# val dataset3
png(filename = "./Figs/PRC_Test3_Signature2.png", width = 2000, height = 2000, res = 300)
prcPlot(datasetObject = ValDataset3, filterObject = New_filter)
dev.off()

# val dataset4
png(filename = "./Figs/PRC_Test4_Signature2.png", width = 2000, height = 2000, res = 300)
prcPlot(datasetObject = ValDataset4, filterObject = New_filter)
dev.off()

# cancer dataset1
png(filename = "./Figs/PRC_Test5_Signature2.png", width = 2000, height = 2000, res = 300)
prcPlot(datasetObject = CanDataset1, filterObject = New_filter2)
dev.off()

# cancer dataset2
png(filename = "./Figs/PRC_Test6_Signature2.png", width = 2000, height = 2000, res = 300)
prcPlot(datasetObject = CanDataset2, filterObject = New_filter2)
dev.off()
##################################

#### Using violin plot

# Val Dataset1 
png(filename = "./Figs/ViolinPlot_Test1_Signature2.png", width = 2000, height = 2000, res = 300)
violinPlot(filterObject = New_filter2, datasetObject = ValDataset1, labelColumn = "DiseaseStatus")
dev.off()

# Val Dataset2
png(filename = "./Figs/ViolinPlot_Test2_Signature2.png", width = 2000, height = 2000, res = 300)
violinPlot(filterObject = New_filter2, datasetObject = ValDataset2, labelColumn = "DiseaseStatus")
dev.off()

# Val Dataset3
png(filename = "./Figs/ViolinPlot_Test3_Signature2.png", width = 2000, height = 2000, res = 300)
violinPlot(filterObject = New_filter2, datasetObject = ValDataset3, labelColumn = "DiseaseStatus")
dev.off()

# Val Dataset4
png(filename = "./Figs/ViolinPlot_Test4_Signature2.png", width = 2000, height = 2000, res = 300)
violinPlot(filterObject = New_filter2, datasetObject = ValDataset4, labelColumn = "DiseaseStatus")
dev.off()

# cancer dataset1
png(filename = "./Figs/ViolinPlot_Test_Cancer1_Signature2.png", width = 2000, height = 2000, res = 300)
violinPlot(filterObject = New_filter2, datasetObject = CanDataset1, labelColumn = "DiseaseStatus")
dev.off()

# cancer dataset2
png(filename = "./Figs/ViolinPlot_Test_Cancer2_Signature2.png", width = 2000, height = 2000, res = 300)
violinPlot(filterObject = New_filter2, datasetObject = CanDataset2, labelColumn = "DiseaseStatus")
dev.off()



###################################################
# Cancer Dataset1
# Violin1 <- violinPlot(filterObject = New_filter2, datasetObject = CanDataset1, labelColumn = "DiseaseStatus")
# library(precrec)
# sscurves1 <- evalmod(scores = canPheno1$score, labels = canPheno1$DiseaseStatus)
# PRC1 <- autoplot(sscurves1, curvetype = c("PRC"))
# ROC1 <- autoplot(sscurves1, curvetype = c("ROC"))
# 
# library(cowplot)
# png(file = "./Figs/CanDataSet1Multi.png", width = 2000, height = 2000)
# plot_grid(ROC1, PRC1, Violin1, 
#           labels = c("A", "B", "C"),
#           ncol = 2, nrow = 2)
# dev.off()
# 


##########################
## Cancer Dataset2
# Violin2 <- violinPlot(filterObject = New_filter2, datasetObject = CanDataset2, labelColumn = "DiseaseStatus")
# sscurves2 <- evalmod(scores = canPheno2$score, labels = canPheno2$DiseaseStatus)
# PRC2 <- autoplot(sscurves2, curvetype = c("PRC"))
# ROC2 <- autoplot(sscurves2, curvetype = c("ROC"))
# 
# png(file = "./Figs/CanDataSet2Multi.png", width = 2000, height = 2000)
# plot_grid(ROC2, PRC2, Violin2, 
#           labels = c("A", "B", "C"),
#           ncol = 2, nrow = 2)
# dev.off()

##########################
## Cell line Dataset1
# Violin_CellLine1 <- violinPlot(filterObject = New_filter2, datasetObject = ValDataset1, labelColumn = "DiseaseStatus")
# sscurves_CellLine1 <- evalmod(scores = Valpheno1$score, labels = Valpheno1$DiseaseStatus)
# PRC_CellLine1 <- autoplot(sscurves_CellLine1, curvetype = c("PRC"))
# ROC_CellLine1 <- autoplot(sscurves_CellLine1, curvetype = c("ROC"))
# 
# png(file = "./Figs/CellLineDataSet1Multi.png", width = 2000, height = 2000)
# plot_grid(ROC_CellLine1, PRC_CellLine1, Violin_CellLine1, 
#           labels = c("A", "B", "C"),
#           ncol = 2, nrow = 2)
# dev.off()


##########################
## Cell line Dataset2
# Violin_CellLine2 <- violinPlot(filterObject = New_filter2, datasetObject = ValDataset2, labelColumn = "DiseaseStatus")
# sscurves_CellLine2 <- evalmod(scores = Valpheno2$score, labels = Valpheno2$DiseaseStatus)
# PRC_CellLine2 <- autoplot(sscurves_CellLine2, curvetype = c("PRC"))
# ROC_CellLine2 <- autoplot(sscurves_CellLine2, curvetype = c("ROC"))
# 
# png(file = "./Figs/CellLineDataSet2Multi.png", width = 2000, height = 2000)
# plot_grid(ROC_CellLine2, PRC_CellLine2, Violin_CellLine2, 
#           labels = c("A", "B", "C"),
#           ncol = 2, nrow = 2)
# dev.off()

##########################
## Cell line Dataset3
# Violin_CellLine3 <- violinPlot(filterObject = New_filter2, datasetObject = ValDataset3, labelColumn = "DiseaseStatus")
# sscurves_CellLine3 <- evalmod(scores = Valpheno3$score, labels = Valpheno3$DiseaseStatus)
# PRC_CellLine3 <- autoplot(sscurves_CellLine3, curvetype = c("PRC"))
# ROC_CellLine3 <- autoplot(sscurves_CellLine3, curvetype = c("ROC"))
# 
# png(file = "./Figs/CellLineDataSet3Multi.png", width = 2000, height = 2000)
# plot_grid(ROC_CellLine3, PRC_CellLine3, Violin_CellLine3, 
#           labels = c("A", "B", "C"),
#           ncol = 2, nrow = 2)
# dev.off()
# 
##########################
## Cell line Dataset4
# Violin_CellLine4 <- violinPlot(filterObject = New_filter2, datasetObject = ValDataset4, labelColumn = "DiseaseStatus")
# sscurves_CellLine4 <- evalmod(scores = Valpheno4$score, labels = Valpheno4$DiseaseStatus)
# PRC_CellLine4 <- autoplot(sscurves_CellLine4, curvetype = c("PRC"))
# ROC_CellLine4 <- autoplot(sscurves_CellLine4, curvetype = c("ROC"))
# 
# png(file = "./Figs/CellLineDataSet4Multi.png", width = 2000, height = 2000)
# plot_grid(ROC_CellLine4, PRC_CellLine4, Violin_CellLine4, 
#           labels = c("A", "B", "C"),
#           ncol = 2, nrow = 2)
# dev.off()



# #####################################################

##### Calculate ROC data

# Val Dataset1
#scoreRs1 <- calculateScore(filterObject = New_filter, datasetObject = valDataSet1)
#rocRes1 <- calculateROC(predictions = scoreRs1, labels = valDataSet1$class)

# Val Dataset2
#scoreRs2 <- calculateScore(filterObject = New_filter, datasetObject = valDataSet2)
#rocRes2 <- calculateROC(predictions = scoreRs2, labels = valDataSet2$class)

##############
#### Calculate a signature score (Z score) and add it to the phenotype table

# Val Dataset1
Valpheno1$score <- calculateScore(filterObject = New_filter2, datasetObject = ValDataset1)

# Val Dataset2
Valpheno2$score <- calculateScore(filterObject = New_filter2, datasetObject = ValDataset2)

# Val Dataset3
Valpheno3$score <- calculateScore(filterObject = New_filter2, datasetObject = ValDataset3)

# Val Dataset4
Valpheno4$score <- calculateScore(filterObject = New_filter2, datasetObject = ValDataset4)

# cancer dataset1
canPheno1$score <- calculateScore(filterObject = New_filter2, datasetObject = CanDataset1)

# cancer dataset2
canPheno2$score <- calculateScore(filterObject = New_filter2, datasetObject = CanDataset2)

#########################################

## ROC curve in the Testing data using pROC

#png(filename = "./Figs/Meta/ROC_Test_pROC.png", width = 2000, height = 2000, res = 300)
#roc(val_pheno1$Metastasis, val_pheno1$score, plot = TRUE, print.auc=TRUE, levels = c("No_Mets", "Mets"), col="blue", lwd=2, grid=TRUE, auc = TRUE, ci = TRUE)
#dev.off()


#############################
##### Predictions

# Val Dataset1
Test_Predictions1 <- ifelse(Valpheno1$score >= Mean_Thr, "case", "control")
confusionMatrix(as.factor(Test_Predictions1), Valpheno1$DiseaseStatus, positive = "case")

mcc(preds = as.factor(Test_Predictions1), actuals = Valpheno1$DiseaseStatus)

# Val Dataset2
Test_Predictions2 <- ifelse(Valpheno2$score >= Mean_Thr, "case", "control")
confusionMatrix(as.factor(Test_Predictions2), Valpheno2$DiseaseStatus, positive = "case")

mcc(preds = as.factor(Test_Predictions2), actuals = Valpheno2$DiseaseStatus)

# Val Dataset3
Test_Predictions3 <- ifelse(Valpheno3$score >= Mean_Thr, "case", "control")
confusionMatrix(as.factor(Test_Predictions3), Valpheno3$DiseaseStatus, positive = "case")

mcc(preds = as.factor(Test_Predictions3), actuals = Valpheno3$DiseaseStatus)

# Val Dataset4
Test_Predictions4 <- ifelse(Valpheno4$score >= Mean_Thr, "case", "control")
confusionMatrix(as.factor(Test_Predictions4), Valpheno4$DiseaseStatus, positive = "case")

mcc(preds = as.factor(Test_Predictions4), actuals = Valpheno4$DiseaseStatus)

# cancer dataset1
Test_Predictions5 <- ifelse(canPheno1$score >= Mean_Thr, "case", "control")
confusionMatrix(as.factor(Test_Predictions5), canPheno1$DiseaseStatus, positive = "case")

mcc(preds = as.factor(Test_Predictions5), actuals = canPheno1$DiseaseStatus)

# cancer dataset2
Test_Predictions6 <- ifelse(canPheno2$score >= Mean_Thr, "case", "control")
confusionMatrix(as.factor(Test_Predictions6), canPheno2$DiseaseStatus, positive = "case")

mcc(preds = as.factor(Test_Predictions6), actuals = canPheno2$DiseaseStatus)



############################################
#############################################
## Big plot1: Discovery performance + validation in the 2 cancer datasets using intersection Signature

# Performance in the 1st cancer dataset
sscurves_CanDataSet1 <- evalmod(scores = canPheno1$score, labels = canPheno1$DiseaseStatus)

# Performance in the 2nd cancer dataset
sscurves_CanDataSet2  <- evalmod(scores = canPheno2$score, labels = canPheno2$DiseaseStatus)

My_Theme = theme(
  axis.title.x = element_text(size = 12),
  axis.text.x = element_text(size = 12),
  axis.title.y = element_text(size = 12),
  axis.text.y = element_text(size = 12),
  plot.title = element_text(size=18)
)

# Pooled AUC in the discovery datasets
P1 <- pooledROCPlot(metaObject = HPylori_metaanalysis, filterObject = New_filter2)

# Violin, AUC, and AUPRC in the 1st cancer dataset
P2 <- violinPlot(filterObject = New_filter2, datasetObject = CanDataset1, labelColumn = "DiseaseStatus") + My_Theme
P3 <- autoplot(sscurves_CanDataSet1, curvetype = c("ROC")) + labs(title = "ROC") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.812"), size =5)
P4 <- autoplot(sscurves_CanDataSet1, curvetype = c("PRC")) +  labs(title = "PRC") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.825"), size = 5)

# Violin, AUC, and AUPRC in the 2nd cancer dataset
P5 <- violinPlot(filterObject = New_filter2, datasetObject = CanDataset2, labelColumn = "DiseaseStatus")+My_Theme
P6 <- autoplot(sscurves_CanDataSet2, curvetype = c("ROC")) + labs(title = "ROC") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.82"), size = 5)
P7 <- autoplot(sscurves_CanDataSet2, curvetype = c("PRC")) + labs(title = "PRC") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.835"), size = 5)

## Plot ALL
tiff(filename = "./Figs/Intersec_Sig.tiff", width = 1500, height = 1000)
P1 + ((P2 / P3 / P4 + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 18))) | (P5 / P6 / P7 + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 18)))) +
  plot_layout(widths = c(0.5, 1)) + 
  plot_annotation(
    title = 'Refined Signature',
    tag_levels = c('A', '1'),
    theme = theme(plot.title = element_text(size = 20, face = "bold"))
  )
dev.off()


#############################################
## Big plot2: Discovery performance + validation in the 4 cell line datasets and the 2 cancer datasets using the original Signature

# Performance in the 1st cell line dataset
sscurves_CellLineDataSet1 <- evalmod(scores = Valpheno1$score, labels = Valpheno1$DiseaseStatus)

# Performance in the 2nd cell line dataset
sscurves_CellLineDataSet2 <- evalmod(scores = Valpheno2$score, labels = Valpheno2$DiseaseStatus)

# Performance in the 3rd cell line dataset
sscurves_CellLineDataSet3 <- evalmod(scores = Valpheno3$score, labels = Valpheno3$DiseaseStatus)

# Performance in the 4th cell line dataset
sscurves_CellLineDataSet4 <- evalmod(scores = Valpheno4$score, labels = Valpheno4$DiseaseStatus)

# Performance in the 1st cancer dataset
sscurves_CanDataSet1 <- evalmod(scores = canPheno1$score, labels = canPheno1$DiseaseStatus)

# Performance in the 2nd cancer dataset
sscurves_CanDataSet2  <- evalmod(scores = canPheno2$score, labels = canPheno2$DiseaseStatus)

My_Theme = theme(
  axis.title.x = element_text(size = 10),
  axis.text.x = element_text(size = 10),
  axis.title.y = element_text(size = 10),
  axis.text.y = element_text(size = 10),
  plot.title = element_text(size=15)
)

# Pooled AUC in the discovery datasets
set.seed(333)
P1 <- pooledROCPlot(metaObject = HPylori_metaanalysis, filterObject = New_filter)

# Violin, AUC, and AUPRC in the 1st cell line dataset
P2 <- violinPlot(filterObject = New_filter, datasetObject = ValDataset1, labelColumn = "DiseaseStatus")+ My_Theme
P3 <- autoplot(sscurves_CellLineDataSet1, curvetype = c("ROC")) + labs(title = "ROC") + annotate("text", x = .65, y = .25, label = paste("AUC = 1.00"), size =3)
P4 <- autoplot(sscurves_CellLineDataSet1, curvetype = c("PRC")) +  labs(title = "PRC") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 1.00"), size = 3)

# Violin, AUC, and AUPRC in the 2nd cell line dataset
P5 <- violinPlot(filterObject = New_filter, datasetObject = ValDataset2, labelColumn = "DiseaseStatus") + My_Theme
P6 <- autoplot(sscurves_CellLineDataSet2, curvetype = c("ROC")) + labs(title = "ROC") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.889"), size =3)
P7 <- autoplot(sscurves_CellLineDataSet2, curvetype = c("PRC")) +  labs(title = "PRC") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.903"), size = 3)

# Violin, AUC, and AUPRC in the 3rd cell line dataset
P8 <- violinPlot(filterObject = New_filter, datasetObject = ValDataset3, labelColumn = "DiseaseStatus") + My_Theme
P9 <- autoplot(sscurves_CellLineDataSet3, curvetype = c("ROC")) + labs(title = "ROC") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.333"), size =3)
P10 <- autoplot(sscurves_CellLineDataSet3, curvetype = c("PRC")) +  labs(title = "PRC") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.372"), size = 3)

# Violin, AUC, and AUPRC in the 4th cell line dataset
P11 <- violinPlot(filterObject = New_filter, datasetObject = ValDataset4, labelColumn = "DiseaseStatus") + My_Theme
P12 <- autoplot(sscurves_CellLineDataSet4, curvetype = c("ROC")) + labs(title = "ROC") + annotate("text", x = .65, y = .25, label = paste("AUC = 1.00"), size =3)
P13 <- autoplot(sscurves_CellLineDataSet4, curvetype = c("PRC")) +  labs(title = "PRC") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 1.00"), size = 3)

# Violin, AUC, and AUPRC in the 1st cancer dataset
P14 <- violinPlot(filterObject = New_filter, datasetObject = CanDataset1, labelColumn = "DiseaseStatus") + My_Theme
P15 <- autoplot(sscurves_CanDataSet1, curvetype = c("ROC")) + labs(title = "ROC") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.818"), size =3)
P16 <- autoplot(sscurves_CanDataSet1, curvetype = c("PRC")) +  labs(title = "PRC") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.81"), size = 3)

# Violin, AUC, and AUPRC in the 2nd cancer dataset
P17 <- violinPlot(filterObject = New_filter, datasetObject = CanDataset2, labelColumn = "DiseaseStatus") + My_Theme
P18 <- autoplot(sscurves_CanDataSet2, curvetype = c("ROC")) + labs(title = "ROC") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.817"), size = 3)
P19 <- autoplot(sscurves_CanDataSet2, curvetype = c("PRC")) + labs(title = "PRC") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.836"), size = 3)

## Plot ALL
tiff(filename = "./Figs/Original_Sig.tiff", width = 1500, height = 1000)
P1 + ((P2 / P3 / P4 + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 20))) | (P5 / P6 / P7+ plot_layout(tag_level = "new")& theme(plot.tag = element_text(size = 20))) |  (P8 / P9 / P10 + plot_layout(tag_level = "new")& theme(plot.tag = element_text(size = 20))) | (P11 / P12 / P13 + plot_layout(tag_level = "new")& theme(plot.tag = element_text(size = 20))) | (P14 / P15 / P16 + plot_layout(tag_level = "new")& theme(plot.tag = element_text(size = 20))) | (P17 / P18 / P19 +plot_layout(tag_level = "new")& theme(plot.tag = element_text(size = 20)))) + 
  plot_layout(widths = c(0.6, 2)) + 
  plot_annotation(
  title = 'Original Signature',
  tag_levels = c('A', '1'),
  theme = theme(plot.title = element_text(size = 20, face = "bold"))
  )
dev.off() 

################################
## Meta-Test
# Meta_Test <- list() 
# Meta_Test$originalData$GSE39919 <- ValDataset1
# Meta_Test$originalData$GSE70394 <- ValDataset2
# Meta_Test$originalData$GSE74577 <- ValDataset3
# Meta_Test$originalData$GSE74492 <- ValDataset4
# 
# Meta_Test <- geneSymbolCorrection(metaObject = Meta_Test)
# 
# #set.seed(333)
# #summaryROCPlot(metaObject = Meta_Test, filterObject = filter)
# 
# set.seed(333)
# png(filename = "./Figs/PooledROC_InitialSig_Testing.png", width = 2000, height = 2000, res = 300)
# pooledROCPlot(metaObject = Meta_Test, filterObject = New_filter)
# dev.off()

###############################################################################
## Gene set enrichment analysis
library(enrichR)
library(clusterProfiler)

# annotate canExpr2 (illumina)
Pos_EntID<- mapIds(org.Hs.eg.db,
                             keys=PositiveGenes,
                             column="ENTREZID",
                             keytype="SYMBOL",
                             multiVals="first")


Neg_EntID<- mapIds(org.Hs.eg.db,
                   keys=NegativeGenes,
                   column="ENTREZID",
                   keytype="SYMBOL",
                   multiVals="first")

My_Theme = theme(
  axis.title.x = element_text(size = 10),
  axis.text.x = element_text(size = 10),
  axis.title.y = element_text(size = 10),
  axis.text.y = element_text(size = 12),
  plot.title = element_text(size=10),
  legend.position="right"
)


## GO BPs
Pos_GO_BPs <- enrichGO(
  Pos_EntID,
  org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = FALSE
)


P1 <- dotplot(Pos_GO_BPs, showCategory=5, title = "GO BPs up-regulated genes") + My_Theme

Neg_GO_BPs <- enrichGO(
  Neg_EntID,
  org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = FALSE
)

P2 <- dotplot(Neg_GO_BPs, showCategory=5, title = "GO BPs down-regulated genes") + My_Theme

##########################3
## GO MFs
Pos_GO_MFs <- enrichGO(
  Pos_EntID,
  org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = FALSE
)


P3 <- dotplot(Pos_GO_MFs, showCategory=5, title = "GO MFs Up-regulated genes") + My_Theme

Neg_GO_MFs <- enrichGO(
  Neg_EntID,
  org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = FALSE
)

P4 <- dotplot(Neg_GO_MFs, showCategory=5, title = "GO MFs down-regulated genes") + My_Theme

#########################
## GO BPs
Pos_GO_CCs <- enrichGO(
  Pos_EntID,
  org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = FALSE
)


P5 <- dotplot(Pos_GO_CCs, showCategory=5, title = "GO CCs Up-regulated genes") + My_Theme

Neg_GO_CCs <- enrichGO(
  Neg_EntID,
  org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = FALSE
)

P6 <- dotplot(Neg_GO_CCs, showCategory=5, title = "GO CCs down-regulated genes") + My_Theme

################################
## KEGG
Pos_KEGG <- enrichKEGG(
  Pos_EntID,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
)

P7 <- dotplot(Pos_KEGG, showCategory=5, title = "KEGG pathways Up-regulated genes") + My_Theme


Neg_KEGG <- enrichKEGG(
  Neg_EntID,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
)


P8 <- dotplot(Neg_KEGG, showCategory=5, title = "KEGG pathways down-regulated genes") + My_Theme


## Plot ALL
tiff(filename = "./Figs/Enrichment.tiff", width = 1500, height = 1000)
((P1 / P2 + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 16))) | (P3 / P4 + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 16))) | (P5  + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 16))) | (P7  + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 16)))) +
  #plot_layout(widths = c(1, 1,1,1,1,1)) + 
  plot_annotation(
    title = 'Enrichment Analysis',
    #tag_levels = c('A', '1'),
    theme = theme(plot.title = element_text(size = 20, face = "bold"))
  )
dev.off()

# dbs <- listEnrichrDbs()
# dbs <- c("GO_Biological_Process_2018", "KEGG_2019_Human", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018")
# Enriched_PositiveGns <- enrichr(genes = PositiveGenes, databases = dbs)
# Enriched_NegativeGns <- enrichr(genes = NegativeGenes, databases = dbs)
# printEnrich(Enriched_PositiveGns, "PosGnsEnrichment.txt" , sep = "\t", columns = c(1:9))
# printEnrich(Enriched_NegativeGns, "NegGnsEnrichment.txt" , sep = "\t", columns = c(1:9))
# 
# Pos_GO_BP <- Enriched_PositiveGns["GO_Biological_Process_2018"]
# Pos_GO_BP <- Pos_GO_BP$GO_Biological_Process_2018
# Pos_GO_BP <- Pos_GO_BP[Pos_GO_BP$P.value <= 0.05, ]
# 
# Neg_GO_BP <- Enriched_NegativeGns["GO_Biological_Process_2018"]
# Neg_GO_BP <- Neg_GO_BP$GO_Biological_Process_2018
# Neg_GO_BP <- Neg_GO_BP[Neg_GO_BP$P.value <= 0.05, ]
# 
# 
# Pos_KEGG <- Enriched_PositiveGns["KEGG_2019_Human"]
# Pos_KEGG <- Pos_KEGG$KEGG_2019_Human
# Pos_KEGG <- Pos_KEGG[Pos_KEGG$P.value <= 0.05, ]
# 
# Neg_KEGG <- Enriched_NegativeGns["KEGG_2019_Human"]
# Neg_KEGG <- Neg_KEGG$KEGG_2019_Human
# Neg_KEGG <- Neg_KEGG[Neg_KEGG$P.value <= 0.05, ]

########################################################################
## Estimating the immune cell proportions based on the gene expression profiles
#immuno_meta <- immunoStatesMeta(metaObject = Prostate_metaanalysis)

## Correct the underlying gene expression data based on the cell proportions
#immuno_meta_corrected <- immunoStatesDecov(metaObject = Prostate_metaanalysis)

#########################################################################
## Using lincsCorrelate to look for a drug with a gene expression profile that reverses the metastasis profile
#lincProstate <- lincsCorrelate(metaObject = Prostate_metaanalysis, filterObject = New_filter, dataset = "CP", direction = "reverse")

####
# val_expr1 <- t(scale(t(val_expr1), center = T, scale = T))
# val_expr2 <- t(scale(t(val_expr2), center = T, scale = T))
# val_expr3 <- t(scale(t(val_expr3), center = T, scale = T))
# val_expr4 <- t(scale(t(val_expr4), center = T, scale = T))
# val_expr5 <- t(scale(t(val_expr5), center = T, scale = T))
# val_expr6 <- t(scale(t(val_expr6), center = T, scale = T))
# 
# 
# 
# AllValExpr <- list(val_expr1, val_expr2, val_expr3, val_expr6)
# CommonGns <- Reduce("intersect", lapply(AllValExpr, rownames))
# 
# AllValExpr <- mapply(x=AllValExpr, FUN=function(x, gns) {
#   x <- x[ gns ,]
# }, MoreArgs=list(gns=CommonGns))
# 
# 
# 
# AllValExpr <- do.call("cbind", AllValExpr)
# AllValExpr <- normalizeBetweenArrays(AllValExpr, method = "quantile")
# 
# 
# val_pheno1 <- as.data.frame(val_pheno1[,"Metastasis"])
# val_pheno2 <- as.data.frame(val_pheno2[,"Metastasis"])
# val_pheno3 <- as.data.frame(val_pheno3[,"Metastasis"])
# val_pheno4 <- as.data.frame(val_pheno4[,"Metastasis"])
# val_pheno5 <- as.data.frame(val_pheno5[,"Metastasis"])
# val_pheno6 <- as.data.frame(val_pheno6[,"Metastasis"])
# 
# colnames(val_pheno1) <- "Metastasis"
# colnames(val_pheno2) <- "Metastasis"
# colnames(val_pheno3) <- "Metastasis"
# colnames(val_pheno4) <- "Metastasis"
# colnames(val_pheno5) <- "Metastasis"
# colnames(val_pheno6) <- "Metastasis"
# 
# 
# AllValPheno <- rbind(val_pheno1, val_pheno2, val_pheno3, val_pheno6)
# 
# all(rownames(AllValPheno) == colnames(AllValExpr))
# rownames(AllValPheno) <- colnames(AllValExpr)
# 
# ValDataSet <- list()
# ValDataSet$expr <- AllValExpr
# ValDataSet$pheno <- AllValPheno
# ValDataSet$keys <- rownames(AllValExpr)
# ValDataSet$formattedName <- "ValDataSet"
# 
# ValDataSet <- classFunction(ValDataSet, column = "Metastasis", diseaseTerms = c("Mets"))
# 
# ##############################################
# 
# ### Now lets examine the performance of our filter on the independent data set
# 
# rocPlot(datasetObject = ValDataSet, filterObject = OncoDx_Filter)
# 
# X <- getSampleLevelGeneData(datasetObject = valDataSet1, geneNames = c(PositiveGenes))
# 
# X <- cleanUpPheno(valDataSet1)