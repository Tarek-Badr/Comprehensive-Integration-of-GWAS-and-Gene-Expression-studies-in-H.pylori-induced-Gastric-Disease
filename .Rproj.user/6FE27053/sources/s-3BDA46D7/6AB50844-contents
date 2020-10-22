###############################################################################################
## Mohamed Omar
## 21/07/2019
## Goal: Discovery and validation of a small gene signature that can predict the metastasis potential in primary prostate cancer
################################################################################################

## Clean work space
rm(list = ls())

## Set the working directory
#setwd("/Volumes/Macintosh/Research/Projects/H_Pylori")


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


#CanDataset3 <- getGEOData(c("GSE79973"))
#CanDataset3 <- CanDataset3$originalData$GSE79973

#CanDataset4 <- getGEOData(c("GSE19826"))
#CanDataset4 <- CanDataset4$originalData$GSE19826

#CanDataset5 <- getGEOData(c("GSE49051"))
#CanDataset5 <- CanDataset5$originalData$GSE49051

#CanDataset5 <- getGEOData(c("GSE13861"))
#CanDataset5 <- CanDataset5$originalData$GSE13861

# CanDataset6 <- getGEOData(c("GSE13911"))
# CanDataset6 <- CanDataset6$originalData$GSE13911

# CanDataset7 <- getGEOData(c("GSE29272"))
# CanDataset7 <- CanDataset7$originalData$GSE29272

#CanDataset8 <- getGEOData(c("GSE29998"))
#CanDataset8 <- CanDataset8$originalData$GSE29998

#CanDataset9 <- getGEOData(c("GSE31811"))
#CanDataset9 <- CanDataset9$originalData$GSE31811

# CanDataset10 <- getGEOData(c("GSE37023"))
# CanDataset11 <- CanDataset10$originalData$GSE37023_GPL96
# CanDataset12 <- CanDataset10$originalData$GSE37023_GPL97
# CanDataset10 <- CanDataset10$originalData$GSE37023_GPL2834

# CanDataset13 <- getGEOData(c("GSE44740"))
# CanDataset13 <- CanDataset13$originalData$GSE44740

# CanDataset14 <- getGEOData(c("GSE64951"))
# CanDataset14 <- CanDataset14$originalData$GSE64951

#CanDataset15 <- getGEOData(c("GSE81948"))
#CanDataset15 <- CanDataset15$originalData$GSE81948

#CanDataset16 <- getGEOData(c("GSE3438"))
#CanDataset16 <- CanDataset15$originalData$GSE3438

#save(CanDataset3, file = "./Data/CanDataset3.rda")
#save(CanDataset4, file = "./Data/CanDataset4.rda")
#save(CanDataset5, file = "./Data/CanDataset5.rda")
#save(CanDataset6, file = "./Data/CanDataset6.rda")
#save(CanDataset7, file = "./Data/CanDataset7.rda")
#save(CanDataset8, file = "./Data/CanDataset8.rda")
#save(CanDataset9, file = "./Data/CanDataset9.rda")
#save(CanDataset10, CanDataset11, CanDataset12, file = "./Data/CanDataset10_11_12.rda")
#save(CanDataset13, file = "./Data/CanDataset13.rda")
#save(CanDataset14, file = "./Data/CanDataset14.rda")

# save(Dataset1, Dataset2, Dataset3, Dataset4, file = "./Data/HPyloriData.rda")
# save(ValDataset1, ValDataset2, ValDataset3, ValDataset4, file = "./Data/HPyloriCellLineData.rda")
#save(CanDataset1, CanDataset3, CanDataset4, CanDataset5, CanDataset6, CanDataset7, CanDataset8, CanDataset9, CanDataset10, CanDataset11, CanDataset12, CanDataset13, CanDataset14, CanDataset15, file = "./Objs/CanDatasetsAll.rda")

#######
# GSE37023_GPL2834 : No Control observations
# GSE44740: From myofibroblast
# GSE64951 : From saliva
## Load Training data sets (4)
load("./Data/HPyloriData.rda")

## Load the validation datasets (cell lines)
load("./Data/HPyloriCellLineData.rda")

## Load the cancer dataset 1
#load("./Data/CanDataSet1.rda")

# load Cancer DataSet2 Pheno
canPheno2 <- read.delim("./Data/EMTAB1440/EMTAB1440Pheno.txt")

# load Cancer DataSet2 Expr
canExpr2 <- read.delim("./Data/EMTAB1440/ExpressionNormalized_v2.txt")


#load("./Data/CanDataset3.rda")
#load("./Data/CanDataset4.rda")
#load("./Data/CanDataset5.rda")

# Load all cancer datasets except canDataset2
load("./Objs/CanDatasetsAll.rda")


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

canPheno3 <- CanDataset3$pheno
canPheno4 <- CanDataset4$pheno

canPheno5 <- CanDataset5$pheno
canPheno6 <- CanDataset6$pheno
canPheno7 <- CanDataset7$pheno
canPheno8 <- CanDataset8$pheno
canPheno9 <- CanDataset9$pheno
canPheno10 <- CanDataset10$pheno
canPheno11 <- CanDataset11$pheno
canPheno12 <- CanDataset12$pheno
canPheno13 <- CanDataset13$pheno
canPheno14 <- CanDataset14$pheno
canPheno15 <- CanDataset15$pheno

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

canExpr3 <- CanDataset3$expr
canExpr4 <- CanDataset4$expr

canExpr5 <- CanDataset5$expr
canExpr6 <- CanDataset6$expr
canExpr7 <- CanDataset7$expr
canExpr8 <- CanDataset8$expr
canExpr9 <- CanDataset9$expr
canExpr10 <- CanDataset10$expr
canExpr11 <- CanDataset11$expr
canExpr12 <- CanDataset12$expr
canExpr13 <- CanDataset13$expr
canExpr14 <- CanDataset14$expr
canExpr15 <- CanDataset15$expr

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

###########################
## CancerExpr3
head(rownames(canExpr3))
rownames(canExpr3) <- CanDataset3$keys
canExpr3 <- canExpr3[!is.na(rownames(canExpr3)), ]
dim(canExpr3)

###########################
## CancerExpr4
head(rownames(canExpr4))
rownames(canExpr4) <- CanDataset4$keys
canExpr4 <- canExpr4[!is.na(rownames(canExpr4)), ]
dim(canExpr4)

###########################
## CancerExpr5
head(rownames(canExpr5))
rownames(canExpr5) <- CanDataset5$keys
canExpr5 <- canExpr5[!is.na(rownames(canExpr5)), ]
dim(canExpr5)

###########################
## CancerExpr6
head(rownames(canExpr6))
rownames(canExpr6) <- CanDataset6$keys
canExpr6 <- canExpr6[!is.na(rownames(canExpr6)), ]
dim(canExpr6)

###########################
## CancerExpr7
head(rownames(canExpr7))
rownames(canExpr7) <- CanDataset7$keys
canExpr7 <- canExpr7[!is.na(rownames(canExpr7)), ]
dim(canExpr7)

###########################
## CancerExpr8
head(rownames(canExpr8))
rownames(canExpr8) <- CanDataset8$keys
canExpr8 <- canExpr8[!is.na(rownames(canExpr8)), ]
dim(canExpr8)

###########################
## CancerExpr9
head(rownames(canExpr9))
rownames(canExpr9) <- CanDataset9$keys
canExpr9 <- canExpr9[!is.na(rownames(canExpr9)), ]
dim(canExpr9)

###########################
## CancerExpr10
head(rownames(canExpr10))
rownames(canExpr10) <- CanDataset10$keys
canExpr10 <- canExpr10[!is.na(rownames(canExpr10)), ]
dim(canExpr10)

###########################
## CancerExpr11
head(rownames(canExpr11))
rownames(canExpr11) <- CanDataset11$keys
canExpr11 <- canExpr11[!is.na(rownames(canExpr11)), ]
dim(canExpr11)

###########################
## CancerExpr12
head(rownames(canExpr12))
rownames(canExpr12) <- CanDataset12$keys
canExpr12 <- canExpr12[!is.na(rownames(canExpr12)), ]
dim(canExpr12)

###########################
## CancerExpr13
head(rownames(canExpr13))
rownames(canExpr13) <- CanDataset13$keys
canExpr13 <- canExpr13[!is.na(rownames(canExpr13)), ]
dim(canExpr13)

###########################
## CancerExpr14
head(rownames(canExpr14))
rownames(canExpr14) <- CanDataset14$keys
canExpr14 <- canExpr14[!is.na(rownames(canExpr14)), ]
dim(canExpr14)

###########################
## CancerExpr15
head(rownames(canExpr15))
rownames(canExpr15) <- CanDataset15$keys
canExpr15 <- canExpr15[!is.na(rownames(canExpr15)), ]
dim(canExpr15)

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


################## 
## Modify cancer pheno3

canPheno3$DiseaseStatus <- canPheno3$`tissue:ch1`
canPheno3$DiseaseStatus[canPheno3$DiseaseStatus == "gastric mucosa"] <- "control"
canPheno3$DiseaseStatus[canPheno3$DiseaseStatus == "gastric adenocarcinoma"] <- "case"

canPheno3$DiseaseStatus <- factor(canPheno3$DiseaseStatus, levels = c("control", "case"))
table(canPheno3$DiseaseStatus)

## Modify sample names to match sample names of valpheno3
all(rownames(canPheno3) == colnames(canExpr3))

## Finally, replace the expression and phenotype data in the dataset with the new modified versions
CanDataset3 <- list()
CanDataset3$expr <- canExpr3
CanDataset3$pheno <- canPheno3
CanDataset3$keys <- rownames(canExpr3)
CanDataset3$formattedName <- "GSE79973"

################## 
## Modify cancer pheno4

canPheno4$DiseaseStatus <- canPheno4$`tissue type:ch1`
canPheno4$DiseaseStatus[canPheno4$DiseaseStatus == "noncancer tissue"] <- "control"
canPheno4$DiseaseStatus[canPheno4$DiseaseStatus == "normal gastric tissue"] <- "control"
canPheno4$DiseaseStatus[canPheno4$DiseaseStatus == "gastric cancer tissue"] <- "case"

canPheno4$DiseaseStatus <- factor(canPheno4$DiseaseStatus, levels = c("control", "case"))
table(canPheno4$DiseaseStatus)

## Modify sample names to match sample names of valpheno4
all(rownames(canPheno4) == colnames(canExpr4))

## Finally, replace the expression and phenotype data in the dataset with the new modified versions
CanDataset4 <- list()
CanDataset4$expr <- canExpr4
CanDataset4$pheno <- canPheno4
CanDataset4$keys <- rownames(canExpr4)
CanDataset4$formattedName <- "GSE19826"

################## 
## Modify cancer pheno5

canPheno5$DiseaseStatus <- canPheno5$characteristics_ch1
canPheno5 <- canPheno5[!(canPheno5$characteristics_ch1 == "GIST"), ]

canPheno5$DiseaseStatus[canPheno5$DiseaseStatus == "normal surrounding gastric tissue"] <- "control"
canPheno5$DiseaseStatus[canPheno5$DiseaseStatus == "gastric adenocarcinoma"] <- "case"

canPheno5$DiseaseStatus <- factor(canPheno5$DiseaseStatus, levels = c("control", "case"))
table(canPheno5$DiseaseStatus)

## Modify sample names to match sample names of valpheno5
canExpr5 <- canExpr5[, colnames(canExpr5) %in% rownames(canPheno5)]
all(rownames(canPheno5) == colnames(canExpr5))

## Finally, replace the expression and phenotype data in the dataset with the new modified versions
CanDataset5 <- list()
CanDataset5$expr <- canExpr5
CanDataset5$pheno <- canPheno5
CanDataset5$keys <- rownames(canExpr5)
CanDataset5$formattedName <- "GSE13861"

################## 
## Modify cancer pheno6

canPheno6$DiseaseStatus <- canPheno6$characteristics_ch1.1

canPheno6$DiseaseStatus[canPheno6$DiseaseStatus == "NORMAL"] <- "control"
canPheno6$DiseaseStatus[canPheno6$DiseaseStatus == "TUMOR"] <- "case"

canPheno6$DiseaseStatus <- factor(canPheno6$DiseaseStatus, levels = c("control", "case"))
table(canPheno6$DiseaseStatus)

## Modify sample names to match sample names of valpheno6
all(rownames(canPheno6) == colnames(canExpr6))

## Finally, replace the expression and phenotype data in the dataset with the new modified versions
CanDataset6$expr <- canExpr6
CanDataset6$pheno <- canPheno6
CanDataset6$keys <- rownames(canExpr6)

################## 
## Modify cancer pheno7

canPheno7$DiseaseStatus <- canPheno7$source_name_ch1

canPheno7$DiseaseStatus[canPheno7$DiseaseStatus == "adjacent normal tissue"] <- "control"
canPheno7$DiseaseStatus[canPheno7$DiseaseStatus == "tumor tissue"] <- "case"

canPheno7$DiseaseStatus <- factor(canPheno7$DiseaseStatus, levels = c("control", "case"))
table(canPheno7$DiseaseStatus)

## Modify sample names to match sample names of valpheno7
all(rownames(canPheno7) == colnames(canExpr7))

## Finally, replace the expression and phenotype data in the dataset with the new modified versions
CanDataset7$expr <- canExpr7
CanDataset7$pheno <- canPheno7
CanDataset7$keys <- rownames(canExpr7)

################## 
## Modify cancer pheno8

canPheno8$DiseaseStatus <- canPheno8$`tissue:ch1`

canPheno8$DiseaseStatus[canPheno8$DiseaseStatus == "Normal"] <- "control"
canPheno8$DiseaseStatus[canPheno8$DiseaseStatus == "Tumor"] <- "case"

canPheno8$DiseaseStatus <- factor(canPheno8$DiseaseStatus, levels = c("control", "case"))
table(canPheno8$DiseaseStatus)

## Modify sample names to match sample names of valpheno8
all(rownames(canPheno8) == colnames(canExpr8))

## Finally, replace the expression and phenotype data in the dataset with the new modified versions
CanDataset8$expr <- canExpr8
CanDataset8$pheno <- canPheno8
CanDataset8$keys <- rownames(canExpr8)

################## 
## Modify cancer pheno9

canPheno9$DiseaseStatus <- canPheno9$`status:ch1`

canPheno9$DiseaseStatus[canPheno9$DiseaseStatus == "normal"] <- "control"
canPheno9$DiseaseStatus[canPheno9$DiseaseStatus == "tumor"] <- "case"

canPheno9$DiseaseStatus <- factor(canPheno9$DiseaseStatus, levels = c("control", "case"))
table(canPheno9$DiseaseStatus)

## Modify sample names to match sample names of valpheno9
all(rownames(canPheno9) == colnames(canExpr9))

## Finally, replace the expression and phenotype data in the dataset with the new modified versions
CanDataset9$expr <- canExpr9
CanDataset9$pheno <- canPheno9
CanDataset9$keys <- rownames(canExpr9)

################## 
## Modify cancer pheno11

canPheno11$DiseaseStatus <- canPheno11$`tissue:ch1`

canPheno11$DiseaseStatus[canPheno11$DiseaseStatus %in% c("Gastric non-malignant", "gastric non-malignant sample")] <- "control"
canPheno11$DiseaseStatus[canPheno11$DiseaseStatus %in% c("gastric tumor", "Gastric tumor")] <- "case"

canPheno11$DiseaseStatus <- factor(canPheno11$DiseaseStatus, levels = c("control", "case"))
table(canPheno11$DiseaseStatus)

## Modify sample names to match sample names of valpheno11
all(rownames(canPheno11) == colnames(canExpr11))

## Finally, replace the expression and phenotype data in the dataset with the new modified versions
CanDataset11$expr <- canExpr11
CanDataset11$pheno <- canPheno11
CanDataset11$keys <- rownames(canExpr11)

########################
## Modify cancer pheno12

canPheno12$DiseaseStatus <- canPheno12$`tissue:ch1`

canPheno12$DiseaseStatus[canPheno12$DiseaseStatus ==  "Gastric non-malignant"] <- "control"
canPheno12$DiseaseStatus[canPheno12$DiseaseStatus == "Gastric tumor"] <- "case"

canPheno12$DiseaseStatus <- factor(canPheno12$DiseaseStatus, levels = c("control", "case"))
table(canPheno12$DiseaseStatus)

## Modify sample names to match sample names of valpheno12
all(rownames(canPheno12) == colnames(canExpr12))

## Finally, replace the expression and phenotype data in the dataset with the new modified versions
CanDataset12$expr <- canExpr12
CanDataset12$pheno <- canPheno12
CanDataset12$keys <- rownames(canExpr12)

########################
## Modify cancer pheno13

canPheno13$DiseaseStatus <- canPheno13$`myofibroblast source:ch1`

canPheno13 <- canPheno13[!(canPheno13$`myofibroblast source:ch1` == "pernicious anaemia patient"), ]
canPheno13$DiseaseStatus[canPheno13$DiseaseStatus %in%  c("normal tissue", "adjacent tissue")] <- "control"
canPheno13$DiseaseStatus[canPheno13$DiseaseStatus == "cancer associated"] <- "case"

canPheno13$DiseaseStatus <- factor(canPheno13$DiseaseStatus, levels = c("control", "case"))
table(canPheno13$DiseaseStatus)

## Modify sample names to match sample names of valpheno13
canExpr13 <- canExpr13[, colnames(canExpr13) %in% rownames(canPheno13)]
all(rownames(canPheno13) == colnames(canExpr13))

## Finally, replace the expression and phenotype data in the dataset with the new modified versions
CanDataset13$expr <- canExpr13
CanDataset13$pheno <- canPheno13
CanDataset13$keys <- rownames(canExpr13)

########################
## Modify cancer pheno14

canPheno14$DiseaseStatus <- canPheno14$`disease state:ch1`

canPheno14$DiseaseStatus[canPheno14$DiseaseStatus == "normal"] <- "control"
canPheno14$DiseaseStatus[canPheno14$DiseaseStatus == "gastric cancer"] <- "case"

canPheno14$DiseaseStatus <- factor(canPheno14$DiseaseStatus, levels = c("control", "case"))
table(canPheno14$DiseaseStatus)

## Modify sample names to match sample names of valpheno14
all(rownames(canPheno14) == colnames(canExpr14))

## Finally, replace the expression and phenotype data in the dataset with the new modified versions
CanDataset14$expr <- canExpr14
CanDataset14$pheno <- canPheno14
CanDataset14$keys <- rownames(canExpr14)

########################
## Modify cancer pheno15

canPheno15$DiseaseStatus <- canPheno15$`tissue:ch1`

canPheno15$DiseaseStatus[canPheno15$DiseaseStatus == "normal"] <- "control"
canPheno15$DiseaseStatus[canPheno15$DiseaseStatus == "tumor"] <- "case"

canPheno15$DiseaseStatus <- factor(canPheno15$DiseaseStatus, levels = c("control", "case"))
table(canPheno15$DiseaseStatus)

## Modify sample names to match sample names of valpheno15
all(rownames(canPheno15) == colnames(canExpr15))

## Finally, replace the expression and phenotype data in the dataset with the new modified versions
CanDataset15$expr <- canExpr15
CanDataset15$pheno <- canPheno15
CanDataset15$keys <- rownames(canExpr15)

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
CanDataset3 <- classFunction(CanDataset3, column = "DiseaseStatus", diseaseTerms = c("case"))
CanDataset4 <- classFunction(CanDataset4, column = "DiseaseStatus", diseaseTerms = c("case"))

CanDataset5 <- classFunction(CanDataset5, column = "DiseaseStatus", diseaseTerms = c("case"))
CanDataset6 <- classFunction(CanDataset6, column = "DiseaseStatus", diseaseTerms = c("case"))
CanDataset7 <- classFunction(CanDataset7, column = "DiseaseStatus", diseaseTerms = c("case"))
CanDataset8 <- classFunction(CanDataset8, column = "DiseaseStatus", diseaseTerms = c("case"))
CanDataset9 <- classFunction(CanDataset9, column = "DiseaseStatus", diseaseTerms = c("case"))
#CanDataset10 <- classFunction(CanDataset10, column = "DiseaseStatus", diseaseTerms = c("case"))
CanDataset11 <- classFunction(CanDataset11, column = "DiseaseStatus", diseaseTerms = c("case"))
CanDataset12 <- classFunction(CanDataset12, column = "DiseaseStatus", diseaseTerms = c("case"))
CanDataset13 <- classFunction(CanDataset13, column = "DiseaseStatus", diseaseTerms = c("case"))
CanDataset14 <- classFunction(CanDataset14, column = "DiseaseStatus", diseaseTerms = c("case"))
CanDataset15 <- classFunction(CanDataset15, column = "DiseaseStatus", diseaseTerms = c("case"))

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
#save(filter, file = "./Objs/filter.rda")

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

## Gene names
PositiveGenes <- New_filter2$posGeneNames
PositiveGenes
NegativeGenes <- New_filter2$negGeneNames
NegativeGenes

PosGenesOrdered <- New_filter2_summary$pos
PosGenesOrdered <- PosGenesOrdered[order(PosGenesOrdered$effectSize, decreasing = T), ]
#PosGenesOrdered <- rownames(PosGenesOrdered)

NegGenesOrdered <- New_filter2_summary$neg
NegGenesOrdered <- NegGenesOrdered[order(NegGenesOrdered$effectSize, decreasing = F), ]
#NegGenesOrdered <- rownames(NegGenesOrdered)

## Save the tables of positive and negative genes
write.table(PosGenesOrdered, file = "./Objs/NewFilter2_Positive_genes.csv", quote = T, sep = "\t", col.names = T, row.names = T)
write.table(NegGenesOrdered, file = "./Objs/NewFilter2_Negative_genes.csv", quote = T, sep = "\t", col.names = T, row.names = T)

## Create a summary ROC curve (Training data sets)
set.seed(333)
png(filename = "./Figs/Pooled_ROC_TrainingDatasets_NewFilter2.png", width = 2000, height = 2000, res = 300)
pooledROCPlot(metaObject = HPylori_metaanalysis, filterObject = New_filter2)
dev.off()

## Effect size visualization in each study
# For the up-regulated genes
pdf(file = "./Figs/Genes_effect_size.pdf", title = "Effect size of up-regulated genes across the discovery data sets", width = 10, height = 15)
par(mfrow=c(5,2))
# Up
forestPlot(metaObject = HPylori_metaanalysis, geneName = "SERPINA3", textColor = "black") 
forestPlot(metaObject = HPylori_metaanalysis, geneName = "XK", textColor = "black")

forestPlot(metaObject = HPylori_metaanalysis, geneName = "CASP1", textColor = "black")
forestPlot(metaObject = HPylori_metaanalysis, geneName = "TPST1", textColor = "black")

forestPlot(metaObject = HPylori_metaanalysis, geneName = "IFNGR1", textColor = "black")
forestPlot(metaObject = HPylori_metaanalysis, geneName = "GUCA2B", textColor = "black")

forestPlot(metaObject = HPylori_metaanalysis, geneName = "TLR8", textColor = "black")
forestPlot(metaObject = HPylori_metaanalysis, geneName = "TDRD3", textColor = "black")

forestPlot(metaObject = HPylori_metaanalysis, geneName = "TNFRSF10B", textColor = "black")
# Down 
forestPlot(metaObject = HPylori_metaanalysis, geneName = "LAPTM4A", textColor = "black")
dev.off()

######
## Tiff more width
tiff(filename = "./Figs/Genes_effect_size.tiff", width = 3000, height = 1500, pointsize = 30)
par(mfrow=c(2,5))
# Up
forestPlot(metaObject = HPylori_metaanalysis, geneName = "SERPINA3", textColor = "black") 
forestPlot(metaObject = HPylori_metaanalysis, geneName = "CASP1", textColor = "black")
forestPlot(metaObject = HPylori_metaanalysis, geneName = "IFNGR1", textColor = "black")
forestPlot(metaObject = HPylori_metaanalysis, geneName = "TLR8", textColor = "black")
forestPlot(metaObject = HPylori_metaanalysis, geneName = "TNFRSF10B", textColor = "black")
# Down
forestPlot(metaObject = HPylori_metaanalysis, geneName = "XK", textColor = "black")
forestPlot(metaObject = HPylori_metaanalysis, geneName = "TPST1", textColor = "black")
forestPlot(metaObject = HPylori_metaanalysis, geneName = "GUCA2B", textColor = "black")
forestPlot(metaObject = HPylori_metaanalysis, geneName = "TDRD3", textColor = "black")
forestPlot(metaObject = HPylori_metaanalysis, geneName = "LAPTM4A", textColor = "black")
dev.off()

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

## Cancer dataset3
png(filename = "./Figs/ROC_Test_Cancer2_Signature2.png", width = 2000, height = 2000, res = 300)
rocPlot(datasetObject = CanDataset3, filterObject = New_filter)
dev.off()

## Cancer dataset4
png(filename = "./Figs/ROC_Test_Cancer2_Signature2.png", width = 2000, height = 2000, res = 300)
rocPlot(datasetObject = CanDataset4, filterObject = New_filter2)
dev.off()

## Cancer dataset5
png(filename = "./Figs/ROC_Test_Cancer2_Signature2.png", width = 2000, height = 2000, res = 300)
rocPlot(datasetObject = CanDataset5, filterObject = New_filter2)
dev.off()

# Cancer dataset6
rocPlot(datasetObject = CanDataset6, filterObject = New_filter2)

# Cancer dataset7
rocPlot(datasetObject = CanDataset7, filterObject = New_filter2)

# Cancer dataset8
rocPlot(datasetObject = CanDataset8, filterObject = New_filter2)

# Cancer dataset9
rocPlot(datasetObject = CanDataset9, filterObject = New_filter2)

# Cancer dataset11
rocPlot(datasetObject = CanDataset11, filterObject = New_filter2)

# Cancer dataset12
rocPlot(datasetObject = CanDataset12, filterObject = New_filter2)

# Cancer dataset13  # myofibroblast
rocPlot(datasetObject = CanDataset13, filterObject = New_filter2)

# Cancer dataset14 # Saliva
rocPlot(datasetObject = CanDataset14, filterObject = New_filter2)

# Cancer dataset15
rocPlot(datasetObject = CanDataset15, filterObject = New_filter2)

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

# cancer dataset3
png(filename = "./Figs/PRC_Test6_Signature2.png", width = 2000, height = 2000, res = 300)
prcPlot(datasetObject = CanDataset3, filterObject = New_filter)
dev.off()

# cancer dataset4
png(filename = "./Figs/PRC_Test6_Signature2.png", width = 2000, height = 2000, res = 300)
prcPlot(datasetObject = CanDataset4, filterObject = New_filter)
dev.off()

# cancer dataset5
png(filename = "./Figs/PRC_Test6_Signature2.png", width = 2000, height = 2000, res = 300)
prcPlot(datasetObject = CanDataset5, filterObject = New_filter2)
dev.off()

# cancer dataset6
prcPlot(datasetObject = CanDataset6, filterObject = New_filter2)

# cancer dataset7
prcPlot(datasetObject = CanDataset7, filterObject = New_filter2)

# cancer dataset8
prcPlot(datasetObject = CanDataset8, filterObject = New_filter2)

# cancer dataset9
prcPlot(datasetObject = CanDataset9, filterObject = New_filter2)

# cancer dataset11
prcPlot(datasetObject = CanDataset11, filterObject = New_filter2)

# cancer dataset12
prcPlot(datasetObject = CanDataset12, filterObject = New_filter2)

# cancer dataset13
prcPlot(datasetObject = CanDataset13, filterObject = New_filter2)

# cancer dataset14
prcPlot(datasetObject = CanDataset14, filterObject = New_filter2)

# cancer dataset15
prcPlot(datasetObject = CanDataset15, filterObject = New_filter2)


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

# cancer dataset3
png(filename = "./Figs/ViolinPlot_Test_Cancer2_Signature2.png", width = 2000, height = 2000, res = 300)
violinPlot(filterObject = New_filter2, datasetObject = CanDataset3, labelColumn = "DiseaseStatus")
dev.off()

# cancer dataset4
png(filename = "./Figs/ViolinPlot_Test_Cancer2_Signature2.png", width = 2000, height = 2000, res = 300)
violinPlot(filterObject = New_filter2, datasetObject = CanDataset4, labelColumn = "DiseaseStatus")
dev.off()

# cancer dataset5
png(filename = "./Figs/ViolinPlot_Test_Cancer2_Signature2.png", width = 2000, height = 2000, res = 300)
violinPlot(filterObject = New_filter2, datasetObject = CanDataset5, labelColumn = "DiseaseStatus")
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

# cancer dataset3
canPheno3$score <- calculateScore(filterObject = New_filter2, datasetObject = CanDataset3)

# cancer dataset4
canPheno4$score <- calculateScore(filterObject = New_filter2, datasetObject = CanDataset4)

# cancer dataset5
canPheno5$score <- calculateScore(filterObject = New_filter2, datasetObject = CanDataset5)

# cancer dataset6
canPheno6$score <- calculateScore(filterObject = New_filter2, datasetObject = CanDataset6)

# cancer dataset7
canPheno7$score <- calculateScore(filterObject = New_filter2, datasetObject = CanDataset7)

# cancer dataset8
canPheno8$score <- calculateScore(filterObject = New_filter2, datasetObject = CanDataset8)

# cancer dataset9
canPheno9$score <- calculateScore(filterObject = New_filter2, datasetObject = CanDataset9)

# cancer dataset11
canPheno11$score <- calculateScore(filterObject = New_filter2, datasetObject = CanDataset11)

# cancer dataset12
canPheno12$score <- calculateScore(filterObject = New_filter2, datasetObject = CanDataset12)

# cancer dataset15
canPheno15$score <- calculateScore(filterObject = New_filter2, datasetObject = CanDataset15)

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


# cancer dataset3
Test_Predictions7 <- ifelse(canPheno3$score >= Mean_Thr, "case", "control")
confusionMatrix(as.factor(Test_Predictions7), canPheno3$DiseaseStatus, positive = "case")

mcc(preds = as.factor(Test_Predictions7), actuals = canPheno3$DiseaseStatus)


# cancer dataset4
Test_Predictions8 <- ifelse(canPheno4$score >= Mean_Thr, "case", "control")
confusionMatrix(as.factor(Test_Predictions8), canPheno4$DiseaseStatus, positive = "case")

mcc(preds = as.factor(Test_Predictions8), actuals = canPheno4$DiseaseStatus)


# cancer dataset5
Test_Predictions9 <- ifelse(canPheno5$score >= Mean_Thr, "case", "control")
confusionMatrix(as.factor(Test_Predictions9), canPheno5$DiseaseStatus, positive = "case")

mcc(preds = as.factor(Test_Predictions9), actuals = canPheno5$DiseaseStatus)


############################################
#############################################
## Big plot1: Discovery performance + validation in the 11 cancer datasets using intersection Signature

# Performance in the 1st cancer dataset
sscurves_CanDataSet1 <- evalmod(scores = canPheno1$score, labels = canPheno1$DiseaseStatus)

# Performance in the 2nd cancer dataset
sscurves_CanDataSet2  <- evalmod(scores = canPheno2$score, labels = canPheno2$DiseaseStatus)

# Performance in the 3rd cancer dataset
sscurves_CanDataSet3  <- evalmod(scores = canPheno3$score, labels = canPheno3$DiseaseStatus)

# Performance in the 4th cancer dataset
sscurves_CanDataSet4  <- evalmod(scores = canPheno4$score, labels = canPheno4$DiseaseStatus)

# Performance in the 5th cancer dataset
sscurves_CanDataSet5  <- evalmod(scores = canPheno5$score, labels = canPheno5$DiseaseStatus)

# Performance in the 6th cancer dataset
sscurves_CanDataSet6  <- evalmod(scores = canPheno6$score, labels = canPheno6$DiseaseStatus)

# Performance in the 7th cancer dataset
sscurves_CanDataSet7  <- evalmod(scores = canPheno7$score, labels = canPheno7$DiseaseStatus)

# Performance in the 8th cancer dataset
sscurves_CanDataSet8  <- evalmod(scores = canPheno8$score, labels = canPheno8$DiseaseStatus)

# Performance in the 9th cancer dataset
sscurves_CanDataSet9  <- evalmod(scores = canPheno9$score, labels = canPheno9$DiseaseStatus)

# Performance in the 11th cancer dataset
sscurves_CanDataSet11  <- evalmod(scores = canPheno11$score, labels = canPheno11$DiseaseStatus)

# Performance in the 12th cancer dataset
sscurves_CanDataSet12  <- evalmod(scores = canPheno12$score, labels = canPheno12$DiseaseStatus)

# Performance in the 15th cancer dataset
sscurves_CanDataSet15  <- evalmod(scores = canPheno15$score, labels = canPheno15$DiseaseStatus)

My_Theme = theme(
  axis.title.x = element_text(size = 12),
  axis.text.x = element_text(size = 12),
  axis.title.y = element_text(size = 12),
  axis.text.y = element_text(size = 12),
  plot.title = element_text(size=18)
)

My_Theme2 = theme(
  axis.title.x = element_text(size = 12),
  axis.text.x = element_text(size = 12),
  axis.title.y = element_text(size = 12),
  axis.text.y = element_text(size = 12),
  plot.title = element_text(size=15)
)

# Pooled AUC in the discovery datasets
Discovery <- pooledROCPlot(metaObject = HPylori_metaanalysis, filterObject = New_filter2)

# Violin, AUC, and AUPRC in the 1st cancer dataset
Can1Violin <- violinPlot(filterObject = New_filter2, datasetObject = CanDataset1, labelColumn = "DiseaseStatus") + My_Theme
Can1ROC <- autoplot(sscurves_CanDataSet1, curvetype = c("ROC")) + labs(title = "ROC") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.812"), size =3)
Can1PRC <- autoplot(sscurves_CanDataSet1, curvetype = c("PRC")) +  labs(title = "PRC") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.83"), size = 3)

# Violin, AUC, and AUPRC in the 2nd cancer dataset
Can2Violin <- violinPlot(filterObject = New_filter2, datasetObject = CanDataset2, labelColumn = "DiseaseStatus")+My_Theme
Can2ROC <- autoplot(sscurves_CanDataSet2, curvetype = c("ROC")) + labs(title = "ROC") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.82"), size = 3)
Can2PRC <- autoplot(sscurves_CanDataSet2, curvetype = c("PRC")) + labs(title = "PRC") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.84"), size = 3)

# Violin, AUC, and AUPRC in the 3rd cancer dataset
Can3Violin <- violinPlot(filterObject = New_filter2, datasetObject = CanDataset3, labelColumn = "DiseaseStatus")+My_Theme
Can3ROC <- autoplot(sscurves_CanDataSet3, curvetype = c("ROC")) + labs(title = "ROC") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.93"), size = 3)
Can3PRC <- autoplot(sscurves_CanDataSet3, curvetype = c("PRC")) + labs(title = "PRC") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.94"), size = 3)

# Violin, AUC, and AUPRC in the 4th cancer dataset
Can4Violin <- violinPlot(filterObject = New_filter2, datasetObject = CanDataset4, labelColumn = "DiseaseStatus")+My_Theme
Can4ROC <- autoplot(sscurves_CanDataSet4, curvetype = c("ROC")) + labs(title = "ROC") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.81"), size = 3)
Can4PRC <- autoplot(sscurves_CanDataSet4, curvetype = c("PRC")) + labs(title = "PRC") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.73"), size = 3)

# Violin, AUC, and AUPRC in the 5th cancer dataset
Can5Violin <- violinPlot(filterObject = New_filter2, datasetObject = CanDataset5, labelColumn = "DiseaseStatus")+My_Theme
Can5ROC <- autoplot(sscurves_CanDataSet5, curvetype = c("ROC")) + labs(title = "ROC") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.71"), size = 3)
Can5PRC <- autoplot(sscurves_CanDataSet5, curvetype = c("PRC")) + labs(title = "PRC") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.88"), size = 3)

# Violin, AUC, and AUPRC in the 6th cancer dataset
Can6Violin <- violinPlot(filterObject = New_filter2, datasetObject = CanDataset6, labelColumn = "DiseaseStatus")+My_Theme2
Can6ROC <- autoplot(sscurves_CanDataSet6, curvetype = c("ROC")) + labs(title = "ROC") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.86"), size = 4)
Can6PRC <- autoplot(sscurves_CanDataSet6, curvetype = c("PRC")) + labs(title = "PRC") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.88"), size = 4)

# Violin, AUC, and AUPRC in the 7th cancer dataset
Can7Violin <- violinPlot(filterObject = New_filter2, datasetObject = CanDataset7, labelColumn = "DiseaseStatus")+My_Theme2
Can7ROC <- autoplot(sscurves_CanDataSet7, curvetype = c("ROC")) + labs(title = "ROC") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.90"), size = 4)
Can7PRC <- autoplot(sscurves_CanDataSet7, curvetype = c("PRC")) + labs(title = "PRC") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.90"), size = 4)

# Violin, AUC, and AUPRC in the 8th cancer dataset
Can8Violin <- violinPlot(filterObject = New_filter2, datasetObject = CanDataset8, labelColumn = "DiseaseStatus")+My_Theme2
Can8ROC <- autoplot(sscurves_CanDataSet8, curvetype = c("ROC")) + labs(title = "ROC") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.86"), size = 4)
Can8PRC <- autoplot(sscurves_CanDataSet8, curvetype = c("PRC")) + labs(title = "PRC") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.85"), size = 4)

# Violin, AUC, and AUPRC in the 9th cancer dataset
Can9Violin <- violinPlot(filterObject = New_filter2, datasetObject = CanDataset9, labelColumn = "DiseaseStatus")+My_Theme2
Can9ROC <- autoplot(sscurves_CanDataSet9, curvetype = c("ROC")) + labs(title = "ROC") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.72"), size = 4)
Can9PRC <- autoplot(sscurves_CanDataSet9, curvetype = c("PRC")) + labs(title = "PRC") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.73"), size = 4)

# Violin, AUC, and AUPRC in the 11th cancer dataset
Can11Violin <- violinPlot(filterObject = New_filter2, datasetObject = CanDataset11, labelColumn = "DiseaseStatus")+My_Theme2
Can11ROC <- autoplot(sscurves_CanDataSet11, curvetype = c("ROC")) + labs(title = "ROC") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.93"), size = 4)
Can11PRC <- autoplot(sscurves_CanDataSet11, curvetype = c("PRC")) + labs(title = "PRC") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.97"), size = 4)

# Violin, AUC, and AUPRC in the 12th cancer dataset
Can12Violin <- violinPlot(filterObject = New_filter2, datasetObject = CanDataset12, labelColumn = "DiseaseStatus")+My_Theme2
Can12ROC <- autoplot(sscurves_CanDataSet12, curvetype = c("ROC")) + labs(title = "ROC") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.76"), size = 4)
Can12PRC <- autoplot(sscurves_CanDataSet12, curvetype = c("PRC")) + labs(title = "PRC") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.72"), size = 4)

# Violin, AUC, and AUPRC in the 15th cancer dataset
Can15Violin <- violinPlot(filterObject = New_filter2, datasetObject = CanDataset15, labelColumn = "DiseaseStatus")+My_Theme2
Can15ROC <- autoplot(sscurves_CanDataSet15, curvetype = c("ROC")) + labs(title = "ROC") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.93"), size = 4)
Can15PRC <- autoplot(sscurves_CanDataSet15, curvetype = c("PRC")) + labs(title = "PRC") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.98"), size = 4)

## Plot the discovery and the first 5 cancer datasets
tiff(filename = "./Figs/Intersec_Sig_5Cancer.tiff", width = 1500, height = 1000)
Discovery +    ((Can1Violin / Can1ROC / Can1PRC + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 12))) | 
                  (Can2Violin / Can2ROC / Can2PRC + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 12))) | 
                  (Can3Violin / Can3ROC / Can3PRC + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 12))) | 
                  (Can4Violin / Can4ROC / Can4PRC + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 12))) |
                  (Can5Violin / Can5ROC / Can5PRC + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 12)))
) +
  plot_layout(widths = c(0.4, 1)) + 
  plot_annotation(
    title = 'Refined Signature performance in the discovery datasets and the fist 5 gastric adenocarcinoma datasets',
    tag_levels = c('A', '1'),
    theme = theme(plot.title = element_text(size = 17, face = "bold"))
  )
dev.off()

########
## Plot the next 7 cancer datasets
tiff(filename = "./Figs/Intersec_Sig_7MoreCancer.tiff", width = 1500, height = 1000)
((Can6Violin / Can6ROC / Can6PRC + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 12))) | 
  (Can7Violin / Can7ROC / Can7PRC + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 12))) | 
  (Can8Violin / Can8ROC / Can8PRC + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 12))) |
  (Can9Violin / Can9ROC / Can9PRC + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 12))) |
  (Can11Violin / Can11ROC / Can11PRC + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 12))) |
  (Can12Violin / Can12ROC / Can12PRC + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 12))) | 
  (Can15Violin / Can15ROC / Can15PRC + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 12)))) +
  #plot_layout(widths = c(0.5, 1)) + 
  plot_annotation(
    title = 'Refined Signature performance in more 7 gastric adenocarcinoma datasets',
    tag_levels = c('A', '1'),
    theme = theme(plot.title = element_text(size = 17, face = "bold"))
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

# Performance in the 3rd cancer dataset
sscurves_CanDataSet3  <- evalmod(scores = canPheno3$score, labels = canPheno3$DiseaseStatus)

# Performance in the 4th cancer dataset
sscurves_CanDataSet4  <- evalmod(scores = canPheno4$score, labels = canPheno4$DiseaseStatus)

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
P15 <- autoplot(sscurves_CanDataSet1, curvetype = c("ROC")) + labs(title = "ROC") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.82"), size =4)
P16 <- autoplot(sscurves_CanDataSet1, curvetype = c("PRC")) +  labs(title = "PRC") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.81"), size = 4)

# Violin, AUC, and AUPRC in the 2nd cancer dataset
P17 <- violinPlot(filterObject = New_filter, datasetObject = CanDataset2, labelColumn = "DiseaseStatus") + My_Theme
P18 <- autoplot(sscurves_CanDataSet2, curvetype = c("ROC")) + labs(title = "ROC") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.82"), size = 4)
P19 <- autoplot(sscurves_CanDataSet2, curvetype = c("PRC")) + labs(title = "PRC") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.84"), size = 4)

# Violin, AUC, and AUPRC in the 3rd cancer dataset
P20 <- violinPlot(filterObject = New_filter, datasetObject = CanDataset3, labelColumn = "DiseaseStatus") + My_Theme
P21 <- autoplot(sscurves_CanDataSet3, curvetype = c("ROC")) + labs(title = "ROC") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.84"), size = 4)
P22 <- autoplot(sscurves_CanDataSet3, curvetype = c("PRC")) + labs(title = "PRC") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.86"), size = 4)

# Violin, AUC, and AUPRC in the 4th cancer dataset
P23 <- violinPlot(filterObject = New_filter, datasetObject = CanDataset4, labelColumn = "DiseaseStatus") + My_Theme
P24 <- autoplot(sscurves_CanDataSet4, curvetype = c("ROC")) + labs(title = "ROC") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.81"), size = 4)
P25 <- autoplot(sscurves_CanDataSet4, curvetype = c("PRC")) + labs(title = "PRC") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.72"), size = 4)

## Plot ALL
tiff(filename = "./Figs/Original_Sig_CancerDatasets.tiff", width = 1500, height = 1000)
P1 +  ((P14 / P15 / P16 + plot_layout(tag_level = "new")& theme(plot.tag = element_text(size = 20))) | (P17 / P18 / P19 +plot_layout(tag_level = "new")& theme(plot.tag = element_text(size = 20))) | (P20 / P21 / P22 +plot_layout(tag_level = "new")& theme(plot.tag = element_text(size = 20))) | (P23 / P24 / P25 +plot_layout(tag_level = "new")& theme(plot.tag = element_text(size = 20)))) + 
  plot_layout(widths = c(0.6, 2)) + 
  plot_annotation(
    title = 'Original Signature performance in the discovery datasets and the 4 gastric adenocarcinoma datasets',
    tag_levels = c('A', '1'),
    theme = theme(plot.title = element_text(size = 20, face = "bold"))
  )
dev.off() 


tiff(filename = "./Figs/Original_SigCellLines.tiff", width = 1500, height = 1000)
((P2 / P3 / P4 + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 20))) | (P5 / P6 / P7+ plot_layout(tag_level = "new")& theme(plot.tag = element_text(size = 20))) |  (P8 / P9 / P10 + plot_layout(tag_level = "new")& theme(plot.tag = element_text(size = 20))) | (P11 / P12 / P13 + plot_layout(tag_level = "new")& theme(plot.tag = element_text(size = 20))))  + 
  plot_annotation(
  title = 'Original Signature performance in the cell lines datasets',
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

##########################################################################
dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018", "KEGG_2019_Human", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018")
Enriched_PositiveGns <- enrichr(genes = PositiveGenes, databases = dbs)
Enriched_NegativeGns <- enrichr(genes = NegativeGenes, databases = dbs)
printEnrich(Enriched_PositiveGns, "PosGnsEnrichment.txt" , sep = "\t", columns = c(1:9))
printEnrich(Enriched_NegativeGns, "NegGnsEnrichment.txt" , sep = "\t", columns = c(1:9))
# 
Pos_GO_BP <- Enriched_PositiveGns["GO_Biological_Process_2018"]
Pos_GO_BP <- Pos_GO_BP$GO_Biological_Process_2018
Pos_GO_BP <- Pos_GO_BP[Pos_GO_BP$P.value <= 0.05, ]

write.table(Pos_GO_BP, file = "./Objs/UpRegulated_GO_BPs.csv", quote = T, sep = "\t", col.names = T, row.names = T)

##############
Neg_GO_BP <- Enriched_NegativeGns["GO_Biological_Process_2018"]
Neg_GO_BP <- Neg_GO_BP$GO_Biological_Process_2018
Neg_GO_BP <- Neg_GO_BP[Neg_GO_BP$P.value <= 0.05, ]

write.table(Neg_GO_BP, file = "./Objs/DownRegulated_GO_BPs.csv", quote = T, sep = "\t", col.names = T, row.names = T)

################
Pos_GO_MFs <- Enriched_PositiveGns["GO_Molecular_Function_2018"]
Pos_GO_MFs <- Pos_GO_MFs$GO_Molecular_Function_2018
Pos_GO_MFs <- Pos_GO_MFs[Pos_GO_MFs$P.value <= 0.05, ]

write.table(Pos_GO_MFs, file = "./Objs/UpRegulated_GO_MFs.csv", quote = T, sep = "\t", col.names = T, row.names = T)

##################
Neg_GO_MFs <- Enriched_NegativeGns["GO_Molecular_Function_2018"]
Neg_GO_MFs <- Neg_GO_MFs$GO_Molecular_Function_2018
Neg_GO_MFs <- Neg_GO_MFs[Neg_GO_MFs$P.value <= 0.05, ]

write.table(Neg_GO_MFs, file = "./Objs/DownRegulated_GO_MFs.csv", quote = T, sep = "\t", col.names = T, row.names = T)

################
Pos_GO_CCs <- Enriched_PositiveGns["GO_Cellular_Component_2018"]
Pos_GO_CCs <- Pos_GO_CCs$GO_Cellular_Component_2018
Pos_GO_CCs <- Pos_GO_CCs[Pos_GO_CCs$P.value <= 0.05, ]

write.table(Pos_GO_CCs, file = "./Objs/UpRegulated_GO_CCs.csv", quote = T, sep = "\t", col.names = T, row.names = T)

##################
Neg_GO_CCs <- Enriched_NegativeGns["GO_Cellular_Component_2018"]
Neg_GO_CCs <- Neg_GO_CCs$GO_Cellular_Component_2018
Neg_GO_CCs <- Neg_GO_CCs[Neg_GO_CCs$P.value <= 0.05, ]

write.table(Neg_GO_CCs, file = "./Objs/DownRegulated_GO_CCs.csv", quote = T, sep = "\t", col.names = T, row.names = T)

################3
Pos_KEGG <- Enriched_PositiveGns["KEGG_2019_Human"]
Pos_KEGG <- Pos_KEGG$KEGG_2019_Human
Pos_KEGG <- Pos_KEGG[Pos_KEGG$P.value <= 0.05, ]

write.table(Pos_KEGG, file = "./Objs/UpRegulated_KEGG.csv", quote = T, sep = "\t", col.names = T, row.names = T)

################## 
Neg_KEGG <- Enriched_NegativeGns["KEGG_2019_Human"]
Neg_KEGG <- Neg_KEGG$KEGG_2019_Human
Neg_KEGG <- Neg_KEGG[Neg_KEGG$P.value <= 0.05, ]

write.table(Neg_KEGG, file = "./Objs/DownRegulated_KEGG.csv", quote = T, sep = "\t", col.names = T, row.names = T)

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