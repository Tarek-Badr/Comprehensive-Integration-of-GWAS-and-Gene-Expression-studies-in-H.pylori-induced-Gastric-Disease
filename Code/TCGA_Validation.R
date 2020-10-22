

# load TCGA Pheno
TCGA_pheno <- read.delim("./Data/TCGA/TCGA_ClinicalData.tsv")

# load TCGA Expr
TCGA_expr <- read.delim("./Data/TCGA/data_RNA_Seq_v2_expression_median.txt")


################3
TCGA_expr <- TCGA_expr[!duplicated(TCGA_expr$Hugo_Symbol), ]
TCGA_expr <- TCGA_expr[!is.na(TCGA_expr$Hugo_Symbol), ]
rownames(TCGA_expr) <- TCGA_expr$Hugo_Symbol
TCGA_expr$Hugo_Symbol <- NULL
TCGA_expr$Entrez_Gene_Id <- NULL
TCGA_expr <- TCGA_expr[-1, ]
TCGA_expr <- TCGA_expr[!is.na(rownames(TCGA_expr)), ]
TCGA_expr <- TCGA_expr[!(rownames(TCGA_expr) == ""), ]

sel <- which(apply(TCGA_expr, 1, function(x) all(is.finite(x)) ))
TCGA_expr <- TCGA_expr[sel, ]


colnames(TCGA_expr) <- gsub("\\.", "-", colnames(TCGA_expr))

# Convert to numeric matrix
COLS <- colnames(TCGA_expr)
ROWS <- rownames(TCGA_expr)
TCGA_expr <- matrix(as.numeric(unlist(TCGA_expr)),nrow=nrow(TCGA_expr))
rownames(TCGA_expr) <- ROWS
colnames(TCGA_expr) <- COLS

TCGA_expr <- log2(TCGA_expr + 1)
boxplot(TCGA_expr[,1:15])


################## 
## Modify TCGA pheno
rownames(TCGA_pheno) <- TCGA_pheno$Sample.ID
TCGA_pheno$DiseaseStatus <- TCGA_pheno$Overall.Survival.Status
TCGA_pheno$DiseaseStatus[TCGA_pheno$DiseaseStatus == "0:LIVING"] <- "control"
TCGA_pheno$DiseaseStatus[TCGA_pheno$DiseaseStatus == "1:DECEASED"] <- "case"

TCGA_pheno$DiseaseStatus <- factor(TCGA_pheno$DiseaseStatus, levels = c("control", "case"))
table(TCGA_pheno$DiseaseStatus)

## Modify sample names to match sample names in the expression
TCGA_pheno <- TCGA_pheno[rownames(TCGA_pheno) %in% colnames(TCGA_expr), ]

TCGA_pheno <- TCGA_pheno[order(rownames(TCGA_pheno)), ]
TCGA_expr <- TCGA_expr[, order(colnames(TCGA_expr)), ]

all(rownames(TCGA_pheno) == colnames(TCGA_expr))

## Finally, replace the expression and phenotype data in the dataset with the new modified versions
TCGA_Dataset <- list()
TCGA_Dataset$expr <- TCGA_expr
TCGA_Dataset$pheno <- TCGA_pheno
TCGA_Dataset$keys <- rownames(TCGA_expr)
TCGA_Dataset$formattedName <- "TCGA"

TCGA_Dataset <- classFunction(TCGA_Dataset, column = "DiseaseStatus", diseaseTerms = c("case"))


rocPlot(datasetObject = TCGA_Dataset, filterObject = New_filter2)






