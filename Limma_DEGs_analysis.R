
rm(list = ls())

library(KEGG.db)
library(Rgraphviz)
library(png)
library(KEGGgraph)
library(dplyr)
library(tibble)
library(gageData)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(sva)
library(bladderbatch)
library(RColorBrewer)
library(pamr)
library(Biobase)
library(GEOquery)
library(limma)
library(genefilter)

setwd("C:/Users/ngs-adm/Desktop/")

data=NULL;eset=NULL;data.exprs=NULL;

HPyloriData <- getGEOData(c("GSE27411", "GSE60427", "GSE60662", "GSE5081"))
CellLinesData <- getGEOData(c("GSE39919", "GSE70394", "GSE74577", "GSE74492"))

#####an Example for a Dataset analysis with Limma

gset <- getGEO("GSE70394", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) det <- grep("GPL6255", attr(gset, "names")) else det <- 1
gsedata <- gset[[det]]
exp_data <- exprs(gsedata)
###########################################################################

fvarLabels(gsedata) <- make.names(fvarLabels(gsedata))

gsmid <- "010101"
sml_GSE70394 <- c()
for (i in 1:nchar(gsmid)) { sml_GSE70394[i] <- substr(gsmid,i,i) }

q_exp_data <- as.numeric(quantile(exp_data, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (q_exp_data[5] > 100) ||
  (q_exp_data[6]-q_exp_data[1] > 50 && q_exp_data[2] > 0) ||
  (q_exp_data[2] > 0 && q_exp_data[2] < 1 && q_exp_data[4] > 1 && q_exp_data[4] < 2)
if (LogC) { exp_data[which(exp_data <= 0)] <- NaN
exprs(gsedata) <- log2(exp_data) }

sml_GSE70394 <- paste("G", sml_GSE70394, sep="")   

samplev <- as.factor(sml_GSE70394) 

gsedata$description <- samplev

design <- model.matrix(~ description + 0, gsedata)

colnames(design) <- levels(samplev)

datafit <- lmFit(gsedata, design=model.matrix(~ description + 0, gsedata))

cont <- makeContrasts(G1-G0, levels=design)

datafit_cont <- contrasts.fit(datafit, cont)

datafit_cont <- eBayes(datafit_cont, 0.01)

#DEGs_GSE70394 <- topTable(datafit_cont, adjust="fdr", coef=2)
DEGs_GSE70394 <- topTable(datafit_cont, adjust="fdr")

DEGs_GSE70394 <- subset(DEGs_GSE70394, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))

write.table(DEGs_GSE70394, file=stdout(), row.names=F, sep="\t")
write.csv(as.data.frame(DEGs_GSE70394), file="GSE70394_DEGs_F.csv")
