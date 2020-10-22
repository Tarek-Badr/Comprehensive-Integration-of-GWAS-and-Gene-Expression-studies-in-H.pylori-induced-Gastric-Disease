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

setwd("C:/Users/ngs-adm/Desktop/")

data=NULL;eset=NULL;data.exprs=NULL;

gset <- getGEO("GSE70394", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) det <- grep("GPL6255", attr(gset, "names")) else det <- 1

gsedata <- gset[[det]]

fvarLabels(gsedata) <- make.names(fvarLabels(gsedata))

gsmid <- "010101"
sml <- c()
for (i in 1:nchar(gsmid)) { sml[i] <- substr(gsmid,i,i) }

exp_data <- exprs(gsedata)
qx <- as.numeric(quantile(exp_data, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { exp_data[which(exp_data <= 0)] <- NaN
exprs(gsedata) <- log2(exp_data) }
sml <- paste("G", sml, sep="")   

samplev <- as.factor(sml)

gsedata$description <- samplev

design <- model.matrix(~ description + 0, gsedata)

colnames(design) <- levels(samplev)

datafit <- lmFit(gsedata, design)

cont <- makeContrasts(G1-G0, levels=design)

datafit2 <- contrasts.fit(datafit, cont)

datafit2 <- eBayes(datafit2, 0.01)

DEGs <- topTable(datafit2, adjust="fdr")

DEGs <- subset(DEGs, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))

write.table(DEGs, file=stdout(), row.names=F, sep="\t")
write.csv(as.data.frame(DEGs), file="GSE70394_d´DEGs.csv")
