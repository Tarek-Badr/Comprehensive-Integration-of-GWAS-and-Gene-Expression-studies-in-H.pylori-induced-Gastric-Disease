
# Pathway Analysis of K313 Mutants

setwd("C:/Users/ngs-adm/Desktop/00-ian_Reanalysis_2020")

DEG_Entriz

FDR001_UP = read.delim("000-FDR001_Up_EntrezID.txt",sep='\t',header = T)
FDR001_dn = read.delim("000-FDR001_Dn_EntrezID.txt",sep='\t',header = T)

DEG_Entriz = read.delim("DEG_Entriz.txt",sep='\t',header = T)

cls <- c(Upregulated="character", Downregulated="character")

read.csv("data.csv", colClasses=cls, stringsAsFactors=FALSE)
d <- read.csv("data.csv",stringsAsFactors=F, na.strings="unknown")

DEG_Ent = read.delim("DEG_Entriz.txt",colClasses = cls, stringsAsFactors=FALSE, sep='\t',header = T)

DEG_Ent_t = t(DEG_Ent)


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("KEGG.db")
BiocManager::install("gage")
BiocManager::install("clusterProfiler")
BiocManager::install("fgsea")
library(KEGG.db)
library(Rgraphviz)
library(png)
library(KEGGgraph)
library(dplyr)
library(tibble)
library(gageData)
library(AnnotationDbi)
library(org.Hs.eg.db)
install.packages("tidyverse")
library(tidyverse)
library(igraph)
library(enrichplot)
library(DOSE)
library(ggplot2)


library(clusterProfiler)
library(gage)
library(pathview)
library(data.table)
library(fgsea)
library(ggplot2)

BiocManager::install("DOSE")

library(org.Mm.eg.db)
keytypes(org.Mm.eg.db)

BiocManager::install("org.Mm.eg.db")
BiocManager::install("org.Hs.eg.db")

data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]

FDR001_UP
gene_Up = FDR001_UP$Entrez_ID
gene_Dn = FDR001_dn$Entrez_ID

ggo <- groupGO(gene     = gene,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)
kk <- enrichKEGG(gene = gene,
                 keyType = "kegg",
                 organism     = 'hsa',
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", universe,
                 maxGSSize = 500, qvalueCutoff = 0.2,
                 use_internal_data = FALSE)


kk_dn <- enrichKEGG(gene = gene_Dn,
                    keyType = "kegg",
                    organism     = 'mmu',
                    pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                    maxGSSize = 500, qvalueCutoff = 0.2,
                    use_internal_data = FALSE)

kk_up <- enrichKEGG(gene = gene_Up,
                    keyType = "kegg",
                    organism     = 'mmu',
                    pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                    maxGSSize = 500, qvalueCutoff = 0.2,
                    use_internal_data = FALSE)

mkk_dn <- enrichMKEGG(gene = gene_Dn,
                   organism = 'mmu')

dp_dn = dotplot(kk_dn)
bp_dn = barplot(kk_dn, showCategory=8)

gseaplot(kk_dn, geneSetID = "mmu04640")


dp_up = dotplot(kk_up)
bp_up = barplot(kk_up, showCategory=8)

dp_all = dp_dn + dp_up + plot_layout(ncol = 1)
bp_all = bp_dn + bp_up + plot_layout(ncol = 1)


####compare cluster

ck <- compareCluster(geneCluster = gcSample, fun = "enrichKEGG")
head(as.data.frame(ck))


xx <- compareCluster(DEG_Ent, fun="enrichKEGG", organism="mmu")
head(as.data.frame(xx))

xx_2 <- compareCluster(DEG_Ent, fun="enrichKEGG", organism="mmu", pvalueCutoff=0.05)
head(as.data.frame(xx_2))

xx_ggo <- compareCluster(DEG_Ent,  fun='groupGO', OrgDb='org.Mm.eg.db')
head(as.data.frame(xx_ggo))

xx_ep<- compareCluster(DEG_Ent,  fun='enrichGO', OrgDb='org.Mm.eg.db')
head(as.data.frame(xx_ep))

xx_DO <- compareCluster(DEG_Ent,  fun='enrichDO')
head(as.data.frame(xx_DO))

dotplot(xx_ep) #GO Enrichment Analysis
dotplot(xx) #KEGG Enrichment Analysis

write.csv(as.data.frame(xx), file="K313_fdr001_FC2_KEGG_pathway.csv")


cnetplot(xx)



