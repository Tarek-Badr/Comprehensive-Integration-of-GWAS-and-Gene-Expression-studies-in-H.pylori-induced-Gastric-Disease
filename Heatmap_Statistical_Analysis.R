getwd()
setwd("C:/Users/Lenovo R2G/Desktop")
getwd()
C:/Users/Lenovo R2G/Desktop/wxs_CRC_genes_ordered

DEGs = read.delim("DEGs_genes_ordered.txt",sep='\t',row.names = 1)
DEGs

library("pheatmap")
library("RColorBrewer")
library(scales)
library(ggplot2)

pheatmap(DEGs)

pheatmap(DEGs, fontsize_number = 2 * fontsize,
         cellwidth = 15, cellheight = 12, scale = "none",
         treeheight_row = 100,
         kmeans_k = NA,
         show_rownames = T, show_colnames = T,
         main = "Full heatmap (avg, eucl, unsc)",
         clustering_method = "average",
         cluster_rows = FALSE, cluster_cols = FALSE,
         clustering_distance_rows = drows1, 
         clustering_distance_cols = dcols1, angle_col = 45) 


cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

DEGs_norm <- t(apply(DEGs, 1, cal_z_score))
pheatmap(DEGs_norm)


pheatmap(DEGs_norm, fontsize_number = 2 * fontsize,
         cellwidth = 15, cellheight = 12, scale = "none",
         treeheight_row = 100,
         kmeans_k = NA,
         show_rownames = T, show_colnames = T,
         main = "Full heatmap (avg, eucl, unsc)",
         clustering_method = "average",
         cluster_rows = FALSE, cluster_cols = FALSE,
         clustering_distance_rows = drows1, 
         clustering_distance_cols = dcols1, angle_col = 45) 
