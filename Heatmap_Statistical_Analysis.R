getwd()
setwd("C:/Users/Lenovo R2G/Desktop")
getwd()
C:/Users/Lenovo R2G/Desktop/wxs_CRC_genes_ordered

#####################################
#Making Heatmaps for Mutations per CLC gene
#####################################

wxs_crc = read.delim("wxs_CRC_genes_ordered.txt",sep='\t',row.names = 1)
wxs_crc


Mutec_Variants = read.table("Mutec_Variants.txt",sep='\t', header = T)
Mutec_Variants


#normalization (as samples size < 30, rld is prefered than vst)

library("pheatmap")
library("RColorBrewer")
library(scales)
library(ggplot2)


################Heatmap of mutations per CRC highly mutated genes
pheatmap(wxs_crc)

pheatmap(wxs_crc, fontsize_number = 2 * fontsize,
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

wxs_crc_norm <- t(apply(wxs_crc, 1, cal_z_score))
pheatmap(wxs_crc_norm)


pheatmap(wxs_crc_norm, fontsize_number = 2 * fontsize,
         cellwidth = 15, cellheight = 12, scale = "none",
         treeheight_row = 100,
         kmeans_k = NA,
         show_rownames = T, show_colnames = T,
         main = "Full heatmap (avg, eucl, unsc)",
         clustering_method = "average",
         cluster_rows = FALSE, cluster_cols = FALSE,
         clustering_distance_rows = drows1, 
         clustering_distance_cols = dcols1, angle_col = 45) 


################Heatmap of mutations according to mutation type

Mutations_HM = read.delim("Mutec_Variants.txt",sep='\t',row.names = 1)

Mutations_HM

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

Mutations_HM_norm <- t(apply(Mutations_HM, 1, cal_z_score))
pheatmap(Mutations_HM_norm)


pheatmap(Mutations_HM_norm, fontsize_number = 2 * fontsize,
         cellwidth = 15, cellheight = 12, scale = "none",
         treeheight_row = 100,
         kmeans_k = NA,
         show_rownames = T, show_colnames = T,
         main = "Full heatmap (avg, eucl, unsc)",
         clustering_method = "average",
         cluster_rows = TRUE, cluster_cols = TRUE) 

####################################
#Adjust Table for plotting and Statistical significance
####################################

Mutations = read.delim("Mutations.txt",sep='\t', header = T)
Mutations

head(Mutations)

df.long <- pivot_longer(Mutations, cols=-1, names_to = "Variants", values_to = "Count")
df.long

df.long$Count <- as.numeric(as.vector(df.long$Count))
ggplot(data=df.long, aes(x=Variants, y=Count, fill=Mouse)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()

df.long$Count <- as.numeric(as.vector(df.long$Count))

##################plotting

P = ggplot(data = df.long, aes(Variants, Count)) +
  geom_boxplot(aes(colour = Mouse))

sp = P + theme(axis.text.x = element_text(face = "bold", color = "#993333", 
                                          size = 12, angle = 45, hjust = 1),
               axis.text.y = element_text(face = "bold", color = "blue", 
                                          size = 12, angle = 45))

###best one##
SD = sp + scale_y_continuous(trans = log2_trans(), breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x)))
SD
#####box Plots
P_box = ggboxplot(df.long, x = "Variants", y = "Count", color = "Mouse", palette = "jco",
                  add = "jitter")

P_box
sp_box = P_box + theme(axis.text.x = element_text(face = "bold", color = "#993333", 
                                          size = 12, angle = 45, hjust = 1),
               axis.text.y = element_text(face = "bold", color = "blue", 
                                          size = 12, angle = 45))
sp_box

SD_box = sp_box + scale_y_continuous(trans = log2_trans(), breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x)))
SD_box





#############################################statistical testing
sp  + annotation_logticks() 


SD_box + stat_compare_means(aes(group = Mouse))
SD_box + stat_compare_means(aes(group = Mouse), label =  "p.signif", label.x = 1.5) #wilcoxon test

##############################################################################################################
#Plotting the vcf files in waterfall plots and mutation load plots
#############################################################################################################














