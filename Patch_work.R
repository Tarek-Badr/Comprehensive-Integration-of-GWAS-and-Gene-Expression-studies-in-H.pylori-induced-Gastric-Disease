rm(list = ls())

setwd("C:/Users/ngs-adm/Desktop")

devtools::install_github("thomasp85/patchwork")
install.packages("patchwork")


par(mfrow=c(3,1))
P1 = BarPlot_Top10_Pathways_GO
P2 = BarPlot_Top10_Pathways_KEGG
P3 = BarPlot_Top10_Pathways_all


P1 + P2 + P3 + plot_layout(ncol = 1)

p1 + p2

(p1 | p2 | p3) /p4
