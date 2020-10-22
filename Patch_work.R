devtools::install_github("thomasp85/patchwork")


install.packages("patchwork")

par(mfrow=c(3,1))
P1 = BarPlot_Top20Family_GTDB
P2 = BarPlot_Top20Genus_GTDB 
P3 = BarPlot_Top20SPP_GTDB


P1 + P2 + P3 + plot_layout(ncol = 1)

p1 + p2

(p1 | p2 | p3) /p4