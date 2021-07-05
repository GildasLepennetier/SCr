library(limma) # For a more efficient implementation of the Wilcoxon Rank Sum Test
library(Seurat)
library(ggplot2) #ggplot
library(ggrepel) #geom_label_repel() not working with Dimplot?
library(dplyr) #%>%
library(tidyr) #separate
library(ggpubr) #grids
library(openxlsx) # read.xlsx
library(pheatmap) #pheatmap
library(patchwork) # wrap_plots
library(gtools) #mixedsort
library(cowplot) # cowplot::plot_grid(gg1, gg2, gg3, nrow = nrow, align = "h", rel_widths = c(1, 1, 1))
library(gt) # great table gt::gt
library(Nebulosa) #plot_density (bioconductor) #https://bioconductor.org/packages/devel/bioc/vignettes/Nebulosa/inst/doc/nebulosa_seurat.html
library(igraph) # used in cellchat
library(CIPR)
library(CellChat) # in case we saved objects: updateCellChat(Object)
library(ggalluvial)
library(harmony)
library(openxlsx)
library(gplots)
library(SCENIC)
	
source("https://raw.githubusercontent.com/GildasLepennetier/R-stuff/master/not_in.R")
source("https://raw.githubusercontent.com/GildasLepennetier/SCr/main/FAM.R")
source("https://raw.githubusercontent.com/GildasLepennetier/SCr/main/ANN.R")
source("https://raw.githubusercontent.com/GildasLepennetier/SCr/main/DOT.R")
source("https://raw.githubusercontent.com/GildasLepennetier/SCr/main/MEXL.R")
source("https://raw.githubusercontent.com/GildasLepennetier/SCr/main/CChat.R")
source("https://raw.githubusercontent.com/GildasLepennetier/SCr/main/SCEN.R")
