library(dplyr)
library(monocle3)
library(here)
library(ggplot2)

cds <- readRDS(here("6month_LMMP.rds"))

pdf(here("plots/supp_figures/UMAP_by_celltype.pdf"))
plot_cells(cds, label_cell_groups = F) + theme(legend.position = NULL)
plot_cells(cds, color_cells_by = "cell_type", group_label_size = 4)   
dev.off()
