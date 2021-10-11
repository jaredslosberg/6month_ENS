library(dplyr)
library(monocle3)
library(here)
library(ggplot2)
library(scales)


cds <- readRDS(here("6month_LMMP.rds"))

pdf(here("plots/supp_figures/UMAP_by_celltype.pdf"))
plot_cells(cds, color_cells_by = "cell_type_factor", label_cell_groups = F) + theme(legend.position = "none")

p <- plot_cells(cds, color_cells_by = "cell_type_factor", label_cell_groups = F) 
p$guides$colour$title <- "Cell Type"
p

plot_cells(cds, color_cells_by = "cell_type_factor", group_label_size = 4)   
dev.off()


