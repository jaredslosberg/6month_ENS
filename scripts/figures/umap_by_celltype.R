#Visualize UMAP colored by cell type
#visualize UMAP colored by batch

library(dplyr)
library(monocle3)
library(here)
library(ggplot2)
library(scales)

source(here("scripts/accessory_functions/monocle_mods.R"))

cds <- readRDS(here("6month_LMMP.rds"))

#by cell type
pdf(here("plots/supp_figures/UMAP_by_celltype.pdf"))
plot_cells(cds, color_cells_by = "cell_type_factor", label_cell_groups = F) + theme(legend.position = "none")

p <- plot_cells(cds, color_cells_by = "cell_type_factor", label_cell_groups = F) 
p$guides$colour$title <- "Cell Type"
p

plot_cells(cds, color_cells_by = "cell_type_factor", group_label_size = 4)   
dev.off()

#by batch
pData(cds) <- pData(cds) %>% as.data.frame() %>%
  mutate(batch_recode = factor(case_when(
    batch == "TH" ~ 1,
    batch == "TL" ~ 2
  ))) %>% DataFrame


q2 <- plot_cells(cds, color_cells_by = "batch_recode", label_cell_groups = F, cell_size = 0.85) 
q2$guides$colour$title <- "Batch"
q2

q3 <- plot_cells_mod(cds, color_cells_by = "batch_recode", label_cell_groups = F, outline_size = 0.95) + facet_wrap(~batch_recode)  
q3$guides$colour$title <- "Batch"
q3

pdf(here("plots/supp_figures/UMAP_by_batch.pdf"), width = 8)
q2
dev.off()

pdf(here("plots/supp_figures/UMAP_by_batch_facet.pdf"), width = 12)
q3
dev.off()


