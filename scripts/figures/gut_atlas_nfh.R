#Visualize the nonfetal & healthy subset of gut cell atlas. 
#Facets by age and their annotated cell types
#Expression of stx3, uchl1, slc17a9

library(tidyverse)
library(monocle3)
library(here)

source(here("scripts/accessory_functions/monocle_mods.R"))

nfh <- readRDS(here("teichman_gut_atlas/data/mesenchymal_gut_nfh_annotated.rds"))
mesenchymal <- readRDS(here("teichman_gut_atlas/data/mesenchymal_gut_annotated.Rds"))


# png(here("plots/supp_figures/mesenchymal_nfh_age_facet.png"), width = 1000, height = 750)

pData(nfh)$pseudo <- as.factor("group")

pdf(here("plots/supp_figures/mesenchymal_nfh_age_facet.pdf"))
#use bigger cell size for png
plot_cells_mod(nfh, label_cell_groups = F, color_cells_by = "Age",  cell_size = .4) +
  facet_wrap(~Age) + 
  scale_color_manual(values = rep("darkblue",11)) +
  theme(legend.position = "none") 
dev.off()

pdf(here("plots/supp_figures/mesenchymal_nfh_celltype_facet.pdf"))
#use bigger cell size for png
plot_cells_mod(nfh, label_cell_groups = F, color_cells_by = "Integrated_05", cell_size = .4) +
  facet_wrap(~Integrated_05) + 
  scale_color_manual(values = rep("darkred",15)) +
  theme(legend.position = "none")
dev.off()

pdf(here("plots/supp_figures/mesenchymal_nfh_uchl1_stx3_slc17a9.pdf"), width = 9)
#use bigger cell size for png
plot_cells(nfh,
           genes = c("SLC17A9","STX3","UCHL1"),
           cell_size = .35,
           label_cell_groups = F) 
dev.off()

pdf(here("plots/supp_figures/mesenchymal_nfh_uchl1_stx3_slc17a9_hand2_pde10a_tubb2b.pdf"), width = 9, height = 6)
#use bigger cell size for png
plot_cells(nfh,
           genes = c("SLC17A9","STX3","UCHL1","TUBB2B" ,"PDE10A","HAND2"),
           cell_size = .35,
           label_cell_groups = F) 
dev.off()

#mesenchymal, non fetal healthy
pdf(here("plots/supp_figures/mesenchymal_nfh_umap.pdf"))
#use bigger cell size for png
plot_cells(nfh, group_label_size = 5)
plot_cells(nfh, label_cell_groups = F)

dev.off()

#whole mesenchymal 
pdf(here("plots/supp_figures/mesenchymal_umap.pdf"))
#use bigger cell size for png
plot_cells(mesenchymal, group_label_size = 5)
plot_cells(mesenchymal, label_cell_groups = F)

dev.off()

pdat <- pData(nfh)  %>% as.data.frame()

pdf(here("plots/mesenchymal_nfh_cluster_by_age"), width = 9)
  ggplot(pdat) + geom_bar(aes(x = clusters, fill= Age), position = "fill")
dev.off()

pdf(here("plots/mesenchymal_nfh_age_by_cluster"), width = 9)
ggplot(pdat) + geom_bar(aes(x = Age, fill= clusters), position = "fill") + ylab("Proportion") + theme_minimal() + 
  theme(panel.grid.major.x = element_blank(),axis.ticks.y = element_line(color = "black"), axis.ticks.length.y = unit(-.15, "cm"))
dev.off()
