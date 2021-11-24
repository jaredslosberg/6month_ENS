#Visualize the nonfetal & healthy subset of gut cell atlas. 
#Facets by age and their annotated cell types
#Expression of stx3, uchl1, slc17a9

library(tidyverse)
library(monocle3)
library(here)

source(here("scripts/accessory_functions/monocle_mods.R"))

nfh <- readRDS(here("teichman_gut_atlas/data/mesenchymal_gut_nfh_annotated.rds"))

# png(here("plots/supp_figures/mesenchymal_nfh_age_facet.png"), width = 1000, height = 750)

pdf(here("plots/supp_figures/mesenchymal_nfh_age_facet.pdf"))
#use bigger cell size for png
plot_cells_mod(nfh, label_cell_groups = F, color_cells_by = "Age", cell_size = .4) +
  facet_wrap(~Age) + 
  theme(legend.position = "none")
dev.off()

pdf(here("plots/supp_figures/mesenchymal_nfh_celltype_facet.pdf"))
#use bigger cell size for png
plot_cells_mod(nfh, label_cell_groups = F, color_cells_by = "Integrated_05", cell_size = .4) +
  facet_wrap(~Integrated_05) + 
  theme(legend.position = "none")
dev.off()

pdf(here("plots/supp_figures/mesenchymal_nfh_uchl1_stx3_slc17a9.pdf"), width = 9)
#use bigger cell size for png
plot_cells(nfh,
           genes = c("SLC17A9","STX3","UCHL1"),
           cell_size = .35) 
dev.off()
