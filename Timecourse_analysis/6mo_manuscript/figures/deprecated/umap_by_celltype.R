#Visualize UMAP colored by cell type
#visualize UMAP colored by batch

library(dplyr)
library(monocle3)
library(here)
library(ggplot2)
library(scales)
library(pals)

source(here("scripts/accessory_functions/monocle_mods.R"))
source(here("scripts/accessory_functions/theme_figure.R"))


cds <- readRDS(here("tc_lmmp_p20.rds"))

pData(cds)$primary_identity <- pData(cds)$primary_identity %>%
  recode("NENs" = "NC-neurons", "neuroglia/neuroendocrine" = "NC-neurons", "Neuroglia" = "NC-glia") %>% 
  as.factor() %>% 
  relevel("MENs", "NC-glia", "NC-neurons")

color_by <- "primary_identity"
color_pal <- unname(pals::glasbey(length(unique(pData(cds)[,color_by]))))
#by cell type
pdf(here("plots/6mo_manuscript/p20_UMAP_by_celltype.pdf"), width = 6, height = 4)

plot_cells(cds, color_cells_by = color_by, label_cell_groups = F) +
  theme(axis.line.x=element_blank(),
      axis.line.y=element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(), 
      axis.text.y = element_blank(), 
      axis.title = element_blank(), 
      legend.key.height = unit(.5, "cm"),
      legend.text = element_text(size = 9)) + 
  scale_color_manual(name = "xxx", values = color_pal) 
dev.off()

#by batch
pData(cds) <- pData(cds) %>% as.data.frame() %>%
  mutate(batch_recode = factor(case_when(
    batch == "0" ~ A,
    batch == "1" ~ B
  ))) %>% DataFrame


q2 <- plot_cells(cds, color_cells_by = "batch_recode", label_cell_groups = F, cell_size = 0.85) 
q2$guides$colour$title <- "Batch"
q2

q3 <- plot_cells_mod(cds, color_cells_by = "batch_recode", label_cell_groups = F, outline_size = 0.95) + facet_wrap(~batch_recode)  
q3$guides$colour$title <- "Batch"
q3

pdf(here("plots/6mo_manusript/p20_UMAP_by_batch.pdf"), width = 8)
q2
dev.off()

pdf(here("plots/6mo_manuscript/p20UMAP_by_batch_facet.pdf"), width = 12)
q3
dev.off()


