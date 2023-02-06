#Plot expression of Ret and Met in p20 TC neuro data


library(monocle3)
library(tidyverse)
library(here)

#at some point this became the p20 subset
cds <- readRDS(here("TC_Neuro.rds"))


colors <- c("#2EA096","#A01515","#C3A459") %>% setNames(c("MENs", "NC-glia", "NC-neurons"))

pData(cds)$'Cell Type' <- pData(cds)$cell_type_aggregate %>%
  recode("NENs" = "NC-neurons", "neuroglia/neuroendocrine" = "NC-neurons", "Neuroglia" = "NC-glia")

pdf(here("plots/6mo_manuscript/p20_neuro_celltype.pdf"), width = 6, height = 6)
plot_cells(cds, color_cells_by = "Cell Type", label_cell_groups = F, cell_size = 0.45) +
  scale_color_manual(values = colors) +
  theme(legend.position = c(0.2,0.8))
dev.off()

pdf(here("plots/6mo_manuscript/p20_neuro_ret_met.pdf"), width = 12, height = 6)
plot_cells(cds, genes = c("Ret", "Met"), label_cell_groups = F) + theme(strip.text = element_text(size = 14))
dev.off()