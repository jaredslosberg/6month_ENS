library(SingleCellExperiment)
library(monocle3)
library(here)
library(dplyr)
library(ggplot2)

cds <- readRDS(here("6month_LMMP.rds"))

cell_types_to_keep <- c("MENs","NENs","Neuroglia")

#if needed, rename cell types for figure
cds_sub <- cds[,pData(cds)$cell_type %in% cell_types_to_keep]
pData(cds_sub)$cell_type <- pData(cds_sub)$cell_type %>%
  recode("NENs" = "NC-neurons", "Neuroglia" = "NC-glia")

pData(cds_sub) <- pData(cds_sub) %>% as.data.frame() %>%
  group_by(cell_type) %>%
  mutate(group_title = paste0(cell_type, " (n = ", n() ,")")) %>%
  ungroup() %>%
  DataFrame()
  
#
genes <- c("Calcb","Aebp1","Cdh3","Cftr","Clic3","Fmo2","Ntf3","Slpi","Smo","Myl7","Met","Il18","Slc17a9","Ret","Nos1","Sox10")
  

vl <- plot_genes_violin(cds_sub[fData(cds)$gene_short_name %in% genes,], group_cells_by = "group_title") + 
  xlab("") + 
  theme(axis.text.x = element_text(size = 14)) + 
  scale_fill_manual(values = rev(wesanderson::wes_palette("GrandBudapest1", n =3)))

pdf(here("plots/supp_figures/violin_expression.pdf"), width = 8, height = 16) 
  vl
dev.off()

vl
