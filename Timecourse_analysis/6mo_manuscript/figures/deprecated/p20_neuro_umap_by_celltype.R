library(monocle3)
library(tidyverse)
library(pals)
library(cowplot)
library(ggplot2)

global_dpi = 150


##Figure 1a ----
cds <- readRDS(here("tc_neuro_p20.rds"))
colors <- c("#2EA096","#A01515","#C3A459") %>% setNames(c("MENs", "NC-glia", "NC-neurons"))

pData(cds)$cell_type_aggregate <- pData(cds)$cell_type_aggregate %>%
  recode("NENs" = "NC-neurons",
         "Neuroglia" = "NC-glia") 

# pData(cds)$figure_celltype <- pData(cds)$figure_celltype %>%
#   recode("NENs" = "NC-neurons", "neuroglia/neuroendocrine" = "NC-neurons", "Neuroglia" = "NC-glia") %>% 
#   factor(levels = c("MENs", "NC-glia", "NC-neurons", "Fibroblasts", "ICC", "SMC", "Macrophage", "Endothelial", "Enterocytes", "Mucosal epith.", "RBC")) 

pData(cds)[,"Cell Type"] <- pData(cds)$cell_type_aggregate
color_by <- "Cell Type"

#by cell type
umap <- plot_cells(cds, color_cells_by = color_by, label_cell_groups = F) +
  theme(axis.line.x=element_blank(),
        axis.line.y=element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title = element_blank(), 
        legend.key.height = unit(.5, "cm"),
        legend.text = element_text(size = 8)) + 
  scale_color_manual(name = "xxx", values = colors) 

#save legend as grob
lgd <- (umap + 
          theme(legend.title.align = 0.5, 
                legend.background = element_rect(fill = "transparent"),
                legend.key = element_rect(fill = "transparent")) + 
          guides(color = guide_legend(ncol = 1, title = "Cell Type",
                                      override.aes = list(size=4)))) %>%
  get_legend()

umap_axes <- ggplot() + theme_minimal() +
  theme(axis.line = element_line(size = .5, arrow = arrow(type='closed', length = unit(6,'pt'))),
        plot.margin = margin(10, 10, 10, 10, "pt")) + 
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.title.x = element_text(vjust = -1, hjust = 0),
        axis.title.y = element_text(vjust = 1, hjust = 0),
        axis.title = element_text(size = 8),
        panel.background = element_rect(fill = "transparent", colour = "transparent"))

proto <- umap + theme(legend.position = "none",  plot.margin = margin(5.5,5.5,5.5,50)) + 
  annotation_custom(lgd, xmin = -20, xmax = -5, ymin = 5, ymax = 11) + 
  annotation_custom(as_grob(umap_axes), xmin = -20, xmax = -10, ymin = -10, ymax = -2.5)

fig <- ggrastr::rasterise(proto, dpi = global_dpi)

pdf(here("plots/6mo_manuscript/p20_neuro_UMAP_by_celltype.pdf"), height = 3, width = 3)
    fig
dev.off()
    

  
  
  
  

