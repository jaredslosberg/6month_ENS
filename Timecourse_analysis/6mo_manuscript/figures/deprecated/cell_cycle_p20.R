library(monocle3)
library(tidyverse)
library(here)
library(tricycle)

source(here("scripts/accessory_functions/monocle_mods.R"))
source(here("scripts/accessory_functions/theme_figure.R"))

cds <- readRDS(here("tc_neuro_p20.rds"))
pData(cds)$cell_type <- pData(cds)$cell_type %>%
  recode("NENs" = "NC-neurons", "neuroglia/neuroendocrine" = "NC-neurons", "Neuroglia" = "NC-glia") %>% 
  as.factor() %>% 
  relevel("MENs", "NC-glia", "NC-neurons")

colors <- c("#2EA096","#A01515","#C3A459") %>% setNames(c("MENs", "NC-glia", "NC-neurons"))
cc_define <- list("lower" = c(0, pi/2), "upper" = c(3*pi/2,2*pi))

#Tricycle density plot----
grDevices::cairo_pdf(here("plots/6mo_manuscript/tricycle_density_p20.pdf"), height = 1.5, width = 4)
#Continuous
tricycle::plot_ccposition_den(cds$tricyclePosition,cds$cell_type,
                              'cell type', type = "linear",bw = 10, palette.v = colors,
                              fig.title = "", line.size = 1) +
# theme_bw(base_size = 14) +
  theme(plot.title = element_blank()) +
  annotate("rect", xmin=cc_define$lower[1], xmax=cc_define$lower[2], ymin=-.025, ymax=.6, alpha=0.1, fill="black") +
  annotate("rect", xmin=cc_define$upper[1], xmax=cc_define$upper[2], ymin=-.025, ymax=.6, alpha=0.1, fill="black") +
    scale_x_continuous(limits = c(0, 2 * pi), breaks = c(0, pi/2, pi, 3 * pi/2, 2*pi),
                     labels = paste0(c(0,0.5, 1, 1.5, 2), "π"), name = "", expand = c(0,0)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0))) +
  theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal")
 #       legend.text = element_text(size = 8), legend.title = element_text(size = 10) ) 
dev.off()


#Embedding ----
pData(cds)$dummy = "All"
hue.colors = c("#2E22EA", "#9E3DFB", "#F86BE2", "#FCCE7B", 
               "#C4E416", "#4BBA0F", "#447D87", "#2C24E9")
hue.n = 500

pdf(here("plots/6mo_manuscript/tricycle_embedding_p20.pdf"), width = 4, height =4)
plot_cells_mod(cds, reduction_method = "tricycleEmbedding",
               color_cells_by = "tricyclePosition", cell_size = 1, label_cell_groups = F, ) + 
  theme_figure() +
  geom_vline(xintercept = 0, linetype= "dashed", size = .3) +
  geom_hline(yintercept = 0, linetype = "dashed", size = .3) + theme(legend.position = "none") +
  ggtitle("") + 
  facet_wrap(~dummy) +
  theme(strip.text.x = element_text(size = 8),
        panel.grid = element_blank()) +
  scale_color_gradientn(name = "tricyclePosition", limits = range(0,2 * pi),
                        breaks = seq(from = 0, to = 2 * pi, length.out = hue.n),
                        colors = hue.colors, guide = "none")
dev.off()

#Embedding faceted by celltype ----
pdf(here("plots/6mo_manuscript/tricycle_embedding_p20_facet.pdf"), width =8 , height = 4)
plot_cells_mod(cds, reduction_method = "tricycleEmbedding",
               color_cells_by = "cell_type", cell_size = 1, label_cell_groups = F, ) +
  theme_figure() +
  facet_wrap(~cell_type) +
  scale_color_manual(values = colors) +
  geom_vline(xintercept = 0, linetype= "dashed", size =.3) +
  geom_hline(yintercept = 0, linetype = "dashed", size = .3) + theme(legend.position = "none") + 
  theme(strip.text.x = element_text(size = 8),
        panel.grid = element_blank()) +
  ggtitle("")
dev.off()

#proportion of cell cycling----
cycling_table <- read.csv(here("results/supp_figures/cycling_proportions_neuro_p20.csv"), row.names = 1) %>%
  rownames_to_column(var = "cell_type")

cycling_table$cell_type <-cycling_table$cell_type %>% 
  recode("NENs" = "NC-neurons", "Neuroglia" = "NC-glia") %>% 
  as.factor() %>% 
  relevel("MENs", "NC-glia", "NC-neurons")

pdf(here("plots/6mo_manuscript/cycling_proportions_neuro_p20.pdf"), width = 4, height = 3)

  ggplot(cycling_table) +
    geom_col(aes(x = cell_type , y = prop_cycling, fill = cell_type)) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(breaks = c(0, .1, .2, .3, .4), labels = paste0(c(0, 10, 20, 30, 40), "%")) +
    theme_bw() +
    theme(panel.grid.major.x  = element_blank(),
          axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 11),
          axis.title.y = element_text(vjust = 2),
          axis.title.x = element_text(vjust = -.5)) + 
    ylab("Proportion of cells cycling") + 
    xlab("Cell Type") +
    theme(legend.position= "none")

dev.off()

