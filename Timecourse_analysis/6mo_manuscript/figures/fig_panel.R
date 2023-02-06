
library(dplyr)
library(monocle3)
library(here)
library(ggplot2)
library(scales)
library(pals)
library(grid)
library(gridExtra)
library(viridis)
library(tidyverse)
library(tricycle)
library(ggrastr)
library(scattermore)
library(ggpubr)
library(cowplot)


source(here("scripts/accessory_functions/monocle_mods.R"))
source(here("scripts/accessory_functions/theme_figure.R"))
source(here("scripts/accessory_functions/pattern_plotting.R"))

global_dpi = 150
var_keep <- c("fig1a", "fig1b", "fig1c", "fig1d", "fig1e", "fig1f" ,"fig1c_key",
              lsf.str(), "var_keep", "global_dpi", "colors")

##Figure 1a ----
cds <- readRDS(here("tc_lmmp_p20.rds"))
colors <- c("#2EA096","#A01515","#C3A459") %>% setNames(c("MENs", "NC-glia", "NC-neurons"))

pData(cds)$figure_celltype <- pData(cds)$figure_celltype %>%
  recode("NENs" = "NC-neurons", "neuroglia/neuroendocrine" = "NC-neurons", "Neuroglia" = "NC-glia") %>% 
  factor(levels = c("MENs", "NC-glia", "NC-neurons", "Fibroblasts", "ICC", "SMC", "Macrophage", "Endothelial", "Enterocytes", "Mucosal epith.", "RBC")) 

pData(cds)[,"Cell Type"] <- pData(cds)$figure_celltype
color_by <- "Cell Type"
color_pal <- unname(pals::glasbey(length(unique(pData(cds)[,color_by]))))

#add mens, nc colors
color_pal[6] <- color_pal[1]
color_pal[1:3] <- colors

#by cell type
umap <- plot_cells(cds, color_cells_by = color_by, label_cell_groups = F) +
  theme(axis.line.x=element_blank(),
        axis.line.y=element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title = element_blank(), 
        legend.key.height = unit(.5, "cm"),
        legend.text = element_text(size = 6)) + 
  scale_color_manual(name = "xxx", values = color_pal) 

#save legend as grob
lgd <- (umap + 
          theme(legend.title.align = 0.5) + 
          guides(color = guide_legend(ncol = 2, title = "Cell Type",
                                           override.aes = list(size=4)))) %>%
  get_legend()

umap_axes <- ggplot() + theme_minimal() +
  theme(axis.line = element_line(size = 1, arrow = arrow(type='closed', length = unit(10,'pt'))),
        plot.margin = margin(10, 10, 10, 10, "pt")) + 
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.title.x = element_text(vjust = -1),
        axis.title.y = element_text(vjust = 1),
        panel.background = element_rect(fill = "transparent", colour = "transparent"))

fig1a_proto <- umap + theme(legend.position = "none",  plot.margin = margin(5.5,10,5.5,5.5)) + 
  annotation_custom(lgd, xmin = 6, xmax = 14, ymin = 6, ymax = 13) + 
  annotation_custom(as_grob(umap_axes), xmin = -13.5, xmax = -8, ymin = -7, ymax = -.5)

fig1a <- ggrastr::rasterise(fig1a_proto, dpi = global_dpi)

rm(list=ls()[! ls() %in% var_keep])

##Figure 1b ----
pattern_cell_weight_filenames <- here("results/6month_projection_weights/tc_in_6mo_LMMP_old_patterns.csv")
pattern_cell_weights <- read.csv(pattern_cell_weight_filenames,row.names = 1) %>%
  tibble::rownames_to_column(var = "cell_id")
colnames(pattern_cell_weights) <- colnames(pattern_cell_weights) %>% str_replace_all(., "X6mo_", "")

pattern_weight_filename <- here("results/6month_projection_weights/NMF_gene_weights.csv")
pattern_gene_weights <- read.csv(pattern_weight_filename, row.names = 1)
n_top_genes_show <- 10

lmmp_filename <- here("tc_lmmp_p20.rds")
lmmp <- readRDS(lmmp_filename)

pData(lmmp) <- pData(lmmp) %>% as.data.frame() %>%
  tibble::rownames_to_column(var = "cell_id") %>%
  left_join(., pattern_cell_weights) %>%
  tibble::column_to_rownames(var = "cell_id") %>%
  DataFrame

npattern <- 50
pattern_colnames <- paste0("cellPattern", 1:npattern)

#declare patterns of interest
MENS_patterns <- c(16,27,32,41)
select_patterns <- c(MENS_patterns)

lmmp_pattern_wt_plots <- lapply(select_patterns,
                                plotCellPatterns, 
                                cds_obj = lmmp,
                                pattern_prefix = "cellPattern",
                                do.clip = c(0.01,0.99),
                                outline_size = 1.5,
                                cell_size = 0.45,
                                outline_color= "black") %>%
  lapply(., function(pl){
    pl <- pl + theme(plot.title = element_text(hjust = 0.5)) 
  })



#gene_tables <-  lapply(paste0("cellPattern",select_patterns), plotGeneWeightTable,
#                       pattern_gene_weights, n_top_genes_show, "Gene weight")



lmmp_pattern_wt_plots_na <- lapply(lmmp_pattern_wt_plots, function(pl){
  pl2 <- pl +
    theme(axis.line.x=element_blank(),
          axis.line.y=element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(), 
          axis.text.y = element_blank(), 
          axis.title = element_blank(), 
          legend.title = element_blank(), 
          legend.key.width = unit(.1, "in"),
          legend.text = element_text(size = 8),
          plot.title = element_text(size = 10), 
          plot.margin = margin(20,5.5,5.5,5.5, "pt"))
  
  return(ggrastr::rasterise(pl2, dpi = global_dpi))
})

fig1b_proto <- do.call(ggarrange, c(lmmp_pattern_wt_plots_na, nrow = 2, ncol = 2))

fig1b <- annotate_figure(fig1b_proto, top = text_grob("P20 LMMP Projected into Learned 6month Patterns", 
                                       color = "black", face = "bold", size = 10))


rm(list=ls()[! ls() %in% var_keep])

#fig1c/d/eTricycle cell cycle plot ----

lgd <- tricycle::circle_scale_legend(addStageLabel = T, y.inner = 0.5,
                                     y.outer = 3.25, text.size = 8, y.text = 4.25) +
  theme(panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        plot.margin = margin(0,4,0,0, "pt"))

lgd$layers[[2]]$aes_params$alpha <- 1
lgd$layers[[2]]$aes_params$size <- 3

ggsave(here("plots/6mo_manuscript/tricycle_estimate_logo.png"), device = "png",
       height = 5, width = 5, units = "in", bg = "transparent")



#tricyle key
img <- png::readPNG(here("plots/6mo_manuscript/tricycle_estimate_logo.png"))
g <- rasterGrob(img, interpolate=TRUE)
# transparent_color <- "#FFFFFF00"
# 
# test <- g
# test$raster <- g$raster %>% str_replace("#FFFFFFFF", "#FFFFFF00") %>% as.matrix(., nrow = 239)



#fig1c
cds <- readRDS(here("tc_neuro_p20.rds"))
pData(cds)$figure_celltype <- pData(cds)$figure_celltype %>%
  recode("NENs" = "NC-neurons", "neuroglia/neuroendocrine" = "NC-neurons", "Neuroglia" = "NC-glia") %>% 
  as.factor() %>% 
  relevel("MENs", "NC-glia", "NC-neurons")


cc_define <- list("lower" = c(0, pi/2), "upper" = c(3*pi/2,2*pi))

#Embedding
pData(cds)$dummy = "All"
hue.colors = c("#2E22EA", "#9E3DFB", "#F86BE2", "#FCCE7B", 
               "#C4E416", "#4BBA0F", "#447D87", "#2C24E9")
hue.n = 500

#define angle for arc
ang = seq(0, pi/4, length.out = 100)
#define radius (e.g. 1/3 = .33)
ang_df <- data.frame(x = cos(ang)/3, y = sin(ang)/3)




p1 <- plot_cells_mod(cds, reduction_method = "tricycleEmbedding",
               color_cells_by = "tricyclePosition", cell_size = 1, label_cell_groups = F, ) + 
  theme_figure() +
  geom_vline(xintercept = 0, linetype= "dashed", size = .3) +
  geom_hline(yintercept = 0, linetype = "dashed", size = .3) + theme(legend.position = "none") +
  ggtitle("") + 
  facet_wrap(~dummy) +
  ylab("") +
  theme(strip.text.x = element_text(size = 8),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_color_gradientn(name = "tricyclePosition", limits = range(0,2 * pi),
                        breaks = seq(from = 0, to = 2 * pi, length.out = hue.n),
                        colors = hue.colors, guide = "none") +
  geom_line(data = ang_df, aes(x = x, y = y)) + 
  geom_path(data = data.frame(x = c(0,.35), y = c(0, .35)), aes(x,y)) +
  geom_path(data = data.frame(x = c(0,.45), y = c(0, 0)), aes(x,y)) + 
  annotate(geom = "text", label = expression(theta), x = .4, y = .125) + 
  #annotate("rect", xmin=0, xmax=.5, ymin=-1.3, ymax=.9, alpha=0.1, fill="black") +
  coord_cartesian(xlim =c(-2.4, 0.5), ylim= c(-1.85,.9), expand = F) +
  #annotate("rect", xmin=-.1, xmax=1, ymin=-2, ymax=-.45, fill="white") +
  annotation_custom(g, xmin = -.9, xmax = .45, ymin = -2.15, ymax= -.2)
  #geom_path(data = data.frame(x = c(-.45,-.45), y = c(-2, -.5)), aes(x,y)) +
  #geom_path(data = data.frame(x = c(-.45,.5), y = c(-.5, -.5)), aes(x,y)) 
  
  # annotation_custom(g, xmin = -1, xmax = 0, ymin = -1, ymax= 0)


fig1c_key <- ggplot() + theme_void() +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

#Embedding faceted by celltype
p2 <- plot_cells_mod(cds, reduction_method = "tricycleEmbedding",
               color_cells_by = "figure_celltype", cell_size = 1, label_cell_groups = F, ) +
  theme_figure() +
  facet_wrap(~figure_celltype, ncol = 1) +
  scale_color_manual(values = colors) +
  geom_vline(xintercept = 0, linetype= "dashed", size =.3) +
  geom_hline(yintercept = 0, linetype = "dashed", size = .3) + theme(legend.position = "none") + 
  theme(strip.text.x = element_text(size = 8),
        panel.grid = element_blank(),
        axis.title.y = element_text(hjust = 1)) +
  ggtitle("") + 
  annotate("rect", xmin=0, xmax=.5, ymin=-1.3, ymax=.9, alpha=0.1, fill="black") +
  coord_cartesian(xlim =c(-2.4, 0.5), ylim= c(-1.3,.9), expand = F)

fig1c <- ggarrange(p1, ggplot() + theme_void(), p2, heights = c(1.3, -.05, 3), ncol = 1, align = "v")


#fig1d expression of genes over cell cycle
genes_sn <- c("Top2a", "Ect2")
cds_sub <- cds[fData(cds)$gene_short_name %in% genes_sn]
rownames(cds_sub) <- fData(cds_sub)$gene_short_name

dat_mat <- assay(cds_sub, "logcounts") %>% as.matrix() %>% t() %>%
  cbind(., as.data.frame(pData(cds_sub)[,c("tricyclePosition", "figure_celltype")]))

exp_plot <- map(genes_sn, function(gn){
  ggplot(dat_mat, aes_string(x = "tricyclePosition", y = gn)) +
    geom_scattermore(aes_string(x = "tricyclePosition", y = gn, color = "figure_celltype"), pointsize = 3.2) +
    facet_wrap(~figure_celltype) +
    geom_smooth(color = "black", alpha = 0.5) +
    ggtitle("") +
    ylab(glue::glue(gn, " (log10-norm exp.)")) +
    scale_x_continuous() + 
    theme_bw() + 
    scale_color_manual(values = colors) +
    theme(legend.position = "none") + 
    scale_x_continuous(limits = c(0, 2 * pi),
                       breaks = c(0,pi/2, pi, 3 * pi/2, 2 * pi),
                       labels = paste0(c(0, 0.5, 1, 1.5, 2), "π"), name = "θ") +
    theme(panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank())
  
})

fig1d <- do.call(ggarrange, c(exp_plot))

#fig1e Tricycle density plot
#Continuous
fig1e <- tricycle::plot_ccposition_den(cds$tricyclePosition,cds$figure_celltype,
                              '', type = "linear",bw = 10, palette.v = colors,
                              fig.title = "", line.size = 1) +
  theme(plot.title = element_blank()) +
  annotate("rect", xmin=cc_define$lower[1], xmax=cc_define$lower[2], ymin=-.025, ymax=.6, alpha=0.1, fill="black") +
  annotate("rect", xmin=cc_define$upper[1], xmax=cc_define$upper[2], ymin=-.025, ymax=.6, alpha=0.1, fill="black") +
  scale_x_continuous(limits = c(0, 2 * pi), breaks = c(0, pi/2, pi, 3 * pi/2, 2*pi),
                     labels = paste0(c(0,0.5, 1, 1.5, 2), "π"), name = "", expand = c(0,0)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0))) +
  theme(legend.position = c(0.5, 0.8),
        legend.direction = "horizontal",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.y = element_text(size = 10)) + 
  annotation_custom(
        grob = textGrob(label = "Non-cycling", hjust = 0, gp = gpar(cex = 1.1, col = "grey50")),
        ymin = .525,      # Vertical position of the textGrob
        ymax =.625,
        xmin = .25,         # Note: The grobs are positioned outside the plot area
        xmax = .35) + 
  annotation_custom(
    grob = textGrob(label = "Cycling", hjust = 0, gp = gpar(cex = 1, col = "black")),
    ymin = .525,      # Vertical position of the textGrob
    ymax =.625,
    xmin = 2.65,         # Note: The grobs are positioned outside the plot area
    xmax = 3.1) +
  annotation_custom(
      grob = textGrob(label = "Non-cycling", hjust = 0, gp = gpar(cex = 1, col = "grey50")),
      ymin = .525,      # Vertical position of the textGrob
      ymax =.625,
      xmin = 5,         # Note: The grobs are positioned outside the plot area
      xmax = 5.2) 
  
#drop dotted line for average
fig1e$layers <- fig1e$layers[-c(2)]
  

#fig1f proportion of cell cycling
cycling_table <- read.csv(here("results/supp_figures/cycling_proportions_neuro_p20.csv"), row.names = 1) %>%
  rownames_to_column(var = "figure_celltype")

cycling_table$figure_celltype <-cycling_table$figure_celltype %>% 
  recode("NENs" = "NC-neurons", "Neuroglia" = "NC-glia") %>% 
  as.factor() %>% 
  relevel("MENs", "NC-glia", "NC-neurons")


fig1f_proto <- ggplot(cycling_table) +
  geom_col(aes(x = figure_celltype , y = prop_cycling, fill = figure_celltype)) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(breaks = c(0, .1, .2, .3, .4), labels = paste0(c(0, 10, 20, 30, 40), "%")) +
  theme_bw() +
  theme(panel.grid.major.x  = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(vjust = 2),
        axis.title.x = element_text(vjust = -.5)) + 
  ylab("Proportion of cells cycling") + 
  xlab("Cell Type") +
  theme(legend.position= "none")

fig1f <- fig1f_proto %>% ggarrange(ggplot() + theme_void(), ., widths = c(0.1, 1), nrow =1) 


###Make fig panel----
blank_space <- ggplot() + theme_void()

umaps <- ggarrange(fig1a, fig1b, labels = c("a", "b"), ncol = 1)
fig1c <- ggarrange(fig1c, labels = "c")

row1 <- ggarrange(umaps, fig1c, widths = c(1, 0.45))
row2 <- ggarrange(fig1d, labels = "d")
row3 <- ggarrange(fig1e,fig1f, widths = c(1, 0.5), labels = c("e", "f"), nrow =1)


cairo_pdf(here("plots/6mo_manuscript/figure1_panel_v4.pdf"), height = 14, width = 10)
ggarrange(row1,  row2, row3, ncol =1, heights = c(1, .33, .33))
dev.off()
