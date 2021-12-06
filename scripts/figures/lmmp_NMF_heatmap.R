#' ---
#' title: "6mo LMMP pattern identification"
#' author: "Jared Slosberg"
#' date: "Nov 1, 2020"
#' ---
#' 

#Generate heatmap for cell pattern weights, annotated by cell type

library(here)
library(monocle3)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(viridis)
source(here("scripts/accessory_functions/pattern_plotting.R"))
source(here("scripts/accessory_functions/monocle_mods.R"))
source(here("scripts/accessory_functions/geneSetEnrichment.R"))


lmmp_6mo <- readRDS(here("6month_LMMP.rds"))

pattern_usages <- read.csv(here("results/NMF/lmmp/old_pattern_run/50dims/pattern_cell_weights.csv"),row.names = 1)

if(all(rownames(pData(lmmp_6mo)) == rownames(pattern_usages))){
  pData(lmmp_6mo) <- cbind(pData(lmmp_6mo), pattern_usages)
}

pData(lmmp_6mo)$cell_type_mod <- pData(lmmp_6mo) %>% as.data.frame %>% dplyr::pull(cell_type) %>%
  recode("Penk+ Fibroblasts" = "Fibroblast-A",
         "Neuroglia" = "NC-glia",
         "NENs" = "NC-neurons",
         "Macrophage-A" = "Macrophage-A",
         "Macrophage-B" = "Macrophage-B",
         "Macrophage-C" = "Macrophage-C",
         "Pdgfra+ Fibroblasts" = "Fibroblast-B",
         "Smooth muscle cells" = "SMC",
         "Smooth muscle cells B" = "SMC-B",
         "Vascular endothelium" = "V Endothelium")

#order by decreasing number of cells
factor_order_by_size <- pData(lmmp_6mo)$cell_type_mod %>%
  table() %>%
  sort(decreasing = T) %>%
  names()

#manual order
factor_order <- c("MENs","NC-glia","NC-neurons","Macrophage-A","Macrophage-B","V Endothelium","SMC","RBC","Fibroblast-B",
"B Lymphocytes","Fibroblast-A","Macrophage-C","T cells","NK cells","SMC-B","Unknown")


#Change this to something like cluster, cell type... whatever cell metadata to group by
group_feature <- "cell_type_mod"

#define color scale
colfun <- circlize::colorRamp2(breaks = quantile(as.matrix(pattern_usages), c(0,.75,.98)),
                               colors = viridis::viridis(3))

#Subset to cell name and feature, group of feature to establish column rows
cells_by_group <- as.factor(pData(lmmp_6mo)[,group_feature])
names(cells_by_group) <- colnames(lmmp_6mo)

#order by UMI
umi_order <- pData(lmmp_6mo)[,"Total_UMIs", drop = F] %>% as.data.frame() %>% arrange(desc(Total_UMIs)) %>% rownames()
cells_by_group <- cells_by_group[umi_order]

#order cells by their cluster (Decreasing size for colors), this object is a named integer array
cells_by_group <- cells_by_group[order(match(cells_by_group,factor_order_by_size))] 
cells_order <- names(cells_by_group)

#Plot heatmap of H (cell weight) matrix
#Assign a color for each group, these match default ggplot colors
my_colors <- scales::hue_pal()(length(unique(cells_by_group)))
names(my_colors) <- unique(cells_by_group)

#order cells by manually specified order, this object is a named integer array
#Colors are assigned by size (i.e. will be consistent with other representations)
cells_by_group <- cells_by_group[order(match(cells_by_group,factor_order))] 
cells_order <- names(cells_by_group)


#Define cell groups to include mark annotations for
cell_type_to_mark <- ""
cell_type_to_mark_idx <- sapply(cell_type_to_mark, function(ct){
  floor(median(which(cells_by_group == ct)))
  })

ha_cellgroup <- HeatmapAnnotation(annot = anno_mark(at = cell_type_to_mark_idx,
                                                    labels = cell_type_to_mark,
                                                    labels_rot = 0),
                                  cell_type = cells_by_group,
                                  col = list(cell_type = my_colors),
                                  show_legend = F)                              
                                  
lgd_boxplot = Legend(labels = factor_order, title = "Cell type",
                     legend_gp = gpar(fill = my_colors[factor_order]))

#We want cells grouped by cluster, as above. We also want to do hierarchical clustering
#on the pre-defined cluster level. Info must be in the same order 
cell_mat <- t(pattern_usages)[,cells_order]
stopifnot(colnames(cell_mat) == cells_order)
stopifnot(colnames(cell_mat) == names(cells_by_group))


#change 'cellPatternX' to 'PatternX'
rownames(cell_mat) <- stringr::str_split_fixed(rownames(cell_mat), "cell", 2)[,2]


pdf(here("plots/supp_figures/NMF_heatmap.pdf"), width = 10)

pl <- ComplexHeatmap::Heatmap(cell_mat, name = "Weight", col = colfun,
                        top_annotation = ha_cellgroup,
                        column_order = cells_order,
                        #cluster_columns = cluster_within_group(cell_mat, groups_to_cluster),
                        show_column_names = F,
                        row_dend_width = unit(1, "in"),
                        row_names_gp = gpar(fontsize = 11))

draw(pl, annotation_legend_list = lgd_boxplot)


#row_order = rownames(decomp$H))
dev.off()



