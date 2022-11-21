##projections were done independently, but should be a way to merge the two technologies into one embedding
##merge the cds's here. Count data might not be comparable but the embedding in projected patterns will be

source(here("./scripts/accessory_functions/monocle_mods.R"))
source(here("./scripts/accessory_functions/patternFeatureCorrelationHeatmap.R"))

library(monocle3)
library(dplyr)
library(SingleCellExperiment)
library(uwot)
library(here)
library(ggplot2)

droplet_tabula_muris <- readRDS(here("./tabula_muris/data/10x_tabula_muris_annotated.rds"))
pData(droplet_tabula_muris)$technology <- "10x"
facs_tabula_muris <- readRDS(here("./tabula_muris/data/smartseq2_tabula_muris_annotated.rds"))
pData(facs_tabula_muris)$technology <- "smartseq2"

droplet_projected_usages <- pData(droplet_tabula_muris) %>%
  as.data.frame() %>% 
  select(starts_with("lmmp_6mo_11_1_pattern")) %>%
  map_df(., function(pattern){
    (pattern - mean(pattern)) / sd(pattern)
  })

facs_projected_usages <- pData(facs_tabula_muris) %>%
  as.data.frame() %>% 
  select(starts_with("lmmp_6mo_11_1_pattern")) %>%
  map_df(., function(pattern){
    (pattern - mean(pattern)) / sd(pattern)
  })


##Merge cell data sets, they have the same genes
merged_tabula_muris <- monocle3::combine_cds(list(droplet_tabula_muris, facs_tabula_muris),
                      cell_names_unique = T, 
                      sample_col_name = "merged_by")

merged_projected_usages <- bind_rows(droplet_projected_usages, facs_projected_usages)

projected_embedding <- uwot::umap(X = merged_projected_usages, n_neighbors = 15)
rownames(projected_embedding) <- rownames(merged_projected_usages)
colnames(projected_embedding) <- paste0("projectedUMAP ", c(1:2))

reducedDim(merged_tabula_muris, "projectedUMAP") <- projected_embedding

pdf(here("tabula_muris/plots/merged_projection_UMAP_normalized.pdf"), width = 15, height = 10)
  plot_cells_mod(merged_tabula_muris, reduction_method = "projectedUMAP", color_cells_by = "technology", label_cell_groups = F)
  plot_cells_mod(merged_tabula_muris, reduction_method = "projectedUMAP", color_cells_by = "tissue", group_label_size = 4)
dev.off()

pdf(here("tabula_muris/plots/merged_projection_UMAP_facet_tissue_normalized.pdf"), width = 15, height = 10)
  plot_cells_mod(merged_tabula_muris, reduction_method = "projectedUMAP", color_cells_by = "technology", label_cell_groups = F) + 
    facet_wrap(vars(tissue))
dev.off()


##correlation of pattern usages
feature <- "technology"
technology_heatmap <- patternFeatureCorrelationHeatmap(cds = merged_tabula_muris,
                                                       cellWeights.df = merged_projected_usages,
                                                       features = feature, 
                                                       pattern_prefix = "lmmp_6mo_11_1_pattern")

pdf(here("./tabula_muris/plots/tabula_muris_technology_feature_correlations.pdf"), width = 10, height = 8)
  technology_heatmap
dev.off()

cor_mat <- cor(merged_projected_usages)
rownames(cor_mat) <- rownames(cor_mat) %>% stringr::str_split_fixed(., "_", 5) %>% .[,5]
colnames(cor_mat) <- colnames(cor_mat) %>% stringr::str_split_fixed(., "_", 5) %>% .[,5]
pdf(here("./tabula_muris/plots/tabula_muris_feature_cocorrelations.pdf"), width = 10, height = 8)
  ComplexHeatmap::Heatmap(cor_mat)
dev.off()

saveRDS(merged_tabula_muris, "./tabula_muris/merged_tabula_muris_annotated.rds")
