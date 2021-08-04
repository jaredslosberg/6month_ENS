library(monocle3)
library(tricycle)
library(SingleCellExperiment)
library(scater)
library(cowplot)
library(ggplot2)
library(here)


tabula_muris <- readRDS(here("tabula_muris/data/10x_tabula_muris_annotated.rds"))

tabula_muris_sub <- tabula_muris[!is.na(fData(tabula_muris)$ensembl_id_trimmed),]
rownames(tabula_muris_sub) <- fData(tabula_muris_sub)[,"ensembl_id_trimmed"]
counts_list <- assays(tabula_muris_sub)
  
# run tricycle workflow ----
tabula_muris_sub <- scuttle::logNormCounts(tabula_muris_sub)

  
  
  
tabula_muris_sub <- tricycle::project_cycle_space(tabula_muris_sub, exprs_values = "logcounts")
  
#plot pca
p1 <- scater::plotReducedDim(tabula_muris_sub, dimred = "tricycleEmbedding") +
    labs(x = "Projected PC1", y = "Projected PC2") +
    ggtitle(sprintf("Projected cell cycle space (n=%d)",
                    ncol(tabula_muris_sub))) + 
    theme_bw(base_size = 14)
  
tabula_muris_sub <- estimate_cycle_position(tabula_muris_sub)
  
#plot pca
p2 <- scater::plotReducedDim(tabula_muris_sub, dimred = "tricycleEmbedding", colour_by = "tricyclePosition") +
    labs(x = "Projected PC1", y = "Projected PC2") +
    ggtitle(sprintf("Projected cell cycle space (n=%d)",
                    ncol(tabula_muris_sub))) + 
    theme_bw(base_size = 14)
  
 
##cylic scale
p <- plot_emb_circle_scale(tabula_muris_sub, dimred = "UMAP", point.size = 3.5, point.alpha = 0.9) +
    theme_bw(base_size = 14)
legend <- circle_scale_legend(text.size = 5, alpha = 0.9)
p3 <- plot_grid(p, legend, ncol = 2, rel_widths = c(1, 0.4))
  
pca_emb_pos <- reducedDim(tabula_muris_sub, "tricycleEmbedding")
tric_pos <- as.data.frame(colData(tabula_muris_sub))[,"tricyclePosition"]
  

  
  p<-plot_emb_circle_scale(tabula_muris_sub, dimred = "UMAP", color_by = "tricyclePosition")
  p4 <- plot_grid(p, legend, ncol = 2, rel_widths = c(1, 0.1))
  
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  dev.off()
  

#Add back to cds object
reducedDim(tabula_muris, "tricycleEmbedding") <- pca_emb_pos
colData(tabula_muris)[,"tricyclePosition"] <- tric_pos
  
saveRDS(tabula_muris, here("tabula_muris/data/10x_tabula_muris_annotated.rds"))

