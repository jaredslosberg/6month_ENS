library(monocle3)
library(tricycle)
library(SingleCellExperiment)
library(scater)
library(cowplot)
library(ggplot2)
library(here)

object_names <- c("TC_LMMP", "TC_MENS")

purrr::map(object_names, function(name){ 


  pdf(here(paste0("./plots/",name,"_tricycle_assignments.pdf")), width = 10, height = 10)
  
  lmmp <- readRDS(here(paste0(name,".rds")))
  rownames(lmmp) <- fData(lmmp)[,"gene_id_trimmed"]
  counts_list <- assays(lmmp)
  
  #create SCE and run tricycle workflow ----
  sce_lmmp <- SingleCellExperiment(assays = counts_list)
  sce_lmmp <- scuttle::logNormCounts(sce_lmmp)
  reducedDim(sce_lmmp, "UMAP") <- reducedDim(lmmp, "UMAP")
  colData(sce_lmmp) <- colData(lmmp)
  
  
  
  sce_lmmp <- tricycle::project_cycle_space(sce_lmmp, exprs_values = "logcounts")
  
  #plot pca
  p1 <- scater::plotReducedDim(sce_lmmp, dimred = "tricycleEmbedding") +
    labs(x = "Projected PC1", y = "Projected PC2") +
    ggtitle(sprintf("Projected cell cycle space (n=%d)",
                    ncol(sce_lmmp))) + 
    theme_bw(base_size = 14)
  
  sce_lmmp <- estimate_cycle_position(sce_lmmp)
  
  #plot pca
  p2 <- scater::plotReducedDim(sce_lmmp, dimred = "tricycleEmbedding", colour_by = "tricyclePosition") +
    labs(x = "Projected PC1", y = "Projected PC2") +
    ggtitle(sprintf("Projected cell cycle space (n=%d)",
                    ncol(sce_lmmp))) + 
    theme_bw(base_size = 14)
  
  # #schwabe
  # sce_lmmp <- estimate_cycle_stage(sce_lmmp, gname.type = 'ENSEMBL', species = 'mouse', exprs_values = "norm_counts")
  # 
  # scater::plotReducedDim(sce_lmmp, dimred = "tricycleEmbedding", colour_by = "CCStage") +
  #   labs(x = "Projected PC1", y = "Projected PC2") +
  #   ggtitle(sprintf("Projected cell cycle space (n=%d)",
  #                   ncol(sce_lmmp))) + 
  #   theme_bw(base_size = 14)
  
  ##cylic scale
  p <- plot_emb_circle_scale(sce_lmmp, dimred = 1, point.size = 3.5, point.alpha = 0.9) +
    theme_bw(base_size = 14)
  legend <- circle_scale_legend(text.size = 5, alpha = 0.9)
  p3 <- plot_grid(p, legend, ncol = 2, rel_widths = c(1, 0.4))
  
  pca_emb_pos <- reducedDim(sce_lmmp, "tricycleEmbedding")
  tric_pos <- as.data.frame(colData(sce_lmmp))[,"tricyclePosition"]
  
  #Add back to cds object
  reducedDim(lmmp, "tricycleEmbedding") <- pca_emb_pos
  pData(lmmp)[,"tricyclePosition"] <- tric_pos
  
  p<-plot_emb_circle_scale(sce_lmmp, dimred = "UMAP", color_by = "tricyclePosition")
  p4 <- plot_grid(p, legend, ncol = 2, rel_widths = c(1, 0.1))

  print(p1)
  print(p2)
  print(p3)
  print(p4)
  dev.off()
  
  
  saveRDS(lmmp, here(paste0(name,".rds")))
})

# ###Back to cds, do followup ----
# pDat <- as.data.frame(pData(lmmp))
# 
# ggplot(pDat) + geom_violin(aes(x = cell_type, y= tricyclePosition, color = cell_type)) +
#   coord_flip()
# 
# dev.off()