library(SingleCellExperiment)
library(monocle3)
library(scrattch.vis)
library(scuttle)
library(here)
library(dplyr)

source("/home/jared/ENS/Timecourse_ENS/scripts/accessory_functions/prep_sparklines_data.R")


cds <- readRDS(here("./6month_LMMP.rds"))
cds <- scuttle::logNormCounts(cds)


cell_types_to_keep <- c("Neuroglia","NENS","MENS")
cds_sub <- cds[,pData(cds)$cell_type %in% cell_types_to_keep]
pData(cds_sub)$cell_type <- pData(cds_sub)$cell_type %>%
  recode("MENS" = "MENs","NENS" = "NC-derived cells", "Neuroglia" = "NC-derived cells")
  

#order cells by decreasing number of genes expressed
idx <- pData(cds_sub)$num_genes_expressed %>% order(., decreasing = F)
cds_sub <- cds_sub[,idx]

#reformats cds so genes are in this order
genes <- c("Cdh3","Clic3","Ntf3","Snap25")
gene_idx <- sapply(genes, function(gene){which(fData(cds_sub)$gene_short_name == gene)})

cds_sub <- cds_sub[gene_idx,]



pData(cds_sub) <- pData(cds_sub) %>% as.data.frame() %>%
  group_by(cell_type) %>%
  mutate(group_title = paste0(cell_type, " (n = ", n() ,")")) %>%
  ungroup() %>%
  DataFrame()
  
colnames(cds_sub) <- paste(pData(cds_sub)$barcode, pData(cds_sub)$sample, sep = ".")

#choose color sets. Previously, neurons were salmon and MENS are green. 
labels <- pData(cds_sub)$group_title %>% unique %>% sort # Glia MENS NENs
colors <- c("#66B200", "#00BFC4")


annotDF <- prep_annotation_matrix(cds_sub, facet_by = "group_title", colors = colors, labels = labels)
matDF <- prep_expression_matrix(cds_sub, norm_method = "none",exprs_assay = "logcounts" , exprs_mat_gene_names = "gene_id", short_gene_names = "gene_short_name")

#Now make sparklines plot
pl <- scrattch.vis::sample_bar_plot(matDF, annotDF,
                                    genes = genes,
                                    grouping = "primary_type",
                                    bg_color ="#f7f7f7",
                                    label_height = 40,
                                    label_type= "angle",
                                    font_size = 12, return_type = "plot")        
pl_titled <- pl +
  theme(plot.title = element_text(color = "black", size = 16, hjust = 0.5))

pdf(here("./plots/supp_figures/selected_genes_sparklines.pdf"), width = 10, height = 6.5)
  pl_titled
dev.off()
