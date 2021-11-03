library(SingleCellExperiment)
library(monocle3)
library(here)
library(dplyr)
library(ggplot2)
library(purrr)
source(here("scripts/accessory_functions/prep_sparklines_data.R"))
library(scrattch.vis)

cds <- readRDS(here("6month_LMMP.rds"))

cell_types_to_keep <- c("MENs","NENs","Neuroglia")

#if needed, rename cell types for figure
cds_sub <- cds[,pData(cds)$cell_type %in% cell_types_to_keep]
colnames(cds_sub) <- paste(pData(cds_sub)$barcode, pData(cds_sub)$sample, sep = ".")

pData(cds_sub)$cell_type <- pData(cds_sub)$cell_type %>%
  recode("NENs" = "NC-neurons", "Neuroglia" = "NC-glia")


pData(cds_sub) <- pData(cds_sub) %>% as.data.frame() %>%
  group_by(cell_type) %>%
  mutate(group_title = paste0(cell_type, " (n = ", n() ,")")) %>%
  ungroup() %>%
  DataFrame()
  
#
genes <- c("Calcb","Aebp1","Cdh3","Cftr","Clic3","Fmo2","Ntf3","Slpi","Smo","Myl7","Met","Il18","Slc17a9","Ret","Nos1","Sox10", 
           "Stmn2","Elavl1","Elavl2","Elavl3","Elavl4","Ncam1","Stx1a","Stx1b","Vsnl1","Hand2", "Pde10a")
genes_list <- shorten_gene_lists(genes, 8)

violins <- map(genes_list, function(gns){
  vl <- plot_genes_violin(cds_sub[fData(cds)$gene_short_name %in% gns,],
                          group_cells_by = "cell_type",
                          pseudocount = 1, log_scale = T, normalize = T,
                          panel_order = gns, min_expr = 0) + 
    xlab("") + 
    theme(axis.text.x = element_text(size = 11)) + 
    scale_fill_manual(values = rev(wesanderson::wes_palette("GrandBudapest1", n =3)))
})

pdf(here("plots/supp_figures/monocle_violin_expression_log.pdf"), width = 8, height = 8) 
  violins
dev.off()

violins <- map(genes_list, function(gns){
  vl <- plot_genes_violin(cds_sub[fData(cds)$gene_short_name %in% gns,],
                          group_cells_by = "cell_type",
                          pseudocount = 0, log_scale = F, normalize = T,
                          panel_order = gns, min_expr = 0) + 
    xlab("") + 
    theme(axis.text.x = element_text(size = 11)) + 
    scale_fill_manual(values = rev(wesanderson::wes_palette("GrandBudapest1", n =3)))
})

pdf(here("plots/supp_figures/monocle_violin_expression.pdf"), width = 8, height = 8) 
violins
dev.off()

##Scratch.vis violins or quasiviolin
labels <- c("MENs","NC-glia","NC-neurons")
colors <- rev(wesanderson::wes_palette("GrandBudapest1", n =3))

colnames(cds_sub) <- paste(pData(cds_sub)$barcode, pData(cds_sub)$sample, sep=".")

dots <- map(genes_list, function(gns){
  print(gns)
  annotDF <- prep_annotation_matrix(cds_sub[fData(cds)$gene_short_name %in% gns,], facet_by = "cell_type", colors = colors, labels = labels)
  matDF <- prep_expression_matrix(cds_sub[fData(cds)$gene_short_name %in% gns,], norm_method = "none",exprs_assay = "logcounts" , exprs_mat_gene_names = "gene_id", short_gene_names = "gene_short_name")
  
  
  scrattch.vis::group_violin_plot(matDF, annotDF, genes = gns,
                                  grouping = "primary_type",
                                  label_height = 15,
                                  font_size = 10, 
                                  max_width = 15)
})

pdf(here("plots/supp_figures/scratch_violin_expression.pdf"), width = 9, height = 8) 
  dots
dev.off()

dots <- map(genes_list, function(gns){
  print(gns)
  
  annotDF <- prep_annotation_matrix(cds_sub[fData(cds)$gene_short_name %in% gns,], 
                                    facet_by = "cell_type", colors = colors, labels = labels)
  matDF <- prep_expression_matrix(cds_sub[fData(cds)$gene_short_name %in% gns,],
                                  norm_method = "none",exprs_assay = "logcounts" ,
                                  exprs_mat_gene_names = "gene_id", short_gene_names = "gene_short_name")
  
  
  scrattch.vis::group_quasirandom_plot(matDF, annotDF, genes = gns,
                                  grouping = "primary_type",
                                  label_height = 15, 
                                  font_size = 10, 
                                  max_width = 15)
})

pdf(here("plots/supp_figures/scratch_quasiviolin_expression.pdf"), width = 9, height = 8) 
dots
dev.off()
