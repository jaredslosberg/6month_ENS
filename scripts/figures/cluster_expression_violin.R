### Visualize expression of IHC-matched and neuronal genes via quasiviolins and sparklines

library(SingleCellExperiment)
library(monocle3)
library(here)
library(dplyr)
library(ggplot2)
library(stringr)
library(purrr)
source(here("scripts/accessory_functions/prep_sparklines_data.R"))



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
genes <- c("Calcb","Aebp1","Cdh3","Cftr","Clic3","Fmo2","Ntf3","Slpi","Smo","Myl7","Met","Il18","Slc17a9","Ret","Nos1","Sox10", 
           "Stmn2","Elavl1","Elavl2","Elavl3","Elavl4","Ncam1","Stx3","Vsnl1","Hand2", "Pde10a")
genes_list <- shorten_gene_lists(genes, 8)  

# vl <- plot_genes_violin(cds_sub[fData(cds)$gene_short_name %in% genes,],
#                         group_cells_by = "group_title",
#                         pseudocount = 1, log_scale = T) + 
#   xlab("") + 
#   theme(axis.text.x = element_text(size = 14)) + 
#   scale_fill_manual(values = rev(wesanderson::wes_palette("GrandBudapest1", n =3)))
# 
# pdf(here("plots/supp_figures/violin_expression.pdf"), width = 8, height = 16) 
#   vl
# dev.off()

###Scratch.vis violins or quasiviolin ----
labels <- c("MENs","NC-glia","NC-neurons")
colors <- rev(wesanderson::wes_palette("GrandBudapest1", n =3))

colnames(cds_sub) <- paste(pData(cds_sub)$barcode, pData(cds_sub)$sample, sep = ".")



  
annotDF <- prep_annotation_matrix(cds_sub[fData(cds)$gene_short_name %in% genes,], 
                                  facet_by = "cell_type", colors = colors, labels = labels)

annotDF_mod <- annotDF
colors_replace <- c("#360000", "#000000","#130B5C") %>% setNames(colors)
annotDF_mod[,"primary_type_color"] <- annotDF_mod[,"primary_type_color"] %>% str_replace_all(.,colors_replace)

matDF <- prep_expression_matrix(cds_sub[fData(cds)$gene_short_name %in% genes,],
                                norm_method = "none",exprs_assay = "logcounts" ,
                                exprs_mat_gene_names = "gene_id", short_gene_names = "gene_short_name")


dots <- scrattch.vis::group_quasirandom_plot(matDF, annotDF, genes = genes,
                                     grouping = "primary_type",
                                     label_height = 5, 
                                     font_size = 10, 
                                     max_width = 15)


pdf(here("plots/supp_figures/quasiviolin_expression.pdf"), width = 9, height = 24) 
dots
dev.off()

#modified color scheme
dots <- scrattch.vis::group_quasirandom_plot(matDF, annotDF_mod, genes = genes,
                                             grouping = "primary_type",
                                             label_height = 5, 
                                             font_size = 10, 
                                             max_width = 15)


pdf(here("plots/supp_figures/quasiviolin_expression_modcolors.pdf"), width = 9, height = 24) 
dots
dev.off()

###Expanded genes lists ----
genes_add <- scan(text = "Prdm8, Stmn2, Gabra1, Gpr88, Plcxd3, Vsnl1, Icam5, Kcnf1, Cck, Asph, Myo5b, Trhde, Camk2b, Lpl, Ttc9, Napb, Cyb561, Cxadr, Unc13a, C1qtnf4, Elavl2, 2010300C02Rik, Sertad4, Ckmt1, Aak1, Erc2, Cdk5r1, Sema3c, Slc4a10, Pnma2, Nptxr, Ablim3, Gabpb2, Nrip3, Sema3a, Runx1t1, Ptprk, Rasgef1a, Egr3, AW551984, Fabp3, Pfkp, Epha3, Kalrn, Nrxn3, Dync1i1, Ntf3, Itpka, Lmo7, Kcnab1, Fgf9, Slc2a3, Pde10a, Cacna2d3, Nrgn, Camta1, Htr1b, St8sia2, Rian, Tmod3, Slco5a1, Isl1, Neto1, Mctp1, Olfm1, Foxp1, Tmsb10, Cdkl2, Grb10, Tubb2b, Snap47, Snap29, Snap23, vamp3, vamp8, stx1a, vamp5, stx8, stx12, stx11, stx1b, napa, smn1, sod1, smpd1, bdnf, ntf5, fig4, Lrrk2, eno2, nr4a2, map2, Nav2, Park7, pomc, Kif1c, Ddr1, Ddr2, Gria3, Pink1",
              sep = ",", what = "character") %>%
  gsub(" ", "", ., fixed = TRUE) %>%
  str_to_title()

genes_not_found <- genes_add[which(genes_add %in% fData(cds)$gene_short_name == F)]
print(genes_not_found)

genes_add <- genes_add[!(genes_add %in%genes_not_found)]
gene_idx <- map(genes_add, function(gene){which(fData(cds_sub)$gene_short_name == gene)}) %>% unlist()

genes_all <- c(genes, genes_add)
genes_all_list <- shorten_gene_lists(genes_all, 8)

dots <- map(genes_all_list, function(gns){
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

pdf(here("plots/supp_figures/quasiviolin_expression_all.pdf"), width = 9, height = 8) 
dots
dev.off()


sparks <- map(genes_all_list, function(gns){
  print(gns)
  
  annotDF <- prep_annotation_matrix(cds_sub[fData(cds)$gene_short_name %in% gns,], 
                                    facet_by = "cell_type", colors = colors, labels = labels)
  matDF <- prep_expression_matrix(cds_sub[fData(cds)$gene_short_name %in% gns,],
                                  norm_method = "none",exprs_assay = "logcounts" ,
                                  exprs_mat_gene_names = "gene_id", short_gene_names = "gene_short_name")
  
  
    scrattch.vis::sample_bar_plot(matDF, annotDF, genes = gns,
                                       grouping = "primary_type",
                                       label_height = 15, 
                                       font_size = 10, 
                                       max_width = 15)
})

pdf(here("plots/supp_figures/sparklines_expression_all.pdf"), width = 9, height = 8) 
sparks
dev.off()