##Visualize MENS + neuronal markers in Vanderbilt data
library(monocle3)
library(scrattch.vis)
library(here)
library(tidyverse)

source(here("scripts/accessory_functions/prep_sparklines_data.R"))


cds <- readRDS(here("vanderbilt/data/cds.vanderbilt.10x_inDrop_annotated.rds"))

pData(cds) <- pData(cds) %>%
  as.data.frame() %>%
  mutate(Group = factor(case_when(
    vanderbilt_clusters == 10 ~ "Putative MENs",
    vanderbilt_clusters != 10 ~ "Enteric neurons"))) %>% 
      DataFrame()

#order cells by decreasing number of genes expressed
idx <- pData(cds)$log10UMI %>% order(., decreasing = F)
cds <- cds[,idx]

pData(cds) <- pData(cds) %>% as.data.frame() %>%
  group_by(cell_type_factor) %>%
  mutate(group_title = paste0(cell_type_factor, " (n = ", n() ,")")) %>%
  ungroup() %>%
  DataFrame()

genes <- c("Calcb","Cdh3","Aebp1","Cftr","Clic3","Fmo2","Ntf3","Slpi","Smo","Il18","Upk3b","Pde10a","Elavl2","Elavl3","Elavl4","Stx3", 
           "Stmn2","Vsnl1","Nos1","Snap25","Malat1")

# genes_add <- scan(text = "Prdm8, Stmn2, Gabra1, Gpr88, Plcxd3, Vsnl1, Icam5, Kcnf1, Cck, Asph, Myo5b, Trhde, Camk2b, Lpl, Ttc9, Napb, Cyb561, Cxadr, Unc13a, C1qtnf4, Elavl2, 2010300C02Rik, Sertad4, Ckmt1, Aak1, Erc2, Cdk5r1, Sema3c, Slc4a10, Pnma2, Nptxr, Ablim3, Gabpb2, Nrip3, Sema3a, Runx1t1, Ptprk, Rasgef1a, Egr3, AW551984, Fabp3, Pfkp, Epha3, Kalrn, Nrxn3, Dync1i1, Ntf3, Itpka, Lmo7, Kcnab1, Fgf9, Slc2a3, Pde10a, Cacna2d3, Nrgn, Camta1, Htr1b, St8sia2, Rian, Tmod3, Slco5a1, Isl1, Neto1, Mctp1, Olfm1, Foxp1, Tmsb10, Cdkl2, Grb10, Tubb2b, Snap47, Snap29, Snap23, vamp3, vamp8, stx1a, vamp5, stx8, stx12, stx11, stx1b, napa, smn1, sod1, smpd1, bdnf, ntf5, fig4, Lrrk2, eno2, nr4a2, map2, Nav2, Park7, pomc, Kif1c, Ddr1, Ddr2, Gria3, Pink1",
#                   sep = ",", what = "character") %>%
#   gsub(" ", "", ., fixed = TRUE) %>%
#   str_to_title()
# 
# genes_not_found <- genes_add[which(genes_add %in% fData(cds)$gene_short_name == F)]
# print(genes_not_found)
# 
# genes_add <- genes_add[!(genes_add %in%genes_not_found)]
# gene_idx <- map(genes_add, function(gene){which(fData(cds)$gene_short_name == gene)}) %>% unlist()
# 
# genes_all <- c(genes, genes_add)


labels <- pData(cds)$Group %>% unique %>% sort %>% as.character()# MENS glia NENs
colors <- wesanderson::wes_palette(name = "GrandBudapest1")[c(1,4)]

annotDF <- prep_annotation_matrix(cds, facet_by = "Group", colors = colors, labels = labels)
matDF <- prep_expression_matrix(cds, norm_method = "none",exprs_assay = "logcounts" , exprs_mat_gene_names = "gene_id", short_gene_names = "gene_short_name")

pl <- scrattch.vis::sample_bar_plot(matDF, annotDF,
                                    genes = genes,
                                    grouping = "primary_type",
                                    bg_color ="#f7f7f7",
                                    label_height = 15,
                                    label_type= "angle",
                                    font_size = 12, return_type = "plot") 

pl_vio <- scrattch.vis::group_quasirandom_plot(matDF, annotDF, genes = genes,
                                               grouping = "primary_type", label_height = 15, font_size = 10, return_type = "plot", max_width = 5)

#Now make sparklines plot
pdf(here("plots/supp_figures/vanderbilt_selected_genes_sparklines.pdf"))
       pl
dev.off()

pdf(here("plots/supp_figures/vanderbilt_selected_genes_quasiviolin.pdf"), width = 7, height = 15)
pl_vio
dev.off()
