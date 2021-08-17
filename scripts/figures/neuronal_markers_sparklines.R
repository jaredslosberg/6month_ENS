library(SingleCellExperiment)
library(monocle3)
library(scrattch.vis)
library(scuttle)
library(here)
library(dplyr)
library(stringr)

source("/home/jared/ENS/Timecourse_ENS/scripts/accessory_functions/prep_sparklines_data.R")

shorten_gene_lists <- function(gene_vec, n){
  gene_list <- list()
  while(length(gene_vec) > n){
    genes <- list(gene_vec[1:n])
    gene_vec <- gene_vec[-c(1:n)]
    gene_list <- c(gene_list, genes)
  }
  gene_list <- c(gene_list, list(gene_vec))
  return(gene_list)
}     

cds <- readRDS(here("./6month_LMMP.rds"))
assay(cds, "logcounts") <- normalized_counts(cds)

cell_types_to_keep <- pData(cds)$cell_type %>% unique()

cds_sub <- cds[,pData(cds)$cell_type %in% cell_types_to_keep]
pData(cds_sub)$cell_type <- pData(cds_sub)$cell_type %>%
  recode("MENS" = "MENs","NENS" = "NC-derived cells", "Neuroglia" = "NC-derived cells")
  

#order cells by decreasing number of genes expressed
idx <- pData(cds_sub)$num_genes_expressed %>% order(., decreasing = F)
cds_sub <- cds_sub[,idx]

#reformats cds so genes are in this order
genes <- scan(text = "Prdm8, Stmn2, Gabra1, Gpr88, Plcxd3, Vsnl1, Icam5, Kcnf1, Cck, Asph, Myo5b, Trhde, Camk2b, Lpl, Ttc9, Napb, Cyb561, Cxadr, Unc13a, C1qtnf4, Elavl2, 2010300C02Rik, Sertad4, Ckmt1, Aak1, Erc2, Cdk5r1, Sema3c, Slc4a10, Pnma2, Nptxr, Ablim3, Gabpb2, Nrip3, Sema3a, Runx1t1, Ptprk, Rasgef1a, Egr3, AW551984, Fabp3, Pfkp, Epha3, Kalrn, Nrxn3, Dync1i1, Ntf3, Itpka, Lmo7, Kcnab1, Fgf9, Slc2a3, Pde10a, Cacna2d3, Nrgn, Camta1, Htr1b, St8sia2, Rian, Tmod3, Slco5a1, Isl1, Neto1, Mctp1, Olfm1, Foxp1, Tmsb10, Cdkl2, Grb10, Tubb2b, Snap47, Snap29, Snap23, vamp3, vamp8, stx1a, vamp5, stx8, stx12, stx11, stx1b, napa, smn1, sod1, smpd1, bdnf, ntf5, fig4, Lrrk2, eno2, nr4a2, map2, Nav2, Park7, pomc, Kif1c, Ddr1, Ddr2, Gria3, Pink1",
     sep = ",", what = "character") %>%
  gsub(" ", "", ., fixed = TRUE) %>%
  str_to_title()

genes_not_found <- genes[which(genes %in% fData(cds)$gene_short_name == F)]
print(genes_not_found)

genes <- genes[!(genes %in%genes_not_found)]
gene_idx <- map(genes, function(gene){which(fData(cds_sub)$gene_short_name == gene)}) %>% unlist()

cds_sub <- cds_sub[gene_idx,]



pData(cds_sub) <- pData(cds_sub) %>% as.data.frame() %>%
  group_by(cell_type) %>%
  mutate(group_title = paste0(cell_type, " (n = ", n() ,")")) %>%
  ungroup() %>%
  DataFrame()
  
colnames(cds_sub) <- paste(pData(cds_sub)$barcode, pData(cds_sub)$sample, sep = ".")

#choose color sets. Previously, neurons were salmon and MENS are green. 
labels <- pData(cds_sub)$group_title %>% unique %>% sort # Glia MENS NENs
colors <- scales::hue_pal()(length(labels))


annotDF <- prep_annotation_matrix(cds_sub, facet_by = "group_title", colors = colors, labels = labels)

matDF <- prep_expression_matrix(cds_sub, norm_method = "none",exprs_assay = "logcounts" , exprs_mat_gene_names = "gene_id", short_gene_names = "gene_short_name")


#make list of genes with a maximum vector length of "max_genes_per_page" so that plot is readable
max_genes_per_page <- 50

if(length(genes) > max_genes_per_page){
  genes_list <- shorten_gene_lists(genes, max_genes_per_page)
}else{
  genes_list <- list(genes_vector)
}

exprs_bar_plots <- purrr::map(genes_list, function(genes_vector_sub){
  pl <- scrattch.vis::sample_bar_plot(matDF, annotDF,
                                      genes = genes_vector_sub,
                                      grouping = "primary_type",
                                      bg_color ="#f7f7f7",
                                      label_height = 17,
                                      font_size = 12,
                                      max_width = 12)        
  pl_titled <- pl + 
    ggtitle("Sparkline expression for neuronal genes") +
    theme(plot.title = element_text(color = "black", size = 16, hjust = 0.5))
})

pdf(here("./plots/supp_figures/neuronal_genes_sparklines.pdf"), width = 8, height = 16)
  exprs_bar_plots
dev.off()
