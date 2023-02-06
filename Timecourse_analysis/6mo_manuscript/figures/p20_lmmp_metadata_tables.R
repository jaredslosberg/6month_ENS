##Create tables for info by cluster for lmmp
##Create tables for NMF pattern definitions and usage
##Create tables for info by cluster for gut cell atlas (nfh subset)


library(tidyverse)
library(monocle3)
library(here)


### lmmp subset metadata----
lmmp <- readRDS(here("tc_lmmp_p20.rds"))

# pData(lmmp) <- pData(lmmp) %>% as.data.frame() %>% 
#   tibble::rownames_to_column(var = "index") %>%
#   left_join(., cell_weights, by = "index") %>%
#   column_to_rownames(var = "index") %>%
#   DataFrame()

pData(lmmp)$cluster <- clusters(lmmp)
lmmp_dat <- pData(lmmp) %>%
  as.data.frame() %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),
    annotated_type = unique(figure_celltype),
    pct_batch1 = round(sum(batch == "1")/n_cells, 3),
    pct_batch2 = round(sum(batch == "0")/n_cells, 3),
    num_genes = round(mean(num_genes_expressed),2),
    mt_ratio = round(mean(mito_ratio), 2),
    umi = round(mean(total_reads),2), 
    umi_sem = round(mean(total_reads)/sqrt(n()), 2)
  )

write.csv(lmmp_dat, here("results/supp_figures/p20_lmmp_metadata_by_cluster.csv"), quote = F, row.names = F)

###Average usage of 6month patterns

cell_weights <- read.csv(here("results/6month_projection_weights/tc_in_6mo_LMMP_old_patterns.csv"))
gene_weights <- read.csv(here("results/6month_projection_weights/pattern_gene_weights.csv"), row.names = 1)

colnames(cell_weights) <- c("cell_id", paste0("sixmo_cellPattern", 1:50))
n_patt <- 50

# pData(lmmp) <- pData(lmmp) %>% as.data.frame() %>% tibble::rownames_to_column(var = "cell_id") %>%
#   left_join(., cell_weights, by = "cell_id") %>%
#   tibble::column_to_rownames(var = "cell_id") %>% DataFrame

#mean and std by cluster 
patt_stats <- pData(lmmp) %>%
  as.data.frame() %>%
  group_by(cluster) %>%
  summarise(across(starts_with("sixmo_cellPattern"), .fns = list(mean = mean, std_dev = sd), .names = "{.col}.{.fn}"))

#add annots
patt_stats <- inner_join(lmmp_dat[,c("cluster","annotated_type", "n_cells")], patt_stats )

#Add ensembl gene_ids 
gene_df <- fData(lmmp) %>% as.data.frame() %>% 
  select(gene_id_trimmed, gene_short_name)
  
gene_weights <- gene_weights %>%
  tibble::rownames_to_column(var = "gene_short_name") %>%
  left_join(gene_df, .)

write.csv(patt_stats, here("results/6month_projection_weights/p20_lmmp_NMF_mean_by_cluster.csv"), quote = F, row.names = F)

#mito and umi for select cell types
lmmp_dat <- pData(lmmp) %>%
  as.data.frame() %>%
  mutate(summ_ct = case_when(
    figure_celltype == "Fibroblasts" ~ "Fibroblasts",
    figure_celltype == "SMC" ~ "SMC",
    figure_celltype == "Endothelial" ~ "Endothelial",
    figure_celltype == "ICC" ~ "ICC",
    figure_celltype == "NC-glia" ~ "Neuroglia",
    figure_celltype == "Macrophage" ~ "Macrophage",
    figure_celltype == "RBC" ~ "RBC", 
    figure_celltype == "NENs" ~ "NENs",
    figure_celltype == "MENs" ~ "MENs",
    figure_celltype == "Enterocytes" ~ "Unknown",
    figure_celltype == "Mucosal epith." ~ "Unknown",
  )) %>%
  select(log10UMI, summ_ct, mito_ratio, batch, num_genes_expressed, total_reads)

lmmp_dat_sum <- lmmp_dat %>%
  group_by(summ_ct) %>%
  summarise(
    n_cells = n(),
    annotated_type = unique(summ_ct),
    pct_batch1 = round(sum(batch == "1")/n_cells, 3),
    pct_batch2 = round(sum(batch == "0")/n_cells, 3),
    num_genes = round(mean(num_genes_expressed),2),
    mt_ratio = round(mean(mito_ratio), 2),
    umi = round(mean(total_reads),2), 
    umi_sem = round(mean(total_reads)/sqrt(n()), 2)
  ) %>% select(-summ_ct) %>%
  relocate(annotated_type)

res1 <- lm(log10UMI ~ summ_ct + 0, lmmp_dat )
anova_res1 <- anova(res1)

lmmp_dat_filt <- lmmp_dat %>% filter(summ_ct %in% c("MENs", "NENs", "Neuroglia"))
res2 <- lm(log10UMI ~ summ_ct + 0, lmmp_dat_filt )
anova_res2 <- anova(res2)

write_csv(lmmp_dat_sum, file = here("results/supp_figures/p20_lmmp_metadata_by_cluster_collapsed.csv"))
write_csv(lmmp_dat, file= here("results/supp_figures/p20_lmmp_metadata_raw.csv"))
