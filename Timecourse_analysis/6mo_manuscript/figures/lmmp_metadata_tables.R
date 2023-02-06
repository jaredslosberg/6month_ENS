##Create tables for info by cluster for lmmp
##Create tables for NMF pattern definitions and usage
##Create tables for info by cluster for gut cell atlas (nfh subset)


library(tidyverse)
library(monocle3)
library(here)


### lmmp subset metadata----
lmmp <- readRDS(here("TC_LMMP.rds"))

# pData(lmmp) <- pData(lmmp) %>% as.data.frame() %>% 
#   tibble::rownames_to_column(var = "index") %>%
#   left_join(., cell_weights, by = "index") %>%
#   column_to_rownames(var = "index") %>%
#   DataFrame()

lmmp_dat <- pData(lmmp) %>%
  as.data.frame() %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),
    annotated_type = unique(cell_type),
    pct_batch1 = round(sum(batch == "1")/n_cells, 3),
    pct_batch2 = round(sum(batch == "0")/n_cells, 3),
    pct_age10 = round(sum(age == "P10")/n_cells, 3),
    pct_age20 = round(sum(age == "P20")/n_cells, 3),
    pct_age60 = round(sum(age == "P60")/n_cells, 3),
    num_genes = round(mean(num_genes_expressed),2),
    mt_ratio = round(mean(mito_ratio), 2),
    umi = round(mean(total_reads),2), 
    umi_sem = round(mean(total_reads)/sqrt(n()), 2)
  )

write.csv(lmmp_dat, here("results/supp_figures/lmmp_metadata_by_cluster.csv"), quote = F, row.names = F)

###Average usage of 6month patterns

cell_weights <- read.csv(here("results/6month_projection_weights/tc_in_6mo_LMMP_old_patterns.csv"))
gene_weights <- read.csv(here("results/6month_projection_weights/pattern_gene_weights.csv"), row.names = 1)

colnames(cell_weights) <- c("cell_id", paste0("sixmo_cellPattern", 1:50))
n_patt <- 50

pData(lmmp) <- pData(lmmp) %>% as.data.frame() %>% tibble::rownames_to_column(var = "cell_id") %>%
  left_join(., cell_weights, by = "cell_id") %>%
  tibble::column_to_rownames(var = "cell_id") %>% DataFrame
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

write.csv(patt_stats, here("results/6month_projection_weights/NMF_mean_by_cluster.csv"), quote = F, row.names = F)
write.csv(gene_weights, here("results/6month_projection_weights/NMF_gene_weights.csv"), quote = F, row.names = F)

  



