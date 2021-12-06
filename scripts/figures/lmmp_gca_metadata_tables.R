##Create tables for info by cluster for lmmp
##Create tables for NMF pattern definitions and usage
##Create tables for info by cluster for gut cell atlas (nfh subset)


library(tidyverse)
library(monocle3)
library(here)


### lmmp metadata----

lmmp <- readRDS(here("6month_LMMP.rds"))

pData(lmmp) <- pData(lmmp) %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "index") %>%
  left_join(., cell_weights, by = "index") %>%
  column_to_rownames(var = "index") %>%
  DataFrame()

lmmp_dat <- pData(lmmp) %>%
  as.data.frame() %>%
  group_by(clusters) %>%
  summarise(
    n_cells = n(),
    annotated_type = unique(cell_type_factor),
    pct_batch1 = round(sum(batch == "TH")/n_cells, 3),
    pct_batch2 = round(sum(batch == "TL")/n_cells, 3),
    num_genes = round(mean(num_genes_expressed),2),
    mt_ratio = round(mean(mito_ratio), 2),
    umi = round(mean(total_reads),2), 
    umi_sem = round(mean(total_reads)/sqrt(n()), 2)
  )

write.csv(lmmp_dat, here("plots/supp_figures/lmmp_metadata_by_cluster.csv"), quote = F, row.names = F)

lmmp_dat_ct <- pData(lmmp) %>%
  as.data.frame() %>%
  filter(cell_type_factor %in% c("MENs","NENs","Neuroglia")) %>%
  group_by(cell_type_factor) %>%
  summarise(
    n_cells = n(),
    annotated_type = unique(cell_type_factor),
    pct_batch1 = round(sum(batch == "TH")/n_cells, 3),
    pct_batch2 = round(sum(batch == "TL")/n_cells, 3),
    num_genes = round(mean(num_genes_expressed),2),
    mt_ratio = round(mean(mito_ratio), 2),
    umi = round(mean(total_reads),2), 
    umi_sem = round(mean(total_reads)/sqrt(n()), 2)
  )

write.csv(lmmp_dat_ct, here("plots/supp_figures/lmmp_metadata_by_celltype.csv"), quote = F, row.names = F)



### lmmp NMF pattern usage by cluster and definitions

cell_weights <- read.csv(here("results/NMF/lmmp/old_pattern_run/50dims/pattern_cell_weights.csv"), row.names = 1) %>% tibble::rownames_to_column(var = "index")
gene_weights <- read.csv(here("results/NMF/lmmp/old_pattern_run/50dims/pattern_gene_weights.csv"), row.names = 1)

n_patt <- 50
#mean and std by cluster 
patt_stats <- pData(lmmp) %>%
  as.data.frame() %>%
  group_by(clusters) %>%
  summarise(across(starts_with("cellPattern"), .fns = list(mean = mean, std_dev = sd), .names = "{.col}.{.fn}"))

#add annots
patt_stats <- inner_join(lmmp_dat[,c("clusters","annotated_type", "n_cells")], patt_stats )

#Add ensembl gene_ids 
gene_df <- fData(lmmp) %>% as.data.frame() %>% 
  select(gene_id_trimmed, gene_short_name)
  
gene_weights <- gene_weights %>%
  tibble::rownames_to_column(var = "gene_short_name") %>%
  left_join(gene_df, .)

write.csv(patt_stats, here("plots/supp_figures/NMF_mean_by_cluster.csv"), quote = F, row.names = F)
write.csv(gene_weights, here("plots/supp_figures/NMF_gene_weights.csv"), quote = F, row.names = F)

  


### Gut cell atlas ----
nfh <- readRDS(here("teichman_gut_atlas/data/mesenchymal_gut_nfh_annotated.rds"))

pData(nfh)$total_reads <- Matrix::colSums(exprs(nfh))

cluster_dat <- pData(nfh) %>%
  as.data.frame() %>%
  mutate(Age_group = str_replace(Age_group, "Adult_MLN","Adult")) %>%
  group_by(clusters) %>%
  summarise(
    n_cells = n(),
    pct_adult = round(sum(Age_group == "Adult")/n_cells, 3),
    pct_pediatric = round(sum(Age_group == "Pediatric")/n_cells, 3),
    mean_pattern16 = round(mean(cellPattern16), 3),
    mean_pattern27 = round(mean(cellPattern27), 3),
    mean_pattern32 = round(mean(cellPattern32),3),
    mean_pattern41 = round(mean(cellPattern41),3),
    num_genes = round(mean(n_genes),2),
    umi = round(mean(total_reads),2)
  )

write.csv(cluster_dat, file = here("plots/supp_figures/teichman_nfh_cluster_metadata.csv"))

#prep for grouping adult vs pediatric

age_df <- pData(nfh) %>%
  as.data.frame() %>%
  mutate(Age_group = factor(case_when(
    Age %in% c(4,6,9,10,12) ~ "Juvenile",
    Age %in% c("20-25","25-30","45-50") ~ "Adult",
    Age %in% c("60-65","65-70","70-75") ~ "Aged Adult"),
    levels = c("Juvenile","Adult","Aged Adult")))

age_counts <-  age_df %>% 
  group_by(Age_group) %>%
  summarise(n_samp = n())

map(levels(age_df$Age_group), function(ag){
  age_df %>% filter(Age_group ==  ag) %>% mutate(sample.name = droplevels(sample.name)) %>% pull(sample.name) %>% nlevels()
}) %>% setNames(levels(age_df$Age_group))
