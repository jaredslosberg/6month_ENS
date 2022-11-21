
library(here)
library(monocle3)
library(tidyverse)

source(here("scripts/accessory_functions/monocle_mods.R"))
source(here("scripts/accessory_functions/pattern_plotting.R"))

lmmp <- readRDS("/data/users/jared/ENS/6mo_LMMP/6month_LMMP.rds")

plot_cells(lmmp)

res <- top_markers(lmmp, group_cells_by = "clusters", genes_to_test_per_group = 50, cores = 12)

pdat <- pData(lmmp)[,c("clusters", "cell_type")] %>% as.data.frame() %>%
  mutate(cell_group = clusters) %>% 
  select(-c("clusters")) %>% unique

res_ct <- res %>% left_join(., pdat)  %>%
  select(c("cell_type" , "cell_group", "gene_id", "gene_short_name","marker_score",
           "mean_expression", "fraction_expressing","specificity" ,"pseudo_R2" ,
           "marker_test_p_value" ,"marker_test_q_value")) %>% 
  arrange(cell_group)

write.csv(res_ct, "/data/users/jared/ENS/6mo_LMMP/results/cluster_markers.csv")
