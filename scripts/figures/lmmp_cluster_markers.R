## For lmmp data, make table with specific genes for each cluster

library(here)
library(tidyverse)
library(monocle3)

lmmp <- readRDS("6month_LMMP.rds")

tm <- monocle3::top_markers(lmmp, group_cells_by = "cell_type_factor", genes_to_test_per_group = 30, cores = 16 )

col_to_keep <- c("cell_group","gene_id","gene_short_name","mean_expression","fraction_expressing","specificity","marker_test_p_value","marker_test_q_value")

tm_filt <- tm[,col_to_keep]

write.csv(tm_filt, here("plots/supp_figures/lmmp_marker_genes_by_celltype.csv"), quote = F, row.names = F)          