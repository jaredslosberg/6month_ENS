library(monocle3)
library(ggplot2)
library(tidyverse)
library(here)
library(ComplexHeatmap)

lmmp <- readRDS(here("6month_LMMP.rds"))
rownames(lmmp) <- fData(lmmp)$gene_short_name

genes <- c("Cdh3")

gr_df <- pData(lmmp) %>% as.data.frame %>% select(cell_type) %>% rownames_to_column(var = "cell_id")

exprs_mat <- assay(lmmp[fData(lmmp)$gene_short_name %in% genes], "logcounts") %>% t() %>% as.matrix() %>% as.data.frame() %>% rownames_to_column(var = "cell_id") %>% 
  inner_join(., gr_df)  %>%
  group_by(cell_type) %>%
  summarise(cdh3 = mean(Cdh3))%>% 
  column_to_rownames(var = "cell_type") %>%
  scale(., center = T, scale = T)

Heatmap(exprs_mat, name = "Z-scored expression")
