---
title: "scDRS followup"
output: html_notebook
---
```{r, message=FALSE}
library(monocle3)
library(ggplot2)
library(tidyverse)
library(here)
library(ComplexHeatmap)
library(circlize)
library(glue)
library(ggpointdensity)

if(!dir.exists(here("DRS/plots"))) dir.create(here("DRS/plots"))
```


```{r}
lmmp <- readRDS(here("TC_LMMP.rds"))

```

```{r}
drs_dir <- here("DRS")
traits_gene_sets_fn <- paste0(drs_dir, "/magma_10kb_1000.74_traits.gs")

traits_gene_sets <- read_tsv(traits_gene_sets_fn)

traits <- traits_gene_sets %>% pull("TRAIT")
```

```{r}
#TODO: Why is this returning only 35537 cells intead of 35657 - almost all RBC - is there some filtering step on genes expressed?
scores <- map(traits, function(trt){
  data.table::fread(here(glue("DRS/out/{trt}.full_score.gz")))
}) %>% setNames(traits)

str(scores[[1]][,1:15])

#long format
scores_df <- scores %>%
  bind_rows(.id = "phenotype")

#pData features to include
feats <- c("log10UMI", "cell_type", "sample")
scores_df <- scores_df %>% left_join(.,
                                     pData(lmmp)[,feats] %>%
                                       as.data.frame() %>%
                                       tibble::rownames_to_column(var = "index")
                                     )



```

```{r}
#for pData
scores_df_sub <- map(traits, function(trt){
  df <- data.frame(scores[[trt]][,"norm_score"]) 
  colnames(df) <- trt
  df
  
}) %>%
  bind_cols %>%
  mutate(index = scores[[1]]$index)

pData(lmmp) <- pData(lmmp) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "index") %>%
  left_join(., scores_df_sub, by = "index") %>%
  tibble::column_to_rownames(var = "index") %>%
  DataFrame()

plot_cells(lmmp, color_cells_by = traits[[1]])
```
```{r}

ggplot(filter(scores_df, phenotype == traits[[1]]), aes(x = log10UMI, y = norm_score)) + geom_point() + ggpointdensity::geom_pointdensity()
ggplot(filter(scores_df, phenotype == traits[[1]]), aes(x = log10UMI, y = pval)) + geom_point() + ggpointdensity::geom_pointdensity()

ggplot(filter(scores_df, phenotype == traits[[1]]), aes(x = log10UMI, y = nlog10_pval)) + geom_point() + ggpointdensity::geom_pointdensity()



```
```{r}
all_boxplots <- ggplot(scores_df) +
  geom_boxplot(aes(x = cell_type, y = norm_score, color = cell_type)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~phenotype, ncol = 2)

pdf(here("DRS/plots/norm_score_boxplots.pdf"), width = 10, height = 220)
  all_boxplots
dev.off()

```

Heatmap for average normalized DRS by cell type
```{r, out.width = "10%" }
mean_score_ct <- map(traits, function(trt){
  df<- scores_df %>% filter(phenotype == trt) %>% group_by(cell_type) %>% summarize(mean_norm_score = mean(norm_score)) 
  colnames(df) <- c("cell_type", trt)
  
  return(df)
}) %>% reduce(left_join, by = "cell_type") %>% tibble::column_to_rownames(var = "cell_type") %>% as.data.frame()

ComplexHeatmap::Heatmap(t(as.matrix(mean_score_ct)))
```



