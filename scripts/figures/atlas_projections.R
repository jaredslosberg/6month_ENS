#Visualize projections of 6mo patterns for our data, teichman atlas, vanderbilt data
#pattern by dataset grid

library(monocle3)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(here)
library(grid)
library(gridExtra)
library(purrr)
library(viridis)

source(here("scripts/accessory_functions/monocle_mods.R"))
source(here("scripts/accessory_functions/pattern_plotting.R"))

pattern_cell_weight_filenames <- here("results/NMF/lmmp/old_pattern_run/50dims/pattern_cell_weights.csv")
pattern_cell_weights <- read.csv(pattern_cell_weight_filenames,row.names = 1)
pattern_weight_filename <- here("results/NMF/lmmp/old_pattern_run/50dims/pattern_gene_weights.csv")
pattern_gene_weights <- read.csv(pattern_weight_filename, row.names = 1)
n_top_genes_show <- 10

lmmp_filename <- here("6month_LMMP.rds")
lmmp <- readRDS(lmmp_filename)
assertthat::assert_that(all(rownames(pattern_cell_weights) == rownames(pData(lmmp))))
pData(lmmp) <- cbind(pData(lmmp), pattern_cell_weights)

teichman_mesenchymal_atlas_annotated_filename <- here("teichman_gut_atlas/data/mesenchymal_gut_annotated.Rds")
teichman_mesenchymal_nfh_annotated_filename <- here("teichman_gut_atlas/data/mesenchymal_gut_nfh_annotated.rds")

vanderbilt_ens_annotated_filename <- here("vanderbilt/data/cds.vanderbilt.10x_inDrop_annotated.rds")

teichman_mesenchymal_atlas <- readRDS(teichman_mesenchymal_atlas_annotated_filename)
teichman_mesenchymal_nfh <- readRDS(teichman_mesenchymal_nfh_annotated_filename)

#prep for grouping adult vs pediatric
pData(teichman_mesenchymal_nfh) <- pData(teichman_mesenchymal_nfh) %>%
  as.data.frame() %>%
  mutate(Age_group = factor(case_when(
    Age %in% c(4,6,9,10,12) ~ "Juvenile",
    Age %in% c("20-25","25-30","45-50") ~ "Adult",
    Age %in% c("60-65","65-70","70-75") ~ "Aged Adult"
  ), levels = c("Juvenile","Adult","Aged Adult"))) %>% DataFrame


nfh_age_usages <- pData(teichman_mesenchymal_nfh) %>%
  as.data.frame() %>%
  group_by(Age) %>%
  summarise(across(paste0("cellPattern", select_patterns), list(mod = mean), .names = "mean_{col}"), group = unique(Age_group)) %>%
  ungroup() %>%
  group_by(group) %>%
  mutate(across(paste0("mean_cellPattern", select_patterns), list(mod = mean), .names = "group_{col}"))


vanderbilt_ens <- readRDS(vanderbilt_ens_annotated_filename)


npattern <- 50
pattern_colnames <- paste0("cellPattern", 1:npattern)

#declare patterns of interest
MENS_patterns <- c(16,27,32,41)
select_patterns <- c(MENS_patterns)

lmmp_pattern_wt_plots <- lapply(select_patterns,
                                plotCellPatterns, 
                                cds_obj = lmmp,
                                do.clip = c(0.01,0.99),
                                outline_size = 1.5,
                                cell_size = 0.45,
                                outline_color= "black") %>%
  lapply(., function(pl){
    pl <- pl + theme(plot.title = element_text(hjust = 0.5))
  })

teichman_pattern_wt_plots <- lapply(select_patterns,
                                    plotCellPatterns,
                                    cds_obj = teichman_mesenchymal_atlas,
                                    do.clip = c(0.01,0.99),
                                    outline_size = 1.5,
                                    cell_size = 0.45,
                                    outline_color= "black") %>% 
  lapply(., function(pl){
    pl <- pl +
      ggtitle("")
  })


teichman_nfh_pattern_wt_plots <- lapply(select_patterns,
                                    plotCellPatterns,
                                    cds_obj = teichman_mesenchymal_nfh,
                                    do.clip = c(0.01,0.99),
                                    outline_size = 1.5,
                                    outline_color= "black",
                                    cell_size = 0.45) %>% 
  lapply(., function(pl){
    pl <- pl +
      ggtitle("")
  })

vanderbilt_pattern_wt_plots <- lapply(select_patterns,
                                      plotCellPatterns,
                                      cds_obj = vanderbilt_ens,
                                      do.clip = c(0.01,0.99),
                                      outline_size = 1.5,
                                      cell_size = 0.45,
                                      outline_color= "black") %>% 
  lapply(., function(pl){
    pl <- pl +
      ggtitle("")
  })

pattern_usage_plots <- lapply(paste0("cellPattern",select_patterns),
                              plotPatternUsageByCondition,
                              teichman_mesenchymal_atlas,
                              bin_by = "Age",
                              color_fun = mean) %>%
  lapply(., function(pl){
    pl <- pl + 
      ggtitle("") +
      theme(axis.text.x = element_text(angle = 90)) + 
      scale_fill_viridis(guide_legend(title="Mean group usage"))
  })

pattern_usage_plots_nfh <- lapply(paste0("cellPattern",select_patterns),
                              plotPatternUsageByCondition,
                              teichman_mesenchymal_nfh,
                              bin_by = "Age",
                              color_fun = mean) %>%
  lapply(., function(pl){
    pl <- pl + 
      ggtitle("") +
      theme(axis.text.x = element_text(angle = 90)) + 
      scale_fill_viridis(guide_legend(title="Mean group usage")) +
      theme_minimal()
  })

pattern_usage_plots_nfh_grouped <- lapply(paste0("cellPattern",select_patterns), function(patt){
  ggplot(nfh_age_usages, aes_string(x = "group", y= paste0("mean_",patt))) +
    geom_col(aes_string(x = "group", y = paste0("group_mean_",patt)),fill = "grey80", color = "black", position = "dodge") +
    geom_point() +
    ggtitle("") + 
    ylab(patt) +
    theme(axis.text.x = element_text(angle = 90)) + 
    theme_minimal()
  })




gene_tables <-  lapply(paste0("cellPattern",select_patterns), plotGeneWeightTable, pattern_gene_weights, n_top_genes_show, "Gene weight")

png(here("plots/supp_figures/atlas_projections_mens_patterns.png"), width = 1800, height = 1300)
do.call(ggarrange,
        c(lmmp_pattern_wt_plots, 
          vanderbilt_pattern_wt_plots,
          teichman_pattern_wt_plots,
          pattern_usage_plots,
          #gene_tables, 
          ncol = 4,
          nrow = 4))
dev.off()

png(here("plots/supp_figures/atlas_projections_nfh_mens_patterns.png"), width = 1800, height = 1300)
do.call(ggarrange,
        c(lmmp_pattern_wt_plots, 
          vanderbilt_pattern_wt_plots,
          teichman_nfh_pattern_wt_plots,
          pattern_usage_plots_nfh, 
          nrow = 4,
          ncol = 4))
dev.off()

png(here("plots/supp_figures/atlas_projections_nfh_binned_mens_patterns.png"), width = 1800, height = 1300)
do.call(ggarrange,
        c(lmmp_pattern_wt_plots, 
          vanderbilt_pattern_wt_plots,
          teichman_nfh_pattern_wt_plots,
          pattern_usage_plots_nfh_grouped, 
          nrow = 4,
          ncol = 4))
dev.off()


###STATS for Teichman atlas over age
patterns_to_match <- paste0("cellPattern",select_patterns)

mean_df <- pData(teichman_mesenchymal_atlas) %>%
  as.data.frame () %>%
  group_by(Age_binned) %>% summarize(across(all_of(patterns_to_match), ~ mean(.x), .names = "{.col}_mean"))

n_df <- pData(teichman_mesenchymal_atlas) %>%
  as.data.frame () %>%
  group_by(Age_binned) %>% summarize(across(all_of(patterns_to_match), ~ length(.x), .names = "{.col}_n"))

sem_df <- pData(teichman_mesenchymal_atlas) %>%
  as.data.frame () %>%
  group_by(Age_binned) %>% summarize(across(all_of(patterns_to_match), ~ sd(.x)/length(.x), .names = "{.col}_sem"))

summ_stats <- cbind(mean_df, n_df, sem_df)
write.csv(summ_stats, here("teichman_gut_atlas/data/pattern_usage_stats_over_age.csv"))

write.csv(nfh_age_usages, here("plots/supp_figures/pattern_usage_stats_over_age_nfh.csv"))
