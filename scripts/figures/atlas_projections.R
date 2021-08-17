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
vanderbilt_ens_annotated_filename <- here("vanderbilt/data/cds.vanderbilt.10x_inDrop_annotated.rds")

teichman_mesenchymal_atlas <- readRDS(teichman_mesenchymal_atlas_annotated_filename)
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
                                legend_title = "Pattern usage")

teichman_pattern_wt_plots <- lapply(select_patterns,
                                    plotCellPatterns,
                                    cds_obj = teichman_mesenchymal_atlas,
                                    do.clip = c(0.01,0.99), 
                                    legend_title = "")

vanderbilt_pattern_wt_plots <- lapply(select_patterns,
                                      plotCellPatterns,
                                      cds_obj = vanderbilt_ens,
                                      do.clip = c(0.01,0.99), 
                                      legend_title = "")

pattern_usage_plots <- lapply(paste0("cellPattern",select_patterns),
                              plotPatternUsageByCondition,
                              teichman_mesenchymal_atlas,
                              bin_by = "Age_binned",
                              fill_by_append = "_mean") 


pattern_usage_plots <- lapply(pattern_usage_plots, function(x){
  x + theme(axis.text.x = element_text(angle = 90)) + 
    scale_fill_viridis(guide_legend(title="Mean group usage"))
})

gene_tables <-  lapply(paste0("cellPattern",select_patterns), plotGeneWeightTable, pattern_gene_weights, n_top_genes_show, "Gene weight")

png(here("plots/supp_figures/atlas_projections.png"), width = 1800, height = 1300)
do.call(ggarrange,
        c(lmmp_pattern_wt_plots, 
          vanderbilt_pattern_wt_plots,
          teichman_pattern_wt_plots,
          pattern_usage_plots, 
          gene_tables, 
          ncol = 4,
          nrow = 5))
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
