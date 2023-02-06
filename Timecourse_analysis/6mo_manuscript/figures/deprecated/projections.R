#Visualize projections of 6mo patterns for our data, teichman atlas, vanderbilt data
#pattern by dataset grid

library(monocle3)
library(ggplot2)
library(ggpubr)
library(here)
library(grid)
library(gridExtra)
library(viridis)
library(tidyverse)

source(here("scripts/accessory_functions/monocle_mods.R"))
source(here("scripts/accessory_functions/pattern_plotting.R"))

pattern_cell_weight_filenames <- here("results/6month_projection_weights/tc_in_6mo_LMMP_old_patterns.csv")
pattern_cell_weights <- read.csv(pattern_cell_weight_filenames,row.names = 1)
colnames(pattern_cell_weights) <- colnames(pattern_cell_weights) %>% str_replace_all(., "X6mo_", "")

pattern_weight_filename <- here("results/6month_projection_weights/NMF_gene_weights.csv")
pattern_gene_weights <- read.csv(pattern_weight_filename, row.names = 1)
n_top_genes_show <- 10

lmmp_filename <- here("tc_lmmp_p20.rds")
lmmp <- readRDS(lmmp_filename)
colnames(pData(lmmp)) <- colnames(pData(lmmp)) %>% str_replace_all(., "sixmo_", "")


npattern <- 50
pattern_colnames <- paste0("cellPattern", 1:npattern)

#declare patterns of interest
MENS_patterns <- c(16,27,32,41)
select_patterns <- c(MENS_patterns)

lmmp_pattern_wt_plots <- lapply(select_patterns,
                                plotCellPatterns, 
                                cds_obj = lmmp,
                                pattern_prefix = "cellPattern",
                                do.clip = c(0.01,0.99),
                                outline_size = 1.5,
                                cell_size = 0.45,
                                outline_color= "black") %>%
  lapply(., function(pl){
    pl <- pl + theme(plot.title = element_text(hjust = 0.5))
  })



#gene_tables <-  lapply(paste0("cellPattern",select_patterns), plotGeneWeightTable,
#                       pattern_gene_weights, n_top_genes_show, "Gene weight")

pdf(here("plots/6mo_manuscript/6mo_pattern_projections_p20.pdf"), width = 24, height = 6)
  do.call(ggarrange, c(lmmp_pattern_wt_plots, ncol = length(select_patterns)))
dev.off()


lmmp_pattern_wt_plots_na <- lapply(lmmp_pattern_wt_plots, function(pl){
  pl2 <- pl +
    theme_figure() + 
    theme(axis.line.x=element_blank(),
          axis.line.y=element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(), 
          axis.text.y = element_blank(), 
          axis.title = element_blank(), 
          legend.title = element_blank(), 
          legend.key.height = unit(1.25, "cm"),
          legend.text = element_text(size = 12)) 
  return(pl2)
})

pdf(here("plots/6mo_manuscript/6mo_pattern_projections_p20_minimal.pdf"), width = 24, height = 6)
  do.call(ggarrange, c(lmmp_pattern_wt_plots_na, ncol = length(select_patterns)))
dev.off()


