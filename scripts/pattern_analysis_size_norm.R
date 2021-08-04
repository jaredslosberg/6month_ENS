#' ---
#' title: "6mo LMMP pattern identification"
#' author: "Jared Slosberg"
#' date: "Nov 1, 2020"
#' ---


library(monocle3)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ComplexHeatmap)
library(htmlwidgets)
library(NNLM)
source("./accessory_functions/pattern_plotting.R")
source("./accessory_functions/monocle_mods.R")
source("./accessory_functions/geneSetEnrichment.R")


lmmp_6mo <- readRDS("../6month_LMMP.rds")
set.seed(42)

lmmp_6mo <- estimate_size_factors(lmmp_6mo)
print(dim(lmmp_6mo))

#use either (size -> log) normalized counts
#normed_exprs <- normalized_counts(lmmp_6mo, norm_method = "log")
#or use just log normalized
normed_exprs <- normalized_counts(lmmp_6mo)

#Choose dimensionality of NMF, how many patterns will be identified
npattern <- 60
system.time(decomp <- nnmf(A=as.matrix(normed_exprs),k = npattern,verbose = F, show.warning = T))
rownames(decomp$H) = paste0("cellPattern",c(1:npattern)) 

#Change this to something like cluster, cell type... whatever cell metadata to group by
group_feature <- "clusters"

#Subset to cell name and feature, group of feature to establish column rows
cells_by_group <- as.integer(pData(lmmp_6mo)[,group_feature])
names(cells_by_group) <- colnames(lmmp_6mo)
#order cells by their cluster, this object is a named integer array
cells_by_group <- cells_by_group[order(cells_by_group)]
cells_order <- names(cells_by_group)


#Plot heatmap of H (cell weight) matrix
#Assign a color for each group, these match default ggplot colors
my_colors <- scales::hue_pal()(length(levels(pData(lmmp_6mo)[,group_feature])))
names(my_colors) <- unique(cells_by_group)
ha_cellgroup <- HeatmapAnnotation(cluster = cells_by_group,
                                  col = list(cluster = my_colors))

#We want cells grouped by cluster, as above. We also want to do hierarchical clustering
#on the pre-defined cluster level. Info must be in the same order 
cell_mat <- decomp$H[,cells_order]
stopifnot(colnames(cell_mat) == cells_order)


pdf("/home/jared/ENS/6mo_LMMP/plots/NMF/lmmp/size_norm/NMF_lmmp_6mo_4-16_cellweights.pdf", width = 10)
plot_cells(lmmp_6mo, color_cells_by = group_feature, group_label_size = 3)

ComplexHeatmap::Heatmap(cell_mat, name = "cellscore",
                        top_annotation = ha_cellgroup,
                        column_order = cells_order,
                        #cluster_columns = cluster_within_group(cell_mat, groups_to_cluster),
                        show_column_names = F)
#row_order = rownames(decomp$H))
dev.off()

##Gene weights associated with patterns
geneWeights.df <- as.data.frame(decomp$W)
colnames(geneWeights.df)<-paste0("cellPattern",c(1:npattern))
#Merge cell weights with pData for plotting pattern scores on UMAP embedding
cellWeights.df <- base::as.data.frame(t(decomp$H))
colnames(cellWeights.df) = paste0("cellPattern",c(1:npattern))


#This adds pattern cell weights to pData, make sure they align
#left_join() would be more specific than this cbind()
add.to.cds <- 1
if(add.to.cds){
  stopifnot(rownames(pData(lmmp_6mo)) == rownames(cellWeights.df))
  pData(lmmp_6mo) <- cbind(pData(lmmp_6mo),cellWeights.df)
}


#Color UMAP embedding by cell weights for each pattern.
#Call function and return a list of ggplot objects, and plot to one page
weighted_emb <- lapply(1:npattern,
                       FUN = plotCellPatterns,
                       cds_obj = lmmp_6mo,
                       red_method = "UMAP",
                       do.clip = c(0.02,.98))
png("/home/jared/ENS/6mo_LMMP/plots/NMF/lmmp/size_norm/lmmp_6mo_4-16_NMF_Patterns.png", height= 2000, width = 2000)
print(do.call(ggarrange, weighted_emb[1:npattern]))
dev.off()


#Look at pattern usage over cell age
plot_over_time <- 1
if(plot_over_time){
  condition_patterns<- as.list(1:npattern)
  myplots <- lapply(paste0("cellPattern",condition_patterns), 
                    FUN= plotPatternUsageByCondition,
                    cds = lmmp_6mo,
                    bin_by = "sample")
  pdf("/home/jared/ENS/6mo_LMMP/plots/NMF/lmmp/size_norm/lmmp_6mo_4-16_patterns_over_sample.pdf", height= 20, width = 20)
  print(do.call(ggarrange, myplots[1:length(condition_patterns)]))
  dev.off()
}


tmp<-as.data.frame(merge(fData(lmmp_6mo)[,c("gene_id","gene_short_name")],geneWeights.df,by=0))
dt<- DT::datatable(tmp[,c("gene_id","gene_short_name",
                     unlist(lapply(1:npattern,function(i){paste0("cellPattern",i)}))
)])

saveWidget(dt, file = "/home/jared/ENS/6mo_LMMP/results/NMF/lmmp/size_norm/lmmp_6mo_4-16_pattern_gene_weights.html")

write.csv(cellWeights.df, "/home/jared/ENS/6mo_LMMP/results/NMF/lmmp/size_norm/lmmp_6mo_4-16_pattern_cell_weights.csv")

write.csv(geneWeights.df, "/home/jared/ENS/6mo_LMMP/results/NMF/lmmp/size_norm/lmmp_6mo_4-16_pattern_gene_weights.csv")

#Run GSEA 
geneWeights.df <- geneWeights.df %>% tibble::rownames_to_column(var = "gene_id")


geneSetEnrichment(gene_weights = geneWeights.df,
                  n_genes = 1000, n_patterns = 50, file_prefix = "/home/jared/ENS/6mo_LMMP/results/GSEA/lmmp/size_norm/lmmp_6mo_4-16_")


#Run pattern correlations
pData(lmmp_6mo)$cycle_phase[is.na(pData(lmmp_6mo)$cycle_phase)] <- "NA"
features <- c("cell_type","cycle_phase", "sample")

heatmaps <- patternFeatureCorrelationHeatmap(lmmp_6mo, cellWeights.df = cellWeights.df, features = features)

pdf("/home/jared/ENS/6mo_LMMP/plots/NMF/lmmp/size_norm/LMMP_4-16_pattern_categorical_feature_correlation.pdf",
    width = 18, height = 10)
heatmaps
dev.off()
