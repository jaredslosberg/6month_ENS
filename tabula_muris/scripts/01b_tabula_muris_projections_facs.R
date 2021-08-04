library(monocle3)
library(SingleCellExperiment)
library(purrr)
library(dplyr)
library(here)
library(ggplot2)
library(projectR)
library(stringr)
library(ggpubr)
library(biomaRt)

#These have pattern plotting on umap with clipped colors, correlation with features and pattern usage
source(here("scripts/accessory_functions/pattern_plotting.R"))
source(here("scripts/accessory_functions/patternFeatureCorrelationHeatmap.R"))

#establish directories for outputs
if(!dir.exists(here("./tabula_muris/plots"))){dir.create(here("./tabula_muris/plots"))}
if(!dir.exists(here("./tabula_muris/data"))){dir.create(here("./tabula_muris/data"))}

sample <- "smartseq2"

#define location to find the processed data
tabula_muris_cds_filepaths <- paste0("/data/users/jared/atlas_processing/tabula_muris/data/",sample,"_tabula_muris.rds")

lmmp <- readRDS(here("6month_LMMP.rds"))

gene_weights <- read.csv(here("results/NMF/lmmp/lmmp_6mo_11-1_pattern_gene_weights.csv"), row.names = 1)
rownames(gene_weights) <- rownames(gene_weights) %>%
  str_split_fixed(., "\\.",2) %>%
  .[,1]
n_patterns <- dim(gene_weights)[2]

cell_weights <- read.csv(here("results/NMF/lmmp/lmmp_6mo_11-1_pattern_cell_weights.csv"), row.names = 1)
    
path <- paste0("/data/users/jared/atlas_processing/tabula_muris/data/",sample,"_tabula_muris.rds")
tabula_muris <- readRDS(path)

#Check to see if there is a relationship between tx length and counts
fData(tabula_muris_filt)$total_reads <- Matrix::rowSums(counts(tabula_muris_filt))
pl <- qplot(x = fData(tabula_muris_filt)$total_reads, y = fData(tabula_muris_filt)$mean_tx_length) +
  ylim(c(0,25000)) + scale_x_log10() + scale_y_log10()

top_hist <- ggplot()+ geom_histogram(aes(fData(tabula_muris_filt)$total_reads)) + scale_x_log10()

right_hist <- ggplot()+ geom_histogram(aes(fData(tabula_muris_filt)$mean_tx_length)) + scale_x_log10() + coord_flip()

empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
        panel.background=element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),           
        axis.title.x=element_blank(), axis.title.y=element_blank())

pl_arranged<- gridExtra::arrangeGrob(top_hist, empty, pl, right_hist,  ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))

png(here("./tabula_muris/plots/smartseq2_tx_length_vs_counts.png"), width = 750, height = 600)
grid.draw(pl_arranged)
dev.off()

#get gene lengths to calculate TPM
ensembl <- useMart("ensembl", host = "http://jul2019.archive.ensembl.org/") #mm 97
ensembl <- useDataset("mmusculus_gene_ensembl",mart=ensembl)

gene_description <- getBM(attributes = c("ensembl_gene_id","external_gene_name","description","chromosome_name", "transcript_length","cds_length"),
                          filters = "ensembl_gene_id", values = fData(tabula_muris)$ensembl_id_trimmed,  mart = ensembl)
  
gene_lengths <- gene_description %>%
  group_by(ensembl_gene_id) %>% 
  summarise(mean_tx_length = mean(transcript_length)) %>%
  rename(ensembl_id_trimmed =ensembl_gene_id)

fData(tabula_muris) <- fData(tabula_muris) %>%
  as.data.frame() %>%
  left_join(., gene_lengths) %>%
  DataFrame()

## This shot out error messages.
# tabula_muris_filt <- tabula_muris[!is.na(fData(tabula_muris)$mean_tx_length),]
# tabula_muris_filt <- scuttle::calculateTPM(tabula_muris_filt, lengths = fData(tabula_muris)$mean_tx_length)

exprs_mat <- normalized_counts(tabula_muris, norm_method = "log") %>% as.matrix()
rownames(exprs_mat) <- fData(tabula_muris)$ensembl_id_trimmed
  
##do the projection
transferred_cell_wts <- t(projectR::projectR(data = exprs_mat,
                                               loadings = as.matrix(gene_weights)))

colnames(transferred_cell_wts) = paste0("cellPattern",1:n_patterns)
pData(tabula_muris) <- cbind(pData(tabula_muris), transferred_cell_wts)
  
##visualize projected pattern usage
##Color UMAP embedding by cell weights for each pattern.
##Call function and return a list of ggplot objects, and plot to one page
weighted_emb <- lapply(1:n_patterns,
                       FUN = plotCellPatterns,
                         cds_obj = tabula_muris,
                         red_method = "UMAP",
                         do.clip = c(0.02,.98))
png(here(paste0("tabula_muris/plots/",sample,"_tabula_muris_in_lmmp_6mo_11_1_patterns.png")), height= 2000, width = 2000)
   print(do.call(ggarrange, weighted_emb))
dev.off()
  
##pattern feature correlations visualized via heatmap
features <- c("sample","cell_ontology_class")
#,"plate.barcode","cell_ontology_class","mouse.sex","mouse.id", "umap.clusters", "num_genes_expressed")
  
heatmaps <- patternFeatureCorrelationHeatmap(tabula_muris, cellWeights.df = transferred_cell_wts, features = features)
  
pdf(here(paste0("tabula_muris/plots/",sample,"_tabula_muris_lmmp_6mo_11_1_patterns_correlations.pdf")), width = 18, height = 10)
    
   heatmaps
dev.off()
  
  
colnames(transferred_cell_wts) = paste0("lmmp_6mo_11_1_pattern",1:n_patterns)
write.csv(transferred_cell_wts, here(paste0("tabula_muris/data/",sample,"_tabula_muris_in_lmmp_6mo_11_1_patterns.csv")))
  
pData(tabula_muris) <- pData(tabula_muris) %>% as.data.frame() %>%
    rename_at(vars(matches('^cellPattern')), ~ str_replace(., 'cellPattern', 'lmmp_6mo_11_1_pattern')) %>% DataFrame()
  
saveRDS(tabula_muris, here(paste0("tabula_muris/data/",sample,"_tabula_muris_annotated.rds")))

##Now look at usages of patterns over annotations such as tissue of origin and annotated cell type

pattern_prefix <- "lmmp_6mo_11_1_pattern"
n_patterns <- 50

usages_annotations <- pData(tabula_muris) %>% as.data.frame()
pl_list <- map(1:n_patterns, function(pattern_no){
  
  pl1 <- ggplot(usages_annotations, aes_string(x= "tissue",
                                               y = paste0(pattern_prefix,pattern_no),
                                               color = "tissue")) +
    geom_violin(position = position_dodge(1)) + 
    geom_boxplot(width = 0.1, position = position_dodge(1)) +
    
    theme(axis.text.x = element_text(angle = 60, hjust=1)) + 
    theme(legend.position = "none")
  
  pl2 <- ggplot(usages_annotations, aes_string(x= "cell_ontology_class",
                                               y = paste0(pattern_prefix,pattern_no),
                                               color = "cell_ontology_class")) +
    geom_boxplot(width = 1, position = position_dodge(1)) +
    
    theme(axis.text.x = element_text(angle = 60, hjust=1)) + 
    theme(legend.position = "none") + 
    coord_flip()
  
  ggpubr::ggarrange(pl1, pl2, ncol = 1, heights = c(1.25,4))
  
})

pdf(here(paste0("tabula_muris/plots/",sample,"_projected_pattern_usages_by_class.pdf")), width = 10, height = 15)
  pl_list
dev.off()

