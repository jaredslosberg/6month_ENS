#A function to do gene set enrichment analysis on the amplitude matrix from NMF

# if(!require(clusterProfiler)){
#   BiocManager::install("clusterProfiler")
# }

library(clusterProfiler)
library(dplyr)
library(org.Mm.eg.db)
library(DOSE)
library(ggplot2)
library(grid)
library(gridExtra)

#Helper function- given a matrix of genes by patterns and a specified pattern, return the decreasing ordered genes
#Such that the first gene in genelist is the most highly weighted
sort_pattern_weights <- function(gw, pattern_col){
  sorted_weights <- gw[,c("gene_id", pattern_col)]
  sorted_weights <- sorted_weights[order(sorted_weights[,pattern_col], decreasing = T),]
  genelist <- as.numeric(sorted_weights[,pattern_col])
  #Trim of transcript ID to get trimmed ensembl id
  names(genelist) <- stringr::str_split_fixed(sorted_weights[,"gene_id"],"\\.",2)[,1]
  return(genelist)
}

#TODO: print top weighted genes for each pattern to make sure patterns are same across comparisons
#TODO: Update the way this is implemented, I've been getting warnings

#gene_weights should be a matrix of genes (row) by patterns (column) with a column titled "gene_id" with ensembl ids 
#ngenes is the number of top weighted genes for each pattern to use for GSEA (1000 has worked ok)
geneSetEnrichment <- function(gene_weights, n_genes, n_patterns = NULL, file_prefix){
  if(is.null(n_patterns)){
    patterns <- dim(gene_weights)[2] - 1
    warning("npatterns not provided, using number of columns minus one")
  }else{
    patterns <- 1:n_patterns
  }
  filename_out <- paste0(file_prefix,"GSEA_patterns_",n_genes,"_genes.pdf")
  pdf(filename_out, width = 12, height = 15)
  
  #For each pattern, run GSEA, print the top terms for CC, BP, MF, then print top 110 terms for BP
  gsea_results <- lapply(1:length(patterns), function(i){
    #pattern should match column names for each dimension
    #MIGHT NEED TO BE CHANGED IF NOT CONSISTENT
    pattern <- paste0("cellPattern",patterns[i])
    #genes is ranked numerical vector with gene ids as names
    genes <- sort_pattern_weights(gene_weights, pattern)
    genes <- genes[1:n_genes]
    
    print("Gene ranks prepped... Running GSEA")
    res <- gseGO(geneList     = genes,
                  OrgDb        = org.Mm.eg.db,
                  ont          = "ALL",
                  keyType      = "ENSEMBL",
                  nPerm        = 1000,
                  minGSSize    = 25,
                  maxGSSize    = 500,
                  pvalueCutoff = 1,
                  verbose      = TRUE)
  
  
  
    print(enrichplot::dotplot(res, showCategory=15, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scales = "free")  + 
      ggtitle(paste0("Subset of GO enriched terms for Pattern ", patterns[i]),
              subtitle = paste0("Enrichment calculated from ", n_genes, " non-zero weights genes")))
    
    #Print top ranked genes for each pattern so that its clear which run of patterns were used
    #Should be names and weights of top-ranked genes for each pattern
    grid.newpage()
    gridExtra::grid.table(as.data.frame(genes[1:20]))
    
    
    summary <- as.data.frame(res) %>%
      dplyr::filter(ONTOLOGY == "BP") %>%
      dplyr::select(ONTOLOGY, Description, p.adjust, qvalues, setSize)
    
    colnames(summary)[colnames(summary) == "Description"] <- paste0("Description - GO Terms for pattern ",patterns[i])
    
   
    sub1 <- summary[1:55,]
    sub2 <- summary[56:110,]
    grid.newpage()
    g1 <- gridExtra::grid.table(sub1,theme=ttheme_minimal(base_size = 8)) 
    grid.newpage()
    print("test")
    g2 <- gridExtra::grid.table(sub2,theme=ttheme_minimal(base_size = 8))
    
    
    
  
    return(res)
  })
  dev.off()
}


###Example of how to run:
if(0){
  gene_weights <- read.csv("/home/jared/ENS/Timecourse_ENS/results/NMF/Neuro/NeuroA_pattern_gene_weights.csv")
  #Depending on how weights were saved, make sure ensembl IDS are in column titled with gene_id
  colnames(gene_weights)[1] <- "gene_id"
  gene_weights <- gene_weights[,1:5]
  
  geneSetEnrichment(gene_weights = gene_weights,
                    n_genes = 1000, n_patterns = 4, file_prefix = "tester")
}

