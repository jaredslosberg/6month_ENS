library(monocle3)
library(dplyr)
library(ggplot2)

#Within Monocle framework, do either LRT or Wald Test
#Provide monocle object, method ("LRT" or "wald"), minimum percent of cells expressing a gene to consider,
#and the linear models to be fit


#Minimum number of cells expressing a gene for the gene to be tested for differential expression
#This will not test for very lowly expressed genes, which could be informative
doDifferentialExpression <- function(cds_obj, method = "LRT", min_pct_cells_exprs, full, reduced, ncores = 1){
  min_pct_exprs <- min_pct_cells_exprs
  min_cells_exprs <- round(dim(cds_obj)[2] * min_pct_exprs)
  
  cds_obj <- detect_genes(cds_obj, min_expr = 0)
  cds_obj <- cds_obj[fData(cds_obj)$num_cells_expressed >= min_cells_exprs,]
  
  print(cds_obj)
  
  gene_fits_full <- fit_models(cds_obj, 
                                model_formula_str = full,
                                expression_family = "negbinomial",
                                cores = ncores)
  if(method == "LRT"){
    
      gene_fits_reduced <- fit_models(cds_obj, model_formula_str = reduced,
                                      expression_family = "negbinomial",
                                      cores = ncores)
      #Likelihood ratio test using a full and reduced model
      res <- compare_models(gene_fits_full, gene_fits_reduced)
      
      # lrt_filt <- lrt %>% filter(q_value < 0.05)
      # print(qplot(lrt_filt$num_cells_expressed))
      # #ensembl ids character vector
      # de_genes <- pull(lrt_filt, gene_id)
  }
  
  #If you want to visualize an embedding just on differentially expressed genes
  do.reembed <- 0
  if(do.reembed){
    #How much variance was explained by PC's before
    p1 <- plot_pc_variance_explained(cds_obj) + ggtitle("PCs and variance on all genes")
    print(p1)
    cds_obj <- estimate_size_factors(cds_obj)
    cds_obj <- preprocess_cds(cds_obj, num_dim = 10, use_genes = de_genes, 
                          verbose=TRUE, method = "PCA")
    #How much variance is explained in the new PCs
    p2 <- plot_pc_variance_explained(cds_obj) + ggtitle("PCs and variance on just differentially expressed genes (LRT)")
    print(p2)
    cds_obj <- reduce_dimension(cds_obj,max_components=2,cores=6,verbose=TRUE, preprocess_method = "PCA")
  }
  
  # plot_cells(cds)
  # plot_cells(cds, color_cells_by = "clusters")
  # plot_cells(cds, color_cells_by = "cycle_phase")
  
  
  
  
  if(method == "Wald"){
  #Wald test
    fit_coefs <- coefficient_table(gene_fits_full)
    
    #fit_coefs <- fit_coefs %>% filter (q_value < 0.05) %>% filter(term != "(Intercept)") 
    #  dplyr::select(gene_short_name, term, q_value, estimate)
    
    res <- fit_coefs
  }
  
  return(res)
    }

  