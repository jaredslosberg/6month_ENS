##functions to take cell data set and prepare expression and annotations matrices for 
##heatmap visualization with scrattch.vis sparklines plots

###Define functions to prepare matrices for sparklines plots ----
prep_annotation_matrix <- function(cds_obj, facet_by, colors = NULL, labels = NULL){
  
  if(is.null(labels)){
    
    labels <- cds_obj %>%
      pData() %>%
      as.data.frame() %>%
      dplyr::select(all_of(facet_by)) %>%
      pull(facet_by) %>%
      unique() %>% as.vector()
  }
  
  if(is.null(colors)){colors <- scales::hue_pal()(length(labels))}
  
  colorsident <- cbind(ident = labels,
                       colors = colors)
  
  # Create annotation data.frame
  anno.df <- as.data.frame(cbind(
    sample_name = colnames(cds_obj),
    primary_type_id = colorsident[match(as.character(pData(cds_obj)[,facet_by]), colorsident[,1]),"ident"],
    primary_type_label = as.character(pData(cds_obj)[,facet_by]),
    primary_type_color = colorsident[match(as.character(pData(cds_obj)[,facet_by]), colorsident[,1]),"colors"]
  ))
  
  return(anno.df)
}

#exprs_mat_gene_names is the column of pData that is the same gene id as the rownames of expression matrix
#short_gene_names is column of pData to include on the plot
prep_expression_matrix <- function(cds, norm_method, exprs_assay = "counts", exprs_mat_gene_names = "gene_id", short_gene_names = "gene_short_name"){
  gene_name_df <- fData(cds)[,short_gene_names, drop = F] %>%
    as.data.frame() %>% 
    tibble::rownames_to_column(var = exprs_mat_gene_names)
  #cells as rows, genes as columns per documentation
  if(norm_method %in% c("log")){
    exprs_mat <- normalized_counts(cds, norm_method = norm_method)
  }
  
  if(norm_method == "none"){
    message(paste0("Normalization method 'none' selected. Using assay ", exprs_assay))
    exprs_mat <- assay(cds, exprs_assay)
  }
  
  exprs_mat <- exprs_mat %>%
    as.matrix() %>% 
    as.data.frame() %>%
    tibble::rownames_to_column(var = exprs_mat_gene_names) %>%
    inner_join(gene_name_df, by = exprs_mat_gene_names) %>%
    tibble::column_to_rownames(var = short_gene_names) %>%
    dplyr::select(-any_of(exprs_mat_gene_names)) %>%
    t() %>% as.data.frame %>%
    tibble::rownames_to_column(var= "sample_name")
  
  return(exprs_mat)
}

shorten_gene_lists <- function(gene_vec, n){
  gene_list <- list()
  
  #while gene set is too long, subset the first n elements and expand list
  while(length(gene_vec) > n){
    genes <- list(gene_vec[1:n])
    gene_vec <- gene_vec[-c(1:n)]
    gene_list <- c(gene_list, genes)
  }
  
  #If there's any genes left (i.e. length of set was not a multiple of n), add them too
  if(!isEmpty(gene_vec)){
    gene_list <- c(gene_list, list(gene_vec))
  }
  return(gene_list)
}     
