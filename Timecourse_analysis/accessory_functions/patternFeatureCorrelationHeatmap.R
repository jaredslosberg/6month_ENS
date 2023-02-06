#### Check correlation of patterns with cell features
#### One hot encoding method.
### For each pattern, compute correlation between two vectors. One is continuous (cell pattern weight)
### and one is binary (yes/no is that cell in that group)

library(ComplexHeatmap)
library(circlize)

#as implemented, features should be a categorical value, can also change for continuous
#cell weights data frame should have cell names on rownames, columns should be cellPattern{#}
patternFeatureCorrelationHeatmap <- function(cds, cellWeights.df, features, corr_method = "pearson", color_scale = NULL, pattern_prefix = "cellPattern", returnData = F){

  #categorical features, one-hot encoding

  corr_hm <- lapply(features, function(feature){
    feature_levels <- sort(unique(pData(cds)[,feature]))

    npattern <- dim(cellWeights.df)[2]
    corr_mat <- matrix(nrow = length(feature_levels), ncol = npattern)


    #loop over each pattern
    for(i in 1:npattern){
      cell_weights <- cellWeights.df[,paste0(pattern_prefix,i)]
      names(cell_weights) <- rownames(cellWeights.df)
      #re-initialize empty correlation vector
      corr_vec <- numeric()
      #loop over each level within the categorical variable
      for(level in feature_levels){

        values <- pData(cds)[,feature]
        one_hot_values <- (values == level)
        names(one_hot_values) <- colnames(cds)
        #Make sure our cells are in the same order
        assertthat::assert_that(sum(names(one_hot_values) == names(cell_weights)) == length(one_hot_values))

        #calculate correlation between two vectors
        level_cor <- cor(cell_weights, one_hot_values, method = corr_method)
        #corr_df is a one-dimensional vector, eg the correlations between all features for a single pattern
        corr_vec <- c(corr_vec, level_cor)
      }
      #Now we have a n_feature_levels x n_pattern matrix with pairwise correlations
      corr_mat[,i] <- corr_vec
    }
    rownames(corr_mat) <- feature_levels
    colnames(corr_mat) <- paste0(pattern_prefix,1:npattern)
    
    if(!returnData){
      if(is.null(color_scale)){
        #default
       color_scale <- circlize::colorRamp2(c(-1,0,1), c("blue","white","red"))
       }
        
      heatmap_figure <- ComplexHeatmap::Heatmap(corr_mat, name = corr_method,
                                                column_title = paste0("correlation between pattern and ",feature, " annotation"),
                                                col = color_scale)
      return(heatmap_figure)
    }else{
      return(corr_mat)
    }
  }) %>% setNames(features)

}





