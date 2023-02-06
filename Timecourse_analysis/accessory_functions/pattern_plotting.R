library(monocle3)
library(viridis)
library(ggplot2)
library(here)
source("/data/users/jared/ENS/Timecourse_ENS/scripts/accessory_functions/monocle_mods.R")

#Define function to create a plot for a given pattern
plotCellPatterns = function (pattern_number, cds_obj, red_method = "UMAP", do.clip = NULL, pattern_prefix = "cellPattern", title_prefix = NULL, ...) {
  
  if(is.null(title_prefix)) title_prefix <- pattern_prefix
  
  if(is.null(do.clip)){
    plot_cells_mod(cds_obj, 
                   reduction_method = red_method, 
                   color_cells_by = paste0(pattern_prefix,pattern_number),
                   ...) 
    #  + ggtitle(paste0("Pattern ",pattern_number)) + 
    #    font("title", size= 40)
  }else{   
    plot_cells_mod(cds_obj,
                   reduction_method = red_method,
                   color_cells_by = paste0(pattern_prefix,pattern_number),
                   ...) +
      scale_color_viridis(limits = quantile(pData(cds_obj)[,paste0(pattern_prefix,pattern_number)], do.clip) , oob = scales::squish,
                          guide_legend(title = paste0(pattern_prefix,pattern_number))) +
      ggtitle(paste0(title_prefix, pattern_number))
  }
}

###Example usage - taken out of context might not run as is
#print to png to minimize file size. clip so that outlier cells do not dominate color scaling
if(0){
weighted_emb <- lapply(1:npattern,
                       FUN = plot_cell_patterns,
                       cds_obj = MENS,
                       red_method <- "UMAP",
                       do.clip <- c(0.02,.98))
png("/home/jared/ENS/Timecourse_ENS/plots/NMF/MENS/MENS_A_NMF_Patterns.png", height= 2000, width = 2000)
print(do.call(ggarrange, weighted_emb[1:length(condition_patterns)]))
dev.off()
}

#Function to plot a continuous variable (eg pattern usage) over a binned variable such as time.
#provide feature (continuous) and bin_by/color_by (categorical, discrete) as a string that matches a column in pData.
plotPatternUsageByCondition <- function(feature, cds, conditions,  bin_by, color_fun=NULL){
  
  pdat <- pData(cds) %>%
    as.data.frame()
  
  if(!is.null(color_fun)){
    #if function provided, operate over feature column
    pdat <- pdat %>%
      group_by(across(bin_by)) %>%
      mutate(across(feature, list(mod=color_fun), .names = "mod_{col}"))
  }
  
  #TODO: assert that each group has one unique value
  
  #plot, and color by modified feature value if provided
  if(!is.null(color_fun)){
    ggplot(pdat, aes_string(x = bin_by, y = feature)) +
      geom_boxplot(aes_string(fill = paste0("mod_", feature))) +
      ggtitle(feature) +
      ylab("cell pattern weight") + 
      xlab(bin_by) + 
      scale_fill_viridis()
  } else{
    ggplot(pdat, aes_string(x = bin_by, y = feature)) +
      geom_boxplot() +
      ggtitle(feature) +
      ylab("cell pattern weight") + 
      xlab(bin_by)
  }
}

###Example usage - taken out of context might not run as is
if(0){
  condition_patterns<- as.list(1:npattern)
  myplots <- lapply(paste0(pattern_prefix,condition_patterns), 
                    FUN= plotPatternUsageByCondition,
                    cds = cds,
                    bin_by = "age")
  pdf("/home/jared/ENS/Timecourse_ENS/plots/NMF/MENS/MENS_patterns_over_age.pdf", height= 20, width = 20)
  print(do.call(ggarrange, myplots[1:length(condition_patterns)]))
  dev.off()
}