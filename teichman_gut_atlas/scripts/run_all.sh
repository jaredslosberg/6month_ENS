#!/bin/bash 

#Run scripts for gut cell atlas analysis that is specific to 6mo_lmmp project


#Projections through the patterns learned in the 6month lmmp scRNA
R -e "rmarkdown::render('gut_atlas_projections.Rmd', output_file= 'gut_atlas_projections.html')"

#TODO: Projection_drivers based analysis - unfinished
R -e "rmarkdown::render('gut_atlas_projection_drivers.Rmd', output_file= 'gut_atlas_projection_drivers.html')"

#Visualize patterns on nfh subsection 
R -e "rmarkdown::render('pattern_vis_nfh.Rmd', output_file= 'pattern_vis_nfh.html')"

#TODO: AUC based projection inference


#Genes specific to clusters identified as MENs-like based on pattern usages
R -e "rmarkdown::render('specific_genes.Rmd', output_file= 'specific_genes.html')"