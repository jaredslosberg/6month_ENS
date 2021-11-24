#!/bin/bash
mkdir -p data/ plots/NMF/ results/ results/GSEA/lmmp results/GSEA/mens results/NMF/MENS results/NMF/lmmp results/MENS

#setwd to project home

#Run kallisto bustools workflow
#Requires building kallistobus.yml
#bash ./preprocessing/scripts/run_all.sh

R -e "rmarkdown::render('scripts/MES_6mo_LMMP_analysis.Rmd', output_format = 'html_document')"

#depends on MES_6mo_LMMP_analysis.Rmd
R -e "rmarkdown::render('scripts/tricycle_analysis.Rmd', output_format = 'html_document')"

#depends on MES_6mo_LMMP_analysis.Rmd
R -e "rmarkdown::render('scripts/pattern_analysis.Rmd', output_format = 'html_document')"

#depends on MES_6mo_LMMP_analysis.Rmd
R -e "rmarkdown::render('scripts/6mo_MENS_followup.Rmd', output_format = 'html_document')"

#Velocity analysis

#Requires building scvelocity.yml
conda activate scVelocity
jupyter nbconvert --clear-output --to notebook --execute ./velocity/scripts/velocity.ipynb


#Figures and supplemental figures
# Rscript --vanilla cluster_expression_violin.R
# Rscript --vanilla lmmp_NMF_heatmap.R
# Rscript --vanilla neuronal_markers_sparklines.R
# Rscript --vanilla suppl_figure_sparklines.R
# Rscript --vanilla atlas_projections.R

#TODO: run external datasets e.g. teichman gut atlas
