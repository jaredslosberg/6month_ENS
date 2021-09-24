#!/bin/bash
mkdir -p data/ plots/ 

#setwd to project home

#Run kallisto bustools workflow
#Requires building kallistobus.yml
bash ./preprocessing/scripts/run_all.sh

R -e "rmarkdown::render('scripts/MES_6mo_LMMP_analysis.Rmd')"

#depends on MES_6mo_LMMP_analysis.Rmd
R -e "rmarkdown::render('scripts/tricycle_analysis.Rmd')"

#depends on MES_6mo_LMMP_analysis.Rmd
R -e "rmarkdown::render('scripts/pattern_analysis.Rmd')"

#depends on pattern_analysis.Rmd
R -e "rmarkdown::render('scripts/NMF_UMI_correlation.Rmd')"

#depends on MES_6mo_LMMP_analysis.Rmd
R -e "rmarkdown::render('scripts/6mo_MENS_followup.Rmd')"

#Velocity analysis

#depends on MES_6mo_LMMP_analysis.Rmd
R -e "rmarkdown::render('scripts/6mo_MENS_followup.Rmd')"

#Requires building scvelocity.yml
conda activate scVelocity
jupyter nbconvert --to notebook --execute ./velocity/scripts/velocity.ipynb