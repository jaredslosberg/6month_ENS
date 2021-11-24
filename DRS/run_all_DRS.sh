#!/bin/bash

cd /data/users/jared/ENS/6mo_LMMP/DRS

mkdir -p out/ plots/

source activate scRNA

#write covariates file, write cds format to h5ad
R --vanilla -e "rmarkdown::render(input = '74_traits_DRS.Rmd', output_file = '74_traits_DRS.html')"

#run scDRS and calculate per cell DRS
python ~/bin/scDRS/compute_score.py \
--h5ad_file ../6month_LMMP.h5ad \
--h5ad_species mouse \
--gs_file ./data/gs_file/magma_10kb_1000.74_traits.gs \
--cov_file ./data/lmmp_covar.tsv \
--gs_species human \
--out_folder ./out \
--flag_filter True \
--flag_raw_count True \
--n_ctrl 1000 \
--flag_return_ctrl_raw_score False \
--flag_return_ctrl_norm_score True \

#plotting, interpretations
# R --vanilla -e "rmarkdown::render(input = '74_traits_DRS_followup.Rmd', output_file = '74_traits_DRS_followup.html')"