#!/bin/bash

# find . -name "*jpeg" | parallel -I% --max-args 1 convert % %.png

mkdir -p ebi_out/ plots/

source activate scRNA

cd /data/users/jared/ENS/6mo_LMMP/DRS

#write covariates file, write cds format to h5ad
# R --vanilla -e "rmarkdown::render(input = '74_traits_DRS.Rmd', output_file = '74_traits_DRS.html')"

#path to gene set files
set_path="/mnt/morbo/Data/Users/jslosberg/ebi_gwas/gwas_subsets"


#For each gene set file in specified directory, run compute_score in parallel

find $set_path -name *.tsv | parallel -I% -j 20 --max-args 1 \
python ~/bin/scDRS/compute_score.py \
--h5ad_file ../6month_LMMP.h5ad \
--h5ad_species mouse \
--gs_file % \
--cov_file ./data/lmmp_covar.tsv \
--gs_species human \
--out_folder ./ebi_out \
--flag_filter True \
--flag_raw_count True \
--n_ctrl 1000 \
--flag_return_ctrl_raw_score False \
--flag_return_ctrl_norm_score True \

#plotting, interpretations
# R --vanilla -e "rmarkdown::render(input = '74_traits_DRS_followup.Rmd', output_file = '74_traits_DRS_followup.html')"
