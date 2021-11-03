#!/bin/bash
source activate scRNA

cd /data/users/jared/ENS/Timecourse_ENS/DRS

python ~/bin/scDRS/compute_score.py \
--h5ad_file ../TC_LMMP.h5ad \
--h5ad_species mouse \
--gs_file magma_10kb_1000.74_traits.gs \
--cov_file lmmp_covar.tsv \
--gs_species human \
--out_folder out \
--flag_filter True \
--flag_raw_count True \
--n_ctrl 1000 \
--flag_return_ctrl_raw_score False \
--flag_return_ctrl_norm_score True \
