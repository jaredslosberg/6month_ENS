#!/bin/bash
#SBATCH --account=jared
#SBATCH --job-name=DRS
#SBATCH --output 74_traits_DRS_ebi_gwas.log
#SBATCH --time=15:00:00
#SBATCH --ntasks=20
#SBATCH --mem=0

##---------------------------------------------------------------------
## SLURM job script to run scDRS
## 12k cells by 30k genes, 20 jobs with 30 gene sets in parallel timed out after 5 hours
##---------------------------------------------------------------------

cd /data/users/jared/ENS/6mo_LMMP/DRS

bash run_all_DRS_ebi_gwas.sh
