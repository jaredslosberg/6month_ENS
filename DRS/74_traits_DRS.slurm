#!/bin/bash
#SBATCH --account=jared
#SBATCH --job-name=DRS
#SBATCH --output 74_traits_DRS.log
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --mem=50gb


##---------------------------------------------------------------------
## SLURM job script to run scDRS
## No idea how long or how much memory for 35k cells x 30k genes, 74 gene sets
##---------------------------------------------------------------------

source activate scRNA

jupyter nbconvert --clear-output --execute --to notebook 74_traits_DRS.ipynb