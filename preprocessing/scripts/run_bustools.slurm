#!/bin/bash
#SBATCH --account=jared
#SBATCH --job-name=Bustools
#SBATCH --output bustools.log
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --mem=40gb
#SBATCH --mail-type=end
#SBATCH --mail-user=jslosbe1@jhu.edu

##---------------------------------------------------------------------
## SLURM job script to run bustools - correct, sort, capture, count reads
## Used 11 Gb of memory over 25 minutes on last run
##---------------------------------------------------------------------
date

#Set up environment
cd /data/users/jared/ENS/6mo_velocity/preprocessing/scripts

bash bustools.sh




echo "slurm script: we got here"


