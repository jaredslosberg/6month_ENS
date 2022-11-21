#!/bin/bash
#SBATCH --account=jared
#SBATCH --job-name=Pseudoalign
#SBATCH --output kallisto_bus_pseudoalign.log
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --mem=70gb
#SBATCH --mail-type=end
#SBATCH --mail-user=jslosbe1@jhu.edu

##---------------------------------------------------------------------
## SLURM job script to run kallisto pseudoalign
## p10 (8 fastqs) run utilized 58 gb of memory over 40 minutes
##---------------------------------------------------------------------

date

##Set up environment
source activate kallistobus

cd /data/users/jared/ENS/6mo_LMMP/preprocessing/scripts
bash kallisto_bus_pseudoalign.sh > kallisto_bus_pseudoalign.log 2>&1

echo "slurm scripts: we got here"


