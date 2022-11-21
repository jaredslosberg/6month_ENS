cd /data/users/jared/ENS/preprocessing/scripts
#cd ./preprocessing/scripts


#Downloads necessary files for index building
Rscript get_velo_files_mm10_10xv2.R > get_velo_files_mm10_10xv2.log 2>&1

#Builds index
bash index_build.sh > index_build.log 2>&1

#run kallisto pseudoalignment
jobid=$(sbatch --parsable run_kallistobus.slurm)

#when pseudoalignment is done, do count quantification with bustools
sbatch --dependency=afterok:$jobid run_bustools.slurm 


