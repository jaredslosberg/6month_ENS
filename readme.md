 # Analysis of 6month LMMP scRNA
Library generated with 10xv2 chemistry
Tissue prepped by longitudinal muscle myenteric plexus (LMMP) from two 6month animals (TH, TL)
 
## Preprocessing to generate counts with kallisto bustools
Scripts in ./preprocessing/
Scripts run in the following order:
get_velo_files_mm10_10xv2.R
  index_build.sh
    kallisto_bus_pseudoalign.sh
      bustools.sh

## Processing from raw UMI counts, QC, dimensionality reduction
./scripts ./data ./results ./plots hold the relevant info
Modular functions held in ./scripts/accessory_functions
Scripts to generate figs and supp figs in .scripts/figures

Scripts run in the following order: 
MES_6mo_LMMP_analysis.Rmd
  pattern_analysis.R
  6mo_MENS_followup.html
  
  
## Analysis of relevant public datasets
tabula_muris/ teichman_gut_atlas/ vanderbilt/
Preprocessing of each datasets elsewhere, e.g. /data/users/jared/atlas_processing
