[![DOI](https://zenodo.org/badge/392759401.svg)](https://zenodo.org/doi/10.5281/zenodo.10278287)

Data analysis repository for Kulkarni et al, [Age-associated changes in lineage composition of the enteric nervous system regulate gut health and disease](https://doi.org/10.7554/eLife.88051.1)

BioRxiv version of our [manuscript](https://www.biorxiv.org/content/10.1101/2020.08.25.262832v1.full)

 # Analysis of 6month LMMP scRNA

Raw FASTQs, count matrices, and metadata are available at GSE156146. Cell metadata available in cell_clusters_and_metadata.csv

Library generated with 10xv2 chemistry
Tissue prepped from longitudinal muscle myenteric plexus (LMMP) from two 6 month mice
 
## Preprocessing to generate counts with kallisto bustools
Scripts in ./preprocessing/
Scripts run in the following order:
get_velo_files_mm10_10xv2.R
  index_build.sh
    kallisto_bus_pseudoalign.sh
      bustools.sh
      
./preprocessing/run_all.sh will run the 4 scripts above. Path to fastqs 
  (in kallisto_bus_pseudoalign.sh) must be correct

## Processing from raw UMI counts, QC, dimensionality reduction
./scripts ./data ./results ./plots hold the relevant info
Modular functions held in ./scripts/accessory_functions

Scripts to generate figures and supplemental figures are available in ./scripts/figures

## Analysis of public datasets
Elmentaite et al: Spliced and unspliced count matrices for mesenchymal subset of the Gut Cell Atlas were obtained from https://www.gutcellatlas.org/

May-Zhang et al: Adult ileal snRNAseq counts matrices et al were obtained from GEO GSE153192


