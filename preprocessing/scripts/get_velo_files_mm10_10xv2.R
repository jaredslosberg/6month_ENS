sessionInfo()
print("R script is running...")




# Install devtools if it's not already installed
#if (!require(devtools)) {
#  install.packages("devtools")
#}
# Install from GitHub
devtools::install_github("BUStools/BUSpaRse")
#devtools::install_github("satijalab/seurat-wrappers")




#if (!require(BiocManager)) {
#  install.packages("BiocManager")
#  BiocManager::install(c("DropletUtils", "BSgenome.Mmusculus.UCSC.mm10",
#                       "AnnotationHub", "SingleR","biomaRt","tidyverse","scales","zeallot","densityClust"))
#}



## Load required libraries
#Once we have installed the necessary packages, we can load them into the R session to begin our analysis.
library(BiocGenerics)
library(AnnotationHub)
library(biomaRt)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BUSpaRse)


# Kallistobus Index Build

#This notebook demonstrates the use of command line tools `kallisto` and `bustools`. Please use `kallisto` >= 0.46, whose binary can be downloaded [here](https://pachterlab.github.io/kallisto/download). Also, please use `bustools` >= 0.39.3, whose binary of `bustools` can be found [here](https://github.com/BUStools/bustools/releases). User interface of `bustools` has changed in version 0.39.3. For version 0.39.2, see earlier git commits of this notebook.

#after you download the binary, you should decompress the file (if it is `tar.gz`) with `tar -xzvf file.tar.gz` in the `bash` terminal, and add the directory containing the binary to `PATH` by `export PATH=$PATH:/foo/bar`, where `/foo/bar` is the directory of interest. Then you can directly invoke the binary on the command line as we will do in this notebook.

## Get gene annotations to prepare a reference index
#In order to map the read sequences to their cognate genes/transcripts, we first must download and prepare a reference transcriptome. Here we are using the mouse reference transcriptome from [Ensembl version 97](http://jul2019.archive.ensembl.org/Mus_musculus/Info/Index).  We download this dataset using a standard Bioconductor 'AnnotationHub' workflow.


# query AnnotationHub for mouse Ensembl annotation
ah <- AnnotationHub() # Connect to the Annotation Hub to query.
record<-query(ah,pattern = c("Ensembl", "97", "Mus musculus", "EnsDb")) # Identify the AnnotationHub record associated with Mouse Ensembl v97


# Once we've identified the AnnotationHub record id `r names(record)` for Mouse Ensembl v97, we can use this to retrieve all of the annotation records for this build.
# Get mouse Ensembl 97 annotation
edb <- ah[[names(record)]] # I need to comment this out to get the material to build currently...Please uncomment if you are running yourself.


## Get and Organize Required Preprocessing Files
#As part of this exercise, we will be performing an RNA Velocity analysis (detailed below).  This requires us to generate two seperate count matrices, one for reads mapping to 'spliced' transcripts and one for reads mapping to 'unspliced' transcripts. To do this dual mapping, we need to extract the information that we need to separately map the reads to cDNA sequences as well as intronic sequences. The `BUSpaRse` library has a function `get_velocity_files()` that will extract and process the necessary information from an `AnnotationHub` object.

#_*Note:* While this approach generates both the spliced and unspliced annotation records, for most of the downstream analysis we will just be using the count matrices associated with the spliced read mappings._

#L=98 for 10xv2 
get_velocity_files(edb, L = 98, Genome = BSgenome.Mmusculus.UCSC.mm10,out_path = "/data/users/jared/ENS/6mo_velocity/bustools_required_files_10xv2_mm/", isoform_action = "separate")


## Build Reference index
#Once we have the files for the cDNA sequences and the intronic sequences, we need to construct a reference kallisto index for the mapping.  This indexing scans through the reference transcriptome and organizes it in a way that makes mapping the read sequences faster and more efficient.
system('wget https://teichlab.github.io/scg_lib_structs/data/737K-august-2016.txt -O /data/users/jared/ENS/6mo_LMMP/preprocessing/bustools_required_files_10xv2_mm/737K-august-2016.txt')
