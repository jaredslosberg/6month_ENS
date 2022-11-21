#!/bin/bash
source activate kallistobus

cd /data/users/jared/ENS/6mo_LMMP/preprocessing
mkdir -p bus_output
cd bus_output

echo "Processing fastqs from TH"
kallisto bus -i ../bustools_required_files_10xv2_mm/mm_cDNA_introns_mm_10xv2.idx -o ./TH \
-x 10xv2 -t 16 \
/mnt/morbo/Data/Sequencing/Raw_Backup/Kulkarni/Kulkarni_10X/fastqs/Kulkarni_LMMP_SK-TH/Kulkarni_SK-TH_MissingLibrary_1_HKJ5JBCXY/SK-TH_S1_L001_R1_001.fastq.gz \
/mnt/morbo/Data/Sequencing/Raw_Backup/Kulkarni/Kulkarni_10X/fastqs/Kulkarni_LMMP_SK-TH/Kulkarni_SK-TH_MissingLibrary_1_HKJ5JBCXY/SK-TH_S1_L001_R2_001.fastq.gz \
/mnt/morbo/Data/Sequencing/Raw_Backup/Kulkarni/Kulkarni_10X/fastqs/Kulkarni_LMMP_SK-TH/Kulkarni_SK-TH_MissingLibrary_1_HKJ5JBCXY/SK-TH_S1_L001_R1_002.fastq.gz \
/mnt/morbo/Data/Sequencing/Raw_Backup/Kulkarni/Kulkarni_10X/fastqs/Kulkarni_LMMP_SK-TH/Kulkarni_SK-TH_MissingLibrary_1_HKJ5JBCXY/SK-TH_S1_L001_R2_002.fastq.gz \
/mnt/morbo/Data/Sequencing/Raw_Backup/Kulkarni/Kulkarni_10X/fastqs/Kulkarni_LMMP_SK-TH/Kulkarni_SK-TH_MissingLibrary_1_HKJ5JBCXY/SK-TH_S1_L001_R1_003.fastq.gz \
/mnt/morbo/Data/Sequencing/Raw_Backup/Kulkarni/Kulkarni_10X/fastqs/Kulkarni_LMMP_SK-TH/Kulkarni_SK-TH_MissingLibrary_1_HKJ5JBCXY/SK-TH_S1_L001_R2_003.fastq.gz \
/mnt/morbo/Data/Sequencing/Raw_Backup/Kulkarni/Kulkarni_10X/fastqs/Kulkarni_LMMP_SK-TH/Kulkarni_SK-TH_MissingLibrary_1_HKJ5JBCXY/SK-TH_S1_L001_R1_004.fastq.gz \
/mnt/morbo/Data/Sequencing/Raw_Backup/Kulkarni/Kulkarni_10X/fastqs/Kulkarni_LMMP_SK-TH/Kulkarni_SK-TH_MissingLibrary_1_HKJ5JBCXY/SK-TH_S1_L001_R2_004.fastq.gz \
/mnt/morbo/Data/Sequencing/Raw_Backup/Kulkarni/Kulkarni_10X/fastqs/Kulkarni_LMMP_SK-TH/Kulkarni_SK-TH_MissingLibrary_1_HKJ5JBCXY/SK-TH_S1_L002_R1_001.fastq.gz \
/mnt/morbo/Data/Sequencing/Raw_Backup/Kulkarni/Kulkarni_10X/fastqs/Kulkarni_LMMP_SK-TH/Kulkarni_SK-TH_MissingLibrary_1_HKJ5JBCXY/SK-TH_S1_L002_R2_001.fastq.gz \
/mnt/morbo/Data/Sequencing/Raw_Backup/Kulkarni/Kulkarni_10X/fastqs/Kulkarni_LMMP_SK-TH/Kulkarni_SK-TH_MissingLibrary_1_HKJ5JBCXY/SK-TH_S1_L002_R1_002.fastq.gz \
/mnt/morbo/Data/Sequencing/Raw_Backup/Kulkarni/Kulkarni_10X/fastqs/Kulkarni_LMMP_SK-TH/Kulkarni_SK-TH_MissingLibrary_1_HKJ5JBCXY/SK-TH_S1_L002_R2_002.fastq.gz \
/mnt/morbo/Data/Sequencing/Raw_Backup/Kulkarni/Kulkarni_10X/fastqs/Kulkarni_LMMP_SK-TH/Kulkarni_SK-TH_MissingLibrary_1_HKJ5JBCXY/SK-TH_S1_L002_R1_003.fastq.gz \
/mnt/morbo/Data/Sequencing/Raw_Backup/Kulkarni/Kulkarni_10X/fastqs/Kulkarni_LMMP_SK-TH/Kulkarni_SK-TH_MissingLibrary_1_HKJ5JBCXY/SK-TH_S1_L002_R2_003.fastq.gz \
/mnt/morbo/Data/Sequencing/Raw_Backup/Kulkarni/Kulkarni_10X/fastqs/Kulkarni_LMMP_SK-TH/Kulkarni_SK-TH_MissingLibrary_1_HKJ5JBCXY/SK-TH_S1_L002_R1_004.fastq.gz \
/mnt/morbo/Data/Sequencing/Raw_Backup/Kulkarni/Kulkarni_10X/fastqs/Kulkarni_LMMP_SK-TH/Kulkarni_SK-TH_MissingLibrary_1_HKJ5JBCXY/SK-TH_S1_L002_R2_004.fastq.gz

echo "Processing fastqs from TL"
kallisto bus -i ../bustools_required_files_10xv2_mm/mm_cDNA_introns_mm_10xv2.idx -o ./TL \
-x 10xv2 -t 16 \
/mnt/morbo/Data/Sequencing/Raw_Backup/Kulkarni/Kulkarni_10X/fastqs/Kulkarni_LMMP_SK-TL/Kulkarni_SK-TL_MissingLibrary_1_HHWFCBCXY/SK-TL_S1_L001_R1_001.fastq.gz \
/mnt/morbo/Data/Sequencing/Raw_Backup/Kulkarni/Kulkarni_10X/fastqs/Kulkarni_LMMP_SK-TL/Kulkarni_SK-TL_MissingLibrary_1_HHWFCBCXY/SK-TL_S1_L001_R2_001.fastq.gz \
/mnt/morbo/Data/Sequencing/Raw_Backup/Kulkarni/Kulkarni_10X/fastqs/Kulkarni_LMMP_SK-TL/Kulkarni_SK-TL_MissingLibrary_1_HHWFCBCXY/SK-TL_S1_L001_R1_002.fastq.gz \
/mnt/morbo/Data/Sequencing/Raw_Backup/Kulkarni/Kulkarni_10X/fastqs/Kulkarni_LMMP_SK-TL/Kulkarni_SK-TL_MissingLibrary_1_HHWFCBCXY/SK-TL_S1_L001_R2_002.fastq.gz \
/mnt/morbo/Data/Sequencing/Raw_Backup/Kulkarni/Kulkarni_10X/fastqs/Kulkarni_LMMP_SK-TL/Kulkarni_SK-TL_MissingLibrary_1_HHWFCBCXY/SK-TL_S1_L001_R1_003.fastq.gz \
/mnt/morbo/Data/Sequencing/Raw_Backup/Kulkarni/Kulkarni_10X/fastqs/Kulkarni_LMMP_SK-TL/Kulkarni_SK-TL_MissingLibrary_1_HHWFCBCXY/SK-TL_S1_L001_R2_003.fastq.gz \
/mnt/morbo/Data/Sequencing/Raw_Backup/Kulkarni/Kulkarni_10X/fastqs/Kulkarni_LMMP_SK-TL/Kulkarni_SK-TL_MissingLibrary_1_HHWFCBCXY/SK-TL_S1_L002_R1_001.fastq.gz \
/mnt/morbo/Data/Sequencing/Raw_Backup/Kulkarni/Kulkarni_10X/fastqs/Kulkarni_LMMP_SK-TL/Kulkarni_SK-TL_MissingLibrary_1_HHWFCBCXY/SK-TL_S1_L002_R2_001.fastq.gz \
/mnt/morbo/Data/Sequencing/Raw_Backup/Kulkarni/Kulkarni_10X/fastqs/Kulkarni_LMMP_SK-TL/Kulkarni_SK-TL_MissingLibrary_1_HHWFCBCXY/SK-TL_S1_L002_R1_002.fastq.gz \
/mnt/morbo/Data/Sequencing/Raw_Backup/Kulkarni/Kulkarni_10X/fastqs/Kulkarni_LMMP_SK-TL/Kulkarni_SK-TL_MissingLibrary_1_HHWFCBCXY/SK-TL_S1_L002_R2_002.fastq.gz \
/mnt/morbo/Data/Sequencing/Raw_Backup/Kulkarni/Kulkarni_10X/fastqs/Kulkarni_LMMP_SK-TL/Kulkarni_SK-TL_MissingLibrary_1_HHWFCBCXY/SK-TL_S1_L002_R1_003.fastq.gz \
/mnt/morbo/Data/Sequencing/Raw_Backup/Kulkarni/Kulkarni_10X/fastqs/Kulkarni_LMMP_SK-TL/Kulkarni_SK-TL_MissingLibrary_1_HHWFCBCXY/SK-TL_S1_L002_R2_003.fastq.gz 

echo "we got here"
