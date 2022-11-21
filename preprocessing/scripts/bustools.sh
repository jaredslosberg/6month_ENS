#!/bin/bash
source activate kallistobus

cd /data/users/jared/ENS/6mo_LMMP/preprocessing/bus_output/TH

echo TH
mkdir cDNA_capture/ introns_capture/ spliced/ unspliced/ tmp/
bustools correct -w ../../bustools_required_files_10xv2_mm/737K-august-2016.txt -p output.bus | bustools sort -o output.correct.sort.bus -t 4 -
bustools capture -s -x -o cDNA_capture/spliced.bus -c ../../bustools_required_files_10xv2_mm/introns_tx_to_capture.txt  -e matrix.ec -t transcripts.txt output.correct.sort.bus
bustools capture -s -x -o introns_capture/unspliced.bus -c ../../bustools_required_files_10xv2_mm/cDNA_tx_to_capture.txt -e matrix.ec -t transcripts.txt output.correct.sort.bus
bustools count -o unspliced/u -g ../../bustools_required_files_10xv2_mm/tr2g.tsv -e matrix.ec -t transcripts.txt --genecounts introns_capture/unspliced.bus 
bustools count -o spliced/s -g ../../bustools_required_files_10xv2_mm/tr2g.tsv -e matrix.ec -t transcripts.txt --genecounts cDNA_capture/spliced.bus



cd /data/users/jared/ENS/6mo_LMMP/preprocessing/bus_output/TL

echo TL
mkdir cDNA_capture/ introns_capture/ spliced/ unspliced/ tmp/
bustools correct -w ../../bustools_required_files_10xv2_mm/10xv2_whitelist.txt -p output.bus | bustools sort -o output.correct.sort.bus -t 4 -
bustools capture -s -x -o cDNA_capture/spliced.bus -c ../../bustools_required_files_10xv2_mm/introns_tx_to_capture.txt  -e matrix.ec -t transcripts.txt output.correct.sort.bus
bustools capture -s -x -o introns_capture/unspliced.bus -c ../../bustools_required_files_10xv2_mm/cDNA_tx_to_capture.txt -e matrix.ec -t transcripts.txt output.correct.sort.bus
bustools count -o unspliced/u -g ../../bustools_required_files_10xv2_mm/tr2g.tsv -e matrix.ec -t transcripts.txt --genecounts introns_capture/unspliced.bus 
bustools count -o spliced/s -g ../../bustools_required_files_10xv2_mm/tr2g.tsv -e matrix.ec -t transcripts.txt --genecounts cDNA_capture/spliced.bus 
