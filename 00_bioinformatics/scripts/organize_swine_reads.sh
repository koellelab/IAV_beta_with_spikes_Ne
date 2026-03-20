#!/bin/bash

#SBATCH --job-name=organize_swine
#SBATCH --partition=swineIAV 
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

while read p; do
  id=$(echo $p | awk -F',' '{print $47}' | sed 's/pig//g' | awk -F'_' '{print $1}')
  date=$(echo $p | awk -F',' '{print $14}')
  acc=$(echo $p | awk -F',' '{print $1}')
  mkdir data/organized_reads/swine_${id}_${date}_${acc}
  scp data/swine/reads/${acc}_1.fastq.gz data/organized_reads/swine_${id}_${date}_${acc}/swine_${id}_${date}_${acc}_R1.fastq.gz
  scp data/swine/reads/${acc}_2.fastq.gz data/organized_reads/swine_${id}_${date}_${acc}/swine_${id}_${date}_${acc}_R2.fastq.gz
  rm -rf data/swine/reads/${acc}_1.fastq.gz 
  rm -rf data/swine/reads/${acc}_2.fastq.gz 
done  < <(grep -f <(awk -F',' '{print $1}' data/vaninsberghe_analysis_files/swine_to_get.csv) data/vaninsberghe_analysis_files/swine_SraRunTable.txt)

