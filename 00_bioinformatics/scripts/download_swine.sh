#!/bin/bash

#SBATCH --job-name=download_mccrone
#SBATCH --partition=swineIAV
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1

while read p; do
  ../sratoolkit.3.1.1-centos_linux64/bin/fasterq-dump $p -O data/swine/reads -e 4
  gzip data/swine/reads/*.fastq
done < <(grep -f data/vaninsberghe_analysis_files/swine_to_get.csv data/vaninsberghe_analysis_files/swine_SraRunTable.txt | awk -F',' '{print $1}')

