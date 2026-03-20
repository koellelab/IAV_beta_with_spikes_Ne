#!/bin/bash

#SBATCH --job-name=download_mccrone
#SBATCH --partition=week-long-cpu
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1

while read p; do
  ~/sratoolkit.3.1.1-centos_linux64/bin/fasterq-dump $p -O data/mccrone/reads -e 4
  gzip data/mccrone/reads/*.fastq
done < <(awk -F',' '{print $1}' <(tail -n +2 data/mccroneIAV_analysis_files/mccrone_SraRunTable_longitudinal.csv))

