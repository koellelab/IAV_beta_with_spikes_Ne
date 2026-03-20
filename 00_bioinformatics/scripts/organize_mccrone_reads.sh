#!/bin/bash

#SBATCH --job-name=organize_mccrone
#SBATCH --partition=day-long-cpu
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1

while read p; do
  array=(${p//,/ })
  new_name=${array[1]}_${array[3]}_${array[0]}
  mkdir data/mccrone/reads/${new_name}
  scp data/mccrone/reads/${array[0]}_1.fastq.gz data/organized_reads/mccrone_${new_name}/${new_name}_R1.fastq.gz
  scp data/mccrone/reads/${array[0]}_2.fastq.gz data/organized_reads/mccrone_${new_name}/${new_name}_R2.fastq.gz
  rm -rf data/mccrone/reads/${array[0]}_1.fastq.gz
  rm -rf data/mccrone/reads/${array[0]}_2.fastq.gz
done <  <(tail -n +2 data/mccroneIAV_analysis_files/mccrone_SraRunTable_longitudinal.csv)

