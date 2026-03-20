#!/bin/bash

# download data
bash scripts/download_mccrone.sh
bash scripts/download_swine.sh

# organize data
bash scripts/organize_mccrone_reads.sh
bash scripts/organize_swine_reads.sh

# assemble genomes
snakemake \
	--snakefile variant_calling_pipeline/assemble_genomes.smk \
	--cores 8  \
	--default-resource tmpdir="'/projects/GII4evolution/mmartin/within_host_ne/tmp'" -p  

# count reads to identify reference to use
bash scripts/count_reads.sh

# call variants
snakemake \
	--snakefile variant_calling_pipeline/call_snps.smk \
	--cores 8  \
	--default-resource tmpdir="'/projects/GII4evolution/mmartin/within_host_ne/tmp'" -p  

# aggregate fasta for position mapping
bash scripts/get_pos_map.sh

# consolidate data
python3 scripts/consolidate_data.py

# filter data
Rscript scripts/process_wide_vcf.R

# compare to JM data
# foramt McCrone no_freq_cut.qual.snv.csv into array format
# https://github.com/lauringlab/Host_level_IAV_evolution/blob/master/data/processed/secondary/no_freq_cut.qual.snv.csv
python3 scripts/format_mccrone_dat.py

# get just used mccrone IDs
cat <(head -n 1 data/all_dates.tsv) \
	<(grep "mccrone" data/all_dates.tsv) \
	> data/all_dates_mccrone.tsv

python3 scripts/compare_var_freqs.py \
	--varFreq1 data/processed_vcf_arr.csv \
	--varFreq2 data/mccroneIAV_analysis_files/no_freq_cut.qual.snv_arr.csv \
	--useIDs data/all_dates_mccrone.tsv \
	--posMap data/genomes/all_aln_map.tsv \
	--labels 'Human IAV iSNV frequencies called by McCrone et al. (2018)' 'Human IAV iSNV frequencies called by our pipeline' 'jt_mm_compare'
