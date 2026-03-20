#!/bin/bash

# count number of reads and get reference genome
rm -rf data/read_counts.tsv
for i in data/IRMA/*/reads/*_R1.fastq.gz; do
	subject=$(echo $i | cut -d"/" -f3 | cut -d"_" -f1)_$(echo $i | cut -d"/" -f3 | cut -d"_" -f2)
	segment=$(echo $i | cut -d"/" -f5 | sed "s/_R1.fastq.gz//g")
	path=$(echo $i | sed "s/_R1.fastq.gz//g")
	read_count=$(gunzip -c $i | wc -l)
	echo -e $subject"\t"$segment"\t"$path"\t"$read_count >> data/read_counts.tsv
done

# sum across segments get highest read sample from the earliest date
awk -F'\t' '{split($3, b, "/"); a[b[3]] += $4}END{for (i in a) print(i"\t"a[i])}' data/read_counts.tsv \
	> data/read_counts_sum.tsv

# for mccrone data
# get list of which have HS in their identifier (indicating household sample)
# this is always the reference

# for swine data
# first get those with the earliest sampling date 
# then get the one with the most reads
cat \
	<(awk -F',' 'FNR==NR{if (substr($(NF-1), 1, 2) == "HS"){a[$1]=1}; next}{
			split($1, c, "\t"); split(c[1], d, "_"); if (d[4] in a) print d[1]"_"d[2]"\t"d[3]"\t"c[1]"\t"c[2]}' \
		data/mccroneIAV_analysis_files/mccrone_SraRunTable.csv \
		data/read_counts_sum.tsv  | \
		sort -nr -k 4 | \
		awk -F'\t' '!seen[$1]++' | \
		sort -k 1) \
	<(awk -F'\t' 'FNR==NR{a[$1] = $2; next}
		{split($1, b, "_"); if (a[b[1]"_"b[2]] == $1) print b[1]"_"b[2]"\t"b[3]"\t"$1"\t"$2}' \
		<(awk -F'\t' '{split($3, a, "_"); split($3, b, "/"); print$1"\t"b[3]"\t"a[3]}' <(grep "swine" data/read_counts.tsv) | sort -k 3 |  awk -F'\t' '!seen[$1]++') \
		<(grep "swine" data/read_counts_sum.tsv) | \
		sort -nr -k 4 | \
		awk -F'\t' '!seen[$1]++' | \
		sort -k 1) \
	> data/references.tsv
