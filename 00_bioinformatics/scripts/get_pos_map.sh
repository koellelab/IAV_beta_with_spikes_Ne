#!/usr/bin/env bash

# aggregate fasta files
# from my analysis
mkdir -p data
mkdir -p data/genomes
for seg in PB2 PB1 PA HA NP NA MP NS; do 
	rm -rf data/genomes/${seg}.fasta
	echo  > data/genomes/${seg}.fasta
	for i in data/IRMA/*/*${seg}*.fasta; do
		id=$(echo $i | cut -d"/" -f4)
		echo ">$id" >> data/genomes/${seg}.fasta
		tail -n +2 $i >> data/genomes/${seg}.fasta
	done
done

# add fasta files from Van Insberghe Plos
# first convert to tsv for easier sorting
awk -F'\t' '{
	if (substr($1, 1,1) == ">")
		if (NR > 1){
			printf "\n%s\t",substr($0,2,length($0)-1)
		}
		else{
			printf "%s\t",substr($0,2,length($0)-1)
		}
	else{
		printf "%s",$0
	}
	}END{printf "\n"}' data/vaninsberghe_analysis_files/all_refs.fasta > tmp.tsv

while read p; do
	seg=$(echo $p | cut -f1 | awk '{print $1}' | cut -d"_" -f2 | sed 's/NEP/NS/g')
	echo $p | sed 's/_/|/g' | awk '{print ">swine_"$1"\n"$2}' >> data/genomes/${seg}.fasta
done < tmp.tsv


# todo only do splitting operation once
rm tmp.tsv
for i in data/mccroneIAV_analysis_files/*.fa; do
	label=$(echo $i | cut -d/ -f3 | cut -d_ -f1 | sed 's/Hong/HK/g')
	awk -F'\t' -v label="$label" '{
	if (substr($1, 1,1) == ">")
		if (NR > 1){
			split(substr($1,2,10), a, " ")
			printf "\n%s\t",label"_"a[1]
		}
		else{
			split(substr($1,2,10), a, " ")
			printf "%s\t",label"_"a[1]
		}
	else{
		printf "%s",$0
	}
	}END{printf "\n"}' $i >> tmp.tsv
done

while read p; do
	seg=$(echo $p | awk '{print $1}' | cut -d"_" -f2 | sed 's/NEP/NS/g')
	echo $p | sed 's/_/|/g' | awk '{print ">mccrone_"$1"\n"$2}' >> data/genomes/${seg}.fasta
done < <(cat tmp.tsv | sed 's/\_NR/\_NA/g')


# align and generate a map between references
for i in PB2 PB1 PA HA NP NA MP NS; do 
	echo $i
	mafft --auto --thread -4 data/genomes/${i}.fasta > data/genomes/${i}_aln.fasta
	python3 scripts/calc_pos_map.py --seqs data/genomes/${i}_aln.fasta
done


# combine all position maps
head -n 1 data/genomes/PB2_aln_map.tsv | awk -F'\t' '{print "seg\t"$0}'> data/genomes/all_aln_map.tsv
for i in PB2 PB1 PA HA NP NA MP NS; do 
	tail -n +2 data/genomes/${i}_aln_map.tsv | awk -F'\t' -v var=$i '{print var"\t"$0}' >> data/genomes/all_aln_map.tsv
done