#!/bin/bash

# This script compiles sequencing quality metrics into one csv file, containing information such as the number of reads taken from the raw fastq files, 
# averaged and median depth of aligned and processed bams and the duplicate percentages. 

batch=$1

DIR_markdup="/groups/wyattgrp/users/amunzur/pipeline/results/metrics/PICARD_markdup_perc/${batch}"
PATH_depth="/groups/wyattgrp/users/amunzur/pipeline/results/metrics/averaged_depth/${batch}/averaged_depths.txt"
DIR_readcounts="/groups/wyattgrp/users/amunzur/pipeline/results/metrics/read_counts/raw/${batch}"
output="/groups/wyattgrp/users/amunzur/pipeline/results/metrics/seq_quality_metrics/${batch}.csv"

paste \
	<(cat ${PATH_depth} | tail -n +2 | sort -n -k 1,1) \
	<(cat ${DIR_markdup}/*.txt | sort -n -k 1,1 | cut -f2) \
	<(cat ${DIR_readcounts}/*.txt | sort -n -k 1,1 | cut -f2 ) |\
	sed 's/\t/,/g' |\
    awk -F , -v OFS=, '$4*=100, $5/=500000' > "${output}_temp" # convert the duplicate fractions to percentage 

cat <(echo -e "Sample_name,Averaged_depth,Median_depth,Duplicates(%),Total_reads(millions)") <(cat "${output}_temp") > ${output}
rm "${output}_temp"