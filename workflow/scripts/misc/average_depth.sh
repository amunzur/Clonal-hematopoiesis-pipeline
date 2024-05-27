#!/bin/bash
# Using the per position depth metrics, this script will go through all depth metrics in a given directory and make one output file 
# consisting of three columns: sample_name, averaged_depth and median depth of all positions.
DIR_depth="/groups/wyattgrp/users/amunzur/pipeline/results/metrics/depth/SSCS2"
PATH_depth="/groups/wyattgrp/users/amunzur/pipeline/results/metrics/averaged_depth/SSCS2/averaged_depths.txt"

for f in ${DIR_depth}/*.txt 
do 
	sample_name=$(basename $f)
	echo $sample_name
	echo "${sample_name%.*}" $(awk '{ total += $3; count++ } END { print total/count }' ${f}) $(cut -f3 ${f} | sort -n | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }') >> ${PATH_depth}
done

# Next part compiles sequencing quality metrics into one csv file, containing information such as the number of reads taken from the raw fastq files, 
# averaged and median depth of aligned and processed bams and the duplicate percentages. 

# DIR_markdup="/groups/wyattgrp/users/amunzur/pipeline/results/metrics/PICARD_markdup_perc/"
# DIR_readcounts="/groups/wyattgrp/users/amunzur/pipeline/results/metrics/read_counts/merged/"
# output="/groups/wyattgrp/users/amunzur/pipeline/results/metrics/seq_quality_metrics/seq_quality_metrics.csv"

# paste \
# 	<(cat ${PATH_depth} | tail -n +2 | sort -n -k 1,1) \
# 	<(cat ${DIR_markdup}/*.txt | sort -n -k 1,1 | cut -f2) \
# 	<(cat ${DIR_readcounts}/*.txt | sort -n -k 1,1 | cut -f2 ) |\
# 	sed 's/\s\+/,/g' |\
#     awk -F , -v OFS=, '$4*=100, $5/=500000' > "${output}_temp" # convert the duplicate fractions to percentage 

# cat <(echo -e "Sample_name,Averaged_depth,Median_depth,Duplicates(%),Total_reads(millions)") <(cat "${output}_temp") > ${output}
# rm "${output}_temp"