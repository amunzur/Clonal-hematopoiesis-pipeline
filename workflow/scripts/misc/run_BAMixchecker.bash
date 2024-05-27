#!/bin/bash
# Runs the BAMixchecker on matched cfDNA and WBC files to ensure they come from the same person. 
# The tool depends on bedtools and GATK, which are both available in the base environment.

PATH_bamixchecker="/home/amunzur/BAMixChecker/BAMixChecker.py"
DIR_bams="/groups/wyattgrp/users/amunzur/pipeline/results/data/bam/SSCS1_final"
PATH_hg38="/groups/wyattgrp/users/jbacon/reference/matti_hg38/hg38_masked.fa"
PATH_bed="/groups/wyattgrp/users/amunzur/pipeline/resources/panel/CHIP/CHIP_targets.bed"
PATH_bam_list="/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/paired_bams_bladder_BAMixchecker.tsv"
DIR_output="/groups/wyattgrp/users/amunzur/pipeline/results/metrics/SNP_matching"

python "${PATH_bamixchecker}" \
 -l ${PATH_bam_list} \
 -p 20 \
 -r "${PATH_hg38}" \
 -o "${DIR_output}" \
 -v hg38 \
 -b "${PATH_bed}"