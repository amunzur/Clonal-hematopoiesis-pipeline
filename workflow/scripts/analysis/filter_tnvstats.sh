#!/bin/bash

# Filter the tnvstats based on variant calling from R

PATH_variants=$1
PATH_tnvstats=$2
identifier=$3 # a string to distinguish the filtered tnvstat names in case we have name clashes

# PATH_variants="/groups/wyattgrp/users/amunzur/pipeline/results/temp/snv.tsv"
# PATH_tnvstats="/groups/wyattgrp/users/amunzur/pipeline/results/metrics/tnvstats/kidney_samples/GU-18-088-cfDNA-UMI-2018Jul23.bam/GU-18-088-cfDNA-UMI-2018Jul23.bam.tnvstat" # one single sample

cat $PATH_variants | grep $(basename ${PATH_tnvstats//.bam.tnvstat/}) | cut -f5 | parallel --env PATH_tnvstats "cat $PATH_tnvstats | grep -w {}" > ${PATH_tnvstats//.tnvstat/.${identifier}.TEMP.tnvstat}
cat <(head -n 1 ${PATH_tnvstats}) <(cat ${PATH_tnvstats//.tnvstat/.${identifier}.TEMP.tnvstat}) > ${PATH_tnvstats//.tnvstat/.${identifier}.filtered.tnvstat} # add the header to the filtered file
rm ${PATH_tnvstats//.tnvstat/.${identifier}.TEMP.tnvstat}



