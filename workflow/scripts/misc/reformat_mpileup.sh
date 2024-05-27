DIR_mpileup="/groups/wyattgrp/users/amunzur/pipeline/results/metrics/mpileup/temp_SSCS1"
DIR_output="/groups/wyattgrp/users/amunzur/pipeline/results/metrics/mpileup/SSCS1"

mkdir -p $DIR_output

for file in $DIR_mpileup/*.mpileup
do
    sample_name=$(basename ${file})
    echo ${file}
    paste <(awk '{print $1$2}' ${file}) <(cat ${file}) > ${DIR_output}/${sample_name}
done

