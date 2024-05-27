# Runs cnvkit on the samples. Need to run conda activate cnvkit.

DIR_bam="/groups/wyattgrp/users/amunzur/pipeline/results/data/bam/cnvkit_bams"
VIP_samples="/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_list_normals.tsv"
PATH_hg38="/groups/wyattgrp/users/jbacon/reference/matti_hg38/hg38_masked.fa"
DIR_results="/groups/wyattgrp/users/amunzur/pipeline/results"

# STEP 3. CONSTRUCT A REFERENCE USING THE cfDNA VIP NORMALS
cnvkit.py reference "${DIR_results}/cnvkit/coverage/SSCS2/"*VIP*cfDNA*coverage.cnn -f ${PATH_hg38} -o "${DIR_results}/cnvkit/references/VIP_cfDNA_reference.cnn" &

# STEP 4. CONSTRUCT A REFERENCE USING THE WBC VIP NORMALS
cnvkit.py reference "${DIR_results}/cnvkit/coverage/SSCS2/"*VIP*WBC*coverage.cnn -f ${PATH_hg38} -o "${DIR_results}/cnvkit/references/VIP_WBC_reference.cnn" &

# For each sample correct for GC bias and
