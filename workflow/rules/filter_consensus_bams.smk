# Generate mapped DCS and SSCS bams from the unmapped consensus bam files.
rule remapBAM:
    input:
        DIR_bams + "/{consensus_type}_uBAM/{wildcard}.bam",
    params:
        PATH_hg38=PATH_hg38,
        sample="{wildcard}",
    output:
        DIR_bams + "/{consensus_type}_mBAM/{wildcard}.bam",
    threads: 12
    run:
        shell(
            "picard SamToFastq -Xmx20G I={input} INCLUDE_NON_PRIMARY_ALIGNMENTS=true F=/dev/stdout INTERLEAVE=true | \
        bwa mem {params.PATH_hg38} /dev/stdin -p -P -Y -t 8 | \
        picard AddOrReplaceReadGroups -Xmx20G I=/dev/stdin O=/dev/stdout RGID=A RGSM={params.sample} RGPL=illumina RGLB=lib1 RGPU=unit1 SORT_ORDER=queryname | \
        picard MergeBamAlignment -Xmx20G UNMAPPED={input} ALIGNED=/dev/stdin O={output} R={params.sample} CLIP_OVERLAPPING_READS=false CLIP_ADAPTERS=false ALIGNER_PROPER_PAIR_FLAGS=true MAX_INSERTIONS_OR_DELETIONS=-1 CREATE_INDEX=true"
        )


# Run ABRA2 and picard's fixmate on the mapped consensus bams
rule process_consensus_mBAM:
    input:
        DIR_bams + "/{consensus_type}_mBAM/{wildcard}.bam",
    params:
        PATH_hg38=PATH_hg38,
        PATH_bed=PATH_bed,
        sample="{wildcard}",
    output:
        DIR_bams + "/{consensus_type}_processed_mBAM/{wildcard}.bam",
    threads: 12
    run:
        shell(
            " java -Xms64G -jar /home/amunzur/anaconda3/envs/snakemake/share/abra2-2.24-1/abra2.jar \
        --in {input} --out /dev/stdout --ref {params.PATH_hg38} --threads {threads} --nosort --mmr=0.1 --cons --no-edge-ci --targets {params.PATH_bed} --tmpdir /groups/wyattgrp/users/amunzur/COMPOST_BIN | \
        picard FixMateInformation -Xmx20G I=/dev/stdin O=/dev/stdout IGNORE_MISSING_MATES=true SORT_ORDER=queryname"
        )


rule FilterConsensusReads_DCS:
    input:
        DIR_bams + "/DCS_processed_mBAM/{wildcard}.bam",
    params:
        PATH_hg38=PATH_hg38,
        PATH_bed=PATH_bed,
        sample="{wildcard}",
        min_reads="2 1 1",
        min_base_quality=45,
    output:
        DIR_bams + "/DCS_filtered_mBAM/{wildcard}.bam",
    threads: 12
    run:
        shell(
            "fgbio FilterConsensusReads -Xmx20G --input={input} --output={output} --ref={params.PATH_hg38} --min-reads={params.min_reads} --min-base-quality={params.min_base_quality} --require-single-strand-agreement=true --sort-order=Coordinate"
        )


# Need to update the min_read here, and change the dir to the outputted bams as well
rule FilterConsensusReads_SSCS:
    input:
        DIR_bams + "/SSCS_processed_mBAM/{wildcard}.bam",
    params:
        PATH_hg38=PATH_hg38,
        PATH_bed=PATH_bed,
        sample="{wildcard}",
        min_reads=3,
        min_base_quality=40,
        max_read_error_rate=0.05,
    output:
        DIR_bams + "/SSCS3_filtered_mBAM/{wildcard}.bam",
    threads: 12
    run:
        shell(
            "fgbio FilterConsensusReads -Xmx20G --input={input} --output={output} --ref={params.PATH_hg38} --min-reads={params.min_reads} \
        --max-read-error-rate={params.max_read_error_rate} --min-base-quality={params.min_base_quality} --reverse-per-base-tags=true --sort-order=Coordinate"
        )
