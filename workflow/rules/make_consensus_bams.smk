rule FastqToBam:
    input:
        R1=DIR_trimmed_fastq + "/{wildcard}_R1.fastq",
        R2=DIR_trimmed_fastq + "/{wildcard}_R2.fastq",
    output:
        DIR_bams + "/uBAM/{wildcard}.bam",
    params:
        sample="{wildcard}",
    run:
        shell(
            "fgbio FastqToBam -Xmx20G \
        --input {input.R1} {input.R2} \
        --output {output} \
        --read-structure 3M2S+T 3M2S+T \
        --umi-tag RX \
        --sample {params.sample} \
        --library lib1 \
        --platform illumina \
        --sort true"
        )


# Generate mapped bam
rule mapBAM:
    input:
        DIR_bams + "/uBAM/{wildcard}.bam",
    params:
        PATH_hg38=PATH_hg38,
    output:
        DIR_bams + "/mBAM/{wildcard}.bam",
    threads: 12
    run:
        shell(
            "picard SamToFastq -Xmx20G I={input} INCLUDE_NON_PF_READS=true INCLUDE_NON_PRIMARY_ALIGNMENTS=true F=/dev/stdout INTERLEAVE=true | \
        bwa mem {params.PATH_hg38} /dev/stdin -p -P -Y -t {threads} | \
        picard AddOrReplaceReadGroups -Xmx20G I=/dev/stdin O=/dev/stdout RGID=A RGSM=$sample RGPL=illumina RGLB=lib1 RGPU=unit1 SORT_ORDER=queryname | \
        picard MergeBamAlignment -Xmx20G UNMAPPED={input} ALIGNED=/dev/stdin O={output} R={params.PATH_hg38} CLIP_OVERLAPPING_READS=false CLIP_ADAPTERS=false ALIGNER_PROPER_PAIR_FLAGS=true MAX_INSERTIONS_OR_DELETIONS=-1 CREATE_INDEX=true"
        )


rule indel_realignment:
    input:
        MAPPED_bam=DIR_bams + "/mBAM/{wildcard}.bam",
        PATH_bed=PATH_bed,
        PATH_hg38=PATH_hg38,
    params:
        min_mapping_quality=20,
    output:
        DIR_bams + "/abra2/{wildcard}.bam",
    threads: 12
    run:
        shell(
            "java -Xms64G -jar /home/amunzur/anaconda3/envs/snakemake/share/abra2-2.24-1/abra2.jar \
        --in {input.SORTED_bam} --out {output} --ref {input.PATH_hg38} --threads {threads} --index \
        --no-edge-ci --nosort --mmr=0.1 --cons --targets {input.PATH_bed} --tmpdir /groups/wyattgrp/users/amunzur/COMPOST_BIN > \
        '/groups/wyattgrp/users/amunzur/pipeline/results/logs_slurm/indel_realignment/{wildcards.wildcard}'"
        )


# tee sends the stdout (the fixmate bam to multiple commands to avoid writing intermediate files)
rule fixmate_and_recalibrate_bases:
    input:
        DIR_bams + "/abra2/{wildcard}.bam",
    output:
        base_scores=DIR_recalibrated_base_scores + "/{wildcard}.table",
        uncollapsed_BAM=DIR_bams + "/uncollapsed_BAM/{wildcard}.bam",
    params:
        PATH_hg38=PATH_hg38,
        PATH_known_indels=PATH_known_indels,
        PATH_gold_std_indels=PATH_gold_std_indels,
        PATH_SNP_db=PATH_SNP_db,
    threads: 12
    run:
        shell(
            "picard -Xmx40g FixMateInformation I={input.ABRA2_bam} O=/dev/stdout IGNORE_MISSING_MATES=true SORT_ORDER=coordinate | \
        tee >(gatk BaseRecalibrator -I /dev/stdin -R {params.PATH_hg38} --known-sites {params.PATH_known_indels} --known-sites {params.PATH_gold_std_indels} --known-sites {params.PATH_SNP_db} -O {output.base_scores}) \
        tee >(gatk ApplyBQSR --input /dev/stdin --output {output.uncollapsed_BAM} --bqsr-recal-file {output.base_scores})"
        )


# Identify reads or read pairs that originate from the same source molecule based on genomic positions and UMI
rule GroupReadsByUmi:
    input:
        DIR_bams + "/uncollapsed_BAM/{wildcard}.bam",  # output of the fixmate_and_recalibrate_bases rule from the process_bams.smk file
    output:
        DCS_uBAM=DIR_bams + "/DCS_uBAM/{wildcard}.bam",
        SSCS_uBAM=DIR_bams + "/SSCS_uBAM/{wildcard}.bam",
    params:
        PATH_hg38=PATH_hg38,
    threads: 12
    run:
        shell(
            "fgbio GroupReadsByUmi -Xmx20G --input={input} --output=/dev/stdout --strategy=paired --edits=1 | \
        tee >(fgbio CallDuplexConsensusReads -Xmx20G --input=/dev/stdin --output={output.DCS_uBAM} --threads={threads} --read-name-prefix=duplex --sort-order=Queryname) \
        tee >(fgbio CallMolecularConsensusReads -Xmx20G --input=/dev/stdin --output={output.SSCS_uBAM} --threads={threads} --min-reads=1 --read-name-prefix=singlex --sort-order=Queryname)"
        )
