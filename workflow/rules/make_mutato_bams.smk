rule FastqToBam:
    input:
        R1=DIR_trimmed_fastq + "/{wildcard}_R1.fastq",
        R2=DIR_trimmed_fastq + "/{wildcard}_R2.fastq",
    output:
        DIR_bams_mutato + "/uBAM/{wildcard}.bam",
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

# Go back to the Fastq format
rule BamtoFastq:
    input:
        DIR_bams_mutato + "/uBAM/{wildcard}.bam",
    params:
        PATH_hg38=PATH_hg38,
    output:
        temp(DIR_fastq + "/consensus_fq/{wildcard}.fastq"),
    threads: 12
    run:
        shell(
            "picard SamToFastq -Xmx20G I={input} F={output} INCLUDE_NON_PF_READS=true INCLUDE_NON_PRIMARY_ALIGNMENTS=true INTERLEAVE=true"
        )

# Generate mapped bam
rule mapBAM:
    input:
        DIR_fastq + "/consensus_fq/{wildcard}.fastq",
    params:
        PATH_hg38=PATH_hg38,
    output:
        temp(DIR_bams_mutato + "/mSAM_raw/{wildcard}.bam")
    threads: 12
    run:
        shell("bowtie -x /groups/wyattgrp/users/amunzur/pipeline/resources/references/index/Homo_sapiens/NCBI/GRCh38/Sequence/Bowtie2Index/genome {input} > {output} ")

rule SAMtoBAM:
    input:
        DIR_bams_mutato + "/mSAM_raw/{wildcard}.bam"
    output:
        DIR_bams_mutato + "/mBAM_raw/{wildcard}.bam"
    threads: 12
    run:
        shell("samtools view -b -h {input} > {output} ")


# Add read groups to the mapped bam file
rule addRG:
    input:
        DIR_bams_mutato + "/mBAM_raw/{wildcard}.bam",
    params:
        PATH_hg38=PATH_hg38,
        sample="{wildcard}"
    output:
        temp(DIR_bams_mutato + "/readGroup/{wildcard}.bam"),
    threads: 12
    run:
        shell(
            "picard AddOrReplaceReadGroups -Xmx20G I={input} O={output} RGID=A RGSM={params.sample} RGPL=illumina RGLB=lib1 RGPU=unit1 SORT_ORDER=queryname"
        )

rule MergeBamAlignment:
    input:
        mBAM=DIR_bams_mutato + "/readGroup/{wildcard}.bam",
        uBAM=DIR_bams_mutato + "/uBAM/{wildcard}.bam"
    params:
        PATH_hg38=PATH_hg38,
    output:
        DIR_bams_mutato + "/mBAM/{wildcard}.bam",
    threads: 12
    run:
        shell(
            "picard MergeBamAlignment -Xmx20G UNMAPPED={input.uBAM} ALIGNED={input.mBAM} O={output} R={params.PATH_hg38} CLIP_OVERLAPPING_READS=false CLIP_ADAPTERS=false ALIGNER_PROPER_PAIR_FLAGS=true MAX_INSERTIONS_OR_DELETIONS=-1 CREATE_INDEX=true"
        )


rule fixmate:
    input:
        DIR_bams_mutato + "/mBAM/{wildcard}.bam",
    output:
        temp(DIR_bams_mutato + "/fixmate/{wildcard}.bam"),
    threads: 12
    run:
        shell(
            "picard -Xmx40g FixMateInformation I={input} O={output} IGNORE_MISSING_MATES=true SORT_ORDER=coordinate"
        )

rule recalibrate_bases:
    input:
        DIR_bams_mutato + "/fixmate/{wildcard}.bam",
    output:
        DIR_recalibrated_base_scores + "/{wildcard}_table",
    params:
        PATH_hg38=PATH_hg38,
        PATH_known_indels=PATH_known_indels,
        PATH_gold_std_indels=PATH_gold_std_indels,
        PATH_SNP_db=PATH_SNP_db,
    threads: 12
    run:
        shell(
            "~/gatk-4.2.0.0/gatk BaseRecalibrator -I {input} -R {params.PATH_hg38} --known-sites {params.PATH_known_indels} --known-sites {params.PATH_gold_std_indels} --known-sites {params.PATH_SNP_db} -O {output}"
        )

rule apply_base_scores:
    input:
        fixmate_BAM=DIR_bams_mutato + "/fixmate/{wildcard}.bam",
        base_scores=DIR_recalibrated_base_scores + "/{wildcard}_table",
    output:
        uncollapsed_BAM =DIR_bams_mutato + "/uncollapsed_BAM/{wildcard}.bam",
    params:
        PATH_hg38=PATH_hg38,
        PATH_known_indels=PATH_known_indels,
        PATH_gold_std_indels=PATH_gold_std_indels,
        PATH_SNP_db=PATH_SNP_db,
    threads: 12
    run:
        shell(
        "~/gatk-4.2.0.0/gatk ApplyBQSR --input {input.fixmate_BAM} --output {output.uncollapsed_BAM} --bqsr-recal-file {input.base_scores}"
        )


# Identify reads or read pairs that originate from the same source molecule based on genomic positions and UMI
rule GroupReadsByUmi:
    input:
        DIR_bams_mutato + "/uncollapsed_BAM/{wildcard}.bam",  # output of the fixmate_and_recalibrate_bases rule from the process_bams.smk file
    output:
        temp(DIR_bams_mutato + "/grouped_umi_BAM/{wildcard}.bam"),
    params:
        PATH_hg38=PATH_hg38,
    threads: 12
    run:
        shell(
            "fgbio GroupReadsByUmi -Xmx20G --input={input} --output={output} --strategy=paired --edits=1"
        )

rule CallDuplexConsensusReads:
    input:
        DIR_bams_mutato + "/grouped_umi_BAM/{wildcard}.bam",  # output of the fixmate_and_recalibrate_bases rule from the process_bams.smk file
    output:
        DIR_bams_mutato + "/DCS_uBAM/{wildcard}.bam",
    threads: 12
    run:
        shell(
        "fgbio CallDuplexConsensusReads -Xmx20G --input={input} --output={output} --threads={threads} --read-name-prefix=duplex --sort-order=Queryname"
        )

rule CallMolecularConsensusReads:
    input:
        DIR_bams_mutato + "/grouped_umi_BAM/{wildcard}.bam",  # output of the fixmate_and_recalibrate_bases rule from the process_bams.smk file
    output:
        DIR_bams_mutato + "/SSCS_uBAM/{wildcard}.bam",
    threads: 12
    run:
        shell(
        "fgbio CallMolecularConsensusReads -Xmx20G --input={input} --output={output} --threads={threads} --min-reads=1 --read-name-prefix=singlex --sort-order=Queryname"
        )

# Generate mapped DCS and SSCS bams from the unmapped consensus bam files.
rule BamtoFastq2:
    input:
        DIR_bams_mutato + "/{consensus_type}_uBAM/{wildcard}.bam",
    output:
        temp(DIR_fastq + "/{consensus_type}_consensus_fq_2/{wildcard}.fastq"),
    threads: 12
    run:
        shell(
            "picard SamToFastq -Xmx20G I={input} INCLUDE_NON_PRIMARY_ALIGNMENTS=true F={output} INTERLEAVE=true"
        )

rule mapBAM2:
    input:
        DIR_fastq + "/{consensus_type}_consensus_fq_2/{wildcard}.fastq",
    params:
        PATH_hg38=PATH_hg38,
    output:
        temp(DIR_bams_mutato + "/{consensus_type}_mBAM_raw_2/{wildcard}.bam"),
    threads: 12
    run:
        shell("bowtie -x /groups/wyattgrp/users/amunzur/pipeline/resources/references/index/Homo_sapiens/NCBI/GRCh38/Sequence/Bowtie2Index/genome {input} > {output} ")

rule addRG2:
    input:
        DIR_bams_mutato + "/{consensus_type}_mBAM_raw_2/{wildcard}.bam",
    params:
        sample="{wildcard}",
    output:
        temp(DIR_bams_mutato + "/{consensus_type}_readGroup_2/{wildcard}.bam"),
    threads: 12
    run:
        shell(
            "picard AddOrReplaceReadGroups -Xmx20G I={input} O={output} RGID=A RGSM={params.sample} RGPL=illumina RGLB=lib1 RGPU=unit1 SORT_ORDER=queryname"
        )


rule MergeBamAlignment2:
    input:
        mBAM=DIR_bams_mutato + "/{consensus_type}_readGroup_2/{wildcard}.bam",
        uBAM=DIR_bams_mutato + "/{consensus_type}_uBAM/{wildcard}.bam",
    params:
        PATH_hg38=PATH_hg38,
    output:
        temp(DIR_bams_mutato + "/{consensus_type}_mBAM/{wildcard}.bam"),
    threads: 12
    run:
        shell(
            "picard MergeBamAlignment -Xmx20G UNMAPPED={input.uBAM} ALIGNED={input.mBAM} O={output} R={params.PATH_hg38} CLIP_OVERLAPPING_READS=false CLIP_ADAPTERS=false ALIGNER_PROPER_PAIR_FLAGS=true MAX_INSERTIONS_OR_DELETIONS=-1 CREATE_INDEX=true"
        )

# Good bams we want to keep before doing any downstream filtering on them, so output shouldn't be temp
rule fixmate2:
    input:
        DIR_bams_mutato + "/{consensus_type}_mBAM/{wildcard}.bam",
    output:
        DIR_bams_mutato + "/{consensus_type}_final/{wildcard}.bam"
    threads: 12
    run:
        shell(
            "picard -Xmx40g FixMateInformation I={input} O={output} IGNORE_MISSING_MATES=true SORT_ORDER=coordinate"
        )


rule FilterConsensusReads_DCS:
    input:
        DIR_bams_mutato + "/DCS_final/{wildcard}.bam",
    params:
        PATH_hg38=PATH_hg38,
        PATH_bed=PATH_bed,
        sample="{wildcard}",
        min_reads="2 1 1",
        min_base_quality=45,
    output:
        temp(DIR_bams_mutato + "/DCS_filtered/{wildcard}.bam"),
    threads: 12
    run:
        shell(
            "samtools sort -n -u {input} | fgbio FilterConsensusReads -Xmx20G --input=/dev/stdin --output={output} --ref={params.PATH_hg38} --min-reads={params.min_reads} --min-base-quality={params.min_base_quality} --require-single-strand-agreement=true --sort-order=Coordinate"
        )

# Need to update the min_read here, and change the dir to the outputted bams as well
rule FilterConsensusReads_SSCS3:
    input:
        DIR_bams_mutato + "/SSCS_final/{wildcard}.bam",
    params:
        PATH_hg38=PATH_hg38,
        PATH_bed=PATH_bed,
        sample="{wildcard}",
        min_reads=3,
        min_base_quality=40,
        max_read_error_rate=0.05,
    output:
        temp(DIR_bams_mutato + "/SSCS3_filtered/{wildcard}.bam"),
    threads: 12
    run:
        shell(
            "samtools sort -n -u {input} | fgbio FilterConsensusReads -Xmx20G --input=/dev/stdin --output={output} --ref={params.PATH_hg38} --min-reads={params.min_reads} \
        --max-read-error-rate={params.max_read_error_rate} --min-base-quality={params.min_base_quality} --reverse-per-base-tags=true --sort-order=Coordinate"
        )

rule subset_to_panel:
    input:
        BAM=DIR_bams_mutato + "/{consensus_type}_filtered/{wildcard}.bam",
        PATH_bed=PATH_bed,
    output:
        DIR_bams_mutato + "/{consensus_type}_panel/{wildcard}.bam"
    threads: 12
    run:
        shell(
            "samtools view -b -h -L {input.PATH_bed} {input.BAM} > {output}"
        )