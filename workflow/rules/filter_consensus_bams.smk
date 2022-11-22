# Generate mapped DCS and SSCS bams from the unmapped consensus bam files.
rule BamtoFastq2:
    input:
        DIR_bams + "/{consensus_type}_uBAM/{wildcard}.bam",
    output:
        temp(DIR_fastq + "/{consensus_type}_consensus_fq_2/{wildcard}.fastq"),
    threads: 12
    run:
        shell(
            "picard SamToFastq -Xmx20G I={input} INCLUDE_NON_PRIMARY_ALIGNMENTS=true INCLUDE_NON_PF_READS=true INTERLEAVE=true F={output}"
        )


rule mapBAM2:
    input:
        DIR_fastq + "/{consensus_type}_consensus_fq_2/{wildcard}.fastq",
    params:
        PATH_hg38=PATH_hg38,
    output:
        temp(DIR_bams + "/{consensus_type}_mBAM_raw_2/{wildcard}.bam"),
    threads: 12
    run:
        shell("bwa mem {params.PATH_hg38} {input} -p -Y -t {threads} > {output}")


rule addRG2:
    input:
        DIR_bams + "/{consensus_type}_mBAM_raw_2/{wildcard}.bam",
    params:
        sample="{wildcard}",
    output:
        temp(DIR_bams + "/{consensus_type}_readGroup_2/{wildcard}.bam"),
    threads: 12
    run:
        shell(
            "picard AddOrReplaceReadGroups -Xmx20G I={input} O={output} RGID=A RGSM={params.sample} RGPL=illumina RGLB=lib1 RGPU=unit1 SORT_ORDER=queryname"
        )


rule MergeBamAlignment2:
    input:
        mBAM=DIR_bams + "/{consensus_type}_readGroup_2/{wildcard}.bam",
        uBAM=DIR_bams + "/{consensus_type}_uBAM/{wildcard}.bam",
    params:
        PATH_hg38=PATH_hg38,
    output:
        temp(DIR_bams + "/{consensus_type}_mBAM/{wildcard}.bam"),
    threads: 12
    run:
        shell(
            "picard MergeBamAlignment -Xmx20G UNMAPPED={input.uBAM} ALIGNED={input.mBAM} O={output} R={params.PATH_hg38} \
                CLIP_OVERLAPPING_READS=false \
                CLIP_ADAPTERS=false \
                ALIGNER_PROPER_PAIR_FLAGS=true \
                INCLUDE_SECONDARY_ALIGNMENTS=true \
                PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
                EXPECTED_ORIENTATIONS=FR \
                MAX_INSERTIONS_OR_DELETIONS=-1 \
                CREATE_INDEX=true"
        )


rule indel_realignment2:
    input:
        MAPPED_bam=DIR_bams + "/{consensus_type}_mBAM/{wildcard}.bam",
        PATH_bed=PATH_bed,
        PATH_hg38=PATH_hg38,
    output:
        temp(DIR_bams + "/{consensus_type}_abra2/{wildcard}.bam"),
    threads: 12
    run:
        shell(
            "java -Xms64G -jar /home/amunzur/anaconda3/envs/snakemake/share/abra2-2.24-1/abra2.jar \
        --in {input.MAPPED_bam} \
        --out {output} \
        --ref {input.PATH_hg38} \
        --threads {threads} \
        --nosort \
        --no-edge-ci \
        --targets {input.PATH_bed} \
        --tmpdir /groups/wyattgrp/users/amunzur/COMPOST_BIN > \
        '/groups/wyattgrp/users/amunzur/pipeline/results/logs_slurm/indel_realignment/{wildcards.wildcard}'"
        )


# Good bams we want to keep before doing any downstream filtering on them, so output shouldn't be temp
rule fixmate2:
    input:
        DIR_bams + "/{consensus_type}_abra2/{wildcard}.bam",
    output:
        DIR_bams + "/{consensus_type}_final/{wildcard}.bam",
    threads: 12
    run:
        shell(
            "picard -Xmx40g FixMateInformation I={input} O={output} SORT_ORDER=coordinate"
        )

rule FilterConsensusReads_SSCS:
    input:
        DIR_bams + "/SSCS_final/{wildcard}.bam",
    params:
        PATH_hg38=PATH_hg38,
        PATH_bed=PATH_bed,
        sample="{wildcard}",
        min_reads=2,
        min_base_quality=40,
        max_read_error_rate=0.025,
    output:
        DIR_bams + "/SSCS1_filtered/{wildcard}.bam",
    threads: 12
    run:
        shell(
            "samtools sort -n {input} | fgbio FilterConsensusReads -Xmx20G \
            --input=/dev/stdin \
            --output={output} \
            --ref={params.PATH_hg38} \
            --min-reads={params.min_reads} \
            --max-read-error-rate={params.max_read_error_rate} \
            --min-base-quality={params.min_base_quality} \
            --reverse-per-base-tags=true \
            --sort-order=Queryname"
        )

rule clip_SSCS:
    input:
        DIR_bams + "/SSCS1_filtered/{wildcard}.bam",
    params:
        PATH_hg38=PATH_hg38,
    output:
        temp(DIR_bams + "/SSCS1_clipped_unsorted/{wildcard}.bam"),
    threads: 12
    run:
        shell(
            "samtools sort -n {input}  | fgbio ClipBam -Xmx20G --ref {params.PATH_hg38} --input /dev/stdin --output {output}  --clip-overlapping-reads=true")

rule sort_clipped:
    input:
        DIR_bams + "/SSCS1_clipped_unsorted/{wildcard}.bam",
    output:
        DIR_bams + "/SSCS1_clipped/{wildcard}.bam",
    threads: 12
    run:
        shell("samtools sort {input} -o {output}")
        
rule index_clipped:
    input:
        DIR_bams + "/SSCS1_clipped/{wildcard}.bam",
    output:
        DIR_bams + "/SSCS1_clipped/{wildcard}.bam.bai",
    threads: 12
    run:
        shell("samtools index {input}")