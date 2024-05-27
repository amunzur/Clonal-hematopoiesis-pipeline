# Generate mapped DCS and SSCS bams from the unmapped consensus bam files.
rule BamtoFastq2:
    input:
        DIR_bams + "/{consensus_type}_uBAM/{wildcard}.bam",
    output:
        temp(DIR_fastq + "/{consensus_type}_consensus_fq_2/{wildcard}.fastq"),
    threads: 12
    run:
        shell(
            "picard SamToFastq -Xmx20G I={input} INTERLEAVE=true F={output}"
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

rule MergeBamAlignment2:
    input:
        mBAM=DIR_bams + "/{consensus_type}_mBAM_raw_2/{wildcard}.bam",
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
                ALIGNER_PROPER_PAIR_FLAGS=true \
                EXPECTED_ORIENTATIONS=FR \
                MAX_INSERTIONS_OR_DELETIONS=-1 \
                VALIDATION_STRINGENCY=SILENT \
                SORT_ORDER=coordinate \
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
        --mad 5000 \
        --targets {input.PATH_bed} \
        --tmpdir /groups/wyattgrp/users/amunzur/COMPOST_BIN > \
        '/groups/wyattgrp/users/amunzur/pipeline/results/logs_slurm/indel_realignment/{wildcards.wildcard}'"
        )


# Good bams we want to keep before doing any downstream filtering on them, so output shouldn't be temp
rule fixmate2:
    input:
        DIR_bams + "/{consensus_type}_abra2/{wildcard}.bam",
    output:
        temp(DIR_bams + "/{consensus_type}_fixmate/{wildcard}.bam"),
    threads: 12
    run:
        shell(
            "picard -Xmx40g FixMateInformation I={input} O={output} VALIDATION_STRINGENCY=SILENT"
        )

rule subset_to_proper_pairs2:
    input:
        DIR_bams + "/{consensus_type}_fixmate/{wildcard}.bam",
    output:
        DIR_bams + "/{consensus_type}_final/{wildcard}.bam",
    threads: 12
    run:
        shell(
            "sambamba sort {input} --sort-picard -F 'proper_pair' -t 12 -o {output}"
        )


rule FilterConsensusReads_SSCS:
    input:
        DIR_bams + "/SSCS_final/{wildcard}.bam",
    params:
        PATH_hg38=PATH_hg38,
        PATH_bed=PATH_bed,
        sample="{wildcard}",
        min_reads=1,
        min_base_quality=30,
        max_read_error_rate=0.025,
        max_no_call_fraction=0.15,
    output:
        DIR_bams + "/SSCS1_filtered/{wildcard}.bam",
    threads: 12
    run:
        shell(
            "fgbio FilterConsensusReads -Xmx20G \
            --input=/dev/stdin \
            --output={output} \
            --ref={params.PATH_hg38} \
            --min-reads={params.min_reads} \
            --max-read-error-rate={params.max_read_error_rate} \
            --max-no-call-fraction={params.max_no_call_fraction} \
            --min-base-quality={params.min_base_quality} \
            --reverse-per-base-tags=true \
            --sort-order=Coordinate"
        )