rule run_mutect2:
    input:
        bam=DIR_bams + "/{consensus_type}_final/{wildcard}.bam",
        index=DIR_bams + "/{consensus_type}_final/{wildcard}.bam.bai",
    output:
        vcf=temp(DIR_results + "/variant_calling_raw/Mutect2/{consensus_type}/{wildcard}.vcf.gz"),
        stats=DIR_results + "/variant_calling_raw/Mutect2/{consensus_type}/{wildcard}.vcf.gz.stats",
        bamout=DIR_results + "/variant_calling_raw/Mutect2/{consensus_type}_bamout/{wildcard}.bam",
    conda:
        "../envs/chip_variantcalling.yaml"
    params:
        PATH_hg38=PATH_hg38,
        PATH_bed=PATH_bed,
    threads: 12
    shell:
        """
        export PATH=/usr/bin:$PATH
        /home/amunzur/gatk-4.2.0.0/gatk Mutect2 \
            --reference {params.PATH_hg38} \
            --intervals {params.PATH_bed} \
            --input {input.bam} \
            --output {output.vcf} \
            --bamout {output.bamout} \
            --force-active true \
            --initial-tumor-lod 0 \
            --tumor-lod-to-emit 0
        """

rule unzip_mutect:
    input:
        DIR_results + "/variant_calling_raw/Mutect2/{consensus_type}/{wildcard}.vcf.gz"
    output:
        temp(DIR_results + "/variant_calling_raw/Mutect2/{consensus_type}/{wildcard}.vcf")
    shell:
        "gunzip {input}"

rule run_VarDict_chip:
    input:
        DIR_bams + "/{consensus_type}_final/{wildcard}.bam",
    output:
        temp(DIR_results + "/variant_calling_raw/Vardict/{consensus_type}/{wildcard}.vcf"),
    params:
        PATH_hg38=PATH_hg38,
        PATH_bed=PATH_bed,
        THRESHOLD_VarFreq="0.001",
        sample_name="{wildcard}",
        min_variant_reads=4,
    threads: 12
    shell:
        "/home/amunzur/VarDictJava/build/install/VarDict/bin/VarDict \
        -G {params.PATH_hg38} \
        -f {params.THRESHOLD_VarFreq} \
        -N {params.sample_name} \
        -r {params.min_variant_reads} \
        -b {input} \
        -k 0 -c 1 -S 2 -E 3 -g 4 {params.PATH_bed} | \
        /home/amunzur/VarDictJava/build/install/VarDict/bin/teststrandbias.R | \
        /home/amunzur/VarDictJava/build/install/VarDict/bin/var2vcf_valid.pl \
        -f {params.THRESHOLD_VarFreq} > {output}"

rule run_freebayes_chip:
    input:
        DIR_bams + "/{consensus_type}_final/{wildcard}.bam",
    output:
        temp(DIR_results + "/variant_calling_raw/freebayes/{consensus_type}/{wildcard}.vcf"),
    conda:
        "../envs/chip_variantcalling.yaml"
    params:
        PATH_hg38=PATH_hg38,
        PATH_bed=PATH_bed,
    threads: 12
    shell:
        "freebayes {input} \
        -f {params.PATH_hg38} \
        -t {params.PATH_bed} \
        --pooled-continuous \
        --min-alternate-fraction 0.001 \
        --min-alternate-count 1 | bcftools norm -m-both -f {params.PATH_hg38} -o {output}"

rule zip_vcf_files_vardict:
    input:
        DIR_results + "/variant_calling_raw/Vardict/{consensus_type}/{wildcard}.vcf",
    output:
        temp(DIR_results + "/variant_calling_raw/Vardict/{consensus_type}/{wildcard}.vcf.gz"),
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bgzip -c {input} > {output}"

rule zip_vcf_files_freebayes:
    input:
        DIR_results + "/variant_calling_raw/freebayes/{consensus_type}/{wildcard}.vcf",
    output:
        temp(DIR_results + "/variant_calling_raw/freebayes/{consensus_type}/{wildcard}.vcf.gz"),
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bgzip -c {input} > {output}"

rule index_vcf_files:
    input:
        DIR_results + "/variant_calling_raw/{variant_caller}/{consensus_type}/{wildcard}.vcf.gz",
    output:
        temp(DIR_results + "/variant_calling_raw/{variant_caller}/{consensus_type}/{wildcard}.vcf.gz.tbi"),
    conda:
        "../envs/bcftools.yaml"
    shell:
        "tabix -p vcf {input}"

rule sort_vcf_chip:
    input:
        vcf = DIR_results + "/variant_calling_raw/{variant_caller}/{consensus_type}/{wildcard}.vcf.gz",
        index = DIR_results + "/variant_calling_raw/{variant_caller}/{consensus_type}/{wildcard}.vcf.gz.tbi"
    output:
        temp(DIR_results + "/variant_calling_sorted/{variant_caller}/{consensus_type}/{wildcard}.vcf.gz"),
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools sort {input.vcf} -Oz -o {output}"

rule index_sorted_vcf_files:
    input:
        DIR_results + "/variant_calling_sorted/{variant_caller}/{consensus_type}/{wildcard}.vcf.gz",
    output:
        temp(DIR_results + "/variant_calling_sorted/{variant_caller}/{consensus_type}/{wildcard}.vcf.gz.tbi"),
    conda:
        "../envs/bcftools.yaml"
    shell:
        "tabix -p vcf {input}"

rule normalize_variants:
    input:
        vcf=DIR_results + "/variant_calling_sorted/{variant_caller}/{consensus_type}/{wildcard}.vcf.gz",
        index=DIR_results + "/variant_calling_sorted/{variant_caller}/{consensus_type}/{wildcard}.vcf.gz.tbi"
    output:
        temp(DIR_results + "/variant_calling_normalized/{variant_caller}/{consensus_type}/{wildcard}.vcf.gz"),
    params:
        PATH_hg38_dict="/groups/wyattgrp/reference/hg38/hg38.fa",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools norm \
            {input.vcf} \
            -m-both \
            -f {params} \
            -o {output}"

rule decompose_blocksubstitutions_chip:
    input:
        DIR_results + "/variant_calling_normalized/{variant_caller}/{consensus_type}/{wildcard}.vcf.gz",
    output:
        DIR_results + "/variant_calling_chip/{variant_caller}/{consensus_type}/{wildcard}.vcf.gz",
    conda:
        "../envs/chip_variantcalling.yaml"
    shell:
        "vt decompose_blocksub {input} -o {output}"