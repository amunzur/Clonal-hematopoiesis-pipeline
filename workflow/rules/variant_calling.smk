rule run_VarScan_snv:
    input:
        DIR_mpileup + "/{consensus_type}/{wildcard}.mpileup",
    output:
        VarScan_snv + "/{consensus_type}/{wildcard}.vcf",
    params:
        min_coverage=8,
        min_reads=2,
        min_avg_base_qual=20,
        min_var_freq=0.001,
        p_value=0.05,
    threads: 12
    shell:
        "java -jar /home/amunzur/VarScan.v2.3.9.jar pileup2snp {input} \
        --min-coverage {params.min_coverage} \
        --min-reads {params.min_reads} \
        --min-avg-qual {params.min_avg_base_qual} \
        --min-var-freq {params.min_var_freq} \
        --p-value {params.p_value} > {output}"


rule run_VarScan_indel:
    input:
        DIR_mpileup + "/{consensus_type}/{wildcard}.mpileup",
    output:
        VarScan_indel + "/{consensus_type}/{wildcard}.vcf",
    params:
        min_coverage=8,
        min_reads=2,
        min_avg_base_qual=30,
        min_var_freq=0.001,
        p_value=0.05,
    threads: 12
    shell:
        "java -jar /home/amunzur/VarScan.v2.3.9.jar pileup2indel {input} \
        --min-coverage {params.min_coverage} \
        --min-reads {params.min_reads} \
        --min-avg-qual {params.min_avg_base_qual} \
        --min-var-freq {params.min_var_freq} \
        --p-value {params.p_value} > {output}"


rule run_VarDict:
    input:
        SC_bam=DIR_bams + "/{consensus_type}_clipped/{wildcard}.bam",
        INDEX = DIR_bams + "/{consensus_type}_clipped/{wildcard}.bam.bai",
    output:
        DIR_Vardict + "/{consensus_type}/{wildcard}.vcf",
    params:
        PATH_hg38=PATH_hg38,
        PATH_bed=PATH_bed,
        THRESHOLD_VarFreq="0.001",
        sample_name="{wildcard}",
    threads: 12
    shell:
        "/home/amunzur/VarDictJava/build/install/VarDict/bin/VarDict -G {params.PATH_hg38} -f {params.THRESHOLD_VarFreq} -u -N {params.sample_name} -b {input.SC_bam} \
        -k 1 -c 1 -S 2 -E 3 -g 4 {params.PATH_bed} | /home/amunzur/VarDictJava/build/install/VarDict/bin/teststrandbias.R | /home/amunzur/VarDictJava/build/install/VarDict/bin/var2vcf_valid.pl \
        -N {params.sample_name} -E -f {params.THRESHOLD_VarFreq} > {output}"

rule run_Mutect:
    input:
        BAM = DIR_bams + "/{consensus_type}_clipped/{wildcard}.bam",
        INDEX = DIR_bams + "/{consensus_type}_clipped/{wildcard}.bam.bai",
    output:
        vcf=DIR_Mutect + "/{consensus_type}/raw/{wildcard}_vcf.gz",
        stats=DIR_Mutect + "/{consensus_type}/raw/{wildcard}_vcf.gz.stats",
        index=DIR_Mutect + "/{consensus_type}/raw/{wildcard}_vcf.gz.tbi",
    params:
        PATH_hg38=PATH_hg38,
        PATH_bed=PATH_bed,
        PATH_PoN=PATH_PoN,
    threads: 12
    shell:
        "/home/amunzur/gatk-4.2.0.0/gatk Mutect2 \
        -R {params.PATH_hg38} \
        -I {input.BAM} \
        -O {output.vcf} \
        -max-mnp-distance 0 \
        --panel-of-normals {params.PATH_PoN} \
        --f1r2-median-mq 10 \
        --f1r2-min-bq 10 \
        --read-filter AllowAllReadsReadFilter \
        --intervals {params.PATH_bed}"

rule filter_Mutect:
    input:
        vcf=DIR_Mutect + "/{consensus_type}/raw/{wildcard}_vcf.gz",
        PATH_hg38=PATH_hg38,
    output:
        DIR_Mutect + "/{consensus_type}/filtered/{wildcard}_vcf.gz",
    threads: 12
    shell:
        "/home/amunzur/gatk-4.2.0.0/gatk FilterMutectCalls \
        -R {input.PATH_hg38} \
        -V {input.vcf} \
        -O {output} \
        --min-median-base-quality 10 \
        --min-median-mapping-quality 10 \
        --min-slippage-length 15 \
        --max-events-in-region 4 \
        --verbosity INFO"

rule unzip_mutect:
    input:
        DIR_Mutect + "/{consensus_type}/filtered/{wildcard}_vcf.gz"
    output:
        DIR_Mutect + "/{consensus_type}/filtered/{wildcard}_vcf"
    shell:
        "gunzip {input}"

# Consensus_type below refers to the consensus for the cfDNA sample
rule GetBaseCounts_Vardict_SSCS1: 
    input: 
        gDNA_vcf=DIR_Vardict + "/SSCS1/{wildcard}.vcf",
        cfDNA_bam=lambda wildcards: get_cfDNA_bam_SSCS1("{wildcard}".format(wildcard=wildcards.wildcard)),
        PATH_hg38=PATH_hg38,
    output: 
        "results/variant_calling/base_counts/Vardict/gDNA_SSCS1_cfDNA_SSCS1" + "/{wildcard}.vcf"
    run: 
        shell("~/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta {input.PATH_hg38} --bam cfDNA_bam:{input.cfDNA_bam} --vcf {input.gDNA_vcf} --output {output}")

rule GetBaseCounts_Vardict_DCS: 
    input: 
        gDNA_vcf=DIR_Vardict + "/SSCS1/{wildcard}.vcf",
        cfDNA_bam=lambda wildcards: get_cfDNA_bam_DCS("{wildcard}".format(wildcard=wildcards.wildcard)),
        PATH_hg38=PATH_hg38,
    output: 
        "results/variant_calling/base_counts/Vardict/gDNA_SSCS1_cfDNA_DCS" + "/{wildcard}.vcf"
    run: 
        shell("~/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta {input.PATH_hg38} --bam cfDNA_bam:{input.cfDNA_bam} --vcf {input.gDNA_vcf} --output {output}")


rule GetBaseCounts_Mutect_SSCS1: 
    input: 
        gDNA_vcf=DIR_Mutect + "/SSCS1/filtered/{wildcard}_vcf",
        cfDNA_bam=lambda wildcards: get_cfDNA_bam_SSCS1("{wildcard}".format(wildcard=wildcards.wildcard)),
        PATH_hg38=PATH_hg38,
    output: 
        "results/variant_calling/base_counts/Mutect2/gDNA_SSCS1_cfDNA_SSCS1" + "/{wildcard}.vcf"
    run: 
        shell("~/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta {input.PATH_hg38} --bam cfDNA_bam:{input.cfDNA_bam} --vcf {input.gDNA_vcf} --output {output}")

rule GetBaseCounts_Mutect_DCS: 
    input: 
        gDNA_vcf=DIR_Mutect + "/SSCS1/filtered/{wildcard}_vcf",
        cfDNA_bam=lambda wildcards: get_cfDNA_bam_DCS("{wildcard}".format(wildcard=wildcards.wildcard)),
        PATH_hg38=PATH_hg38,
    output: 
        "results/variant_calling/base_counts/Mutect2/gDNA_SSCS1_cfDNA_DCS" + "/{wildcard}.vcf"
    run: 
        shell("~/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta {input.PATH_hg38} --bam cfDNA_bam:{input.cfDNA_bam} --vcf {input.gDNA_vcf} --output {output}")