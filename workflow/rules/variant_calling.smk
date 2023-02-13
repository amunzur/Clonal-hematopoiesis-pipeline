

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