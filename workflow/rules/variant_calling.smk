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


# Modify the VarScan2 snv output in such a way that ANNOVAR can handle it.
rule make_ANNOVAR_snv_input:
    input:
        VarScan_snv + "/{consensus_type}/{wildcard}.vcf",
    output:
        ANNOVAR_snv_input + "/{consensus_type}/{wildcard}_anno.tsv",
    shell:
        "paste <(cat {input} | cut -f1,2)  <(cat {input} | cut -f2,3,19) > {output}"


# Modify the VarScan2 indel output in such a way that ANNOVAR can handle it.
rule make_ANNOVAR_indel_input:
    input:
        VarScan_indel + "/{consensus_type}/{wildcard}.vcf",
    output:
        ANNOVAR_indel_input + "/{consensus_type}/{wildcard}_anno.tsv",
    conda:
        "../envs/r_env_v2.yaml"
    shell:
        "Rscript --silent --slave /groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/analysis/make_anno_input_indel.R\
        --PATH_VarScan_indel {input} \
        --ANNOVAR_indel_input {output}"


# Annotate the SNVs
rule run_ANNOVAR_snv:
    input:
        ANNOVAR_snv_input + "/{consensus_type}/{wildcard}_anno.tsv",
    output:
        ANNOVAR_snv_output + "/{consensus_type}/{wildcard}.hg38_multianno.txt",
    params:
        actual_output_file=ANNOVAR_snv_output + "/{wildcard}",
    threads: 12
    shell:
        "perl /groups/wyattgrp/software/annovar/annovar/table_annovar.pl {input} /groups/wyattgrp/software/annovar/annovar/humandb/ \
        -buildver hg38 \
        -out {params.actual_output_file} \
        -remove \
        -protocol refGene,knownGene,avsnp147,exac03,cosmic70,clinvar_20170130,kaviar_20150923,gnomad_exome,dbnsfp33a,dbscsnv11,hrcr1,mcap,revel \
        -operation g,g,f,f,f,f,f,f,f,f,f,f,f \
        -nastring ."


# Annotate the indels
rule run_ANNOVAR_indel:
    input:
        make_ANNOVAR_indel_input=ANNOVAR_indel_input
        + "/{consensus_type}/{wildcard}_anno.tsv",
        VarScan_indel=VarScan_indel + "/{consensus_type}/{wildcard}.vcf",
    output:
        ANNOVAR_indel_output=ANNOVAR_indel_output
        + "/{consensus_type}/{wildcard}.hg38_multianno.txt",
    params:
        actual_output_file=ANNOVAR_indel_output + "/{consensus_type}/{wildcard}",
        temp_output_file=ANNOVAR_indel_output + "/{consensus_type}/{wildcard}_temp",
    threads: 12
    shell:
        "perl /groups/wyattgrp/software/annovar/annovar/table_annovar.pl {input.make_ANNOVAR_indel_input} /groups/wyattgrp/software/annovar/annovar/humandb/ \
        -buildver hg38 \
        -out {params.actual_output_file} \
        -remove \
        -protocol refGene,knownGene,avsnp147,exac03,cosmic70,clinvar_20170130,kaviar_20150923,gnomad_exome,dbnsfp33a,dbscsnv11,hrcr1,mcap,revel \
        -operation g,g,f,f,f,f,f,f,f,f,f,f,f \
        -nastring ."


rule run_VarDict:
    input:
        SC_bam=DIR_bams + "/{consensus_type}_filtered/{wildcard}.bam",
        SC_bam_index=DIR_bams + "/{consensus_type}_filtered/{wildcard}.bam.bai",  # helps make sure the bam was indexed before variant calling
    output:
        DIR_Vardict + "/{consensus_type}/{wildcard}.vcf",
    params:
        PATH_hg38=PATH_hg38,
        PATH_bed=PATH_bed,
        THRESHOLD_VarFreq="0.001",
        sample_name="{wildcard}",
    threads: 12
    shell:
        "/home/amunzur/VarDictJava/build/install/VarDict/bin/VarDict -G {params.PATH_hg38} -f {params.THRESHOLD_VarFreq} -N {params.sample_name} -b {input.SC_bam} \
        -k 0 -c 1 -S 2 -E 3 -g 4 {params.PATH_bed} | /home/amunzur/VarDictJava/build/install/VarDict/bin/teststrandbias.R | /home/amunzur/VarDictJava/build/install/VarDict/bin/var2vcf_valid.pl \
        -N {params.sample_name} -E -f {params.THRESHOLD_VarFreq} > {output}"


rule make_ANNOVAR_Vardict_input:
    input:
        DIR_Vardict + "/{consensus_type}/{wildcard}.vcf",
    output:
        ANNOVAR_Vardict_input + "/{consensus_type}/{wildcard}_anno.tsv",
    conda:
        "../envs/r_env_v2.yaml"
    shell:
        "Rscript --silent --slave /groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/analysis/make_anno_input_vardict.R\
        --PATH_Vardict_output {input} \
        --PATH_ANNOVAR_input {output}"


# Runs directly on the vcf file, no need for manual conversion
rule run_ANNOVAR_vardict:
    input:
        DIR_Vardict + "/{consensus_type}/{wildcard}.vcf",
    output:
        ANNOVAR_Vardict_output + "/{consensus_type}/{wildcard}.hg38_multianno.txt",
    params:
        actual_output_file=ANNOVAR_Vardict_output + "/{consensus_type}/{wildcard}",
    shell:
        "perl /groups/wyattgrp/software/annovar/annovar/table_annovar.pl {input} /groups/wyattgrp/software/annovar/annovar/humandb/ \
        -vcfinput \
        -buildver hg38 \
        -out {params.actual_output_file} \
        -remove \
        -protocol refGene,knownGene,avsnp147,exac03,cosmic70,clinvar_20170130,kaviar_20150923,gnomad_exome,dbnsfp33a,dbscsnv11,hrcr1,mcap,revel \
        -operation g,g,f,f,f,f,f,f,f,f,f,f,f \
        -nastring ."


rule run_Mutect:
    input:
        DIR_bams + "/{consensus_type}_filtered/{wildcard}.bam",
    output:
        vcf=DIR_Mutect + "/{consensus_type}/raw/{wildcard}_vcf.gz",
        stats=DIR_Mutect + "/{consensus_type}/raw/{wildcard}_vcf.gz.stats",
        index=DIR_Mutect + "/{consensus_type}/raw/{wildcard}_vcf.gz.tbi",
    params:
        PATH_hg38=PATH_hg38,
    threads: 12
    shell:
        "/home/amunzur/gatk-4.2.0.0/gatk Mutect2 \
        -R {params.PATH_hg38} \
        -I {input} \
        -O {output.vcf} \
        -max-mnp-distance 0 \
        -max-af 0.00001 \
        --f1r2-median-mq 10 \
        --f1r2-min-bq 10 \
        --read-filter AllowAllReadsReadFilter"


rule filter_Mutect:
    input:
        vcf=DIR_Mutect + "/{consensus_type}/raw/{wildcard}_vcf.gz",
        PATH_hg38=PATH_hg38,
    output:
        DIR_Mutect + "/{consensus_type}_filtered/{wildcard}_vcf.gz",
    threads: 12
    shell:
        "/home/amunzur/gatk-4.2.0.0/gatk FilterMutectCalls \
        -R {input.PATH_hg38} \
        -V {input.vcf} \
        -O {output} \
        --min-median-base-quality 10 \
        --min-median-mapping-quality 10 \
        --min-slippage-length 15 \
        --verbosity INFO"

# Given a vcf from a cfDNA vcf file, get read counts from the WBC bams in the positions in the vcf file
rule GetBaseCounts_SSCS3: 
    input: 
        cfDNA_vcf=DIR_Vardict + "/SSCS3/{wildcard}.vcf",
        WBC_bam=lambda wildcards: get_WBC_bam_SSCS3("{wildcard}".format(wildcard=wildcards.wildcard)),
        PATH_hg38=PATH_hg38,
    wildcard_constraints:
        wildcard="/cfDNA/"
    output: 
        "results/variant_calling/base_counts/SSCS3" + "/{wildcard}" + ".vcf"
    run: 
        shell("~/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta {input.PATH_hg38} --bam {input.WBC_bam} --vcf {input.cfDNA_vcf} --output {output}")



        
