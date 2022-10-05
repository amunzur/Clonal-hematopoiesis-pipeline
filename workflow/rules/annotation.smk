##########################
# VARSCAN
##########################
rule make_ANNOVAR_snv_input:
    input:
        VarScan_snv + "/{consensus_type}/{wildcard}.vcf",
    output:
        ANNOVAR_snv_input + "/{consensus_type}/{wildcard}_anno.tsv",
    shell:
        "paste <(cat {input} | cut -f1,2)  <(cat {input} | cut -f2,3,19) > {output}"

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

##########################
# VARDICT
##########################

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

##########################
# MUTECT2
##########################
# Runs directly on the vcf file, no need for manual conversion
rule run_ANNOVAR_mutect2:
    input:
        DIR_Mutect + "/{consensus_type}_filtered/{wildcard}_vcf.gz",
    output:
        ANNOVAR_mutect_output + "/{consensus_type}/{wildcard}.hg38_multianno.txt",
    params:
        actual_output_file=ANNOVAR_mutect_output + "/{consensus_type}/{wildcard}",
    shell:
        "perl /groups/wyattgrp/software/annovar/annovar/table_annovar.pl {input} /groups/wyattgrp/software/annovar/annovar/humandb/ \
        -vcfinput \
        -buildver hg38 \
        -out {params.actual_output_file} \
        -remove \
        -protocol refGene,knownGene,avsnp147,exac03,cosmic70,clinvar_20170130,kaviar_20150923,gnomad_exome,dbnsfp33a,dbscsnv11,hrcr1,mcap,revel \
        -operation g,g,f,f,f,f,f,f,f,f,f,f,f \
        -nastring ."