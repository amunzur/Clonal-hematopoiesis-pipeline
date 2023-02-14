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

