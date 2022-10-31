rule run_fastqc_merged:
    input:
        DIR_merged_fastq + "/{wildcard}.fastq",
    output:
        output_zip=DIR_merged_fastqc + "/{wildcard}_fastqc.zip",
        output_html=DIR_merged_fastqc + "/{wildcard}_fastqc.html",
    threads: 5
    shell:
        "/home/amunzur/FastQC/fastqc {input} --outdir=`dirname {output.output_zip}`"


rule fastq_read_counts:
    input:
        DIR_merged_fastq + "/{wildcard}_R1.fastq",
    params:
        sample_name="{wildcard}",
    output:
        DIR_readcounts_metrics + "/raw/{wildcard}.txt",
    shell:
        "paste <(echo {params}) <(echo $(cat {input}|wc -l)/4|bc) > {output}"


# mask low quality bases in fastq files
rule mask_fastq:
    input:
        DIR_merged_fastq + "/{wildcard}.fastq",
    output:
        DIR_masked_fastq + "/{wildcard}_masked.fastq",
    params:
        min_base_quality=20,
    run:
        shell(
            "/groups/wyattgrp/users/amunzur/gene_panel_pipeline/dependencies/fasta mask by quality {input} {params.min_base_quality} > {output}"
        )


rule trim_fastq:
    input:
        R1=DIR_masked_fastq + "/{wildcard}_R1_masked.fastq",
        R2=DIR_masked_fastq + "/{wildcard}_R2_masked.fastq",
        PATH_adapters=PATH_adapters,
    output:
        R1=DIR_trimmed_fastq + "/{wildcard}_R1.fastq",
        R2=DIR_trimmed_fastq + "/{wildcard}_R2.fastq",
        html_report="results/reports/fastp/{wildcard}.html",
        json_report="results/reports/fastp/{wildcard}.json",
    threads: 12
    params:
        minimum_read_length=50,
    run:
        shell(
            "fastp -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} \
        --length_required {params.minimum_read_length} \
        --correction \
        --length_required 32 \
        --html {output.html_report} \
        --json {output.json_report} \
        --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        --thread {threads}"
        )


rule run_fastqc_trimmed:
    input:
        R1=DIR_trimmed_fastq + "/{wildcard}_R1.fastq",
        R2=DIR_trimmed_fastq + "/{wildcard}_R2.fastq",
    output:
        output_zip_R1=DIR_trimmed_fastqc + "/{wildcard}_R1_fastqc.zip",
        output_html_R1=DIR_trimmed_fastqc + "/{wildcard}_R1_fastqc.html",
        output_zip_R2=DIR_trimmed_fastqc + "/{wildcard}_R2_fastqc.zip",
        output_html_R2=DIR_trimmed_fastqc + "/{wildcard}_R2_fastqc.html",
    threads: 5
    shell:
        """
        /home/amunzur/FastQC/fastqc {input.R1} --outdir=`dirname {output.output_zip_R1}`
        /home/amunzur/FastQC/fastqc {input.R2} --outdir=`dirname {output.output_zip_R2}` 
        """
