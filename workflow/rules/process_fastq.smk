
rule run_fastqc_merged:
    input:
        DIR_fastq + "/merged/{wildcard}",
    output:
        output_zip= DIR_fastq + "/merged/{wildcard}_fastqc.zip",
        output_html= DIR_fastq + "/merged/{wildcard}_fastqc.html",
    threads: 5
    shell:
        "/home/amunzur/FastQC/fastqc {input} --outdir=`dirname {output.output_zip}`"

rule mask_fastq:
    input:
        DIR_fastq + "/merged/{wildcard}.fq.gz",
    output:
        temp(DIR_fastq + "/masked/{wildcard}.fq.gz"),
    params:
        min_base_quality=20,
    run:
        shell(
            "zcat {input} | /groups/wyattgrp/users/amunzur/gene_panel_pipeline/dependencies/fasta mask by quality - {params.min_base_quality} | gzip > {output}"
        )

rule trim_fastq:
    input:
        R1= DIR_fastq + "/masked/{wildcard}_1.fq.gz",
        R2= DIR_fastq + "/masked/{wildcard}_2.fq.gz",
    output:
        R1=temp(DIR_fastq + "/trimmed/{wildcard}_1.fq.gz"),
        R2=temp(DIR_fastq + "/trimmed/{wildcard}_2.fq.gz"),
        html_report="results/reports/fastp/{wildcard}.html",
        json_report="results/reports/fastp/{wildcard}.json",
    threads: 12
    params:
        minimum_read_length=32,
    run:
        shell(
            "fastp -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} \
        --dont_eval_duplication \
        --cut_tail \
        --length_required {params.minimum_read_length} \
        --html {output.html_report} \
        --json {output.json_report} \
        --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        --thread {threads}"
        )

rule run_fastqc_trimmed:
    input:
        DIR_fastq + "/trimmed/{wildcard}",
    output:
        output_zip="results/reports/fastqc/trimmed/{wildcard}_fastqc.zip",
        output_html="results/reports/fastqc/trimmed/{wildcard}_fastqc.html",
    threads: 5
    shell:
        """
        /home/amunzur/FastQC/fastqc {input} --outdir=`dirname {output}`
        """