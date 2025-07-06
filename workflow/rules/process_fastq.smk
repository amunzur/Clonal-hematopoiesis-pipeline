
rule run_fastqc_merged:
    input:
        DIR_merged_fastq + "/{wildcard}",
    output:
        output_zip=DIR_merged_fastqc + "/{wildcard}_fastqc.zip",
        output_html=DIR_merged_fastqc + "/{wildcard}_fastqc.html",
    threads: 5
    shell:
        "/home/amunzur/FastQC/fastqc {input} --outdir=`dirname {output.output_zip}`"

rule fastq_read_counts:
    input:
        DIR_merged_fastq + "/{wildcard}_1.fq.gz",
    params:
        sample_name="{wildcard}",
    output:
        DIR_readcounts_metrics + "/raw/{wildcard}.txt",
    shell:
        "paste <(echo {params}) <(echo $(cat {input}|wc -l)/4|bc) > {output}"

rule trim_fastq:
    input:
        R1=DIR_merged_fastq + "/{wildcard}_1.fq.gz",
        R2=DIR_merged_fastq + "/{wildcard}_2.fq.gz",
    output:
        R1=temp(DIR_trimmed_fastq + "/{wildcard}_1.fq.gz"),
        R2=temp(DIR_trimmed_fastq + "/{wildcard}_2.fq.gz"),
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
        DIR_trimmed_fastq + "/{wildcard}",
    output:
        output_zip=DIR_trimmed_fastqc + "/{wildcard}_fastqc.zip",
        output_html=DIR_trimmed_fastqc + "/{wildcard}_fastqc.html",
    threads: 5
    shell:
        """
        /home/amunzur/FastQC/fastqc {input} --outdir=`dirname {output}`
        """
