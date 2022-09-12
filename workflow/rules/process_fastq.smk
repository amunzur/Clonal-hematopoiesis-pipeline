rule run_fastqc_merged:
	input:
		DIR_merged_fastq + "/{wildcard}.fastq"
	output:
		output_zip = DIR_merged_fastqc + "/{wildcard}_fastqc.zip",
		output_html = DIR_merged_fastqc + "/{wildcard}_fastqc.html"
	threads: 5
	shell:
		"/home/amunzur/FastQC/fastqc {input} --outdir=`dirname {output.output_zip}`"

rule fastq_read_counts:
	input:
		DIR_merged_fastq + "/{wildcard}_R1.fastq"
	params:
		sample_name = "{wildcard}"
	output:
		DIR_readcounts_metrics + "/raw/{wildcard}.txt"
	shell:
		"paste <(echo {params}) <(echo $(cat {input}|wc -l)/4|bc) > {output}"
		
# mask low quality bases in fastq files
rule mask_fastq:
	input: 
		DIR_merged_fastq + "/{wildcard}.fastq"
	output:
		DIR_masked_fastq + "/{wildcard}_masked.fastq"
	params: 
		min_base_quality = 20
	run:
		shell("/groups/wyattgrp/users/amunzur/gene_panel_pipeline/dependencies/fasta mask by quality {input} {params.min_base_quality} > {output}")

# extracts from the 5' end
rule extract_UMIs:
	input: 
		R1 = DIR_masked_fastq + "/{wildcard}_R1_masked.fastq",
		R2 = DIR_masked_fastq + "/{wildcard}_R2_masked.fastq"
	output:
		R1_extracted = DIR_extracted_fastq + "/{wildcard}_R1_extracted.fastq",
		R2_extracted = DIR_extracted_fastq + "/{wildcard}_R2_extracted.fastq"
	conda: 
		"../envs/umi_tools.yaml" # two dots here because it starts in the workflow/rules directory. To go to envs we need to jump up one directory first.
	shell:
		"umi_tools extract \
							--ignore-read-pair-suffixes \
							--bc-pattern=NNNCC \
							--bc-pattern2=NNNCC \
		                    --extract-method=string \
		                    --stdin={input.R1} \
		                    --stdout={output.R1_extracted} \
		                    --read2-in={input.R2} \
		                    --read2-out={output.R2_extracted}"

rule trim_fastq: 
	input: 
		R1_extracted = DIR_extracted_fastq + "/{wildcard}_R1_extracted.fastq",
		R2_extracted = DIR_extracted_fastq + "/{wildcard}_R2_extracted.fastq",
		PATH_adapters = PATH_adapters
	output:
		read1_P = DIR_trimmed_fastq + "/{wildcard}_R1_P.fq",
		read1_U = DIR_trimmed_fastq + "/{wildcard}_R1_U.fq",
		read2_P = DIR_trimmed_fastq + "/{wildcard}_R2_P.fq",
		read2_U = DIR_trimmed_fastq + "/{wildcard}_R2_U.fq"
	threads: 12
	params: 
		sliding_window_size = 4,
		sliding_window_quality = 20,
		minimum_read_length = 25,
		five_prime_clipping = 5
	run:
		shell('trimmomatic PE -threads {threads} \
			{input.R1_extracted} \
			{input.R2_extracted} \
			{output.read1_P} \ 
			{output.read1_U} \ 
			{output.read2_P} \ 
			{output.read2_U} \
			SLIDINGWINDOW:{params.sliding_window_size}:{params.sliding_window_quality} \
			MINLEN:{params.minimum_read_length} \
			ILLUMINACLIP:{input.PATH_adapters}:2:40:15 \
			HEADCROP:{PARAMS.five_prime_clipping}' )
