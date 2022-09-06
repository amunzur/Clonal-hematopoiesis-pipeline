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
		DIR_merged_fastq + "/{cohort_wildcard}/{wildcard}_R1.fastq"
	params:
		sample_name = "{wildcard}"
	output:
		DIR_readcounts_metrics + "/raw/{cohort_wildcard}/{wildcard}.txt"
	shell:
		"paste <(echo {params}) <(echo $(cat {input}|wc -l)/4|bc) > {output}"
		
# mask low quality bases in fastq files
rule mask_fastq:
	input: 
		DIR_merged_fastq + "/{cohort_wildcard}/{wildcard}.fastq"
	output:
		DIR_masked_fastq + "/{cohort_wildcard}/{wildcard}_masked.fastq"
	params: 
		min_base_quality = 20
	run:
		shell("/groups/wyattgrp/users/amunzur/gene_panel_pipeline/dependencies/fasta mask by quality {input} {params.min_base_quality} > {output}")

# extracts from the 5' end
rule extract_UMIs:
	input: 
		R1 = DIR_masked_fastq + "/{cohort_wildcard}/{wildcard}_R1_masked.fastq",
		R2 = DIR_masked_fastq + "/{cohort_wildcard}/{wildcard}_R2_masked.fastq"
	output:
		R1_extracted = DIR_extracted_fastq + "/{cohort_wildcard}/{wildcard}_R1_extracted.fastq",
		R2_extracted = DIR_extracted_fastq + "/{cohort_wildcard}/{wildcard}_R2_extracted.fastq"
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

# Trim fastq files using trim galore and run fastqc on them. R1 and R2 fastq files should be provided at the same time. 
rule trim_fastq: 
	input: 
		R1_extracted = DIR_extracted_fastq + "/{cohort_wildcard}/{wildcard}_R1_extracted.fastq",
		R2_extracted = DIR_extracted_fastq + "/{cohort_wildcard}/{wildcard}_R2_extracted.fastq"
	output:
		pair1_trimmed = DIR_trimmed_fastq + "/{cohort_wildcard}/{wildcard}_R1_extracted_val_1.fq",
		pair2_trimmed = DIR_trimmed_fastq + "/{cohort_wildcard}/{wildcard}_R2_extracted_val_2.fq",
		pair1_trimmed_fastqc_html = DIR_trimmed_fastqc + "/{cohort_wildcard}/{wildcard}_R1_extracted_val_1_fastqc.html",
		pair2_trimmed_fastqc_html = DIR_trimmed_fastqc + "/{cohort_wildcard}/{wildcard}_R2_extracted_val_2_fastqc.html",
		pair1_trimmed_fastqc_zip = DIR_trimmed_fastqc + "/{cohort_wildcard}/{wildcard}_R1_extracted_val_1_fastqc.zip",
		pair2_trimmed_fastqc_zip = DIR_trimmed_fastqc + "/{cohort_wildcard}/{wildcard}_R2_extracted_val_2_fastqc.zip"
	params: 
		min_base_quality = 20, 
		clip_R1 = 0, 
		clip_R2 = 0, 
		three_prime_clip_R1 = 5, 
		three_prime_clip_R2 = 5
	run:
		shell('trim_galore --quality {params.min_base_quality} \
			--phred33 \
			--fastqc \
			--illumina \
			--paired \
			--dont_gzip \
			--three_prime_clip_R1 {params.three_prime_clip_R1} \
			--three_prime_clip_R2 {params.three_prime_clip_R2} \
			--output_dir `dirname {output.pair1_trimmed}` \
			--fastqc_args "--nogroup --outdir `dirname {output.pair1_trimmed_fastqc_html}`" \
			{input.R1_extracted} \
			{input.R2_extracted}')