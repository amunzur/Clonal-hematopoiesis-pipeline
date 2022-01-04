rule PDF_to_PNG: 
	input: 
		DIR_insertsize_figures + "{cohort_wildcard}/{wildcard}.pdf"
	output:
		DIR_insertsize_figures_PNG + "{cohort_wildcard}/{wildcard}.png"
	threads: 3
	shell:
		'bash /groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/analysis/pdf_to_png.sh {input} {output}'

# Compute depth at each position
rule run_depth: 
	input: 
		BAM = DIR_bams + "/{cohort_wildcard}/SC_penalty/{wildcard}.bam",
		PATH_bed = PATH_bed 
	output:
		DIR_depth_metrics + "/{cohort_wildcard}/{wildcard}.txt"
	threads: 12
	shell:
		"samtools depth -b {input.PATH_bed} {input.BAM} > {output}"

# Compute depth at each position for UMI deduplicated bams
rule run_depth_dedup: 
	input: 
		BAM = DIR_bams + "/{cohort_wildcard}/dedup/{wildcard}.bam",
		PATH_bed = PATH_bed 
	output:
		DIR_depth_metrics_dedup + "/{cohort_wildcard}/{wildcard}.txt"
	threads: 12
	shell:
		"samtools depth -b {input.PATH_bed} {input.BAM} > {output}"

rule run_mpileup: 
	input: 
		BAM = DIR_bams + "/{cohort_wildcard}/SC_penalty/{wildcard}.bam", 
		BAM_index = DIR_bams + "/{cohort_wildcard}/SC_penalty/{wildcard}.bam.bai",
		PATH_hg38 = PATH_hg38
	output:
		DIR_mpileup + "/{cohort_wildcard}/{wildcard}.mpileup"
	threads: 12
	shell:
		"samtools mpileup -f {input.PATH_hg38} {input.BAM} -o {output}"

# extract duplicate perc from the PICARD markdup bams
rule extract_duplicate_percentage: 
	input: 		
		DIR_markdup_metrics + "/{cohort_wildcard}/{wildcard}.txt"
	params: 
		sample_name = "{wildcard}"
	output: 
		DIR_markdup_perc_metrics + "/{cohort_wildcard}/{wildcard}.txt"
	shell: 
		"paste <(echo {params}) <(head {input} | grep Library | cut -f9) > {output}"

# extract duplicate perc from the UMI deduplicated bams
# rule extract_duplicate_percentage_dedup: 
# 	input: 		
# 		DIR_markdup_metrics + "/{cohort_wildcard}/{wildcard}.txt"
# 	params: 
# 		sample_name = "{wildcard}"
# 	output: 
# 		DIR_markdup_perc_metrics + "/{cohort_wildcard}/{wildcard}.txt"
# 	shell: 
# 		"paste <(echo {params}) <(head {input} | grep Library | cut -f9) > {output}"


# run complexity estimates after deduplication
rule PICARD_complexity_dedup:
	input: 
		DEDUP_bam = DIR_bams + "/{cohort_wildcard}/dedup/{wildcard}.bam"
	output: 
		complexity_metrics = DIR_complexity_metrics + "/{cohort_wildcard}/dedup/{wildcard}.txt"
	threads: 12
	shell:
		"picard -Xmx40g EstimateLibraryComplexity \
					I={input} \
					O={output}"

rule run_insert_size: 
	input: 
		DIR_bams + "/{cohort_wildcard}/readGroup/{wildcard}.bam"
	output:
		metrics = DIR_insertsize_metrics + "/{cohort_wildcard}/{wildcard}.txt",
		figures = DIR_insertsize_figures + "/{cohort_wildcard}/{wildcard}.pdf"
	threads: 12
	shell:
		"picard CollectInsertSizeMetrics \
			I={input} \
			O={output.metrics} \
			H={output.figures} \
			M=0.5"
