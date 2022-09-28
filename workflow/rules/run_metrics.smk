# Compute depth at each position
# consensus type should SSCS3 or DCS.
rule run_depth: 
	input: 
		BAM = DIR_bams + "/{consensus_type}_filtered/{wildcard}.bam",
		PATH_bed = PATH_bed 
	output:
		DIR_depth_metrics + "{consensus_type}/{wildcard}.txt"
	threads: 12
	shell:
		"samtools depth -b {input.PATH_bed} {input.BAM} > {output}"

rule run_mpileup: 
	input: 
		BAM = DIR_bams + "/{consensus_type}/{wildcard}.bam", 
		BAM_index = DIR_bams + "/{consensus_type}/{wildcard}.bam.bai",
		PATH_hg38 = PATH_hg38
	output:
		DIR_mpileup + "/{consensus_type}/{wildcard}.mpileup"
	threads: 12
	shell:
		"samtools mpileup -f {input.PATH_hg38} {input.BAM} -o {output}"

# extract duplicate perc from the PICARD markdup bams
rule extract_duplicate_percentage: 
	input: 		
		DIR_markdup_metrics + "/{wildcard}.txt"
	params: 
		sample_name = "{wildcard}"
	output: 
		DIR_markdup_perc_metrics + "/{wildcard}.txt"
	shell: 
		"paste <(echo {params}) <(head {input} | grep Library | cut -f9) > {output}"

rule run_insert_size: 
	input: 
		DIR_bams + "/readGroup/{wildcard}.bam"
	output:
		metrics = DIR_insertsize_metrics + "/{wildcard}.txt",
		figures = DIR_insertsize_figures + "/{wildcard}.pdf"
	threads: 12
	shell:
		"picard CollectInsertSizeMetrics \
			I={input} \
			O={output.metrics} \
			H={output.figures} \
			M=0.5"