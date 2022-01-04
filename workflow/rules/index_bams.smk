rule index_sorted_bams: 
	input: 
		DIR_bams + "/{cohort_wildcard}/sorted/{wildcard}.bam"
	output:
		DIR_bams + "/{cohort_wildcard}/sorted/{wildcard}.bam.bai"
	threads: 3
	shell:
		"samtools index {input}"

rule index_readGroup_bams: 
	input: 
		DIR_bams + "/{cohort_wildcard}/readGroup/{wildcard}.bam"
	output:
		DIR_bams + "/{cohort_wildcard}/readGroup/{wildcard}.bam.bai"
	threads: 3
	shell:
		"samtools index {input}"
		
rule index_SC_penalty_bams: 
	input: 
		DIR_bams + "/{cohort_wildcard}/SC_penalty/{wildcard}.bam"
	output:
		DIR_bams + "/{cohort_wildcard}/SC_penalty/{wildcard}.bam.bai"
	threads: 3
	shell:
		"samtools index {input}"
