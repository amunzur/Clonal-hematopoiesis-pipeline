# align, sort, remove dups from trimmed fastq files. 
rule align_sort:
	input:
		pair1 = DIR_trimmed_fastq + "/{wildcard}_R1_extracted_val_1.fq",
		pair2 = DIR_trimmed_fastq + "/{wildcard}_R2_extracted_val_2.fq",
		PATH_hg38 = PATH_hg38,
		PATH_bed = PATH_bed
	params: 
		min_mapping_quality = 20,
		bitwise_flag = 12, # remove if the read and the mate is unmapped
	output:
		DIR_bams + "/sorted/{wildcard}.bam"
	threads: 12
	shell:
		"bwa mem -t {threads} {input.PATH_hg38} {input.pair1} {input.pair2} | \
		samtools view -h -q {params.min_mapping_quality} -F {params.bitwise_flag} - | \
		samtools sort -o {output}"

rule PICARD_fixmate:
	input: 
		SORTED_bam = DIR_bams + "/sorted/{wildcard}.bam",
		SORTED_bam_index = DIR_bams + "/sorted/{wildcard}.bam.bai" # to make sure this doesnt run before we index the aligned bams
	output: 
		DIR_bams + "/fixmate/{wildcard}.bam"
	threads: 12
	shell:
		"picard -Xmx40g FixMateInformation \
					I={input.SORTED_bam} \
					O={output}"

# Mark and remove duplicates after aligning and sorting. 
rule mark_duplicates_PICARD: 
	input:
		DIR_bams + "/fixmate/{wildcard}.bam"
	output:
		MARKDUP_bam = DIR_bams + "/markdup/{wildcard}.bam",
		MARKDUP_metrics = DIR_markdup_metrics + "/{wildcard}.txt"
	threads: 12
	shell:
		"picard -Xmx40g MarkDuplicates I={input} O={output.MARKDUP_bam} M={output.MARKDUP_metrics}"

# Deduplication using UMI tools
rule dedup_UMITOOLS: 
	input:
		SORTED_bam = DIR_bams + "/sorted/{wildcard}.bam",
		SORTED_bam_idx = DIR_bams + "/sorted/{wildcard}.bam.bai"
	params:
		DEDUP_metrics_general = DIR_dedup_metrics + "/{wildcard}"
	output:
		DEDUP_bam = DIR_bams + "/dedup/{wildcard}.bam",
		DEDUP_metrics1 = DIR_dedup_metrics + "/{wildcard}_edit_distance.tsv",
		DEDUP_metrics2 = DIR_dedup_metrics + "/{wildcard}_per_umi_per_position.tsv",
		DEDUP_metrics3 = DIR_dedup_metrics + "/{wildcard}_per_umi.tsv"
	threads: 12
	conda: 
		"../envs/umi_tools.yaml" # two dots here because it starts in the workflow/rules directory. To go to envs we need to jump up one directory first.
	shell:
		"umi_tools dedup -I {input} --output-stats={params.DEDUP_metrics_general} -S {output.DEDUP_bam}"

# runs matti's sam_mark_duplicates tool on the sorted bams
rule dedup_sam_mark_duplicates: 
	input:
		DIR_bams + "/sorted/{wildcard}.bam"
	output: 
		DIR_bams + "/sam_mark_duplicates/{wildcard}.bam"
	threads: 12
	shell: 
		"sam_mark_duplicates {input} > {output}"

# Add read groups after removing duplicates
rule add_read_groups_PICARD: 
	input:
		DIR_bams + "/markdup/{wildcard}.bam"	
	params: 
		rglb = "library",
		rgpl = "ILLUMINA",
		rgpu = "unit",
		rgsm = "sample"
	output:
		DIR_bams + "/readGroup/{wildcard}.bam"
	threads: 12
	shell:
		"picard -Xmx40g AddOrReplaceReadGroups I={input} O={output} RGID=1 RGLB={params.rglb} RGPL={params.rgpl} RGPU={params.rgpu} RGSM={params.rgsm}"

rule filter_clipped_reads:
	input: 
		BAM = DIR_bams + "/readGroup/{wildcard}.bam",
		PATH_hg38 = PATH_hg38
	output:
		BAM = DIR_bams + "/SC_penalty/{wildcard}.bam",
		# BAM_index = DIR_bams + "/SC_penalty/{wildcard}.bam.bai"
	threads: 12
	params:
		clipping_threshold = 10
	shell:
		"samtools view -h {input.BAM} | /home/amunzur/samclip --ref {input.PATH_hg38} --max {params.clipping_threshold} | samtools view -bh > {output.BAM}"
		# "samtools index {output.BAM}"
