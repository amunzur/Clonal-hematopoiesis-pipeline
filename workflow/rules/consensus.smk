rule FastqToBam:
    input:
        R1=DIR_trimmed_fastq + "/{wildcard}_R1.fastq",
        R2=DIR_trimmed_fastq + "/{wildcard}_R2.fastq",
    output:
        DIR_bams + "/uBAM/{wildcard}.bam",
    params:
        sample="{wildcard}",
    shell:
        "fgbio FastqToBam -Xmx20G \
        --input {input.R1} {input.R2} \
        --output {output} \
        --read-structure 3M2S+T 3M2S+T \
        --umi-tag RX \
        --sample {params.sample} \
        --library lib1 \
        --platform illumina \
        --sort true"


# Generate mapped bam
rule mapBAM:
	input:
		DIR_bams + "/uBAM/{wildcard}.bam"
	params:
		PATH_hg38 = PATH_hg38,
	output:
		DIR_bams + "/mBAM/{wildcard}.bam"
	threads: 12
	shell:
		"picard SamToFastq -Xmx20G I={input} INCLUDE_NON_PF_READS=true INCLUDE_NON_PRIMARY_ALIGNMENTS=true F=/dev/stdout INTERLEAVE=true | \
        bwa mem {params.PATH_hg38} /dev/stdin -p -P -Y -t {threads} | \
        picard AddOrReplaceReadGroups -Xmx20G I=/dev/stdin O=/dev/stdout RGID=A RGSM=$sample RGPL=illumina RGLB=lib1 RGPU=unit1 SORT_ORDER=queryname | \
        picard MergeBamAlignment -Xmx20G UNMAPPED={input} ALIGNED=/dev/stdin O={output} R={params.PATH_hg38} CLIP_OVERLAPPING_READS=false CLIP_ADAPTERS=false ALIGNER_PROPER_PAIR_FLAGS=true MAX_INSERTIONS_OR_DELETIONS=-1 CREATE_INDEX=true"

