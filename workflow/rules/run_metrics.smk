rule run_depth: 
	input: 
		BAM = DIR_bams + "/{consensus_type}_filtered/{wildcard}.bam",
		PATH_bed = PATH_bed 
	output:
		DIR_depth_metrics + "/{consensus_type}/{wildcard}.txt"
	threads: 12
	shell:
		"samtools sort {input.BAM} | samtools depth -b {input.PATH_bed} /dev/stdin > {output}"

rule run_mpileup: 
	input: 
		BAM = DIR_bams + "/{consensus_type}_clipped/{wildcard}.bam", 
		BAM_index = DIR_bams + "/{consensus_type}_clipped/{wildcard}.bam.bai",
		PATH_hg38 = PATH_hg38
	output:
		DIR_mpileup + "/temp_{consensus_type}/{wildcard}.mpileup"
	threads: 12
	shell:
		"samtools mpileup -A -f {input.PATH_hg38} {input.BAM} -o {output}"

rule modify_mpileup: 
	input: 
		DIR_mpileup + "/temp_{consensus_type}/{wildcard}.mpileup"
	output:
		DIR_mpileup + "/{consensus_type}/{wildcard}.mpileup"
	shell:
		"paste <(awk '{{print $1$2}}' {input}) <(cat {input}) > {output}"

rule run_insert_size: 
	input: 
		DIR_bams + "/{consensus_type}_filtered/{wildcard}.bam"
	output:
		metrics = DIR_insertsize_metrics + "{consensus_type}/{wildcard}.txt",
		figures = DIR_insertsize_figures + "{consensus_type}/{wildcard}.pdf"
	threads: 12
	shell:
		"picard CollectInsertSizeMetrics \
			I={input} \
			O={output.metrics} \
			H={output.figures} \
			M=0.5"

rule run_read_counts: 
	input: 
		DIR_merged_fastq + "/{wildcard}_R1.fastq"
	params:
		sample_name = "{wildcard}"
	output:
		DIR_readcounts_metrics + "/merged/{wildcard}.txt",
	shell:
		"paste <(echo {params}) <(echo $(cat {input}|wc -l)/4|bc) > {output}"

# rule run_tnvstats:
# 	input:
# 		BAM_tumor_path=DIR_bams + "/SSCS1_clipped/{wildcard}.bam",
# 		BAM_normal_path=
# 		BAM_normal_name=
# 				BAM_index=DIR_bams + "/SSCS1_clipped/{wildcard}.bam.bai",

# 	params:
# 	output:
# 	conda:
#         "../envs/pysamstat_v2.yaml"
# 	shell:
# 		'script_dir="/groups/wyattgrp/users/amunzur/gene_panel_pipeline/scripts";
# 		mkdir -p /groups/wyattgrp/users/amunzur/chip_project/tnvstats/kidney_samples; # output dir to place all tnvstats
# 		configfile="/groups/wyattgrp/users/amunzur/gene_panel_pipeline/scripts/config.txt"; # conda envs and stuff

# 		source ${configfile};
# 		source ${conda_profile_path};
# 		conda activate ${conda_pysam_env};

# 		while read name_tumor path_tumor name_normal path_normal
# 		do
# 		    echo $name_tumor
# 		    echo $name_normal
# 		    outputdir="/groups/wyattgrp/users/amunzur/chip_project/tnvstats/kidney_samples/${name_tumor}";

# 		    # tbam=$(readlink -ve ${path_tumor})
# 		    # nbam=$(readlink -ve ${path_normal})
# 		    printf "bash ${script_dir}/make_tnvstat_file.bash ${path_tumor} ${path_normal} ${name_tumor} ${name_normal} ${outputdir} ${configfile}\n";
# 		    sbatch ${script_dir}/make_tnvstat_file.bash ${path_tumor} ${path_normal} ${name_tumor} ${name_normal} ${outputdir} ${configfile};
# 		done < /groups/wyattgrp/users/amunzur/pipeline/results/metrics/tnvstats/kidney_samples/bam_subset.txt'