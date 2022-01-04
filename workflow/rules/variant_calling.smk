rule run_VarScan_snv: 
	input: 
		DIR_mpileup + "/{cohort_wildcard}/{wildcard}.mpileup"
	output: 
		VarScan_snv + "/{cohort_wildcard}/{wildcard}.vcf"
	params:
		min_coverage = 8,
		min_reads = 2, 
		min_avg_base_qual = 20, 
		min_var_freq = 0.001,
		p_value = 0.05
	threads: 12
	shell:
		"java -jar /home/amunzur/VarScan.v2.3.9.jar pileup2snp {input} \
			--min-coverage {params.min_coverage} \
			--min-reads {params.min_reads} \
			--min-avg-qual {params.min_avg_base_qual} \
			--min-var-freq {params.min_var_freq} \
			--p-value {params.p_value} > {output}"

rule run_VarScan_indel: 
	input: 
		DIR_mpileup + "/{cohort_wildcard}/{wildcard}.mpileup"
	output: 
		VarScan_indel + "/{cohort_wildcard}/{wildcard}.vcf"
	params:
		min_coverage = 8,
		min_reads = 2, 
		min_avg_base_qual = 30, 
		min_var_freq = 0.001,
		p_value = 0.05
	threads: 12
	shell:
		"java -jar /home/amunzur/VarScan.v2.3.9.jar pileup2indel {input} \
			--min-coverage {params.min_coverage} \
			--min-reads {params.min_reads} \
			--min-avg-qual {params.min_avg_base_qual} \
			--min-var-freq {params.min_var_freq} \
			--p-value {params.p_value} > {output}"

# Modify the VarScan2 snv output in such a way that ANNOVAR can handle it.
rule make_ANNOVAR_snv_input: 
	input: 
		VarScan_snv + "/{cohort_wildcard}/{wildcard}.vcf"
	output: 
		ANNOVAR_snv_input + "/{cohort_wildcard}/{wildcard}_anno.tsv"
	shell:
		'paste <(cat {input} | cut -f1,2)  <(cat {input} | cut -f2,3,19) > {output}'

# Modify the VarScan2 indel output in such a way that ANNOVAR can handle it.
rule make_ANNOVAR_indel_input: 
	input: 
		VarScan_indel + "/{cohort_wildcard}/{wildcard}.vcf"
	output: 
		ANNOVAR_indel_input + "/{cohort_wildcard}/{wildcard}_anno.tsv"
	conda: 
		"envs/r_env_v2.yaml"
	shell:
		'Rscript --silent --slave /groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/analysis/make_anno_input_indel.R\
			--PATH_VarScan_indel {input} \
			--ANNOVAR_indel_input {output}'

# Annotate the SNVs
rule run_ANNOVAR_snv: 
	input: 
		ANNOVAR_snv_input + "/{cohort_wildcard}/{wildcard}_anno.tsv"
	output: 
		ANNOVAR_snv_output + "/{cohort_wildcard}/{wildcard}.hg38_multianno.txt"
	params: 
		actual_output_file = ANNOVAR_snv_output + "/{cohort_wildcard}/{wildcard}"
	threads: 12
	shell:
		'perl /groups/wyattgrp/software/annovar/annovar/table_annovar.pl {input} /groups/wyattgrp/software/annovar/annovar/humandb/ \
			-buildver hg38 \
			-out {params.actual_output_file} \
			-remove \
			-protocol refGene,knownGene,avsnp147,exac03,cosmic70,clinvar_20170130,kaviar_20150923,gnomad_exome,dbnsfp33a,dbscsnv11,hrcr1,mcap,revel \
			-operation g,g,f,f,f,f,f,f,f,f,f,f,f \
			-nastring .'

# Annotate the indels, replace the ALT field in the ANNOVAR outputs with more information
rule run_ANNOVAR_indel: 
	input: 
		make_ANNOVAR_indel_input = ANNOVAR_indel_input + "/{cohort_wildcard}/{wildcard}_anno.tsv",
		VarScan_indel = VarScan_indel + "/{cohort_wildcard}/{wildcard}.vcf"
	output: 
		ANNOVAR_indel_output = ANNOVAR_indel_output + "/{cohort_wildcard}/{wildcard}.hg38_multianno.txt"
	params: 
		actual_output_file = ANNOVAR_indel_output + "/{cohort_wildcard}/{wildcard}",
		temp_output_file = ANNOVAR_indel_output + "/{cohort_wildcard}/{wildcard}_temp"
	threads: 12
	shell:
		'perl /groups/wyattgrp/software/annovar/annovar/table_annovar.pl {input.make_ANNOVAR_indel_input} /groups/wyattgrp/software/annovar/annovar/humandb/ \
			-buildver hg38 \
			-out {params.actual_output_file} \
			-remove \
			-protocol refGene,knownGene,avsnp147,exac03,cosmic70,clinvar_20170130,kaviar_20150923,gnomad_exome,dbnsfp33a,dbscsnv11,hrcr1,mcap,revel \
			-operation g,g,f,f,f,f,f,f,f,f,f,f,f \
			-nastring . && paste <(cut -f1-4 {output.ANNOVAR_indel_output}) <(cut -f19 {input.VarScan_indel}) <(cut -f6- {output.ANNOVAR_indel_output}) > {params.temp_output_file} && \
			rm {output.ANNOVAR_indel_output} && \
			mv {params.temp_output_file} {output.ANNOVAR_indel_output} && \
			perl -pi -e s/VarAllele/Alt/g {output.ANNOVAR_indel_output}'

rule run_VarDict: 
	input: 
		SC_bam = DIR_bams + "/{cohort_wildcard}/SC_penalty/{wildcard}.bam",
		SC_bam_index = DIR_bams + "/{cohort_wildcard}/SC_penalty/{wildcard}.bam.bai" # helps make sure the bam was indexed before variant calling
	output: 
		DIR_Vardict + "/{cohort_wildcard}/{wildcard}.vcf" # snv and indel together
	params: 
		PATH_hg38 = PATH_hg38, 
		PATH_bed = PATH_bed,
		THRESHOLD_VarFreq = "0.001", 
		sample_name = "{wildcard}"
	threads: 12
	shell:
		"/home/amunzur/VarDictJava/build/install/VarDict/bin/VarDict -G {params.PATH_hg38} -f {params.THRESHOLD_VarFreq} -N {params.sample_name} -b {input.SC_bam} \
		-c 1 -S 2 -E 3 -g 4 {params.PATH_bed} | /home/amunzur/VarDictJava/build/install/VarDict/bin/teststrandbias.R | /home/amunzur/VarDictJava/build/install/VarDict/bin/var2vcf_valid.pl \
		-N {params.sample_name} -E -f {params.THRESHOLD_VarFreq} > {output}"

rule make_ANNOVAR_Vardict_input: 
	input: 
		DIR_Vardict + "/{cohort_wildcard}/{wildcard}.vcf"
	output: 
		ANNOVAR_Vardict_input + "/{cohort_wildcard}/{wildcard}_anno.tsv"
	conda: 
		"envs/r_env_v2.yaml"
	shell:
		'Rscript --silent --slave /groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/analysis/make_anno_input_vardict.R\
			--PATH_Vardict_output {input} \
			--PATH_ANNOVAR_input {output}'

rule run_ANNOVAR_vardict: 
	input: 
		ANNOVAR_Vardict_input + "/{cohort_wildcard}/{wildcard}_anno.tsv"
	output: 
		ANNOVAR_Vardict_output + "/{cohort_wildcard}/{wildcard}.hg38_multianno.txt"
	params: 
		actual_output_file = ANNOVAR_Vardict_output + "/{cohort_wildcard}/{wildcard}"		
	shell:
		'perl /groups/wyattgrp/software/annovar/annovar/table_annovar.pl {input} /groups/wyattgrp/software/annovar/annovar/humandb/ \
			-buildver hg38 \
			-out {params.actual_output_file} \
			-remove \
			-protocol refGene,knownGene,avsnp147,exac03,cosmic70,clinvar_20170130,kaviar_20150923,gnomad_exome,dbnsfp33a,dbscsnv11,hrcr1,mcap,revel \
			-operation g,g,f,f,f,f,f,f,f,f,f,f,f \
			-nastring .'

rule reformat_vardict_results: 
	input: 
		DIR_Vardict + "/{cohort_wildcard}/{wildcard}.vcf"
	output: 
		DIR_Vardict + "/{cohort_wildcard}_reformatted/{wildcard}.tsv"
	conda: 
		"envs/r_env_v2.yaml"
	shell:
		'Rscript --silent --slave /groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/analysis/reformat_vardict.R\
			--PATH_Vardict_output {input} \
			--PATH_Vardict_reformatted {output}'

# rule make_IGV_snapshot_script:
# 	input: 
# 		DIR_Vardict + "/{cohort_wildcard}/{wildcard}.vcf"
# 	output: 
# 		DIR_Vardict + "/{cohort_wildcard}_reformatted/{wildcard}.tsv"
# 	shell:
# 		'Rscript --silent --slave /groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/analysis/reformat_vardict.R\
# 			--PATH_Vardict_output {input} \
# 			--PATH_Vardict_reformatted {output}'
