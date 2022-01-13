rule Haplotype_Caller:
	input:
		SC_bam = DIR_bams + "/{cohort_wildcard}/SC_penalty/{wildcard}.bam",
		PATH_hg38 = PATH_hg38,
		PATH_bed = PATH_bed
	output: 
		DIR_Haplotype + "/GVCFs/{cohort_wildcard}/{wildcard}.g.vcf.gz"
	priority: 1 # This allows haplotype caller to run before Create_Genomics_DB
	threads: 12
	shell:
		'/home/amunzur/gatk-4.2.0.0/gatk --java-options "-Xmx4g" HaplotypeCaller \
			--reference {input.PATH_hg38} \
			--input {input.SC_bam} \
			--output {output} \
			-ERC GVCF \
			--intervals {input.PATH_bed} \
			--read-filter AllowAllReadsReadFilter \
			--native-pair-hmm-threads {threads}'

# wildcard below gets replaced with the cohort name, like batch3, 4 etc 
rule Create_Genomics_DB:
	input:
		PATH_sample_map = DIR_sample_maps + "/{wildcard}.txt",
		PATH_bed = PATH_bed,
		DIR_temp = DIR_temp
	output: 
		PATH_DB = directory(DIR_Haplotype + "/Genomics_DB/{wildcard}")
	params:
		batch_size = 54
	threads: 16
	shell:
		'/home/amunzur/gatk-4.2.0.0/gatk --java-options "-Xmx20g -Xms20g" GenomicsDBImport \
			--genomicsdb-workspace-path {output} \
			--batch-size {params.batch_size} \
			--intervals {input.PATH_bed} \
			--sample-name-map {input.PATH_sample_map} \
			--tmp-dir {input.DIR_temp} \
			--read-filter AllowAllReadsReadFilter \
			--reader-threads {threads}'

# wildcard below gets replaced with the cohort name, like batch3, 4 etc 
rule Genotype_vcf:
	input:
		PATH_DB = DIR_Haplotype + "/Genomics_DB/{wildcard}",
		PATH_bed = PATH_bed,
		PATH_hg38 = PATH_hg38
	output: 
		DIR_Haplotype + "/GenotypeGVCFs/{wildcard}_genotyped.vcf.gz"
	params:
		batch_size = 54
	threads: 16
	shell:
		'/home/amunzur/gatk-4.2.0.0/gatk --java-options "-Xmx4g" GenotypeGVCFs \
			--reference {input.PATH_hg38} \
			--variant "gendb://{input.PATH_DB}" \
			--intervals {input.PATH_bed} \
			--read-filter AllowAllReadsReadFilter \
			--output {output}'

rule run_ANNOVAR_GATK: 
	input: 
		DIR_Haplotype + "/GenotypeGVCFs/{wildcard}_genotyped.vcf.gz"
	output: 
		ANNOVAR_GATK_output + "/{wildcard}/anno.hg38_multianno.txt"
	params: 
		actual_output_file = ANNOVAR_GATK_output + "/{wildcard}"		
	shell:
		'perl /groups/wyattgrp/software/annovar/annovar/table_annovar.pl {input} /groups/wyattgrp/software/annovar/annovar/humandb/ \
			-vcfinput \
			-buildver hg38 \
			-out {params.actual_output_file} \
			-remove \
			-protocol refGene,knownGene,avsnp147,exac03,cosmic70,clinvar_20170130,kaviar_20150923,gnomad_exome,dbnsfp33a,dbscsnv11,hrcr1,mcap,revel \
			-operation g,g,f,f,f,f,f,f,f,f,f,f,f \
			-nastring .'
