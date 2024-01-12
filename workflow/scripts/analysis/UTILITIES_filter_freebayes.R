

filter_somatic <- function(vars, min_alt_reads_cfDNA, min_depth_gDNA, max_VAF, min_VAF, min_test_ref_ratio, min_VAF_bg_ratio, max_SBF_cfDNA) {
	vars <- vars %>% 
		mutate(Status = "Somatic", 
			   test_ref_ratio = cfDNA_AF/gDNA_AF) %>%
		filter(Func.refGene == "exonic",
			   Effects != "synonymous", 
			   cfDNA_VD >= min_alt_reads_cfDNA, # alt reads in cfDNA
			   gDNA_DP >= min_depth_gDNA, # total depth at gDNA
			   cfDNA_AF < max_VAF, # max vaf
			   cfDNA_AF > min_VAF, # min vaf
			   test_ref_ratio > min_test_ref_ratio,
			   VAF_bg_ratio_cfDNA >= min_VAF_bg_ratio, # bg error rate
			   cfDNA_SBF < max_SBF_cfDNA) # strand bias 
		
	return(vars)
}


MAINsomatic <- function(max_VAF, 
						min_VAF,
						max_SBF_cfDNA,
						min_alt_reads_cfDNA,
						min_depth_gDNA,
						min_test_ref_ratio,
						min_VAF_bg_ratio, 
						DIR_variant_tables_somatic,
						DIR_mpileup,
						DIR_mpileup_filtered_somatic,
						bg,
						PATH_before_filtering,
						PATH_after_filtering){

	vars <- parse_anno_output(DIR_variant_tables_chip)
	vars <- add_patient_information(vars, PATH_cancer_type)
	vars <- add_bg_error_rate(vars, bg)
	vars <- add_AAchange_effect(vars, "freebayes")
	vars <- evaluate_strand_bias(vars)
	write_csv(vars, PATH_before_filtering)

	vars <- filter_somatic(vars, min_alt_reads_cfDNA, min_depth_gDNA, max_VAF, min_VAF, min_test_ref_ratio, min_VAF_bg_ratio, max_SBF_cfDNA)
	vars <- add_N_fraction(vars, DIR_mpileup, DIR_mpileup_filtered_chip, force = FALSE, file_pattern = "*cfDNA*")
	vars <- add_N_fraction(vars, DIR_mpileup, DIR_mpileup_filtered_chip, force = TRUE, file_pattern = "*gDNA*")
	write_csv(vars, PATH_after_filtering)
	vars$Variant_caller <- "Freebayes"

 	return(vars)
}


MAINchip <- function(	max_VAF, 
						min_VAF,
						max_SBF_cfDNA,
						min_alt_reads_cfDNA,
						min_depth_gDNA,
						min_test_ref_ratio,
						min_VAF_bg_ratio, 
						DIR_variant_tables_chip,
						DIR_mpileup,
						DIR_mpileup_filtered_somatic,
						bg,
						PATH_before_filtering,
						PATH_after_filtering){

	vars <- parse_anno_output(DIR_variant_tables_chip, PATH_sample_list = PATH_sample_list)
	vars <- add_patient_information(vars, PATH_sample_information)
	vars <- add_bg_error_rate(vars, bg)
	vars <- add_AAchange_effect(vars, "freebayes")
	vars <- evaluate_strand_bias(vars)
	write_csv(vars, PATH_before_filtering)

	vars <- filter_variants_chip(vars, min_alt_reads, min_depth, min_VAF_low, max_VAF_low, min_VAF_high, max_VAF_high, min_VAF_bg_ratio, PATH_blacklist, blacklist = TRUE)
	vars <- add_N_fraction(vars, DIR_mpileup, DIR_mpileup_filtered_chip, force = FALSE)
	vars$Variant_caller <- "freebayes"
	write_csv(vars, PATH_after_filtering)
	
	vars <- combine_tumor_wbc(vars)
	vars <- filter(vars, Strand_bias_fishers_n != TRUE & Strand_bias_fishers_t != TRUE)
	write_csv(vars, PATH_final)

 	return(vars)
}