# do a merge based on column to combine metadata 
parse_anno_output <- function(DIR_variant_tables, mode, PATH_sample_list = NULL) {
	if (is.null(PATH_sample_list)) {
		if (mode == "somatic") {
			anno_df_list <- lapply(as.list(list.files(DIR_variant_tables, full.names = TRUE, pattern = ".tsv$")), return_anno_output_vardict_somatic)
		} else if (mode == "chip") {
			anno_df_list <- lapply(as.list(list.files(DIR_variant_tables, full.names = TRUE, pattern = ".tsv$")), return_anno_output_vardict_chip)
		} 
	} else {
		# We subset to a certain group of samples.
		sample_df <- as.data.frame(read_delim(PATH_sample_list, delim = "\t"))
		samples <- c(sample_df$cfDNA, sample_df$WBC)

		# List files in the directory matching the samples
		files_to_load <- list.files(DIR_variant_tables, full.names = TRUE, pattern = ".tsv$") # all files in the dir
		files_to_load <- files_to_load[sapply(files_to_load, function(file) any(sapply(samples, grepl, file)))] # choose a subset based on the provided samples file
		
		if (mode == "somatic") {
			anno_df_list <- lapply(as.list(files_to_load), return_anno_output_vardict_somatic)
		} else if (mode == "chip") {
			anno_df_list <- lapply(as.list(files_to_load), return_anno_output_vardict_chip)
		}
	}
	anno_df <- as.data.frame(do.call(rbind, anno_df_list))
	return(anno_df)
}

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
						PATH_before_filtering_chip,
						PATH_after_filtering){

	vars <- parse_anno_output(DIR_variant_tables_somatic, "somatic")
	vars <- add_patient_information(vars, PATH_sample_information)
	vars <- add_bg_error_rate(vars, bg)
	vars <- add_AAchange_effect(vars)
	write_csv(vars, PATH_before_filtering_chip)

	vars <- filter_variants_chip_or_germline("chip", vars, min_alt_reads, min_depth, min_VAF_low, max_VAF_low, min_VAF_high, max_VAF_high, min_VAF_bg_ratio, PATH_blacklist = PATH_blacklist)
	vars <- add_N_fraction(vars, DIR_mpileup, DIR_mpileup_filtered_chip, force = TRUE, file_pattern = "*cfDNA*")
	vars <- add_N_fraction(vars, DIR_mpileup, DIR_mpileup_filtered_chip, force = TRUE, file_pattern = "*gDNA*")
	vars$Variant_caller <- "Vardict"
	write_csv(vars, PATH_after_filtering)
	

 	return(vars)
}

MAINchip <- function(	min_VAF_low, 
						max_VAF_low,
						min_VAF_high,
						max_VAF_high,
						max_SBF_cfDNA,
						min_alt_reads,
						min_depth_gDNA,
						min_test_ref_ratio,
						min_VAF_bg_ratio, 
						DIR_variant_tables_chip,
						DIR_mpileup,
						DIR_mpileup_filtered_somatic,
						bg,
						PATH_sample_information,
						PATH_blacklist,
						PATH_before_filtering,
						PATH_after_filtering){

	vars <- parse_anno_output(DIR_variant_tables_chip, "chip")
	vars <- add_patient_information(vars, PATH_sample_information)
	vars <- add_bg_error_rate(vars, bg)
	vars <- add_AAchange_effect(vars, "Vardict")
	vars <- evaluate_strand_bias(vars)
	write_csv(vars, PATH_before_filtering)

	vars <- filter_variants_chip_or_germline("chip", vars, min_alt_reads, min_depth, min_VAF_low, max_VAF_low, min_VAF_high, max_VAF_high, min_VAF_bg_ratio, PATH_blacklist, blacklist = TRUE)
	vars <- add_N_fraction(vars, DIR_mpileup, DIR_mpileup_filtered_chip, force = TRUE)
	vars$Variant_caller <- "Vardict"
	write_csv(vars, PATH_after_filtering)
	
	vars <- combine_tumor_wbc(vars)
	vars <- filter(vars, Strand_bias_fishers_n != TRUE & Strand_bias_fishers_t != TRUE)
	write_csv(vars, PATH_final)

 	return(vars)
}