return_anno_output_vardict <- function(PATH_variant_table) {

	df <- as.data.frame(read_delim(PATH_variant_table, delim = "\t"))
	colnames(df) <- gsub("^.*_gDNA_.*\\.", "gDNA_", colnames(df)) 
	colnames(df) <- gsub("^.*_cfDNA_.*\\.", "cfDNA_", colnames(df)) 

	# Add patient id
	df$Patient_id <- sub(".*\\/([^\\/]+)_cfDNA.*", "\\1", PATH_variant_table)
	df <- df[, c("Patient_id", colnames(df)[-which(names(df) == "Patient_id")])] # move to the beginning

	# Indicate insertion, deletion or SNV
	df <- mutate(df, TYPE = case_when(
					nchar(REF) > nchar(ALT) ~ "Deletion", 
					nchar(REF) < nchar(ALT) ~ "Insertion",
					nchar(REF) == nchar(ALT) ~ "SNV", 
					TRUE ~ "Error")) %>%
		  mutate(test_ref_ratio = cfDNA_AF/gDNA_AF) %>%
		  separate(cfDNA_ALD, into = c("cfDNA_ALT_Fw", "cfDNA_ALT_Rv"), sep = ",", remove = TRUE) %>%
		  separate(cfDNA_RD, into = c("cfDNA_REF_Fw", "cfDNA_REF_Rv"), sep = ",", remove = TRUE) %>%
		  separate(gDNA_ALD, into = c("gDNA_ALT_Fw", "gDNA_ALT_Rv"), sep = ",", remove = TRUE) %>%
		  separate(gDNA_RD, into = c("gDNA_REF_Fw", "gDNA_REF_Rv"), sep = ",", remove = TRUE)


	return(df)
}

# do a merge based on column to combine metadata 
parse_anno_output <- function(DIR_variant_tables) {

	anno_df_list <- lapply(as.list(list.files(DIR_variant_tables, full.names = TRUE, pattern = ".tsv$")), return_anno_output_vardict)
	anno_df <- as.data.frame(do.call(rbind, anno_df_list))
	return(anno_df)
}

add_bg_error_rate <- function(vars, bg) {

	# modify the vars df
	vars <- vars %>% mutate(ERROR_TYPE = paste0("mean_error", vars$REF, "to", vars$ALT))
	
	# add the deletions 
	vars$ERROR_TYPE[grepl("Deletion", vars$TYPE)] <- "mean_errordel"
	vars$ERROR_TYPE[grepl("Insertion", vars$TYPE)] <- "mean_errorins"
	
	# modify the bg error rate df
	# bg <- read_delim(PATH_bg, delim = "\t")
	bg <- gather(bg, "error_type", "error_rate", starts_with("mean_error"))
	names(bg) <- c("CHROM", "POS", "REF", "ERROR_TYPE", "ERROR_RATE")

	vars <- left_join(vars, bg) # to add the error rates from bg to variants df
	vars <- mutate(vars, 
				   VAF_bg_ratio_cfDNA = cfDNA_AF/ERROR_RATE, 
				   VAF_bg_ratio_gDNA = gDNA_AF/ERROR_RATE)
	return(vars)
}

evaluate_strand_bias <- function(combined) {
	combined <- combined %>%
		mutate(ALT_Fw_total = ALT_Fw_t + ALT_Fw_n, 
			   ALT_Rv_total = ALT_Rv_t + ALT_Rv_n) %>%
		filter(ALT_Fw_total > 0 & ALT_Rv_total > 0)

	return(combined)
	
	}

filter_somatic <- function(vars, min_alt_reads_cfDNA, min_ref_reads_gDNA, max_VAF, min_VAF, min_test_ref_ratio, min_VAF_bg_ratio, max_SBF_cfDNA) {
	vars <- vars %>% 
		mutate(Status = "Somatic", 
			   test_ref_ratio = cfDNA_AF/gDNA_AF) %>%
		filter(Func.refGene == "exonic",
			   Effects != "synonymous", 
			   cfDNA_VD >= min_alt_reads_cfDNA, # alt reads in cfDNA
			   gDNA_DP >= min_ref_reads_gDNA, # total depth at gDNA
			   cfDNA_AF < max_VAF, # max vaf
			   cfDNA_AF > min_VAF, # min vaf
			   test_ref_ratio > min_test_ref_ratio,
			   VAF_bg_ratio_cfDNA >= min_VAF_bg_ratio, # bg error rate
			   cfDNA_SBF < max_SBF_cfDNA) # strand bias 
		
	return(vars)
}


# main function to run everything
MAIN <- function(
					max_VAF, 
					min_VAF,
					max_SBF_cfDNA,
					min_alt_reads_cfDNA,
					min_ref_reads_gDNA,
					min_test_ref_ratio,
					min_VAF_bg_ratio, 
					DIR_variant_tables,
					DIR_mpileup,
					DIR_mpileup_filtered_somatic,
					bg,
					PATH_before_filtering,
					PATH_after_filtering){

	vars <- parse_anno_output(DIR_variant_tables)
	vars <- add_bg_error_rate(vars, bg)
	vars <- add_AAchange_effect(vars)
	write_csv(vars, PATH_before_filtering)

	vars <- filter_somatic(vars, min_alt_reads_cfDNA, min_ref_reads_gDNA, max_VAF, min_VAF, min_test_ref_ratio, min_VAF_bg_ratio, max_SBF_cfDNA)
	vars <- add_N_fraction(vars, DIR_mpileup, DIR_mpileup_filtered_somatic, force = FALSE, file_pattern = "*cfDNA*")
	vars <- add_N_fraction(vars, DIR_mpileup, DIR_mpileup_filtered_somatic, force = TRUE, file_pattern = "*gDNA*")
	write_csv(vars, PATH_after_filtering)

 	return(vars)

}