return_anno_output_mutect_chip <- function(PATH_variant_table) {

	print(PATH_variant_table)
	df <- as.data.frame(read.delim(PATH_variant_table, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE))
	colnames(df) <- gsub(paste0(file_path_sans_ext(basename(PATH_variant_table)), "."), "", colnames(df)) 

	df <- df %>%
		separate(col = SB, sep = ",", into = c("Ref_forward", "Ref_reverse", "Alt_forward", "Alt_reverse"), remove = TRUE) %>%
		select(-EVENTLENGTH) %>% # really didnt need to include this column in the variant tables
		mutate(Sample = gsub(".tsv", "", basename(PATH_variant_table)), 
			   Sample_type = str_extract(basename(PATH_variant_table), "cfDNA|WBC"),
			   Sample_name = file_path_sans_ext(basename(PATH_variant_table)),
			   Date_collected = str_split(Sample_name, "[-_]")  %>% sapply(tail, 1), 
			   nchar_ref = nchar(REF), 
			   nchar_alt = nchar(ALT), 
			   ALT = ifelse(grepl(",", ALT), sub(",.*", "", ALT), ALT)) %>% # Select the first element if ALT contains commas
		mutate(across(c(Alt_forward, Alt_reverse, Ref_forward, Ref_reverse), as.numeric)) %>%
		mutate(TYPE = case_when(nchar_ref > nchar_alt ~ "Deletion", 
			   					   nchar_ref < nchar_alt ~ "Insertion", 
								   nchar_ref == nchar_alt ~ "SNV"), 
			   VAF = (Alt_forward+Alt_reverse)/(Ref_forward+Ref_reverse+Alt_forward+Alt_reverse), 
			   Depth = Ref_forward+Ref_reverse+Alt_forward+Alt_reverse) %>%
		select(-nchar_ref, -nchar_alt) %>%
		select(Sample_name, Sample_type, Date_collected, CHROM, POS, REF, ALT, TYPE, VAF, Depth, Alt_forward, Ref_forward, Alt_reverse, Ref_reverse, Func.refGene, Gene.refGene, ExonicFunc.refGene, 
				 AAChange.refGene, cosmic97_coding, avsnp150, CLNALLELEID, CLNSIG)

	names(df) <- c("Sample_name", "Sample_type", "Date_collected", "Chrom", "Position", "Ref", "Alt", "Type", "VAF", "Depth", "Alt_forward", "Ref_forward", "Alt_reverse", "Ref_reverse", "Function", "Gene", "Consequence", 
				   "AAchange", "cosmic97_coding", "avsnp150", "CLNALLELEID", "CLNSIG")
	return(df)

}

return_anno_output_mutect_somatic <- function(PATH_variant_table) {

	print(PATH_variant_table)
	df <- as.data.frame(read.delim(PATH_variant_table, stringsAsFactors = FALSE, check.names = FALSE))
	df$Sample_name_n = str_split(colnames(df)[grep("WBC", colnames(df))][1], "\\.")[[1]][1]
	colnames(df) <- gsub("^.*[_-]WBC[-_].*\\.", "WBC_", colnames(df)) 
	colnames(df) <- gsub("^.*[_-]cfDNA[-_].*\\.", "cfDNA_", colnames(df)) 

	df <- df %>%
		separate(col = cfDNA_SB, sep = ",", into = c("Ref_forward_t", "Ref_reverse_t", "Alt_forward_t", "Alt_reverse_t"), remove = TRUE) %>%
		separate(col = WBC_SB, sep = ",", into = c("Ref_forward_n", "Ref_reverse_n", "Alt_forward_n", "Alt_reverse_n"), remove = TRUE) %>%
		mutate(Patient_id = gsub("GU-", "", str_split(file_path_sans_ext(basename(PATH_variant_table)), "_")[[1]][1]),
			   Sample_name_t = file_path_sans_ext(basename(PATH_variant_table)),
			   Date_collected = tail(str_split(file_path_sans_ext(basename(PATH_variant_table)), "[-_]")[[1]], n = 1), 
			   TYPE = case_when(
							nchar(REF) > nchar(ALT) ~ "Deletion", 
							nchar(REF) < nchar(ALT) ~ "Insertion",
							nchar(REF) == nchar(ALT) ~ "SNV", 
							TRUE ~ "Error"), 
			   across(c(Alt_forward_t, Alt_reverse_t, Ref_forward_t, Ref_reverse_t, Alt_forward_n, Alt_reverse_n, Ref_forward_n, Ref_reverse_n), as.numeric), 
			   VAF_t = (Alt_forward_t + Alt_reverse_t)/(Ref_forward_t + Ref_reverse_t + Alt_forward_t + Alt_reverse_t), 
			   VAF_n = (Alt_forward_n + Alt_reverse_n)/(Ref_forward_n + Ref_reverse_n + Alt_forward_n + Alt_reverse_n), 
			   Depth_t = Ref_forward_t + Ref_reverse_t + Alt_forward_t + Alt_reverse_t, 
			   Depth_n = Ref_forward_n + Ref_reverse_n + Alt_forward_n + Alt_reverse_n, 
			   tumor_to_normal_VAF_ratio = VAF_t/VAF_n) %>%
		select(Patient_id, Sample_name_t, Date_collected, CHROM, POS, REF, ALT, TYPE, VAF_t, Depth_t, Alt_forward_t, Ref_forward_t, Alt_reverse_t, Ref_reverse_t, Func.refGene, Gene.refGene, ExonicFunc.refGene, 
		AAChange.refGene, cosmic97_coding, avsnp150, CLNALLELEID, CLNSIG, Sample_name_n, VAF_n, Depth_n, Alt_forward_n, Ref_forward_n, Alt_reverse_n, Ref_reverse_n, tumor_to_normal_VAF_ratio)

	names(df) <- c("Patient_id", "Sample_name_t", "Date_collected", "Chrom", "Position", "Ref", "Alt", "Type", "VAF_t", "Depth_t", "Alt_forward_t", "Ref_forward_t", "Alt_reverse_t", "Ref_reverse_t", "Function", "Gene", "Consequence", 
	"AAchange", "cosmic97_coding", "avsnp150", "CLNALLELEID", "CLNSIG", "Sample_name_n", "VAF_n", "Depth_n", "Alt_forward_n", "Ref_forward_n", "Alt_reverse_n", "Ref_reverse_n", "tumor_to_normal_VAF_ratio")

	return(df)
}


# do a merge based on column to combine metadata 
parse_anno_output <- function(DIR_variant_tables, mode, PATH_sample_list = NULL) {
	# No subsetting, all files in the directory are used.
	if (is.null(PATH_sample_list)) {
		if (mode == "somatic") {
			anno_df_list <- lapply(as.list(list.files(DIR_variant_tables, full.names = TRUE, pattern = ".tsv$")), return_anno_output_mutect_somatic)
		} else if (mode == "chip") {
			anno_df_list <- lapply(as.list(list.files(DIR_variant_tables, full.names = TRUE, pattern = ".tsv$")), return_anno_output_mutect_chip)
		}
	} else {
		# We subset to a certain group of samples.
		sample_df <- as.data.frame(read_delim(PATH_sample_list, delim = "\t"))
		samples <- c(sample_df$cfDNA, sample_df$WBC)

		# List files in the directory matching the samples
		files_to_load <- list.files(DIR_variant_tables, full.names = TRUE, pattern = ".tsv$") # all files in the dir
		files_to_load <- files_to_load[sapply(files_to_load, function(file) any(sapply(samples, grepl, file)))] # choose a subset based on the provided samples file
		
		if (mode == "somatic") {
			anno_df_list <- lapply(as.list(files_to_load), return_anno_output_mutect_somatic)
		} else if (mode == "chip") {
			anno_df_list <- lapply(as.list(files_to_load), return_anno_output_mutect_chip)
		}
	}
	anno_df <- as.data.frame(do.call(rbind, anno_df_list))
	return(anno_df)
}

# main function to run everything
MAIN <- function(
					VALUE_Func_refGene, 
					THRESHOLD_VarFreq, 
					THRESHOLD_alt_reads, 
					THRESHOLD_VAF_bg_ratio, 
					DIR_variant_tables_chip,
					DIR_mpileup,
					DIR_mpileup_filtered,
					bg,
					PATH_cancer_type,
					PATH_panel_genes,
					PATH_bed,
					DIR_depth_metrics,
					PATH_collective_depth_metrics,
					DIR_temp,
					variant_caller,
					PATH_blacklist, 
					PATH_before_filtering, 
					PATH_after_filtering, 
					PATH_final){

	variant_caller = variant_caller # just so the function doesn't complain about us not using this argument
	vars <- parse_anno_output(DIR_variant_tables_chip)
	vars <- add_patient_information(vars, PATH_sample_information)
	vars <- add_bg_error_rate(vars, bg)
	vars <- add_AAchange_effect(vars, "Mutect2")
	vars <- evaluate_strand_bias(vars)
	write_csv(vars, PATH_before_filtering)

	vars <- filter_variants_chip(vars, min_alt_reads, min_depth, min_VAF_low, max_VAF_low, min_VAF_high, max_VAF_high, min_VAF_bg_ratio, PATH_blacklist, blacklist = TRUE)
	vars <- add_N_fraction(vars, DIR_mpileup, DIR_mpileup_filtered_chip, force = FALSE)
	vars$Variant_caller <- "Mutect2"
	write_csv(vars, PATH_after_filtering)
	
	vars <- combine_tumor_wbc(vars)
	vars <- filter(vars, Strand_bias_fishers_n != TRUE & Strand_bias_fishers_t != TRUE)
	write_csv(vars, PATH_final)

	return(combined)

}
