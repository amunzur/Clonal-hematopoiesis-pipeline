# reformat such that each variant is associated with a sample, vaf etc
parse_variant_info <- function(PATH_vcf_table){
	df <- read_delim(PATH_vcf_table, delim = "\t")
	df_AD <- df[, grep("AD", names(df))] # slice of df with allele count info
	variant_info <- df[, c(1, 2, 3, 4, 5)]	# info about variants: chrom, pos etc
	df <- as.data.frame(cbind(variant_info, df_AD)) # put back together 
	# names(df) <- gsub("\\.GT|\\.AD|\\.DP|\\.GQ|\\.PL", "", names(df)) # remove the ".AD" at the end of the column names

	df <- df %>%
			gather("sample_names", "read_depth", ends_with(".AD")) %>%
			mutate(
				read_depth = sub(",", "_", read_depth), 
				ALT = gsub("\\*", 0, ALT), 
				sample_names = gsub("\\.AD", "", sample_names)) %>%
			separate(
				col = read_depth, 
				into = c("Ref_reads", "Alt_reads"), sep = "_") %>%
			separate_rows(
				ALT, 
				Alt_reads) %>%
			mutate(
				Alt_reads = as.numeric(Alt_reads), 
				Ref_reads = as.numeric(Ref_reads), 
				VAF = Alt_reads/(Alt_reads + Ref_reads))

	return(df)
}

# This function is needed to change the format of REF and ALT fields as used in the ANNOVAR. This helps doing a merge later on
# with the annovar results and the sample information.
annovariaze_sample_info <- function(sample_info) {

	df <- sample_info %>%
		mutate(
			nchar_REF = nchar(REF), 
			nchar_ALT = nchar(ALT), 
			nchar_dif = nchar_REF - nchar_ALT, 
			TYPE = case_when(
				nchar_dif == 0 ~ "snv", 
				nchar_dif < 0 ~ "insertion",
				nchar_dif > 0 ~ "deletion"),
			POS_new = case_when(
				nchar_REF == 1 ~ POS, # SNP 
				nchar_REF > 1 ~ POS + 1, # more than 1 base at REF, add 1 to pos
				TRUE ~ 9999), 
			END_new = case_when(
				nchar_REF == 1 ~ POS, # SNP
				nchar_REF > 1 ~ POS + nchar_REF - 1, # more than 
				TRUE ~ 9999), 
			REF = case_when(
				TYPE == "snv" ~ REF, 
				TYPE != "snv" ~ substring(REF, 2)),
			ALT = case_when(
				TYPE == "snv" ~ ALT, 
				TYPE != "snv" ~ substring(ALT, 2)), 
			REF = case_when(
				nchar(REF) == 0 ~ "-",
				nchar(REF) != 0 ~ REF),
			ALT = case_when(
				nchar(ALT) == 0 ~ "-",
				nchar(ALT) != 0 ~ ALT)) %>%
		select(CHROM, POS, REF, ALT, TYPE, sample_names, Ref_reads, Alt_reads, VAF) %>%
		rename(
			Sample_name = sample_names, 
			variant = TYPE)

	return(df)

}

# Going back to normal reference and alt alleles after annovarizing them
normalize_ref_alt <- function(variant_df, sample_info){

	# sample_info has the ref and alt alleles in the format we want, we will do a merge on that
	# change the "annovarized ref and alt"s back to the format that is in accordance with vardict and varscan. This will help with merging later on.
	names(sample_info) <- c("Chrom", "Position", "Ref", "Alt", "Variant", "Sample_name", "Ref_reads", "Alt_reads", "VAF")
	
	variant_df <- variant_df[, -which(names(variant_df) %in% c("Ref", "Alt"))]
	variant_df <- variant_df %>% 
						left_join(sample_info) %>%
						select(
							"Sample_name", 
							"Sample_type", 
							"Patient_ID", 
							"Cohort_name", 
							"Chrom", 
							"Position", 
							"Ref", 
							"Alt", 
							"VAF", 
							"Ref_reads", 
							"Alt_reads", 
							"Function", 
							"Gene", 
							"AAchange",
							"Protein_annotation", 
							"Effects", 
							"ExAC_ALL", 
							"Variant", 
							"Error_rate", 
							"VAF_bg_ratio", 
							"Total_reads", 
							"Duplicate", 
							"Depth")

	return(variant_df)
}

main <- function(
	PATH_vcf_table, 
	PATH_ANNOVAR_output,
	cohort_name, 
	THRESHOLD_ExAC_ALL, 
	VALUE_Func_refGene, 
	THRESHOLD_VarFreq, 
	THRESHOLD_VAF_bg_ratio,
	bg, 
	PATH_bets_somatic, 
	PATH_bets_germline, 
	PATH_panel_genes, 
	PATH_bed, 
	DIR_depth_metrics, 
	PATH_collective_depth_metrics, 
	DIR_finland_bams, 
	DIR_temp, 
	DIR_tnvstats, 
	PATH_filter_tnvstats_script,
	variant_caller){
		
	variant_caller = variant_caller

	sample_info <- parse_variant_info(PATH_vcf_table) # each row is a sample, each variant is evaluated in each sample
	sample_info <- annovariaze_sample_info(sample_info)

	annovar <- read_delim(PATH_ANNOVAR_output, delim = "\t")
	combined <- left_join(
					annovar, 
					sample_info, 
					by = c(
						"Chr" = "CHROM", 
						"Start" = "POS", 
						"Ref" = "REF", 
						"Alt" = "ALT")) %>%
				mutate(Start = as.character(Start))

	combined <- add_patient_id(combined, cohort_name)
	combined <- add_bg_error_rate(combined, bg) # background error rate
	combined <- add_AAchange_effect(combined) # protein annot and the effects
	combined$Cohort_name <- cohort_name
	combined <- add_sample_type(combined)

	combined <- combined %>%
						mutate(Total_reads = Ref_reads + Alt_reads, 
								VAF_bg_ratio = VAF/error_rate, 
								ExAC_ALL = replace_na(as.numeric(ExAC_ALL), 0), 
								gnomAD_exome_ALL = replace_na(as.numeric(gnomAD_exome_ALL), 0)) %>%
						select(Sample_name, Sample_type, patient_id, Cohort_name, Chr, Start, Ref, Alt, VAF, Ref_reads, Alt_reads, Func.refGene, Gene.refGene, AAChange.refGene, Protein_annotation, Effects, ExAC_ALL, variant, error_rate, VAF_bg_ratio, Total_reads) %>%
						filter(ExAC_ALL <= THRESHOLD_ExAC_ALL, 
								Func.refGene != VALUE_Func_refGene,
								VAF_bg_ratio >= THRESHOLD_VAF_bg_ratio) 

	# a common naming convention i will be sticking to from now on
	names(combined) <- c("Sample_name", "Sample_type", "Patient_ID", "Cohort_name", "Chrom", "Position", "Ref", "Alt", "VAF", "Ref_reads", "Alt_reads", "Function", "Gene", "AAchange", "Protein_annotation", "Effects", "ExAC_ALL", "Variant", "Error_rate", "VAF_bg_ratio", "Total_reads")

	idx <- grep("splicing", combined$Function)
	combined$Effects[idx] <- "splicing"

	# add a new col indicating if the variant is duplicated or not
	combined <- identify_duplicates(combined)

	# now filtering based on vaf, read support and depth etc. 
	combined <- combined %>%
				filter((Total_reads >= 1000 & VAF >= 0.005) | 
						(Total_reads <= 1000 & Alt_reads >= 5), 
						VAF < THRESHOLD_VarFreq)

	combined <- subset_to_panel(PATH_bed, combined) # subset to panel
	combined <- add_depth(DIR_depth_metrics, PATH_collective_depth_metrics, combined) # add depth information at these positions 
	combined <- compare_with_bets(PATH_bets_somatic, PATH_bets_germline, PATH_panel_genes, combined) # adds three new columns

	# add an extra col for alerting the user if the variant isn't found, despite gene being in the bets
	combined <- combined %>% mutate(Status = case_when(
										(In_germline_bets == FALSE & In_panel == TRUE) ~ "ALERT", 
										(In_germline_bets == TRUE & In_panel == TRUE) ~ "Great",
										(In_germline_bets == TRUE & In_panel == FALSE) ~ "Error",
										TRUE ~ "OK"), 
								Position = as.numeric(Position))

	if (cohort_name == "new_chip_panel"){

		# add the sample ID from the finland bams
		finland_sample_IDs <- gsub(".bam", "", grep("*.bam$", list.files(DIR_finland_bams), value = TRUE))
		combined$Sample_name_finland <- unlist(lapply(as.list(gsub("GUBB", "GU", combined$Patient_ID)), function(x) grep(x, finland_sample_IDs, value = TRUE)))
		combined <- combined %>% relocate(Sample_name_finland, .after = Sample_name)

		# write_csv(combined, "/groups/wyattgrp/users/amunzur/pipeline/results/temp/vardict_combined.csv")
		# combined <- read_csv("/groups/wyattgrp/users/amunzur/pipeline/results/temp/vardict_combined.csv")

		message("Filtering tnvstats right now.")
		filtered_tnvstats <- filter_tnvstats_by_variants(
								variants_df = combined, 
								DIR_tnvstats = DIR_tnvstats, 
								PATH_temp = file.path(DIR_temp, variant_caller), 
								PATH_filter_tnvstats_script = PATH_filter_tnvstats_script, 
								identifier = variant_caller)

		combined <- add_finland_readcounts(combined, filtered_tnvstats, variant_caller) }

	# Check for duplicated rows just before returning the object.
	combined <- check_duplicated_rows(combined, TRUE)
	combined <- normalize_ref_alt(combined, sample_info)

	return(combined)}

combine_and_save <- function(variants, PATH_SAVE_chip_variants){

	dir.create(dirname(PATH_SAVE_chip_variants))

	write_csv(variants, PATH_SAVE_chip_variants) # snv + indel, csv
	write_delim(variants, gsub(".csv", ".tsv", PATH_SAVE_chip_variants), delim = "\t") # snv + indel, tsv
}