return_vardict_output <- function(DIR_vardict) {

	df_main <- as.data.frame(read_delim(DIR_vardict, delim = "\t")) 
	df <- df_main %>% 
		mutate(Sample_name = gsub(".tsv", "", basename(DIR_vardict))) %>%
		filter(!str_detect(Alt, 'N')) 	# Drop variants if the called variant has an "N" in it. 

	return(df)

}

# do a merge based on column to combine metadata 
combine_anno_vardict <- function(DIR_vardict, DIR_ANNOVAR) {

	# determine which dir to scan based on the cariant_type given

	anno_df_list <- lapply(as.list(list.files(DIR_ANNOVAR, full.names = TRUE, pattern = "\\.hg38_multianno.txt$")), return_anno_output)
	anno_df <- as.data.frame(do.call(rbind, anno_df_list)) %>%
				mutate(Sample_name = gsub(".hg38_multianno.txt", "", Sample_name), 
						Start = as.character(Start)) 

	vardict_df_list <- lapply(as.list(list.files(DIR_vardict, full.names = TRUE, pattern = "\\.tsv$")), return_vardict_output)
	vardict_df <- do.call(rbind, vardict_df_list) %>%
				mutate(Start = as.character(Start)) 

	combined <- left_join(vardict_df, anno_df, by = c("Sample_name", "Chr", "Start")) %>%
					mutate(Ref = vardict_df$Ref, 
						Alt = vardict_df$Alt,
						ExAC_ALL = replace_na(as.numeric(ExAC_ALL), 0), 
						gnomAD_exome_ALL = replace_na(as.numeric(gnomAD_exome_ALL), 0))

	combined$Ref.x <- NULL
	combined$Alt.x <- NULL
	combined$Ref.y <- NULL
	combined$Alt.y <- NULL

	return(combined)

}

add_bg_error_rate <- function(variants_df, bg) {

	# modify the vars df
	variants_df <- variants_df %>% mutate(error_type = paste0("mean_error", variants_df$Ref, "to", variants_df$Alt))
	
	# add the deletions 
	idx <- which(variants_df$variant == "deletion")
	variants_df$error_type[idx] <- "mean_errordel"
	
	# add the insertions
	idx <- which(variants_df$variant == "insertion")
	variants_df$error_type[idx] <- "mean_errorins"

	# modify the bg error rate df
	# bg <- read_delim(PATH_bg, delim = "\t")
	bg$chrom <- paste0("chr", bg$chrom)
	bg <- gather(bg, "error_type", "error_rate", starts_with("mean_error"))

	# now if we have a deletion, we want to recover the pos + 1 in the bg file, this is actual place where deletion begins 
	var_pos <- variants_df$Start
	deletion_idx <- which(variants_df$variant == "deletion")
	
	if (length(deletion_idx) > 0) {
		var_pos[deletion_idx] <- as.numeric(var_pos[deletion_idx]) + 1 # add 1 to the position of all the deletions
		variants_df$Start <- var_pos}

	# this merge adds an exta col with the error rate 
	bg$pos <- as.character(bg$pos)
	variants_df <- left_join(variants_df, bg, by = c("Chr" = "chrom", "Start" = "pos", "error_type" = "error_type")) # to add the error rates from bg to variants df

	# subtract 1 from the positions after merging
	if (length(deletion_idx) > 0) {
		deletion_idx <- which(variants_df$variant == "deletion")
		variants_df$Start[deletion_idx] <- as.numeric(variants_df$Start[deletion_idx]) - 1}

	return(variants_df)
}

# main function to run everything
MAIN <- function(cohort_name,
					THRESHOLD_ExAC_ALL, 
					VALUE_Func_refGene, 
					THRESHOLD_VarFreq, 
					THRESHOLD_Reads2, 
					THRESHOLD_VAF_bg_ratio, 
					DIR_vardict, 
					DIR_ANNOVAR,
					bg,
					PATH_bets_somatic,
					PATH_bets_germline,
					PATH_PCa_panel_2017,
					PATH_bed,
					DIR_depth_metrics,
					PATH_collective_depth_metrics,
					DIR_finland_bams, 
					DIR_temp,
					DIR_tnvstats,
					PATH_filter_tnvstats_script, 
					variant_caller){

	variant_caller = variant_caller # just so the function doesn't complain about us not using this argument
	combined <- combine_anno_vardict(DIR_vardict, DIR_ANNOVAR) # add annovar annots to the varscan outputs

	# add patient id
	x <- str_split(combined$Sample_name, "-") # split the sample name
	combined$patient_id <- paste(lapply(x, "[", 1), lapply(x, "[", 2), lapply(x, "[", 3), sep = "-") # paste 2nd and 3rd elements to create the sample name 

	combined <- add_bg_error_rate(combined, bg) # background error rate
	combined <- add_AAchange_effect(combined) # protein annot and the effects

	combined <- combined %>%
						mutate(Total_reads = Reads1 + Reads2, 
								VAF_bg_ratio = VarFreq/error_rate) %>%
						select(Sample_name, patient_id, Chr, Start, Ref, Alt, VarFreq, Reads1, Reads2, StrandBias_Fisher_pVal, StrandBias_OddsRatio, REF_Fw, REF_Rv, ALT_Fw, ALT_Rv, Func.refGene, Gene.refGene, AAChange.refGene, Protein_annotation, Effects, ExAC_ALL, variant, error_rate, VAF_bg_ratio, Total_reads) %>%
						filter(ExAC_ALL <= THRESHOLD_ExAC_ALL, 
								Func.refGene != VALUE_Func_refGene,
								VAF_bg_ratio >= THRESHOLD_VAF_bg_ratio, # vaf should be at least 15 times more than the bg error rate
								StrandBias_Fisher_pVal > 0.05) 

	# a common naming convention i will be sticking to from now on
	names(combined) <- c("Sample_name", "Patient_ID", "Chrom", "Position", "Ref", "Alt", "VAF", "Ref_reads", "Alt_reads", "StrandBias_Fisher_pVal", "StrandBias_OddsRatio", "REF_Fw", "REF_Rv", "ALT_Fw", "ALT_Rv", "Function", "Gene", "AAchange", "Protein_annotation", "Effects", "ExAC_ALL", "Variant", "Error_rate", "VAF_bg_ratio", "Total_reads")

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
	combined <- compare_with_bets(PATH_bets_somatic, PATH_bets_germline, PATH_PCa_panel_2017, combined) # adds three new columns

	# add an extra col for alerting the user if the variant isn't found, despite gene being in the bets
	combined <- combined %>% mutate(Status = case_when(
										(In_germline_bets == FALSE & In_2017_PCa == TRUE) ~ "ALERT", 
										(In_germline_bets == TRUE & In_2017_PCa == TRUE) ~ "Great",
										(In_germline_bets == TRUE & In_2017_PCa == FALSE) ~ "Error",
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

		combined <- add_finland_readcounts(combined, filtered_tnvstats, variant_caller)

	}

	return(combined)

} # end of function

combine_and_save <- function(variants, PATH_validated_variants, PATH_SAVE_chip_variants){

	validated_vars <- compare_with_jacks_figure(PATH_validated_variants, variants)
	# variants$detected <- variants$Position %in% validated_vars$Position # add a new col to show if the same var was detected in jacks figure

	# PATH_SAVE_chip_variants <- "/groups/wyattgrp/users/amunzur/chip_project/variant_lists/chip_variants.csv"
	write_csv(variants, PATH_SAVE_chip_variants) # snv + indel, csv
	write_delim(variants, gsub(".csv", ".tsv", PATH_SAVE_chip_variants), delim = "\t") # snv + indel, tsv

	write_csv(validated_vars, gsub("chip_variants.csv", "validated_variants.csv", PATH_SAVE_chip_variants)) # jack df as csv
	write_delim(validated_vars, gsub("chip_variants.csv", "validated_variants.tsv", PATH_SAVE_chip_variants), delim = "\t") # jack df as tsv
}