return_anno_output <- function(PATH_ANNOVAR_output, PATH_ANNOVAR_input) {

	PATH_ANNOVAR_input <- gsub("outputs", "inputs", PATH_ANNOVAR_output)
    PATH_ANNOVAR_input <- gsub(".hg38_multianno.txt", "_anno.tsv", PATH_ANNOVAR_input)

    df_main <- read_delim(PATH_ANNOVAR_output, delim = "\t") %>%
			   select(Chr, Start, Ref, Alt, Func.refGene, Gene.refGene, Func.knownGene, Gene.knownGene, AAChange.refGene, ExAC_ALL, gnomAD_exome_ALL) %>%
               mutate_at(c("Start", "ExAC_ALL", "gnomAD_exome_ALL"), as.numeric) %>%
			   replace_na(list(ExAC_ALL = 0, gnomAD_exome_ALL = 0)) %>%
               slice(-1)

    df_supp <- read_delim(PATH_ANNOVAR_input, "\t") %>%
               select(Chrom, Position, Ref, Alt, sample_names_n, alt_reads_n, total_reads_n, mapq_n, VAF_n, sample_names_t, alt_reads_t, total_reads_t, mapq_t, VAF_t, Type) %>%
               rename(Chr = Chrom, Start = Position, variant = Type)

	df_main <- left_join(df_main, df_supp, by = c("Chr", "Start", "Ref", "Alt"))
    df_main$variant <- tolower(df_main$variant)

	return(df_main)

}

# do a merge based on column to combine metadata 
parse_anno_output <- function(DIR_ANNOVAR) {

	anno_df_list <- lapply(as.list(list.files(DIR_ANNOVAR, full.names = TRUE, pattern = "GUBB.+hg38_multianno.txt$")), return_anno_output)
	anno_df <- as.data.frame(do.call(rbind, anno_df_list)) %>%
				mutate(Start = as.character(Start)) 

	return(anno_df)

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

# based on the cohort_name, add the patient id
add_patient_id <- function(variant_df){

	x <- str_split(variant_df$sample_names_t, "_gDNA|_WBC|_cfDNA")
	variant_df$patient_id <- unlist(lapply(x, "[", 1))

	return(variant_df)

}

# main function to run everything
MAIN <- function(
					VALUE_Func_refGene, 
					THRESHOLD_VarFreq, 
					THRESHOLD_Reads2, 
					THRESHOLD_VAF_bg_ratio, 
					DIR_ANNOVAR,
					bg,
					PATH_panel_genes,
					PATH_bed,
					DIR_depth_metrics,
					PATH_collective_depth_metrics,
					DIR_temp,
					DIR_tnvstats,
					PATH_filter_tnvstats_script, 
					variant_caller){

	variant_caller = variant_caller # just so the function doesn't complain about us not using this argument
	combined <- parse_anno_output(DIR_ANNOVAR)
	combined <- add_patient_id(combined)
	combined <- add_bg_error_rate(combined, bg)
	combined <- add_AAchange_effect(combined)
  
	write_csv(combined, "/groups/wyattgrp/users/amunzur/COMPOST_BIN/mutect_before_filtering.csv")

	combined <- combined %>%
						mutate(VAF_bg_ratio = VAF_n/error_rate, Ref_reads_n = total_reads_n - alt_reads_n, Ref_reads_t = total_reads_t - alt_reads_t) %>%
						select(sample_names_n, patient_id, Chr, Start, Ref, Alt, VAF_n, Ref_reads_n, alt_reads_n, Func.refGene, Gene.refGene, AAChange.refGene, Protein_annotation, Effects, ExAC_ALL, variant, error_rate, VAF_bg_ratio, total_reads_n, VAF_t, Ref_reads_t, alt_reads_t) %>%
						filter(Func.refGene != VALUE_Func_refGene,
								VAF_bg_ratio >= THRESHOLD_VAF_bg_ratio) 
	names(combined) <- c("Sample_name", "Patient_ID", "Chrom", "Position", "Ref", "Alt", "VAF_n", "Ref_reads_n", "Alt_reads_n", "Function", "Gene", "AAchange", "Protein_annotation", "Effects", "ExAC_ALL", "Variant", "Error_rate", "VAF_bg_ratio", "Total_reads", "VAF_t", "Ref_reads_t", "Alt_reads_t")

	combined <- combined %>%
				filter((Total_reads >= 1000 & VAF_n >= 0.005) | 
						(Total_reads <= 1000 & Alt_reads_n >= 5), 
						VAF_n < THRESHOLD_VarFreq) %>%
				mutate(Position = as.character(Position))

	combined <- add_depth(DIR_depth_metrics, PATH_collective_depth_metrics, combined) # add depth information at these positions 
	combined <- subset_to_panel(PATH_bed, combined)
	write_csv(combined, "/groups/wyattgrp/users/amunzur/COMPOST_BIN/mutato_after_filtering.csv")
	
	# combine with the basecounts df 
	basecounts_df <- parse_basecount_vcf(DIR_basecounts)
	combined <- left_join(combined, basecounts_df, by = c("Chrom", "Position", "Ref", "Alt", "Sample_name"))
	names(combined) <- gsub("\\.x", "_n", names(combined))
    names(combined) <- gsub("\\.y", "_t", names(combined))

	# final filtering 
	combined <- filter(combined, Alt_reads_t > 0)

	return(combined)

} # end of function

combine_and_save <- function(variants, PATH_SAVE_chip_variants){

	dir.create(dirname(PATH_SAVE_chip_variants))

	write_csv(variants, PATH_SAVE_chip_variants) # snv + indel, csv
	write_delim(variants, gsub(".csv", ".tsv", PATH_SAVE_chip_variants), delim = "\t") # snv + indel, tsv
}