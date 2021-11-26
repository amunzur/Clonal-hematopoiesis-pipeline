return_varscan_output <- function(PATH_varscan) {

	df_main <- as.data.frame(read_delim(PATH_varscan, delim = "\t")) 
	df <- df_main %>%
			mutate(Sample_name = gsub(".vcf", "", basename(PATH_varscan))) %>%
			rename(Chr = Chrom, Start = Position, Alt = VarAllele) %>%
			filter(!str_detect(Alt, 'N')) 	# Drop variants if the called variant has an "N" in it. 

	return(df)

}

# do a merge based on column to combine metadata 
combine_anno_varscan <- function(DIR_varscan_snv, DIR_varscan_indel, ANNOVAR_snv_output, ANNOVAR_indel_output, variant_type) {

	# determine which dir to scan based on the cariant_type given
	if (variant_type == "snv") {
		DIR_varscan <- DIR_varscan_snv
		DIR_ANNOVAR <- ANNOVAR_snv_output
	} else {DIR_varscan <- DIR_varscan_indel
		DIR_ANNOVAR <- ANNOVAR_indel_output}

	anno_df_list <- lapply(as.list(list.files(DIR_ANNOVAR, full.names = TRUE, pattern = "\\.hg38_multianno.txt$")), return_anno_output)
	anno_df <- as.data.frame(do.call(rbind, anno_df_list)) %>%
				mutate(Sample_name = gsub(".hg38_multianno.txt", "", Sample_name), 
						Start = as.character(Start)) 

	varscan_df_list <- lapply(as.list(list.files(DIR_varscan, full.names = TRUE, pattern = "\\.vcf$")), return_varscan_output)
	varscan_df <- do.call(rbind, varscan_df_list) %>%
				mutate(Start = as.character(Start)) 

	# if (variant_type == "indel"){
	# 	varscan_df$Alt <- gsub("\\+", "", varscan_df$Alt)
	# 	varscan_df$Alt[grep("-", varscan_df$Alt)] <- "-"}

	combined <- left_join(varscan_df, anno_df, by = c("Sample_name", "Chr", "Start", "Ref", "Alt")) %>%
					mutate(ExAC_ALL = as.numeric(gsub("\\.", 0, ExAC_ALL)), 
						gnomAD_exome_ALL = as.numeric(gsub("\\.", 0, gnomAD_exome_ALL)))

	return(combined)

}

add_bg_error_rate <- function(variants_df, bg) {

	# modify the vars df
	variants_df <- variants_df %>% mutate(error_type = paste0("mean_error", variants_df$Ref, "to", variants_df$Alt))
	variants_df$error_type[grep("\\+", variants_df$error_type)] <- "mean_errorins"
	variants_df$error_type[grep("\\-", variants_df$error_type)] <- "mean_errordel"

	# modify the bg error rate df
	bg$chrom <- paste0("chr", bg$chrom)
	bg <- gather(bg, "error_type", "error_rate", starts_with("mean_error"))

	# now if we have a deletion, we want to recover the pos + 1 in the bg file, this is actual place where deletion begins 
	original_pos <- variants_df$Start

	var_pos <- variants_df$Start
	deletion_idx <- which(variants_df$variant == "deletion")
	
	if (length(deletion_idx) > 0) {
		var_pos[deletion_idx] <- as.numeric(var_pos[deletion_idx]) + 1 # add 1 to the position of all the deletions
		variants_df$Start <- var_pos}

	bg <- subset(bg, pos %in% variants_df$Start) # subset to the positions we have in the variants so that the df is smaller and more managable

	# this merge adds an exta col with the error rate 
	bg$pos <- as.character(bg$pos)
	variants_df <- left_join(variants_df, bg, by = c("Start" = "pos", "error_type" = "error_type"))

	# the merge above adds two extra cols we dont want, remove those 
	variants_df$chrom <- NULL
	variants_df$ref <- NULL

	# subtract 1 from the positions after merging
	if (length(deletion_idx) > 0) {
		deletion_idx <- which(variants_df$variant == "deletion")
		variants_df$Start[deletion_idx] <- as.numeric(variants_df$Start[deletion_idx]) - 1}

	return(variants_df)
}

# Do Fisher's test to calculate the strand bias. This function should be applied across all rows with an apply function. 
evaluate_strandBias <- function(variants_df){

	message("Evaluating strand bias.")

	minifunction <- function(some_row) {
		
		dat <- matrix(c(as.numeric(some_row[15]), 
						as.numeric(some_row[17]), 
						as.numeric(some_row[16]), 
						as.numeric(some_row[18])), nrow = 2, ncol = 2, byrow = FALSE)

		dimnames(dat) <- list(c("Ref", "Alt"), c("Forward", "Reverse"))

		fishers_pval <- fisher.test(dat)$p.value
		odds_ratio <- oddsratio(dat + 1)$measure[2]

		return(list(fishers_pval = fishers_pval, odds_ratio = odds_ratio))}

	values_list <- apply(variants_df, 1, minifunction) # contains both the pval from fishers test and the odds ratio
	variants_df$StrandBias_Fisher_pVal <- unname(unlist(values_list)[seq(1, length(unlist(values_list)), 2)]) # select every odd element starting from 0 - fishers exact
	variants_df$StrandBias_OddsRatio <- unname(unlist(values_list)[seq(2, length(unlist(values_list)), 2)]) # select every even element starting from 2 - odds ratio

	return(variants_df)

}

# main function to run everything
MAIN <- function(THRESHOLD_ExAC_ALL, 
					VALUE_Func_refGene, 
					THRESHOLD_VarFreq, 
					THRESHOLD_Reads2, 
					THRESHOLD_VAF_bg_ratio, 
					DIR_varscan_snv, 
					DIR_varscan_indel, 
					ANNOVAR_snv_output, 
					ANNOVAR_indel_output,
					bg,
					PATH_bets,
					PATH_bed,
					DIR_depth_metrics,
					PATH_collective_depth_metrics,
					DIR_finland_bams,
					DIR_temp,
					DIR_tnvstats,
					PATH_filter_tnvstats_script,
					variant_type,
					variant_caller){

	combined <- combine_anno_varscan(DIR_varscan_snv, DIR_varscan_indel, ANNOVAR_snv_output, ANNOVAR_indel_output, variant_type) # add annovar annots to the varscan outputs

	# add a new col to show if we have an indel or an snv
	if (variant_type == "indel") {combined$variant <- ifelse(unlist(str_match(combined$Alt, "\\-|\\+")[, 1]) == "-", "deletion", "insertion")} else {
		combined$variant <- "snv"
	}

	# add patient id
	x <- str_split(combined$Sample_name, "-") # split the sample name
	combined$patient_id <- paste(lapply(x, "[", 1), lapply(x, "[", 2), lapply(x, "[", 3), sep = "-") # paste 2nd and 3rd elements to create the sample name 

	combined <- add_bg_error_rate(combined, bg) # background error rate

	# some filtering based on some criteria
	combined_not_intronic <- combined %>%
						mutate(VarFreq = as.numeric(gsub("%", "", VarFreq))/100, 
								Total_reads = Reads1 + Reads2, 
								VAF_bg_ratio = VarFreq/error_rate) %>%
						filter(ExAC_ALL <= THRESHOLD_ExAC_ALL, 
								Func.refGene != VALUE_Func_refGene,
								VAF_bg_ratio >= THRESHOLD_VAF_bg_ratio) # vaf should be at least 15 times more than the bg error rate
	
	combined_not_intronic <- evaluate_strandBias(combined_not_intronic) # add strand bias
	combined_not_intronic <- add_AAchange_effect(combined_not_intronic) # Add protein annot and the effect of the mutations as two separate columns 

	combined_not_intronic <- combined_not_intronic %>%
						filter(StrandBias_Fisher_pVal > 0.05) %>%
						select(Sample_name, patient_id, Chr, Start, Ref, Alt, VarFreq, Reads1, Reads2, StrandBias_Fisher_pVal, StrandBias_OddsRatio, Reads1Plus, Reads1Minus, Reads2Plus, Reads2Minus, Func.refGene, Gene.refGene, AAChange.refGene, Protein_annotation, Effects, ExAC_ALL, variant, error_rate, VAF_bg_ratio, Total_reads)

	names(combined_not_intronic) <- c("Sample_name", "Patient_ID", "Chrom", "Position", "Ref", "Alt", "VAF", "Ref_reads", "Alt_reads", "StrandBias_Fisher_pVal", "StrandBias_OddsRatio", "REF_Fw", "REF_Rv", "ALT_Fw", "ALT_Rv", "Function", "Gene", "AAchange", "Protein_annotation", "Effects", "ExAC_ALL", "Variant", "Error_rate", "VAF_bg_ratio", "Total_reads")

	# add for splicing variants, make sure both the "Function" and "Effects" column have the string splicing
	idx <- grep("splicing", combined_not_intronic$Function)
	combined_not_intronic$Effects[idx] <- "splicing"

	# remove duplicated vars 
	dedup <- distinct(combined_not_intronic, Chrom, Position, Ref, Alt, .keep_all = TRUE) # remove duplicated variants
	
	# now filtering based on vaf, read support and depth etc. 
	dedup <- dedup %>%
				filter((Total_reads >= 1000 & VAF >= 0.005) | 
						(Total_reads <= 1000 & Alt_reads >= 5), 
						VAF < THRESHOLD_VarFreq)

	dedup <- subset_to_panel(PATH_bed, dedup) # subset to panel 
	dedup <- add_depth(DIR_depth_metrics, PATH_collective_depth_metrics, dedup) # add depth information at these positions
	dedup <- compare_with_bets(PATH_bets, dedup) # indicate whether the variant has been found previously and whether the gene where the variant is in exists in the bets table.

	# add an extra col for alerting the user if the variant isn't found, despite gene being in the bets
	dedup <- dedup %>% mutate(Status = case_when(
										(detected == FALSE & Gene_in_bets == TRUE) ~ "ALERT", 
										(detected == TRUE & Gene_in_bets == TRUE) ~ "Great",
										(detected == TRUE & Gene_in_bets == FALSE) ~ "Error",
										TRUE ~ "OK"), 
								Position = as.numeric(Position))

	# add the sample ID from the finland bams
	finland_sample_IDs <- gsub(".bam", "", grep("*.bam$", list.files(DIR_finland_bams), value = TRUE))
	dedup$Sample_name_finland <- unlist(lapply(as.list(gsub("GUBB", "GU", dedup$Patient_ID)), function(x) grep(x, finland_sample_IDs, value = TRUE)))
	dedup <- dedup %>% relocate(Sample_name_finland, .after = Sample_name)

	# write_csv(dedup, "/groups/wyattgrp/users/amunzur/pipeline/results/temp/snv_varscan_dedup.csv")
	# write_csv(dedup, "/groups/wyattgrp/users/amunzur/pipeline/results/temp/indel_varscan_dedup.csv")
	# dedup <- read_csv("/groups/wyattgrp/users/amunzur/pipeline/results/temp/snv_varscan_dedup.csv")
	# dedup <- read_csv("/groups/wyattgrp/users/amunzur/pipeline/results/temp/indel_varscan_dedup.csv")

	message("Filtering tnvstats right now.")
	filtered_tnvstats <- filter_tnvstats_by_variants(
							variants_df = dedup, 
							DIR_tnvstats = DIR_tnvstats, 
							PATH_temp = file.path(DIR_temp, paste(variant_caller, paste0(variant_type, ".tsv"), sep = "_")), 
							PATH_filter_tnvstats_script = PATH_filter_tnvstats_script, 
							identifier = paste(variant_caller, variant_type, sep = "_"))

	dedup <- add_finland_readcounts(dedup, filtered_tnvstats)

	return(dedup)
}

combine_and_save <- function(snv, indel, PATH_validated_variants, PATH_SAVE_chip_variants){

	variants_chip <- as.data.frame(rbind(snv, indel))

	# variants_documentation <- add_documentation(THRESHOLD_ExAC_ALL, VALUE_Func_refGene, THRESHOLD_VarFreq, THRESHOLD_Reads2, THRESHOLD_VAF_bg_ratio)
	# comment(variants_chip) <- variants_documentation

	# PATH_to_jack <- "/groups/wyattgrp/users/amunzur/chip_project/validated_kidney_variants/chip_muts_locations.tsv" # compare with previously validated variants
	validated_vars <- compare_with_jacks_figure(PATH_validated_variants, variants_chip)
	# variants_chip$detected <- variants_chip$Position %in% validated_vars$Position # add a new col to show if the same var was detected in jacks figure

	# PATH_SAVE_chip_variants <- "/groups/wyattgrp/users/amunzur/chip_project/variant_lists/chip_variants.csv"
	write_csv(variants_chip, PATH_SAVE_chip_variants) # snv + indel, csv
	write_delim(variants_chip, gsub(".csv", ".tsv", PATH_SAVE_chip_variants), delim = "\t") # snv + indel, tsv

	write_csv(validated_vars, gsub("chip_variants.csv", "validated_variants.csv", PATH_SAVE_chip_variants)) # jack df as csv
	write_delim(validated_vars, gsub("chip_variants.csv", "validated_variants.tsv", PATH_SAVE_chip_variants), delim = "\t") # jack df as tsv

	return(variants_chip)
}