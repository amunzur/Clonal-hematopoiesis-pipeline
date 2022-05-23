return_varscan_output <- function(PATH_varscan) {

	df_main <- as.data.frame(suppress_messages(read_delim(PATH_varscan, delim = "\t"))) 
	df <- df_main %>%
			mutate(Sample_name = gsub(".vcf", "", basename(PATH_varscan))) %>%
			rename(Chr = Chrom, Start = Position, Alt = VarAllele)
			
	return(df)

}

# This helps put the varscan output for indels in the same format as vardict so that they are consistent. 
# This also helps when merging the annovar outputs with the varscan indels
# Very similar to the make anno indel function.
modify_varscan_output <- function(varscan_df, variant_type) {

	if (variant_type == "indel") {
		
		varscan_df <- varscan_df %>%
		select(-Ref) %>% # we remove this, to regenerate it based on the Cons column
		rename(Ref = Cons) %>%
		mutate(
			variant = case_when(
				grepl("-", Ref) ~ "deletion", 
				grepl("+", Ref) ~ "insertion",
				TRUE ~ "FUCK"), # a little moment of frustration
			Start = ifelse(variant == "deletion", as.numeric(Start) + 1, as.numeric(Start)), # if DELETION, add 1 to the position. INS positions stay the way they are.
			Ref = basename(as.character(Ref)), # to remove everything before the slash
			Ref = gsub("[[:punct:]]", "", Ref), # consider everything but the punctuation, this one removes + and - signs from the ALT allele
			Alt = Ref,
			Alt = ifelse(variant == "deletion", "-", Alt), 
			Ref = ifelse(variant == "insertion", "-", Ref), 
			End = Start + nchar(Ref) - 1, 
			Start = as.character(Start))

	} else {

		varscan_df$Start <- as.character(varscan_df$Start)

	}

	return(varscan_df)
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
	varscan_df <- do.call(rbind, varscan_df_list)

	varscan_df <- modify_varscan_output(varscan_df, variant_type)

	combined <- left_join(varscan_df, anno_df, by = c("Sample_name", "Chr", "Start", "Ref", "Alt")) %>%
					mutate(ExAC_ALL = replace_na(as.numeric(ExAC_ALL), 0), 
						gnomAD_exome_ALL = replace_na(as.numeric(gnomAD_exome_ALL), 0))
					
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

	var_pos <- variants_df$Start
	deletion_idx <- which(variants_df$variant == "deletion")
	
	# if (length(deletion_idx) > 0) {
	# 	var_pos[deletion_idx] <- as.numeric(var_pos[deletion_idx]) + 1 # add 1 to the position of all the deletions
	# 	variants_df$Start <- var_pos}

	# now subset the bg based on the chr and position from the variants df 
	bg$pos <- as.character(bg$pos)
	# bg <- left_join(variants_df[, c("Chr", "Start")], bg, by = c("Chr" = "chrom", "Start" = "pos")) # to filter bg with the chrom and pos of interest
	variants_df <- left_join(variants_df, bg, by = c("Chr" = "chrom", "Start" = "pos", "error_type" = "error_type")) # to add the error rates from bg to variants df

	# subtract 1 from the positions after merging
	if (length(deletion_idx) > 0) {
		deletion_idx <- which(variants_df$variant == "deletion")
		variants_df$Start[deletion_idx] <- as.numeric(variants_df$Start[deletion_idx]) - 1}

	return(variants_df)
}

# Do Fisher's test to calculate the strand bias. This function should be applied across all rows with an apply function. 
evaluate_strandBias <- function(variants_df){

	message("Evaluating strand bias.")

	col_ids_list <- grep(c("Reads1Plus|Reads2Plus|Reads1Minus|Reads2Minus"), names(variants_df))

	minifunction <- function(some_row, col_ids_list) {
		
		dat <- matrix(c(as.numeric(some_row[col_ids_list[1]]), #Reads1Plus
						as.numeric(some_row[col_ids_list[2]]), #Reads2Plus
						as.numeric(some_row[col_ids_list[3]]), #Reads1Minus
						as.numeric(some_row[col_ids_list[4]])),#Reads2Minus
						nrow = 2, ncol = 2, byrow = FALSE)

		dimnames(dat) <- list(c("Ref", "Alt"), c("Forward", "Reverse"))

		fishers_pval <- fisher.test(dat)$p.value
		odds_ratio <- oddsratio(dat + 1)$measure[2]

		return(list(fishers_pval = fishers_pval, odds_ratio = odds_ratio))}

	values_list <- apply(variants_df, 1, minifunction, col_ids_list) # contains both the pval from fishers test and the odds ratio
	variants_df$StrandBias_Fisher_pVal <- unname(unlist(values_list)[seq(1, length(unlist(values_list)), 2)]) # select every odd element starting from 0 - fishers exact
	variants_df$StrandBias_OddsRatio <- unname(unlist(values_list)[seq(2, length(unlist(values_list)), 2)]) # select every even element starting from 2 - odds ratio

	return(variants_df)

}

# main function to run everything
MAIN <- function(cohort_name, 
					THRESHOLD_ExAC_ALL, 
					VALUE_Func_refGene, 
					THRESHOLD_VarFreq, 
					THRESHOLD_Reads2, 
					THRESHOLD_VAF_bg_ratio, 
					DIR_varscan_snv, 
					DIR_varscan_indel, 
					ANNOVAR_snv_output, 
					ANNOVAR_indel_output,
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
					var_type,
					variant_caller){

	variant_caller = variant_caller # so that the function doesnt complain this argument isn't used
	combined <- combine_anno_varscan(DIR_varscan_snv, DIR_varscan_indel, ANNOVAR_snv_output, ANNOVAR_indel_output, var_type) # add annovar annots to the varscan outputs

	# add a new col to show if we have an snv, for indels we add this annotation in the modify_varscan_output function
	if (var_type == "snv") {combined$variant <- "snv"} 

	combined <- add_patient_id(combined, cohort_name)
	combined <- add_bg_error_rate(combined, bg) 
	combined <- add_sample_type(combined)
	combined$Cohort_name <- cohort_name

	combined <- combined %>%
						mutate(VarFreq = as.numeric(gsub("%", "", VarFreq))/100, 
								Total_reads = Reads1 + Reads2, 
								VAF_bg_ratio = VarFreq/error_rate) %>%
						filter(ExAC_ALL <= THRESHOLD_ExAC_ALL, 
								Func.refGene != VALUE_Func_refGene,
								VAF_bg_ratio >= THRESHOLD_VAF_bg_ratio) # vaf should be at least 15 times more than the bg error rate
	
	combined <- evaluate_strandBias(combined) # add strand bias
	combined <- add_AAchange_effect(combined) # Add protein annot and the effect of the mutations as two separate columns 

	combined <- combined %>%
						filter(StrandBias_Fisher_pVal > 0.05) %>%
						select(Sample_name, Sample_type, patient_id, Cohort_name, Chr, Start, Ref, Alt, VarFreq, Reads1, Reads2, StrandBias_Fisher_pVal, StrandBias_OddsRatio, Reads1Plus, Reads1Minus, Reads2Plus, Reads2Minus, Func.refGene, Gene.refGene, AAChange.refGene, Protein_annotation, Effects, ExAC_ALL, variant, error_rate, VAF_bg_ratio, Total_reads)

	names(combined) <- c("Sample_name", "Sample_type", "Patient_ID", "Cohort_name", "Chrom", "Position", "Ref", "Alt", "VAF", "Ref_reads", "Alt_reads", "StrandBias_Fisher_pVal", "StrandBias_OddsRatio", "REF_Fw", "REF_Rv", "ALT_Fw", "ALT_Rv", "Function", "Gene", "AAchange", "Protein_annotation", "Effects", "ExAC_ALL", "Variant", "Error_rate", "VAF_bg_ratio", "Total_reads")

	# add for splicing variants, make sure both the "Function" and "Effects" column have the string splicing
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
	combined <- compare_with_bets(PATH_bets_somatic, PATH_bets_germline, PATH_panel_genes, combined) # adds 3 new columns

	# add an extra col for alerting the user if the variant isn't found, despite gene being in the bets
	combined <- combined %>% mutate(Status = case_when(
										(In_germline_bets == FALSE & In_panel == TRUE) ~ "ALERT", 
										(In_germline_bets == TRUE & In_panel == TRUE) ~ "Great",
										(In_germline_bets == TRUE & In_panel == FALSE) ~ "Error",
										TRUE ~ "OK"), 
								Position = as.numeric(Position))

	# Check for duplicated rows just before returning the object.
	combined <- check_duplicated_rows(combined, TRUE)

	return(combined)
} # end of function

combine_and_save <- function(snv, indel, PATH_validated_variants, PATH_SAVE_chip_variants){

	variants_chip <- as.data.frame(rbind(snv, indel))
	validated_vars <- compare_with_jacks_figure(PATH_validated_variants, variants_chip)

	dir.create(dirname(PATH_SAVE_chip_variants))
	
	write_csv(variants_chip, PATH_SAVE_chip_variants) # snv + indel, csv
	write_delim(variants_chip, gsub(".csv", ".tsv", PATH_SAVE_chip_variants), delim = "\t") # snv + indel, tsv

	write_csv(validated_vars, gsub("chip_variants.csv", "validated_variants.csv", PATH_SAVE_chip_variants)) # jack df as csv
	write_delim(validated_vars, gsub("chip_variants.csv", "validated_variants.tsv", PATH_SAVE_chip_variants), delim = "\t") # jack df as tsv

	return(variants_chip)
}