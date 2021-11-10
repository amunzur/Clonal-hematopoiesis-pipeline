# This R script filters and plots VarScan2 output. 

library(tidyverse)
library(stringr)

# DIR_figures <- "/groups/wyattgrp/users/amunzur/chip_project/figures/main_figures/new_chip_panel_CP"
# PATH_SAVE_chip_variants <- "/groups/wyattgrp/users/amunzur/chip_project/VarScan2_results/WBC_only/new_chip_panel_CP/snv/exonic.csv"

# filter anno output from all samples based on exac score and combine into one big df
return_anno_output <- function(DIR_ANNOVAR) {

	df_main <- as.data.frame(read_delim(DIR_ANNOVAR, delim = "\t"))

	df <- df_main %>% 
		mutate(Sample_name = gsub(".hg38_multianno.txt", "", basename(DIR_ANNOVAR))) %>%
		select(Sample_name, Chr, Start, Ref, Alt, Func.refGene, Gene.refGene, Func.knownGene, Gene.knownGene, AAChange.refGene, ExAC_ALL, gnomAD_exome_ALL) 

	return(df)

}

# filter anno output from all samples based on exac score and combine into one big df
return_vardict_output <- function(DIR_vardict) {

	df_main <- as.data.frame(read_delim(DIR_vardict, delim = "\t")) 
	df <- df_main %>% mutate(Sample_name = gsub(".tsv", "", basename(DIR_vardict)))

	return(df)

}

# do a merge based on column to combine metadata 
combine_anno_varscan <- function(DIR_vardict, DIR_ANNOVAR) {

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
						ExAC_ALL = as.numeric(gsub(".", 0, ExAC_ALL)), 
						gnomAD_exome_ALL = as.numeric(gsub(".", 0, gnomAD_exome_ALL)))

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

# compare with the bets list, add an extra col to show if the vars we called also appear in the bets
compare_with_bets <- function(PATH_bets, variant_df){

	bets <- as.data.frame(read_delim(PATH_bets, delim = "\t"))
	variant_df$detected <- variant_df$Position %in% bets$POSITION

	# add a new col to indicate if the gene of interest in each row is found in the bets. 
	variant_df$Gene_in_bets <- variant_df$Gene %in% bets$GENE

	return(variant_df)
}

# this function compares variants to variants jack found. need a path to jack
compare_with_jacks_figure <- function(PATH_validated_variants, variant_df){

	validated_vars <- read.delim(PATH_validated_variants, header = FALSE, sep = "\t")
	names(validated_vars) <- c("Sample", "Gene", "Chrom", "Position", "Ref", "Alt")

	validated_vars$detected <- validated_vars$Position %in% variant_df$Position

	return(validated_vars)
}

# add the depth at each variant position
add_depth <- function(DIR_depth_metrics, variant_df) {

	file_name_list <- list.files(DIR_depth_metrics, full.names = TRUE)
	depth_file_list <- lapply(file_name_list, read.delim, header = FALSE) # load the related files to a list
	depth_file <- do.call(rbind, depth_file_list) # combine all depth files into one large one

	sample_names_repeated <- mapply(rep, unlist(lapply(file_name_list, basename)), unlist(lapply(depth_file_list, nrow)))
	sample_names_repeated <- unlist(lapply(sample_names_repeated, function(x) gsub(".txt", "", x)))
	depth_file$Sample_name <- sample_names_repeated # add the sample names so that we can do a join based on sample names with the variant df later on
	names(depth_file) <- c("Chrom", "Position", "Depth", "Sample_name")
	depth_file$Position <- as.character(depth_file$Position) # add the sample names so that we can do a join based on sample names with the variant df later on


	combined <- left_join(variant_df, depth_file, by = c("Sample_name", "Chrom", "Position"))

	return(combined)

}

# add_mutation_type <- function(variant_df){

# 	Mutation_type <- as.list(str_match(unlist(variant_df$AAchange), "p.[A-z][0-9][0-9[a-zA-Z0-9]][0-9[a-zA-Z0-9]][a-zA-Z0-9]*"))

# lapply(Mutation_type, substr, 3, 3)





# }

subset_to_panel <- function(PATH_bed, variant_df) {
	
	bed <- as.data.frame(read.delim(PATH_bed, header = FALSE, sep = "\t"))

	if (ncol(bed) == 4) {names(bed) <- c("chrom", "start", "stop", "gene")} else {names(bed) <- c("chrom", "start", "stop")}

	to_keep <- list() # list of positions to remove

	print("Subsetting to the panel.")
	i <- 1 
	print(i)
	while (i <= dim(variant_df)[1]){
		chrom_subsetted <- variant_df[i, 3] # pick the chrom we are at 
		location <- variant_df[i, 4] # pick the location we are at 
		bed_subsetted <- bed %>% filter(chrom == chrom_subsetted) # subset the bed by chrom
		message(c("position:", i))

		j <- 1
		while(j <= dim(bed_subsetted)[1]) {
			start_pos <- bed_subsetted[j, 2]
			end_pos <- bed_subsetted[j, 3]

			if (all(location >= start_pos, location <= end_pos)) {
				to_keep <- append(to_keep, i) # saving the row index of muts we are keeping
				break 
				print("Position in panel.")
				} else {j <- j + 1} # if the location isn't in the panel, go check out the next position in the bed file.
				# print(c(j, "j"))

		} # end of inner while loop

	i <- i + 1 # next identified variant

	} # end of outer while loop - looping through identified variants

	variant_df <- variant_df[unlist(to_keep), ]

	return(variant_df)

} # end of function

# This function generates a string of explanations to add to the final df as comments. 
add_documentation <- function(THRESHOLD_ExAC_ALL, VALUE_Func_refGene, THRESHOLD_VarFreq, THRESHOLD_Reads2, THRESHOLD_VAF_bg_ratio){

	docs <- paste(paste0("Exac score: less than or equal to ", THRESHOLD_ExAC_ALL),
				paste0("Function: not ", VALUE_Func_refGene),
				paste0("VAF: less than or equal to:", THRESHOLD_VarFreq),
				paste0("Number of mutant reads:", THRESHOLD_Reads2), 
				paste0("VAF_bg_ratio", THRESHOLD_VAF_bg_ratio), 
				sep = "\n")

	return(docs)

}

# main function to run everything
MAIN <- function(THRESHOLD_ExAC_ALL, 
					VALUE_Func_refGene, 
					THRESHOLD_VarFreq, 
					THRESHOLD_Reads2, 
					THRESHOLD_VAF_bg_ratio, 
					DIR_vardict, 
					DIR_ANNOVAR,
					bg,
					PATH_bets,
					PATH_bed,
					DIR_depth_metrics,
					variant_type){

	combined <- combine_anno_varscan(DIR_vardict, DIR_ANNOVAR) # add annovar annots to the varscan outputs

	# add patient id
	x <- str_split(combined$Sample_name, "-") # split the sample name
	combined$patient_id <- paste(lapply(x, "[", 1), lapply(x, "[", 2), lapply(x, "[", 3), sep = "-") # paste 2nd and 3rd elements to create the sample name 

	combined <- add_bg_error_rate(combined, bg) # background error rate
	combined_not_intronic <- combined %>%
						mutate(Total_reads = Reads1 + Reads2, 
								VAF_bg_ratio = VarFreq/error_rate) %>%
						select(Sample_name, patient_id, Chr, Start, Ref, Alt, VarFreq, Reads1, Reads2, Func.refGene, Gene.refGene, AAChange.refGene, ExAC_ALL, variant, error_rate, VAF_bg_ratio, Total_reads) %>%
						filter(ExAC_ALL <= THRESHOLD_ExAC_ALL, 
								Func.refGene != VALUE_Func_refGene,
								VAF_bg_ratio >= THRESHOLD_VAF_bg_ratio) # vaf should be at least 15 times more than the bg error rate

	# a common naming convention i will be sticking to from now on
	names(combined_not_intronic) <- c("Sample_name", "Patient_ID", "Chrom", "Position", "Ref", "Alt", "VAF", "Ref_reads", "Alt_reads", "Function", "Gene", "AAchange", "ExAC_ALL", "Variant", "Error_rate", "VAF_bg_ratio", "Total_reads")

	dedup <- distinct(combined_not_intronic, Chrom, Position, Ref, Alt, .keep_all = TRUE) # remove duplicated variants

	# now filtering based on vaf, read support and depth etc. 
	dedup <- dedup %>%
				filter((Total_reads >= 1000 & VAF >= 0.005) | 
						(Total_reads <= 1000 & Alt_reads >= 5), 
						VAF < THRESHOLD_VarFreq)

	dedup <- subset_to_panel(PATH_bed, dedup) # subset to panel
	dedup <- add_depth(DIR_depth_metrics, dedup) # add depth information at these positions 
	dedup <- compare_with_bets(PATH_bets, dedup)

	# add an extra col for alerting the user if the variant isn't found, despite gene being in the bets
	dedup <- dedup %>% mutate(Status = case_when(
										(detected == FALSE & Gene_in_bets == TRUE) ~ "ALERT", 
										(detected == TRUE & Gene_in_bets == TRUE) ~ "Great",
										(detected == TRUE & Gene_in_bets == FALSE) ~ "Error",
										TRUE ~ "OK"))

	return(dedup)

}

combine_and_save <- function(variants, PATH_validated_variants, PATH_SAVE_chip_variants){

	# make a separate df for igv because the indel positions need to be adjusted before we can take snapshots 
	idx <- which(variants$Variant == "deletion")
	if (length(idx) > 0) {
		variants_igv <- variants %>% mutate(Position = as.numeric(Position))
		variants_igv$Position[idx] <- variants_igv$Position[idx] + 1} else {variants_igv <- variants}

	validated_vars <- compare_with_jacks_figure(PATH_validated_variants, variants)
	# variants$detected <- variants$Position %in% validated_vars$Position # add a new col to show if the same var was detected in jacks figure

	# PATH_SAVE_chip_variants <- "/groups/wyattgrp/users/amunzur/chip_project/variant_lists/chip_variants.csv"
	write_csv(variants, PATH_SAVE_chip_variants) # snv + indel, csv
	write_delim(variants, gsub(".csv", ".tsv", PATH_SAVE_chip_variants), delim = "\t") # snv + indel, tsv
	write_delim(variants_igv, gsub("chip_variants.csv", "chip_variants_igv.tsv", PATH_SAVE_chip_variants), delim = "\t") # snv + indel for igv

	write_csv(validated_vars, gsub("chip_variants.csv", "validated_variants.csv", PATH_SAVE_chip_variants)) # jack df as csv
	write_delim(validated_vars, gsub("chip_variants.csv", "validated_variants.tsv", PATH_SAVE_chip_variants), delim = "\t") # jack df as tsv
}

THRESHOLD_ExAC_ALL <- 0.005
VALUE_Func_refGene <- "intronic"
THRESHOLD_VarFreq <- 0.30
THRESHOLD_Reads2 <- 5
THRESHOLD_VAF_bg_ratio <- 10
DIR_vardict <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/Vardict/new_chip_panel_reformatted"
DIR_ANNOVAR <- "/groups/wyattgrp/users/amunzur/pipeline/results/data/annovar_outputs/vardict/new_chip_panel"
PATH_bg <- "/groups/wyattgrp/users/amunzur/pipeline/resources/bg_error_rate/bg_error.tsv"
PATH_bets <- "/groups/wyattgrp/users/amunzur/pipeline/resources/betastasis/CLEANED_mutations_no_germline_filter.tsv"
PATH_bed  <- "/groups/wyattgrp/users/amunzur/pipeline/resources/panel/1000012543_CHIP_Design_selection_results_Version2/capture_targets.bed"
DIR_depth_metrics <- "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/depth/new_chip_panel"

PATH_validated_variants <- "/groups/wyattgrp/users/amunzur/pipeline/resources/validated_variants/chip_muts_locations.tsv"
PATH_SAVE_chip_variants <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/Vardict/finalized/chip_variants.csv"

bg <- read_delim(PATH_bg, delim = "\t")

variants <- MAIN(THRESHOLD_ExAC_ALL, 
				VALUE_Func_refGene,
				THRESHOLD_VarFreq,
				THRESHOLD_Reads2,
				THRESHOLD_VAF_bg_ratio,
				DIR_vardict,
				DIR_ANNOVAR,
				bg,
				PATH_bets,
				PATH_bed,
				DIR_depth_metrics)

combine_and_save(variants,
				PATH_validated_variants, 
				PATH_SAVE_chip_variants)

# variants_df <- read_csv(PATH_SAVE_chip_variants)
# write_csv(variants, "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/Vardict/finalized/vars.csv") # snv + indel, csv

# snv <- MAIN(THRESHOLD_ExAC_ALL = THRESHOLD_ExAC_ALL, 
# 				VALUE_Func_refGene = VALUE_Func_refGene, 
# 				THRESHOLD_VarFreq = THRESHOLD_VarFreq, 
# 				THRESHOLD_Reads2 = THRESHOLD_Reads2, 
# 				THRESHOLD_VAF_bg_ratio = THRESHOLD_VAF_bg_ratio, 
# 				DIR_varscan_snv = DIR_varscan_snv,
# 				DIR_varscan_indel = DIR_varscan_indel,
# 				ANNOVAR_snv_output = ANNOVAR_snv_output,
# 				ANNOVAR_indel_output = ANNOVAR_indel_output,
# 				PATH_bg = PATH_bg,
# 				PATH_bets = PATH_bets, 
# 				PATH_bed = PATH_bed,
# 				DIR_depth_metrics = DIR_depth_metrics, 
# 				"snv")

# indel <- MAIN(THRESHOLD_ExAC_ALL = THRESHOLD_ExAC_ALL, 
# 				VALUE_Func_refGene = VALUE_Func_refGene, 
# 				THRESHOLD_VarFreq = THRESHOLD_VarFreq, 
# 				THRESHOLD_Reads2 = THRESHOLD_Reads2, 
# 				THRESHOLD_VAF_bg_ratio = THRESHOLD_VAF_bg_ratio, 
# 				DIR_varscan_snv = DIR_varscan_snv,
# 				DIR_varscan_indel = DIR_varscan_indel,
# 				ANNOVAR_snv_output = ANNOVAR_snv_output,
# 				ANNOVAR_indel_output = ANNOVAR_indel_output,
# 				PATH_bg = PATH_bg,
# 				PATH_bets = PATH_bets, 
# 				PATH_bed = PATH_bed,
# 				DIR_depth_metrics = DIR_depth_metrics, 
# 				"indel")

