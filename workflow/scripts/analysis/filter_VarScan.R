# This R script filters and plots VarScan2 output. 

library(tidyverse)
library(stringr)
library(argparse)

parser <- ArgumentParser(description = "Filter the variants called by VarScan based on various parameters")

# passed from Snakemake into the script
parser$add_argument('--THRESHOLD_ExAC_ALL', metavar='VARIABLE', type='double', help="Variants with ExAC values more than this will be filtered out, they might be germline variants.")
parser$add_argument('--VALUE_Func_refGene', metavar='VARIABLE', type='character', help="Whether to exclude intronic variants or not.")
parser$add_argument('--THRESHOLD_VarFreq', metavar='VARIABLE', type='double', help="Minimum VAF we allow. If it above, will be filtered out. Give a fraction, not percentage.")
parser$add_argument('--THRESHOLD_Reads2', metavar='VARIABLE', type='double', help="Minimum number of mutant reads. If a variant has less, will be filtered out.")
parser$add_argument('--THRESHOLD_VAF_bg_ratio', metavar='VARIABLE', type='double', help="The ratio of VAF and the bg error rate at that position must be more than this.")
parser$add_argument('--DIR_varscan', metavar='DIRECTORY', type='character', help="Directory where the VarScan outputs are located.")
parser$add_argument('--PATH_bg', metavar='FILE', type='character', help="{Path to the background error rate.}")
parser$add_argument('--PATH_bets', metavar='FILE', type='character', help="Path to the betastasis table from this cohort.")
parser$add_argument('--PATH_bed', metavar='FILE', type='character', help="Path to the bed file for the panel.")
parser$add_argument('--DIR_depth_metrics', metavar='DIRECTORY', type='character', help="Dir where we keep the depths at individivual base positions.")
parser$add_argument('--PATH_SAVE_chip_variants', metavar='FILE', type='character', help="Path to save the cleaned up chip variants to after filtering.")
parser$add_argument('--PATH_validated_variants', metavar='FILE', type='character', help="Path to save the validated variants file with information about whether we detected them all or not.")

args <- parser$parse_args()

# DIR_figures <- "/groups/wyattgrp/users/amunzur/chip_project/figures/main_figures/new_chip_panel_CP"
# PATH_SAVE_chip_variants <- "/groups/wyattgrp/users/amunzur/chip_project/VarScan2_results/WBC_only/new_chip_panel_CP/snv/exonic.csv"

# filter anno output from all samples based on exac score and combine into one big df
return_anno_output <- function(PATH_annovar) {

	df_main <- as.data.frame(read_delim(PATH_annovar, delim = "\t"))

	df <- df_main %>% 
		mutate(Sample_name = gsub(".hg38_multianno.txt", "", basename(PATH_annovar))) %>%
		select(Sample_name, Chr, Start, Ref, Alt, Func.refGene, Gene.refGene, Func.knownGene, Gene.knownGene, AAChange.refGene, ExAC_ALL, gnomAD_exome_ALL) 

	return(df)

}

# filter anno output from all samples based on exac score and combine into one big df
return_varscan_output <- function(PATH_varscan) {

	print(PATH_varscan)

	df_main <- as.data.frame(read_delim(PATH_varscan, delim = "\t")) 
	# names(df_main) <- c("Chr", "Start", "Ref", "Cons", "Reads1", "Reads2", "VarFreq", "Strands1", "Strands2", "Qual1", 
						# "Qual2", "Pvalue", "MapQual1", "MapQual2", "Reads1Plus", "Reads1Minus", "Reads2Plus", "Reads2Minus", "Alt", "Sample_name")

	df <- df_main %>%
			mutate(Sample_name = gsub(".vcf", "", basename(PATH_varscan))) %>%
			rename(Chr = Chrom, Start = Position, Alt = VarAllele)

	return(df)

}

# do a merge based on column to combine metadata 
combine_anno_varscan <- function(DIR_varscan_snv, DIR_varscan_indel, ANNOVAR_snv_output, ANNOVAR_indel_output, variant_type) {

	# determine which dir to scan based on the cariant_type given
	if (variant_type == "snv") {
		DIR_varscan <- DIR_varscan_snv
		DIR_ANNOVAR <- ANNOVAR_snv_output} 
	else {DIR_varscan <- DIR_varscan_indel
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
					mutate(ExAC_ALL = as.numeric(gsub(".", 0, ExAC_ALL)), 
						gnomAD_exome_ALL = as.numeric(gsub(".", 0, gnomAD_exome_ALL)))

	return(combined)

}

add_bg_error_rate <- function(variants_df, PATH_bg) {

	# modify the vars df
	variants_df <- variants_df %>% mutate(error_type = paste0("mean_error", variants_df$Ref, "to", variants_df$Alt))
	variants_df$error_type[grep("\\+", variants_df$error_type)] <- "mean_errorins"
	variants_df$error_type[grep("\\-", variants_df$error_type)] <- "mean_errordel"

	# modify the bg error rate df
	bg <- read_delim(PATH_bg, delim = "\t")
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
	variant_df$Position_in_bets <- variant_df$Position %in% bets$POSITION

	return(variant_df)
}

# this function compares variants to variants jack found. need a path to jack
compare_with_jacks_figure <- function(PATH_to_jack, variant_df){

	validated_vars <- read.delim(PATH_to_jack, header = FALSE)
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
					DIR_varscan_snv, 
					DIR_varscan_indel, 
					ANNOVAR_snv_output, 
					ANNOVAR_indel_output,
					PATH_bg,
					PATH_bets,
					PATH_bed,
					DIR_depth_metrics,
					variant_type){

	combined <- combine_anno_varscan(DIR_varscan_snv, DIR_varscan_indel, ANNOVAR_snv_output, ANNOVAR_indel_output, variant_type) # add annovar annots to the varscan outputs

	# add a new col to show if we have an indel or an snv
	if (variant_type == "indel") {combined$variant <- ifelse(unlist(str_match(combined$Alt, "\\-|\\+")[, 1]) == "-", "deletion", "insertion")} else {
		combined$variant <- "snv"
	}

	# add patient id
	x <- str_split(combined$Sample_name, "-") # split the sample name
	combined$patient_id <- paste(lapply(x, "[", 1), lapply(x, "[", 2), lapply(x, "[", 3), sep = "-") # paste 2nd and 3rd elements to create the sample name 

	combined <- add_bg_error_rate(combined, PATH_bg) # background error rate
	combined_not_intronic <- combined %>%
						mutate(VarFreq = as.numeric(gsub("%", "", VarFreq))/100, 
								Total_reads = Reads1 + Reads2, 
								VAF_bg_ratio = VarFreq/error_rate) %>%
						select(Sample_name, patient_id, Chr, Start, Ref, Alt, VarFreq, Reads1, Reads2, Func.refGene, Gene.refGene, AAChange.refGene, ExAC_ALL, variant, error_rate, VAF_bg_ratio) %>%
						filter(ExAC_ALL <= THRESHOLD_ExAC_ALL, 
								Func.refGene != VALUE_Func_refGene,
								VarFreq <= THRESHOLD_VarFreq, 
								Reads2 >= THRESHOLD_Reads2, 
								VAF_bg_ratio >= THRESHOLD_VAF_bg_ratio) # vaf should be at least 15 times more than the bg error rate

	# a common naming convention i will be sticking to from now on
	names(combined_not_intronic) <- c("Sample_name", "Patient_ID", "Chrom", "Position", "Ref", "Alt", "VAF", "Ref_reads", "Alt_reads", "Function", "Gene", "AAchange", "ExAC_ALL", "Variant", "Error_rate", "VAF_bg_ratio")
	dedup <- compare_with_bets(PATH_bets, combined_not_intronic)
	dedup <- distinct(combined_not_intronic, Chrom, Position, Ref, Alt, .keep_all = TRUE) # remove duplicated variants
	
	dedup <- subset_to_panel(PATH_bed, dedup) # subset to panel 
	dedup <- add_depth(DIR_depth_metrics, dedup) # add depth information at these positions
	
	return(dedup)

}

combine_and_save <- function(snv, indel, PATH_validated_variants, PATH_SAVE_chip_variants){

	indel_igv <- indel %>% mutate(Position = Position + 1) 	# make a separate df just for igv

	variants_chip <- as.data.frame(rbind(snv, indel))
	variants_chip_igv <- rbind(snv, indel_igv)

	# variants_documentation <- add_documentation(THRESHOLD_ExAC_ALL, VALUE_Func_refGene, THRESHOLD_VarFreq, THRESHOLD_Reads2, THRESHOLD_VAF_bg_ratio)
	# comment(variants_chip) <- variants_documentation

	# PATH_to_jack <- "/groups/wyattgrp/users/amunzur/chip_project/validated_kidney_variants/chip_muts_locations.tsv" # compare with previously validated variants
	validated_vars <- compare_with_jacks_figure(PATH_validated_variants, variants_chip)
	variants_chip$detected <- variants_chip$Position %in% validated_vars$Position # add a new col to show if the same var was detected in jacks figure

	# PATH_SAVE_chip_variants <- "/groups/wyattgrp/users/amunzur/chip_project/variant_lists/chip_variants.csv"
	write_csv(variants_chip, PATH_SAVE_chip_variants) # snv + indel, csv
	write_delim(variants_chip, gsub(".csv", ".tsv", PATH_SAVE_chip_variants), delim = "\t") # snv + indel, tsv
	write_delim(variants_chip_igv, gsub("chip_variants.csv", "chip_variants_igv.tsv", PATH_SAVE_chip_variants), delim = "\t") # snv + indel for igv

	write_csv(validated_vars, gsub("chip_variants.csv", "validated_variants.csv", PATH_SAVE_chip_variants)) # jack df as csv
	write_delim(validated_vars, gsub("chip_variants.csv", "validated_variants.tsv", PATH_SAVE_chip_variants), delim = "\t") # jack df as tsv
}

snv <- MAIN(THRESHOLD_ExAC_ALL = args$THRESHOLD_ExAC_ALL, 
				VALUE_Func_refGene = args$VALUE_Func_refGene, 
				THRESHOLD_VarFreq = args$THRESHOLD_VarFreq, 
				THRESHOLD_Reads2 = args$THRESHOLD_Reads2, 
				THRESHOLD_VAF_bg_ratio = args$THRESHOLD_VAF_bg_ratio, 
				DIR_varscan_snv = args$DIR_varscan_snv,
				DIR_varscan_indel = args$DIR_varscan_indel,
				ANNOVAR_snv_output = args$ANNOVAR_snv_output,
				ANNOVAR_indel_output = args$ANNOVAR_indel_output,
				PATH_bg = args$PATH_bg,
				PATH_bets = args$PATH_bets, 
				PATH_bed = args$PATH_bed,
				DIR_depth_metrics = args$DIR_depth_metrics, 
				"snv")

indel <- MAIN(THRESHOLD_ExAC_ALL = args$THRESHOLD_ExAC_ALL, 
				VALUE_Func_refGene = args$VALUE_Func_refGene, 
				THRESHOLD_VarFreq = args$THRESHOLD_VarFreq, 
				THRESHOLD_Reads2 = args$THRESHOLD_Reads2, 
				THRESHOLD_VAF_bg_ratio = args$THRESHOLD_VAF_bg_ratio, 
				DIR_varscan_snv = args$DIR_varscan_snv,
				DIR_varscan_indel = args$DIR_varscan_indel,
				ANNOVAR_snv_output = args$ANNOVAR_snv_output,
				ANNOVAR_indel_output = args$ANNOVAR_indel_output,
				PATH_bg = args$PATH_bg,
				PATH_bets = args$PATH_bets, 
				PATH_bed = args$PATH_bed,
				DIR_depth_metrics = args$DIR_depth_metrics, 
				"indel")

combine_and_save(snv = snv,
					indel = indel, 
					PATH_validated_variants = args$PATH_validated_variants, 
					PATH_SAVE_chip_variants = args$PATH_SAVE_chip_variants)


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

