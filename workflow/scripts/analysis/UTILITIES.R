return_anno_output <- function(DIR_ANNOVAR) {

	df_main <- as.data.frame(read_delim(DIR_ANNOVAR, delim = "\t"))

	df <- df_main %>% 
		mutate(Sample_name = gsub(".hg38_multianno.txt", "", basename(DIR_ANNOVAR))) %>%
		select(Sample_name, Chr, Start, Ref, Alt, Func.refGene, Gene.refGene, Func.knownGene, Gene.knownGene, AAChange.refGene, ExAC_ALL, gnomAD_exome_ALL) 

	return(df)

}

# based on the cohort_name, add the patient id
add_patient_id <- function(variant_df, cohort_name){

	if (cohort_name == "new_chip_panel"){

		x <- str_split(variant_df$Sample_name, "-") # split the sample name
		variant_df$patient_id <- paste(lapply(x, "[", 1), lapply(x, "[", 2), lapply(x, "[", 3), sep = "-") # paste 2nd and 3rd elements to create the sample name 

	} else {

		x <- str_split(variant_df$Sample_name, "_gDNA|_WBC|_cfDNA")
		variant_df$patient_id <- unlist(lapply(x, "[", 1))
	}

	return(variant_df)

}

# add a new column to indicate sample type, WBC or tumor. 
add_sample_type <- function(variant_df){

	x <- str_split(variant_df$Sample_name, "_") # splitted strings as a list
	y <- unlist(lapply(x, "[", 2)) # gDNA or cfDNA

	my_map <- c("gDNA" = "WBC", "gdna" = "WBC", "cfDNA" = "Tumor")
	variant_df$Sample_type <- unname(my_map[y])

	return(variant_df)

}

# compare with the somatic bets and germline bets list, add an extra col to show if the vars we called also appear in either one of them bets
compare_with_bets <- function(PATH_bets_somatic, PATH_bets_germline, PATH_panel_genes, variant_df){

	bets_somatic <- as.data.frame(read_delim(PATH_bets_somatic, delim = "\t"))
	bets_germline <- as.data.frame(read_delim(PATH_bets_germline, delim = "\t"))

	variant_df$In_germline_bets <- variant_df$Position %in% bets_germline$POSITION
	variant_df$In_somatic_bets <- variant_df$Position %in% bets_somatic$POSITION

	# add a new col to indicate if the gene of interest in each row was also found in the 2017 PCa panel design 	
	panel <- read_csv(PATH_panel_genes)
	variant_df$In_panel <- variant_df$Gene %in% panel$GENE

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
add_depth <- function(DIR_depth_metrics, PATH_collective_depth_metrics, variant_df) {

	file_name_list <- list.files(DIR_depth_metrics, full.names = TRUE)
	depth_file_list <- lapply(file_name_list, read.delim, header = FALSE) # load the related files to a list
	depth_file <- do.call(rbind, depth_file_list) # combine all depth files into one large one

	sample_names_repeated <- mapply(rep, unlist(lapply(file_name_list, basename)), unlist(lapply(depth_file_list, nrow)))
	sample_names_repeated <- unlist(lapply(sample_names_repeated, function(x) gsub(".txt", "", x)))
	depth_file$Sample_name <- sample_names_repeated # add the sample names so that we can do a join based on sample names with the variant df later on
	names(depth_file) <- c("Chrom", "Position", "Depth", "Sample_name")
	depth_file$Position <- as.character(depth_file$Position) # add the sample names so that we can do a join based on sample names with the variant df later on

	combined <- left_join(variant_df, depth_file, by = c("Sample_name", "Chrom", "Position"))

	# Add the median depth across all positions
	depth_file <- as.data.frame(read_delim(PATH_collective_depth_metrics, delim = "\t"))
	combined <- left_join(combined, depth_file, by = "Sample_name")

	return(combined)

}

subset_to_panel <- function(PATH_bed, variant_df) {
	
	variant_df$Chrom <- as.character(variant_df$Chrom)

	bed <- as.data.frame(read.delim(PATH_bed, header = FALSE, sep = "\t"))

	if (ncol(bed) == 4) {names(bed) <- c("chrom", "start", "stop", "gene")} else {names(bed) <- c("chrom", "start", "stop")}
	bed$chrom <- as.character(bed$chrom)

	to_keep <- list() # list of positions to remove

	message("Started subsetting to the panel.")
	i <- 1 
	while (i <= dim(variant_df)[1]){
		chrom_subsetted <- as.character(variant_df[i, 5]) # pick the chrom we are at 
		location <- as.character(variant_df[i, 6]) # pick the location we are at 
		bed_subsetted <- bed %>% dplyr::filter(chrom == chrom_subsetted) # subset the bed by chrom

		j <- 1
		while(j <= dim(bed_subsetted)[1]) {
			start_pos <- bed_subsetted[j, 2]
			end_pos <- bed_subsetted[j, 3]

			if (all(location >= start_pos, location <= end_pos)) {
				to_keep <- append(to_keep, i) # saving the row index of muts we are keeping
				break 
				} else {j <- j + 1} # if the location isn't in the panel, go check out the next position in the bed file.
				# print(c(j, "j"))

		} # end of inner while loop

	i <- i + 1 # next identified variant

	} # end of outer while loop - looping through identified variants
	
	to_print <- paste("Out of", nrow(variant_df), "variants", length(to_keep), "is retained as a part of the panel.")
	message(to_print)
	variant_df <- variant_df[unlist(to_keep), ]

	return(variant_df)

} # end of function

add_AAchange_effect <- function(variants_df){

	# Extract all annotations that start with "p." from each variant.
	p_list <- lapply(str_split(variants_df$AAChange.refGene, ":"), function(x) grep("p.", x, value = TRUE))
	p_list <- lapply(p_list, function(x) unique(gsub(",.*", "", x))) # gene name follows the annotation, remove it and subset to unique values
	
	add_effect <- function(pattern, effect_name, p_list){
		
		effects <- lapply(p_list, function(x) grepl(pattern, x)) # if it contains the string "fs" anywhere
		effects <- lapply(effects, function(my_vector) case_when(my_vector == TRUE ~ effect_name))

		return(effects)
	}

	# creating a list of first and second AA will help identify missense and synonymous mutations 
	first_AA <- lapply(p_list, function(x) substr(x, 3, 3)) # first AA
	second_AA <- lapply(p_list, function(x) substr(x, nchar(x), nchar(x))) # second AA that first one changes into
	second_AA <- lapply(second_AA, function(x) grep("s|X", x, invert = TRUE, value = TRUE)) # exclude the frameshift

	# missense mutations 
	idx <- mapply(function(x, y) {x != y}, first_AA, second_AA)
	missense_list <- lapply(idx, function(my_vector) case_when(my_vector == TRUE ~ "missense"))

	# synonymous mutations
	idx <- mapply(function(x, y) {x == y}, first_AA, second_AA)
	synonymous_list <- lapply(idx, function(my_vector) case_when(my_vector == TRUE ~ "synonymous"))

	# frameshift, stop gain and start loss
	fs_list <- add_effect("fs", "frameshift", p_list)
	stop_gain_list <- add_effect("X$", "stop_gain", p_list)
	start_loss_list <- add_effect("^p.M.+[a-zA-Z]$", "start_loss", p_list) # Methionine becomes something else

	# now combine all mutational effects into one list
	effects <- mapply(c, fs_list, stop_gain_list, start_loss_list, synonymous_list, missense_list, SIMPLIFY=TRUE)
	effects <- lapply(effects, function(x) x[!is.na(x)]) # remove all NAs

	# convert all elements with a length 0 to NA
	effects <- purrr::modify_if(effects, ~ length(.) == 0, ~ NA_character_)
	p_list <- purrr::modify_if(p_list, ~ length(.) == 0, ~ NA_character_)

	# collapse into single strings separated by :
	effects <- unlist(lapply(effects, function(x) paste(x, collapse = ":"))) # concatenate the strings from the same mutation with a colon 

	# add two new cols to the df
	variants_df$Protein_annotation <- unlist(lapply(p_list, function(x) paste(x, collapse = ":")))
	variants_df$Effects <- unlist(lapply(effects, function(x) paste(x, collapse = ":")))

	return(variants_df)

} # end of function

filter_tnvstats_by_variants <- function(variants_df, DIR_tnvstats, PATH_temp, PATH_filter_tnvstats_script, identifier){

	# Add 1 to indels so that the merge happens correctly
	variants_df <- variants_df %>% 
	mutate(
		Position = case_when(
			Variant == "deletion" ~ Position + 1, 
			TRUE ~ Position))

	write_delim(variants_df, PATH_temp, delim = "\t") # save a temp file before we add the read counts, depth and vaf from finland bams

	PATHS_tnvstats <- as.list(unique(grep(paste(variants_df$Sample_name_finland, collapse = "|"), list.files(DIR_tnvstats, full.names = TRUE, recursive = TRUE), value = TRUE))) # complete file paths to tnvstats
	PATHS_tnvstats <- unique(grep(".bam.tnvstat$", PATHS_tnvstats, value = TRUE)) # make sure we filter for the original tnvstat files

	# Run system commands for each tnvstats
	message("Filtering tnvstats in through system commands.")
	for (PATH_tnvstat in PATHS_tnvstats) {
		message(basename(PATH_tnvstat))
		system(paste("bash", PATH_filter_tnvstats_script, PATH_temp, PATH_tnvstat, identifier))
	}

	PATHS_filtered_tnvstats <- unique(
		grep(paste(variants_df$Sample_name_finland, collapse = "|"), 
			list.files(DIR_tnvstats, full.names = TRUE, recursive = TRUE), value = TRUE)) # complete file paths to tnvstats
	
	PATHS_filtered_tnvstats <- grep(identifier, PATHS_filtered_tnvstats, value = TRUE) # grep for the correct identifier to ensure we grab the right filtered file
	filtered_tnvstats <- lapply(PATHS_filtered_tnvstats, function(some_path) read_delim(some_path, delim = "\t")) # load all tnvstats
	filtered_tnvstats <- do.call(rbind, filtered_tnvstats) %>% select(-contains("_n"), -tumor_max_value, -gc, -target)

	# fix some formatting issues
	filtered_tnvstats$base_tumor <- gsub("TRUE", "T", filtered_tnvstats$base_tumor)
	filtered_tnvstats$ref <- gsub("TRUE", "T", filtered_tnvstats$ref)
	filtered_tnvstats$sample_t <- gsub(".bam", "", filtered_tnvstats$sample_t)

	return(filtered_tnvstats)

}

add_finland_readcounts <- function(variants_df, filtered_tnvstats, variant_caller){
	# Add the number of reads supporting the variants in Matti's bams. This helps justify that the variant is real, but we weren't able to detect it.
	
	modify_pos <- function(variants_df, variant_caller, add_or_subtract) {

		if (add_or_subtract == "add") {

				variants_df <- variants_df %>% mutate(
						Position = case_when(
							Variant == "deletion" ~ Position + 1, 
							TRUE ~ Position))

		} else {

				variants_df <- variants_df %>% mutate(
						Position = case_when(
							Variant == "deletion" ~ Position - 1, 
							TRUE ~ Position))

		} # end of if loop

		return(variants_df)

	} # end of function

	# Add 1 to the pos for indels because varcallers report one base too early
	variants_df <- modify_pos(variants_df, variant_caller, "add")
	
	variants_df_subsetted <- variants_df %>% 
			select(Chrom, Position, Ref, Alt, Sample_name_finland, Variant)

	# go through each variant and identify the correct position in the tnvstats file
	finland_df <- left_join(variants_df_subsetted, filtered_tnvstats, by = c("Chrom" = "chrom", "Position" = "pos", "Sample_name_finland" = "sample_t")) # this allows for exact matches
	
	# Choose the corresponding mutant reads based on what the alt is
	finland_df <- finland_df %>%
		mutate(
		F_tumor_alt_reads = case_when(
			Variant == "snv" & Alt == "A" ~ A_t, 
			Variant == "snv" & Alt == "C" ~ C_t, 
			Variant == "snv" & Alt == "T" ~ T_t, 
			Variant == "snv" & Alt == "G" ~ G_t,
			Variant == "deletion" ~ deletions_t, 
			Variant == "insertion" ~ insertions_t,
			TRUE ~ NA_real_), 
		F_tumor_vaf = case_when(
			Variant == "snv" & Alt == "A" ~ AAF_t, 
			Variant == "snv" & Alt == "C" ~ CAF_t, 
			Variant == "snv" & Alt == "T" ~ TAF_t, 
			Variant == "snv" & Alt == "G" ~ GAF_t, 
			Variant == "deletion" ~ deletions_t/reads_all_t, 
			Variant == "insertion" ~ insertions_t/reads_all_t,
			TRUE ~ NA_real_)) %>%
		select(Chrom, Position, Ref, Sample_name_finland, reads_all_t, F_tumor_alt_reads, F_tumor_vaf, N_t)

	# combine with the variants_df
	names(finland_df) <- c("Chrom", "Position", "Ref", "Sample_name_finland", "F_tumor_depth", "F_tumor_alt_reads", "F_tumor_vaf", "F_N_counts") # some renaming so that the merge happens with ease
	variants_df <- left_join(variants_df, finland_df)

	# revert back to the original position values 
	variants_df <- modify_pos(variants_df, variant_caller, "subtract")

	return(variants_df)
}

# Based on the annotation, identify which variants are duplicates without filtering them out.
identify_duplicates <- function(variants_df) {

	dups <- variants_df %>% 
			get_dupes(AAchange) %>%
			select(names(variants_df))

	dups$dupe_count <- NULL # drop an extra col that the function above adds
	non_dups <- setdiff(variants_df, dups)

	dups$Duplicate <- TRUE
	non_dups$Duplicate <- FALSE

	variants_df <- rbind(non_dups, dups)

	return(variants_df)
}

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

# Checks for any duplicated variants in the finalized output. Helpful to notice any errors.
check_duplicated_rows <- function(variants_df, drop_dups) {

	idx <- which(duplicated(variants_df))
	
	if (length(idx) > 0){
		warning_message <- paste("Duplicated rows detected at", idx)
		message(warning_message)

		if (drop_dups == TRUE){
			variants_df <- variants_df[!duplicated(variants_df), ]
		}
	} else {
		message("No duplicates. All good!")
	}

	return(variants_df)
}

