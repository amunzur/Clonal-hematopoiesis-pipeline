return_anno_output <- function(DIR_ANNOVAR) {

	df_main <- as.data.frame(read_delim(DIR_ANNOVAR, delim = "\t"))

	df <- df_main %>% 
		mutate(Sample_name = gsub(".hg38_multianno.txt", "", basename(DIR_ANNOVAR))) %>%
		select(Sample_name, Chr, Start, Ref, Alt, Func.refGene, Gene.refGene, Func.knownGene, Gene.knownGene, AAChange.refGene, ExAC_ALL, gnomAD_exome_ALL) 

	return(df)

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
				} else {j <- j + 1} # if the location isn't in the panel, go check out the next position in the bed file.
				# print(c(j, "j"))

		} # end of inner while loop

	i <- i + 1 # next identified variant

	} # end of outer while loop - looping through identified variants

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

add_finland_readcounts <- function(variants_df, DIR_tnvstats, PATH_temp, PATH_filter_tnvstats_script, identifier){
	# Add the number of reads supporting the variants in Matti's bams. This helps justify that the variant is real, but we weren't able to detect it.

	write_delim(dedup, PATH_temp, delim = "\t") # save a temp file before we add the read counts, depth and vaf from finland bams

	DIR_tnvstats <- "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/tnvstats/kidney_samples"
	PATH_temp <- "/groups/wyattgrp/users/amunzur/pipeline/results/temp/snv.tsv"
	PATH_filter_tnvstats_script <- "/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/analysis/filter_tnvstats.sh"
	PATHS_tnvstats <- as.list(unique(grep(paste(variants_df$Sample_name_finland, collapse = "|"), list.files(DIR_tnvstats, full.names = TRUE, recursive = TRUE), value = TRUE))) # complete file paths to tnvstats

	# Make a list of system commands to run
	LIST_system_command <- list()
	for (PATH_tnvstat in PATHS_tnvstats) {
		system_command <- paste("bash", PATH_filter_tnvstats_script, PATH_temp, PATH_tnvstat, identifier)
		LIST_system_command <- append(LIST_system_command, system_command)
	}

	# Run them in a loop to generate smaller tnvstats files that only contain the positions of interest
	for (system_command in LIST_system_command) {system(system_command)}

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





