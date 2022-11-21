return_anno_output <- function(PATH_annovar) {

	df_main <- as.data.frame(read_delim(PATH_annovar, delim = "\t"))

	df <- df_main %>% 
		mutate(Sample_name = gsub(".hg38_multianno.txt", "", basename(PATH_annovar))) %>%
		select(Sample_name, Chr, Start, Ref, Alt, Func.refGene, Gene.refGene, Func.knownGene, Gene.knownGene, AAChange.refGene, ExAC_ALL, gnomAD_exome_ALL) 

	return(df)

}

# add a new column to indicate sample type, WBC or tumor. 
add_sample_type <- function(variant_df){

	x <- str_split(variant_df$Sample_name, "_") # splitted strings as a list
	y <- unlist(lapply(x, "[", 2)) # gDNA or cfDNA

	my_map <- c("gDNA" = "WBC", "gdna" = "WBC", "cfDNA" = "Tumor")
	variant_df$Sample_type <- unname(my_map[y])

	return(variant_df)

}

# add the depth at each variant position
add_depth <- function(DIR_depth_metrics, PATH_collective_depth_metrics, variant_df) {

	file_name_list <- list.files(DIR_depth_metrics, full.names = TRUE, pattern = "\\.txt$")
	depth_file_list <- lapply(file_name_list, read.delim, header = FALSE) # load the related files to a list
	depth_file <- do.call(rbind, depth_file_list) # combine all depth files into one large one

	sample_names_repeated <- mapply(rep, unlist(lapply(file_name_list, basename)), unlist(lapply(depth_file_list, nrow)))
	sample_names_repeated <- unlist(lapply(sample_names_repeated, function(x) gsub(".txt", "", x)))
	depth_file$Sample_name <- sample_names_repeated # add the sample names so that we can do a join based on sample names with the variant df later on
	names(depth_file) <- c("Chrom", "Position", "Depth", "Sample_name")
	depth_file$Position <- as.character(depth_file$Position) # add the sample names so that we can do a join based on sample names with the variant df later on
	depth_file$Chrom <- as.character(depth_file$Chrom)
	combined <- left_join(variant_df, depth_file, by = c("Sample_name", "Chrom", "Position"))

	# Add the median depth across all positions
	depth_file <- as.data.frame(read_delim(PATH_collective_depth_metrics, delim = " "))
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
		chrom_subsetted <- as.character(variant_df[i, "Chrom"]) # pick the chrom we are at 
		position <- as.character(variant_df[i, "Position"]) # pick the genomic we are at 
		bed_subsetted <- bed %>% dplyr::filter(chrom == chrom_subsetted) # subset the bed by chrom

		j <- 1
		while(j <= dim(bed_subsetted)[1]) {
			start_pos <- bed_subsetted[j, 2]
			end_pos <- bed_subsetted[j, 3]

			if (all(position >= start_pos, position <= end_pos)) {
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

	add_effect <- function(pattern, effect_name, p_list){
		
		effects <- lapply(p_list, function(x) grepl(pattern, x)) # if it contains the string "fs" anywhere
		effects <- lapply(effects, function(my_vector) case_when(my_vector == TRUE ~ effect_name))

		return(effects)
	}

	# Extract all annotations that start with "p." from each variant.
	p_list <- lapply(str_split(variants_df$AAChange.refGene, ",|:"), function(x) grep("^p.", x, value = TRUE))
	p_list <- lapply(p_list, function(x) unique(gsub(",.*", "", x))) # gene name follows the annotation, remove it and subset to unique values
	
	# creating a list of first and second AA will help identify missense and synonymous mutations 
	first_AA <- lapply(p_list, function(x) substr(x, 3, 3)) # first AA
	first_AA[grep("delins|del$", p_list)] <- NA # exclude nonfs deletions and complex delins
	second_AA <- lapply(p_list, function(x) substr(x, nchar(x), nchar(x))) # second AA that first one changes into
	second_AA <- lapply(second_AA, function(x) grep("s|X", x, invert = TRUE, value = TRUE)) # exclude the frameshift

	# missense mutations 
	idx <- mapply(function(x, y) {x != y}, first_AA, second_AA)
	missense_list <- lapply(idx, function(my_vector) case_when(my_vector == TRUE ~ "missense"))

	# synonymous mutations
	idx <- mapply(function(x, y) {x == y}, first_AA, second_AA)
	synonymous_list <- lapply(idx, function(my_vector) case_when(my_vector == TRUE ~ "synonymous"))

	# complicated deletions and insertions together
	idx <- grepl("delins", p_list)
	delins_list <- lapply(idx, function(my_vector) case_when(my_vector == TRUE ~ "non_frameshift_delins"))

	# non frameshift deletion
	idx <- grepl("del$", p_list)
	nonfsdel_list <- lapply(idx, function(my_vector) case_when(my_vector == TRUE ~ "nonframeshift_deletion"))

	# non frameshift insertion
	idx <- grepl("ins$", p_list)
	nonfsins_list <- lapply(idx, function(my_vector) case_when(my_vector == TRUE ~ "nonframeshift_insertion"))

	# frameshift, stop gain and start loss
	fs_list <- add_effect("fs", "frameshift", p_list)
	stop_gain_list <- add_effect("X$", "stop_gain", p_list)
	start_loss_list <- add_effect("^p.M1.+[a-zA-Z]$", "start_loss", p_list) # Methionine1 becomes something else

	# now combine all mutational effects into one list
	effects <- mapply(c, fs_list, stop_gain_list, start_loss_list, synonymous_list, missense_list, delins_list, nonfsdel_list, nonfsins_list, SIMPLIFY=TRUE)
	effects <- lapply(effects, function(x) x[!is.na(x)]) # remove all NAs

	# convert all elements with a length 0 to NA
	effects <- purrr::modify_if(effects, ~ length(.) == 0, ~ NA_character_)
	p_list <- purrr::modify_if(p_list, ~ length(.) == 0, ~ NA_character_)

	# Choosing the longest spicing variant from the protein_annotation column
	idx_list <- lapply(p_list, function(some_vector) which.max(str_extract(some_vector, "[[:digit:]]+")))
	idx_list[lengths(idx_list) == 0] <- NA # if annotation is NA, set the index to 1
	idx_list <- unlist(idx_list)

	# add two new cols to the df
	variants_df$Protein_annotation <- unlist(mapply("[", p_list, idx_list))
	variants_df$Effects <- unlist(mapply("[", effects, idx_list))
	variants_df$Effects[grep("splicing", variants_df$Func.refGene)] <- "splicing"

	return(variants_df)

} # end of function

# based on the cohort_name, add the patient id
add_patient_id <- function(variant_df){

	x <- str_split(variant_df$Sample_name, "_gDNA|_WBC|_cfDNA")
	variant_df$patient_id <- unlist(lapply(x, "[", 1))

	return(variant_df)

}

# Removes variants if they are found in more than n_times (appears multiple times)
# For the remaining, add a new column to indicate if the variant occurs more than once
find_and_filter_duplicated_variants <- function(variants_df, n_times) {

	# Removing
	tab <- table(variants_df$AAchange)
	variants_df <- variants_df[variants_df$AAchange %in% names(tab[tab < n_times]), ]

	# Adding a new column to indicate duplicates
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

process_basecounts_vcf <- function(PATH_basecounts) {
	# PATH_basecounts <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/base_counts/DCS/GUBB-18-029_gDNA_Baseline_IDT_2018Apr18.vcf"

	vcf <- read.table(PATH_basecounts, stringsAsFactors = FALSE, col.names = c("Chrom", "Position", "ID", "Ref", "Alt", "Qual", "Filter", "Info", "Format", "Values")) %>%
			separate(col = Values, sep = ":", into = c("DP", "Ref_reads", "Alt_reads", "VAF", "DPP", "DPN", "RDP", "RDN", "ADP", "ADN")) %>%
			mutate(Sample_name = gsub(".vcf", "", basename(PATH_basecounts))) %>%
			select(Sample_name, Chrom, Position, Ref, Alt, Ref_reads, Alt_reads, VAF) %>%
			mutate(Type = case_when(
								nchar(Ref) > nchar(Alt) ~ "Deletion", 
								nchar(Ref) < nchar(Alt) ~ "Insertion",
								nchar(Ref) == nchar(Alt) ~ "SNV", 
								TRUE ~ "Error")) %>%
			mutate(Position = as.double(Position)) %>%
			mutate(Position = case_when(
								Type == "Deletion" ~ as.double(Position) + 1, 
								TRUE ~ Position)) %>%
			mutate(Ref = ifelse(Type == "Deletion", sub(".", "", Ref), Ref), 
				   Alt = ifelse(Type == "Deletion", sub(".", "-", Alt), Alt), 
				   Ref = ifelse(Type == "Insertion", sub(".", "-", Ref), Ref), 
				   Alt = ifelse(Type == "Insertion", sub(".", "", Alt), Alt)) %>%
			mutate(Position = as.character(Position))
	return(vcf)

}

parse_basecount_vcf <- function(DIR_basecounts){

	vcf_list <- lapply(as.list(list.files(DIR_basecounts, full.names = TRUE, pattern = "gDNA.+vcf$")), process_basecounts_vcf) # only look for cfDNA files
	vcf_df <- as.data.frame(do.call(rbind, vcf_list)) %>% mutate(Position = as.character(Position))
}

add_cfDNA_information <- function(variants_df, DIR_basecounts_DCS, DIR_basecounts_SSCS1) {

	# combine with the DCS - only keep the vars with read support in the cfDNA DCS file
	basecounts_df_DCS <- parse_basecount_vcf(DIR_basecounts_DCS)
	variants_df <- left_join(variants_df, basecounts_df_DCS, by = c("Chrom", "Position", "Ref", "Alt", "Sample_name"))
	names(variants_df) <- gsub("\\.x", "_n", names(variants_df))
    names(variants_df) <- gsub("\\.y", "_DCS_t", names(variants_df))
	variants_df <- filter(variants_df, Alt_reads_DCS_t > 0)

	# combine with SSCS1 - to obtain the cfDNA VAF
	basecounts_df_SSCS1 <- parse_basecount_vcf(DIR_basecounts_SSCS1)
	names(basecounts_df_SSCS1)[c(6, 7, 8)] <- c("Ref_reads_SSCS1_t", "Alt_reads_SSCS1_t", "VAF_SSCS1_t")
	variants_df <- left_join(variants_df, basecounts_df_SSCS1)

	return(variants_df)
}

combine_tumor_wbc <- function(variants_df){

    variants_df <- variants_df %>% mutate(Sample_type = str_extract(Sample_name, "cfDNA|gDNA"))
	tumor <- variants_df %>% filter(Sample_type == "cfDNA") %>% select(-Sample_type)
	wbc <- variants_df %>% filter(Sample_type == "gDNA") 

    combined <- inner_join(tumor, wbc, by = c("Patient_ID", "Chrom", "Position", "Ref", "Alt", "Function", "Gene", "AAchange", "Protein_annotation", "Effects", "ExAC_ALL", "Variant"))
    names(combined) <- gsub("\\.x", "_t", names(combined)) 
    names(combined) <- gsub("\\.y", "_n", names(combined)) 

    combined <- combined %>%
            mutate(tumor_wbc_vaf_ratio = round((VAF_t / VAF_n), 2), 
                    tumor_wbc_depth_ratio = round((Depth_t / Depth_n), 2))

    return(combined)
}

blacklist_variants <- function(variants_df, PATH_blacklist) {
	df_blacklist <- as.data.frame(read_csv(PATH_blacklist)) %>% mutate(Position = as.character(Position))
	combined <- anti_join(variants_df, df_blacklist)
}

filter_pileup <- function(PATH_mpileup, DIR_mpileup_filtered, variants_df) {

	df <- variants_df %>% 
		  filter(Sample_name == gsub(".mpileup", "", basename(PATH_mpileup))) %>% 
		  select(Chrom, Position) %>% 
		  distinct(.keep_all = TRUE) %>%
		  unite(combined_names, Chrom, Position, sep = "[[:blank:]]")

	message("Grepping ", gsub(".mpileup", "", basename(PATH_mpileup)))
	if (!file.exists(file.path(DIR_mpileup_filtered, basename(PATH_mpileup)))) {
		system(paste0("grep -E ", "'", (paste(df$combined_names, collapse = "|")), "'", " ", PATH_mpileup, " > ", file.path(DIR_mpileup_filtered, basename(PATH_mpileup))))
	} else {
		message("Filtered pileup exists, skipping this sample.")
	}
	
	mpileup <- read_delim(file.path(DIR_mpileup_filtered, basename(PATH_mpileup)), delim = "\t", col_names = c("Chrom", "Position", "Ref", "Depth", "Read_bases", "Read_quality")) %>% select(-Read_quality)

	mpileup <- mpileup %>% 
			   mutate(N_bases = str_count(Read_bases, "N|n"),
			   		  N_fraction = as.numeric(N_bases) / as.numeric(Depth), 
					  Sample_name = gsub(".mpileup", "", basename(PATH_mpileup))) %>%

	return(mpileup)
}

add_N_fraction <- function(variants_df, DIR_mpileup, DIR_mpileup_filtered) {

	idx <- grep(paste(unique(variants_df$Sample_name), collapse = "|"), list.files(DIR_mpileup, full.names = TRUE))
	files <- list.files(DIR_mpileup, full.names = TRUE)[idx]

	mpileup_list <- lapply(files, filter_pileup, DIR_mpileup_filtered = DIR_mpileup_filtered, variants_df = variants_df)
	mpileup <- do.call(rbind, mpileup_list) %>% 
				select(Sample_name, Chrom, Position, N_fraction)

	combined <- left_join(variants_df, mpileup) 

}