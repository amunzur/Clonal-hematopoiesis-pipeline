# Some callers will reports REF and ALT with two bases like TG T, meaning G was deleted. 
# However this causes in some other downstream tools, so here I convert all to a single bp notation.
# So TG T would become G -.
# Assuming there is a column named TYPE that has the type information like SNV, Deletion and Insertion.
standardize_ref_and_alt <- function(vars) {
	vars <- vars %>%
		mutate(REF = 
				case_when(
					Variant == "Deletion" ~ sub('.', '', Ref),
					Variant == "Insertion" ~ "-",
					TRUE ~ Ref
						 ), 
			   ALT = 
				case_when(
					Variant == "Deletion" ~ "-",
					Variant == "Insertion" ~ sub('.', '', Alt),
					TRUE ~ Alt
						 ))
	return(vars)
}

add_bg_error_rate <- function(variants_df, bg) {

	# modify the vars df
	variants_df <- variants_df %>% mutate(error_type = paste0("mean_error", variants_df$Ref, "to", variants_df$Alt))
	variants_df <- variants_df %>% mutate(Position = ifelse(Type == "Deletion", Position + 1, Position))

	# add the deletions 
	idx <- which(tolower(variants_df$Type) == "deletion")
	variants_df$error_type[idx] <- "mean_errordel"
	
	# add the insertions
	idx <- which(tolower(variants_df$Type) == "insertion")
	variants_df$error_type[idx] <- "mean_errorins"

	# modify the bg error rate df
	bg <- gather(bg, "error_type", "error_rate", starts_with("mean_error"))
	bg$POS <- as.numeric(bg$POS)
	bg <- select(bg, -REF)
	
	# this merge adds an exta col with the error rate 
	variants_df <- left_join(variants_df, bg, by = c("Chrom" = "CHROM", "Position" = "POS", "error_type" = "error_type")) # to add the error rates from bg to variants df
	
	if ("VAF" %in% colnames(variants_df)) {
		variants_df <- mutate(variants_df, VAF_bg_ratio=VAF/error_rate)
	} else if ("VAF_t" %in% colnames(variants_df)) {
		variants_df <- mutate(variants_df, VAF_bg_ratio=VAF_t/error_rate)
	}
	
	return(variants_df)
}

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

filter_somatic_variants <- function(vars, min_alt_reads_t, max_alt_reads_n, min_depth_n, min_VAF_low, min_tumor_to_normal_vaf_ratio, min_VAF_bg_ratio, PATH_blacklist, blacklist = TRUE) {
	vars$VAF_bg_ratio <- na.fill(vars$VAF_bg_ratio, fill = 99999)
	vars <- vars %>%
		mutate(Status = "SOMATIC") %>%
		filter(Alt_forward_t + Alt_reverse_t >= min_alt_reads_t, 
			   VAF_t >= min_VAF_low, 
			   Depth_n >= min_depth_n, 
			   (tumor_to_normal_VAF_ratio >= min_tumor_to_normal_vaf_ratio | Alt_forward_n + Alt_reverse_n <= max_alt_reads_n), 
			   VAF_bg_ratio > min_VAF_bg_ratio, 
			   Strand_bias_fishers != TRUE, 
			   ifelse(Gene == "TERT", Function %in% c("exonic", "splicing", "upstream"), Function %in% c("exonic", "splicing")),
			   is.na(Effects) | Effects != "synonymous")	

	if (blacklist == TRUE) {
		blacklist_df_protein_annot <- as.data.frame(read_csv(PATH_blacklist)) %>% select(Chrom, Gene, Protein_annotation) %>% filter(!is.na(Protein_annotation))
		blacklist_df_position <- as.data.frame(read_csv(PATH_blacklist)) %>% select(Chrom, Gene, Position) %>% filter(!is.na(Position))

		vars <- anti_join(vars, blacklist_df_protein_annot)
		vars <- anti_join(vars, blacklist_df_position)
	}
	return(vars)
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

# Based on the number of forward and reverse alt and ref reads, calculates strand bias using odds ratio and fisher's exact test
evaluate_strand_bias <- function(variants_df) {
	if (all(c("Alt_forward_t", "Alt_reverse_t", "Ref_forward_t", "Ref_reverse_t") %in% colnames(variants_df))) {
		suffix <- "_t"
	} else (
		suffix <- ""
	)
	
	f <- function(row) {
    # Create a 2x2 contingency table for the current row
    contingency_table <- matrix(c(row[paste0("Alt_forward", suffix)], row[paste0("Alt_reverse", suffix)],
                                  row[paste0("Ref_forward", suffix)], row[paste0("Ref_reverse", suffix)]),
                                nrow = 2,
                                dimnames = list(c("Alt", "Ref"),
                                                c("Forward", "Reverse")))
	
	if (any(contingency_table == 0)) {
      p.value <- 0.05 # Set a default p-value indicating potential strand bias
    } else {
      p.value <- fisher.test(contingency_table)$p.value
    }
    return(p.value<=0.05) # if TRUE, this means there is strand bias.
  }

	variants_df$Strand_bias_fishers <- apply(variants_df[, c(paste0("Alt_forward", suffix), paste0("Ref_forward", suffix), paste0("Alt_reverse", suffix), paste0("Ref_reverse", suffix))], 1, f)
	return(variants_df)
}

evaluate_strand_bias2 <- function(vars) {
	if (all(c("Alt_forward_t", "Alt_reverse_t", "Ref_forward_t", "Ref_reverse_t") %in% colnames(vars))) {
		suffix <- "_t"
	} else (
		suffix <- ""
	)

	vars$Alt_Reads_total <- vars[paste0("Alt_forward", suffix)] + vars[paste0("Alt_reverse", suffix)]
	condition1 <- vars[paste0("Alt_forward", suffix)] >= 0.95*vars["Alt_Reads_total"]
	condition2 <- vars[paste0("Alt_forward", suffix)] <= 0.05*vars["Alt_Reads_total"]
	vars$Strand_bias_fishers <- c(condition1 | condition2)
	vars <- select(vars, -Alt_Reads_total)

	return(vars)
}

add_AAchange_effect <- function(vars, variant_caller){

	add_effect <- function(pattern, effect_name, p_list){
		
		effects <- lapply(p_list, function(x) grepl(pattern, x)) # if it contains the string "fs" anywhere
		effects <- lapply(effects, function(my_vector) case_when(my_vector == TRUE ~ effect_name))

		return(effects)
	}

	# Extract all annotations that start with "p." from each variant.
	p_list <- lapply(str_split(vars$AAchange, ",|:"), function(x) grep("^p.", x, value = TRUE))
	p_list <- lapply(p_list, function(x) unique(gsub(",.*", "", x))) # gene name follows the annotation, remove it and subset to unique values
	
	# creating a list of first and second AA will help identify missense and synonymous mutations 
	first_AA <- lapply(p_list, function(x) substr(x, 3, 3)) # first AA
	first_AA[grep("delins|del|ins", p_list)] <- NA # exclude nonfs deletions and complex delins
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
	vars$Protein_annotation <- unlist(mapply("[", p_list, idx_list))
	vars$Effects <- unlist(mapply("[", effects, idx_list))
	vars$Effects[grep("splicing", vars$Function)] <- "splicing"
	
	# move the newly added columsn to the beginning of the dataframe
	# alt_col_index <- which(names(vars) == "Alt")
	# vars <- vars %>% 
	# 		select(Patient_id, Sample, Chrom, Position, Ref, Alt, Protein_annotation, Effects, everything()[-(1:(alt_col_index+2))]) %>%
	# if (variant_caller == "Mutect2") {
	# 	vars <- select(vars, Patient_id, Sample_name, Sample_type, Chrom, Position, Ref, Alt, Type, VAF, Function, Gene, Consequence, AAchange, Protein_annotation, Effects, Depth, Ref_forward, Ref_reverse, Alt_forward, Alt_reverse, cosmic97_coding, avsnp150, CLNALLELEID, CLNSIG, error_type, error_rate)
	# } else if (variant_caller == "Freebayes") {
	# 	vars <- select(vars, Patient_id, Sample_name, Sample_type, Chrom, Position, Ref, Alt, Type, VAF, Function, Gene, Consequence, AAchange, Protein_annotation, Effects, Ref_forward, Ref_reverse, Alt_forward, Alt_reverse, cosmic97_coding, avsnp150, CLNALLELEID, CLNSIG, error_type, error_rate)
	# } else if (variant_caller == "Vardict") {
	# 	vars <- select(vars, Patient_id, Sample_name, Sample_type, Chrom, Position, Ref, Alt, Type, VAF, Function, Gene, Consequence, AAchange, Protein_annotation, Effects, Depth, Ref_forward, Ref_reverse, Alt_forward, Alt_reverse, cosmic97_coding, avsnp150, CLNALLELEID, CLNSIG, error_type, error_rate)	
	# }
	
	return(vars)

} # end of function

# based on the cohort_name, add the patient id
add_patient_information <- function(variant_df, PATH_sample_information){

	x <- str_split(variant_df$Sample_name, "_gDNA|_WBC|_cfDNA")
	variant_df$Patient_id <- unlist(lapply(x, "[", 1))
	variant_df$Patient_id <- gsub("GU-|GUBB-", "", variant_df$Patient_id)

	sample_info <- read_delim(PATH_sample_information, delim = "\t", col_names = FALSE)
	names(sample_info) <- c("Patient_id", "Date_collected", "Diagnosis", "Timepoint")
	variant_df <- left_join(variant_df, sample_info)

	variant_df <- variant_df[, c("Patient_id", "Sample_name", "Sample_type", "Date_collected", "Timepoint",  "Diagnosis", colnames(variant_df)[-which(names(variant_df) %in% c("Patient_id", "Sample_name", "Sample_type", "Date_collected", "Timepoint",  "Diagnosis"))])] # move to the beginning
	return(variant_df)

}

add_patient_information_somatic <- function(variant_df, PATH_sample_information){

	sample_info <- read_delim(PATH_sample_information, delim = "\t", col_names = FALSE)
	names(sample_info) <- c("Patient_id", "Date_collected", "Diagnosis", "Timepoint")
	variant_df <- left_join(variant_df, sample_info)

	variant_df <- variant_df[, c("Patient_id", "Sample_name_t", "Date_collected", "Timepoint",  "Diagnosis", colnames(variant_df)[-which(names(variant_df) %in% c("Patient_id", "Sample_name_t", "Date_collected", "Timepoint",  "Diagnosis"))])] # move to the beginning
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

combine_tumor_wbc <- function(vars, match_cfDNA_and_WBC_dates = TRUE){
	# Same collection date set to TRUE will require matching entries for Timepoint and Date_collected columns.

	tumor <- vars %>% filter(Sample_type == "cfDNA") %>% select(-Sample_type)
	wbc <- vars %>% filter(Sample_type == "WBC") %>% select(-Sample_type)

	if (match_cfDNA_and_WBC_dates == TRUE ) {
    	combined <- inner_join(tumor, wbc, by = c("Patient_id", "Date_collected", "Timepoint", "Diagnosis", "Chrom", "Position", "Ref", "Alt", "Type", "Function", "Gene", "Consequence", "AAchange", "Protein_annotation", "Effects", "cosmic97_coding", "avsnp150", "CLNALLELEID", "CLNSIG", "Variant_caller"))
	} else {
    	combined <- inner_join(tumor, wbc, by = c("Patient_id", "Diagnosis", "Chrom", "Position", "Ref", "Alt", "Type", "Function", "Gene", "Consequence", "AAchange", "Protein_annotation", "Effects", "cosmic97_coding", "avsnp150", "CLNALLELEID", "CLNSIG", "Variant_caller"))
	}
    names(combined) <- gsub("\\.x", "_t", names(combined)) 
    names(combined) <- gsub("\\.y", "_n", names(combined)) 

    combined <- combined %>%
            mutate(tumor_wbc_vaf_ratio = round((VAF_t / VAF_n), 2), 
                    tumor_wbc_depth_ratio = round((Depth_t / Depth_n), 2))
	
	combined <- combined[, c(setdiff(names(combined), "Variant_caller"), "Variant_caller")] # move the variant caller column to the end of the df

    return(combined)
}

blacklist_variants <- function(vars, PATH_blacklist) {
	df_blacklist <- as.data.frame(read_csv(PATH_blacklist)) %>% mutate(Position = as.character(Position))
	combined <- anti_join(vars, df_blacklist)
}

filter_pileup <- function(PATH_mpileup, DIR_mpileup_filtered, vars, force) {

	PATH_mpileup_filtered <- file.path(DIR_mpileup_filtered, basename(PATH_mpileup))
	if (!file.exists(DIR_mpileup_filtered)) {dir.create(DIR_mpileup_filtered, recursive = TRUE, showWarnings = FALSE)} 

	if ("Sample_name_t" %in% colnames(vars)) {
		df <- vars %>% filter(Sample_name_t == file_path_sans_ext(basename(PATH_mpileup)))
	} else {
		df <- vars %>% filter(Sample_name == file_path_sans_ext(basename(PATH_mpileup)))
	}
	df <- df %>% 
		  select(Chrom, Position) %>% 
		  distinct(.keep_all = TRUE) %>%
		  unite(combined_names, Chrom, Position, sep = "[[:space:]]")

	# Check if 'df' has more than 5000 rows
	if (!file.exists(PATH_mpileup_filtered) | force) {
		if (nrow(df) > 5000) {
  			split_dfs <- split(df, ceiling(seq_along(df$combined_names) / 5000))

  			# Run the system command on each smaller dataframe
  			for (i in seq_along(split_dfs)) {
    			current_df <- split_dfs[[i]]
    			current_file <- file.path(DIR_mpileup_filtered, paste0(basename(PATH_mpileup), "_part", i))
    			current_grep_command <- paste0("grep -E '", paste(current_df$combined_names, collapse = "|"), "' ", PATH_mpileup, " > ", current_file)

    			# Execute the system command
    			message("Grepping ", gsub(".mpileup", "", basename(PATH_mpileup)), " - Part ", i)
    			system(current_grep_command)
  			}

  			# Combine the existing files into one file
  			combined_file <- file.path(DIR_mpileup_filtered, basename(PATH_mpileup))
  			combined_command <- paste0("cat ", file.path(DIR_mpileup_filtered, basename(PATH_mpileup)), "_part* > ", combined_file)

  			# Execute the system command to combine files
  			system(combined_command)
			} else {
  			# Execute the system command directly without splitting
  			message("Grepping ", gsub(".mpileup", "", basename(PATH_mpileup)))
  			system(paste0("grep -E '", paste(df$combined_names, collapse = "|"), "' ", PATH_mpileup, " > ", file.path(DIR_mpileup_filtered, basename(PATH_mpileup))))
			}
	}
	print(file.path(DIR_mpileup_filtered, basename(PATH_mpileup)))
	mpileup <- read_delim(file.path(DIR_mpileup_filtered, basename(PATH_mpileup)), delim = "\t", col_names = c("Chrom", "Position", "Ref", "Depth", "Read_bases", "Read_quality")) %>% select(-Read_quality)
	
	if ("Sample_name_t" %in% colnames(vars)) {
		mpileup <- mpileup %>% mutate(Sample_name_t = file_path_sans_ext(basename(PATH_mpileup)))
	} else {
		mpileup <- mpileup %>% mutate(Sample_name = file_path_sans_ext(basename(PATH_mpileup)))
	}

	mpileup <- mpileup %>% 
			   mutate(N_bases = str_count(Read_bases, "N|n"),
			   		  N_fraction = as.numeric(N_bases) / as.numeric(Depth))
	return(mpileup)
}

add_N_fraction <- function(vars, DIR_mpileup, DIR_mpileup_filtered, force) {

	if ("Sample_name_t" %in% colnames(vars)) {
		idx <- grep(paste(unique(vars$Sample_name_t), collapse = "|"), list.files(DIR_mpileup, full.names = TRUE))
	} else {
		idx <- grep(paste(unique(vars$Sample_name), collapse = "|"), list.files(DIR_mpileup, full.names = TRUE))
	}
	
	files <- list.files(DIR_mpileup, full.names = TRUE)[idx]
	mpileup_list <- lapply(files, filter_pileup, DIR_mpileup_filtered = DIR_mpileup_filtered, vars = vars, force = force)

	if ("Sample_name_t" %in% colnames(vars)) {
		mpileup <- do.call(rbind, mpileup_list) %>% select(Sample_name_t, Chrom, Position, N_fraction, N_bases)
	} else {
		mpileup <- do.call(rbind, mpileup_list) %>% select(Sample_name, Chrom, Position, N_fraction, N_bases)
	}	
	combined <- left_join(vars, mpileup)

	return(combined)

}

filter_variants_chip_or_germline <- function(mode, vars, min_alt_reads, min_depth, min_VAF_low, max_VAF_low, min_VAF_high, max_VAF_high, min_VAF_bg_ratio, PATH_blacklist, blacklist = TRUE, filter_by_min_depth = TRUE) {
	vars$VAF_bg_ratio <- na.fill(vars$VAF_bg_ratio, fill = 9999999)
	# When gnomad is undefined replace it with 0.
	vars$gnomad40_exome_AF <- as.numeric(vars$gnomad40_exome_AF)
	mask <- is.na(vars$gnomad40_exome_AF)
	vars$gnomad40_exome_AF[mask] <- 0
	
	vars <- vars %>% 
			mutate(Status = toupper(mode)) %>% 
			filter(Alt_forward + Alt_reverse >= min_alt_reads,
			       VAF_bg_ratio >= min_VAF_bg_ratio, 
				   gnomad40_exome_AF < 0.0005)
	
	if (filter_by_min_depth) {
		vars <- vars %>% 
			filter(Alt_forward + Ref_forward + Alt_reverse + Ref_reverse >= min_depth)
	} 
	
	if (toupper(mode) == "CHIP") {
		vars <- vars %>%
			filter((VAF >= min_VAF_low & VAF <= max_VAF_low) | (VAF >= min_VAF_high & VAF <= max_VAF_high), 
					ifelse(Gene == "TERT", Function %in% c("exonic", "splicing", "upstream"), Function %in% c("exonic", "splicing")),
			   		is.na(Effects) | Effects != "synonymous")	
	} else if (toupper(mode) == "GERMLINE"){
		vars <- vars %>%
			filter((VAF > min_VAF_low & VAF < max_VAF_low) | (VAF > min_VAF_high & VAF <= max_VAF_high), 
					Sample_type == "WBC")
	}

	if (blacklist == TRUE) {

		blacklist_df_protein_annot <- as.data.frame(read_csv(PATH_blacklist)) %>% select(Chrom, Gene, Protein_annotation) %>% filter(!is.na(Protein_annotation))
		blacklist_df_position <- as.data.frame(read_csv(PATH_blacklist)) %>% select(Chrom, Gene, Position) %>% filter(!is.na(Position))

		vars <- anti_join(vars, blacklist_df_protein_annot)
		vars <- anti_join(vars, blacklist_df_position)

		# result_df <- anti_join(df, exclude_ranges, by = 'chromosome') %>%
  		# 			 filter(end < start_exclude | start > end_exclude)
	}
	return(vars)
}

# FUNCTIONS TO LOAD VARIANT TABLES FROM ANNOVAR
# do a merge based on column to combine metadata 
parse_anno_output <- function(DIR_variant_tables, mode, variant_caller, PATH_sample_list = NULL) {

	if (toupper(variant_caller) == "VARDICT") {
		somatic_func = get("return_anno_output_vardict_somatic")
		chip_func = get("return_anno_output_vardict_chip")
	} else if (toupper(variant_caller) == "MUTECT2") {
		somatic_func = get("return_anno_output_mutect_somatic")
		chip_func = get("return_anno_output_mutect_chip")
	} else if (toupper(variant_caller) == "FREEBAYES") {
		# somatic_func = get("return_anno_output_freebayes_somatic")
		chip_func = get("return_anno_output_freebayes_chip")
	} else {
		stop("Invalid variant caller specified.")}

	if (is.null(PATH_sample_list)) {
		if (mode == "somatic") {
			anno_df_list <- lapply(as.list(list.files(DIR_variant_tables, full.names = TRUE, pattern = ".tsv$")), somatic_func)
		} else if (mode == "chip") {
			anno_df_list <- lapply(as.list(list.files(DIR_variant_tables, full.names = TRUE, pattern = ".tsv$")), chip_func)
		} 
	} else {
		# We subset to a certain group of samples.
		sample_df <- as.data.frame(read_delim(PATH_sample_list, delim = "\t"))
		samples <- c(sample_df$cfDNA, sample_df$WBC)

		# List files in the directory matching the samples
		files_to_load <- list.files(DIR_variant_tables, full.names = TRUE, pattern = ".tsv$") # all files in the dir
		files_to_load <- files_to_load[sapply(files_to_load, function(file) any(sapply(samples, grepl, file)))] # choose a subset based on the provided samples file
		
		if (mode == "somatic") {
			anno_df_list <- lapply(as.list(files_to_load), somatic_func)
		} else if (mode == "chip") {
			anno_df_list <- lapply(as.list(files_to_load), chip_func)
		}
	}
	anno_df <- as.data.frame(do.call(rbind, anno_df_list))
	return(anno_df)
}

return_anno_output_mutect_chip <- function(PATH_variant_table) {

	print(PATH_variant_table)
	df <- as.data.frame(read.delim(PATH_variant_table, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE))
	colnames(df) <- gsub(paste0(file_path_sans_ext(basename(PATH_variant_table)), "."), "", colnames(df)) 

	df <- df %>%
		separate(col = SB, sep = ",", into = c("Ref_forward", "Ref_reverse", "Alt_forward", "Alt_reverse"), remove = TRUE) %>%
		select(-EVENTLENGTH) %>% # really didnt need to include this column in the variant tables
		mutate(Sample = gsub(".tsv", "", basename(PATH_variant_table)), 
			   Sample_type = str_extract(basename(PATH_variant_table), "cfDNA|WBC"),
			   Sample_name = file_path_sans_ext(basename(PATH_variant_table)),
			   Date_collected = str_split(Sample_name, "[-_]")  %>% sapply(tail, 1), 
			   nchar_ref = nchar(REF), 
			   nchar_alt = nchar(ALT), 
			   ALT = ifelse(grepl(",", ALT), sub(",.*", "", ALT), ALT)) %>% # Select the first element if ALT contains commas
		mutate(across(c(Alt_forward, Alt_reverse, Ref_forward, Ref_reverse), as.numeric)) %>%
		mutate(TYPE = case_when(nchar_ref > nchar_alt ~ "Deletion", 
			   					   nchar_ref < nchar_alt ~ "Insertion", 
								   nchar_ref == nchar_alt ~ "SNV"), 
			   VAF = (Alt_forward+Alt_reverse)/(Ref_forward+Ref_reverse+Alt_forward+Alt_reverse), 
			   Depth = Ref_forward+Ref_reverse+Alt_forward+Alt_reverse) %>%
		select(-nchar_ref, -nchar_alt) %>%
		select(Sample_name, Sample_type, Date_collected, CHROM, POS, REF, ALT, TYPE, VAF, Depth, Alt_forward, Ref_forward, Alt_reverse, Ref_reverse, Func.refGene, Gene.refGene, ExonicFunc.refGene, 
				 AAChange.refGene, cosmic97_coding, avsnp150, gnomad40_exome_AF, CLNALLELEID, CLNSIG)

	names(df) <- c("Sample_name", "Sample_type", "Date_collected", "Chrom", "Position", "Ref", "Alt", "Type", "VAF", "Depth", "Alt_forward", "Ref_forward", "Alt_reverse", "Ref_reverse", "Function", "Gene", "Consequence", 
				   "AAchange", "cosmic97_coding", "avsnp150", "gnomad40_exome_AF", "CLNALLELEID", "CLNSIG")
	return(df)

}

return_anno_output_mutect_somatic <- function(PATH_variant_table) {

	print(PATH_variant_table)
	df <- as.data.frame(read.delim(PATH_variant_table, stringsAsFactors = FALSE, check.names = FALSE))
	df$Sample_name_n = str_split(colnames(df)[grep("WBC", colnames(df))][1], "\\.")[[1]][1]
	colnames(df) <- gsub("^.*[_-]WBC[-_].*\\.", "WBC_", colnames(df)) 
	colnames(df) <- gsub("^.*[_-]cfDNA[-_].*\\.", "cfDNA_", colnames(df)) 
	colnames(df) <- gsub("^.*[_-]FiT[-_].*\\.", "cfDNA_", colnames(df)) 

	df <- df %>%
		separate(col = cfDNA_SB, sep = ",", into = c("Ref_forward_t", "Ref_reverse_t", "Alt_forward_t", "Alt_reverse_t"), remove = TRUE) %>%
		separate(col = WBC_SB, sep = ",", into = c("Ref_forward_n", "Ref_reverse_n", "Alt_forward_n", "Alt_reverse_n"), remove = TRUE) %>%
		mutate(Patient_id = gsub("GU-", "", str_split(file_path_sans_ext(basename(PATH_variant_table)), "_")[[1]][1]),
			   Sample_name_t = file_path_sans_ext(basename(PATH_variant_table)),
			   Date_collected = tail(str_split(file_path_sans_ext(basename(PATH_variant_table)), "[-_]")[[1]], n = 1), 
			   TYPE = case_when(
							nchar(REF) > nchar(ALT) ~ "Deletion", 
							nchar(REF) < nchar(ALT) ~ "Insertion",
							nchar(REF) == nchar(ALT) ~ "SNV", 
							TRUE ~ "Error"), 
			   across(c(Alt_forward_t, Alt_reverse_t, Ref_forward_t, Ref_reverse_t, Alt_forward_n, Alt_reverse_n, Ref_forward_n, Ref_reverse_n), as.numeric), 
			   VAF_t = (Alt_forward_t + Alt_reverse_t)/(Ref_forward_t + Ref_reverse_t + Alt_forward_t + Alt_reverse_t), 
			   VAF_n = (Alt_forward_n + Alt_reverse_n)/(Ref_forward_n + Ref_reverse_n + Alt_forward_n + Alt_reverse_n), 
			   Depth_t = Ref_forward_t + Ref_reverse_t + Alt_forward_t + Alt_reverse_t, 
			   Depth_n = Ref_forward_n + Ref_reverse_n + Alt_forward_n + Alt_reverse_n, 
			   tumor_to_normal_VAF_ratio = VAF_t/VAF_n) %>%
		select(Patient_id, Sample_name_t, Date_collected, CHROM, POS, REF, ALT, TYPE, VAF_t, Depth_t, Alt_forward_t, Ref_forward_t, Alt_reverse_t, Ref_reverse_t, Func.refGene, Gene.refGene, ExonicFunc.refGene, 
		AAChange.refGene, cosmic97_coding, avsnp150, gnomad40_exome_AF, CLNALLELEID, CLNSIG, Sample_name_n, VAF_n, Depth_n, Alt_forward_n, Ref_forward_n, Alt_reverse_n, Ref_reverse_n, tumor_to_normal_VAF_ratio)

	names(df) <- c("Patient_id", "Sample_name_t", "Date_collected", "Chrom", "Position", "Ref", "Alt", "Type", "VAF_t", "Depth_t", "Alt_forward_t", "Ref_forward_t", "Alt_reverse_t", "Ref_reverse_t", "Function", "Gene", "Consequence", 
	"AAchange", "cosmic97_coding", "avsnp150", "gnomad40_exome_AF", "CLNALLELEID", "CLNSIG", "Sample_name_n", "VAF_n", "Depth_n", "Alt_forward_n", "Ref_forward_n", "Alt_reverse_n", "Ref_reverse_n", "tumor_to_normal_VAF_ratio")

	return(df)
}

return_anno_output_freebayes_chip <- function(PATH_variant_table) {
	print(PATH_variant_table)
	df <- as.data.frame(read.delim(PATH_variant_table, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE))
	colnames(df) <- gsub(paste0(file_path_sans_ext(basename(PATH_variant_table)), "."), "", colnames(df)) 
	df <- df[, -ncol(df)] # remove the last col, it just an empty column

	df <- df %>%
		rename(Ref_forward = SRF, Ref_reverse = SRR, Alt_forward = SAF, Alt_reverse = SAR) %>%
		mutate_at(c("POS", "Ref_forward", "Ref_reverse", "Alt_forward", "Alt_reverse"), as.numeric) %>%
		mutate(Sample = gsub(".tsv", "", basename(PATH_variant_table)), 
			   Sample_name = file_path_sans_ext(basename(PATH_variant_table)),
			   Sample_type = str_extract(basename(PATH_variant_table), "cfDNA|WBC"),
			   Date_collected = str_split(Sample_name, "[-_]")  %>% sapply(tail, 1),
			   nchar_ref = nchar(REF), 
			   nchar_alt = nchar(ALT), 
			   ALT = ifelse(grepl(",", ALT), sub(",.*", "", ALT), ALT)) %>% # Select the first element if ALT contains commas
		mutate(TYPE = case_when(nchar_ref > nchar_alt ~ "Deletion", 
			   					   nchar_ref < nchar_alt ~ "Insertion", 
								   nchar_ref == nchar_alt ~ "SNV"), 
			   VAF = (Alt_forward+Alt_reverse)/(Ref_forward+Ref_reverse+Alt_forward+Alt_reverse), 
			   Depth = Ref_forward+Ref_reverse+Alt_forward+Alt_reverse) %>%
		select(-nchar_ref, -nchar_alt) %>%
		select(Sample_name, Sample_type, Date_collected, CHROM, POS, REF, ALT, TYPE, VAF, Depth, Alt_forward, Ref_forward, Alt_reverse, Ref_reverse, Func.refGene, Gene.refGene, ExonicFunc.refGene, 
				 AAChange.refGene, cosmic97_coding, avsnp150, gnomad40_exome_AF, CLNALLELEID, CLNSIG)

		names(df) <- c("Sample_name", "Sample_type", "Date_collected", "Chrom", "Position", "Ref", "Alt", "Type", "VAF", "Depth", "Alt_forward", "Ref_forward", "Alt_reverse", "Ref_reverse", "Function", "Gene", "Consequence", 
		"AAchange", "cosmic97_coding", "avsnp150", "gnomad40_exome_AF", "CLNALLELEID", "CLNSIG")
	
	return(df)
}

return_anno_output_vardict_somatic <- function(PATH_variant_table) {

	df <- as.data.frame(read.delim(PATH_variant_table, stringsAsFactors = FALSE, check.names = FALSE))
	df$Sample_name_n = str_split(colnames(df)[grep("WBC", colnames(df))][1], "\\.")[[1]][1]
	colnames(df) <- gsub("^.*[_-]WBC[-_].*\\.", "WBC_", colnames(df)) 
	colnames(df) <- gsub("^.*[_-]cfDNA[-_].*\\.", "cfDNA_", colnames(df)) 
	colnames(df) <- gsub("^.*[_-]FiT[-_].*\\.", "cfDNA_", colnames(df)) 

	# Indicate insertion, deletion or SNV
	df <- df %>% mutate(Patient_id = gsub("GU-", "", str_split(file_path_sans_ext(basename(PATH_variant_table)), "_")[[1]][1]),
						Sample_name_t = file_path_sans_ext(basename(PATH_variant_table)),
						Date_collected = tail(str_split(file_path_sans_ext(basename(PATH_variant_table)), "[-_]")[[1]], n = 1), 
						TYPE = case_when(
							nchar(REF) > nchar(ALT) ~ "Deletion", 
							nchar(REF) < nchar(ALT) ~ "Insertion",
							nchar(REF) == nchar(ALT) ~ "SNV", 
							TRUE ~ "Error"), 
		   		 		VAF_t = cfDNA_VD/cfDNA_DP, 
						VAF_n = WBC_VD/WBC_DP, 
						tumor_to_normal_VAF_ratio = VAF_t/VAF_n) %>%
				 filter(STATUS %in% c("StrongSomatic", "LikelySomatic")) %>%
				 rename(Depth_t = cfDNA_DP, Depth_n = WBC_DP) %>%
				 select(-cfDNA_VD, -WBC_VD, -cfDNA_AF, -WBC_AF) %>%
		  		 separate(cfDNA_ALD, into = c("Alt_forward_t", "Alt_reverse_t"), sep = ",", remove = TRUE) %>%
		  		 separate(cfDNA_RD, into = c("Ref_forward_t", "Ref_reverse_t"), sep = ",", remove = TRUE) %>%
		  		 separate(WBC_ALD, into = c("Alt_forward_n", "Alt_reverse_n"), sep = ",", remove = TRUE) %>%
		  		 separate(WBC_RD, into = c("Ref_forward_n", "Ref_reverse_n"), sep = ",", remove = TRUE) %>%
				 select(Patient_id, Sample_name_t, Date_collected, STATUS, CHROM, POS, REF, ALT, TYPE, VAF_t, Depth_t, Alt_forward_t, Ref_forward_t, Alt_reverse_t, Ref_reverse_t, Func.refGene, Gene.refGene, ExonicFunc.refGene, 
				 AAChange.refGene, cosmic97_coding, avsnp150, gnomad40_exome_AF, CLNALLELEID, CLNSIG, Sample_name_n, VAF_n, Depth_n, Alt_forward_n, Ref_forward_n, Alt_reverse_n, Ref_reverse_n, tumor_to_normal_VAF_ratio) %>%
				 mutate(across(c(VAF_t, Depth_t, Alt_forward_t, Alt_reverse_t, Ref_forward_t, Ref_reverse_t, VAF_n, Depth_n, Alt_forward_n, Alt_reverse_n, Ref_forward_n, Ref_reverse_n), as.numeric))

	names(df) <- c("Patient_id", "Sample_name_t", "Date_collected", "Status", "Chrom", "Position", "Ref", "Alt", "Type", "VAF_t", "Depth_t", "Alt_forward_t", "Ref_forward_t", "Alt_reverse_t", "Ref_reverse_t", "Function", "Gene", "Consequence", 
	"AAchange", "cosmic97_coding", "avsnp150", "gnomad40_exome_AF", "CLNALLELEID", "CLNSIG", "Sample_name_n", "VAF_n", "Depth_n", "Alt_forward_n", "Ref_forward_n", "Alt_reverse_n", "Ref_reverse_n", "tumor_to_normal_VAF_ratio")

	return(df)
}

return_anno_output_vardict_chip <- function(PATH_variant_table) {
	print(PATH_variant_table)
	df <- as.data.frame(read.delim(PATH_variant_table, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE))
	colnames(df) <- gsub(paste0(file_path_sans_ext(basename(PATH_variant_table)), "."), "", colnames(df)) 

	# Indicate insertion, deletion or SNV
	df <- df %>% mutate(Sample_type = str_extract(basename(PATH_variant_table), "cfDNA|WBC"),
						Sample_name = file_path_sans_ext(basename(PATH_variant_table)),
						Date_collected = str_split(Sample_name, "[-_]")  %>% sapply(tail, 1), 
						TYPE = case_when(
							nchar(REF) > nchar(ALT) ~ "Deletion", 
							nchar(REF) < nchar(ALT) ~ "Insertion",
							nchar(REF) == nchar(ALT) ~ "SNV", 
							TRUE ~ "Error"), 
		   		 		VAF = VD/DP) %>%
				 rename(Depth = DP) %>%
				 select(-VD, -AF) %>%
		  		 separate(ALD, into = c("Alt_forward", "Alt_reverse"), sep = ",", remove = TRUE) %>%
		  		 separate(RD, into = c("Ref_forward", "Ref_reverse"), sep = ",", remove = TRUE) %>%
				 select(Sample_name, Sample_type, Date_collected, CHROM, POS, REF, ALT, TYPE, VAF, Depth, Alt_forward, Ref_forward, Alt_reverse, Ref_reverse, Func.refGene, Gene.refGene, ExonicFunc.refGene, 
				 AAChange.refGene, cosmic97_coding, avsnp150, gnomad40_exome_AF, CLNALLELEID, CLNSIG) %>%
				 mutate(across(c(VAF, Depth, Alt_forward, Alt_reverse, Ref_forward, Ref_reverse), as.numeric))

	names(df) <- c("Sample_name", "Sample_type", "Date_collected", "Chrom", "Position", "Ref", "Alt", "Type", "VAF", "Depth", "Alt_forward", "Ref_forward", "Alt_reverse", "Ref_reverse", "Function", "Gene", "Consequence", 
	"AAchange", "cosmic97_coding", "avsnp150", "gnomad40_exome_AF", "CLNALLELEID", "CLNSIG")
		  
	return(df)
}



