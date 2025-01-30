read_variants_df <- function(PATH, variant_caller, consensus_type){
	df <- read_csv(PATH)
	return(df)
}

identify_varcaller <- function(combined, varcaller_df, varcaller_name){

	# select cols to avoid duplication after merge
	varcaller_df <- varcaller_df %>%
						select(Sample_name_t, Sample_name_n, Chrom, Position, Ref, Alt, Gene, Protein_annotation) %>%
						mutate(Variant_caller = TRUE)

	names(varcaller_df)[ncol(varcaller_df)] <- paste(varcaller_name)
	combined <- left_join(combined, varcaller_df)

	return(combined)
}

identify_varcaller_WBC_only <- function(combined, varcaller_df, varcaller_name){

	# select cols to avoid duplication after merge
	varcaller_df <- varcaller_df %>%
						select(Sample_name, Chrom, Position, Ref, Alt, Gene, Protein_annotation) %>%
						mutate(Variant_caller = TRUE)

	names(varcaller_df)[ncol(varcaller_df)] <- paste(varcaller_name)
	combined <- left_join(combined, varcaller_df)

	return(combined)
}



combine_variant_callers <- function(df_list, additional_args = list()) {
  # Combine all input data frames
  combined <- dplyr::bind_rows(df_list) %>%
    dplyr::distinct(Sample_name_t, Sample_name_n, Chrom, Position, Ref, Alt, Gene, Protein_annotation, .keep_all = TRUE)

  # Identify the variant caller for each row
  caller_names <- list()
  for (df in df_list) {
    caller_name <- unique(df$Variant_caller)
    combined <- identify_varcaller(combined, df, caller_name)
	caller_names <- append(caller_names, caller_name)
  }

  # Further modifications
  counts_df <- as.matrix(combined[, unique(unlist(caller_names))])
  combined <- select(combined, -contains(unique(unlist(caller_names))))

  counts_df <- replace_na(counts_df, FALSE)
  counts_vector <- as.vector(rowCounts(counts_df, value = TRUE))  # Number of variant callers that called the variant

  # Add counts to the combined data frame
  combined <- cbind(combined, as.data.frame(counts_df))
  combined$n_callers <- counts_vector

  return(combined)
}

combine_variant_callers_WBC_only_calls <- function(df_list, additional_args = list()) {
  # Combine all input data frames
  combined <- dplyr::bind_rows(df_list) %>%
    dplyr::distinct(Sample_name, Chrom, Position, Ref, Alt, Gene, Protein_annotation, .keep_all = TRUE)

  # Identify the variant caller for each row
  caller_names <- list()
  for (df in df_list) {
    caller_name <- unique(df$Variant_caller)
    combined <- identify_varcaller_WBC_only(combined, df, caller_name)
	caller_names <- append(caller_names, caller_name)
  }

  # Further modifications
  counts_df <- as.matrix(combined[, unique(unlist(caller_names))])
  combined <- select(combined, -contains(unique(unlist(caller_names))))

  counts_df <- replace_na(counts_df, FALSE)
  counts_vector <- as.vector(rowCounts(counts_df, value = TRUE))  # Number of variant callers that called the variant

  # Add counts to the combined data frame
  combined <- cbind(combined, as.data.frame(counts_df))
  combined$n_callers <- counts_vector

  return(combined)
}



# adds the abs path to the bam files to make it easier to run IGV snapshots
add_bam_path <- function(variant_df, DIR_bams){

    PATH_tumor_bam <- file.path(DIR_bams, paste0(variant_df$Sample_name_t, ".bam"))
    PATH_wbc_bam <- file.path(DIR_bams, paste0(variant_df$Sample_name_n, ".bam"))

    variant_df <- variant_df %>% 
        mutate(Path_bam_t = PATH_tumor_bam, Path_bam_n = PATH_wbc_bam) %>% 
        relocate(Path_bam_t, .after = Sample_name_t) %>%
        relocate(Path_bam_n, .after = Sample_name_n)

    return(variant_df)
}

add_bam_path_WBC_only <- function(variant_df, DIR_bams){

    PATH_wbc_bam <- file.path(DIR_bams, paste0(variant_df$Sample_name, ".bam"))

    variant_df <- variant_df %>% 
        mutate(Path_bam_n = PATH_wbc_bam) %>% 
        relocate(Path_bam_n, .after = Sample_name)

    return(variant_df)
}