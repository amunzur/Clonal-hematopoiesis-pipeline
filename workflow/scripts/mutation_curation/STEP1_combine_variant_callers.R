library(tidyverse)
library(stringr)
library(matrixStats)

# This script combines the outcome of multiple variant callers: VarScan and Vardict. 
# It counts how many variant callers identified each given variant. 
# ca r_env_v2

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

# COMBINE ALL BLADDER AND KIDNEY IN ONE FILE
consensus <- "SSCS2"
type <- "chip" # chip or somatic
variant_callers <- c("Mutect2", "Vardict", "freebayes")
path_output <- sprintf("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/%s_%s.csv", type, consensus)
DIR_bams <- sprintf("/groups/wyattgrp/users/amunzur/pipeline/results/data/bam/%s_final", consensus)
df_list <- list()  # Initialize an empty list

for (variant_caller in variant_callers) {
  file_path <- sprintf("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/%s/finalized/%s/%s_final.csv", variant_caller, consensus, toupper(type))
  if (file.exists(file_path)) {
    df = read_csv(file_path)
    df$Position <- as.numeric(df$Position)
    double_cols <- sapply(df, function(x) inherits(x, "numeric") && identical(typeof(x), "double"))
    df[double_cols] <- lapply(df[double_cols], as.numeric)
    df_list <- append(df_list, list(df))
  } else {
    warning(sprintf("File not found: %s", file_path))
  }
}

combined_df <- combine_variant_callers(df_list)
combined_df$Consensus <- consensus
combined_df <- add_bam_path(combined_df, DIR_bams)
combined_df <- filter(combined_df, n_callers > 1)

to_exclude = c("20-323", "18-439", "18-500", "19-005", "19-097", "20-231", "21-187", "21-302", "22-320", "19-334", "23-098", "23-414", "22-563", "23-018", "23-036")
to_exclude_ctDNA = c("20-313", "21-184", "21-430") # oxidative damage

if (type == "chip") {
  filtered_df <- combined_df %>% filter(!Patient_id %in% to_exclude)
} else if (type == "somatic") {
  exclude_combined = c(to_exclude, to_exclude_ctDNA)
  filtered_df <- combined_df %>% filter(!Patient_id %in% exclude_combined)
}

write_csv(filtered_df, path_output)

# 