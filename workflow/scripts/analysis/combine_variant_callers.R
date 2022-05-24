library(tidyverse)
library(stringr)
library(matrixStats)

# This script combines the outcome of multiple variant callers: VarScan and Vardict. 
# It counts how many variant callers identified each given variant. 
# TO DO: Will add GATK later, after I manage to run it on the first batch of samples (new_chip_panel)

process_dfs <- function(path_df, varcaller_name){

	cols_to_subset <- c(
		"Sample_name", 
		"Patient_ID",
		"Cohort_name", 
		"Chrom", 
		"Position", 
		"Ref", 
		"Alt", 
		"VAF", 
		"Ref_reads", 
		"Alt_reads", 
		"Function", 
		"Gene", 
		"AAchange", 
		"Protein_annotation", 
		"Effects", 
		"ExAC_ALL", 
		"Variant", 
		"Error_rate", 
		"VAF_bg_ratio", 
		"Total_reads", 
		"Duplicate", 
		"Depth"
	)
	
	df <- read_csv(path_df)
	df <- df[, cols_to_subset]
	df$variant_caller <- varcaller_name

	return(df)

}

cohort_name <- "new_chip_panel"
varscan_PATH <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/VarScan2/finalized", cohort_name, "chip_variants.csv")
vardict_PATH <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/Vardict/finalized", cohort_name, "chip_variants.csv")

varscan <- process_dfs(varscan_PATH, "varscan")
vardict <- process_dfs(vardict_PATH, "vardict")

df <- rbind(varscan, vardict)
df[which(duplicated(df$AAchange)), "variant_caller"] <- "Both"

# to save
PATH_to_save_csv <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/combined", cohort_name, "combined_vars.csv")
PATH_to_save_tsv <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/combined", cohort_name, "combined_vars.tsv")

dir.create(dirname(PATH_to_save_csv))
write_csv(df, PATH_to_save_csv)
write_delim(df, PATH_to_save_tsv, delim = "\t")