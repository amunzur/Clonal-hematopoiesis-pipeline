library(tidyverse)
library(stringr)
library(matrixStats)

# This script combines the outcome of multiple variant callers: VarScan and Vardict. 
# It counts how many variant callers identified each given variant. 
# TO DO: Will add GATK later, after I manage to run it on the first batch of samples (new_chip_panel)

# To help with the joins later on
identify_varcaller <- function(combined, varcaller_df, varcaller_name){

	# select cols to avoid duplication after merge
	varcaller_df <- varcaller_df %>%
					select(
						Sample_name,
						Sample_type, 
						Patient_ID, 
						Cohort_name, 
						Chrom, 
						Position, 
						Ref, 
						Alt, 
						Gene) %>%
					mutate(variant_caller = TRUE)

	names(varcaller_df)[ncol(varcaller_df)] <- paste(varcaller_name)
	combined <- left_join(combined, varcaller_df)

	return(combined)
}


cohort_name <- "batch5"
varscan_PATH <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/VarScan2/finalized", cohort_name, "chip_variants.csv")
vardict_PATH <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/Vardict/finalized", cohort_name, "chip_variants.csv")

varscan <- read_csv(varscan_PATH) %>% mutate(variant_caller = "varscan")
vardict <- read_csv(vardict_PATH) %>% mutate(variant_caller = "vardict")

# rbind and drop duplicates, this helps curate a list of unique variants
combined <- rbind(varscan, vardict)
combined <- combined %>%
			distinct(
				Sample_name,
				Sample_type, 
				Patient_ID, 
				Cohort_name, 
				Chrom, 
				Position, 
				Ref, 
				Alt, 
				Gene,
				.keep_all= TRUE) %>%
			select(-variant_caller)

combined <- identify_varcaller(combined, varscan, "varscan")
combined <- identify_varcaller(combined, vardict, "vardict")

# add 1 to all deletion positions identified by varscan
combined <- combined %>%
				mutate(
					Position = case_when(
						grepl("-", Alt) & varscan == TRUE ~ Position + 1, 
						TRUE ~ Position))

# further modifications
counts_df <- as.matrix(combined[, (ncol(combined)-1):ncol(combined)])
combined <- select(combined, -varscan, -vardict)

counts_df <- replace_na(counts_df, FALSE)
counts_vector <- as.vector(rowCounts(counts_df, value = TRUE)) # number of variant callers that called the variant

combined <- cbind(combined, as.data.frame(counts_df))
combined$n_callers <- counts_vector

# to save
PATH_to_save_csv <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/combined", cohort_name, "combined_vars.csv")
PATH_to_save_tsv <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/combined", cohort_name, "combined_vars.tsv")

dir.create(dirname(PATH_to_save_csv))
write_csv(combined, PATH_to_save_csv)
write_delim(combined, PATH_to_save_tsv, delim = "\t")