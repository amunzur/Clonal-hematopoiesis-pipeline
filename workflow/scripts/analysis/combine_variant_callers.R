library(tidyverse)
library(stringr)
library(matrixStats)

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

cohort_name <- "new_chip_panel"
varscan_PATH <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/VarScan2/finalized", cohort_name, "chip_variants.csv")
vardict_PATH <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/Vardict/finalized", cohort_name, "chip_variants.csv")
gatk_PATH <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/Haplotype_caller/finalized", cohort_name, "chip_variants.csv")

varscan <- read_csv(varscan_PATH)
vardict <- read_csv(vardict_PATH)
gatk <- read_csv(gatk_PATH)

varscan$variant_caller <- "Varscan"
vardict$variant_caller <- "Vardict"
gatk$variant_caller <- "gatk"

# subset to cols in gatk, that one doesn't have the extra columns
varscan <- varscan[, colnames(gatk)]
vardict <- vardict[, colnames(gatk)]

# rbind and drop duplicates, this helps curate a list of unique variants
combined <- rbind(varscan, vardict, gatk)
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
combined <- identify_varcaller(combined, gatk, "gatk")

# further modifications
counts_df <- as.matrix(combined[, (ncol(combined)-2):ncol(combined)])
combined <- select(combined, -varscan, -vardict, -gatk)

counts_df <- replace_na(counts_df, FALSE)
counts_vector <- as.vector(rowCounts(counts_df, value = TRUE)) # number of variant callers that called the variant

combined <- cbind(combined, as.data.frame(counts_df))
combined$n_callers <- counts_vector

# to save
PATH_to_save_csv <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/combined", cohort_name, "combined.csv")
PATH_to_save_tsv <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/combined", cohort_name, "combined.tsv")

dir.create(dirname(PATH_to_save_csv))
write_csv(combined, PATH_to_save_csv)
write_delim(combined, PATH_to_save_tsv, delim = "\t")