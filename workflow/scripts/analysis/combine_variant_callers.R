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
					select(Sample_name, Sample_type, Patient_ID, Cohort_name, Chrom, Position, Ref, Alt, Gene) %>%
					mutate(variant_caller = TRUE)

	names(varcaller_df)[ncol(varcaller_df)] <- paste(varcaller_name)
	combined <- left_join(combined, varcaller_df)

	return(combined)
}

combine_variant_callers <- function(cohort_name, PATH_varscan, PATH_vardict){

	varscan <- read_csv(varscan_PATH) %>% mutate(variant_caller = "varscan")
	vardict <- read_csv(vardict_PATH) %>% mutate(variant_caller = "vardict")

	# rbind and drop duplicates, this helps curate a list of unique variants
	combined <- rbind(varscan, vardict)
	combined <- combined %>% distinct(Sample_name,Sample_type, Patient_ID, Cohort_name, Chrom, Position, Ref, Alt, Gene, .keep_all= TRUE) %>%
				select(-variant_caller)

	combined <- identify_varcaller(combined, varscan, "varscan")
	combined <- identify_varcaller(combined, vardict, "vardict")

	# add 1 to all deletion positions identified by varscan
	combined <- combined %>% mutate(Position = case_when(grepl("-", Alt) & varscan == TRUE ~ Position + 1, TRUE ~ Position))

	# further modifications
	counts_df <- as.matrix(combined[, (ncol(combined)-1):ncol(combined)])
	combined <- select(combined, -varscan, -vardict)

	counts_df <- replace_na(counts_df, FALSE)
	counts_vector <- as.vector(rowCounts(counts_df, value = TRUE)) # number of variant callers that called the variant

	combined <- cbind(combined, as.data.frame(counts_df))
	combined$n_callers <- counts_vector

	dir.create(dirname(PATH_to_save_csv))
	write_csv(combined, PATH_to_save_csv)
	write_delim(combined, gsub("csv", "tsv", PATH_to_save_csv), delim = "\t")

}

# Makes one master df that contains only the variants found in both tumor and wbc samples
combine_tumor_wbc <- function(DIR_vars){

    file_paths <- as.list(grep("combined_vars.csv", list.files(DIR_vars, recursive = TRUE, full.names = TRUE), value = TRUE))

    vars <- as.data.frame(do.call(rbind, lapply(file_paths, read_csv)))
    
    tumor <- vars %>% filter(Sample_type == "Tumor")
    wbc <- vars %>% filter(Sample_type == "WBC") 

    combined <- inner_join(tumor, wbc, by = c("Patient_ID", "Chrom", "Position", "Ref", "Alt", "Function", "Gene", "AAchange", "Protein_annotation", "Effects", "ExAC_ALL", "Variant"))
    names(combined) <- gsub("\\.x", "_t", names(combined)) 
    names(combined) <- gsub("\\.y", "_n", names(combined)) 

    combined <- combined %>%
            mutate(tumor_wbc_vaf_ratio = round((VAF_t / VAF_n), 2), 
                    tumor_wbc_depth_ratio = round((Depth_t / Depth_n), 2))

    return(combined)
}

# adds the abs path to the bam files to make it easier to run IGV snapshots
add_bam_path <- function(variant_df, DIR_bams){

    PATH_tumor_bam <- file.path(DIR_bams, variant_df$Cohort_name_t, "SC_penalty", paste0(variant_df$Sample_name_t, ".bam"))
    PATH_wbc_bam <- file.path(DIR_bams, variant_df$Cohort_name_n, "SC_penalty", paste0(variant_df$Sample_name_n, ".bam"))

    variant_df <- variant_df %>% 
        mutate(Path_bam_t = PATH_tumor_bam, Path_bam_n = PATH_wbc_bam) %>% 
        relocate(Path_bam_t, .after = Sample_name_t) %>%
        relocate(Path_bam_n, .after = Sample_name_n)

    return(variant_df)
}


# COMBINE VAR CALLERS
cohort_name <- "batch5"
varscan_PATH <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/VarScan2/finalized/chip_variants.csv")
vardict_PATH <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/Vardict/finalized/chip_variants.csv")
PATH_to_save_csv <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/combined/combined_vars.csv")

combine_variant_callers(PATH_varscan, PATH_vardict, PATH_to_save_csv)

# COMBINE TUMOR AND WBC
DIR_vars <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/combined"
DIR_bams <- "/groups/wyattgrp/users/amunzur/pipeline/results/data/bam"
combined <- combine_tumor_wbc(DIR_vars)
combined <- add_bam_path(combined, DIR_bams)

DIR_output <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/combined"
write_csv(combined, file.path(DIR_output, "tumor_wbc.csv"))
write_delim(combined, file.path(DIR_output, "tumor_wbc.tsv"), delim = "\t")
