args <- commandArgs(trailingOnly=TRUE) # right now only arg needed is the name of the cohort 

# consensus <- "SSCS2"
# variant_caller <- "Mutect2"
# type <- "chip"

consensus <- NULL
variant_caller <- NULL
type <- NULL
forced_rerun <- TRUE  # Default value for forced_rerun

if (length(args) == 0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
} else if (length(args) > 1) {
  consensus <- args[1]
  variant_caller <- args[2]
  type <- args[3]

  # Check if noforce is provided by the user
  if ("--noforce" %in% args) {
    forced_rerun <- FALSE
  } else {
    # Check if forced_rerun is provided by the user
    if ("--forced_rerun" %in% args) {
      forced_rerun <- TRUE
    }
  }
}


log_file <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/logs_slurm/misc_logs", paste0(paste(consensus, variant_caller, type, sep = "_"), ".txt"))
sink(log_file, append = FALSE, split = TRUE)

cat("Running with the following arguments:", consensus, variant_caller, type, "\n")


library(tidyverse)
library(stringr)
library(janitor)
library(pkgcond)
library(epitools)
library(tools)
library(zoo)

setwd("/groups/wyattgrp/users/amunzur/pipeline/results")
source("/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/analysis/UTILITIES.R")

PATH_sample_list <- sprintf("/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/paired_samples.tsv") # must be paired
DIR_mpileup <- sprintf("metrics/mpileup/%s", consensus)
PATH_bg <- "/groups/wyattgrp/users/amunzur/therapy_chip/resources/bg_error/error_rate/error_rates.tsv"
PATH_sample_information <- "../resources/sample_lists/sample_information.tsv"
PATH_blacklist <- "../resources/validated_variants/blacklisted_variants.csv"

chroms <- paste0("chr", c(as.character(1:22), "X", "Y"))
bg <- read_delim(PATH_bg, delim = "\t", col_types = cols(CHROM = col_character())) %>%
	  filter(CHROM %in% chroms) %>%
	  select(-X1)

#########################################################################################
# CHIP
min_alt_reads <- 5
min_depth <- 200
min_VAF_low <- 0.0025
max_VAF_low <- 0.45
min_VAF_high <- 0.55
max_VAF_high <- 0.90
min_VAF_bg_ratio <- 10

DIR_variant_tables_chip <- sprintf("data/variant_tables/%s/%s/%s", tolower(type), variant_caller, consensus)
DIR_mpileup_filtered_chip <- sprintf("metrics/mpileup_filtered_%s/chip/%s", consensus, variant_caller)
PATH_before_filtering <- sprintf("variant_calling/%s/finalized/%s/CHIP_before_filtering.csv", variant_caller, consensus)
PATH_after_filtering <- sprintf("variant_calling/%s/finalized/%s/CHIP_after_filtering.csv", variant_caller, consensus)
PATH_final_chip <- sprintf("variant_calling/%s/finalized/%s/CHIP_final.csv", variant_caller, consensus)

#########################################################################################
# GERMLINE
min_VAF_low_GERMLINE <- 0.45
max_VAF_low_GERMLINE <- 0.55
min_VAF_high_GERMLINE <- 0.95
max_VAF_high_GERMLINE <- 1

PATH_germline <- sprintf("variant_calling/%s_germline.csv", consensus)

#########################################################################################
# SOMATIC 
min_tumor_to_normal_vaf_ratio <- 5
min_alt_reads_t <- 5
max_alt_reads_n <- 2
min_depth_n <- 25
min_VAF_low_somatic <- 0.01
min_VAF_bg_ratio <- 10

DIR_variant_tables_somatic <- sprintf("data/variant_tables/%s/%s/%s", tolower(type), variant_caller, consensus)
DIR_mpileup_filtered_somatic <- sprintf("metrics/mpileup_filtered_%s/somatic/%s", consensus, variant_caller)
PATH_before_filtering_somatic <- sprintf("variant_calling/%s/finalized/%s/SOMATIC_before_filtering.csv", variant_caller, consensus)
PATH_final_somatic <- sprintf("variant_calling/%s/finalized/%s/SOMATIC_final.csv", variant_caller, consensus)

if (toupper(type) == "CHIP"){
	if (forced_rerun) {
		vars <- parse_anno_output(DIR_variant_tables_chip, "chip", variant_caller, PATH_sample_list = PATH_sample_list)
		vars <- add_patient_information(vars, PATH_sample_information)
		vars <- add_bg_error_rate(vars, bg)
		vars <- add_AAchange_effect(vars, variant_caller)
		vars <- evaluate_strand_bias2(vars)
		write_csv(vars, PATH_before_filtering)
	} else {
		vars <- read_csv(PATH_before_filtering)
		# vars <- evaluate_strand_bias2(vars)
		# write_csv(vars, PATH_before_filtering)
	}
	vars <- filter_variants_chip_or_germline("chip", vars, min_alt_reads, min_depth, min_VAF_low, max_VAF_low, min_VAF_high, max_VAF_high, min_VAF_bg_ratio, PATH_blacklist, blacklist = TRUE)
	# vars <- add_N_fraction(vars, DIR_mpileup, DIR_mpileup_filtered_chip, force = FALSE)
	vars$Variant_caller <- variant_caller
	write_csv(vars, PATH_after_filtering)

	vars <- combine_tumor_wbc(vars)
	# vars <- vars %>% filter(Gene == "DNMT3A",)
	vars <- filter(vars, Strand_bias_fishers_n != TRUE & Strand_bias_fishers_t != TRUE)
	
	cat("Saving the final list to", PATH_final_chip)
	write_csv(vars, PATH_final_chip)

} else if (toupper(type) == "SOMATIC") {
	if (forced_rerun) {
		vars <- parse_anno_output(DIR_variant_tables_somatic, "somatic", variant_caller, PATH_sample_list = PATH_sample_list)
		vars <- add_patient_information_somatic(vars, PATH_sample_information)
		vars <- add_bg_error_rate(vars, bg)
		vars <- add_AAchange_effect(vars, variant_caller)
		vars <- evaluate_strand_bias2(vars)
		write_csv(vars, PATH_before_filtering_somatic)
	} else {
		vars <- read_csv(PATH_before_filtering_somatic)
		# vars <- evaluate_strand_bias2(vars)
		# write_csv(vars, PATH_before_filtering_somatic)
	}
	vars <- filter_somatic_variants(vars, min_alt_reads_t, max_alt_reads_n, min_depth_n, min_VAF_low_somatic, min_tumor_to_normal_vaf_ratio, min_VAF_bg_ratio, PATH_blacklist, blacklist = TRUE)
	table(vars$Gene)
	# vars <- add_N_fraction(vars, DIR_mpileup, DIR_mpileup_filtered_somatic, force = FALSE)
	vars$Variant_caller <- variant_caller
	cat("Saving the final list to", PATH_final_somatic)
	write_csv(vars, PATH_final_somatic)

 } else if (toupper(type) == "GERMLINE") {
	germline <- read_csv(PATH_before_filtering)
	germline <- filter_variants_chip_or_germline("germline", germline, min_alt_reads, min_depth, min_VAF_low_GERMLINE, max_VAF_low_GERMLINE, min_VAF_high_GERMLINE, max_VAF_high_GERMLINE, min_VAF_bg_ratio, PATH_blacklist, blacklist = FALSE)
	germline$Variant_caller <- variant_caller
	germline <- filter(germline, !Strand_bias_fishers, CLNSIG %in% c("Pathogenic/Likely_pathogenic", "Pathogenic"))
	cat("Saving the final list to", PATH_germline)
	write_csv(germline, PATH_germline)
 }

sink()
