args <- commandArgs(trailingOnly=TRUE) # right now only arg needed is the name of the cohort 

library(tidyverse)
library(stringr)
library(janitor)
library(pkgcond)
library(epitools)
library(tools)

PATH_utilities_file_mutect <- "/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/analysis/UTILITIES_filter_gatk.R"
PATH_utilities_file <- "/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/analysis/UTILITIES.R"
source(PATH_utilities_file_mutect) 
source(PATH_utilities_file)

setwd("/groups/wyattgrp/users/amunzur/pipeline/results")

diagnosis <- "kidney"
variant_caller <- "Mutect2"
consensus <- "SSCS2"

PATH_sample_list <- sprintf("/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/paired_samples_%s.tsv", tolower(diagnosis)) # must be paired
# DIR_variant_tables_chip <- sprintf("data/variant_tables/chip/%s/test", variant_caller)
DIR_mpileup <- sprintf("metrics/mpileup/%s", consensus)
PATH_bg <- "/groups/wyattgrp/users/jbacon/err_rates/new_chip/error_rate/chip_error_rates.tsv"
PATH_sample_information <- "../resources/sample_lists/sample_information.tsv"
PATH_blacklist <- "../resources/validated_variants/blacklisted_variants.csv"

chroms <- paste0("chr", c(as.character(1:22), "X", "Y"))
bg <- read_delim(PATH_bg, delim = "\t", col_types = cols(CHROM = col_character())) %>%
	  filter(CHROM %in% chroms) %>%
	  select(-X1)

################################
# CHIP
min_alt_reads <- 5
min_depth <- 200
min_VAF_low <- 0.001
max_VAF_low <- 0.45
min_VAF_high <- 0.55
max_VAF_high <- 0.90
min_VAF_bg_ratio <- 10

type <- "chip"

DIR_variant_tables_chip <- sprintf("data/variant_tables/%s/%s/%s", tolower(type), variant_caller, consensus)
DIR_mpileup_filtered_chip <- sprintf("metrics/mpileup_filtered_%s/chip/%s", consensus, variant_caller)
PATH_before_filtering <- sprintf("variant_calling/%s/finalized/%s/%s_before_filtering.csv", variant_caller, consensus, toupper(diagnosis))
PATH_after_filtering <- sprintf("variant_calling/%s/finalized/%s/%s_after_filtering.csv", variant_caller, consensus, toupper(diagnosis))
PATH_final <- sprintf("variant_calling/%s/finalized/%s/%s_final.csv", variant_caller, consensus, toupper(diagnosis))

vars <- parse_anno_output(DIR_variant_tables_chip, "chip", PATH_sample_list = PATH_sample_list) # if PATH_sample_list is null, loads all files in the dir.
vars <- add_patient_information(vars, PATH_sample_information)
vars <- add_bg_error_rate(vars, bg)
vars <- add_AAchange_effect(vars, variant_caller)
vars <- evaluate_strand_bias(vars)
write_csv(vars, PATH_before_filtering)

vars <- filter_variants_chip_or_germline("chip", vars, min_alt_reads, min_depth, min_VAF_low, max_VAF_low, min_VAF_high, max_VAF_high, min_VAF_bg_ratio, PATH_blacklist, blacklist = TRUE)
vars <- add_N_fraction(vars, DIR_mpileup, DIR_mpileup_filtered_chip, force = FALSE)
vars$Variant_caller <- variant_caller
write_csv(vars, PATH_after_filtering)
	
vars <- combine_tumor_wbc(vars)
vars <- filter(vars, Strand_bias_fishers_n != TRUE & Strand_bias_fishers_t != TRUE)
write_csv(vars, PATH_final)

#########################################################################################
# SOMATIC MUTATIONS
min_tumor_to_normal_vaf_ratio <- 5
min_alt_reads_t <- 5
max_alt_reads_n <- 2
min_depth_n <- 25
min_VAF_low <- 0.1
min_VAF_bg_ratio <- 10

type <- "somatic"

DIR_variant_tables_somatic <- sprintf("data/variant_tables/%s/%s/%s", tolower(type), variant_caller, consensus)
DIR_mpileup_filtered_somatic <- sprintf("metrics/mpileup_filtered_%s/somatic/%s", consensus, variant_caller)
PATH_before_filtering_somatic <- sprintf("variant_calling/%s/finalized/%s/%s_before_filtering_somatic.csv", variant_caller, consensus, toupper(diagnosis))
PATH_final_somatic <- sprintf("variant_calling/%s/finalized/%s/%s_final_somatic.csv", variant_caller, consensus, toupper(diagnosis))

vars <- parse_anno_output(DIR_variant_tables_somatic, "somatic", variant_caller, PATH_sample_list = PATH_sample_list)
vars <- add_patient_information_somatic(vars, PATH_sample_information)
vars <- add_bg_error_rate(vars, bg)
vars <- add_AAchange_effect(vars, "Mutect2")
vars <- evaluate_strand_bias(vars)
write_csv(vars, PATH_before_filtering_somatic)

vars <- filter_somatic_variants(vars, min_alt_reads_t, max_alt_reads_n, min_depth_n, min_VAF_low, min_tumor_to_normal_vaf_ratio, min_VAF_bg_ratio, PATH_blacklist, blacklist = TRUE)
vars <- add_N_fraction(vars, DIR_mpileup, DIR_mpileup_filtered_somatic, force = FALSE)
vars$Variant_caller <- "Mutect2"
write_csv(vars, PATH_final_somatic)