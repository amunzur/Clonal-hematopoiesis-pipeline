args <- commandArgs(trailingOnly=TRUE) # right now only arg needed is the name of the cohort 

variant_caller <- "Lofreq"

library(tidyverse)
library(stringr)
library(janitor)
library(pkgcond)

PATH_utilities_variant_caller <- file.path("/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/analysis", paste0("UTILITIES_filter_", variant_caller, ".R"))
PATH_utilities_file <- "/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/analysis/UTILITIES.R"
source(PATH_utilities_variant_caller) 
source(PATH_utilities_file)

setwd("/groups/wyattgrp/users/amunzur/pipeline/results")


################################
# SOMATIC
################################
max_VAF <- 0.40
min_VAF <- 0.01
max_SBF_cfDNA <- 0.0001
min_alt_reads_cfDNA <- 3
min_depth_gDNA <- 3
min_VAF_bg_ratio <- 15
min_test_ref_ratio <- 5
DIR_variant_tables_somatic <- file.path("data/variant_tables/somatic", variant_caller, "SSCS1")
DIR_variant_tables_chip <- file.path("data/variant_tables/chip", variant_caller, "SSCS1")
DIR_mpileup <- "metrics/mpileup/SSCS1"
DIR_mpileup_filtered_somatic <- file.path("metrics/mpileup_filtered/somatic", variant_caller)
DIR_mpileup_filtered_chip <- file.path("metrics/mpileup_filtered/chip", variant_caller)
PATH_bg <- "/groups/wyattgrp/users/jbacon/err_rates/chip/bgerror/chip_errorrate.tsv"

PATH_before_filtering_somatic <- file.path("variant_calling_final", variant_caller, "somatic_before_filtering.csv")
PATH_after_filtering_somatic <- file.path("variant_calling_final", variant_caller, "somatic_after_filtering.csv")

################################
# CHIP
################################
min_alt_reads <- 3
min_depth <- 20
max_VAF <- 0.40

PATH_before_filtering_chip <- file.path("variant_calling_final", variant_caller, "chip_before_filtering_CHIP.csv")
PATH_after_filtering_chip <- file.path("variant_calling_final", variant_caller, "chip_after_filtering_CHIP.csv")

bg <- read_delim(PATH_bg, delim = "\t")

vars_somatic <- MAINsomatic(min_alt_reads_cfDNA,
						min_VAF_bg_ratio, 
						DIR_variant_tables_somatic,
						bg,
						PATH_before_filtering_somatic,
						PATH_after_filtering_somatic)
vars_chip <- MAINchip(min_alt_reads,
						min_depth,
						max_VAF,
						min_VAF_bg_ratio,
						bg,
						PATH_before_filtering_chip,
						PATH_after_filtering_chip)