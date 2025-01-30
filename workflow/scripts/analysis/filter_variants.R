library(argparse)

parser <- ArgumentParser(description = 'Filter variants')

# Add arguments
parser$add_argument('--consensus', type = 'character', required = TRUE, help = 'Consensus (e.g. SSCS2)')
parser$add_argument('--caller', type = 'character', required = TRUE, help = 'Caller name (freebayes, Mutect2 or Vardict)')
parser$add_argument('--type', type = 'character', required = TRUE, help = 'Type (chip or somatic)')
parser$add_argument('--combine_tumor_and_wbc', type = 'character', required = TRUE, help = 'Boolean')
parser$add_argument('--min_alt_reads', type = 'integer', required = TRUE, help = 'Mininum number of alt reads')
parser$add_argument('--match_cfDNA_and_WBC_dates', type = 'character', required = TRUE, help = 'When combining cfDNA and WBC, do the dates need to match exactly?')
parser$add_argument('--forced_rerun', type = 'character', required = TRUE, help = 'Boolean.')
parser$add_argument('--working_dir', type = 'character', required = TRUE, help = 'Working directory where outputs will be saved.')
parser$add_argument('--exclude_nonpathogenic_germline', action = 'store_true', help = 'Flag to set exclude_nonpathogenic_germline to TRUE if provided.')
parser$add_argument('--do_not_impose_vaf_upper_limit', action = 'store_true', help = 'Flag to set do_not_impose_vaf_upper_limit to TRUE if provided. Does not filter based on VAF in CH mutation calling.')

args <- parser$parse_args()

# Convert string arguments to logical

# Access the arguments
consensus <- args$consensus
variant_caller <- args$caller
type <- args$type
min_alt_reads <- args$min_alt_reads
match_cfDNA_and_WBC_dates <- args$match_cfDNA_and_WBC_dates
working_dir <- args$working_dir
combine_tumor_and_wbc <- tolower(args$combine_tumor_and_wbc) == 'true'
forced_rerun <- tolower(args$forced_rerun) == 'true'
exclude_nonpathogenic_germline <- args$exclude_nonpathogenic_germline
do_not_impose_vaf_upper_limit <- args$do_not_impose_vaf_upper_limit

# Print arguments (for testing purposes)
cat("Consensus:", consensus, "\n")
cat("Caller:", variant_caller, "\n")
cat("Type:", type, "\n")
cat("combine_tumor_and_wbc:", combine_tumor_and_wbc, "\n")
cat("min_alt_reads:", min_alt_reads, "\n")
cat("match_cfDNA_and_WBC_dates:", match_cfDNA_and_WBC_dates, "\n")
cat("forced_rerun:", forced_rerun, "\n")
cat("working_dir:", working_dir, "\n")

# Minor fixes to account for different spellings:
if (toupper(variant_caller) == "FREEBAYES") {
	variant_caller <- "freebayes" 
} else if (toupper(variant_caller) == "MUTECT2") {
	variant_caller <- "Mutect2" 
} else if (toupper(variant_caller) == "VARDICT") {
	variant_caller <- "Vardict"
} 

# Example run:
# Rscript filter_variants.R \
#             --consensus SSCS2 \
#             --caller Mutect2 \
#             --type chip \
#             --combine_tumor_and_wbc FALSE \
#             --min_alt_reads 5 \
#			  --match_cfDNA_and_WBC_dates FALSE \
#             --forced_rerun TRUE \
#             --working_dir "/groups/wyattgrp/users/amunzur/pipeline/results/wbc_downsampling/depth_1500"

library(tidyverse)
library(stringr)
library(janitor)
library(pkgcond)
library(epitools)
library(tools)
library(zoo)

setwd(working_dir)
source("/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/analysis/UTILITIES.R")

PATH_sample_list <- sprintf("../resources/sample_lists/paired_samples.tsv") # must be paired
DIR_mpileup <- sprintf("metrics/mpileup/%s", consensus)
PATH_bg <- "/groups/wyattgrp/users/jbacon/err_rates/chip_panel_CHIP/error_rate/error_rates.tsv"
PATH_bg_ctDNA <- "/groups/wyattgrp/users/amunzur/pipeline/resources/bg_error_rate/bgerror/error_rates_ctDNA.tsv"

PATH_sample_information <- "../resources/sample_lists/sample_information.tsv"
PATH_blacklist <- "/groups/wyattgrp/users/amunzur/pipeline/resources/validated_variants/blacklisted_variants.csv"

chroms <- paste0("chr", c(as.character(1:22), "X", "Y"))
bg <- read_delim(PATH_bg, delim = "\t", col_types = cols(CHROM = col_character())) %>%
	  filter(CHROM %in% chroms)
	  
bg_ctDNA <- read_delim(PATH_bg_ctDNA, delim = "\t", col_types = cols(chrom = col_character())) %>%
	  rename(CHROM = chrom, POS = pos, REF = ref) %>%
	  filter(CHROM %in% chroms)

#########################################################################################
# CHIP
min_alt_reads <- min_alt_reads
min_depth <- 200
min_VAF_bg_ratio <- 10
min_VAF_low <- 0.0025
max_VAF_low <- 0.45
min_VAF_high <- 0.55
max_VAF_high <- 0.90

if (do_not_impose_vaf_upper_limit) {
	max_VAF_low <- 1
	path_suffix <- "_no_VAF_UL"
} else {
	path_suffix <- ""
}

DIR_variant_tables_chip <- sprintf("data/variant_tables/%s/%s/%s", tolower(type), variant_caller, consensus)
DIR_mpileup_filtered_chip <- sprintf("metrics/mpileup_filtered_%s/chip/%s", consensus, variant_caller)
PATH_before_filtering <- sprintf("variant_calling/%s/finalized/%s/CHIP_before_filtering.csv", variant_caller, consensus, min_alt_reads)
PATH_after_filtering <- sprintf("variant_calling/%s/finalized/%s/min_alt_reads_%s_CHIP_after_filtering%s.csv", variant_caller, consensus, min_alt_reads, path_suffix)
PATH_final_chip <- sprintf("variant_calling/%s/finalized/%s/min_alt_reads_%s_CHIP_final%s.csv", variant_caller, consensus, min_alt_reads, path_suffix)

if (!dir.exists(sprintf("variant_calling/%s/finalized/%s/", variant_caller, consensus))) {
  dir.create(sprintf("variant_calling/%s/finalized/%s/", variant_caller, consensus), recursive = TRUE)
}
#########################################################################################
# GERMLINE
min_VAF_low_GERMLINE <- 0.45
max_VAF_low_GERMLINE <- 0.55
min_VAF_high_GERMLINE <- 0.95
max_VAF_high_GERMLINE <- 1

PATH_germline <- sprintf("variant_calling/%s_germline.csv", consensus)

if (!dir.exists("variant_calling")) {
  dir.create("variant_calling", recursive = TRUE)
}

#########################################################################################
# SOMATIC 
min_tumor_to_normal_vaf_ratio <- 5
min_alt_reads_t <- 5
max_alt_reads_n <- 2
min_depth_n <- 25
min_VAF_low_somatic <- 0.0025
min_VAF_bg_ratio <- 10

DIR_variant_tables_somatic <- sprintf("data/variant_tables/%s/%s/%s", tolower(type), variant_caller, consensus)
DIR_mpileup_filtered_somatic <- sprintf("metrics/mpileup_filtered_%s/somatic/%s", consensus, variant_caller)
PATH_before_filtering_somatic <- sprintf("variant_calling/%s/finalized/%s/SOMATIC_before_filtering.csv", variant_caller, consensus)
PATH_final_somatic <- sprintf("variant_calling/%s/finalized/%s/SOMATIC_final.csv", variant_caller, consensus)

if (!dir.exists(sprintf("variant_calling/%s/finalized/%s", variant_caller, consensus))) {
  dir.create(sprintf("variant_calling/%s/finalized/%s", variant_caller, consensus), recursive = TRUE)
}

if (toupper(type) == "CHIP"){
	if (forced_rerun || !file.exists(PATH_before_filtering)) {
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
	vars <- filter_variants_chip_or_germline("chip", vars, min_alt_reads, min_depth, min_VAF_low, max_VAF_low, min_VAF_high, max_VAF_high, min_VAF_bg_ratio, PATH_blacklist, blacklist = TRUE, filter_by_min_depth = FALSE)
	# vars <- add_N_fraction(vars, DIR_mpileup, DIR_mpileup_filtered_chip, force = FALSE)
	vars$Variant_caller <- variant_caller
	write_csv(vars, PATH_after_filtering)
	
	if (combine_tumor_and_wbc) {
		vars <- combine_tumor_wbc(vars, match_cfDNA_and_WBC_dates)
		vars <- filter(vars, Strand_bias_fishers_n != TRUE & Strand_bias_fishers_t != TRUE)
	} else {
		vars <- filter(vars, Strand_bias_fishers != TRUE)
	}
	
	cat("Saving the final list to", PATH_final_chip)
	write_csv(vars, PATH_final_chip)

} else if (toupper(type) == "SOMATIC") {
	if (forced_rerun || !file.exists(PATH_before_filtering)) {
		vars <- parse_anno_output(DIR_variant_tables_somatic, "somatic", variant_caller, PATH_sample_list = PATH_sample_list)
		vars <- add_patient_information_somatic(vars, PATH_sample_information)
		vars <- add_bg_error_rate(vars, bg_ctDNA)
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
	germline <- filter(germline, !Strand_bias_fishers)
	if (exclude_nonpathogenic_germline) {
		germline <- filter(germline, CLNSIG %in% c("Pathogenic/Likely_pathogenic", "Pathogenic"))
	} 
	cat("Saving the final list to", PATH_germline)
	write_csv(germline, PATH_germline)
 }