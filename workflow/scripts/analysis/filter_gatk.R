# This R script filters and plots VarScan2 output. 

args <- commandArgs(trailingOnly=TRUE) # right now only arg needed is the name of the cohort 

library(tidyverse)
library(stringr)
library(janitor)
library(pkgcond)

variant_caller <- "Mutect2"

VALUE_Func_refGene <- "intronic"
THRESHOLD_VarFreq <- 0.40
THRESHOLD_Reads2 <- 1
THRESHOLD_VAF_bg_ratio <- 10
DIR_ANNOVAR <- "/groups/wyattgrp/users/amunzur/pipeline/results/data/annovar_outputs/Mutect2/SSCS1"
DIR_basecounts_DCS <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/base_counts/Mutect2/gDNA_SSCS1_cfDNA_DCS"
DIR_basecounts_SSCS1 <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/base_counts/Mutect2/gDNA_SSCS1_cfDNA_SSCS1"
DIR_mpileup <- "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/mpileup/temp_SSCS1"
DIR_mpileup_filtered <- "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/mpileup_filtered/Mutect"

PATH_bg <- "/groups/wyattgrp/users/amunzur/pipeline/resources/bg_error_rate/bg_error.tsv"

PATH_panel_genes <- "/groups/wyattgrp/users/amunzur/pipeline/resources/panel/kidney/genes.tsv"

PATH_bed  <- "/groups/wyattgrp/users/amunzur/pipeline/resources/panel/1000012543_CHIP_Design_selection_results_Version2/capture_targets.bed"
DIR_depth_metrics <- "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/depth/SSCS1"
PATH_collective_depth_metrics <- "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/averaged_depth/SSCS1/averaged_depths.txt"
DIR_temp <- "/groups/wyattgrp/users/amunzur/pipeline/results/temp"

PATH_blacklist <- "/groups/wyattgrp/users/amunzur/pipeline/resources/validated_variants/blacklisted_variants.csv"

PATH_before_filtering <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/Mutect2/finalized/SSCS1/before_filtering.csv"
PATH_after_filtering <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/Mutect2/finalized/SSCS1/after_filtering.csv"
PATH_final <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/Mutect2/finalized/SSCS1/chip_variants.csv"

PATH_utilities_file_mutect <- "/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/analysis/UTILITIES_filter_gatk.R"
PATH_utilities_file <- "/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/analysis/UTILITIES.R"

source(PATH_utilities_file_mutect) # functions specificic to vardict
source(PATH_utilities_file) # functions shared between vardict and varscan

chroms <- c(as.character(1:22), "X", "Y")
bg <- read_delim(PATH_bg, delim = "\t", col_types = cols(chrom = col_character())) %>%
	  filter(chrom %in% chroms)

variants <- MAIN(
				VALUE_Func_refGene,
				THRESHOLD_VarFreq,
				THRESHOLD_Reads2,
				THRESHOLD_VAF_bg_ratio,
				DIR_ANNOVAR,
				DIR_mpileup,
				DIR_mpileup_filtered,
				bg,
				PATH_panel_genes,
				PATH_bed,
				DIR_depth_metrics,
				PATH_collective_depth_metrics, 
				DIR_temp,
				variant_caller,
				PATH_blacklist, 
				PATH_before_filtering,
				PATH_after_filtering,
				PATH_final)