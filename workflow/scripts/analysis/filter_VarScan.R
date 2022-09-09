# This R script filters and plots VarScan2 output. 

args <- commandArgs(trailingOnly=TRUE) # right now only arg needed is the name of the cohort 

library(tidyverse)
library(stringr)
library(epitools)
library(janitor)
library(pkgcond)

variant_caller <- "varscan"

THRESHOLD_ExAC_ALL <- 0.005
VALUE_Func_refGene <- "intronic"
THRESHOLD_VarFreq <- 0.40
THRESHOLD_Reads2 <- 5
THRESHOLD_VAF_bg_ratio <- 10
DIR_varscan_snv <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/VarScan2/snv")
DIR_varscan_indel <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/VarScan2/indel")
DIR_annovar_snv <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/data/annovar_outputs/snv")
DIR_annovar_indel <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/data/annovar_outputs/indel")
PATH_bg <- "/groups/wyattgrp/users/amunzur/pipeline/resources/bg_error_rate/bg_error.tsv"

PATH_bets_somatic <- "/groups/wyattgrp/users/amunzur/pipeline/resources/betastasis/CLEANED_mutations_kidney_cancer_somatic.tsv"
PATH_bets_germline <- "/groups/wyattgrp/users/amunzur/pipeline/resources/betastasis/CLEANED_mutations_no_germline_filter.tsv"
PATH_panel_genes <- "/groups/wyattgrp/users/amunzur/pipeline/resources/panel/kidney/genes.tsv"

PATH_bed  <- "/groups/wyattgrp/users/amunzur/pipeline/resources/panel/1000012543_CHIP_Design_selection_results_Version2/capture_targets.bed"
DIR_depth_metrics <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/metrics/depth")
PATH_collective_depth_metrics <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/metrics/averaged_depth/PICARD_markdup/averaged_depths.txt")
DIR_tnvstats <- "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/tnvstats/kidney_samples"
DIR_temp <- "/groups/wyattgrp/users/amunzur/pipeline/results/temp"
PATH_filter_tnvstats_script <- "/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/analysis/filter_tnvstats.sh"

PATH_validated_variants <- "/groups/wyattgrp/users/amunzur/pipeline/resources/validated_variants/chip_muts_locations.tsv"
PATH_SAVE_chip_variants <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/VarScan2/finalized/chip_variants.csv")
DIR_finland_bams <- "/groups/wyattgrp/data/bam/kidney"

PATH_utilities_file_varscan <- "/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/analysis/UTILITIES_filter_varscan.R" # has the functions we use here 
PATH_utilities_file <- "/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/analysis/UTILITIES.R"

source(PATH_utilities_file_varscan) # functions specific to varscan
source(PATH_utilities_file) # functions shared between vardict and varscan

bg <- read_delim(PATH_bg, delim = "\t")

DIR_annovar <- DIR_annovar_snv
DIR_varscan <- DIR_varscan_snv

snv <- MAIN(
				THRESHOLD_ExAC_ALL, 
				VALUE_Func_refGene, 
				THRESHOLD_VarFreq, 
				THRESHOLD_Reads2, 
				THRESHOLD_VAF_bg_ratio, 
				DIR_varscan_snv, 
				DIR_annovar_snv, 
				bg,
				PATH_panel_genes,
				PATH_bed,
				DIR_depth_metrics,
				PATH_collective_depth_metrics,
				DIR_temp,
				DIR_tnvstats,
				PATH_filter_tnvstats_script,
				"snv",
				variant_caller)

indel <- MAIN(
				THRESHOLD_ExAC_ALL, 
				VALUE_Func_refGene, 
				THRESHOLD_VarFreq, 
				THRESHOLD_Reads2, 
				THRESHOLD_VAF_bg_ratio, 
				DIR_varscan_indel,
				DIR_annovar_indel,
				bg,
				PATH_panel_genes,
				PATH_bed,
				DIR_depth_metrics,
				PATH_collective_depth_metrics,
				DIR_temp,
				DIR_tnvstats,
				PATH_filter_tnvstats_script,
				"indel", 
				variant_caller)

variants_chip <- combine_and_save(snv,
						indel, 
						PATH_validated_variants, 
						PATH_SAVE_chip_variants)
