# This R script filters and plots VarScan2 output. 

args <- commandArgs(trailingOnly=TRUE) # right now only arg needed is the name of the cohort 

library(tidyverse)
library(stringr)
library(janitor)
library(pkgcond)

VALUE_Func_refGene <- "intronic"
THRESHOLD_VarFreq <- 0.40
THRESHOLD_Reads2 <- 1
THRESHOLD_VAF_bg_ratio <- 10
DIR_ANNOVAR <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/data/annovar_outputs/mutato")
PATH_bg <- "/groups/wyattgrp/users/amunzur/pipeline/resources/bg_error_rate/bg_error.tsv"

PATH_panel_genes <- "/groups/wyattgrp/users/amunzur/pipeline/resources/panel/kidney/genes.tsv"

PATH_bed  <- "/groups/wyattgrp/users/amunzur/pipeline/resources/panel/1000012543_CHIP_Design_selection_results_Version2/capture_targets.bed"
DIR_depth_metrics <- "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/depth/SSCS3"
PATH_collective_depth_metrics <- "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/averaged_depth/SSCS3/averaged_depths.txt"
DIR_tnvstats <- "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/tnvstats/kidney_samples"
DIR_temp <- "/groups/wyattgrp/users/amunzur/pipeline/results/temp"
PATH_filter_tnvstats_script <- "/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/analysis/filter_tnvstats.sh"

PATH_SAVE_chip_variants <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/Mutato/finalized/chip_variants_consensus.csv")

# PATH_utilities_file_vardict <- "/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/analysis/UTILITIES_filter_vardict.R"
PATH_utilities_file <- "/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/analysis/UTILITIES.R"

source(PATH_utilities_file_vardict) # functions specificic to vardict
source(PATH_utilities_file) # functions shared between vardict and varscan

bg <- read_delim(PATH_bg, delim = "\t") # background error rate file

variants <- MAIN(
				VALUE_Func_refGene,
				THRESHOLD_VarFreq,
				THRESHOLD_Reads2,
				THRESHOLD_VAF_bg_ratio,
				DIR_ANNOVAR,
				bg,
				PATH_panel_genes,
				PATH_bed,
				DIR_depth_metrics,
				PATH_collective_depth_metrics, 
				DIR_temp,
				DIR_tnvstats,
				PATH_filter_tnvstats_script,
				variant_caller)

combine_and_save(variants,
				PATH_SAVE_chip_variants)
