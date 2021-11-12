# This R script filters and plots VarScan2 output. 

library(tidyverse)
library(stringr)

THRESHOLD_ExAC_ALL <- 0.005
VALUE_Func_refGene <- "intronic"
THRESHOLD_VarFreq <- 0.30
THRESHOLD_Reads2 <- 5
THRESHOLD_VAF_bg_ratio <- 10
DIR_varscan_snv <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/VarScan2/snv/new_chip_panel"
DIR_varscan_indel <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/VarScan2/indel/new_chip_panel"
ANNOVAR_snv_output <- "/groups/wyattgrp/users/amunzur/pipeline/results/data/annovar_outputs/snv/new_chip_panel"
ANNOVAR_indel_output <- "/groups/wyattgrp/users/amunzur/pipeline/results/data/annovar_outputs/indel/new_chip_panel"
PATH_bg <- "/groups/wyattgrp/users/amunzur/pipeline/resources/bg_error_rate/bg_error.tsv"
PATH_bets <- "/groups/wyattgrp/users/amunzur/pipeline/resources/betastasis/CLEANED_mutations_no_germline_filter.tsv"
PATH_bed  <- "/groups/wyattgrp/users/amunzur/pipeline/resources/panel/1000012543_CHIP_Design_selection_results_Version2/capture_targets.bed"
DIR_depth_metrics <- "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/depth/new_chip_panel"

PATH_validated_variants <- "/groups/wyattgrp/users/amunzur/pipeline/resources/validated_variants/chip_muts_locations.tsv"
PATH_SAVE_chip_variants <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/VarScan2/finalized/chip_variants.csv"
DIR_finland_bams <- "/groups/wyattgrp/data/bam/kidney"

PATH_utilities_file_varscan <- "/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/analysis/UTILITIES_filter_varscan.R" # has the functions we use here 
PATH_utilities_file <- "/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/analysis/UTILITIES.R"

source(PATH_utilities_file_varscan) # functions specific to varscan
source(PATH_utilities_file) # functions shared between vardict and varscan


bg <- read_delim(PATH_bg, delim = "\t")

snv <- MAIN(THRESHOLD_ExAC_ALL, 
				VALUE_Func_refGene,
				THRESHOLD_VarFreq,
				THRESHOLD_Reads2,
				THRESHOLD_VAF_bg_ratio,
				DIR_varscan_snv,
				DIR_varscan_indel,
				ANNOVAR_snv_output,
				ANNOVAR_indel_output,
				bg,
				PATH_bets,
				PATH_bed,
				DIR_depth_metrics,
				DIR_finland_bams,
				"snv")

indel <- MAIN(THRESHOLD_ExAC_ALL, 
				VALUE_Func_refGene, 
				THRESHOLD_VarFreq, 
				THRESHOLD_Reads2, 
				THRESHOLD_VAF_bg_ratio, 
				DIR_varscan_snv,
				DIR_varscan_indel,
				ANNOVAR_snv_output,
				ANNOVAR_indel_output,
				bg,
				PATH_bets,
				PATH_bed,
				DIR_depth_metrics,
				DIR_finland_bams,
				"indel")

variants_chip <- combine_and_save(snv,
						indel, 
						PATH_validated_variants, 
						PATH_SAVE_chip_variants)