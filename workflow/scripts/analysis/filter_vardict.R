# This R script filters and plots VarScan2 output. 

library(tidyverse)
library(stringr)

THRESHOLD_ExAC_ALL <- 0.005
VALUE_Func_refGene <- "intronic"
THRESHOLD_VarFreq <- 0.30
THRESHOLD_Reads2 <- 5
THRESHOLD_VAF_bg_ratio <- 10
DIR_vardict <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/Vardict/new_chip_panel_reformatted"
DIR_ANNOVAR <- "/groups/wyattgrp/users/amunzur/pipeline/results/data/annovar_outputs/vardict/new_chip_panel"
PATH_bg <- "/groups/wyattgrp/users/amunzur/pipeline/resources/bg_error_rate/bg_error.tsv"
PATH_bets <- "/groups/wyattgrp/users/amunzur/pipeline/resources/betastasis/CLEANED_mutations_no_germline_filter.tsv"
PATH_bed  <- "/groups/wyattgrp/users/amunzur/pipeline/resources/panel/1000012543_CHIP_Design_selection_results_Version2/capture_targets.bed"
DIR_depth_metrics <- "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/depth/new_chip_panel"
PATH_collective_depth_metrics <- "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/averaged_depth/new_chip_panel/averaged_depths.txt"

PATH_validated_variants <- "/groups/wyattgrp/users/amunzur/pipeline/resources/validated_variants/chip_muts_locations.tsv"
PATH_SAVE_chip_variants <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/Vardict/finalized/chip_variants.csv"
DIR_finland_bams <- "/groups/wyattgrp/data/bam/kidney"

PATH_utilities_file_vardict <- "/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/analysis/UTILITIES_filter_vardict.R"
PATH_utilities_file <- "/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/analysis/UTILITIES.R"

source(PATH_utilities_file_vardict) # functions specificic to vardict
source(PATH_utilities_file) # functions shared between vardict and varscan

bg <- read_delim(PATH_bg, delim = "\t") # background error rate file

variants <- MAIN(THRESHOLD_ExAC_ALL, 
				VALUE_Func_refGene,
				THRESHOLD_VarFreq,
				THRESHOLD_Reads2,
				THRESHOLD_VAF_bg_ratio,
				DIR_vardict,
				DIR_ANNOVAR,
				bg,
				PATH_bets,
				PATH_bed,
				DIR_depth_metrics,
				PATH_collective_depth_metrics, 
				DIR_finland_bams)

combine_and_save(variants,
				PATH_validated_variants, 
				PATH_SAVE_chip_variants)

# variants_df <- read_csv(PATH_SAVE_chip_variants)
# write_csv(variants, "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/Vardict/finalized/vars.csv") # snv + indel, csv

# snv <- MAIN(THRESHOLD_ExAC_ALL = THRESHOLD_ExAC_ALL, 
# 				VALUE_Func_refGene = VALUE_Func_refGene, 
# 				THRESHOLD_VarFreq = THRESHOLD_VarFreq, 
# 				THRESHOLD_Reads2 = THRESHOLD_Reads2, 
# 				THRESHOLD_VAF_bg_ratio = THRESHOLD_VAF_bg_ratio, 
# 				DIR_varscan_snv = DIR_varscan_snv,
# 				DIR_varscan_indel = DIR_varscan_indel,
# 				ANNOVAR_snv_output = ANNOVAR_snv_output,
# 				ANNOVAR_indel_output = ANNOVAR_indel_output,
# 				PATH_bg = PATH_bg,
# 				PATH_bets = PATH_bets, 
# 				PATH_bed = PATH_bed,
# 				DIR_depth_metrics = DIR_depth_metrics, 
# 				"snv")

# indel <- MAIN(THRESHOLD_ExAC_ALL = THRESHOLD_ExAC_ALL, 
# 				VALUE_Func_refGene = VALUE_Func_refGene, 
# 				THRESHOLD_VarFreq = THRESHOLD_VarFreq, 
# 				THRESHOLD_Reads2 = THRESHOLD_Reads2, 
# 				THRESHOLD_VAF_bg_ratio = THRESHOLD_VAF_bg_ratio, 
# 				DIR_varscan_snv = DIR_varscan_snv,
# 				DIR_varscan_indel = DIR_varscan_indel,
# 				ANNOVAR_snv_output = ANNOVAR_snv_output,
# 				ANNOVAR_indel_output = ANNOVAR_indel_output,
# 				PATH_bg = PATH_bg,
# 				PATH_bets = PATH_bets, 
# 				PATH_bed = PATH_bed,
# 				DIR_depth_metrics = DIR_depth_metrics, 
# 				"indel")

