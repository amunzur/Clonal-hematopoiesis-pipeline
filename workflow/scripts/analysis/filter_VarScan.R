# This R script filters and plots VarScan2 output. 

library(tidyverse)
library(stringr)
library(epitools)
library(janitor)
library(pkgcond)

cohort_name <- "batch2"
variant_caller <- "varscan"

THRESHOLD_ExAC_ALL <- 0.005
VALUE_Func_refGene <- "intronic"
THRESHOLD_VarFreq <- 0.40
THRESHOLD_Reads2 <- 5
THRESHOLD_VAF_bg_ratio <- 10
DIR_varscan_snv <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/VarScan2/snv", cohort_name)
DIR_varscan_indel <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/VarScan2/indel", cohort_name)
ANNOVAR_snv_output <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/data/annovar_outputs/snv", cohort_name)
ANNOVAR_indel_output <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/data/annovar_outputs/indel", cohort_name)
PATH_bg <- "/groups/wyattgrp/users/amunzur/pipeline/resources/bg_error_rate/bg_error.tsv"

PATH_bets_somatic <- "/groups/wyattgrp/users/amunzur/pipeline/resources/betastasis/CLEANED_mutations_kidney_cancer_somatic.tsv"
PATH_bets_germline <- "/groups/wyattgrp/users/amunzur/pipeline/resources/betastasis/CLEANED_mutations_no_germline_filter.tsv"
PATH_panel_genes <- "/groups/wyattgrp/users/amunzur/pipeline/resources/panel/kidney/genes.tsv"

PATH_bed  <- "/groups/wyattgrp/users/amunzur/pipeline/resources/panel/1000012543_CHIP_Design_selection_results_Version2/capture_targets.bed"
DIR_depth_metrics <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/metrics/depth", cohort_name)
PATH_collective_depth_metrics <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/metrics/averaged_depth/PICARD_markdup", cohort_name, "averaged_depths.txt")
DIR_tnvstats <- "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/tnvstats/kidney_samples"
DIR_temp <- "/groups/wyattgrp/users/amunzur/pipeline/results/temp"
PATH_filter_tnvstats_script <- "/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/analysis/filter_tnvstats.sh"

PATH_validated_variants <- "/groups/wyattgrp/users/amunzur/pipeline/resources/validated_variants/chip_muts_locations.tsv"
PATH_SAVE_chip_variants <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/VarScan2/finalized", cohort_name, "chip_variants.csv")
DIR_finland_bams <- "/groups/wyattgrp/data/bam/kidney"

PATH_utilities_file_varscan <- "/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/analysis/UTILITIES_filter_varscan.R" # has the functions we use here 
PATH_utilities_file <- "/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/analysis/UTILITIES.R"

source(PATH_utilities_file_varscan) # functions specific to varscan
source(PATH_utilities_file) # functions shared between vardict and varscan

bg <- read_delim(PATH_bg, delim = "\t")

variant_type <- "snv"
snv <- MAIN(cohort_name, 
				THRESHOLD_ExAC_ALL, 
				VALUE_Func_refGene,
				THRESHOLD_VarFreq,
				THRESHOLD_Reads2,
				THRESHOLD_VAF_bg_ratio,
				DIR_varscan_snv,
				DIR_varscan_indel,
				ANNOVAR_snv_output,
				ANNOVAR_indel_output,
				bg,
				PATH_bets_somatic,
				PATH_bets_germline,
				PATH_panel_genes,
				PATH_bed,
				DIR_depth_metrics,
				PATH_collective_depth_metrics,
				DIR_finland_bams,
				DIR_temp,
				DIR_tnvstats,
				PATH_filter_tnvstats_script,
				variant_type, 
				variant_caller)

variant_type <- "indel"
indel <- MAIN(cohort_name, 
				THRESHOLD_ExAC_ALL, 
				VALUE_Func_refGene, 
				THRESHOLD_VarFreq, 
				THRESHOLD_Reads2, 
				THRESHOLD_VAF_bg_ratio, 
				DIR_varscan_snv,
				DIR_varscan_indel,
				ANNOVAR_snv_output,
				ANNOVAR_indel_output,
				bg,
				PATH_bets_somatic,
				PATH_bets_germline,
				PATH_panel_genes,
				PATH_bed,
				DIR_depth_metrics,
				PATH_collective_depth_metrics,
				DIR_finland_bams,
				DIR_temp,
				DIR_tnvstats,
				PATH_filter_tnvstats_script,
				variant_type, 
				variant_caller)

variants_chip <- combine_and_save(snv,
						indel, 
						PATH_validated_variants, 
						PATH_SAVE_chip_variants)
