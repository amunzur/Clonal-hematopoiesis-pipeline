library(argparse)

parser <- ArgumentParser(description = 'Combine variant callers')

# Add arguments
parser$add_argument('--depth', type = 'character', required = TRUE, help = '')
parser$add_argument('--min_alt_reads', type = 'character', required = TRUE, help = '')

args <- parser$parse_args()

# Convert string arguments to logical

# Access the arguments
depth <- as.numeric(args$depth)
min_alt_reads <- as.numeric(args$min_alt_reads)

############################################################################

library(tidyverse)
library(stringr)
library(matrixStats)

# This script combines the outcome of multiple variant callers: VarScan and Vardict. 
# It counts how many variant callers identified each given variant. 
# ca r_env_v2

source("/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/mutation_curation/combine_variant_callers_UTILITIES.R")

# COMBINE ALL BLADDER AND KIDNEY IN ONE FILE
consensus <- "SSCS2"
type <- "chip" # chip or somatic
variant_callers <- c("Mutect2", "Vardict")
# depth <- "depth_1500"
# min_alt_reads <- "min_alt_reads_5"

# path_output <- sprintf("/groups/wyattgrp/users/amunzur/pipeline/results/wbc_downsampling/depth_%s/variant_calling/min_alt_reads_%s_%s.csv", depth, min_alt_reads, type)
path_output <- sprintf("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/%s_SSCS2.csv", type)
DIR_bams <- "/groups/wyattgrp/users/amunzur/pipeline/results/data/bam/SSCS2_final"
df_list <- list()  # Initialize an empty list

for (variant_caller in variant_callers) {
  # file_path <- sprintf("/groups/wyattgrp/users/amunzur/pipeline/results/wbc_downsampling/depth_%s/variant_calling/%s/finalized/%s/min_alt_reads_%s_%s_final.csv", depth, variant_caller, consensus, min_alt_reads, toupper(type))
    if (type == "chip") {
      print("LOL")
      file_path <- sprintf("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/%s/finalized/SSCS2/min_alt_reads_5_CHIP_final.csv", variant_caller)
    } else {
      file_path <- sprintf("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/%s/finalized/SSCS2/SOMATIC_final.csv", variant_caller)
    }
  if (file.exists(file_path)) {
    df = read_csv(file_path)
    df$Position <- as.numeric(df$Position)
    double_cols <- sapply(df, function(x) inherits(x, "numeric") && identical(typeof(x), "double"))
    df[double_cols] <- lapply(df[double_cols], as.numeric)
    df_list <- append(df_list, list(df))
  } else {
    warning(sprintf("File not found: %s", file_path))
  }
}

# combined_df <- combine_variant_callers_WBC_only_calls(df_list)
combined_df <- combine_variant_callers(df_list)

combined_df$Consensus <- consensus
combined_df <- add_bam_path(combined_df, DIR_bams)
combined_df <- filter(combined_df, n_callers > 1)

# There are some samples with oxidative damage, for these ones only include mutations with VAF > 10%.
# These samples are 20-313, 21-184, 21-430, 20-265, but 20-313 and 20-265 are already dropped due to insufficient coverage.
# filtered_df <- combined_df %>%
#   filter(!(Patient_id %in% c("21-284", "21-430")) | (Patient_id %in% c("21-284", "21-430") & VAF_t >= 15))

# to_exclude = c("20-323", "18-439", "18-500", "19-005", "19-097", "20-231", "21-187", "21-302", "22-320", "19-334", "23-098", "23-414", "22-563", "23-018", "23-036")
to_exclude_ctDNA = c("20-313", "21-184", "21-430", "20-265") # oxidative damage

if (type == "somatic") {
  exclude_combined = to_exclude_ctDNA
  combined_df <- combined_df %>% filter(!Patient_id %in% exclude_combined)
}

write_csv(combined_df, path_output)