library(argparse)

parser <- ArgumentParser(description = 'Combines variant calls from the three callers.')

# Add arguments
# parser$add_argument('--depth', type = 'character', required = TRUE, help = '')
# parser$add_argument('--min_alt_reads', type = 'character', required = TRUE, help = '')
parser$add_argument('--dir_working', type = 'character', required = TRUE, help = 'Working directory. Not results dir.')
parser$add_argument('--mutation_type', type = 'character', required = TRUE, help = 'chip or somatic, case does not matter')
parser$add_argument('--consensus', type = 'character', required = FALSE, help = '')
parser$add_argument('--DIR_bams', type = 'character', required = TRUE, help = 'Location of the bams.')
parser$add_argument('--path_output', type = 'character', required = FALSE, help = 'Output file will be saved here.')
parser$add_argument('--WBC_only', action = 'store_true', help = 'Are these WBC only calls?')
parser$add_argument('--input_file_keyword', type = 'character', required = FALSE, help = 'If the input files contain a keyword, specify it here.')

args <- parser$parse_args()
print(args)

# Rscript /path/to/STEP2_combine_variant_callers.R \
# --dir_working /groups/wyattgrp/users/amunzur/hla_pipeline/ \
# --mutation_type somatic \
# --DIR_bams /groups/wyattgrp/users/amunzur/hla_pipeline/results/data/bam/sorted \
# --path_output /groups/wyattgrp/users/amunzur/hla_pipeline/results/variant_calling/ctdna_mutations.csv \
# --WBC_only \
# --input_file_keyword WBConly

# Access the arguments
# depth <- as.numeric(args$depth)
# min_alt_reads <- as.numeric(args$min_alt_reads)
dir_working <- args$dir_working
mutation_type <- toupper(args$mutation_type)
DIR_bams <- args$DIR_bams
WBC_only <- args$WBC_only
path_output <- args$path_output
input_file_keyword <- args$input_file_keyword

if (is.null(args$consensus)) {
  consensus <- ""
} else {
  consensus <- args$consensus
}

print(input_file_keyword)
############################################################################

library(tidyverse)
library(stringr)
library(matrixStats)

# It counts how many variant callers identified each given variant. 
# ca r_env_v2

path_utils <- sprintf("%s/src/combine_variant_callers_UTILITIES.R", dir_working)
source(path_utils)

if (mutation_type == "CHIP") {
  variant_callers <- c("Mutect2", "Vardict", "freebayes")
} else {
  variant_callers <- c("Mutect2", "Vardict")
}

# If path_output is not provided, set a default value
if (is.null(path_output)) {
  path_output <- sprintf("%s/results/variant_calling/%s_SSCS2.csv", dir_working, mutation_type)
}

if (is.null(input_file_keyword)) {
  input_file_keyword <- ""
}

# DIR_bams <- "/groups/wyattgrp/users/amunzur/lu_chip/results/bam/SSCS2_final"
df_list <- list()  # Initialize an empty list

for (variant_caller in variant_callers) {
  # file_path <- sprintf("/groups/wyattgrp/users/amunzur/pipeline/results/wbc_downsampling/depth_%s/variant_calling/%s/finalized/%s/min_alt_reads_%s_%s_final.csv", depth, variant_caller, consensus, min_alt_reads, toupper(type))
    if (mutation_type == "CHIP") {
      file_path <- sprintf("%s/results/variant_calling/%s/finalized/%s/min_alt_reads_5_CHIP_final%s.csv", dir_working, variant_caller, consensus, input_file_keyword)
    } else {
      file_path <- sprintf("%sresults/variant_calling/%s/finalized/%s/SOMATIC_final%s.csv", dir_working, variant_caller, consensus, input_file_keyword)
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

if (WBC_only) {
  combined_df <- combine_variant_callers_WBC_only_calls(df_list)
  combined_df$Consensus <- consensus
  combined_df <- add_bam_path_WBC_only(combined_df, DIR_bams)
} else {
  combined_df <- combine_variant_callers(df_list)
  combined_df$Consensus <- consensus
  combined_df <- add_bam_path(combined_df, DIR_bams)
}

combined_df <- filter(combined_df, n_callers > 1)

write_csv(combined_df, path_output)
print(path_output)