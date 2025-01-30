"""
What would our CH calls look like if we did CH calling only on the WBC?
"""

library(tidyverse)
library(stringr)
library(matrixStats)

source("/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/mutation_curation/combine_variant_callers_UTILITIES.R")

library(data.table)

var_callers <- c("freebayes", "Mutect2", "Vardict")
dir_vars <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling"
path_sample_info <- "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"
DIR_bams <- "/groups/wyattgrp/users/amunzur/pipeline/results/data/bam/SSCS2_final"
chip_muts_to_exclude <- "/groups/wyattgrp/users/amunzur/pipeline/resources/validated_variants/chip_to_exclude_IGV.csv"

sample_df <- unique(fread(path_sample_info, sep="\t", header=TRUE, col.names=c("Patient_id", "Date", "Diagnosis", "Timepoint"))[, c("Patient_id"), with=FALSE])

df_list <- list()

for (var_caller in var_callers) {
  path_calls <- file.path(dir_vars, var_caller, "finalized/SSCS2", "CHIP_after_filtering.csv")
  df <- read.csv(file = path_calls) %>%
    filter(Sample_type == "WBC",
           !Strand_bias_fishers,
           (VAF >= 0.0050 & VAF <= 0.45) | (VAF >= 0.55 & VAF <= 0.9),
           VAF_bg_ratio > 20, 
           Timepoint == "Baseline") %>%
    mutate(Variant_caller = var_caller)

  df$Position <- as.numeric(df$Position)
  double_cols <- sapply(df, function(x) inherits(x, "numeric") && identical(typeof(x), "double"))
  df[double_cols] <- lapply(df[double_cols], as.numeric)

  df <- merge(df, sample_df, by="Patient_id", all.x=TRUE)
  df_list[[var_caller]] <- df
}

combined_df <- combine_variant_callers_WBC_only_calls(df_list)
combined_df$Consensus <- "SSCS2"
combined_df <- add_bam_path_WBC_only(combined_df, DIR_bams)
combined_df <- filter(combined_df, n_callers > 1)
combined_df <- filter(combined_df, Timepoint == "Baseline")


# exclude blacklisted events
to_exclude <- read.csv(chip_muts_to_exclude)

# Purge events that happen in more than 7 people
pair_counts <- combined_df %>%
  group_by(Gene, Protein_annotation) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  mutate(pair_concat = paste(Gene, Protein_annotation, sep = "|"))

# Filter pairs that appear more than 7 times
pairs_to_remove <- pair_counts %>%
  filter(count > 5) %>%
  pull(pair_concat)

# Remove rows with pairs that appear more than 7 times
filtered_df <- combined_df %>%
  filter(!(paste(Gene, Protein_annotation, sep = "|") %in% pairs_to_remove))

write.csv(filtered_df, file = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/chip_WBC_only.csv", row.names = FALSE)







