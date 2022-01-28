# This function compares the variants called in tumor and normal samples. 
# Before running this script, first all the variant calling results must be compiled in one csv file for both tumor and normal sample.

library(tidyverse)

combine_tumor_wbc <- function(variant_caller, DIR_output, DIR_vars){

    file_paths <- as.list(grep(variant_caller, list.files(DIR_vars, recursive = TRUE, full.names = TRUE), value = TRUE))[1:4]

    vars <- as.data.frame(do.call(rbind, lapply(file_paths, read_csv)))
    tumor <- vars %>%
            filter(Sample_type == "Tumor") %>%
            select(Patient_ID, Chrom, Position, Ref, Alt, Function, Gene, AAchange, Protein_annotation, Effects, ExAC_ALL, Variant, Sample_type, Sample_name, Cohort_name, VAF, Ref_reads, Alt_reads, Duplicate, Depth, variant_caller)
    wbc <- vars %>% 
            filter(Sample_type == "WBC") %>%
            select(Patient_ID, Chrom, Position, Ref, Alt, Function, Gene, AAchange, Protein_annotation, Effects, ExAC_ALL, Variant, Sample_type, Sample_name, Cohort_name, VAF, Ref_reads, Alt_reads, Duplicate, Depth, variant_caller)

    combined <- inner_join(tumor, wbc, by = c("Patient_ID", "Chrom", "Position", "Ref", "Alt", "Function", "Gene", "AAchange", "Protein_annotation", "Effects", "ExAC_ALL", "Variant"))
    combined <- inner_join(tumor, wbc, by = c("Patient_ID", "Chrom", "Position", "Ref", "Alt"))
    names(combined) <- gsub("\\.x", "_t", names(combined)) 
    names(combined) <- gsub("\\.y", "_n", names(combined)) 

    combined <- combined %>%
            mutate(tumor_wbc_vaf_ratio = round((VAF_t / VAF_n), 2), 
                    tumor_wbc_depth_ratio = round((Depth_t / Depth_n), 2))

    write_csv(combined, file.path(DIR_output, paste(variant_caller, "tumor_wbc.csv", sep = "_")))
    write_delim(combined, file.path(DIR_output, paste(variant_caller, "tumor_wbc.tsv", sep = "_")))

}

DIR_output <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/VarScan2/tumor_wbc"
DIR_vars <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/combined")
combine_tumor_wbc("varscan", DIR_output, DIR_vars)

DIR_output <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/Vardict/tumor_wbc"
combine_tumor_wbc("vardict", DIR_output, DIR_vars)

####################################################################################
# The function below works with the most recent combined file that has info about all variant callers in one file
combine_tumor_wbc <- function(DIR_vars){

    file_paths <- as.list(grep("combined.csv", list.files(DIR_vars, recursive = TRUE, full.names = TRUE), value = TRUE))[1:4]

    vars <- as.data.frame(do.call(rbind, lapply(file_paths, read_csv)))
    tumor <- vars %>%
            filter(Sample_type == "Tumor") %>%
            select(Patient_ID, Chrom, Position, Ref, Alt, Function, Gene, AAchange, Protein_annotation, Effects, ExAC_ALL, Variant, Sample_type, Sample_name, Cohort_name, VAF, Ref_reads, Alt_reads, Duplicate, Depth, varscan, vardict, gatk, n_callers)
    wbc <- vars %>% 
            filter(Sample_type == "WBC") %>%
            select(Patient_ID, Chrom, Position, Ref, Alt, Function, Gene, AAchange, Protein_annotation, Effects, ExAC_ALL, Variant, Sample_type, Sample_name, Cohort_name, VAF, Ref_reads, Alt_reads, Duplicate, Depth, varscan, vardict, gatk, n_callers)

    combined <- inner_join(tumor, wbc, by = c("Patient_ID", "Chrom", "Position", "Ref", "Alt", "Function", "Gene", "AAchange", "Protein_annotation", "Effects", "ExAC_ALL", "Variant"))
    names(combined) <- gsub("\\.x", "_t", names(combined)) 
    names(combined) <- gsub("\\.y", "_n", names(combined)) 

    combined <- combined %>%
            mutate(tumor_wbc_vaf_ratio = round((VAF_t / VAF_n), 2), 
                    tumor_wbc_depth_ratio = round((Depth_t / Depth_n), 2))

    return(combined)
}

# adds the abs path to the bam files to make it easier to run IGV snapshots
add_bam_path <- function(variant_df, DIR_bams){

        PATH_tumor_bam <- file.path(DIR_bams, variant_df$Cohort_name_t, "SC_penalty", paste0(variant_df$Sample_name_t, ".bam"))
        PATH_wbc_bam <- file.path(DIR_bams, variant_df$Cohort_name_n, "SC_penalty", paste0(variant_df$Sample_name_n, ".bam"))

        variant_df <- variant_df %>% 
                                mutate(
                                    Path_bam_t = PATH_tumor_bam,
                                    Path_bam_n = PATH_wbc_bam) %>%
                                relocate(
                                    Path_bam_t, .after = Sample_name_t) %>%
								relocate(
									Path_bam_n, .after = Sample_name_n)
        return(variant_df)
}

DIR_vars <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/combined"
DIR_bams <- "/groups/wyattgrp/users/amunzur/pipeline/results/data/bam"
combined <- combine_tumor_wbc(DIR_vars)
combined <- add_bam_path(combined, DIR_bams)

DIR_output <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/combined"
write_csv(combined, file.path(DIR_output, "tumor_wbc.csv"))
write_delim(combined, file.path(DIR_output, "tumor_wbc.tsv"), delim = "\t")
