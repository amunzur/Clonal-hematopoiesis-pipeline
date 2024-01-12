library(tidyverse)
library(gsubfn)

PATH_mutato <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/Mutato/finalized/chip_bg20.vcf"
PATH_final <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/Mutato/finalized/chip_bg20_processed.csv"    

df <- as.data.frame(read.delim(PATH_mutato)) %>%
    rename(Chrom = CHROM,
           Position = POSITION,
           Ref = REF, 
           Alt = ALT,
           error_rate = NOTES, 
           Gene = GENE) %>%
    mutate(across(starts_with(c("GUBB", "VIP")), as.character)) %>%
    gather("sample_names", "variants", starts_with(c("GUBB", "VIP"))) %>%
    mutate(sample_type = str_extract(sample_names, "cfDNA|gDNA"), 
           sample_names = gsub("[.]", "-", sample_names), 
           patient_id = unlist(lapply(str_split(sample_names, "_gDNA|_WBC|_cfDNA"), "[", 1))) %>%
    relocate(patient_id, .before = EFFECT) %>%
    filter(grepl("[*]", variants)) %>%
    separate(variants, c("alt_reads", "total_reads", "mapq", "sidedness", "perc_strand", "perc_other_strand"), remove = TRUE) %>%
    mutate(VAF = as.numeric(alt_reads)/as.numeric(total_reads))

cfdna <- filter(df, sample_type == "cfDNA")
gdna <- filter(df, sample_type == "gDNA")

combined <- inner_join(gdna, cfdna, suffix = c("_n", "_t"), by = c("Chrom", "Position", "Ref", "Alt", "Gene", "EFFECT", "patient_id", "error_rate"))

write_csv(combined, PATH_final)
