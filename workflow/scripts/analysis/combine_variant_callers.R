library(tidyverse)
library(stringr)

cohort_name <- "batch5"
varscan_PATH <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/VarScan2/finalized", cohort_name, "chip_variants.csv")
vardict_PATH <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/Vardict/finalized", cohort_name, "chip_variants.csv")

varscan <- read_csv(varscan_PATH)
vardict <- read_csv(vardict_PATH)

varscan$variant_caller <- "Varscan"
vardict$variant_caller <- "Vardict"

# variants identified by both tools 
vardict$Position <- as.character(vardict$Position)
varscan$Position <- as.character(varscan$Position)
both <- inner_join(varscan, vardict, by = c("Sample_name", "Chrom", "Position"))
# both$variant_caller <- "Both"

both_varscan <- both[, -grep("\\.y", names(both))]
names(both_varscan) <- names(varscan)
varscan_only <- anti_join(varscan, vardict, by = c("Sample_name", "Chrom", "Position"))
varscan_only$variant_caller <- "Varscan"
both_varscan$variant_caller <- "Both"

both_vardict <- both[, -grep("\\.x", names(both))]
names(both_vardict) <- gsub("\\.y", "", names(both_vardict))
both_vardict <- select(both_vardict, names(both_varscan))

# if dealing with batch1, consider a few more steps: 
if (cohort_name == "new_chip_panel") {

	both_vardict <- both_vardict %>% relocate(Sample_name_finland.y, .after = Sample_name)
	both_vardict <- both_vardict %>% relocate(Patient_ID.y, .after = Sample_name_finland.y)

} else {

	# both_vardict <- both_vardict %>% relocate(Patient_ID.y, .after = Sample_name)

}

names(both_vardict) <- names(vardict)
vardict_only <- anti_join(vardict, varscan, by = c("Sample_name", "Chrom", "Position"))
vardict_only$variant_caller <- "Vardict"
both_vardict$variant_caller <- "Both"

varscan <- rbind(both_varscan, varscan_only)
vardict <- rbind(both_vardict, vardict_only)
combined <- rbind(both_varscan, varscan_only, vardict_only)

varscan_tosave <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/combined", cohort_name, "varscan.csv")
vardict_tosave <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/combined", cohort_name, "vardict.csv")
combined_tosave <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/combined", cohort_name, "combined.csv")

dir.create(dirname(varscan_tosave))
dir.create(dirname(vardict_tosave))
dir.create(dirname(combined_tosave))

write_csv(varscan, varscan_tosave)
write_csv(vardict, vardict_tosave)
write_csv(combined, combined_tosave)