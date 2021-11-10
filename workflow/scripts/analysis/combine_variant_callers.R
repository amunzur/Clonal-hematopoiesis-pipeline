library(tidyverse)
library(stringr)

varscan_PATH <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/VarScan2/finalized/chip_variants.csv"
vardict_PATH <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/Vardict/finalized/chip_variants.csv"

varscan <- read_csv(varscan_PATH)
vardict <- read_csv(vardict_PATH)

varscan$variant_caller <- "Varscan"
vardict$variant_caller <- "Vardict"

# variants identified by both tools 
both <- inner_join(varscan, vardict, by = c("Sample_name", "Chrom", "Position"))
# both$variant_caller <- "Both"

both_varscan <- both[, -grep("\\.y", names(both))]
names(both_varscan) <- names(varscan)
varscan_only <- anti_join(varscan, vardict, by = c("Sample_name", "Chrom", "Position"))
varscan_only$variant_caller <- "Varscan"
both_varscan$variant_caller <- "Both"

both_vardict <- both[, -grep("\\.x", names(both))]
both_vardict %>% relocate(Patient_ID.y, .after = Sample_name)
names(both_vardict) <- names(vardict)
vardict_only <- anti_join(vardict, varscan, by = c("Sample_name", "Chrom", "Position"))
vardict_only$variant_caller <- "Vardict"
both_vardict$variant_caller <- "Both"

varscan <- rbind(both_varscan, varscan_only)
vardict <- rbind(both_vardict, vardict_only)

# find duplicated rows
# idx <- duplicated(combined[, c("Patient_ID", "Chrom", "Position")])

varscan_tosave <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/combined/varscan.csv"
vardict_tosave <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/combined/vardict.csv"

write_csv(varscan, varscan_tosave)
write_csv(vardict, vardict_tosave)

# dedup <- combined %>% distinct(Patient_ID, Chrom, Position, .keep_all = FALSE)
# write_csv(dedup, "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/VarScan2/finalized/combined.csv")