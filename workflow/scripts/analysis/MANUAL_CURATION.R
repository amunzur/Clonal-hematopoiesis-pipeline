# After generating a list of potential variants, this script helps further narrow them down based on manual curation

library(tidyverse)

excluded_annotations <- c(
    "SETDB1:NM_001243491:exon9:c.1180dupG:p.G393fs", # long repeat
    "BRCA2:NM_000059:exon10:c.1806dupA:p.G602fs", 
    "BRCA1:NM_007297:exon9:c.1820dupA:p.K607fs,BRCA1:NM_007294:exon10:c.1961dupA:p.K654fs,BRCA1:NM_007300:exon10:c.1961dupA:p.K654fs", 
    "DNMT3A:NM_001320893:exon14:c.1799_1801del:p.600_601del,DNMT3A:NM_153759:exon15:c.1688_1690del:p.563_564del,DNMT3A:NM_022552:exon19:c.2255_2257del:p.752_753del,DNMT3A:NM_175629:exon19:c.2255_2257del:p.752_753del", 
    "ASXL1:NM_015338:exon12:c.1927dupG:p.G642fs", 
    "BRCA2:NM_000059:exon23:c.9090dupA:p.T3030fs", 
    "BRCA1:NM_007297:exon9:c.1820dupA:p.K607fs,BRCA1:NM_007294:exon10:c.1961dupA:p.K654fs,BRCA1:NM_007300:exon10:c.1961dupA:p.K654fs", 
    "ASXL1:NM_015338:exon12:c.1927dupG:p.G642fs", 
    "BRCA2:NM_000059:exon10:c.1806delA:p.G602fs", 
    "CUX1:NM_001202543:exon18:c.2373delC:p.N791fs,CUX1:NM_181552:exon18:c.2340delC:p.N780fs", 
    "TET2:NM_001127208:exon3:c.1961delA:p.Q654fs,TET2:NM_017628:exon3:c.1961delA:p.Q654fs", 
    "IDH1:NM_001282386:exon6:c.638delA:p.N213fs,IDH1:NM_001282387:exon6:c.638delA:p.N213fs,IDH1:NM_005896:exon6:c.638delA:p.N213fs", 
    "DNMT3A:NM_001320893:exon10:c.1255delG:p.A419fs,DNMT3A:NM_153759:exon11:c.1144delG:p.A382fs,DNMT3A:NM_022552:exon15:c.1711delG:p.A571fs,DNMT3A:NM_175629:exon15:c.1711delG:p.A571fs", 
    "PPM1D:NM_003620:exon1:c.298delT:p.F100fs")

df_path <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/combined/tumor_wbc.csv"
output_path <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/combined/curated.csv"
df_main <- read_csv(df_path)

df <- df_main %>%
    filter(
        Function == "exonic", 
        !(AAchange %in% excluded_annotations))

write_csv(df, output_path)


