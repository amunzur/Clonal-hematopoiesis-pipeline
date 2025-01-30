"""
From the final list, divides the VAFs of mutations in X chromosome genes by 2, in men only.
"""

import pandas as pd
import numpy as np
import os
import re

path_chip_muts = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/chip_unadjusted_vaf.csv"
path_ctdna_muts = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/somatic_unadjusted_vaf.csv"

path_chip_muts_adjusted = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/chip.csv"
path_ctdna_muts_adjusted = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/somatic.csv"

path_bladder_clinical = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/bladder/clinical_data.csv"
path_kidney_clinical = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/Supplementary tables - Clinical data - mRCC.tsv"

ch_muts = pd.read_csv(path_chip_muts)
ctdna_muts = pd.read_csv(path_ctdna_muts)
bladder_clin = pd.read_csv(path_bladder_clinical)
bladder_clin = bladder_clin[bladder_clin["First sample?"]].reset_index(drop = True)
kidney_clin = pd.read_csv(path_kidney_clinical, sep = "\t").rename(columns = {"GUBB ID": "Patient_id"})

# Make sex df
sex_df = pd.concat([bladder_clin[["Patient_id", "Sex"]], kidney_clin[["Patient_id", "Sex"]]]).reset_index(drop = True)
additional_pts = {
    "17-051": "Female",
    "17-114": "Male",
    "17-118": "Male",
    "18-129": "Male",
    "20-300": "Male",
    "21-390": "Male",
    "21-391": "Male",
    "21-420": "Male",
}

additional_pts_df = pd.DataFrame(list(additional_pts.items()), columns=["Patient_id", "Sex"])

# Concatenate the additional patients DataFrame with the existing sex_df DataFrame
sex_df = pd.concat([sex_df, additional_pts_df]).reset_index(drop=True)

# Genes on the X chrom
genes_list = ["BCORL1", "BCOR", "BRCC3", "ZRSR2", "AR", "KDM6A", "STAG2"]

def adjust_vafs(muts_df, sex_df, genes_list):
    """
    Given a muts_df, divides the VAF in X chromosome genes in men.
    """
    muts_df = muts_df.merge(sex_df, how = "left")
    
    # Define the condition for adjusting VAF
    condition = (muts_df["Sex"] == "Male") & (muts_df["Gene"].isin(genes_list))
    muts_df.loc[condition, "VAF_t"] /= 2
    muts_df.loc[condition, "VAF_n"] /= 2
    
    del muts_df["Sex"]
    
    return(muts_df)

chip_adjusted = adjust_vafs(ch_muts, sex_df, genes_list)
ctdna_adjusted = adjust_vafs(ctdna_muts, sex_df, genes_list)

chip_adjusted.to_csv(path_chip_muts_adjusted, index = False)
ctdna_adjusted.to_csv(path_ctdna_muts_adjusted, index = False)