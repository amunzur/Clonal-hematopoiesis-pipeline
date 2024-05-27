"""
This script explores the false positives we would call as CH if we didn't have the matched WBC.
We use the variant calls after filtering but before combining with WBC.
"""

path_true_chip = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/chip_SSCS2_curated_complete.csv"
curated_chip = pd.read_csv(path_true_chip)

gene_list = ["TP53", "ATM", "BRCA1", "BRCA2", "CHEK2", "PTEN"]
for gene in gene_list: 
    for diagnosis, cohort_size in zip(["Kidney", "Bladder"], [239, 125]):
        df = curated_chip[(curated_chip["Gene"] == gene) & (curated_chip["VAF_t"] > 1) & (curated_chip["Diagnosis"] == diagnosis)]
        n_pts_miscall = len(df["Patient_id"].unique())
        perc_miscall = round((n_pts_miscall/cohort_size)*100, 2)
        print(f"{gene}, {diagnosis}, {n_pts_miscall}, {perc_miscall}")

