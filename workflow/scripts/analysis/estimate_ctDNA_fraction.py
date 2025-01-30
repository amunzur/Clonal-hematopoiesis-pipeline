"""
Given a list of somatic mutations calculate the ctDNA fraction of the samples.
"""
# Estimate cfDNA using mutations
import pandas as pd
import numpy as np
import os
from scipy.stats import binom

def get_adj_vaf(vaf, depth):
    """
    Adjust the mutation VAF using binomial distribution.
    """
    alt = vaf*depth
    for p in np.arange(0,1,0.005):
        dist = binom(depth, p)
        if dist.cdf(alt) < 0.95:
            vaf = p; break;
    return(vaf)

DIR_working = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling"
PATH_mutations_curated = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/somatic.csv"
PATH_sample_dict = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"
PATH_kidney_samples = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_list_kidney.tsv"
PATH_bladder_samples = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_list_bladder.tsv"

PATH_mutation_ctfractions = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/ctDNA_fractions.csv" # output path
PATH_bladder_ctDNA = "/groups/wyattgrp/users/amunzur/pipeline/resources/ct_fractions/ERBB2 Supplemental Tables - Table S1 - ctDNA cohort.csv"


DIR_cnvkit = "/groups/wyattgrp/users/amunzur/pipeline/results/cnvkit/final_result"

sample_dict = pd.read_csv(PATH_sample_dict, "\t", names =["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
muts = pd.read_csv(PATH_mutations_curated)
muts = muts[muts["Dependent"] == False]

# Gnenerate the median log2 ratio per gene per sample
mylist = []
for file in os.listdir(DIR_cnvkit):
    sample_name = file.replace(".tsv", "")
    print(sample_name)
    median_log2 = pd.read_csv(os.path.join(DIR_cnvkit, file), sep = "\t")[["gene", "log2"]]
    median_log2["median"] = median_log2.groupby("gene").transform("median")
    median_log2 = median_log2[["gene", "median"]].drop_duplicates().reset_index(drop = True)
    median_log2["Sample_name_t"] = sample_name
    mylist.append(median_log2)

gene_medians = pd.concat(mylist).reset_index(drop = True).rename(columns = {"gene": "Gene", "median": "Median log2"})
muts = muts.merge(gene_medians, how = "left")
muts["VAF_t"] = muts["VAF_t"]/100

muts = muts[(muts['Median log2'] > -0.3) & (muts['Median log2'] < 0.3) &
            (muts['Depth_t'] > 40) &
            ((muts['Diagnosis'] == "Kidney") | (muts['Diagnosis'] == "Bladder")) &
            (~muts['Chrom'].str.contains("X", na = False)) &
            (~muts['Chrom'].str.contains("Y", na = False)) & 
            (~muts["Patient_id"].isin(['20-313', '21-430', '21-184']))]

muts['Adj_allele_frequency'] = muts.apply(lambda row: get_adj_vaf(row['VAF_t'], row['Depth_t']), axis = 1)
muts["Max_Sample_adj_VAF"] = muts.groupby('Sample_name_t')["Adj_allele_frequency"].transform('max')
muts = muts[muts["Adj_allele_frequency"] == muts["Max_Sample_adj_VAF"]]
muts['row_max'] = muts[['Depth_t']].max(axis=1) #To break any ties for max VAF!
muts = muts.loc[muts.groupby('Sample_name_t')['row_max'].idxmax()]

muts['Mutation_ctDNA_fraction'] = 2/(1+1/(muts['Max_Sample_adj_VAF']))
muts = muts.merge(sample_dict, how = "outer") # some pts have no mutations available to call CT fractions, this retains them and allows a helpful comparison between CN based and mut based estimates across the whole cohort.
muts = muts[["Patient_id", "Diagnosis", "Sample_name_t", "Date_collected", "VAF_t", "Depth_t", "Gene", "Consequence", "Mutation_ctDNA_fraction"]]
muts = muts[~muts["Patient_id"].str.contains("Normal")] # remove VIP normals
muts["Mutation_ctDNA_fraction"] = muts["Mutation_ctDNA_fraction"].fillna(0)
muts["ctDNA_status"] = muts["Mutation_ctDNA_fraction"].apply(lambda x: "ctDNA negative" if x == 0 else "ctDNA positive")
muts['Date_collected'] = pd.to_datetime(muts['Date_collected'], format='%Y%b%d')

# add the sample name to the ctneg samples
kidney_samps = pd.read_csv(PATH_kidney_samples, sep = "\t").rename(columns = {"sample_names": "Sample_name_t"})
bladder_samps = pd.read_csv(PATH_bladder_samples, sep = "\t").rename(columns = {"sample_names": "Sample_name_t"})
all_samps = pd.concat([kidney_samps, bladder_samps]).reset_index(drop = True)

# all_samps = pd.read_csv(PATH_all_samples, sep = "\t").rename(columns = {"sample_names": "Sample_name_t"})
all_samps["Patient_id"] = all_samps["Sample_name_t"].str.split("_").str[0]
all_samps["Patient_id"] = all_samps["Patient_id"].str.replace("GU-", "")
all_samps["Date_collected"] = all_samps["Sample_name_t"].str.split("[-_]").str[-1]
cfDNA_samps = all_samps[all_samps["Sample_name_t"].str.contains("cfDNA")]
cfDNA_samps["Date_collected"] = pd.to_datetime(cfDNA_samps['Date_collected'], format='%Y%b%d')
del muts["Sample_name_t"]
muts = muts.merge(cfDNA_samps, how = "left")

# Dropping these two, since they had oxidative damage 
muts = muts[muts["Patient_id"] != "21-430"]
muts.to_csv("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/ctDNA_fractions.csv", index = False)


# Add one sample manually. This library was dropped from the CH cohort because it had oxidative damage. 
# Mutation ctDNA estimation is from the bladder cohort.
# muts = muts[muts["Patient_id"] != "21-430"]
# dict_sample_1 = {
#     "Patient_id": "21-430",
#     "Diagnosis": "Bladder",
#     "Sample_name_t": "GU-21-430_cfDNA-Baseline-IDT-2021Nov10"
#     "Date_collected": "2021-11-10",
#     "VAF_t": 0.15,
#     "Depth_t": 2234,
#     "Gene": "PIK3CA",
#     "Consequence": "Missense",
#     "Mutation_ctDNA_fraction": 0.26,
#     "ctDNA_status": "ctDNA positive",
# }

# dict_sample_2 = {
#     "Patient_id": "21-430",
#     "Diagnosis": "Bladder",
#     "Sample_name_t": "GU-21-430_cfDNA-OnTreatment-IDT-2023Jan06"
#     "Date_collected": "2023-01-06",
#     "VAF_t": 0.11,
#     "Depth_t": 3128,
#     "Gene": "PIK3CA",
#     "Consequence": "Missense",
#     "Mutation_ctDNA_fraction": 0.24,
#     "ctDNA_status": "ctDNA positive",
# }

# dict_sample_3 = {
#     "Patient_id": "21-184",
#     "Diagnosis": "Bladder",
#     "Sample_name_t": "GU-21-184_cfDNA-Baseline-IDT-2021Mar12"
#     "Date_collected": "2021-03-12",
#     "VAF_t": 0.11,
#     "Depth_t": 3128,
#     "Gene": "FGFR3",
#     "Consequence": "Missense",
#     "Mutation_ctDNA_fraction": 0.24,
#     "ctDNA_status": "ctDNA positive",
# }

# .bam

# Add the ctDNA fraction estimations from the bladder project to compare
bladder_ctDNA = (
    pd.read_csv(PATH_bladder_ctDNA)[["Patient ID", "Collection date", "ctDNA fraction", "Targeted DNA sequencing panel"]]
    .rename(columns={"Patient ID": "Patient_id", "Collection date": "Date_collected"})
    .assign(Date_collected=lambda x: pd.to_datetime(x['Date_collected'], format='%d-%b-%Y'))
    .assign(Date_collected=lambda x: x['Date_collected'].dt.strftime('%Y%b%d'))
    .query('Patient_id.str.startswith("GU-")', engine='python')  # Filter rows where Patient_id starts with 'GU-'
    .assign(Patient_id=lambda x: x['Patient_id'].str.replace('GU-', ''))
)

bladder_ctDNA['Date_collected'] = pd.to_datetime(bladder_ctDNA['Date_collected'], format='%Y%b%d')

# All patients that have a ctDNA estimate done in the bladder cohort will have that estimate.
sample_dict["Date_collected"] = pd.to_datetime(sample_dict["Date_collected"], format='%Y%b%d')
bladder_panel_estimates = bladder_ctDNA.merge(sample_dict, on = ["Patient_id", "Date_collected"])
bladder_panel_estimates = bladder_panel_estimates[["Patient_id", "Date_collected", "Diagnosis", "Timepoint", "ctDNA fraction", "Targeted DNA sequencing panel"]]

# Patients for which we will use the CHIP panel estimate
chip_panel_estimates = muts[muts["Diagnosis"] == "Bladder"].merge(bladder_ctDNA, indicator = True, on = ["Patient_id", "Date_collected"], how = "outer")
chip_panel_estimates = chip_panel_estimates[chip_panel_estimates["_merge"] == "left_only"].reset_index(drop = True)
chip_panel_estimates = chip_panel_estimates[["Patient_id", "Diagnosis", "Date_collected", "Mutation_ctDNA_fraction"]]
chip_panel_estimates["Targeted DNA sequencing panel"] = "CHIP"
chip_panel_estimates = chip_panel_estimates.rename(columns = {"Mutation_ctDNA_fraction": "ctDNA fraction"})

kidney_ctdna_fr = muts[muts["Diagnosis"] == "Kidney"]
kidney_ctdna_fr = kidney_ctdna_fr.merge(sample_dict, how = "inner")
kidney_ctdna_fr = kidney_ctdna_fr[["Patient_id", "Date_collected", "Diagnosis", "Timepoint", "Mutation_ctDNA_fraction"]].rename(columns = {"Mutation_ctDNA_fraction": "ctDNA fraction"})
kidney_ctdna_fr["Targeted DNA sequencing panel"] = "CHIP"

# Final file on ctDNa fractions
combined = pd.concat([bladder_panel_estimates, chip_panel_estimates, kidney_ctdna_fr]).reset_index(drop = True)
combined.to_csv("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/ctDNA_fractions.csv", index = False)



# combined = muts[["Patient_id", "Date_collected", "Diagnosis", "Mutation_ctDNA_fraction"]].merge(bladder_ctDNA[["Patient_id", "Date_collected", "ctDNA_fraction_bladder_project"]], how = "inner", on = ["Patient_id", "Date_collected"])
# combined.to_csv("/groups/wyattgrp/users/amunzur/pipeline/results/misc/ERBB2_ctDNA_frac_difference.csv", index = False)




xx = muts.merge(bladder_ctDNA, on = ["Patient_id", "Date_collected"], how = "inner")
xx = xx[["Patient_id", "Date_collected", "Mutation_ctDNA_fraction", "ctDNA fraction"]]
xx = xx.rename(columns = {"Mutation_ctDNA_fraction": "CHIP panel", "ctDNA fraction": "Bladder panel"})
xx.to_csv("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/bladder_vs_chip_panel_ct_estimates.csv",index = False)