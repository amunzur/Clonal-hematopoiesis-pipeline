import pandas as pd
import numpy as np
import os 
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from scipy.stats import fisher_exact
from scipy.stats import fisher_exact
from scipy.stats import mannwhitneyu
from scipy.stats import wilcoxon
from scipy.stats import ttest_rel

mpl.rcParams['font.size'] = 10
mpl.rcParams['text.color'] = 'k'
mpl.rcParams['legend.fontsize'] = 10
mpl.rcParams['legend.handletextpad'] = '0.8'
mpl.rcParams['legend.labelspacing'] = '0.4'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10
mpl.rcParams['axes.labelsize'] = 10

def get_baseline_and_ot_vafs(muts_df, gene, pts_with_ot, return_difference = False): 
    """
    Given a pt id and a gene name, follows all the mutations in that gene and returns their vaf in all samples of that patient.
    If return_difference = True, then returns the VAF difference in the two timepoints.
    """
    # Determine if this is a CHIP df or a ctDNA mutation df
    if "Status" in muts_df.columns:
        col = "VAF_t"
    else:
        col = "VAF_n"
    
    muts_df = muts_df[(muts_df["Patient_id"].isin(pts_with_ot)) & (muts_df["Gene"] == gene)]
    results_list = []
    for pt in pts_with_ot:
        pt_df = muts_df[muts_df["Patient_id"] == pt]
        # Go through each mutation in the gene
        for i, group in pt_df.groupby("Protein_annotation"): 
            mut_name = group["Protein_annotation"].unique()[0]
            group["Date_collected"] = pd.to_datetime(group["Date_collected"], format = '%Y%b%d')
            group = group.sort_values(by = "Date_collected")
            vafs_list = group[col].tolist()[0:2] # only showing the extra timepoint after the baseline
            if return_difference:
                vafs_list = calculate_distances(vafs_list)
            results_list.append(vafs_list)
    return(results_list)

def calculate_distances(lst):
    """
    computes the distances (differences) between consecutive elements in a given list.
    """
    distances = [lst[i+1] - lst[i] for i in range(len(lst) - 1)]
    return distances

def flatten(nested_list):
    flat_list = []
    for item in nested_list:
        if isinstance(item, list):
            flat_list.extend(flatten(item))  # Recursively flatten if the item is a list
        else:
            flat_list.append(item)
    return flat_list

PATH_sample_information = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"
path_chip = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/chip.csv"
path_somatic = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/somatic.csv"

sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
pts_with_ot = sample_info[sample_info["Timepoint"] == "During treatment"]["Patient_id"].unique().tolist()

chip = pd.read_csv(path_chip)
somatic = pd.read_csv(path_somatic)

###############################################################
# Mean change from baseline to OT
genes_list = chip["Gene"].unique()
baseline_vafs_list = []
ot_vafs_list = []

for gene in genes_list:
    print(gene)
    mutations = get_baseline_and_ot_vafs(chip, gene, pts_with_ot, return_difference = False) # get the vafs all mutations in a given gene from baseline to timepoint2
    baseline_vafs = flatten([i[0] for i in mutations])
    ot_vafs = flatten([i[1] for i in mutations])
    
    # add to list
    baseline_vafs_list.extend(baseline_vafs)
    ot_vafs_list.extend(ot_vafs)

np.mean(baseline_vafs_list)
np.mean(ot_vafs_list)
np.mean(ot_vafs_list) - np.mean(baseline_vafs_list)

stat, p_value = wilcoxon(baseline_vafs_list, ot_vafs_list)

###############################################################
# CHIP
genes_list = ["DNMT3A", "TET2", "ASXL1", "TP53", "PPM1D", "CHEK2", "ATM", "BRCA1", "BRCA2"]
dta_genes_muts_chip = []
ddr_genes_muts_chip = []

for gene in genes_list:
    print(gene)
    mutations = get_baseline_and_ot_vafs(chip, gene, pts_with_ot, return_difference = False) # get the vafs all mutations in a given gene from baseline to timepoint2
    mutations_vaf_diff = get_baseline_and_ot_vafs(chip, gene, pts_with_ot, return_difference = True) # get the vaf differences from baseline to the next OT
    mutations_vaf_diff = flatten(mutations_vaf_diff)
    
    # Append to list for doing stats later
    if gene in ["DNMT3A", "TET2", "ASXL1"]:
        dta_genes_muts_chip.extend(mutations_vaf_diff)
    elif gene in ["TP53", "PPM1D", "CHEK2"]:
        ddr_genes_muts_chip.extend(mutations_vaf_diff)


np.var(dta_genes_muts_chip)
np.var(ddr_genes_muts_chip)

np.mean(dta_genes_muts_chip + ddr_genes_muts_chip)

np.median(dta_genes_muts_chip)
np.median(ddr_genes_muts_chip)
mannwhitneyu(dta_genes_muts_chip, ddr_genes_muts_chip)

# Genes most commonly mutated in OT samples
chip_ot = chip[(chip["Timepoint"] == "During treatment") & (chip["Dependent"] == False)]
n_patient_OT = chip_ot["Patient_id"].unique().shape[0] # number of that are CH+ at any OT sample

chip_ot_genes = chip_ot[["Patient_id", "Gene"]].drop_duplicates().sort_values(by = "Patient_id").reset_index(drop = True)
chip_ot_genes = chip_ot_genes["Gene"].value_counts().reset_index().rename(columns = {"index": "Gene", "Gene": "Counts"}) 
chip_ot_genes["Patient prevalence"] = round(chip_ot_genes["Counts"]/n_patient_OT*100, 2)

# Did PPM1D mutation all increase in VAF in patients whose OT samples was taken after platinum chemo?
path_clinical = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/bladder/clinical_data.csv"
clin = pd.read_csv(path_clinical)

clin_platinum = clin[(clin["First sample?"] == False) & (clin["Treatment"].isin(["Platinum chemotherapy", "Combination"]))]
platinum_samples = clin_platinum[["Patient_id", "Date_collected"]] # OT samples taken after platinum chemo
platinum_samples["Date_collected"] = pd.to_datetime(platinum_samples["Date_collected"])

chip["Date_collected"] = pd.to_datetime(chip["Date_collected"], format = '%Y%b%d')
platinum_samples = chip[["Patient_id", "Date_collected", "Timepoint", "VAF_n", "Gene"]].merge(platinum_samples).sort_values(by = "Gene")

# Checking if the pts with OT samples also had DNMT3A mutations at baseline
chip[(chip["Patient_id"].isin(pts_with_ot)) & (chip["Timepoint"] == "Baseline")][["Patient_id", "Gene"]].drop_duplicates()["Gene"].value_counts()
