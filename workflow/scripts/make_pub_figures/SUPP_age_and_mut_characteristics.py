"""
Makes 3 plots exploring age and :
1. fraction of the group that is CH+ 
2. Median CH vaf in the age group
3. Median number of mutations per patient in the group
"""

import pandas as pd
import numpy as np
import os 
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

path_chip = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/chip.csv"
path_sample_information = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"
path_kidney_clin = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/Supplementary tables - Clinical data - mRCC.tsv"
path_bladder_clin = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/bladder/clinical_data.csv"

# LOAD CHIP DATASETS
all_vars_chip = pd.read_csv(path_chip)
baseline = all_vars_chip[(all_vars_chip["Dependent"] == False) & (all_vars_chip["Timepoint"] == "Baseline")]

kidney_clin = pd.read_csv(path_kidney_clin, sep = "\t")[["GUBB ID", "Age at GUBB draw"]].rename(columns= {"GUBB ID": "Patient_id", "Age at GUBB draw": "Age"})
kidney_clin["Diagnosis"] = "mRCC"
kidney_clin["Color"] = "orangered"

bladder_clin = pd.read_csv(path_bladder_clin)
bladder_clin = bladder_clin[bladder_clin["First sample?"] == True]
bladder_clin = bladder_clin[["Patient_id", "Age at blood draw"]].rename(columns = {"Age at blood draw": "Age"})
bladder_clin["Diagnosis"] = "mUC"
bladder_clin["Color"] = "deepskyblue"

# Prepare age data
age_df = pd.concat([kidney_clin, bladder_clin]).reset_index(drop = True)
bins = [20, 40, 50, 60, 70, 80, 90]
labels = [f'{bins[i]}-{bins[i+1]}' for i in range(len(bins) - 1)]
age_df['Age_Bin'] = pd.cut(age_df['Age'], bins=bins, labels=labels, right=False)
age_df['Age_Bin'] = pd.Categorical(age_df['Age_Bin'], categories=labels, ordered=True)

fig = plt.figure(figsize=(4, 8))
gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 1], hspace = 0.1, wspace = 0.1) # outer most gs with 3 rows
ax0 = plt.subplot(gs[0])
ax1 = plt.subplot(gs[1])
ax2 = plt.subplot(gs[2])

def age_and_n_muts(age_df, baseline, ax):
    """
    Correlation between age and the number of mutations per patient.
    """
    n_muts = baseline["Patient_id"].value_counts().reset_index().rename(columns = {"index": "Patient_id", "Patient_id": "n_muts"})
    merged = n_muts.merge(age_df, how = "inner")
    merged = merged.sort_values(by='Age_Bin') # Sort by Age_Bin to ensure proper order
    merged['Age_Bin_Num'] = merged['Age_Bin'].cat.codes
    
    # Fit and plot second-degree polynomial for mRCC (kidney)
    kidney_data = merged[merged["Diagnosis"] == "mRCC"]
    if not kidney_data.empty:
        kidney_fit = np.polyfit(kidney_data["Age_Bin_Num"], kidney_data["n_muts"], 2)
        kidney_poly = np.poly1d(kidney_fit)
        ax.plot(kidney_data["Age_Bin"], kidney_poly(kidney_data["Age_Bin_Num"]), color="red", label="Kidney (mRCC)")
    
    # Fit and plot second-degree polynomial for mUC (bladder)
    bladder_data = merged[merged["Diagnosis"] == "mUC"]
    if not bladder_data.empty:
        bladder_fit = np.polyfit(bladder_data["Age_Bin_Num"], bladder_data["n_muts"], 2)
        bladder_poly = np.poly1d(bladder_fit)
        ax.plot(bladder_data["Age_Bin"], bladder_poly(bladder_data["Age_Bin_Num"]), color="blue", label="Bladder (mUC)")
    
    ax.scatter(merged["Age_Bin"], merged["n_muts"], color=merged["Color"], s = 2)
    ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
    ax.set_yticks([0, 4, 8, 12, 16, 20])
    ax.set_yticklabels(["0", "4", "8", "12", "16", "20"])
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylabel("Number of mutations\nper patient")
    return(ax)


def age_and_median_vaf(age_df, baseline, ax):
    """
    Correlation between age and the number of mutations per patient.
    """
    median_vaf = baseline.groupby("Patient_id")["VAF_n"].median().reset_index().rename(columns = {"median vaf": "Patient_id", "VAF_n": "median vaf"})
    merged = median_vaf.merge(age_df, how = "inner")
    merged = merged.sort_values(by='Age_Bin') # Sort by Age_Bin to ensure proper order
    merged['Age_Bin_Num'] = merged['Age_Bin'].cat.codes
    
    # Fit and plot second-degree polynomial for mRCC (kidney)
    kidney_data = merged[merged["Diagnosis"] == "mRCC"]
    if not kidney_data.empty:
        kidney_fit = np.polyfit(kidney_data["Age_Bin_Num"], kidney_data["median vaf"], 2)
        kidney_poly = np.poly1d(kidney_fit)
        ax.plot(kidney_data["Age_Bin"], kidney_poly(kidney_data["Age_Bin_Num"]), color="red", label="Kidney (mRCC)")
    
    # Fit and plot second-degree polynomial for mUC (bladder)
    bladder_data = merged[merged["Diagnosis"] == "mUC"]
    if not bladder_data.empty:
        bladder_fit = np.polyfit(bladder_data["Age_Bin_Num"], bladder_data["median vaf"], 2)
        bladder_poly = np.poly1d(bladder_fit)
        ax.plot(bladder_data["Age_Bin"], bladder_poly(bladder_data["Age_Bin_Num"]), color="blue", label="Bladder (mUC)")
    
    ax.scatter(merged["Age_Bin"], merged["median vaf"], color=merged["Color"], s = 2)
    ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
    ax.set_yticks([0, 4, 8, 12, 16, 20])
    ax.set_yticklabels(["0", "4", "8", "12", "16", "20"])
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylabel("Number of mutations\nper patient")
    return(ax)



ax0 = age_and_n_muts(age_df, baseline, ax0)
ax1 = age_and_median_vaf(age_df, baseline, ax1)

gs.tight_layout(fig)
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/age_and_ch_characteristics.png")
    
    

