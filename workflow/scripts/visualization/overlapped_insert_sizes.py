"""
Makes insert size diagrams
"""
import pandas as pd
import numpy as np
import os
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm

def convert_to_fraction(path_to_insert_size):
    """
    Given an insert size df with absolute values convert it to a df with fractions.
    """
    df = pd.read_csv(path_to_insert_size, sep = "\t", skiprows = 10)
    total_count = df["All_Reads.fr_count"].sum()
    df["Fraction"] = df["All_Reads.fr_count"]/total_count
    
    return(df)

def plot_overlapping_insert_size_histograms(patient_names, insert_size_paths, ax, alpha, color):
    """
    To a given ax plots the insert size histogram for one 
    """
    for pt in patient_names:
        print(f"Working on {pt}")
        
        path = [f for f in insert_size_paths if pt in f and "cfDNA" in f][0]
        df_fractions = convert_to_fraction(path)
        ax.bar(df_fractions["insert_size"], df_fractions["Fraction"], color=color, alpha = alpha, width = 0.1, edgecolor = "None")
        
        ax.set_xlim(0, 1000)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_ylabel("Fraction of insert sizes")
        ax.set_xlabel("Insert size")    
    return(ax)

# SOURCE FUNCTIONS
source_functions = os.path.join("/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/visualization/UTILITIES_make_chip_plots.py")
with open(source_functions, 'r') as file:
    script_code = file.read()

exec(script_code)

DIR_insert_sizes = "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/insert_size/SSCS2"
path_sample_information = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"

# LOAD DATASETS
all_vars_somatic = pd.read_csv(os.path.join("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/somatic.csv"))
base_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Dependent"] == False)].reset_index(drop = True)
ctdna_status = annotate_mutation_status(base_somatic, "Both", path_sample_information, annotate_what = "ctDNA")

sample_info = pd.read_csv(path_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
insert_size_paths = [os.path.join(DIR_insert_sizes, f) for f in os.listdir(DIR_insert_sizes)]

# Determine patient groups
ctdna_pos_kidney_patients = ctdna_status[(ctdna_status["Diagnosis"] == "Kidney") & (ctdna_status["ctDNA status"] == "Positive")]["Patient_id"].unique()
ctdna_neg_kidney_patients = ctdna_status[(ctdna_status["Diagnosis"] == "Kidney") & (ctdna_status["ctDNA status"] == "Negative")]["Patient_id"].unique()
ctdna_pos_bladder_patients = ctdna_status[(ctdna_status["Diagnosis"] == "Bladder") & (ctdna_status["ctDNA status"] == "Positive")]["Patient_id"].unique()
ctdna_neg_bladder_patients = ctdna_status[(ctdna_status["Diagnosis"] == "Bladder") & (ctdna_status["ctDNA status"] == "Negative")]["Patient_id"].unique()

# PLOTTING
# 1. CTDNA POSITIVE SAMPLES
# Makes overlapping plots for kidney and bladder. 
fig_ctpos = plt.figure(figsize=(8, 4))
outer_gs_ctpos = gridspec.GridSpec(1, 2, width_ratios=[1, 1], hspace = 0.3, wspace = 0.3) # outer most gs with 3 rows

ax_kidney_ct_pos = plt.subplot(outer_gs_ctpos[0])
ax_bladder_ct_pos = plt.subplot(outer_gs_ctpos[1], sharey = ax_kidney_ct_pos, sharex = ax_kidney_ct_pos)

ax_kidney_ct_pos = plot_overlapping_insert_size_histograms(ctdna_pos_kidney_patients, insert_size_paths, ax = ax_kidney_ct_pos, alpha = 0.2, color = "orangered")
ax_bladder_ct_pos = plot_overlapping_insert_size_histograms(ctdna_pos_bladder_patients, insert_size_paths, ax = ax_bladder_ct_pos, alpha = 0.2, color = "deepskyblue")

ax_kidney_ct_pos.set_title("ctDNA+ kidney")
ax_bladder_ct_pos.set_title("ctDNA+ bladder")
fig_ctpos.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/insert_size/insert_size_ctDNA_POS.pdf")

# 2. CTDNA NEGATIVE SAMPLES
# Makes overlapping plots for kidney and bladder. 
fig_ctneg = plt.figure(figsize=(8, 4))
outer_gs_ctneg = gridspec.GridSpec(1, 2, width_ratios=[1, 1], hspace = 0.3, wspace = 0.3) # outer most gs with 3 rows

ax_kidney_ct_neg = plt.subplot(outer_gs_ctneg[0])
ax_bladder_ct_neg = plt.subplot(outer_gs_ctneg[1], sharey = ax_kidney_ct_neg, sharex = ax_kidney_ct_neg)

ax_kidney_ct_neg = plot_overlapping_insert_size_histograms(ctdna_neg_kidney_patients, insert_size_paths, ax = ax_kidney_ct_neg, alpha = 0.05, color = "orangered")
ax_bladder_ct_neg = plot_overlapping_insert_size_histograms(ctdna_neg_bladder_patients, insert_size_paths, ax = ax_bladder_ct_neg, alpha = 0.05, color = "deepskyblue")

ax_kidney_ct_neg.set_title("ctDNA- kidney")
ax_bladder_ct_neg.set_title("ctDNA- bladder")
fig_ctneg.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/insert_size/insert_size_ctDNA_NEG.pdf")

# 3. SELECT KIDNEY SAMPLES WE RAN ON THE TAPE STATION
# Plots samples individually
pt_names = ["21-041", "23-085", "20-038", "23-492", "20-394", "23-499", "20-221", "22-526", "21-353", "23-536", "17-094", "19-436"]
colors = ["Orangered", "Orangered", "Orangered", "Orangered", "Orangered", "Orangered", "Orangered", "Orangered", "Orangered", "Orangered", "Orangered", "Deepskyblue"]
# colors = cm.tab20(np.linspace(0, 1, len(pt_names)))

# Set up the GridSpec layout
n_cols = 3  # Number of columns in the grid
n_rows = 4  # Calculate required rows
fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4, n_rows * 3), sharex=True, sharey=True)
axes = axes.flatten() # Flatten the axes array for easier indexing

# Loop over patients and plot each one in the appropriate subplot
for i, (pt, color) in enumerate(zip(pt_names, colors)):
    ax = axes[i]
    print(pt)
        
    # Plotting
    path = [f for f in insert_size_paths if pt in f and "cfDNA" in f][0]
    df_fractions = convert_to_fraction(path)
    ax.bar(df_fractions["insert_size"], df_fractions["Fraction"], color=color, width = 0.1, edgecolor = "None")
    
    # Line through mode
    max_val = df_fractions[df_fractions["Fraction"] == df_fractions["Fraction"].max()]["insert_size"].reset_index(drop = True)[0]
    ax.axvline(x=max_val, linestyle='--', color='black', lw = 0.5)
    ax.text(max_val+10, ax.get_ylim()[1]*0.8, max_val, ha='left', va='center', color="black", fontsize=10)
    
    # Aes
    ax.set_title(pt)
    ax.set_xlim(0, 700)
    ax.spines[["top", "right"]].set_visible(False)
    
    if i in [9, 10, 11]: 
        ax.set_xlabel("Insert size")
    if i in [0, 3, 6, 9]:
        ax.set_ylabel("Fraction of insert sizes")

fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/insert_size/tape_station_samples.pdf")

# 4. GROUPING SAMPLES BASED ON METS STATUS
# Makes overlapping plots based on mets status.
path_clinical_data = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/Supplementary tables - Clinical data - mRCC.csv"
clinical_data= pd.read_csv(path_clinical_data)

