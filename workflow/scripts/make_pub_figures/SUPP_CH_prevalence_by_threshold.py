"""
Shows the prevalence of CH in the cohort across various minimum detection thresholds.
"""

import pandas as pd
import numpy as np
import os 
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
from scipy.stats import fisher_exact
from scipy.stats import fisher_exact # needs to stay twice

mpl.rcParams['font.size'] = 10

DIR_working = "/groups/wyattgrp/users/amunzur/pipeline"
PATH_sample_information = os.path.join(DIR_working, "resources/sample_lists/sample_information.tsv")
figure_dir = os.path.join(DIR_working, "results/figures/amazing_figures")
source_functions = os.path.join(DIR_working, "workflow/scripts/visualization/UTILITIES_make_chip_plots.py")
sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])

color_dict = {"Bladder": "deepskyblue", "Kidney": "orangered"}

with open(source_functions, 'r') as file:
    script_code = file.read()

exec(script_code)

# LOAD CHIP DATASETS
all_vars_chip = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/chip.csv"))
all_vars_chip = all_vars_chip[all_vars_chip["Dependent"] == False]
base_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
prog_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)
prog_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)

baseline = pd.concat([base_kidney_chip, base_bladder_chip]).reset_index(drop = True)
baseline_somatic = pd.concat([base_kidney_somatic, base_bladder_somatic]).reset_index(drop = True)

def plotting(all_vars_chip, ax, PATH_sample_information):
    sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
    n_total = sample_info["Patient_id"].unique().shape[0]
    
    thresholds = [0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2]
    
    # Determine genes lists
    all_genes = all_vars_chip["Gene"].unique()
    DTA_genes = ["DNMT3A", "TET2", "ASXL1"]
    DDR_genes = ["TP53", "BRCA1", "BRCA2", "ARID1A", "ATM", "CHEK2"]
    genes_list = [all_genes, DTA_genes, DDR_genes]
    
    colors = ["black", "hotpink", "limegreen"]
    
    for genes, color in zip(genes_list, colors):
        perc_pos_list = []  # to store percentage values for each threshold
        for threshold in thresholds:
            ch_calls = all_vars_chip[(all_vars_chip["Gene"].isin(genes)) & (all_vars_chip["VAF_n"] >= threshold)]
            n_pts_with_ch = ch_calls["Patient_id"].unique().shape[0]
            perc_pos = n_pts_with_ch / n_total * 100
            perc_pos_list.append(perc_pos)
            ax.scatter(threshold, perc_pos, color=color, s=50)
        
        # Plot lines connecting the dots
        ax.plot(thresholds, perc_pos_list, color=color, linestyle='-', linewidth=2)
    
    ax.set_yticks([0, 10, 20, 30, 40, 50, 60, 70, 80])
    ax.set_yticklabels(["0", "10", "20", "30", "40", "50", "60", "70", "80"], fontsize = 10)
    ax.set_ylabel("% CH+ in the cohort", fontsize = 10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlabel("Minimum CH detection threshold %", fontsize = 10)
    ax.tick_params(axis='x', labelsize = 10)
    
    # legend
    legend_labels = ["All genes", "DTA genes", "DDR genes"]
    legend_colors = ["black", "hotpink", "limegreen"]
    legend_handles = [plt.Line2D([0], [0], marker='o', color=color, label=label, markersize=3, linestyle='') for color, label in zip(legend_colors, legend_labels)]
    ax.legend(handles=legend_handles, loc="upper right", frameon=False, handletextpad=0.3)
    
    return(ax)

fig = plt.figure(figsize=(8, 6))
gs = gridspec.GridSpec(2, 1, height_ratios=[0.7, 1], hspace = 0.2, wspace = 0) # outer most gs with 3 rows
ax_threshold = plt.subplot(gs[0])
ax_threshold = plotting(all_vars_chip, ax_threshold, PATH_sample_information)

# Now variance
n_pts_kidney = base_kidney_chip["Patient_id"].unique().size
n_pts_bladder = base_bladder_chip["Patient_id"].unique().size

gs_variance = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=[n_pts_kidney, n_pts_bladder], hspace = 0.5, wspace = 0.03, subplot_spec=gs[1])
kidney_ax = plt.subplot(gs_variance[0])
bladder_ax = plt.subplot(gs_variance[1])
kidney_ax, kidney_twin_ax = plot_variance_per_pt(baseline[baseline["Diagnosis"] == "Kidney"], kidney_ax, scatter_color = "orangered", bar_color = "#FFEDE6")
bladder_ax, bladder_twin_ax = plot_variance_per_pt(baseline[baseline["Diagnosis"] == "Bladder"], bladder_ax, scatter_color = "deepskyblue", bar_color = "#E8F9FF")

# minor aes
kidney_twin_ax.set_ylabel("")
kidney_twin_ax.set_yticks([])
kidney_twin_ax.set_yticklabels([])
kidney_twin_ax.spines["right"].set_visible(False)
kidney_ax.spines["right"].set_visible(False)
kidney_ax.set_xlabel("mRCC")

bladder_ax.set_ylabel("")
bladder_ax.set_yticks([])
bladder_ax.set_yticklabels([])
bladder_ax.spines["left"].set_visible(False)
bladder_twin_ax.spines["left"].set_visible(False)
bladder_ax.set_xlabel("mUC")

fig.tight_layout()
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/SUPP_CH_prevalence_threshold.png")
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/SUPP_CH_prevalence_threshold.pdf")