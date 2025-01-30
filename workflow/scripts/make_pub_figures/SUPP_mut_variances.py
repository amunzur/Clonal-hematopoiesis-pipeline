"""
For each patient with more than one CH mutation, plots the variance.
"""

import pandas as pd
import numpy as np
import os 
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from itertools import cycle
from scipy.stats import ttest_ind
import seaborn as sns

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

DIR_working = "/groups/wyattgrp/users/amunzur/pipeline"
PATH_sample_information = os.path.join(DIR_working, "resources/sample_lists/sample_information.tsv")
figure_dir = os.path.join(DIR_working, "results/figures/pub_figures")
PATH_kidney_clinical = os.path.join(DIR_working, "resources/clinical_data/RCC clinical data - mRCC clinical Data.csv")
PATH_clinical_bladder = os.path.join(DIR_working, "resources/clinical_data/bladder/clinical_data.csv")

# LOAD CHIP DATASETS
all_vars_chip = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/chip.csv"))
baseline = all_vars_chip[(all_vars_chip["Dependent"] == False) & (all_vars_chip["Timepoint"] == "Baseline")]

kidney_age = pd.read_csv(PATH_kidney_clinical)[["Patient_id", "Age at GUBB draw"]].rename(columns = {"Age at GUBB draw": "Age"}).assign(Diagnosis="Kidney")
bladder_age = pd.read_csv(PATH_clinical_bladder)[["Patient_id", "Age at blood draw"]].rename(columns = {"Age at blood draw": "Age"}).assign(Diagnosis="Bladder")
age_df = pd.concat([kidney_age, bladder_age]).reset_index(drop=True)

def plot_variance_per_pt(baseline, ax, scatter_color, bar_color):
    baseline.loc[baseline["Gene"] == "ZRSR2", "VAF_n"] = baseline.loc[baseline["Gene"] == "ZRSR2", "VAF_n"]/2  
    
    # Calculate coefficient of variance per pt
    grouped = baseline.groupby("Patient_id")["VAF_n"]
    cv_per_pt = (grouped.std() / grouped.mean()).reset_index(name='CV')
    
    cv_per_pt.columns = ["Patient_id", "CV"]
    cv_per_pt.sort_values(by = "CV", ascending = False, inplace = True)
    cv_per_pt = cv_per_pt.dropna(subset = "CV").reset_index(drop = True)
    cv_per_pt = cv_per_pt.reset_index()
    
    # Pull the WBC VAF of all mutations in these patients
    vaf_df = baseline[["Patient_id", "VAF_n"]]
    vaf_df = vaf_df[vaf_df["Patient_id"].isin(cv_per_pt["Patient_id"])].reset_index(drop = True)
    
    # merge
    cv_per_pt = cv_per_pt.merge(vaf_df)
    
    # plotting
    ax_twin = ax.twinx()
    for i, group in cv_per_pt.groupby("Patient_id"):
        var = group["CV"].unique()[0]
        x_pos = group["index"].unique()[0]
        ax.bar(x_pos, var, color=bar_color, edgecolor = "None")
        ax_twin.scatter(np.repeat(x_pos, len(group)), group["VAF_n"], s = 2, edgecolor = None, color = scatter_color, alpha = 0.7, zorder = 100)
    
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_zorder(1000)
    ax.set_ylabel("Coefficient of variance")
    ax_twin.set_ylim((0, 50))
    ax_twin.set_xlim()
    
    ax_twin.spines["top"].set_visible(False)
    ax_twin.set_ylabel("WBC VAF%")
    ax_twin.spines["bottom"].set_zorder(1000)
    
    # Set x limits
    xmin = cv_per_pt["index"].min() - 1
    xmax = cv_per_pt["index"].max() + 1
    ax.set_xlim((xmax, xmin))
    ax_twin.set_xlim((xmax, xmin))
    ax_twin.set_xlabel("Baseline samples with >1 CH mutation")
    
    return(ax, ax_twin)

fig = plt.figure(figsize=(8, 3))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1], hspace = 0, wspace = 0) # outer most gs with 3 rows

# FIGURE 1B.
kidney_ax = plt.subplot(gs[0])
bladder_ax = plt.subplot(gs[1])

kidney_ax, kidney_twin_ax = plot_variance_per_pt(baseline, kidney_ax, scatter_color = "orangered", bar_color = "#ffd2c2")
bladder_ax, bladder_twin_ax = plot_variance_per_pt(baseline, bladder_ax, scatter_color = "deepskyblue", bar_color = "#d1f4ff")

# minor aes
kidney_twin_ax.set_ylabel("")
kidney_twin_ax.set_yticks([])
kidney_twin_ax.set_yticklabels([])
kidney_twin_ax.spines["right"].set_visible(False)
kidney_ax.spines["right"].set_visible(False)
kidney_ax.set_xlabel("mRCC samples with >1 CH mutation")

bladder_ax.set_ylabel("")
bladder_ax.set_yticks([])
bladder_ax.set_yticklabels([])
bladder_ax.spines["left"].set_visible(False)
bladder_twin_ax.spines["left"].set_visible(False)
bladder_ax.set_xlabel("mUC samples with >1 CH mutation")

gs.tight_layout(fig)
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/mut_variance.png")
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/mut_variance.pdf")

