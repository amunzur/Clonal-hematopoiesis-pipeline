"""
This script provides some stats and figures on the SSCS1 vs SSCS2 comparison.
"""

import os 
import pandas as pd
import numpy as np
import re
from collections import Counter
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
from math import floor



DIR_UMI_metrics = "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/umi_metrics"
DIR_bladder_fastq = "/groups/wyattgrp/users/amunzur/pipeline/results/data/fastq/merged/bladder" # helps get a list of bladder patients
DIR_kidney_fastq = "/groups/wyattgrp/users/amunzur/pipeline/results/data/fastq/merged/kidney"
DIR_figures = "/groups/wyattgrp/users/amunzur/pipeline/results/figures/family_size"
PATH_depth = "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/seq_quality_metrics/seq_quality_metrics.tsv"

bladder_patients = pd.read_csv("/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_list_bladder.tsv")
kidney_patients = pd.read_csv("/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_list_kidney.tsv")

all_umi_metrics_files = os.listdir(DIR_UMI_metrics)

def get_family_df(all_umi_metrics_files, patients, family_size): 
    """
    Get the umi metrics files from bladder samples, and return one df containing the requested family size from all bladder samples. Used for histogram plotting.
    """
    df_list = []
    for file in all_umi_metrics_files:
        substring = re.split(r"_cfDNA|_WBC", file)[0]
        if any(substring in patient_id for patient_id in patients["sample_names"]):
            df = pd.read_csv(os.path.join(DIR_UMI_metrics, file), sep="\t")
            df["Sample"] = file.split(".")[0]
            df_list.append(df)
    
    df = pd.concat(df_list).reset_index(drop = True)
    family_df = df[df["family_size"] == family_size].reset_index(drop = True)
    family_df["type"] = family_df["Sample"].str.extract("(WBC|gDNA|cfDNA)", expand = False)
    family_df.loc[family_df["type"] == "gDNA", "type"] = "WBC" # replace gDNA with WBC        
    return(family_df)

def make_hist(values, ax, colors, legend = False):
    """
    Plot the histogram of the fraction of reads with a given family size to the given axis.
    """
    ax.hist(values, bins=15, color = colors)
    # First dataset (cfDNA)
    median_value = values[0].median()
    ax.axvline(median_value, color=colors[0], linestyle='dashed', linewidth=1, label='Med cfDNA')
    ax.text(median_value + 0.001, ax.get_ylim()[1] * 0.9, f'Med = {median_value:.3f}', color=colors[0])
    # second dataset (WBC)
    median_value = values[1].median()
    ax.axvline(median_value, color=colors[1], linestyle='dashed', linewidth=1, label='Med WBC')
    ax.text(median_value + 0.01, ax.get_ylim()[1] * 0.95, f'Med = {median_value:.3f}', color=colors[1])
    # Aesthetics
    ax.set_yticks(list(range(0, floor(ax.get_ylim()[1]), 10)))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylabel("Number of samples")
    ax.set_xlabel("Fraction of all families")
    if legend: 
        handles = [Rectangle((0,0),1,1, color=c, ec="k") for c in colors]
        labels= ["cfDNA", "WBC"]
        ax.legend(handles, labels, bbox_to_anchor=(0.7, 0.6), loc="upper center", borderaxespad=0.)
    return(ax)

# FIGURE 1
# Histogram showing the percentage of reads with a family size of 1, 2, 3 and 4
fig = plt.figure(figsize = (7, 7))
gs = gridspec.GridSpec(3, 2, width_ratios=[1, 1], height_ratios=[1, 1, 1])

for number in [0, 1, 2, 3, 4, 5]: 
    ax = fig.add_subplot(gs[number]) # gene CN and mutation status
    ax.set_title(f"Family size = {number + 1}")    
    family_df = get_family_df(all_umi_metrics_files, kidney_patients, family_size = number + 1)
    family_df_cfDNA = family_df[family_df["type"] == "cfDNA"].reset_index(drop = True)
    family_df_WBC = family_df[family_df["type"] == "WBC"].reset_index(drop = True)
    if number == 0: 
        legend = True
    else: 
        legend = False
    make_hist([family_df_cfDNA["fraction"], family_df_WBC["fraction"]], ax, colors = ["black", "grey"], legend = legend)

gs.tight_layout(fig)
fig.savefig(os.path.join(DIR_figures, "kidney_read_family_histograms.png"))
fig.savefig(os.path.join(DIR_figures, "kidney_read_family_histograms.pdf"))

# FIGURE 2
# SOME STATS on family size = 1
family1 = get_family_df(all_umi_metrics_files, kidney_patients, family_size = 1).sort_values(by = "fraction")
family1["type"] = family1["Sample"].str.extract("(WBC|gDNA|cfDNA)", expand = False)
family1.loc[family1["type"] == "gDNA", "type"] = "WBC" # replace gDNA with WBC
family1_filter = family1[family1["fraction"] > 0.3]
family1_filter["type"].value_counts()

# Exploring the distribution of samples with a fraction higher than the median, which is 10%.
fig, ax = plt.subplots(figsize=(3, 5))
ax.bar(cfDNA_dict.keys(), cfDNA_dict.values(), color='gray', edgecolor='black', width = 0.6)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlabel('')
ax.set_ylabel('Number of samples')
ax.set_xticklabels(["Baseline cfDNA\nsamples", "On treatment cfDNA\nsamples"])
ax.set_title('Bladder cfDNA samples')

fig.tight_layout()
fig.savefig(os.path.join(DIR_figures, "baseline_and_OT_bladder_samples.png"))