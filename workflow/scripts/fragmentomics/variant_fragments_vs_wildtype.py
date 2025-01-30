"""
For each ctDNA and CH mutation this script asks the question, just based on fragment length is it possible to distingush CH 
from normal DNA?
"""
# ca snakemake

import pandas as pd
import numpy as np
import os
import ast
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import matplotlib.gridspec as gridspec
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

# conda activate pysam_environment
# DIR_working = "/groups/wyattgrp/users/amunzur/pipeline"
# PATH_sample_information = os.path.join(DIR_working, "resources/sample_lists/sample_information.tsv")
DIR_output = "/groups/wyattgrp/users/amunzur/pipeline/results/fragmentomics"

# sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
# sample_info = sample_info[sample_info["Diagnosis"] != "Healthy"]

# all_vars_chip = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/chip.csv"))
# all_vars_chip = sample_info.merge(all_vars_chip, how = "inner")
# all_vars_chip = all_vars_chip[all_vars_chip["Timepoint"] == "Baseline"].reset_index(drop = True)

# all_vars_ctDNA = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/somatic.csv"))
# all_vars_ctDNA = sample_info.merge(all_vars_ctDNA, how = "inner")
# all_vars_ctDNA = all_vars_ctDNA[all_vars_ctDNA["Timepoint"] == "Baseline"].reset_index(drop = True)

# Now for each CH mutation, ask the question: Can we distinguish CH from WT fragments based on fragment length only? 


p_value_list = []
for i, row in df_ctDNA.iterrows():
    ctDNA_fragments = ast.literal_eval(row["ctDNA fragments"])
    wt_fragments = ast.literal_eval(row["WT fragments"])
    statistic, p_value = mannwhitneyu(ctDNA_fragments, wt_fragments)
    p_value_list.append(p_value)

_, corrected_p_values_ctDNA, _, _ = multipletests(p_value_list, method='fdr_bh')

len(corrected_p_values_ctDNA[corrected_p_values_ctDNA<0.05])



def compare_mutated_to_WT(df_path, genes = None):
    """
    Given a muts_df, with CH and WT fragments added as two additional columns
    """
    df = pd.read_csv(df_path)
    
    if genes is not None:
        if isinstance(genes, str): 
            df = df[df["Gene"] == genes]
        elif isinstance(genes, list): 
            df = df[df["Gene"].isin(genes)]
    
    if "ctDNA fragments" in df.columns:
        col_to_use = "ctDNA fragments"
    else: 
        col_to_use = "CH fragments"
    
    # Go through each mut and compare the WT to mutated
    p_value_list = []
    for i, row in df.iterrows():
        mutated_fragments = ast.literal_eval(row[col_to_use])
        wt_fragments = ast.literal_eval(row["WT fragments"])
        statistic, p_value = mannwhitneyu(mutated_fragments, wt_fragments)
        p_value_list.append(p_value)
    
    _, corrected_p_values, _, _ = multipletests(p_value_list, method='fdr_bh')
    
    # add to the df
    df["MWU p value"] = p_value_list
    df["corrected MWU p value"] = corrected_p_values
    
    n_sig = len(corrected_p_values[corrected_p_values<0.05])
    
    print(f"There are {n_sig} significant p values.")
    return(df)

def plot_mutated_median_vs_wt_median(df, ax, title, set_y_label = True, add_legend = True):
    if "ctDNA fragments" in df.columns:
        col_to_use = "ctDNA fragments"
    else: 
        col_to_use = "CH fragments"
    
    color_dict = {True: "cornflowerblue", False: "black"}
    zorder_dict = {True: 50, False: 49}
    df["sig"] = df["corrected MWU p value"] < 0.05
    df["color"] = df["sig"].map(color_dict)
    df["zorder"] = df["sig"].map(zorder_dict)
    
    for i, row in df.iterrows():
        mutated_fragments = ast.literal_eval(row[col_to_use])
        wt_fragments = ast.literal_eval(row["WT fragments"])
        
        mutated_median = np.median(mutated_fragments)
        wt_median = np.median(wt_fragments)
        color = row["color"]
        
        ax.scatter(mutated_median, wt_median, alpha=0.7, color=color, s=4, edgecolor = None, zorder = row["zorder"])
    
    # set ax lims
    ax_max = max(max(ax.get_xlim()), max(ax.get_ylim()))
    ax_min = min(max(ax.get_xlim()), min(ax.get_ylim())) 
    ax.set_xlim((ax_min-20, ax_max))
    ax.set_ylim((ax_min-10, ax_max))
    
    # Aes
    if set_y_label: 
        ax.set_ylabel('Median length of wildtype fragments', fontsize = 9)
    ax.set_title(title, fontsize = 9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel("Median length of mutated fragments")
    
    if add_legend:
        legend_colors = ["cornflowerblue", "black"]
        legend_labels = ["p < 0.05", "p > 0.05"]
        legend_handles = [plt.Line2D([0], [0], color=color, marker='o', label=label, linestyle='None', markersize=5) for color, label in zip(legend_colors, legend_labels)]
        ax.legend(handles=legend_handles, loc="upper right", frameon=False, handletextpad=0.2)
    
    return(ax)

def plot_gs(path_save, fig_title = None):
    fig = plt.figure(figsize=(6, 2.5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])
    ax1 = fig.add_subplot(gs[0])
    ax1 = plot_mutated_median_vs_wt_median(df_ctDNA, ax=ax1, title = "ctDNA mutations", add_legend = False)
    ax2 = fig.add_subplot(gs[1])
    ax2 = plot_mutated_median_vs_wt_median(df_CH, ax=ax2, title = "CH mutations", set_y_label = False)
    
    fig.text(0.5, 0.02, 'Mutated fragments', ha='center', fontsize = 9)
    gs.tight_layout(fig)
    
    if fig_title is not None:
        fig.suptitle(fig_title, fontsize = 9)
    
    fig.savefig(path_save)

# ctDNA_nsig = df_ctDNA[df_ctDNA["corrected MWU p value"] < 0.05].shape[0]
# ctDNA_nsig/df_ctDNA.shape[0]*100

# CH_nsig = df_CH[df_CH["corrected MWU p value"] < 0.05].shape[0]
# CH_nsig/df_CH.shape[0]*100

ctDNA_path = os.path.join(DIR_output, "ctDNA.csv")
df_ctDNA = compare_mutated_to_WT(ctDNA_path)
CH_path = os.path.join(DIR_output, "CH.csv")
df_CH = compare_mutated_to_WT(CH_path)
plot_gs("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/fragmentomics/mutated_vs_WT_scatter.png")

ctDNA_path = os.path.join(DIR_output, "ctDNA.csv")
df_ctDNA = compare_mutated_to_WT(ctDNA_path, genes = "TP53")
CH_path = os.path.join(DIR_output, "CH.csv")
df_CH = compare_mutated_to_WT(CH_path, genes = "TP53")
plot_gs("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/fragmentomics/mutated_vs_WT_scatter_TP53.png", fig_title = "TP53")

ctDNA_path = os.path.join(DIR_output, "ctDNA.csv")
df_ctDNA = compare_mutated_to_WT(ctDNA_path, genes = ["TP53", "BRCA2", "BRCA1", "CHEK2", "ATM"])
CH_path = os.path.join(DIR_output, "CH.csv")
df_CH = compare_mutated_to_WT(CH_path, genes = ["TP53", "BRCA2", "BRCA1", "CHEK2", "ATM"])
plot_gs("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/fragmentomics/mutated_vs_WT_scatter_DDR.png", fig_title = "DDR genes")