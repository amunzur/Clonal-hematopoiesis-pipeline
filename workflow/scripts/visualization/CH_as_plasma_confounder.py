"""
Explores various plotting options and figures regarding CH as a confounder of plasma ctDNA analysis.
"""
import pandas as pd
import numpy as np
import os 
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import linregress
import matplotlib.gridspec as gridspec
import seaborn as sns


def bin_ct_fractions(PATH_mutation_ctfractions, all_vars_chip):
    """
    For each sample add the ctDNA fraction and place the same into the appropriate bin.
    """
    ct_df = pd.read_csv(PATH_mutation_ctfractions)[["Patient_id", "Diagnosis", "Date_collected", "Mutation_ctDNA_fraction"]]
    ct_df['Date_collected'] = pd.to_datetime(ct_df['Date_collected']).dt.strftime('%Y%b%d')
    all_vars_chip = all_vars_chip.merge(ct_df, how = "left", on = ["Patient_id", "Diagnosis", "Date_collected"])
    all_vars_chip["Mutation_ctDNA_fraction"] = all_vars_chip["Mutation_ctDNA_fraction"]*100
    
    bins = [0, 2, 30, 100]
    labels = ['0-2', '2-30', '30-100']
    all_vars_chip['Mutation_ctDNA_bin'] = pd.cut(all_vars_chip['Mutation_ctDNA_fraction'], bins=bins, labels=labels, right=False)
    
    return(all_vars_chip)

def correlate_VAFt_and_VAFn(all_vars_chip):
    """
    Correlates VAF_t and VAF_n in df_muts. Returns the pearson R
    """
    labels = ['0-2', '2-30', '30-100'] # labels for ctDNA groupings
    
    # Correlate WBC VAf with cfDNA VAF in each bin
    correlations = {}
    for label in labels:
        bin_data = all_vars_chip[all_vars_chip['Mutation_ctDNA_bin'] == label] # Filter the DataFrame for the current bin
        correlation = bin_data['VAF_t'].corr(bin_data['VAF_n']) # Calculate Pearson correlation between 'VAF_t' and 'VAF_n' in the current bin
        correlations[label] = correlation # Store the correlation in the dictionary
    return(correlations)

def plot_ctDNA_bins_and_VAF_correlations(correlations, all_vars_chip, dir_figures, ext = "png"):    
    # Set up the plot
    bin_labels = ['0-2', '2-30', '30-100']
    bin_numeric = [1, 2, 3]  # Numeric representation of the bins
    corr_values = [correlations[label] for label in bin_labels]
    
    fig = plt.figure(figsize=(3, 3))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
    
    ax1 = plt.subplot(gs[0]) # dotplot
    ax2 = plt.subplot(gs[1]) # boxplot
    
    # Perform linear regression
    slope, intercept, r_value, p_value, std_err = linregress(bin_numeric, corr_values)
    
    # Plot the connected line plot for the pearson
    ax1.scatter(labels, correlations.values(), s=70, color = "black")
    ax1.plot(bin_numeric, intercept + slope * np.array(bin_numeric), 'r', label=f'y={intercept:.2f}+{slope:.2f}x\n(p={p_value:.3f})')
    ax1.set_ylabel("Pearson R")
    ax1.set_xlabel("")
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.set_yticks([1, 0.98, 0.96, 0.94, 0.92])
    ax1.set_yticklabels(["1", "0.98", "0.96", "0.94", "0.92"])
    ax1.set_xticks([])
    ax1.tick_params(axis="both", direction="out", which="both", left=True, bottom=False, colors='k')
    
    # Plot the boxplots with scatters
    sns.boxplot(x='Mutation_ctDNA_bin', y="Mutation_ctDNA_fraction", data=all_vars_chip, order=['0-2', '2-30', '30-100'], ax=ax2, color='black', medianprops={'color': 'black'}, fliersize = 0, boxprops={"alpha": .2, "edgecolor": 'black'})
    sns.stripplot(x='Mutation_ctDNA_bin', y="Mutation_ctDNA_fraction", data=all_vars_chip, order=['0-2', '2-30', '30-100'], ax=ax2, jitter=True, size = 3, color='black')
    ax2.set_ylabel("ctDNA%")
    ax2.set_xlabel("Binned ctDNA fractions")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.set_yticks([0, 25, 50, 75, 100])
    ax2.set_yticklabels(["0", "25", "50", "75", "100"])
    ax2.tick_params(axis="both", direction="out", which="both", left=True, bottom=True, colors='k')
    fig.tight_layout()
    fig.savefig(os.path.join(dir_figures, f"binned_VAF_correlation_plot.{ext}"))

def plot_wbc_and_cfDNA_vaf_correlation(all_vars_chip, dir_figures, ext = "pdf", vaf_n_limit = None):
    """
    For each ctDNA bin plot the WBC and cfDNA vaf scatter showing how well they are correlated.
    vaf_limit is the highest VAF that will be allowed. Give an integer not fraction.
    """
    if vaf_n_limit is not None:
        all_vars_chip = all_vars_chip[all_vars_chip["VAF_n"] <= vaf_n_limit]
    fig = plt.figure(figsize=(7, 2.8))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1])
    
    color_dict = {"Bladder": "deepskyblue", "Kidney": "orangered"}
    
    bin_labels = ['0-2', '2-30', '30-100']
    for i, label in enumerate(bin_labels):
        df = all_vars_chip[all_vars_chip["Mutation_ctDNA_bin"] == label]
        ax = plt.subplot(gs[i])
        sns.scatterplot(x="VAF_n", y="VAF_t", data=df, ax=ax, s=9, hue="Diagnosis", edgecolor = None, palette = color_dict, legend=False)
        sns.regplot(x="VAF_n", y="VAF_t", data=df, ax=ax, scatter = False, line_kws={'color': 'black', 'linewidth': 1, "linestyle": "--"})
        spearman_corr = df[["VAF_n", "VAF_t"]].corr(method="spearman").iloc[0, 1]
        ax.text(0.05, 0.85, f'Spearman\'s\ncorrelation: {spearman_corr:.2f}', transform=ax.transAxes, fontsize=8, color='black')
        ax.set_title(f"ctDNA between {label}%", fontsize = 8, pad = -5)
        # ax.set_xlim((0, 80))
        # ax.set_ylim((0, 80))
        if i == 0:
            ax.set_xlabel("")
            ax.set_ylabel("Plasma cfDNA VAF %")
        elif i == 1:
            ax.set_xlabel("WBC VAF %")
            ax.set_ylabel("")
        else:
            ax.set_xlabel("")
            ax.set_ylabel("")
        ax.spines[["top", "right"]].set_visible(False)
    if vaf_n_limit is not None:
        fig.suptitle(f"Only including mutations with WBC VAF > {vaf_n_limit}%", fontsize = 10)
        gs.tight_layout(fig)
        fig.savefig(os.path.join(dir_figures, f"binned_VAF_correlation_plot_scatter_VAF{vaf_n_limit}.{ext}"))
    else:
        fig.suptitle(f"Including all CH mutations", fontsize = 10)
        gs.tight_layout(fig)
        fig.savefig(os.path.join(dir_figures, f"binned_VAF_correlation_plot_scatter.{ext}"))

DIR_working = "/groups/wyattgrp/users/amunzur/pipeline"
PATH_sample_information = os.path.join(DIR_working, "resources/sample_lists/sample_information.tsv")
PATH_mutation_ctfractions = os.path.join(DIR_working, "results/variant_calling/ctDNA_fractions.csv")
DIR_figures = os.path.join(DIR_working, "results/figures/amazing_figures")

all_vars_chip = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/chip.csv"))
sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
all_vars_chip = sample_info.merge(all_vars_chip, how = "inner")
all_vars_chip = all_vars_chip[all_vars_chip["Dependent"] == False][['Patient_id', 'Diagnosis', 'Timepoint', 'Date_collected', 'Chrom', 'Position', 'Ref', 'Alt', 'Type', 'VAF_t', 'VAF_n']]

# add ct fraction
all_vars_chip = bin_ct_fractions(PATH_mutation_ctfractions, all_vars_chip)

# correlate ctDNA with CH
correlations = correlate_VAFt_and_VAFn(all_vars_chip)

# plot the correlations
plot_ctDNA_bins_and_VAF_correlations(correlations, all_vars_chip, DIR_figures, ext = "png")

# Plot individual scatters and correlations depending on the ctDNA fraction
plot_wbc_and_cfDNA_vaf_correlation(all_vars_chip, dir_figures, ext = "pdf")
plot_wbc_and_cfDNA_vaf_correlation(all_vars_chip, dir_figures, ext = "pdf", vaf_n_limit = 2)
plot_wbc_and_cfDNA_vaf_correlation(all_vars_chip, dir_figures, ext = "pdf", vaf_n_limit = 10)
plot_wbc_and_cfDNA_vaf_correlation(all_vars_chip, dir_figures, ext = "pdf", vaf_n_limit = 20)