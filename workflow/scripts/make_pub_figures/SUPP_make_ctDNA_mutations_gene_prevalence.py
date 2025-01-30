"""
Simple bar chart. Kidney and bladder separate. 
Shows the percent mutated in the cohort, swarm plot indicates the VAF.
"""

import pandas as pd
import numpy as np
import os 
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from scipy.stats import mannwhitneyu

PATH_sample_information = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"
PATH_ctDNA_fraction = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/ctDNA_fractions.csv"

# LOAD SOMATIC DATASETS
all_vars_somatic = pd.read_csv(os.path.join("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/somatic.csv"))
all_vars_somatic = all_vars_somatic[~all_vars_somatic["Patient_id"].isin(["20-313", "21-184", "21-430"])] # exclude some samples due to oxidative damage
all_vars_somatic = all_vars_somatic[all_vars_somatic["Dependent"] == False]
base_kidney_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Kidney")].reset_index(drop = True)
prog_kidney_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "During treatment") & (all_vars_somatic["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Bladder")].reset_index(drop = True)
prog_bladder_somatic = all_vars_somatic[(all_vars_somatic["Tibase_bladder_somaticmepoint"] == "During treatment") & (all_vars_somatic["Diagnosis"] == "Bladder")].reset_index(drop = True)

baseline_ctDNA = pd.concat([base_kidney_somatic, ]).reset_index(drop = True)

def make_percent_mutated_barchart(muts_df, bar_color, ax, ymax):
    diagnosis = muts_df["Diagnosis"].unique()[0]
    denom = muts_df["Patient_id"].unique().shape[0] # Using ctDNA+ patients as denominator
    
    # subset and go from counts to fractions
    genes_df = muts_df[["Patient_id", "Gene"]].drop_duplicates().reset_index(drop = True)["Gene"].value_counts().reset_index().rename(columns = {"index": "Gene", "Gene": "Counts"}).reset_index()
    genes_df["frac"] = (genes_df["Counts"]/denom)*100
    genes_df = genes_df.iloc[0:9, ] # We show only top 10 genes
    vaf = muts_df[["Gene", "VAF_t"]].reset_index(drop = True).merge(genes_df[["Gene", "index"]])    
    
    # plotting
    # fig, ax = plt.subplots(figsize = (5.6, 2.8))
    ax2 = ax.twinx()
    ax.bar(genes_df["index"], genes_df["frac"], color=bar_color, width = 0.8, edgecolor = "None")
    # plot the vafs for kidney
    for i, row in vaf.iterrows():
        jitter = np.random.uniform(-0.2, 0.2, 1)
        ax2.scatter(row["index"]+jitter, row["VAF_t"], color="black", s = 2, alpha = 0.6)    
    
    # aesthetics
    ax.spines["top"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    
    ax.set_xticks(genes_df["index"].tolist())
    ax.set_xticklabels(genes_df["Gene"].tolist(), rotation = 90, fontstyle = "italic")
    ax.tick_params(axis = "x", which='both', length=0)
    ax.set_ylabel("% of ctDNA+ patients")
    ax.set_ylim((0, ymax))
    ax.set_xlim((-0.6, 8.6))
    ax2.set_ylabel("Plasma VAF%")
    ax2.set_ylim((0, 100))
    ax2.set_yticks([0, 20, 40, 60, 80, 100])
    ax2.set_yticklabels(["0", "20", "40", "60", "80", "100"])
    ax2.tick_params(axis='x', pad=1)
    
    # Annotate the numbers for each bar
    ypos = ymax-1
    for i, row in genes_df.iterrows():
        ax.text(row["index"], ypos, row["Counts"], fontsize=6, color='black',  ha='center', va='center')
    
    return(ax, ax2)

def make_min_vaf_ctDNA_prevalence_barchart(muts_df, PATH_sample_information, ax, bar_color):
    """
    Barchart showing the prevalence of ctDNA mutations at different VAF thresholds.
    """
    diagnosis = muts_df["Diagnosis"].unique()[0]
    if diagnosis == "Bladder":
        n_total = 113 # excluding oxidative damage samples
        ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8])
        ax.set_yticklabels(["0", "20", "40", "60", "80"])
        # ax.set_title("mUC")
    else:
        n_total = 186
        ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
        ax.set_yticklabels(["0", "10", "20", "30", "40", "50"])
        # ax.set_title("mRCC")
    
    value1 = muts_df["Patient_id"].unique().shape[0] # As is
    value2 = muts_df[muts_df["VAF_t"] >=1]["Patient_id"].unique().shape[0] # Only muts with VAF > 1%
    value3 = muts_df[muts_df["VAF_t"] >=10]["Patient_id"].unique().shape[0] # Only muts with VAF > 10%
    ax.bar([0, 1, 2], [value1/n_total, value2/n_total, value3/n_total], color=bar_color, width = 0.8, edgecolor = "None")
    
    # Annotate ns
    ypos = ax.get_ylim()[1]-1
    for x_value, txt in zip([0, 1, 2], [value1, value2, value3]):
        ypos = txt/n_total+0.01
        ax.text(x_value, ypos, str(txt), fontsize=6, color='black',  ha='center', va='bottom')
    # Aes
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["≥0.5%", "≥1%", "≥10%"], fontsize = 7, rotation = 45)
    ax.set_ylabel("% of ctDNA+ patients")
    ax.set_xlabel("Minimum ctDNA VAF%")
    return(ax)  

def make_pie_chart(muts_df, PATH_sample_information, ax):
    diagnosis = muts_df["Diagnosis"].unique()[0]
    annotated = annotate_mutation_status(muts_df, diagnosis, PATH_sample_information, annotate_what = "ctDNA")
    annotated = annotated[annotated["Timepoint"] == "Baseline"].reset_index(drop = True)
    
    pos_counts = annotated[annotated["ctDNA status"] == "Positive"].shape[0]
    neg_counts = annotated[annotated["ctDNA status"] == "Negative"].shape[0]
    
    labels = ['', '']
    sizes = [pos_counts, neg_counts]
    if diagnosis == "Bladder": 
        colors = ['deepskyblue', '#b1ebff']
    else:
        colors = ['orangered', '#ff8d62']
    ax.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, colors=colors, radius = 1.2)
    # ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    return(ax)

def plot_ctDNA_vafs(all_vars_somatic, ax):
    """
    Plot the ctDNA vafs in log from RCC and mUC into the same plot and compare.
    """
    muc_vafs = all_vars_somatic[all_vars_somatic["Diagnosis"] == "Bladder"].reset_index(drop = True)["VAF_t"]
    mrcc_vafs = all_vars_somatic[all_vars_somatic["Diagnosis"] == "Kidney"].reset_index(drop = True)["VAF_t"]
    
    # Log them
    muc_vafs = np.log10(muc_vafs)
    mrcc_vafs = np.log10(mrcc_vafs)
    
    muc_color = "deepskyblue"
    mrcc_color = "orangered"
    
    # Add jitter
    muc_xvalues = np.random.uniform(-0.3, 0.3, muc_vafs.shape[0])
    mrcc_xvalues = np.random.uniform(-0.3, 0.3, mrcc_vafs.shape[0])+1
    
    # Plotting
    ax.scatter(muc_xvalues, muc_vafs, color=muc_color, s = 3)
    ax.scatter(mrcc_xvalues, mrcc_vafs, color=mrcc_color, s = 3)
    
    # Add overlaying boxplot
    muc_boxprops_dict = dict(facecolor=(143/255, 215/255, 239/255, 0), edgecolor='black', linewidth = 1)
    mrcc_boxprops_dict = dict(facecolor=(239/255, 169/255, 143/255, 0), edgecolor='black', linewidth = 1)
    flierprops_dict = dict(marker='o', markersize=5, markeredgecolor='black', linestyle='None')
    whiskerprops_dict =dict(color='black')
    medianprops_dict = dict(color='black')
    capprops_dict = dict(color='black')
    ax.boxplot(muc_vafs, positions = [0], flierprops = flierprops_dict, boxprops = muc_boxprops_dict, medianprops = medianprops_dict, capprops = capprops_dict, widths = 0.4, showfliers = False, patch_artist = True, zorder = 100)
    ax.boxplot(mrcc_vafs, positions = [1], flierprops = flierprops_dict, boxprops = mrcc_boxprops_dict, medianprops = medianprops_dict, capprops = capprops_dict, widths = 0.4, showfliers = False, patch_artist = True, zorder = 100)
    
    # Aes
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylabel("ctDNA VAF%")
    ax.set_ylim(np.log10(0.5), np.log10(100))
    ax.set_yticks([np.log10(0.5), np.log10(1), np.log10(2), np.log10(10), np.log10(50), np.log10(100)])
    ax.set_yticklabels(["0.5", "1", "2", "10", "50", "100"])
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["mUC", "mRCC"])
    ax.set_xlabel("Cancer type")
    
    # Add p values from MWU
    stats, pvalue = mannwhitneyu(muc_vafs, mrcc_vafs)
    return(ax)

source_functions = "/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/visualization/UTILITIES_make_chip_plots.py"

with open(source_functions, 'r') as file:
    script_code = file.read()

exec(script_code)

fig = plt.figure(figsize=(8, 7))
outermost_gs = gridspec.GridSpec(2, 1, height_ratios=[0.5, 1], hspace=0.9, wspace=0.7)  # Adjusted hspace
outer_gs = gridspec.GridSpecFromSubplotSpec(2, 3, width_ratios=[0.35, 1, 1], subplot_spec=outermost_gs[1], wspace=0.7, hspace = 0.7)

# outer_gs = gridspec.GridSpec(3, 3, height_ratios=[1, 1, 1], width_ratios = [1, 1, 1], hspace=0.7, wspace=0.7)  # Adjusted hspace

gs_vafs = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=[1, 1], subplot_spec=outermost_gs[0], wspace=0.3, hspace = 0.3)
ax_vafs = plt.subplot(gs_vafs[0])

ax_bladder_prevalence = plt.subplot(outer_gs[0, 0])
ax_bladder_nmuts = plt.subplot(outer_gs[0, 1])
ax_bladder_genes = plt.subplot(outer_gs[0, 2])

ax_kidney_prevalence = plt.subplot(outer_gs[1, 0])
ax_kidney_nmuts = plt.subplot(outer_gs[1, 1])
ax_kidney_genes = plt.subplot(outer_gs[1, 2])

ax_vafs = plot_ctDNA_vafs(all_vars_somatic, ax_vafs)

color_dict = {"Bladder": "deepskyblue", "Kidney": "orangered"}

ax_kidney_genes, ax_kidney_twin = make_percent_mutated_barchart(base_kidney_somatic, color_dict["Kidney"], ax_kidney_genes, ymax = 45)
ax_bladder_genes, ax_bladder_twin = make_percent_mutated_barchart(base_bladder_somatic, color_dict["Bladder"], ax_bladder_genes, ymax = 60)
# ax_kidney_genes.set_ylim((0, 45))
ax_kidney_genes.set_yticks([0, 10, 20, 30, 40])
ax_kidney_genes.set_yticklabels(["0", "10", "20", "30", "40"])
ax_kidney_twin.set_ylim((0, 50))
ax_kidney_twin.set_yticks([0, 10, 20, 30, 40, 50])
ax_kidney_twin.set_yticklabels(["0", "10", "20", "30", "40", "50"])

# ax_bladder_genes.set_ylim((0, 55))
ax_bladder_genes.set_yticks([0, 10, 20, 30, 40, 50])
ax_bladder_genes.set_yticklabels(["0", "10", "20", "30", "40", "50"])

ax_kidney_nmuts, ax_kidney_nmuts_twin = plot_per_patient_counts(base_kidney_somatic, ax = ax_kidney_nmuts, bar_color = color_dict["Kidney"], superimpose_swarm = True)
ax_bladder_nmuts, ax_bladder_nmuts_twin = plot_per_patient_counts(base_bladder_somatic, ax = ax_bladder_nmuts, bar_color = color_dict["Bladder"], superimpose_swarm = True)

ax_kidney_nmuts.set_xlabel("Number of ctDNA mutations")
ax_kidney_nmuts.set_ylim((0, 32))
ax_kidney_nmuts.set_yticks((0, 8, 16, 24, 32))
ax_kidney_nmuts.set_yticklabels(("0", "8", "16", "24", "32"))
ax_kidney_nmuts.tick_params(axis='x', labelsize=9)  # Set the fontsize for x-tick labels
ax_kidney_nmuts.set_title("", loc = "left")
ax_kidney_nmuts.set_title("", loc = "center")
ax_kidney_nmuts_twin.set_ylim(0, 50)
ax_kidney_nmuts_twin.set_yticks([0, 10, 20, 30, 40, 50])
ax_kidney_nmuts_twin.set_yticklabels(["0", "10", "20", "30", "40", "50"])

ax_bladder_nmuts.set_xlabel("Number of ctDNA mutations")
ax_bladder_nmuts.set_ylim((0, 20))
ax_bladder_nmuts.set_yticks((0, 5, 10, 15, 20))
ax_bladder_nmuts.set_yticklabels(("0", "5", "10", "15", "20"))
ax_bladder_nmuts.set_title("", loc = "left")
ax_bladder_nmuts.set_title("", loc = "center")
ax_bladder_nmuts_twin.set_ylim((0, 100))
ax_bladder_nmuts_twin.set_yticks([0, 20, 40, 60, 80, 100])
ax_bladder_nmuts_twin.set_yticklabels(["0", "20", "40", "60", "80", "100"])


# Prevalence plots
ax_bladder_prevalence = make_min_vaf_ctDNA_prevalence_barchart(base_bladder_somatic, PATH_sample_information, ax_bladder_prevalence, bar_color = color_dict["Bladder"])
ax_kidney_prevalence = make_min_vaf_ctDNA_prevalence_barchart(base_kidney_somatic, PATH_sample_information, ax_kidney_prevalence, bar_color = color_dict["Kidney"])

ax_bladder_prevalence.set_yticks([0, 0.2, 0.40, 0.60, 0.80])
ax_bladder_prevalence.set_yticklabels(["0", "20", "40", "60", "80"])

# ax_kidney_pie = make_pie_chart(base_kidney_somatic, PATH_sample_information, ax_kidney_pie)
# ax_bladder_pie = make_pie_chart(base_bladder_somatic, PATH_sample_information, ax_bladder_pie)

fig.text(0.01, 0.98, 'A', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)

fig.text(0.01, 0.66, 'B', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
fig.text(0.23, 0.66, 'C', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
fig.text(0.61, 0.66, 'D', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)

fig.text(0.01, 0.33, 'E', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
fig.text(0.23, 0.33, 'F', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
fig.text(0.61, 0.33, 'G', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)

outermost_gs.tight_layout(fig)
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/SUPP_ctDNA_mutation_prevalence.png")
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/SUPP_ctDNA_mutation_prevalence.pdf")

# Stat for the figure
mannwhitneyu(base_kidney_somatic["VAF_t"], base_bladder_somatic["VAF_t"])
np.median(base_kidney_somatic["VAF_t"])
np.median(base_bladder_somatic["VAF_t"])

np.percentile(base_kidney_somatic["VAF_t"], 25)
np.percentile(base_kidney_somatic["VAF_t"], 75)

np.percentile(base_bladder_somatic["VAF_t"], 25)
np.percentile(base_bladder_somatic["VAF_t"], 75)