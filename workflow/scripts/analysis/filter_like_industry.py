"""
How many false calls do we have when we try to filter out the data like the industry does?
"""

import pandas as pd
import numpy as np
import os 
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def clean_up_df(muts_df):
    """
    Subsetting the df to the relevant cols before analysis.
    """
    muts_df = muts_df[(muts_df["Dependent"] == False) & (muts_df["Timepoint"] == "Baseline")].reset_index(drop = True)
    
    if "Status" in muts_df.columns:
        variant_type = "ctDNA"
    else:
        variant_type = "CH"
    
    muts_df_subsetted = muts_df[["Patient_id", "Gene", "VAF_t", "Consequence"]]
    muts_df_subsetted["Variant_type"] = variant_type
    
    return(muts_df_subsetted)

def calculate_sensitivity(filtered_muts_df, ctdna_vars, chip_vars):
    # From ground truth
    n_chip_total = chip_vars.shape[0]
    n_ctdna_total = ctdna_vars.shape[0]
    
    # From filtered df
    n_ch_called = filtered_muts_df[filtered_muts_df["Variant_type"] == "CH"].shape[0]
    n_ctdna_called = filtered_muts_df[filtered_muts_df["Variant_type"] == "ctDNA"].shape[0]
    
    # Proportion of true ctDNA variants recovered
    sensitivity = round((n_ctdna_called/n_ctdna_total)*100, 1)
    
    # Proportion of all recovered variants that are genuinely from ctDNA
    specificity = round(n_ctdna_called/(n_ch_called+n_ctdna_called)*100, 1)
    
    mydict = {"Sensitivity": sensitivity, "Specificity": specificity}
    
    return(mydict)

def calculate_percentage_CH_kept(filtered_muts_df, chip_vars):
    """
    Calculates the fraction of the CH mutations kept in a set of curated mutations.
    """
    n_ch_total = chip_vars.shape[0]
    n_ch_df = filtered_muts_df[filtered_muts_df["Variant_type"] == "CH"].shape[0]
    frac = round((n_ch_df/n_ch_total)*100, 1)
    
    return(frac)

def calculate_retention_fraction(muts_df_filtered, muts_df_total, gene, mut_type):
    """
    Given a gene name calculates the fraction of ctDNA mutations we correctly identified.
    Provide CH or ctDNA for mut_type.
    """
    muts_df_filtered_subset = muts_df_filtered[muts_df_filtered["Variant_type"] == mut_type]
    muts_df_total_subset = muts_df_total[muts_df_total["Variant_type"] == mut_type]
    
    n_total_gene = muts_df_total_subset[muts_df_total_subset["Gene"] == gene].shape[0]
    n_identified = muts_df_filtered_subset[muts_df_filtered_subset["Gene"] == gene].shape[0]
    
    if n_total_gene == 0:
        identification_fraction = np.nan
    else:
        identification_fraction = round(n_identified/n_total_gene, 3)
    
    result_dict = {"Gene": gene, "n total": n_total_gene, "n retained": n_identified, "retention fraction": identification_fraction}
    
    return(result_dict)

def make_stacked_bar(df, ax):
    """
    For each gene shows the fraction of correctly identified variants vs the fraction missed
    """
    ax.bar(df.index, df['identification_fraction'], color='black')
    ax.bar(df.index, df['missed_fraction'], bottom = df["identification_fraction"], color='grey')
    ax.set_xticks(df.index)
    ax.set_xticklabels(df["Gene"], rotation = 90)
    
    # Annotate numbers
    for xpos, row in df.iterrows():
        ax.text(xpos, row["identification_fraction"]-0.05, row["n identified"], ha='center', va='center', color = "white", fontsize = 5)
        ax.text(xpos, 0.95, row["n_missed"], ha='center', va='center', color = "black", fontsize = 5)
    
    # Add legend
    legend_colors = ["black", "gray"]
    legend_labels = ["Identified", "missed"]
    legend_handles = [plt.Line2D([0], [0], marker='s', color=color, label=label, markersize=3, linestyle='') for color, label in zip(legend_colors, legend_labels)]
    ax.legend(handles=legend_handles, loc="lower right", frameon=True, fontsize = 8, handletextpad=0.1)
    
    return (ax)

def make_variant_retention_barchart(results, ax, color_dict):
    """
    For each gene shows the fraction and number of variants retained and missed, for both CH and ctDNA variants
    """
    # Plot ctDNA stats
    bar_height = 0.8
    ax.barh(results.index, results["n retained_ctDNA"], height=bar_height, color=color_dict["retention_color"], alpha=0.8)
    ax.barh(results.index, results["n missed_ctDNA"], left = results["n retained_ctDNA"], height=bar_height, color=color_dict["missed_color"], alpha=0.8)
    
    # Plot CH stats
    ax.barh(results.index, results["n retained_CH"], left = 70, height=bar_height, color=color_dict["retention_color"], alpha=0.8)
    ax.barh(results.index, results["n missed_ctDNA"], left = 70+results["n retained_CH"], height=bar_height, color=color_dict["missed_color"], alpha=0.8)
    
    # Aes
    ax.axvline(70, color='black', linestyle='--', linewidth = 0.25)
    ax.set_yticks(results.index)
    ax.set_yticklabels(results["Gene"], fontsize = 8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylim((results.index.min() - bar_height*0.75, results.index.max() + bar_height*0.75))    
    ax.set_xticks(0, 20, 60, 80, )
    return(ax)

path_ctDNA_muts = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/somatic.csv"
path_chip_muts = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/chip.csv"

all_vars_ctDNA = pd.read_csv(path_ctDNA_muts)
ctdna_vars = clean_up_df(all_vars_ctDNA)

all_vars_chip = pd.read_csv(path_chip_muts)
chip_vars = clean_up_df(all_vars_chip)

combined = pd.concat([ctdna_vars, chip_vars]).reset_index(drop = True)

# Now filtering based on a few criteria.
# 1. Excluding DTA genes
combined_filtered = combined[~combined["Gene"].isin(["DNMT3A", "TET2", "ASXL1"])]

# 2. Exclude all variants less than 1% in VAF
combined_filtered = combined_filtered[combined_filtered["VAF_t"] > 1]

# 3. Excluding variants that are <25% of the sample's maximum VAF (i.e. a surrogate for tumor fraction), i.e. highly "subclonal" variants.
# Compute sample maximum VAF
subclonal_df = combined.groupby("Patient_id")["VAF_t"].max().reset_index()
subclonal_df["Subclonal VAF threshold"] = subclonal_df["VAF_t"]*0.25
del subclonal_df["VAF_t"]
combined_filtered = combined_filtered.merge(subclonal_df, how = "left")
combined_filtered = combined_filtered[combined_filtered["VAF_t"] > combined_filtered["Subclonal VAF threshold"]]

# Combining results
ctDNA_results_list = []
CH_results_list = []
for gene in ctdna_vars["Gene"].unique(): 
    ctDNA_retention = calculate_retention_fraction(combined_filtered, combined, gene = gene, mut_type = "ctDNA")
    CH_retention = calculate_retention_fraction(combined_filtered, combined, gene = gene, mut_type = "CH")
    
    # Append to list
    ctDNA_results_list.append(ctDNA_retention)
    CH_results_list.append(CH_retention)

# Make a ctDNA df
ctDNA_df = pd.DataFrame(ctDNA_results_list)
ctDNA_df = ctDNA_df.sort_values(by = "n total").reset_index(drop = True)
ctDNA_df["missed fraction"] = 1 - ctDNA_df["retention fraction"]
ctDNA_df["n missed"] = ctDNA_df["n total"] - ctDNA_df["n retained"]

# Make a CH df
CH_df = pd.DataFrame(CH_results_list)
CH_df = CH_df.sort_values(by = "n total").reset_index(drop = True)
CH_df["missed fraction"] = 1 - CH_df["retention fraction"]
CH_df["n missed"] = CH_df["n total"] - CH_df["n retained"]

results = ctDNA_df.merge(CH_df, on = "Gene", suffixes = ["_ctDNA", "_CH"])

# Specificity and sensiticity calculations
df = results[["n retained_ctDNA", "n retained_CH"]].sum().reset_index()
n_ctdna_called = df[df["index"] == "n retained_ctDNA"][0][0]
n_ch_called = df[df["index"] == "n retained_CH"][0][1]
n_ctdna_total = 564
n_chip_total = 607
sensitivity = round((n_ctdna_called/n_ctdna_total)*100, 1)
specificity = round((n_ctdna_called/(564+607))*100, 1)


# Plotting
fig, ax = plt.subplots(figsize=(8, 3))
ax = make_stacked_bar(df, ax)
fig.tight_layout()
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/filter_like_industry.png")
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/filter_like_industry.pdf")

fig, ax = plt.subplots(figsize=(3, 5))

ax = make_variant_retention_barchart(results, ax, color_dict = {"retention_color": "#6B3724", "missed_color": "#CC8166"})
fig.tight_layout()
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/filter_like_industry_retention.pdf")

# STRATEGY 2. APPLY THE FILTERING TECHNIQUES INDIVIDUALLY.
# For each approach we calculate sensitivity and specificity.
# sensitivity: Proportion of true ctDNA variants we detected.
# specificity: Proportion of all recovered variants that are genuinely from ctDNA

ctdna_vars = clean_up_df(all_vars_ctDNA)
chip_vars = clean_up_df(all_vars_chip)

combined = pd.concat([ctdna_vars, chip_vars]).reset_index(drop = True)

# Now filtering based on a few criteria.
# 1. Excluding DTA genes
dta_filtered = combined[~combined["Gene"].isin(["DNMT3A", "TET2", "ASXL1"])]
dta_filtering_dict = calculate_sensitivity(dta_filtered, ctdna_vars, chip_vars)
dta_ch = calculate_percentage_CH_kept(dta_filtered, chip_vars)

# 2. Exclude all variants less than 1% in VAF
vaf_filtered = combined[combined["VAF_t"] > 1]
vaf_filtering_dict = calculate_sensitivity(vaf_filtered, ctdna_vars, chip_vars)
vaf_ch =calculate_percentage_CH_kept(vaf_filtered, chip_vars)
    
# 3. Excluding variants that are <25% of the sample's maximum VAF (i.e. a surrogate for tumor fraction), i.e. highly "subclonal" variants.
# Determine patients with more than one mutation, this will be applied to those only
n_muts_df = combined.value_counts("Patient_id").reset_index().rename(columns = {0: "n_muts"})
max_vaf_df = combined.groupby("Patient_id")["VAF_t"].max().reset_index()
max_vaf_df["Subclonal VAF threshold"] = max_vaf_df["VAF_t"]*0.25
del max_vaf_df["VAF_t"]
max_vaf_df = max_vaf_df.merge(n_muts_df)

combined_merged = combined.merge(max_vaf_df, how = "left")
combined_filtered = combined_merged[
    (combined_merged['n_muts'] > 1) & (combined_merged['VAF_t'] > combined_merged['Subclonal VAF threshold']) |
    (combined_merged['n_muts'] == 1)
]
subclonal_dict = calculate_sensitivity(combined_filtered, ctdna_vars, chip_vars)
subclonal_ch = calculate_percentage_CH_kept(combined_filtered, chip_vars)

# and now plot these as barchart
def make_sensitivity_specificity_bar(subclonal_dict, vaf_filtering_dict, dta_filtering_dict, ax):
    
    sensitivity_color = "black"
    specificity_color = "gray"
    for i, mydict in enumerate([subclonal_dict, vaf_filtering_dict, dta_filtering_dict]):
        ax.bar(i-0.1, mydict["Sensitivity"], width=0.12, color=sensitivity_color)
        ax.bar(i+0.1, mydict["Specificity"], width=0.12, color=specificity_color)
    
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["subclonal", "vaf", "dta"])
    ax.spines[["top", "right", "bottom"]].set_visible(False)
    ax.set_yticks((0, 50, 100))
    return(ax)

def make_sensitivity_specificity_bar(subclonal_dict, vaf_filtering_dict, dta_filtering_dict, ax):
    
    sensitivity_color = "black"
    specificity_color = "gray"
    for i, mydict in enumerate([subclonal_dict, vaf_filtering_dict, dta_filtering_dict]):
        ax.bar(i-0.1, mydict["Sensitivity"], width=0.12, color=sensitivity_color)
        ax.bar(i+0.1, mydict["Specificity"], width=0.12, color=specificity_color)
    
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["subclonal", "vaf", "dta"])
    ax.spines[["top", "right", "bottom"]].set_visible(False)
    ax.set_yticks((0, 50, 100))
    return(ax)

def make_piecharts(subclonal_dict, vaf_filtering_dict, dta_filtering_dict, ax):
    # Define the labels and colors for the pie charts
    labels = ['Sensitivity', 'Specificity']
    colors = ['black', 'gray']
    
    # Dictionaries to loop through
    dicts = [subclonal_dict, vaf_filtering_dict, dta_filtering_dict]
    
    # Create pie charts for each dictionary
    for i, mydict in enumerate(dicts):
        # Sensitivity pie chart
        ax[i, 0].pie([mydict["Sensitivity"], 100 - mydict["Sensitivity"]], labels=[labels[0], ''], colors=[colors[0], 'white'], startangle=90, autopct='%1.1f%%',  wedgeprops={'edgecolor': 'black'})
        ax[i, 0].set_title(f"{['Subclonal', 'VAF', 'DTA'][i]} Sensitivity")
        
        # Specificity pie chart
        ax[i, 1].pie([mydict["Specificity"], 100 - mydict["Specificity"]], labels=[labels[1], ''], colors=[colors[1], 'white'], startangle=90, autopct='%1.1f%%', wedgeprops={'edgecolor': 'black'})
        ax[i, 1].set_title(f"{['Subclonal', 'VAF', 'DTA'][i]} Specificity")

    
    

fig, ax = plt.subplots(figsize=(2, 2))
ax = make_sensitivity_specificity_bar(subclonal_dict, vaf_filtering_dict, dta_filtering_dict,ax)
fig.tight_layout()
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/filter_like_industry.png")
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/filter_like_industry.pdf")

# Create a figure with gridspec (3 rows, 2 columns)
fig, ax = plt.subplots(3, 2, figsize=(5, 5), gridspec_kw={'wspace': 0.5})
# Call the function with the axis array
make_piecharts(subclonal_dict, vaf_filtering_dict, dta_filtering_dict, ax)
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/sensitivity_specificity_pies.pdf")
