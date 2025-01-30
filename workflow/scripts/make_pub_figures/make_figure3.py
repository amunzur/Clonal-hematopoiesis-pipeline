"""
Makes fragmentomics figures.
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
from lifelines import KaplanMeierFitter
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.cm import get_cmap 
from datetime import datetime
import matplotlib.cm as cm
from lifelines.statistics import logrank_test
from scipy.stats import fisher_exact
from lifelines import KaplanMeierFitter
import matplotlib.ticker as ticker
import upsetplot
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportions_ztest
from scipy.stats import fisher_exact
from scipy.stats import gaussian_kde

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

def compare_to_groundtruth(chip, df_depth):
    """
    Determine how many of the mutations could be detected in the df.
    """
    df_depth = df_depth[["Patient_id", "Sample_name_n", 'Date_collected', 'Timepoint', 'Diagnosis', 
        'Chrom', 'Position', 'Ref', 'Alt', 'Gene', 'Protein_annotation']]
    chip["Date_collected"] = pd.to_datetime(chip["Date_collected"], format='%Y%b%d')
    df_depth["Date_collected"] = pd.to_datetime(df_depth["Date_collected"], format='%Y%b%d')
    merged = chip.merge(df_depth, how = "left", indicator = True)
    merged_df = merged["_merge"].value_counts(normalize = True).reset_index()
    perc_value = merged_df[merged_df["index"] == "both"].reset_index()["_merge"][0]
    return(perc_value)

def run_comparison(chip_df):
    results_dict_1 = {}
    for depth in [1500, 1400, 1300, 1200, 1100, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100]:
        print(f"Working on depth {depth}")
        results_dict_2 = {}
        for min_reads in [2, 3, 4, 5]:
            df_depth = return_chip_df(f"depth_{depth}", min_reads)  # Assuming this function returns a dataframe
            perc_value = compare_to_groundtruth(chip_df, df_depth)  # Your comparison function
            results_dict_2[min_reads] = perc_value
        results_dict_1[depth] = results_dict_2
    return results_dict_1


def plotting_wbc_downsampling(results_dict_1, ax, title, set_ylabel = True):
    """
    Plots the percentage of mutations caught for each depth value.
    """
    df = pd.DataFrame(results_dict_1)
    xticks_raw = df.columns.tolist()
    xticks_fraction = [round(x/2290, 2) for x in xticks_raw]
    xticks_combined = [f'{raw}\n{fraction}' for fraction, raw in zip(xticks_fraction, xticks_raw)]
    
    # Write the WBC depth as fraction of the 
    df_long = df.reset_index().melt(id_vars='index', var_name='Depth', value_name='Value')
    
    colors_dict = {5: "#264653", 4: "#2a9d8f", 3: "#e9c46a", 2: "#f4a261"}
    df_long["color"] = df_long["index"].map(colors_dict)
    
    ax.scatter(df_long["Depth"], df_long["Value"], color=df_long["color"], zorder = 100, s = 30)
    
    for i, group in df_long.groupby("index"):
        color = group["color"].unique()[0]
        ax.plot(group["Depth"], group["Value"], linestyle='--', color = color, linewidth = 1)
        
    ax.invert_xaxis()
    ax.set_xlim(1550, 50)
    ax.set_xticks(range(1500, 0, -100))
    # ticklabels = [f"{tick}X" for tick in ax.get_xticks()]
    ax.set_xticklabels(xticks_combined, ha = "center", fontsize = 7)
    ax.set_xlabel("Downsampled WBC sequencing depth")
    
    ax.set_ylim((0, 1.05))
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax.set_yticklabels(["0", "20", "40", "60", "80", "100"])
    if set_ylabel:
        ax.set_ylabel("% of CH calls identified")
    
    # add horizontal lines through ticks
    for tick in ax.get_yticks():
        ax.hlines(y=tick, xmin=0.0, xmax=1.0, color='b')
    
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(title)
    return(ax)

def plotting_wbc_downsampling_make_legend(ax_wbc_legend):
    """
    Makes legend.
    """
    colors_dict = {"≥5 reads": "#264653", "≥4 reads": "#2a9d8f", "≥3 reads": "#e9c46a", "≥2 reads": "#f4a261"}
    legend_colors = colors_dict.values()
    legend_labels = colors_dict.keys()
    legend_handles = [plt.Line2D([0], [0], marker='o', color=color, label=label, markersize=5, linestyle='') for color, label in zip(legend_colors, legend_labels)]
    ax_wbc_legend.legend(handles=legend_handles, loc="center", frameon=False, fontsize = 8, handletextpad=0.3, ncol = len(legend_colors))
    ax_wbc_legend.spines[["top", "right", "bottom", "left"]].set_visible(False)
    ax_wbc_legend.set_yticks([])
    ax_wbc_legend.set_xticks([])
    return(ax_wbc_legend)

def return_chip_df(depth_name, min_alt_reads):
    """
    Given a depth name, return the relevant chip output file.
    """
    path_chip = os.path.join("/groups/wyattgrp/users/amunzur/pipeline/results/wbc_downsampling", depth_name, f"variant_calling/min_alt_reads_{min_alt_reads}_chip.csv")
    print(f"Loading {path_chip}")
    df = pd.read_csv(path_chip)
    df = df.rename(columns = {"Sample_name": "Sample_name_n"})
    df.loc[df["Gene"] == "U2AF1\\x3bU2AF1L5", "Gene"] = "U2AF1"
    return(df)

#############################
# Fragmentomics functions
def load_fragmentomics_datafiles(path_ch_fr, path_ctDNA_fr):
    """
    Loads the raw dfs.
    """
    path_fragmentomics_functions = "/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/fragmentomics/source_functions.py"
    with open(path_fragmentomics_functions, 'r') as file:
        script_code = file.read()
    
    exec(script_code)
    
    # path_to_df = os.path.join(DIR_fragment_counts, "CH.csv")
    df_ch = MWU_mutated_to_WT_fragments(path_ch_fr, mutated_col_name = "CH fragments", WT_col_name = "WT fragments") # 18 sig
    df_ch['CH fragments'] = df_ch['CH fragments'].apply(ast.literal_eval)
    df_ch['WT fragments'] = df_ch['WT fragments'].apply(ast.literal_eval)
    
    # ctDNA vs WT fragments
    # path_to_df = os.path.join(DIR_fragment_counts, "ctDNA.csv")
    df_ctDNA = MWU_mutated_to_WT_fragments(path_ctDNA_fr, mutated_col_name = "ctDNA fragments", WT_col_name = "WT fragments") # 189 sig
    df_ctDNA['ctDNA fragments'] = df_ctDNA['ctDNA fragments'].apply(ast.literal_eval)
    df_ctDNA['WT fragments'] = df_ctDNA['WT fragments'].apply(ast.literal_eval)
    
    CH_mutant_fragments = list(chain.from_iterable(df_ch['CH fragments']))
    CH_WT_fragments = list(chain.from_iterable(df_ch['WT fragments']))
    ctDNA_mutant_fragments = list(chain.from_iterable(df_ctDNA['ctDNA fragments']))
    ctDNA_WT_fragments = list(chain.from_iterable(df_ctDNA['WT fragments']))
    
    result_dict = {
        "CH_mutant_fragments": CH_mutant_fragments,
        "CH_WT_fragments": CH_WT_fragments,
        "ctDNA_mutant_fragments": ctDNA_mutant_fragments,
        "ctDNA_WT_fragments": ctDNA_WT_fragments
    }
    
    return(result_dict)

def make_kde(result_dict, smoothing = 0.01):
    """
    Generates KDE values from the dict.
    As input, uses the output of the load_fragmentomics_datafiles function.
    """
    CH_mutant_fragments = result_dict["CH_mutant_fragments"]
    CH_WT_fragments = result_dict["CH_WT_fragments"]
    ctDNA_mutant_fragments = result_dict["ctDNA_mutant_fragments"]
    ctDNA_WT_fragments = result_dict["ctDNA_WT_fragments"]
    
    # generate kde values and plot
    kde1 = gaussian_kde(CH_mutant_fragments, bw_method=smoothing)
    kde2 = gaussian_kde(CH_WT_fragments, bw_method=smoothing)
    kde3 = gaussian_kde(ctDNA_mutant_fragments, bw_method=smoothing)
    kde4 = gaussian_kde(ctDNA_WT_fragments, bw_method=smoothing)
    
    x_range = np.linspace(50, 400, 1000)
    
    kde1_values = kde1(x_range)
    kde2_values = kde2(x_range)
    kde3_values = kde3(x_range) 
    kde4_values = kde4(x_range)
    
    kde_dict = {
        "CH_mutant_fragments": kde1_values,
        "CH_WT_fragments": kde2_values,
        "ctDNA_mutant_fragments": kde3_values,
        "ctDNA_WT_fragments": kde4_values,
    }
    
    return(kde_dict)

def plot_fragmentomics_kde(kde_dict, ax1, ax2):
    """
    Plots the KDEs. Uses the output of the make_kde function.
    """    
    KDE_CH_mutant_fragments = kde_dict["CH_mutant_fragments"]
    KDE_CH_WT_fragments = kde_dict["CH_WT_fragments"]
    KDE_ctDNA_mutant_fragments = kde_dict["ctDNA_mutant_fragments"]
    KDE_ctDNA_WT_fragments = kde_dict["ctDNA_WT_fragments"]
    
    # plotting
    x_range = np.linspace(50, 400, 1000)
    ax1.plot(x_range, KDE_CH_WT_fragments, color = "black")
    ax1.plot(x_range, KDE_CH_mutant_fragments, color = "red")
    ax2.plot(x_range, KDE_ctDNA_WT_fragments, color = "black")
    ax2.plot(x_range, KDE_ctDNA_mutant_fragments, color = "forestgreen")
    
    for ax in [ax1, ax2]:
        ax.spines[["right", "top"]].set_visible(False)
        ax.set_xlabel("Fragment length")
        ax.set_xticks([100, 200, 300, 400])
        ax.set_xticklabels(["100", "200", "300", "400"])
    
    ax1.set_title("CH mutations", loc = "center")
    ax2.set_title("ctDNA mutations", loc = "center")
    ax2.set_ylabel("")
    ax1.set_ylabel("Density")
    ax2.tick_params(axis='y', labelleft=False)
    ax2.spines["left"].set_visible(False)
    ax2.tick_params(axis='y', left=False, labelleft=False)
    
    # Add legend to the right ax only
    legend_colors = ["red", "black"]
    legend_labels = ["CH", "WT"]
    legend_handles = [plt.Line2D([0], [0], color=color, label=label, linestyle='-') for color, label in zip(legend_colors, legend_labels)]
    ax1.legend(handles=legend_handles, loc="upper right", frameon=False, fontsize = 8, handletextpad=0.3)
    
    legend_colors = ["forestgreen", "black"]
    legend_labels = ["ctDNA", "WT"]
    legend_handles = [plt.Line2D([0], [0], color=color, label=label, linestyle='-') for color, label in zip(legend_colors, legend_labels)]
    ax2.legend(handles=legend_handles, loc="upper right", frameon=False, fontsize = 8, handletextpad=0.3)
    
    return(ax1, ax2)    

DIR_working = "/groups/wyattgrp/users/amunzur/pipeline"
PATH_sample_information = os.path.join(DIR_working, "resources/sample_lists/sample_information.tsv")
PATH_gene_categories = os.path.join(DIR_working, "resources/panel/chip_panel_gene_categories.tsv")
PATH_treatment_landscape = os.path.join(DIR_working, "resources/clinical_data/bladder/treatment.csv")
PATH_mutation_ctfractions = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/ctDNA_fractions.csv"
figure_dir = os.path.join(DIR_working, "results/figures/pub_figures")
source_functions = os.path.join(DIR_working, "workflow/scripts/visualization/UTILITIES_make_chip_plots.py")
path_gene_exons = os.path.join(DIR_working, "resources/references/panel_genes_exons_refseq.tsv")
sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])

frag_source_functions = "/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/fragmentomics/source_functions.py"

color_dict = {"Bladder": "deepskyblue", "Kidney": "orangered"}

with open(source_functions, 'r') as file:
    script_code = file.read()

exec(script_code)

with open(frag_source_functions, 'r') as file:
    script_code = file.read()

exec(script_code)

# LOAD CHIP DATASETS
all_vars_chip = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/chip.csv"))
all_vars_chip = all_vars_chip[(all_vars_chip["Dependent"] == False) & (all_vars_chip["Timepoint"] == "Baseline")]
base_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
prog_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)
prog_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)

# LOAD SOMATIC DATASETS
all_vars_somatic = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/somatic.csv"))
all_vars_somatic = all_vars_somatic[(all_vars_somatic["Dependent"] == False) & (all_vars_somatic["Timepoint"] == "Baseline")]
all_vars_somatic = all_vars_somatic[~all_vars_somatic["Patient_id"].isin(["20-313", "21-184", "21-430"])] # exclude some samples due to oxidative damage
base_kidney_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Kidney")].reset_index(drop = True)
prog_kidney_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "During treatment") & (all_vars_somatic["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Bladder")].reset_index(drop = True)
prog_bladder_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "During treatment") & (all_vars_somatic["Diagnosis"] == "Bladder")].reset_index(drop = True)

baseline = pd.concat([base_kidney_chip, base_bladder_chip]).reset_index(drop = True)
baseline_somatic = pd.concat([base_kidney_somatic, base_bladder_somatic]).reset_index(drop = True)

fig = plt.figure(figsize=(8, 11))
outer_gs = gridspec.GridSpec(4, 1, height_ratios=[1, 1, 1.5, 1], hspace = 0.2, wspace = 1) # outer most gs with 3 rows
inner_gs0 = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios = [1, 2], wspace = 0.2, hspace = 0.2, subplot_spec=outer_gs[0])
inner_gs1 = gridspec.GridSpecFromSubplotSpec(1, 3, width_ratios=[1, 1, 1], wspace=0.05, hspace = 0.05, subplot_spec=outer_gs[1])
inner_gs2 = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios = [1, 0.05], subplot_spec=outer_gs[2], wspace=0.1, hspace = 1.2)
inner_inner_gs2 = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios = [1, 1], subplot_spec=inner_gs2[0], wspace=0.2, hspace = 0.2)
inner_gs3 = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios = [0.6, 1], subplot_spec=outer_gs[3], wspace=0.1, hspace = 0.1)

# FIGURE 3A. Overlap of VAF distribution in cfDNA and CH mutations
ax0 = plt.subplot(inner_gs0[0])
path_sample_information = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"
ax0 = plot_ctDNA_and_CH_VAF_box(ch_df = baseline, ctdna_df = baseline_somatic, PATH_sample_information = PATH_sample_information, ax = ax0, ax_title = "", hide_xticks = False)
# ax0 = plot_ctDNA_and_CH_VAF_split_violin(PATH_sample_information, baseline, baseline_somatic, ax0, ax_title = "", hide_xticks = False, show_legend = True)

# FIGURE 3B. Fragment distribution
fragmentomics_gs = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios = [1, 1], wspace = 0.05, hspace = 0.05, subplot_spec=inner_gs0[1])
path_ch_fr = "/groups/wyattgrp/users/amunzur/pipeline/results/fragmentomics/CH.csv"
path_ctDNA_fr = "/groups/wyattgrp/users/amunzur/pipeline/results/fragmentomics/ctDNA.csv"
# result_dict = load_fragmentomics_datafiles(path_ch_fr, path_ctDNA_fr)
# kde_dict = make_kde(result_dict)
ax1 = plt.subplot(fragmentomics_gs[0])
ax2 = plt.subplot(fragmentomics_gs[1], sharey = ax1)
ax1, ax2 = plot_fragmentomics_kde(kde_dict, ax1, ax2)

# FIGURE 2C - Fragmentomics scatter
ax3 = plt.subplot(inner_gs1[1])
ax4 = plt.subplot(inner_gs1[2], sharey = ax3, sharex = ax3)
# df_ctDNA = compare_mutated_to_WT(df_path = "/groups/wyattgrp/users/amunzur/pipeline/results/fragmentomics/ctDNA.csv")
# df_CH = compare_mutated_to_WT(df_path = "/groups/wyattgrp/users/amunzur/pipeline/results/fragmentomics/CH.csv")
ax3, ax3_inset = plot_mutated_median_vs_wt_median(df_CH, ax3, title = "CH mutations", set_y_label = True, add_legend = False, add_pie = True)
ax4, ax4_inset = plot_mutated_median_vs_wt_median(df_ctDNA, ax4, title = "ctDNA mutations", set_y_label = False, add_legend = False, add_pie = True)
ax4.tick_params(axis='y', left=False, labelleft = False)
ax4.spines["left"].set_visible(False)
ax4.set_ylim((150, 325))
ax4.set_yticks([150, 200, 250, 300])
ax4.set_xticks([150, 250, 350])
ax4.tick_params(axis='y', labelleft=False)

# # The relationship between WBC and cfDNA vaf correlation and ctDNA fraction
# gs_corr = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios = [1, 1], subplot_spec=inner_gs2[2])
# ax_corr_dot = plt.subplot(gs_corr[0])
# ax_corr_box = plt.subplot(gs_corr[1], sharex = ax_corr_dot)
# ax_corr_dot, ax_corr_box = plot_ctDNA_bins_and_VAF_correlations(PATH_sample_information, ch_df = baseline, ctdna_df = baseline_somatic, ax1 = ax_corr_dot, ax2 = ax_corr_box)

# WBC downsampling
ax5 = plt.subplot(inner_inner_gs2[0])
ax6 = plt.subplot(inner_inner_gs2[1])
ax_wbc_legend = plt.subplot(inner_gs2[1])

chip = baseline[["Patient_id", "Sample_name_n", 'Date_collected', 'Timepoint', 'Diagnosis', 'Chrom', 'Position', 'Ref', 'Alt', 'Gene', 'AAchange', 'VAF_n']]
chip_1_perc = baseline[baseline["VAF_n"] > 1]

# results_dict_chip = run_comparison(chip) # Run comparison for chip
# results_dict_chip_1_perc = run_comparison(chip_1_perc) # Run comparison for chip_1_perc

ax5 = plotting_wbc_downsampling(results_dict_chip, ax5, title = "All CH variants")
ax6 = plotting_wbc_downsampling(results_dict_chip_1_perc, ax6, title = "CH variants with VAF > 1%", set_ylabel = False)
ax_wbc_legend = plotting_wbc_downsampling_make_legend(ax_wbc_legend)

# Filtering the variants step by step
ax7 = plt.subplot(inner_gs3[1]) # cartoon

industry_results =  filter_like_industry(all_vars_chip, all_vars_somatic)
industry_results = industry_results[industry_results["n total_ctDNA"] > 4].reset_index(drop = True)
ax7 = make_variant_retention_barchart(industry_results, ax = ax7, color_dict = {
        "True positive": "forestgreen",
        "False negative": (0.13, 0.55, 0.13, 0.3),  # forestgreen with alpha 0.3
        "False positive": "red",
        "True negative": (1.0, 0, 0, 0.3)  # red with alpha 0.3
})

# Now add the 6 piecharts
gs_pies = gridspec.GridSpecFromSubplotSpec(3, 2, width_ratios = [1, 1], height_ratios = [1, 1, 1], wspace = 0.1, hspace = 0.1, subplot_spec=inner_gs3[0])
add_sensitivity_specificity_pies(gs_pies, all_vars_chip, all_vars_somatic)




fig.text(0.01, 0.99, 'a', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
fig.text(0.35, 0.99, 'b', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
fig.text(0.01, 0.69, 'c', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
fig.text(0.01, 0.38, 'd', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)

outer_gs.tight_layout(fig)
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/figure3_panel.png")
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/figure3_panel.pdf")