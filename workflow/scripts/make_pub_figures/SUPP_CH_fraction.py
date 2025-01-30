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
from adjustText import adjust_text
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

def vafs_of_mutations_same_gene_same_person(all_vars_chip, all_vars_ctDNA, ax_ch, ax_ctdna):
    """
    Determines patients with both CH and ctDNA mutations in the same gene and plots their VAF distribution.
    """
    ch = all_vars_chip[["Patient_id", "Diagnosis", "Gene", "VAF_n"]]
    ctdna = all_vars_ctDNA[["Patient_id", "Diagnosis", "Gene", "VAF_t"]]
    merged = ch.merge(ctdna)    
    
    merged["VAF_t"] = np.log10(merged["VAF_t"])
    merged["VAF_n"] = np.log10(merged["VAF_n"])
    
    gene_order = ["TP53", "CHEK2", "ATM", "KMT2D", "SETD2", "STAG2", "ARID1A", "DNMT3A", "TET2"]
    merged['Gene'] = pd.Categorical(merged['Gene'], categories=gene_order, ordered=True)
    merged = merged.sort_values('Gene').reset_index(drop = True)
    # sorted_groups = merged.groupby("Gene").apply(lambda x: x["VAF_n"].max()).sort_values(ascending=False).index
    # Generate the x ticks list
    merged["xtick_mask"] = merged["Gene"].astype(str) + "_" + merged["Patient_id"]
    merged["xtick"] = pd.factorize(merged['xtick_mask'])[0]
    
    for i, row in merged.iterrows():
        patient_id = row["Patient_id"]
        xpos = row["xtick"]
        VAF_n = row["VAF_n"]
        VAF_t = row["VAF_t"] * -1
        ax_ch.scatter(xpos, VAF_n, color="red", s=6, zorder = 100)
        ax_ctdna.scatter(xpos, VAF_t, color="forestgreen", s=6, zorder = 100)
    
    ax_ch.set_ylim((np.log10(0.2), np.log10(45)))
    ax_ctdna.set_ylim((-np.log10(45), -np.log10(0.25)))
    
    xticklabel_list = merged["xtick_mask"].apply(lambda x: x.split("_")[1]).unique()
    
    # Add the alternating rectangles to indicate genes
    # color = "gray"
    # for gene, group in merged.groupby("Gene"):
    #     xpos_min = group["xtick"].min()
    #     xpos_max = group["xtick"].max()
    #     width_rectangle = xpos_max-xpos_min
    #     if width_rectangle == 0:
    #         width_rectangle = 1
    #     print(f"{gene} - {width_rectangle}")
    #     if color is "gray": 
    #         color = "black"
    #     else:
    #         color = "gray"
    #     ax_ch.bar(xpos_min, np.log10(45)-np.log10(0.25), bottom = np.log10(0.25), align = "edge", width = width_rectangle, color = color)
    #     ax_ctdna.bar(xpos_min, -np.log10(45)+np.log10(0.25), bottom = -np.log10(0.25), align = "edge", width = width_rectangle, color = color)
    
    ax_ch.spines[["top", "right"]].set_visible(False)
    ax_ch.set_xticks(np.arange(0, merged["xtick"].max()+1))
    ax_ch.set_xticklabels(xticklabel_list, rotation = 90)
    ax_ch.set_yticks([np.log10(0.25), np.log10(1), np.log10(2), np.log10(10), np.log10(45)])
    ax_ch.set_yticklabels(["0.25", "1", "2", "10", "45"])
    ax_ch.tick_params(axis='x', bottom=False)
    ax_ch.set_ylabel("WBC VAF% of\nCH mutations")   
    
    ax_ctdna.spines[["bottom", "right"]].set_visible(False)
    ax_ctdna.set_yticklabels(["0.25", "1", "2", "10", "45"])
    ax_ctdna.set_yticks([-np.log10(0.25), -np.log10(1), -np.log10(2), -np.log10(10), -np.log10(45)])
    ax_ctdna.tick_params(axis='x', bottom=False, labelbottom=False)
    ax_ctdna.set_ylabel("cfDNA VAF% of\ntumor mutations")
    
    return(ax_ch, ax_ctdna)

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
PATH_gene_categories = os.path.join(DIR_working, "resources/panel/chip_panel_gene_categories.tsv")
PATH_bladder_clinical_data = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/bladder/clinical_data.csv"
PATH_kidney_clinical_data = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/RCC clinical data - mRCC clinical Data.csv"
PATH_treatment_landscape = os.path.join(DIR_working, "resources/clinical_data/bladder/treatment.csv")
figure_dir = os.path.join(DIR_working, "results/figures/pub_figures")
source_functions = os.path.join(DIR_working, "workflow/scripts/visualization/UTILITIES_make_chip_plots.py")
sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])

color_dict = {"Bladder": "deepskyblue", "Kidney": "orangered"}

with open(source_functions, 'r') as file:
    script_code = file.read()

exec(script_code)

# LOAD CHIP DATASETS
# ctDNA mutations
all_vars_ctDNA = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/somatic.csv"))
all_vars_ctDNA = all_vars_ctDNA[(all_vars_ctDNA["Dependent"] == False) & (all_vars_ctDNA["Timepoint"] == "Baseline")]

# CH mutations
all_vars_chip = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/chip.csv"))
all_vars_chip = all_vars_chip[(all_vars_chip["Dependent"] == False) & (all_vars_chip["Timepoint"] == "Baseline")]
all_vars_chip_vaf = all_vars_chip[all_vars_chip["VAF_n"]>=1]

fig = plt.figure(figsize=(8, 6))
outer_gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1], hspace = 0.8) # outer most gs with 3 rows
inner_gs0 = gridspec.GridSpecFromSubplotSpec(3, 1, height_ratios = [0.008, 1, 1], hspace = 1.1, subplot_spec=outer_gs[0])
inner_inner_gs0 = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios = [1, 0.3], wspace=0.2, hspace = 0.2, subplot_spec=inner_gs0[2])

ax0 = fig.add_subplot(inner_gs0[1])
ax1 = fig.add_subplot(inner_inner_gs0[0], sharex=ax0)
ax2 = fig.add_subplot(inner_inner_gs0[1], sharex=ax0)

ax0, ax1, ax2 = plot_fraction_of_CH_calls_with_bars(df_somatic_main = all_vars_ctDNA, df_ch_main = all_vars_chip_vaf, ax0 = ax0, ax1 = ax1, ax2 = ax2, stacked = False, yval_max = -110) # plot_fraction_of_CH_calls was the old function, still available in the source functions.

# Gene vafs plot
inner_gs_vafs = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios = [1, 1], wspace = 0.2, subplot_spec=outer_gs[1])
ax3 = fig.add_subplot(inner_gs_vafs[1])
baseline_filtered = all_vars_chip[~all_vars_chip["Gene"].isin(["DNMT3A", "TET2"])]
baseline_somatic_filtered = all_vars_ctDNA[~all_vars_ctDNA["Gene"].isin(["DNMT3A", "TET2"])]
ax3 = plot_ctDNA_and_CH_VAF_box(ch_df = baseline_filtered, ctdna_df = baseline_somatic_filtered, PATH_sample_information = PATH_sample_information, ax = ax3, ax_title = "", hide_xticks = False)

# ks_2samp(baseline_filtered["VAF_t"], baseline_somatic_filtered["VAF_t"])

# Patients with CH and ctDNA mutations in the same genes
inner_gs_vafs_same_gene_vafs = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios = [1, 1], hspace = 0.6, subplot_spec=inner_gs_vafs[0])
ax_ch = fig.add_subplot(inner_gs_vafs_same_gene_vafs[0])
ax_ctdna = fig.add_subplot(inner_gs_vafs_same_gene_vafs[1], sharex = ax_ch)
ax_ch, ax_ctdna = vafs_of_mutations_same_gene_same_person(all_vars_chip, all_vars_ctDNA, ax_ch, ax_ctdna)

fig.text(0.01, 0.99, 'A', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
fig.text(0.01, 0.55, 'B', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
fig.text(0.5, 0.55, 'C', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)

outer_gs.tight_layout(fig)
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/SUPP_CH_fraction.png")
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/SUPP_CH_fraction.pdf")


