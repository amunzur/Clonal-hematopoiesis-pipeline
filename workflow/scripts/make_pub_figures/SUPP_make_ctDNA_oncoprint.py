import os
import pandas as pd
import numpy as np
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
from scipy import stats
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.font_manager import FontProperties
from matplotlib.colors import Normalize, to_rgba
from matplotlib.patches import Rectangle
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import seaborn as sns
import matplotlib.colors as mcolors

mut_dict = {
    "frameshift indel": '#FFC907',
    "nonframeshift indel": '#a9a9a9',
    "missense": '#79B443',
    "stopgain": '#BD4398',
    "splicing": "darkorange"}

mpl.rcParams['font.size'] = 8
mpl.rcParams['text.color'] = 'k'
mpl.rcParams['legend.fontsize'] = 10
mpl.rcParams['legend.handletextpad'] = '0.8'
mpl.rcParams['legend.labelspacing'] = '0.4'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['axes.linewidth'] = 0.5
mpl.rcParams['xtick.labelsize'] = 7
mpl.rcParams['ytick.labelsize'] = 7
mpl.rcParams['axes.labelsize'] = 10
mpl.rcParams['xtick.major.size'] = 3
mpl.rcParams['ytick.major.size'] = 3

mpl.rcParams['legend.fontsize'] = 6
mpl.rcParams['legend.handletextpad'] = '0.7'
mpl.rcParams['legend.labelspacing'] = '0.2'
plt.rcParams['legend.handlelength'] = 0.4
plt.rcParams['legend.handleheight'] = 0.70

# LOAD DATASETS
DIR_working = "/groups/wyattgrp/users/amunzur/pipeline"
PATH_sample_information = os.path.join(DIR_working, "resources/sample_lists/sample_information.tsv")
PATH_mutation_ctfractions = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/ctDNA_fractions.csv"
PATH_bladder_clinical_data = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/bladder/clinical_data.csv"
PATH_kidney_clinical_data = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/RCC clinical data - mRCC clinical Data.csv"

# ctDNA mutations
all_vars_ctDNA = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/somatic.csv"))
all_vars_ctDNA = all_vars_ctDNA[all_vars_ctDNA["Dependent"] == False]
all_vars_ctDNA["Consequence"] = all_vars_ctDNA["Consequence"].replace({"Frameshift deletion":"Frameshift indel", "Frameshift insertion":"Frameshift indel", "Nonframeshift deletion":"Nonframeshift indel", "Nonframeshift insertion":"Nonframeshift indel"})
bladder_ctDNA_main = all_vars_ctDNA[(all_vars_ctDNA["Diagnosis"] == "Bladder") & (all_vars_ctDNA["Timepoint"] == "Baseline")]
kidney_ctDNA_main = all_vars_ctDNA[(all_vars_ctDNA["Diagnosis"] == "Kidney") & (all_vars_ctDNA["Timepoint"] == "Baseline")]

# CH mutations
all_vars_chip = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/chip.csv"))
all_vars_chip = all_vars_chip[all_vars_chip["Dependent"] == False]
all_vars_chip["Consequence"] = all_vars_chip["Consequence"].replace({"Frameshift deletion":"Frameshift indel", "Frameshift insertion":"Frameshift indel", "Nonframeshift deletion":"Nonframeshift indel", "Nonframeshift insertion":"Nonframeshift indel"})
bladder_chip_main = all_vars_chip[(all_vars_chip["Diagnosis"] == "Bladder") & (all_vars_chip["Timepoint"] == "Baseline")]
kidney_chip_main = all_vars_chip[(all_vars_chip["Diagnosis"] == "Kidney") & (all_vars_chip["Timepoint"] == "Baseline")]

# source functions
source_functions = "/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/make_pub_figures/make_oncoprint_functions.py"
with open(source_functions, 'r') as file:
    script_code = file.read()

exec(script_code)

# OPs for the ctDNA mutations
##############################################################
# BLADDER
# PATH_gene_categories_bladder_somatic = "/groups/wyattgrp/users/amunzur/pipeline/resources/panel/for_bladder_somatic_op.tsv"
bladder_OP_genes_dict = {
    "DNMT3A":17, "TET2":16, "ASXL1":15, "JAK2":14, "ATM":13, "BRCA2":12, "PPM1D":11, "CHEK2":10,
    "TP53":9, "ARID1A":8,"ERBB2":7,"STAG2":6,"RB1":5, "FGFR3":4,"PIK3CA":3,"KMT2D":2, "KDM6A":1, "TERT":0
    }

fig_width = 8.5
fig_height = 11
tc_height = 2
mutcounts_height = 2
sex_height = 0.5
age_height = 0.5
main_height = len(bladder_OP_genes_dict.keys())*2

fig = plt.figure(figsize = (fig_width, fig_height))
gs = mpl.gridspec.GridSpec(ncols = 2, nrows = 7, height_ratios = [tc_height, mutcounts_height, sex_height, age_height, main_height, 2, 3], width_ratios = [1, 0.1], hspace = 0, wspace = 0)

# Add and space out subplots
ax_tc = fig.add_subplot(gs[0, 0]) 
ax_mut_counts_patient = fig.add_subplot(gs[1, 0], sharex = ax_tc) 
ax_main = fig.add_subplot(gs[4, 0], sharex = ax_tc)
ax_mut_counts = fig.add_subplot(gs[4, 1], sharey = ax_main) 
ax_sex = fig.add_subplot(gs[2, 0], sharex = ax_main) 
ax_age = fig.add_subplot(gs[3, 0], sharex = ax_main) 

bar_height = 0.9
bar_width = 0.8

# Generate dfs
[bladder_ctDNA, gene_order_df, chip_fr, CTDNA_mut_counts_df, CTDNA_vaf_df, CTDNA_pt_counts, pts_ctdna_NO_nan, pts_ctdna_nan, DF_Samples_enumerated, mut_dict_kidney, sex_df, age_df, combined_age_df] = prepare_datasets_for_plotting(bladder_ctDNA_main, mut_dict, bladder_OP_genes_dict, PATH_mutation_ctfractions, PATH_bladder_clinical_data, PATH_kidney_clinical_data)
[bladder_ch, gene_order_df, chip_fr, CH_mut_counts_df, CH_vaf_df, CH_pt_counts, pts_ctdna_NO_nan, CH_pts_ctdna_nan, DF_Samples_enumerated, mut_dict_kidney, sex_df, age_df, combined_age_df] = prepare_datasets_for_plotting(bladder_chip_main, mut_dict, bladder_OP_genes_dict, PATH_mutation_ctfractions, PATH_bladder_clinical_data, PATH_kidney_clinical_data)

# Samples will be ordered based on the CH DF samples enumerated, need to add that to all the dfs
bladder_ctDNA = bladder_ctDNA.drop("Samples_enumerated", axis = 1).merge(DF_Samples_enumerated)
chip_fr = chip_fr.drop("Samples_enumerated", axis = 1).merge(DF_Samples_enumerated)
CTDNA_pt_counts = CTDNA_pt_counts.drop("Samples_enumerated", axis = 1).merge(DF_Samples_enumerated)
pts_ctdna_NO_nan = pts_ctdna_NO_nan.drop("Samples_enumerated", axis = 1).merge(DF_Samples_enumerated)
pts_ctdna_nan = pts_ctdna_nan.drop("Samples_enumerated", axis = 1).merge(DF_Samples_enumerated)
sex_df = sex_df.drop("Samples_enumerated", axis = 1).merge(DF_Samples_enumerated)
age_df = age_df.drop("Samples_enumerated", axis = 1).merge(DF_Samples_enumerated)

ax_main = plot_oncoprint_combined(ax_main, bladder_ctDNA, bladder_ch, bar_height, bar_width, genes = list(bladder_OP_genes_dict.keys()))
[ax_mut_counts, ax_secondary] = plot_mut_counts_per_gene(mut_counts_df, vaf_df, ax_mut_counts, diagnosis = "Bladder")
ax_mut_counts_patient = plot_mut_counts_per_patient(ax_mut_counts_patient, pt_counts, bar_width)
[ax_sex, ax_age] = plot_sex_and_age(ax_sex, ax_age, sex_df, age_df, bar_width, bar_height = 9.5)
ax_tc = plot_ctdna_fractions(ax_tc, pts_ctdna_NO_nan, pts_ctdna_nan, bar_width)

# Minor tweaks
ax_mut_counts.set_xlim((0, 50))
ax_mut_counts.set_xticks([0, 25, 50])
ax_mut_counts.set_xticklabels(["0", "25", "50"])

ax_mut_counts_patient.set_ylim((0, 10))
ax_mut_counts_patient.set_yticks([0, 5, 10])
ax_mut_counts_patient.set_yticklabels(["0", "5", "10"])

# ax_mut_counts = add_legend(mut_dict_kidney, bbox_location = (2, 0.05), loc = 'lower right', ax = ax_mut_counts, age_df = age_df_bladder)
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/OP_bladder_ctDNA.pdf")
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/OP_bladder_ctDNA.png")


##############################################################
# RUNNING ONCOPRINT FOR KIDNEY SAMPLES
PATH_gene_categories_kidney_somatic = "/groups/wyattgrp/users/amunzur/pipeline/resources/panel/for_kidney_somatic_op.tsv"

fig_width = 8.5
fig_height = 5
tc_height = 2
mutcounts_height = 2
sex_height = 0.5
age_height = 0.5
main_height = all_vars_chip["Gene"].unique().shape[0]

fig = plt.figure(figsize = (fig_width, fig_height))
gs = mpl.gridspec.GridSpec(ncols = 2, nrows = 7, height_ratios = [tc_height, mutcounts_height, sex_height, age_height, main_height, 2, 3], width_ratios = [1, 0.1], hspace = 0, wspace = 0)

# Add and space out subplots
ax_tc = fig.add_subplot(gs[0, 0]) 
ax_mut_counts_patient = fig.add_subplot(gs[1, 0], sharex = ax_tc) 
ax_main = fig.add_subplot(gs[4, 0], sharex = ax_tc)
ax_mut_counts = fig.add_subplot(gs[4, 1], sharey = ax_main) 
ax_sex = fig.add_subplot(gs[2, 0], sharex = ax_main) 
ax_age = fig.add_subplot(gs[3, 0], sharex = ax_main) 
ax_legend = fig.add_subplot(gs[6, 0], sharex = ax_main) 

bar_height = 0.9
bar_width = 0.8

[kidney_ctDNA, gene_order_df, chip_fr, mut_counts_df, vaf_df, pt_counts, pts_ctdna_NO_nan, pts_ctdna_nan, DF_Samples_enumerated, mut_dict_kidney, sex_df, age_df_kidney, combined_age_df] = prepare_datasets_for_plotting(kidney_ctDNA_main, mut_dict, PATH_gene_categories_kidney_somatic, PATH_mutation_ctfractions, PATH_bladder_clinical_data, PATH_kidney_clinical_data)
ax_main = plot_oncoprint(ax_main, kidney_ctDNA, bar_height, bar_width)
ax_main = beautify_OP_ax(ax_main, gene_order_df, DF_Samples_enumerated, xlabel = "Baseline mRCC samples")
[ax_mut_counts, ax_secondary] = plot_mut_counts_per_gene(mut_counts_df, vaf_df, ax_mut_counts, diagnosis = "Kidney")
ax_mut_counts_patient = plot_mut_counts_per_patient(ax_mut_counts_patient, pt_counts, bar_width)
[ax_sex, ax_age] = plot_sex_and_age(ax_sex, ax_age, sex_df, age_df_kidney, bar_width, bar_height = bar_height)
ax_tc = plot_ctdna_fractions(ax_tc, pts_ctdna_NO_nan, pts_ctdna_nan, bar_width)
ax_legend = add_legend(mut_dict_kidney, ax = ax_legend, age_df = combined_age_df, fig = fig)

# Minor tweaks
ax_mut_counts.set_xlim((0, 40))
ax_mut_counts.set_xticks([0, 20, 40])
ax_mut_counts.set_xticklabels(["0", "20", "40"])

ax_secondary.set_xlim((0, 50))
ax_secondary.set_xticks([0, 25, 50])
ax_secondary.set_xticklabels(["0", "25", "50"])

ax_mut_counts_patient.set_ylim((0, 6))
ax_mut_counts_patient.set_yticks([0, 3, 6])
ax_mut_counts_patient.set_yticklabels(["0", "3", "6"])

# gs.tight_layout(fig)
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/OP_kidney_cfDNA.pdf")
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/OP_kidney_cfDNA.png")