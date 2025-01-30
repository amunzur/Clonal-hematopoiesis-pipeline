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
all_vars_ctDNA = all_vars_ctDNA[(all_vars_ctDNA["Dependent"] == False) & (all_vars_ctDNA["Timepoint"] == "Baseline")]
all_vars_ctDNA["Consequence"] = all_vars_ctDNA["Consequence"].replace({"Frameshift deletion":"Frameshift indel", "Frameshift insertion":"Frameshift indel", "Nonframeshift deletion":"Nonframeshift indel", "Nonframeshift insertion":"Nonframeshift indel"})
bladder_ctDNA_main = all_vars_ctDNA[all_vars_ctDNA["Diagnosis"] == "Bladder"]
kidney_ctDNA_main = all_vars_ctDNA[all_vars_ctDNA["Diagnosis"] == "Kidney"]

# CH mutations
all_vars_chip = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/chip.csv"))
all_vars_chip = all_vars_chip[(all_vars_chip["Dependent"] == False) & (all_vars_chip["Timepoint"] == "Baseline")]
all_vars_chip["Consequence"] = all_vars_chip["Consequence"].replace({"Frameshift deletion":"Frameshift indel", "Frameshift insertion":"Frameshift indel", "Nonframeshift deletion":"Nonframeshift indel", "Nonframeshift insertion":"Nonframeshift indel"})
bladder_chip_main = all_vars_chip[all_vars_chip["Diagnosis"] == "Bladder"]
kidney_chip_main = all_vars_chip[all_vars_chip["Diagnosis"] == "Kidney"]

# source functions
source_functions = "/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/make_pub_figures/make_oncoprint_functions.py"
with open(source_functions, 'r') as file:
    script_code = file.read()

exec(script_code)

# OPs for the ctDNA mutations
##############################################################
# BLADDER
bladder_genes = [
    'DNMT3A', 'TET2', 'ASXL1', 'JAK2', 'ATM', 'BRCA2', 'PPM1D', 'CHEK2', 'TP53', 
    'ARID1A', 'ERBB2', 'STAG2', 'RB1', 'FGFR3', 'PIK3CA', 'KMT2D', 'KDM6A', 'TERT'
    ]

bladder_ctDNA_main = bladder_ctDNA_main[bladder_ctDNA_main["Gene"].isin(bladder_genes)].reset_index(drop = True)
bladder_chip_main = bladder_chip_main[bladder_chip_main["Gene"].isin(bladder_genes)].reset_index(drop = True)

# Compute gene and sample ordering
gene_order_df         = assign_mutation_values(genes = genes, muts_ch = bladder_chip_main, muts_ctDNA = bladder_ctDNA_main, bar_height=1)
ctdna_fr              = return_ctDNA_fractions(PATH_sample_information, PATH_mutation_ctfractions, bladder_ctDNA_main, bladder_chip_main)

# Sample ordering
ctdna_pt_counts = calculate_patient_counts(ctdna_muts)
ch_pt_counts = calculate_patient_counts(ch_muts)
combined_counts = ctdna_pt_counts.merge(ch_pt_counts, how = "outer", on = "Patient_id").fillna(0)
combined_counts["Sum count"] = combined_counts["Count_x"] + combined_counts["Count_y"]
combined_counts = combined_counts.sort_values(by = "Sum count", ascending = False)

samples_enumerated = combined_counts.reset_index(drop = True).reset_index()[["index", "Patient_id"]].rename(columns = {"index": "Samples_enumerated"})
ctdna_fr = ctdna_fr.merge(samples_enumerated, how = "left")

# Generate ctDNA DFs
ctdna_muts          = calculate_mutation_burden_per_gene(bladder_ctDNA_main)
ctdna_muts          = assign_mutations_to_colors(ctdna_muts, mut_dict).merge(samples_enumerated, how = "left").merge(gene_order_df[["Gene", "ctDNA position"]], how = "left").rename(columns = {"ctDNA position": "Order on oncoprint"})
# ctdna_mut_counts    = get_mutation_counts_per_gene(ctdna_muts)[0].merge(gene_order_df[["Gene", "ctDNA position"]], how = "left").rename(columns = {"ctDNA position": "Order on oncoprint"}).set_index("Order on oncoprint").drop("Gene", axis = 1)
ctdna_vaf           = get_mutation_counts_per_gene(ctdna_muts)[1].merge(gene_order_df[["Gene", "ctDNA position"]], how = "left").rename(columns = {"ctDNA position": "Order on oncoprint"}).set_index("Order on oncoprint").drop("Gene", axis = 1)
ctdna_pt_counts     = calculate_patient_counts(ctdna_muts).merge(samples_enumerated, how = "left")

# Generate CH DFs
ch_muts         = calculate_mutation_burden_per_gene(bladder_chip_main)
ch_muts         = assign_mutations_to_colors(ch_muts, mut_dict).merge(samples_enumerated, how = "left").merge(samples_enumerated, how = "left").merge(gene_order_df[["Gene", "CH position"]], how = "left").rename(columns = {"CH position": "Order on oncoprint"})
# ch_mut_counts   = get_mutation_counts_per_gene(ch_muts)[0].merge(gene_order_df[["Gene", "CH position"]], how = "left").rename(columns = {"CH position": "Order on oncoprint"}).set_index("Order on oncoprint").drop("Gene", axis = 1)
ch_vaf          = get_mutation_counts_per_gene(ch_muts)[1].merge(gene_order_df[["Gene", "CH position"]], how = "left").rename(columns = {"CH position": "Order on oncoprint"}).set_index("Order on oncoprint").drop("Gene", axis = 1)
ch_pt_counts    = calculate_patient_counts(ch_muts).merge(samples_enumerated, how = "left")
ch_fr           = calculate_ch_fraction(ch_muts).merge(samples_enumerated, how = "left")

age_df, sex_df = get_sex_and_age(PATH_bladder_clinical_data)
age_df = age_df.merge(samples_enumerated, how = "left")
sex_df = sex_df.merge(samples_enumerated, how = "left")

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
ax_legend = fig.add_subplot(gs[6, 0], sharex = ax_main) 

bar_height = 0.8
bar_width = 0.8

ax_main = plot_oncoprint_combined(ax_main, gene_order_df, ctdna_muts, ch_muts, bar_height, bar_width)
ax_main = beautify_OP_ax(ax_main, gene_order_df, samples_enumerated, xlabel = "Baseline mUC samples")
[ax_mut_counts, ax_secondary] = plot_mut_counts_per_gene(ctdna_mut_counts, ctdna_vaf, ax_mut_counts, mut_dict, bar_height = bar_height, diagnosis = "Bladder") # Plot ctDNA mutations
[ax_mut_counts, ax_secondary] = plot_mut_counts_per_gene(ch_mut_counts, ch_vaf, ax_mut_counts, mut_dict, bar_height = bar_height, diagnosis = "Bladder") # Plot CH mutations
ax_mut_counts_patient = plot_mut_counts_per_patient_CH_and_ctDNA_combined(ax_mut_counts_patient, ctdna_pt_counts, ch_pt_counts, bar_width)
[ax_sex, ax_age] = plot_sex_and_age(ax_sex, ax_age, sex_df, age_df, bar_width, bar_height = 9.5)
ax_tc = plot_ctdna_fractions(ax_tc, ctdna_fr, bar_width)
ax_legend = add_legend(mut_dict, ax = ax_legend, age_df = age_df, fig = fig, ct_and_ctdna_combined = True)

fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/OP_bladder.pdf")
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/OP_bladder.png")


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