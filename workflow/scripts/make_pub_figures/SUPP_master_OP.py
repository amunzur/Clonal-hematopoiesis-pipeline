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
    "multiple": "#00FFFF",
    "missense": '#79B443',
    "stopgain": '#BD4398',
    "frameshift indel": '#FFC907',
    "nonframeshift indel": '#a9a9a9',
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
all_vars_ctDNA["Consequence"] = all_vars_ctDNA["Consequence"].str.lower()
all_vars_ctDNA["Consequence"] = all_vars_ctDNA["Consequence"].replace(
    {"frameshift deletion":"frameshift indel", 
     "frameshift insertion":"frameshift indel", 
     "nonframeshift deletion":"nonframeshift indel", 
     "nonframeshift insertion":"nonframeshift indel", 
     "startloss": "missense"})

# CH mutations
all_vars_chip = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/chip.csv"))
all_vars_chip = all_vars_chip[(all_vars_chip["Dependent"] == False) & (all_vars_chip["Timepoint"] == "Baseline")]
all_vars_chip["Consequence"] = all_vars_chip["Consequence"].replace({"Frameshift deletion":"Frameshift indel", "Frameshift insertion":"Frameshift indel", "Nonframeshift deletion":"Nonframeshift indel", "Nonframeshift insertion":"Nonframeshift indel"})
all_vars_chip["Consequence"] = all_vars_chip["Consequence"].str.lower()
all_vars_chip["Consequence"] = all_vars_chip["Consequence"].replace(
    {"frameshift deletion":"frameshift indel", 
     "frameshift insertion":"frameshift indel", 
     "nonframeshift deletion":"nonframeshift indel", 
     "nonframeshift insertion":"nonframeshift indel", 
     "startloss": "missense"})

# source functions
source_functions = "/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/make_pub_figures/make_oncoprint_functions.py"
with open(source_functions, 'r') as file:
    script_code = file.read()

exec(script_code)

# OPs for the ctDNA mutations
##############################################################
# Get clinical data
bladder_age_df, bladder_sex_df = get_sex_and_age(PATH_bladder_clinical_data)
bladder_diagnosis_df = bladder_sex_df[["Patient_id"]].assign(Diagnosis="mUC", diagnosis_color="deepskyblue")

kidney_age_df, kidney_sex_df = get_sex_and_age(PATH_kidney_clinical_data)
kidney_diagnosis_df = kidney_sex_df[["Patient_id"]].assign(Diagnosis="mRCC", diagnosis_color="orangered")

age_df = pd.concat([bladder_age_df, kidney_age_df])
sex_df = pd.concat([bladder_sex_df, kidney_sex_df])
diagnosis_df = pd.concat([bladder_diagnosis_df, kidney_diagnosis_df])

# Compute gene ordering
genes_to_include = [
    'DNMT3A', 'TET2', 'ASXL1', 'JAK2', 'ATM', 'BRCA2', 'PPM1D', 'CHEK2', 'TP53', 
    'ARID1A', 'ERBB2', 'STAG2', 'RB1', 'FGFR3', 'PIK3CA', 'KMT2D', 'KDM6A', 'TERT',
    "BAP1", "MTOR", "PBRM1", "SETD2", "VHL"
    ]

ctDNA_main = all_vars_ctDNA[all_vars_ctDNA["Gene"].isin(genes_to_include)].reset_index(drop = True)
chip_main = all_vars_chip[all_vars_chip["Gene"].isin(genes_to_include)].reset_index(drop = True)
gene_order_df = assign_mutation_values(genes = genes_to_include, muts_ch = chip_main, muts_ctDNA = ctDNA_main, bar_height=1, omit_missing_row = True)

# ctDNA dfs
ctdna_plotting_df = generate_plotting_df(ctDNA_main, mut_dict, gene_order_df) # makes up bulk of the plot
ctdna_vaf = get_mutation_counts_per_gene(ctDNA_main)[1].merge(gene_order_df[["Gene", "ctDNA position"]], how = "left").rename(columns = {"ctDNA position": "Order on oncoprint"}).set_index("Order on oncoprint").drop("Gene", axis = 1)
ctdna_mut_counts = get_mutation_counts_per_gene(ctDNA_main)[0].merge(gene_order_df[["Gene", "ctDNA position"]], how = "left").rename(columns = {"ctDNA position": "Order on oncoprint"}).set_index("Order on oncoprint").drop("Gene", axis = 1)
ch_pt_counts = calculate_patient_counts(chip_main)

# CH dfs
ch_plotting_df = generate_plotting_df(chip_main, mut_dict, gene_order_df) # makes up bulk of the plot
ch_vaf = get_mutation_counts_per_gene(chip_main)[1].merge(gene_order_df[["Gene", "CH position"]], how = "left").rename(columns = {"CH position": "Order on oncoprint"}).set_index("Order on oncoprint").drop("Gene", axis = 1)
ch_mut_counts   = get_mutation_counts_per_gene(chip_main)[0].merge(gene_order_df[["Gene", "CH position"]], how = "left").rename(columns = {"CH position": "Order on oncoprint"}).set_index("Order on oncoprint").drop("Gene", axis = 1)
ctdna_pt_counts = calculate_patient_counts(ctDNA_main)

# Get samples enumerated
# samples_enumerated = ch_plotting_df.merge(diagnosis_df).merge(ch_pt_counts)
# samples_enumerated["Diagnosis"] = pd.Categorical(samples_enumerated["Diagnosis"], categories=["mUC", "mRCC"], ordered=True)
# samples_enumerated["Gene"] = pd.Categorical(samples_enumerated["Gene"], categories=genes_to_include, ordered=True)
# samples_enumerated["Consequence"] = pd.Categorical(samples_enumerated["Consequence"], categories=mut_dict.keys(), ordered=True)
# samples_enumerated = samples_enumerated.sort_values(by=["Diagnosis", "Gene", "Consequence", "Count"], ascending=[True, True, True, False]).reset_index(drop = True).reset_index()
# samples_enumerated = samples_enumerated[["Patient_id", "index"]].rename(columns = {"index": "Samples_enumerated"})


ch_plotting_df_renamed = ch_plotting_df.rename(columns=lambda x: 'ch_' + x if x != 'Patient_id' else x)
ctdna_plotting_df_renamed = ctdna_plotting_df.rename(columns=lambda x: 'ctdna_' + x if x != 'Patient_id' else x)
samples_enumerated = ch_plotting_df_renamed.merge(ctdna_plotting_df_renamed, how = "outer").merge(diagnosis_df).merge(ch_pt_counts)
samples_enumerated = samples_enumerated[[
    "Patient_id", 
    'ch_Gene',
    'ch_Consequence',
    'ch_Counts',
    'ctdna_Gene', 
    'ctdna_Consequence',
    'ctdna_Counts',
    'Diagnosis'
]]

samples_enumerated["Diagnosis"] = pd.Categorical(samples_enumerated["Diagnosis"], categories=["mUC", "mRCC"], ordered=True)
samples_enumerated["ch_Gene"] = pd.Categorical(samples_enumerated["ch_Gene"], categories=genes_to_include, ordered=True)
samples_enumerated["ch_Consequence"] = pd.Categorical(samples_enumerated["ch_Consequence"], categories=mut_dict.keys(), ordered=True)
samples_enumerated = samples_enumerated.sort_values(
    by=["Diagnosis", "ch_Gene", "ch_Consequence", "ch_Counts", 'ctdna_Gene', 'ctdna_Consequence', 'ctdna_Counts'], 
    ascending=[True, True, True, False, True, True, False]).reset_index(drop = True).reset_index()
samples_enumerated = samples_enumerated[["Patient_id", "index"]].rename(columns = {"index": "Samples_enumerated"})

age_df = age_df.merge(samples_enumerated, how = "left")
sex_df = sex_df.merge(samples_enumerated, how = "left")
diagnosis_df = diagnosis_df.merge(samples_enumerated, how = "left")

ctdna_plotting_df = ctdna_plotting_df.merge(samples_enumerated, how = "left")
ch_plotting_df = ch_plotting_df.merge(samples_enumerated, how = "left")

ch_pt_counts = ch_pt_counts.merge(samples_enumerated, how = "left")
ctdna_pt_counts = ctdna_pt_counts.merge(samples_enumerated, how = "left")

fig_width = 8
fig_height = 5
tc_height = 2
mutcounts_height = 2
sex_height = 0.5
age_height = 0.5
diagnosis_height = 0.5
main_height = len(genes_to_include)*2

fig = plt.figure(figsize = (fig_width, fig_height))
gs = mpl.gridspec.GridSpec(ncols = 2, nrows = 9, height_ratios = [tc_height, tc_height, mutcounts_height, sex_height, age_height, diagnosis_height, main_height, 2, 3], width_ratios = [1, 0.1], hspace = 0, wspace = 0)

# Add and space out subplots
ax_ctDNA_vaf = fig.add_subplot(gs[0, 0]) #Plotting the highest VAF per patient
ax_CH_vaf = fig.add_subplot(gs[1, 0]) #Plotting the highest VAF per patient

ax_mut_counts_patient = fig.add_subplot(gs[2, 0]) 
ax_main = fig.add_subplot(gs[6, 0], sharex = ax_mut_counts_patient)
ax_mut_counts = fig.add_subplot(gs[6, 1], sharey = ax_main) 
ax_sex = fig.add_subplot(gs[3, 0], sharex = ax_main) 
ax_age = fig.add_subplot(gs[4, 0], sharex = ax_main) 
ax_diagnosis = fig.add_subplot(gs[5, 0], sharex = ax_main) 
ax_legend = fig.add_subplot(gs[8, 0], sharex = ax_main) 

bar_height = 0.5
bar_width = 1

ax_ctDNA_vaf = plot_highest_VAF_per_patient(ax_ctDNA_vaf, ctDNA_main, samples_enumerated, bar_width = bar_width)
ax_CH_vaf = plot_highest_VAF_per_patient(ax_CH_vaf, chip_main, samples_enumerated, bar_width = bar_width)
ax_main = plot_oncoprint_compressed(ax_main, gene_order_df, ctdna_plotting_df, ch_plotting_df, bar_height, bar_width)
ax_main = beautify_OP_ax(ax_main, gene_order_df, samples_enumerated, xlabel = "Baseline samples", fontsize = 6)
[ax_mut_counts, ax_secondary] = plot_mut_counts_per_gene(ctdna_mut_counts, ctdna_vaf, ax_mut_counts, mut_dict, bar_height = bar_height, diagnosis = "Both", logscale = True, vaf_col = "VAF_t") # Plot ctDNA mutations
[ax_mut_counts, ax_secondary] = plot_mut_counts_per_gene(ch_mut_counts, ch_vaf, ax_mut_counts, mut_dict, bar_height = bar_height, diagnosis = "Both", logscale = True, vaf_col = "VAF_n") # Plot CH mutations
ax_mut_counts_patient = plot_mut_counts_per_patient_CH_and_ctDNA_combined(ax_mut_counts_patient, ctdna_pt_counts, ch_pt_counts, bar_width)
[ax_sex, ax_age, ax_diagnosis] = plot_sex_age_diagnosis(ax_sex, ax_age, ax_diagnosis, sex_df, age_df, diagnosis_df, bar_width, bar_height = 9.5)
# ax_main = add_gene_group_annotations(ax_main, bar_height)

ax_legend = add_legend(mut_dict, ax = ax_legend, age_df = age_df, fig = fig, ct_and_ctdna_combined = True)

fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/OP_compressed.pdf")
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/OP_compressed.png")