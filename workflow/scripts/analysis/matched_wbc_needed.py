"""
What perc of variants we would m
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
from adjustText import adjust_text
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.cm import get_cmap
from datetime import datetime
import matplotlib.cm as cm
from lifelines.statistics import logrank_test
from scipy.stats import fisher_exact
from lifelines import KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
import matplotlib.ticker as ticker
import upsetplot
from scipy.stats import mannwhitneyu
from lifelines import CoxPHFitter
import matplotlib.dates as mdates

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
PATH_kidney_clinical = os.path.join(DIR_working, "resources/clinical_data/Kidney/mRCC clinical Data_clean.csv")
PATH_clinical_bladder = os.path.join(DIR_working, "resources/clinical_data/Bladder_enrollment.csv")
PATH_treatment_landscape = os.path.join(DIR_working, "resources/clinical_data/bladder/treatment.csv")
PATH_mutation_ctfractions = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/ctDNA_fractions.csv"
figure_dir = os.path.join(DIR_working, "results/figures/amazing_figures")
source_functions = os.path.join(DIR_working, "workflow/scripts/visualization/UTILITIES_make_chip_plots.py")

color_dict = {"Bladder": "deepskyblue", "Kidney": "orangered"}

with open(source_functions, 'r') as file:
    script_code = file.read()

exec(script_code)

# LOAD CHIP DATASETS
all_vars_chip = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/chip.csv"))
sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
all_vars_chip = sample_info.merge(all_vars_chip, how = "inner")
all_vars_chip = all_vars_chip[all_vars_chip["Dependent"] == False]

base_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
prog_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)
prog_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)

#################################################################
# Analysis 1. What percentage of variants would we miss if we didn't have the matched WBC in the progression samples?
#################################################################

genes = ["TP53", "ATM", "CHEK2", "PPM1D", "ERBB2", "TERT", "ASXL1", "KMT2D"]
base_bladder_chip = base_bladder_chip[base_bladder_chip["Gene"].isin(["TP53", "ATM", "CHEK2", "PPM1D", "ERBB2", "TERT", "ASXL1", "KMT2D"])].reset_index(drop = True)
prog_bladder_chip = prog_bladder_chip[prog_bladder_chip["Gene"].isin(["TP53", "ATM", "CHEK2", "PPM1D", "ERBB2", "TERT", "ASXL1", "KMT2D"])].reset_index(drop = True)

base_bladder = base_bladder_chip[["Patient_id", "Timepoint", 'Chrom', 'Position', 'Ref', 'Alt', "Gene", "Protein_annotation", "Consequence", "VAF_n", "VAF_t"]]
ot_bladder = prog_bladder_chip[["Patient_id", "Timepoint", 'Chrom', 'Position', 'Ref', 'Alt', "Gene", "Protein_annotation", "Consequence", "VAF_n", "VAF_t"]]
merged = ot_bladder.merge(base_bladder, on = ["Patient_id", 'Chrom', 'Position', 'Ref', 'Alt', "Gene", "Protein_annotation", "Consequence"], indicator = True, how = "left")

missed = merged[merged["_merge"] == "left_only"] # we would be mischaracterizing these as treatment related without having the progression WBCs - essentially these are the CH mutatoin that emerged later on during treatment with no evidence at baseline 
x = missed.drop_duplicates(["Patient_id", "Gene", "Protein_annotation"])["Gene"].value_counts().reset_index().rename(columns = {"index": "Gene", "Gene": "Counts"}).reset_index()
missed = missed.merge(x, how = "left") # helps with plotting

# plotting
fig, ax = plt.subplots(figsize=(3.5, 3))
ax2 = ax.twinx()
ax.bar(x["Gene"], x["Counts"], color='#4fb948', edgecolor = None, linewidth=0.0)
# plot the scatter
for i, row in missed.iterrows():
    ax2.scatter(np.random.uniform(-0.02, 0.2, 1) + row['index'], row["VAF_t_x"], s = 10, color = "black", zorder = 100)

ax.spines["top"].set_visible(False)
ax.set_xticklabels(x["Gene"].tolist(), fontsize = 9)
ax2.spines["top"].set_visible(False)

ax.set_ylabel("Mutation count")
ax2.set_ylabel("Plasma VAF%")
ax2.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
ax2.set_yticklabels(["0", "0.2", "0.4", "0.6", "0.8", "1.0"])

fig.tight_layout()
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/bladder/CH_wbc_prog.png")