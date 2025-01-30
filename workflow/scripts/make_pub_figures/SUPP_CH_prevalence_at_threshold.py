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
all_vars_chip = all_vars_chip[all_vars_chip["VAF_n"]>=1]

fig = plt.figure(figsize=(8, 3.14))
inner_gs0 = gridspec.GridSpec(3, 1, height_ratios = [0.008, 1, 1], hspace = 0.3)
inner_inner_gs0 = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios = [1, 0.3], wspace=0.2, hspace = 0.2, subplot_spec=inner_gs0[2])

ax0 = fig.add_subplot(inner_gs0[1])
ax1 = fig.add_subplot(inner_inner_gs0[0], sharex=ax0)
ax2 = fig.add_subplot(inner_inner_gs0[1], sharex=ax0)

ax0, ax1, ax2 = plot_fraction_of_CH_calls_with_bars(df_somatic_main = all_vars_ctDNA, df_ch_main = all_vars_chip, ax0 = ax0, ax1 = ax1, ax2 = ax2, stacked = False, yval_max = -110) # plot_fraction_of_CH_calls was the old function, still available in the source functions.

inner_gs0.tight_layout(fig)
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/SUPP_CH_fraction.png")
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/SUPP_CH_fraction.pdf")