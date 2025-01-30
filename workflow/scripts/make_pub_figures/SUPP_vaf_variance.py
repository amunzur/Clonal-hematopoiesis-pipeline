"""
Similar to make_chip_plots.py. Generates paper figures.
"""

import pandas as pd
import numpy as np
import os 
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

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
figure_dir = os.path.join(DIR_working, "results/figures/pub_figures")
source_functions = os.path.join(DIR_working, "workflow/scripts/visualization/UTILITIES_make_chip_plots.py")

with open(source_functions, 'r') as file:
    script_code = file.read()

exec(script_code)

# LOAD CHIP DATASETS
all_vars_chip = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/chip.csv"))
all_vars_chip = all_vars_chip[all_vars_chip["Dependent"] == False]
base_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
prog_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)
prog_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)

baseline = pd.concat([base_kidney_chip, base_bladder_chip]).reset_index(drop = True)

pt_counts = baseline["Patient_id"].value_counts().reset_index()
pt_counts.columns = ["Patient_id", "nmuts"]
pt_counts = pt_counts[pt_counts["nmuts"] > 1].reset_index(drop = True)

baseline = baseline[baseline["Patient_id"].isin(pt_counts["Patient_id"])]
baseline_somatic = pd.concat([base_kidney_somatic, base_bladder_somatic]).reset_index(drop = True)

fig = plt.figure(figsize=(7, 2.5))
gs = gridspec.GridSpec(1, 2, width_ratios=[1.2, 1], hspace = 0, wspace = 0) # outer most gs with 3 rows

kidney_ax = plt.subplot(gs[0])
bladder_ax = plt.subplot(gs[1])
kidney_ax, kidney_twin_ax = plot_variance_per_pt(baseline[baseline["Diagnosis"] == "Kidney"], kidney_ax, scatter_color = "orangered", bar_color = "#FFEDE6")
bladder_ax, bladder_twin_ax = plot_variance_per_pt(baseline[baseline["Diagnosis"] == "Bladder"], bladder_ax, scatter_color = "deepskyblue", bar_color = "#E8F9FF")

# minor aes
kidney_twin_ax.set_ylabel("")
kidney_twin_ax.set_yticks([])
kidney_twin_ax.set_yticklabels([])
kidney_twin_ax.spines["right"].set_visible(False)
kidney_ax.spines["right"].set_visible(False)
kidney_ax.set_xlabel("mRCC")

bladder_ax.set_ylabel("")
bladder_ax.set_yticks([])
bladder_ax.set_yticklabels([])
bladder_ax.spines["left"].set_visible(False)
bladder_twin_ax.spines["left"].set_visible(False)
bladder_ax.set_xlabel("mUC")

gs.tight_layout(fig)
fig.savefig(os.path.join(figure_dir, "SUPP_vaf_variance.png"))
fig.savefig(os.path.join(figure_dir, "SUPP_vaf_variance.pdf"))