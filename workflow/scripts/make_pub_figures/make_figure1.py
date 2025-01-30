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
import statsmodels.api as sm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

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
PATH_kidney_clinical = os.path.join(DIR_working, "resources/clinical_data/Supplementary tables - Clinical data - mRCC.csv")
PATH_clinical_bladder = os.path.join(DIR_working, "resources/clinical_data/Bladder_enrollment.csv")
PATH_treatment_landscape = os.path.join(DIR_working, "resources/clinical_data/bladder/treatment.csv")
PATH_mutation_ctfractions = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/ctDNA_fractions.csv"

figure_dir = os.path.join(DIR_working, "results/figures/pub_figures")
source_functions = os.path.join(DIR_working, "workflow/scripts/visualization/UTILITIES_make_chip_plots.py")
path_gene_exons = os.path.join(DIR_working, "resources/references/panel_genes_exons_refseq.tsv")
sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
color_dict = {"Bladder": "deepskyblue", "Kidney": "orangered"}

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

# LOAD SOMATIC DATASETS
all_vars_somatic = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/somatic.csv"))
all_vars_somatic = all_vars_somatic[~all_vars_somatic["Patient_id"].isin(["20-313", "21-184", "21-430"])] # exclude some samples due to oxidative damage
base_kidney_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Kidney")].reset_index(drop = True)
prog_kidney_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "During treatment") & (all_vars_somatic["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Bladder")].reset_index(drop = True)
prog_bladder_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "During treatment") & (all_vars_somatic["Diagnosis"] == "Bladder")].reset_index(drop = True)

PATH_clinical_bladder = os.path.join(DIR_working, "resources/clinical_data/bladder/clinical_data.csv")
kidney_age = pd.read_csv(PATH_kidney_clinical)[["Patient_id", "Age at GUBB draw"]].rename(columns = {"Age at GUBB draw": "Age"}).assign(Diagnosis="Kidney")
# kidney_age.loc[kidney_age["Age"] == -2, "Age"] = 68

bladder_age = pd.read_csv(PATH_clinical_bladder)
bladder_age = bladder_age[bladder_age["First sample?"] == True][["Patient_id", "Age at blood draw"]].rename(columns = {"Age at blood draw": "Age"}).assign(Diagnosis="Bladder")
age_df = pd.concat([kidney_age, bladder_age]).reset_index(drop=True).drop_duplicates()
baseline = pd.concat([base_kidney_chip, base_bladder_chip]).reset_index(drop = True)
baseline_somatic = pd.concat([base_kidney_somatic, base_bladder_somatic]).reset_index(drop = True)

kidney_clin = pd.read_csv(PATH_kidney_clinical)
bladder_clin = pd.read_csv(PATH_clinical_bladder)
kidney_pts = sample_info[(sample_info["Diagnosis"] == "Kidney") & (sample_info["Timepoint"] == "Baseline")]["Patient_id"].tolist()
bladder_pts = sample_info[(sample_info["Diagnosis"] == "Bladder") & (sample_info["Timepoint"] == "Baseline")]["Patient_id"].tolist()
kidney_clin = kidney_clin[kidney_clin["Patient_id"].isin(kidney_pts)]
bladder_clin = bladder_clin[bladder_clin["Patient_id"].isin(bladder_pts)]

with open(source_functions, 'r') as file:
    script_code = file.read()

exec(script_code)

fig = plt.figure(figsize=(8, 9))
outer_gs = gridspec.GridSpec(3, 1, height_ratios=[0.8, 1.7, 1], hspace = 0, wspace = 0) # outer most gs with 3 rows
inner_gs1 = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=[1, 1], subplot_spec=outer_gs[1], wspace=0.15, hspace = 0.15)
inner_gs3 = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=[1, 1], subplot_spec=outer_gs[2], wspace=0.35, hspace = 0.3)
# inner_inner_gs3 = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=[1, 0.25], subplot_spec=inner_gs3[0], wspace=0.2, hspace = 0.3)

# inner_inner_gs3 = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios=[1, 1], subplot_spec=inner_gs3[0], hspace = 0.2)

inner_inner_left_gs1 = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios=[1, 1], subplot_spec=inner_gs1[0], wspace=0.5, hspace = 0.5)
inner_inner_left_top_gs1 = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=[1, 1], subplot_spec=inner_inner_left_gs1[0], wspace=0.5, hspace = 0.3)

# FIGURE 1C - vaf correlation
ax2 = plt.subplot(inner_inner_left_top_gs1[1])
ax2 = plot_vaf_scatter(baseline, ax2, annotate_genes = False)
outer_gs.tight_layout(fig)
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/figure1.png")

# FIGURE 1B
ax0 = plt.subplot(inner_inner_left_top_gs1[0])
ax0 = plot_ch_presence_absence_bars(baseline, PATH_sample_information, ax0, "CHIP")

# FIGURE 1D - number of muts
ax3 = plt.subplot(inner_inner_left_gs1[1])
ax3, ax3_twin = plot_per_patient_counts_grouped_bar_chart(baseline, PATH_sample_information, ax3)

# Gene counts grouped, old function we were using was plot_gene_counts_grouped.
gs_genes = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=[1, 1], subplot_spec=inner_gs1[1], wspace=0.4)
ax_genes_main = plt.subplot(gs_genes[0])
ax_genes_main.set_yticks([])
ax_genes_main.set_yticklabels([])
ax_genes_main.spines[["top", "left", "right", "bottom"]].set_visible(False)
ax_genes_main.tick_params(axis='x', labeltop = False)
ax_genes_main.tick_params(axis='y', labelleft=False, left = False, right = False, labelright = False)
ax_genes = ax_genes_main.twinx()
ax_swarm = plt.subplot(gs_genes[1], sharey = ax_genes)

gene_list = ['JAK2', "SF3B1", 'BRCC3', 'RAD21', 'CBL', 'STAG2', 'KMT2D', 'GNAS', 'CHEK2', 'ATM', 'TP53', 'ASXL1', 'PPM1D', 'TET2', 'DNMT3A']
ax_genes, gene_order = plot_gene_counts_grouped_single_ax(baseline, ax_genes, gene_list = gene_list, horizontal = True)
ax_swarm = plot_vafs_swarm(baseline, ax_swarm, revert_vafs = True, gene_list = gene_list, gene_order = gene_order, horizontal = True, add_violins = True)
ax_swarm.set_ylim(ax_genes.get_ylim())
ax_genes.set_xlabel("% mutated in CH+ pts")

# FIGURE 1G - CH and age box plots
ax6 = plt.subplot(inner_gs3[0])
ax6 = age_vs_CH_presence_grouped_ch_thresholds(df_CH= baseline, age_df= age_df, PATH_sample_information= PATH_sample_information, ax=ax6, fontsize = 10)

# FIGURE 1E - line plot and accompanying boxplots
ax4 = plt.subplot(inner_gs3[1]) # lineplot
ax4_forest_inset = inset_axes(ax4, width="25%", height="25%", loc='upper left')
ax4 = CH_presence_line_thresholds(baseline, age_df, PATH_sample_information, ax4, axforest = ax4_forest_inset, fontsize = 10)

fig.text(0.02, 0.99, 'a', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
fig.text(0.02, 0.79, 'b', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
fig.text(0.27, 0.79, 'c', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
fig.text(0.53, 0.79, 'd', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
fig.text(0.02, 0.56, 'e', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
fig.text(0.02, 0.33, 'f', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
fig.text(0.51, 0.33, 'g', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)

outer_gs.tight_layout(fig)
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/figure1.png")
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/figure1.pdf")