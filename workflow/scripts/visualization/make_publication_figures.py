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
PATH_kidney_clinical = os.path.join(DIR_working, "resources/clinical_data/RCC clinical data - mRCC clinical Data.csv")
PATH_clinical_bladder = os.path.join(DIR_working, "resources/clinical_data/Bladder_enrollment.csv")
PATH_treatment_landscape = os.path.join(DIR_working, "resources/clinical_data/bladder/treatment.csv")
PATH_mutation_ctfractions = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/ctDNA_fractions.csv"
figure_dir = os.path.join(DIR_working, "results/figures/amazing_figures")
source_functions = os.path.join(DIR_working, "workflow/scripts/visualization/UTILITIES_make_chip_plots.py")
path_gene_exons = os.path.join(DIR_working, "resources/references/panel_genes_exons_refseq.tsv")

color_dict = {"Bladder": "deepskyblue", "Kidney": "orangered"}

with open(source_functions, 'r') as file:
    script_code = file.read()

exec(script_code)

# LOAD CHIP DATASETS
all_vars_chip = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/chip.csv"))
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

PATH_kidney_clinical = os.path.join(DIR_working, "resources/clinical_data/Kidney/mRCC clinical Data_clean.csv")
PATH_clinical_bladder = os.path.join(DIR_working, "resources/clinical_data/Bladder_enrollment.csv")

kidney_age = pd.read_csv(PATH_kidney_clinical)[["Patient_id", "Age at GUBB draw"]].rename(columns = {"Age at GUBB draw": "Age"}).assign(Diagnosis="Kidney")
bladder_age = pd.read_csv(PATH_clinical_bladder)[["Patient_id", "Age at baseline blood draw"]].rename(columns = {"Age at baseline blood draw": "Age"}).assign(Diagnosis="Bladder")
age_df = pd.concat([kidney_age, bladder_age]).reset_index(drop=True)

baseline = pd.concat([base_kidney_chip, base_bladder_chip]).reset_index(drop = True)
baseline_somatic = pd.concat([base_kidney_somatic, base_bladder_somatic]).reset_index(drop = True)

kidney_clinical = pd.read_csv(PATH_kidney_clinical)
bladder_clinical = pd.read_csv(PATH_clinical_bladder)

# START PLOTTING
# PLOT 1. Number of CH mutations and number of patients plots. One for kidney and one for bladder. They are stacked on top of each other.
fig = plt.figure(figsize=(2.5, 2.7))
gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1], hspace = 0.4)
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])
ax1 = plot_per_patient_counts(base_kidney_chip, figure_dir, f"patient_mutation_counts_baseline_kidney.pdf", colorby = None, ax = ax1, bar_color = color_dict["Kidney"]) # number of patients and number of mutations
ax1.set_xlabel("")
ax1.set_xticklabels(ax1.get_xticklabels(), fontsize = 9)
ax1.set_yticks([0, 10, 20, 30, 40, 50])
ax1.set_yticklabels(["0", "10", "20", "30", "40", "50"])
ax2 = plot_per_patient_counts(base_bladder_chip, figure_dir, f"patient_mutation_counts_baseline_bladder.pdf", colorby = None, ax = ax2, bar_color = color_dict["Bladder"]) # number of patients and number of mutations
gs.tight_layout(fig)
fig.savefig(os.path.join(figure_dir, "patient_mutation_counts.pdf"))

# PLOT 2. The proportion of the study cohort that is CH+ based on age.
plot_age_plots_stacked(pd.concat([base_kidney_chip, base_bladder_chip]), age_df, PATH_sample_information, figure_dir, 10, color_dict, ext = "pdf")

# PLOT 3. Grouped boxplots for CH+ and CH-. Separated by kidney - bladder as well.
age_vs_CH_presence_grouped(pd.concat([base_kidney_chip, base_bladder_chip]), age_df, "Age", PATH_sample_information, "Age and CH status", figure_dir, gene = None, min_WBC_VAF_threshold = None, test_to_use = "MWU", fontsize = 10)

# PLOT 4. 2 Upset plots. 
generate_upset_plot(baseline, os.path.join(figure_dir, "TP53_PPM1D_ATM_upset.pdf"), color_dict, genes = ["TP53", "PPM1D", "ATM"])
generate_upset_plot(baseline, os.path.join(figure_dir, "DTA_upset.pdf"), color_dict, genes = ["DNMT3A", "TET2", "ASXL1"])

# PLOT 5. Pie charts for CH presence absence in the two cohorts.
plot_pie_for_ch_presence_in_df(base_kidney_chip, PATH_sample_information, figure_dir, f"Kidney", f"kidney_pie_cohort_baseline.pdf", color_dict = ["orangered", "#f9b7a7"])
plot_pie_for_ch_presence_in_df(base_bladder_chip, PATH_sample_information, figure_dir, f"Bladder", f"bladder_pie_cohort_baseline.pdf", color_dict = ["deepskyblue", "#b2e4f7"])

# PLOT 6. Scatter plot for vafs
plot_vaf_scatter(all_vars_chip, figure_dir, f"baseline_vaf_correlation.pdf", color_dict, colorby = "Diagnosis", fontsize = 10)

# PLOT 7. Gene counts colored by mut type
# plot_gene_counts(baseline, figure_dir, "gene_counts.pdf", "Bladder and kidney", unique = True)
# plot_gene_counts(base_kidney_chip, figure_dir, "kidney_gene_counts.pdf", "Kidney", unique = True)
# plot_gene_counts(base_bladder_chip, figure_dizr, "bladder_counts_mirror.pdf", kidney_clinical, bladder_clinical)

gene_list = ["DNMT3A", "TET2", "PPM1D", "ASXL1", "TP53", "ATM", "CHEK2", "GNAS", "KMT2D", "STAG2", "CBL", "RAD21", "BRCC3", "JAK2", "SF3B1"]
gene_colors = 
plot_gene_counts_grouped(all_vars_chip, figure_dir, "gene_counts_grouped.pdf", kidney_clinical, bladder_clinical, gene_list = gene_list)







# PLOT 8. 
plot_gene_counts_percent_miscalled(baseline_somatic, baseline, True, figure_dir, "ctDNA_CH_perc_1_perc.png")
plot_gene_counts_percent_miscalled(baseline_somatic, baseline, False, figure_dir, "ctDNA_CH_perc.png")

# PLOT 8. DETECTION LIMIT
fig = plt.figure(figsize=(2.9, 2.7))
gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1], hspace = 0, wspace = 0)
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])
ax1 = plot_detection_limit(baseline, PATH_sample_information, os.path.join(figure_dir, "lod_plot.pdf"), ax = ax1, legendon = False)
ax2 = plot_detection_limit_per_gene(baseline, genes_list = ["TP53", "ATM", "CHEK2", "BRCA1", "BRCA2", "CHEK2"], ax = ax2, fontsize = 10)
gs.tight_layout(fig)
fig.savefig(os.path.join(figure_dir, "lod_plot.pdf"))

# PLOT 8. Are pts with immune AEs are significantly older than pts with no AEs? Just for kidney.
PATH_clin_data = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/Kidney/mRCC clinical Data_clean.csv"
clin = pd.read_csv(PATH_clin_data)[["Patient_id", "irAE", "Age at GUBB draw"]]
clin = clin[~pd.isnull(clin["irAE"])]
clin["Grade 3 or 4 irAE"] = clin["irAE"] > 0
# clin.loc[clin["irAE bool"] == True, "irAE bool"] = "Present"
# clin.loc[clin["irAE bool"] == False, "irAE bool"] = "Absent"
plot_boxplot(clin, "Age at GUBB draw", "irAE bool", ["Black", "Black"], ["Black", "Black"], figure_dir, "irAE_and_ages")

clin = clin.merge(data, how = "left")
data['CH status'] = data['CH status'].astype('category')

# Perform one-hot encoding for 'CH status'


clin = clin[~pd.isnull(clin["irAE"])]
X =clin[['CH status', 'Age at GUBB draw']] # Dropping first column to avoid multicollinearity
# clin["Grade 3 or 4 irAE"] = clin["irAE"] > 0

# X['Age at GUBB draw'] = clin['Age at GUBB draw']  # Add 'Age at GUBB draw' column
y = clin['Grade 3 or 4 irAE']  # Dependent variable

# Add constant for intercept term
X = sm.add_constant(X)

# Fit logistic regression model
logit_model = sm.Logit(y, X)

# Obtain results
logit_result = logit_model.fit()
print(logit_result.summary())




clin['Age_Positive_Interact'] = clin['Age at GUBB draw'] * (clin['CH status'] == 'Positive')
clin["CH status dummy"] = pd.get_dummies(clin['CH status'], drop_first=True)  
clin = clin.dropna()

X = clin[['Age at GUBB draw', 'CH status dummy', 'Age_Positive_Interact']]  # Predictor variables
X = sm.add_constant(X)  # Add constant for intercept
y = clin['Grade 3 or 4 irAE']  # Outcome variable

# Fit logistic regression model with interaction term
logit_model = sm.Logit(y, X)
logit_result = logit_model.fit()

# Print model summary
print(logit_result.summary())


x = data[["Patient_id", "CH status", "irAE bool"]]
y = x.merge(clin[["Patient_id", "Age at GUBB draw"]], how = "left")

data['CH status'] = data['CH status'].astype('category')

data = y 