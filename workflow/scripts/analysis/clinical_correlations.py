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

############################################### OR analysis
# OR and age, and presence absence of CH mutations in select genes
PATH_kidney_clinical = os.path.join(DIR_working, "resources/clinical_data/Kidney/mRCC clinical Data_clean.csv")
PATH_clinical_bladder = os.path.join(DIR_working, "resources/clinical_data/Bladder_enrollment.csv")
kidney_clinical = pd.read_csv(PATH_kidney_clinical)
bladder_clinical = pd.read_csv(PATH_clinical_bladder)

kidney_age = pd.read_csv(PATH_kidney_clinical)[["Patient_id", "Age at GUBB draw"]].rename(columns = {"Age at GUBB draw": "Age"}).assign(Diagnosis="Kidney")
bladder_age = pd.read_csv(PATH_clinical_bladder)[["Patient_id", "Age at baseline blood draw"]].rename(columns = {"Age at baseline blood draw": "Age"}).assign(Diagnosis="Bladder")
age_df = pd.concat([kidney_age, bladder_age]).reset_index(drop=True)

gene_list = ["DNMT3A", "TET2", "PPM1D", "ASXL1", "TP53", "ATM", "CHEK2", "GNAS", "KMT2D", "STAG2", "CBL", "RAD21", "BRCC3", "JAK2", "SF3B1"]

results_kidney = calculate_OR(base_kidney_chip, kidney_clinical, gene_list, "Age at GUBB draw")
plot_OR_forest_plots(results_kidney, figure_dir, "CH_and_age_OR_kidney.pdf", plot_title = "Age at blood draw (K)") # plot the ORs
results_kidney = calculate_OR(base_kidney_chip, kidney_clinical, gene_list, "Age at GUBB draw", vaf_cutoff = 1)
plot_OR_forest_plots(results_kidney, figure_dir, "CH_and_age_OR_kidney_1_perc.pdf", plot_title = "Age at blood draw (K)") # plot the ORs

results_bladder = calculate_OR(base_bladder_chip, bladder_clinical, gene_list, "Age at baseline blood draw")
plot_OR_forest_plots(results_bladder, figure_dir, "CH_and_age_OR_bladder.pdf", plot_title = "Age at blood draw (B)") # plot the ORs
results_bladder = calculate_OR(base_bladder_chip, bladder_clinical, gene_list, "Age at baseline blood draw", vaf_cutoff = 1)
plot_OR_forest_plots(results_bladder, figure_dir, "CH_and_age_OR_bladder_1_perc.pdf", plot_title = "Age at blood draw (B)") # plot the ORs

combined_age_df = pd.concat([kidney_clinical[["Patient_id", "Age at GUBB draw"]], bladder_clinical[["Patient_id", "Age at baseline blood draw"]].rename(columns = {"Age at baseline blood draw": "Age at GUBB draw"})])
results_combined = calculate_OR(all_vars_chip, combined_age_df, gene_list, "Age at GUBB draw", vaf_cutoff = 1)
plot_OR_forest_plots(results_combined, figure_dir, "CH_and_age_OR_combined_1_perc.pdf", plot_title = "Age at blood draw (B+K)") # plot the ORs

# CH and treatment history
PATH_treatment_landscape_bladder = os.path.join(DIR_working, "resources/clinical_data/bladder/treatment.csv")
bladder_treatment = pd.read_csv(PATH_treatment_landscape_bladder)
bladder_treatment = bladder_treatment[bladder_treatment["Treatment"] != "Palliative / best supportive care"].reset_index(drop = True)[["Patient_id", "Treatment"]]
bladder_treatment_or = calculate_OR_two_categorical(base_bladder_chip, bladder_treatment, gene_list, "Age at GUBB draw", vaf_cutoff = 1)
plot_OR_forest_plots(results_combined, figure_dir, "CH_and_age_OR_combined_1_perc.pdf", plot_title = "Age at blood draw (B+K)") # plot the ORs





