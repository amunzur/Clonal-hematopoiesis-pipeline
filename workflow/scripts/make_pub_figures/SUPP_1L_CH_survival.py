"""
Makes a supp figure looking into the impact of ctDNA mutations on overall survival, subsetting to patients starting 1L treatment
"""

import pandas as pd
import numpy as np
import os 
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from itertools import cycle
# from scipy.stats import ttest_ind
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
from lifelines import CoxPHFitter
from lifelines.plotting import add_at_risk_counts

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
PATH_clinical_bladder = os.path.join(DIR_working, "resources/clinical_data/bladder/clinical_data.csv")
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

with open("/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/visualization/survival_analysis/kidney_functions.py", 'r') as file:
    script_code = file.read()

exec(script_code)

# LOAD CHIP DATASETS
all_vars_chip = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/chip.csv"))
all_vars_chip = all_vars_chip[all_vars_chip["Dependent"] == False]
base_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)

all_pts = pd.read_csv(PATH_sample_information, sep="\t", names=["Patient_id", "Date", "Diagnosis", "Timepoint"])
kidney_pts = all_pts[all_pts["Diagnosis"] == "Kidney"]["Patient_id"].unique()
bladder_pts = all_pts[all_pts["Diagnosis"] == "Bladder"]["Patient_id"].unique()

all_vars_somatic = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/somatic.csv"))
all_vars_somatic = all_vars_somatic[all_vars_somatic["Dependent"] == False]
all_vars_somatic = all_vars_somatic[~all_vars_somatic["Patient_id"].isin(["20-313", "21-184", "21-430"])] # exclude some samples due to oxidative damage
base_kidney_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Bladder")].reset_index(drop = True)

############################ SET UP FIGURE
fig = plt.figure(figsize=(8, 9))
gs = gridspec.GridSpec(3, 2, width_ratios = [1, 1], height_ratios = [1, 1, 1], hspace = 0, wspace = 0)

# CH KMs
CHIP_OS_bladder = plt.subplot(gs[0, 0])
CHIP_OS_kidney = plt.subplot(gs[0, 1])

CH_1L_OS_bladder = plt.subplot(gs[1, 0])
CH_1L_OS_kidney = plt.subplot(gs[1, 1])

CHIP_1L_OS_bladder = plt.subplot(gs[2, 0])
CHIP_1L_OS_kidney = plt.subplot(gs[2, 1])

# CHIP_1L_PFS_bladder = plt.subplot(gs[4, 0])
# CHIP_1L_PFS_kidney = plt.subplot(gs[4, 1])

# Step 1. Prepare the surv DF and the multivar df for bladder.
clin_df_bladder = pd.read_csv(PATH_clinical_bladder)
clin_df_bladder = clin_df_bladder[clin_df_bladder["First sample?"] == True]

sex_and_age_df_bladder = clin_df_bladder[["Patient_id", "Sex", "Age at blood draw"]]
smoking_status = clin_df_bladder[["Patient_id", "Previous smoking history"]]
ctDNA_status_bladder = prepare_ctDNA_status(base_bladder_somatic, "Bladder", PATH_sample_information)
ch_status_bladder = prepare_chip_status(base_bladder_chip, "Bladder", PATH_sample_information)
chip_status_bladder = prepare_chip_status(base_bladder_chip[base_bladder_chip["VAF_n"] >= 2], "Bladder", PATH_sample_information)
surv_df_bladder = prepare_survival_data(clin_df_bladder)
surv_df_bladder["Death"] = surv_df_bladder["Death"].map({True: True, False: False, 'True': True, 'False': False, 'Lost to follow-up': False}) # Map values to boolean, treating 'Lost to follow-up' as False

# Step 2. Prepare the surv DF and the multivar df for kidney.
clin_df_kidney = pd.read_csv(PATH_kidney_clinical)
clin_df_kidney = clin_df_kidney.rename(columns = {"GUBB ID": "Patient_id"})
clin_df_kidney["Date of GUBB draw"] = pd.to_datetime(clin_df_kidney["Date of GUBB draw"])
clin_df_kidney["Date of start of 1L systemic therapy"] = pd.to_datetime(clin_df_kidney["Date of start of 1L systemic therapy"])
clin_df_kidney["GUBB draw minus date start"] = abs(clin_df_kidney["Date of GUBB draw"] - clin_df_kidney["Date of start of 1L systemic therapy"]).dt.days

sex_and_age_df_kidney, cci_df, mets_df, imdc_df, subtype_df, surv_df_kidney = prepare_clinical_data(clin_df_kidney)
surv_df_kidney = surv_df_kidney[surv_df_kidney["Patient_id"].isin(kidney_pts)]
ctDNA_status_kidney = prepare_ctDNA_status(base_kidney_somatic, "Kidney", PATH_sample_information)
ch_status_kidney = prepare_chip_status(base_kidney_chip, "Kidney", PATH_sample_information)
chip_status_kidney = prepare_chip_status(base_kidney_chip[base_kidney_chip["VAF_n"] >=2], "Kidney", PATH_sample_information)

###########################################################################
# PANEL A. CH 1L OS - mUC
pts_1L = clin_df_bladder[(clin_df_bladder["Event"] == "1L") & (clin_df_bladder["Patient_id"] != "23-088")]["Patient_id"].reset_index(drop = True)
multivars_df_binarized_bladder = prepare_multivars_df_bladder(sex_and_age_df_bladder, smoking_status, ctDNA_status_bladder, ch_status_bladder)
df = surv_df_bladder.merge(multivars_df_binarized_bladder, how = "left")
pts_1L_bladder = clin_df_bladder[(clin_df_bladder["Event"] == "1L") & (clin_df_bladder["Patient_id"] != "23-088")]["Patient_id"].reset_index(drop = True)
df = df[df["Patient_id"].isin(pts_1L_bladder)].reset_index(drop = True)
del df["Patient_id"]
del df["ctDNA+"]

CH_1L_OS_bladder = make_survival_curve(
    df, 
    "CH+", 
    event_col = "Death",
    duration_col = "OS from cfDNA collection (mo)",
    output_path = None, 
    plot_title = "", 
    print_HRs = False, 
    add_legend = False,
    pos_curve_label_positions = [0.2, 0.4], 
    neg_curve_label_positions =  [0.4, 0.7], 
    print_survival = False,
    xlabel = "OS from cfDNA collection (mo)",
    ax = CH_1L_OS_bladder)

# Calculate median survival to print on the plot
muc_ctDNA_dict = return_median_survival(df, stratify_by="CH+", event_col="Death", duration_col="OS from cfDNA collection (mo)")
muc_cox_dict = run_cox_proportional_hazards(df, stratify_by="CH+", event_col="Death", duration_col="OS from cfDNA collection (mo)")

CH_1L_OS_bladder.set_title("mUC")
CH_1L_OS_bladder.text(0.93, 1, "Median", transform=CH_1L_OS_bladder.transAxes, fontsize=8, ha='center', va='center')
CH_1L_OS_bladder.text(0.93, 0.93, "mo", transform=CH_1L_OS_bladder.transAxes, fontsize=8, ha='center', va='center', fontstyle = "italic")
CH_1L_OS_bladder.text(0.90, 0.85, 'CH(+)', transform=CH_1L_OS_bladder.transAxes, fontsize=8, ha='right', va='center')
CH_1L_OS_bladder.text(0.90, 0.78, 'CH(-)', transform=CH_1L_OS_bladder.transAxes, fontsize=8, ha='right', va='center')
CH_1L_OS_bladder.text(0.95, 0.85, muc_ctDNA_dict["Positive"], transform=CH_1L_OS_bladder.transAxes, fontsize=8, ha='center', va='center')
CH_1L_OS_bladder.text(0.95, 0.78, muc_ctDNA_dict["Negative"], transform=CH_1L_OS_bladder.transAxes, fontsize=8, ha='center', va='center')

# Cox results
hr = round(muc_cox_dict["HR"], 2)
ci_lower = round(muc_cox_dict["CI_lower"], 2)
ci_upper = round(muc_cox_dict["CI_upper"], 2)
p = round(float(muc_cox_dict['p']), 3)
CH_1L_OS_bladder.text(0.8, 0.6, f"HR {hr} (95%, {ci_lower}-{ci_upper})", transform=CH_1L_OS_bladder.transAxes, fontsize=8, ha='center', va='center')
CH_1L_OS_bladder.text(1.02, 0.53, f"p={p}", transform=CH_1L_OS_bladder.transAxes, fontsize=8, ha='right', va='center')

###########################################################################
# PANEL B. CH 1L OS - mRCC
multivars_df_binarized_kidney = prepare_multivars_df(sex_and_age_df_kidney, cci_df, mets_df, ctDNA_status_kidney, subtype_df, ch_status_kidney, imdc_df)
df = surv_df_kidney.merge(multivars_df_binarized_kidney, how = "left")
pts_1L_kidney = clin_df_kidney[clin_df_kidney["GUBB draw minus date start"] < 21]["Patient_id"].reset_index(drop = True)
df = df[df["Patient_id"].isin(pts_1L_kidney)].reset_index(drop = True)
del df["Patient_id"]
del df["ctDNA+"]

CH_1L_OS_kidney = make_survival_curve(
    df, 
    "CH+",
    event_col = "Death",
    duration_col = "OS from cfDNA collection (mo)",
    output_path = None,
    plot_title = "", 
    print_HRs = True,
    pos_curve_label_positions = [0.5, 0.6], 
    neg_curve_label_positions =  [0.27, 0.35], 
    add_legend = False,
    print_survival = False,
    xlabel = "OS from cfDNA collection (mo)",
    ax = CH_1L_OS_kidney)

CH_1L_OS_kidney.set_title("mRCC")

mrcc_ctDNA_dict = return_median_survival(df, stratify_by="CH+", event_col="Death", duration_col="OS from cfDNA collection (mo)")
mrcc_cox_dict = run_cox_proportional_hazards(df, duration_col="OS from cfDNA collection (mo)", event_col="Death", stratify_by="CH+", ax = None)

# Log rank test
CH_1L_OS_kidney.set_title("mRCC")
CH_1L_OS_kidney.text(0.93, 1, "Median", transform=CH_1L_OS_kidney.transAxes, fontsize=8, ha='center', va='center')
CH_1L_OS_kidney.text(0.93, 0.93, "mo", transform=CH_1L_OS_kidney.transAxes, fontsize=8, ha='center', va='center', fontstyle = "italic")
CH_1L_OS_kidney.text(0.90, 0.85, 'CH(+)', transform=CH_1L_OS_kidney.transAxes, fontsize=8, ha='right', va='center')
CH_1L_OS_kidney.text(0.90, 0.78, 'CH(-)', transform=CH_1L_OS_kidney.transAxes, fontsize=8, ha='right', va='center')
CH_1L_OS_kidney.text(0.95, 0.85, mrcc_ctDNA_dict["Positive"], transform=CH_1L_OS_kidney.transAxes, fontsize=8, ha='center', va='center')
CH_1L_OS_kidney.text(0.95, 0.78, mrcc_ctDNA_dict["Negative"], transform=CH_1L_OS_kidney.transAxes, fontsize=8, ha='center', va='center')

# Cox results
hr = round(mrcc_cox_dict["HR"], 2)
ci_lower = round(mrcc_cox_dict["CI_lower"], 2)
ci_upper = round(mrcc_cox_dict["CI_upper"], 2)
p = round(float(mrcc_cox_dict['p']), 3)
CH_1L_OS_kidney.text(0.8, 0.6, f"HR {hr} (95%, {ci_lower}-{ci_upper})", transform=CH_1L_OS_kidney.transAxes, fontsize=8, ha='center', va='center')
CH_1L_OS_kidney.text(1.02, 0.53, f"p={p}", transform=CH_1L_OS_kidney.transAxes, fontsize=8, ha='right', va='center')

# ###########################################################################
# # PANEL C. CH 1L PFS - mUC
# pts_1L = clin_df_bladder[(clin_df_bladder["Event"] == "1L") & (clin_df_bladder["Patient_id"] != "23-088")]["Patient_id"].reset_index(drop = True)
# multivars_df_binarized_bladder = prepare_multivars_df_bladder(sex_and_age_df_bladder, smoking_status, ctDNA_status_bladder, ch_status_bladder)
# df = surv_df_bladder.merge(multivars_df_binarized_bladder, how = "left")
# df = df[df["Patient_id"].isin(pts_1L_bladder)].reset_index(drop = True)
# df = df[~df["PFS (mo)"].isna()]
# del df["Patient_id"]
# del df["ctDNA+"]

# CH_1L_PFS_bladder = make_survival_curve(
#     df, 
#     "CH+", 
#     event_col = "Progression",
#     duration_col = "PFS (mo)",
#     output_path = None, 
#     plot_title = "", 
#     print_HRs = False, 
#     add_legend = False,
#     pos_curve_label_positions = [0.2, 0.4], 
#     neg_curve_label_positions =  [0.4, 0.7], 
#     print_survival = False,
#     xlabel = "PFS (mo)",
#     ax = CH_1L_PFS_bladder)

# # Calculate median survival to print on the plot
# muc_ctDNA_dict = return_median_survival(df, stratify_by="CH+", event_col="Progression", duration_col="PFS (mo)")
# muc_cox_dict = run_cox_proportional_hazards(df, stratify_by="CH+", event_col="Progression", duration_col="PFS (mo)")

# CH_1L_PFS_bladder.set_title("mUC")
# CH_1L_PFS_bladder.text(0.93, 1, "Median", transform=CH_1L_PFS_bladder.transAxes, fontsize=8, ha='center', va='center')
# CH_1L_PFS_bladder.text(0.93, 0.93, "mo", transform=CH_1L_PFS_bladder.transAxes, fontsize=8, ha='center', va='center', fontstyle = "italic")
# CH_1L_PFS_bladder.text(0.90, 0.85, 'CH(+)', transform=CH_1L_PFS_bladder.transAxes, fontsize=8, ha='right', va='center')
# CH_1L_PFS_bladder.text(0.90, 0.78, 'CH(-)', transform=CH_1L_PFS_bladder.transAxes, fontsize=8, ha='right', va='center')
# CH_1L_PFS_bladder.text(0.95, 0.85, muc_ctDNA_dict["Positive"], transform=CH_1L_PFS_bladder.transAxes, fontsize=8, ha='center', va='center')
# CH_1L_PFS_bladder.text(0.95, 0.78, muc_ctDNA_dict["Negative"], transform=CH_1L_PFS_bladder.transAxes, fontsize=8, ha='center', va='center')

# # Cox results
# hr = round(muc_cox_dict["HR"], 2)
# ci_lower = round(muc_cox_dict["CI_lower"], 2)
# ci_upper = round(muc_cox_dict["CI_upper"], 2)
# p = round(float(muc_cox_dict['p']), 3)
# CH_1L_PFS_bladder.text(0.8, 0.6, f"HR {hr} (95%, {ci_lower}-{ci_upper})", transform=CH_1L_PFS_bladder.transAxes, fontsize=8, ha='center', va='center')
# CH_1L_PFS_bladder.text(1.02, 0.53, f"p={p}", transform=CH_1L_PFS_bladder.transAxes, fontsize=8, ha='right', va='center')

# ###########################################################################
# # PANEL D. CH 1L PFS - mRCC
# multivars_df_binarized_kidney = prepare_multivars_df(sex_and_age_df_kidney, cci_df, mets_df, ctDNA_status_kidney, subtype_df, ch_status_kidney, imdc_df)
# df = surv_df_kidney.merge(multivars_df_binarized_kidney, how = "left")
# pts_1L_kidney = clin_df_kidney[clin_df_kidney["GUBB draw minus date start"] < 21]["Patient_id"].reset_index(drop = True)
# df = df[df["Patient_id"].isin(pts_1L_kidney)].reset_index(drop = True)
# df = df[~df["PFS (mo)"].isna()]
# del df["Patient_id"]
# del df["ctDNA+"]

# CH_1L_PFS_kidney = make_survival_curve(
#     df, 
#     "CH+",
#     event_col = "Progression",
#     duration_col = "PFS (mo)",
#     output_path = None,
#     plot_title = "", 
#     print_HRs = True,
#     pos_curve_label_positions = [0.5, 0.6], 
#     neg_curve_label_positions =  [0.27, 0.35], 
#     add_legend = False,
#     print_survival = False,
#     xlabel = "PFS (mo)",
#     ax = CH_1L_PFS_kidney)

# CH_1L_PFS_kidney.set_title("mRCC")

# mrcc_ctDNA_dict = return_median_survival(df, stratify_by="CH+", event_col="Progression", duration_col="PFS (mo)")
# mrcc_cox_dict = run_cox_proportional_hazards(df, duration_col="PFS (mo)", event_col="Progression", stratify_by="CH+", ax = None)

# # Log rank test
# CH_1L_PFS_kidney.set_title("mRCC")
# CH_1L_PFS_kidney.text(0.93, 1, "Median", transform=CH_1L_PFS_kidney.transAxes, fontsize=8, ha='center', va='center')
# CH_1L_PFS_kidney.text(0.93, 0.93, "mo", transform=CH_1L_PFS_kidney.transAxes, fontsize=8, ha='center', va='center', fontstyle = "italic")
# CH_1L_PFS_kidney.text(0.90, 0.85, 'CH(+)', transform=CH_1L_PFS_kidney.transAxes, fontsize=8, ha='right', va='center')
# CH_1L_PFS_kidney.text(0.90, 0.78, 'CH(-)', transform=CH_1L_PFS_kidney.transAxes, fontsize=8, ha='right', va='center')
# CH_1L_PFS_kidney.text(0.95, 0.85, mrcc_ctDNA_dict["Positive"], transform=CH_1L_PFS_kidney.transAxes, fontsize=8, ha='center', va='center')
# CH_1L_PFS_kidney.text(0.95, 0.78, mrcc_ctDNA_dict["Negative"], transform=CH_1L_PFS_kidney.transAxes, fontsize=8, ha='center', va='center')

# # Cox results
# hr = round(mrcc_cox_dict["HR"], 2)
# ci_lower = round(mrcc_cox_dict["CI_lower"], 2)
# ci_upper = round(mrcc_cox_dict["CI_upper"], 2)
# p = round(float(mrcc_cox_dict['p']), 3)
# CH_1L_PFS_kidney.text(0.8, 0.6, f"HR {hr} (95%, {ci_lower}-{ci_upper})", transform=CH_1L_PFS_kidney.transAxes, fontsize=8, ha='center', va='center')
# CH_1L_PFS_kidney.text(1.02, 0.53, f"p={p}", transform=CH_1L_PFS_kidney.transAxes, fontsize=8, ha='right', va='center')



























# CHIP PLOTS
###########################################################################
# PANEL E. CHIP OS - mUC
multivars_df_binarized_bladder = prepare_multivars_df_bladder(sex_and_age_df_bladder, smoking_status, ctDNA_status_bladder, chip_status_bladder)
df = surv_df_bladder.merge(multivars_df_binarized_bladder, how = "left")
del df["Patient_id"]
del df["ctDNA+"]

CHIP_OS_bladder = make_survival_curve(
    df, 
    "CH+", 
    event_col = "Death",
    duration_col = "OS from cfDNA collection (mo)",
    output_path = None, 
    plot_title = "", 
    print_HRs = False, 
    add_legend = False,
    pos_curve_label_positions = [0.2, 0.4], 
    neg_curve_label_positions =  [0.4, 0.7], 
    print_survival = False,
    xlabel = "OS from cfDNA collection (mo)",
    ax = CHIP_OS_bladder)

# Calculate median survival to print on the plot
muc_ctDNA_dict = return_median_survival(df, stratify_by="CH+", event_col="Death", duration_col="OS from cfDNA collection (mo)")
muc_cox_dict = run_cox_proportional_hazards(df, stratify_by="CH+", event_col="Death", duration_col="OS from cfDNA collection (mo)")

CHIP_OS_bladder.set_title("mUC")
CHIP_OS_bladder.text(0.93, 1, "Median", transform=CHIP_OS_bladder.transAxes, fontsize=8, ha='center', va='center')
CHIP_OS_bladder.text(0.93, 0.93, "mo", transform=CHIP_OS_bladder.transAxes, fontsize=8, ha='center', va='center', fontstyle = "italic")
CHIP_OS_bladder.text(0.90, 0.85, 'CH(+)', transform=CHIP_OS_bladder.transAxes, fontsize=8, ha='right', va='center')
CHIP_OS_bladder.text(0.90, 0.78, 'CH(-)', transform=CHIP_OS_bladder.transAxes, fontsize=8, ha='right', va='center')
CHIP_OS_bladder.text(0.95, 0.85, muc_ctDNA_dict["Positive"], transform=CHIP_OS_bladder.transAxes, fontsize=8, ha='center', va='center')
CHIP_OS_bladder.text(0.95, 0.78, muc_ctDNA_dict["Negative"], transform=CHIP_OS_bladder.transAxes, fontsize=8, ha='center', va='center')

# Cox results
hr = round(muc_cox_dict["HR"], 2)
ci_lower = round(muc_cox_dict["CI_lower"], 2)
ci_upper = round(muc_cox_dict["CI_upper"], 2)
p = round(float(muc_cox_dict['p']), 3)
CHIP_OS_bladder.text(0.8, 0.6, f"HR {hr} (95%, {ci_lower}-{ci_upper})", transform=CHIP_OS_bladder.transAxes, fontsize=8, ha='center', va='center')
CHIP_OS_bladder.text(1.02, 0.53, f"p={p}", transform=CHIP_OS_bladder.transAxes, fontsize=8, ha='right', va='center')

###########################################################################
# PANEL F. CHIP OS RCC
multivars_df_binarized_kidney = prepare_multivars_df(sex_and_age_df_kidney, cci_df, mets_df, ctDNA_status_kidney, subtype_df, chip_status_kidney, imdc_df)
df = surv_df_kidney.merge(multivars_df_binarized_kidney, how = "left")
del df["Patient_id"]
del df["ctDNA+"]

CHIP_OS_kidney = make_survival_curve(
    df, 
    "CH+",
    event_col = "Death",
    duration_col = "OS from cfDNA collection (mo)",
    output_path = None,
    plot_title = "", 
    print_HRs = True,
    pos_curve_label_positions = [0.5, 0.6], 
    neg_curve_label_positions =  [0.27, 0.35], 
    add_legend = False,
    print_survival = False,
    xlabel = "OS from cfDNA collection (mo)",
    ax = CHIP_OS_kidney)

CHIP_OS_kidney.set_title("mRCC")

mrcc_ctDNA_dict = return_median_survival(df, stratify_by="CH+", event_col="Death", duration_col="OS from cfDNA collection (mo)")
mrcc_cox_dict = run_cox_proportional_hazards(df, duration_col="OS from cfDNA collection (mo)", event_col="Death", stratify_by="CH+", ax = None)

# Log rank test
CHIP_OS_kidney.set_title("mRCC")
CHIP_OS_kidney.text(0.93, 1, "Median", transform=CHIP_OS_kidney.transAxes, fontsize=8, ha='center', va='center')
CHIP_OS_kidney.text(0.93, 0.93, "mo", transform=CHIP_OS_kidney.transAxes, fontsize=8, ha='center', va='center', fontstyle = "italic")
CHIP_OS_kidney.text(0.90, 0.85, 'CH(+)', transform=CHIP_OS_kidney.transAxes, fontsize=8, ha='right', va='center')
CHIP_OS_kidney.text(0.90, 0.78, 'CH(-)', transform=CHIP_OS_kidney.transAxes, fontsize=8, ha='right', va='center')
CHIP_OS_kidney.text(0.95, 0.85, mrcc_ctDNA_dict["Positive"], transform=CHIP_OS_kidney.transAxes, fontsize=8, ha='center', va='center')
CHIP_OS_kidney.text(0.95, 0.78, mrcc_ctDNA_dict["Negative"], transform=CHIP_OS_kidney.transAxes, fontsize=8, ha='center', va='center')

# Cox results
hr = round(mrcc_cox_dict["HR"], 2)
ci_lower = round(mrcc_cox_dict["CI_lower"], 2)
ci_upper = round(mrcc_cox_dict["CI_upper"], 2)
p = round(float(mrcc_cox_dict['p']), 3)
CHIP_OS_kidney.text(0.8, 0.6, f"HR {hr} (95%, {ci_lower}-{ci_upper})", transform=CHIP_OS_kidney.transAxes, fontsize=8, ha='center', va='center')
CHIP_OS_kidney.text(1.02, 0.53, f"p={p}", transform=CHIP_OS_kidney.transAxes, fontsize=8, ha='right', va='center')








###########################################################################
# PANEL E. CHIP OS - mUC in 1L
pts_1L = clin_df_bladder[(clin_df_bladder["Event"] == "1L") & (clin_df_bladder["Patient_id"] != "23-088")]["Patient_id"].reset_index(drop = True)
multivars_df_binarized_bladder = prepare_multivars_df_bladder(sex_and_age_df_bladder, smoking_status, ctDNA_status_bladder, chip_status_bladder)
df = surv_df_bladder.merge(multivars_df_binarized_bladder, how = "left")
df = df[df["Patient_id"].isin(pts_1L_bladder)].reset_index(drop = True)
del df["Patient_id"]
del df["ctDNA+"]

CHIP_1L_OS_bladder = make_survival_curve(
    df, 
    "CH+", 
    event_col = "Death",
    duration_col = "OS from cfDNA collection (mo)",
    output_path = None, 
    plot_title = "", 
    print_HRs = False, 
    add_legend = False,
    pos_curve_label_positions = [0.2, 0.4], 
    neg_curve_label_positions =  [0.4, 0.7], 
    print_survival = False,
    xlabel = "OS from cfDNA collection (mo)",
    ax = CHIP_1L_OS_bladder)

# Calculate median survival to print on the plot
muc_ctDNA_dict = return_median_survival(df, stratify_by="CH+", event_col="Death", duration_col="OS from cfDNA collection (mo)")
muc_cox_dict = run_cox_proportional_hazards(df, stratify_by="CH+", event_col="Death", duration_col="OS from cfDNA collection (mo)")

CHIP_1L_OS_bladder.set_title("mUC")
CHIP_1L_OS_bladder.text(0.93, 1, "Median", transform=CHIP_1L_OS_bladder.transAxes, fontsize=8, ha='center', va='center')
CHIP_1L_OS_bladder.text(0.93, 0.93, "mo", transform=CHIP_1L_OS_bladder.transAxes, fontsize=8, ha='center', va='center', fontstyle = "italic")
CHIP_1L_OS_bladder.text(0.90, 0.85, 'CH(+)', transform=CHIP_1L_OS_bladder.transAxes, fontsize=8, ha='right', va='center')
CHIP_1L_OS_bladder.text(0.90, 0.78, 'CH(-)', transform=CHIP_1L_OS_bladder.transAxes, fontsize=8, ha='right', va='center')
CHIP_1L_OS_bladder.text(0.95, 0.85, muc_ctDNA_dict["Positive"], transform=CHIP_1L_OS_bladder.transAxes, fontsize=8, ha='center', va='center')
CHIP_1L_OS_bladder.text(0.95, 0.78, muc_ctDNA_dict["Negative"], transform=CHIP_1L_OS_bladder.transAxes, fontsize=8, ha='center', va='center')

# Cox results
hr = round(muc_cox_dict["HR"], 2)
ci_lower = round(muc_cox_dict["CI_lower"], 2)
ci_upper = round(muc_cox_dict["CI_upper"], 2)
p = round(float(muc_cox_dict['p']), 3)
CHIP_1L_OS_bladder.text(0.8, 0.6, f"HR {hr} (95%, {ci_lower}-{ci_upper})", transform=CHIP_1L_OS_bladder.transAxes, fontsize=8, ha='center', va='center')
CHIP_1L_OS_bladder.text(1.02, 0.53, f"p={p}", transform=CHIP_1L_OS_bladder.transAxes, fontsize=8, ha='right', va='center')

###########################################################################
# PANEL F. CHIP OS RCC 1L
multivars_df_binarized_kidney = prepare_multivars_df(sex_and_age_df_kidney, cci_df, mets_df, ctDNA_status_kidney, subtype_df, chip_status_kidney, imdc_df)
df = surv_df_kidney.merge(multivars_df_binarized_kidney, how = "left")
pts_1L_kidney = clin_df_kidney[clin_df_kidney["GUBB draw minus date start"] < 21]["Patient_id"].reset_index(drop = True)
df = df[df["Patient_id"].isin(pts_1L_kidney)].reset_index(drop = True)
del df["Patient_id"]
del df["ctDNA+"]

CHIP_1L_OS_kidney = make_survival_curve(
    df, 
    "CH+",
    event_col = "Death",
    duration_col = "OS from cfDNA collection (mo)",
    output_path = None,
    plot_title = "", 
    print_HRs = True,
    pos_curve_label_positions = [0.5, 0.6], 
    neg_curve_label_positions =  [0.27, 0.35], 
    add_legend = False,
    print_survival = False,
    xlabel = "OS from cfDNA collection (mo)",
    ax = CHIP_1L_OS_kidney)

CHIP_1L_OS_kidney.set_title("mRCC")

mrcc_ctDNA_dict = return_median_survival(df, stratify_by="CH+", event_col="Death", duration_col="OS from cfDNA collection (mo)")
mrcc_cox_dict = run_cox_proportional_hazards(df, duration_col="OS from cfDNA collection (mo)", event_col="Death", stratify_by="CH+", ax = None)

# Log rank test
CHIP_1L_OS_kidney.set_title("mRCC")
CHIP_1L_OS_kidney.text(0.93, 1, "Median", transform=CHIP_1L_OS_kidney.transAxes, fontsize=8, ha='center', va='center')
CHIP_1L_OS_kidney.text(0.93, 0.93, "mo", transform=CHIP_1L_OS_kidney.transAxes, fontsize=8, ha='center', va='center', fontstyle = "italic")
CHIP_1L_OS_kidney.text(0.90, 0.85, 'CH(+)', transform=CHIP_1L_OS_kidney.transAxes, fontsize=8, ha='right', va='center')
CHIP_1L_OS_kidney.text(0.90, 0.78, 'CH(-)', transform=CHIP_1L_OS_kidney.transAxes, fontsize=8, ha='right', va='center')
CHIP_1L_OS_kidney.text(0.95, 0.85, mrcc_ctDNA_dict["Positive"], transform=CHIP_1L_OS_kidney.transAxes, fontsize=8, ha='center', va='center')
CHIP_1L_OS_kidney.text(0.95, 0.78, mrcc_ctDNA_dict["Negative"], transform=CHIP_1L_OS_kidney.transAxes, fontsize=8, ha='center', va='center')

# Cox results
hr = round(mrcc_cox_dict["HR"], 2)
ci_lower = round(mrcc_cox_dict["CI_lower"], 2)
ci_upper = round(mrcc_cox_dict["CI_upper"], 2)
p = round(float(mrcc_cox_dict['p']), 3)
CHIP_1L_OS_kidney.text(0.8, 0.6, f"HR {hr} (95%, {ci_lower}-{ci_upper})", transform=CHIP_1L_OS_kidney.transAxes, fontsize=8, ha='center', va='center')
CHIP_1L_OS_kidney.text(1.02, 0.53, f"p={p}", transform=CHIP_1L_OS_kidney.transAxes, fontsize=8, ha='right', va='center')








# ###########################################################################
# # PANEL CHIP 1L PFS
# ###########################################################################
# # PANEL C. CH 1L PFS - mUC
# pts_1L = clin_df_bladder[(clin_df_bladder["Event"] == "1L") & (clin_df_bladder["Patient_id"] != "23-088")]["Patient_id"].reset_index(drop = True)
# multivars_df_binarized_bladder = prepare_multivars_df_bladder(sex_and_age_df_bladder, smoking_status, ctDNA_status_bladder, chip_status_bladder)
# df = surv_df_bladder.merge(multivars_df_binarized_bladder, how = "left")
# df = df[df["Patient_id"].isin(pts_1L_bladder)].reset_index(drop = True)
# df = df[~df["PFS (mo)"].isna()]
# del df["Patient_id"]
# del df["ctDNA+"]

# CHIP_1L_PFS_bladder = make_survival_curve(
#     df, 
#     "CH+", 
#     event_col = "Progression",
#     duration_col = "PFS (mo)",
#     output_path = None, 
#     plot_title = "", 
#     print_HRs = False, 
#     add_legend = False,
#     pos_curve_label_positions = [0.2, 0.4], 
#     neg_curve_label_positions =  [0.4, 0.7], 
#     print_survival = False,
#     xlabel = "PFS (mo)",
#     ax = CHIP_1L_PFS_bladder)

# # Calculate median survival to print on the plot
# muc_ctDNA_dict = return_median_survival(df, stratify_by="CH+", event_col="Progression", duration_col="PFS (mo)")
# muc_cox_dict = run_cox_proportional_hazards(df, stratify_by="CH+", event_col="Progression", duration_col="PFS (mo)")

# CHIP_1L_PFS_bladder.set_title("mUC")
# CHIP_1L_PFS_bladder.text(0.93, 1, "Median", transform=CHIP_1L_PFS_bladder.transAxes, fontsize=8, ha='center', va='center')
# CHIP_1L_PFS_bladder.text(0.93, 0.93, "mo", transform=CHIP_1L_PFS_bladder.transAxes, fontsize=8, ha='center', va='center', fontstyle = "italic")
# CHIP_1L_PFS_bladder.text(0.90, 0.85, 'CH(+)', transform=CHIP_1L_PFS_bladder.transAxes, fontsize=8, ha='right', va='center')
# CHIP_1L_PFS_bladder.text(0.90, 0.78, 'CH(-)', transform=CHIP_1L_PFS_bladder.transAxes, fontsize=8, ha='right', va='center')
# CHIP_1L_PFS_bladder.text(0.95, 0.85, muc_ctDNA_dict["Positive"], transform=CHIP_1L_PFS_bladder.transAxes, fontsize=8, ha='center', va='center')
# CHIP_1L_PFS_bladder.text(0.95, 0.78, muc_ctDNA_dict["Negative"], transform=CHIP_1L_PFS_bladder.transAxes, fontsize=8, ha='center', va='center')

# # Cox results
# hr = round(muc_cox_dict["HR"], 2)
# ci_lower = round(muc_cox_dict["CI_lower"], 2)
# ci_upper = round(muc_cox_dict["CI_upper"], 2)
# p = round(float(muc_cox_dict['p']), 3)
# CHIP_1L_PFS_bladder.text(0.8, 0.6, f"HR {hr} (95%, {ci_lower}-{ci_upper})", transform=CHIP_1L_PFS_bladder.transAxes, fontsize=8, ha='center', va='center')
# CHIP_1L_PFS_bladder.text(1.02, 0.53, f"p={p}", transform=CHIP_1L_PFS_bladder.transAxes, fontsize=8, ha='right', va='center')

# ###########################################################################
# # PANEL D. CH 1L PFS - mRCC
# multivars_df_binarized_kidney = prepare_multivars_df(sex_and_age_df_kidney, cci_df, mets_df, ctDNA_status_kidney, subtype_df, chip_status_kidney, imdc_df)
# df = surv_df_kidney.merge(multivars_df_binarized_kidney, how = "left")
# pts_1L_kidney = clin_df_kidney[clin_df_kidney["GUBB draw minus date start"] < 21]["Patient_id"].reset_index(drop = True)
# df = df[df["Patient_id"].isin(pts_1L_kidney)].reset_index(drop = True)
# df = df[~df["PFS (mo)"].isna()]
# del df["Patient_id"]
# del df["ctDNA+"]

# CHIP_1L_PFS_kidney = make_survival_curve(
#     df, 
#     "CH+",
#     event_col = "Progression",
#     duration_col = "PFS (mo)",
#     output_path = None,
#     plot_title = "", 
#     print_HRs = True,
#     pos_curve_label_positions = [0.5, 0.6], 
#     neg_curve_label_positions =  [0.27, 0.35], 
#     add_legend = False,
#     print_survival = False,
#     xlabel = "PFS (mo)",
#     ax = CHIP_1L_PFS_kidney)

# CHIP_1L_PFS_kidney.set_title("mRCC")

# mrcc_ctDNA_dict = return_median_survival(df, stratify_by="CH+", event_col="Progression", duration_col="PFS (mo)")
# mrcc_cox_dict = run_cox_proportional_hazards(df, duration_col="PFS (mo)", event_col="Progression", stratify_by="CH+", ax = None)

# # Log rank test
# CHIP_1L_PFS_kidney.set_title("mRCC")
# CHIP_1L_PFS_kidney.text(0.93, 1, "Median", transform=CHIP_1L_PFS_kidney.transAxes, fontsize=8, ha='center', va='center')
# CHIP_1L_PFS_kidney.text(0.93, 0.93, "mo", transform=CHIP_1L_PFS_kidney.transAxes, fontsize=8, ha='center', va='center', fontstyle = "italic")
# CHIP_1L_PFS_kidney.text(0.90, 0.85, 'CH(+)', transform=CHIP_1L_PFS_kidney.transAxes, fontsize=8, ha='right', va='center')
# CHIP_1L_PFS_kidney.text(0.90, 0.78, 'CH(-)', transform=CHIP_1L_PFS_kidney.transAxes, fontsize=8, ha='right', va='center')
# CHIP_1L_PFS_kidney.text(0.95, 0.85, mrcc_ctDNA_dict["Positive"], transform=CHIP_1L_PFS_kidney.transAxes, fontsize=8, ha='center', va='center')
# CHIP_1L_PFS_kidney.text(0.95, 0.78, mrcc_ctDNA_dict["Negative"], transform=CHIP_1L_PFS_kidney.transAxes, fontsize=8, ha='center', va='center')

# # Cox results
# hr = round(mrcc_cox_dict["HR"], 2)
# ci_lower = round(mrcc_cox_dict["CI_lower"], 2)
# ci_upper = round(mrcc_cox_dict["CI_upper"], 2)
# p = round(float(mrcc_cox_dict['p']), 3)
# CHIP_1L_PFS_kidney.text(0.8, 0.6, f"HR {hr} (95%, {ci_lower}-{ci_upper})", transform=CHIP_1L_PFS_kidney.transAxes, fontsize=8, ha='center', va='center')
# CHIP_1L_PFS_kidney.text(1.02, 0.53, f"p={p}", transform=CHIP_1L_PFS_kidney.transAxes, fontsize=8, ha='right', va='center')

gs.tight_layout(fig)
fig.savefig(os.path.join(figure_dir, "SUPP_CH_KM.png"))
fig.savefig(os.path.join(figure_dir, "SUPP_CH_KM.pdf"))