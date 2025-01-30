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

# LOAD SOMATIC DATASETS
all_vars_somatic = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/somatic.csv"))
all_vars_somatic = all_vars_somatic[all_vars_somatic["Dependent"] == False]
all_vars_somatic = all_vars_somatic[~all_vars_somatic["Patient_id"].isin(["20-313", "21-184", "21-430"])] # exclude some samples due to oxidative damage
base_kidney_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Bladder")].reset_index(drop = True)

all_pts = pd.read_csv(PATH_sample_information, sep="\t", names=["Patient_id", "Date", "Diagnosis", "Timepoint"])
kidney_pts = all_pts[all_pts["Diagnosis"] == "Kidney"]["Patient_id"].unique()
bladder_pts = all_pts[all_pts["Diagnosis"] == "Bladder"]["Patient_id"].unique()

############################ SET UP FIGURE
fig = plt.figure(figsize=(8, 8))
gs = gridspec.GridSpec(3, 2, width_ratios = [1, 1], height_ratios = [1, 1, 1], hspace = 0, wspace = 0)

# ctDNA KMs
os_1L_bladder_ax = plt.subplot(gs[0, 0])
os_1L_kidney_ax = plt.subplot(gs[0, 1])

os_with_ch_bladder_ax = plt.subplot(gs[1, 0])
os_with_ch_kidney_ax = plt.subplot(gs[1, 1])

# Step 1. Prepare the surv DF and the multivar df for bladder.
clin_df_bladder = pd.read_csv(PATH_clinical_bladder)
clin_df_bladder = clin_df_bladder[clin_df_bladder["First sample?"] == True]

sex_and_age_df_bladder = clin_df_bladder[["Patient_id", "Sex", "Age at blood draw"]]
smoking_status = clin_df_bladder[["Patient_id", "Previous smoking history"]]
ctDNA_status_bladder = prepare_ctDNA_status(base_bladder_somatic, "Bladder", PATH_sample_information)
chip_status_bladder = prepare_chip_status(base_bladder_chip, "Bladder", PATH_sample_information)
surv_df_bladder = prepare_survival_data(clin_df_bladder)
surv_df_bladder["Death"] = surv_df_bladder["Death"].map({True: True, False: False, 'True': True, 'False': False, 'Lost to follow-up': False}) # Map values to boolean, treating 'Lost to follow-up' as False
multivars_df_binarized_bladder = prepare_multivars_df_bladder(sex_and_age_df_bladder, smoking_status, ctDNA_status_bladder, chip_status_bladder)

# Step 2. Prepare the surv DF and the multivar df for kidney.
clin_df_kidney = pd.read_csv(PATH_kidney_clinical)
clin_df_kidney = clin_df_kidney.rename(columns = {"GUBB ID": "Patient_id"})
clin_df_kidney["Date of GUBB draw"] = pd.to_datetime(clin_df_kidney["Date of GUBB draw"])
clin_df_kidney["Date of start of 1L systemic therapy"] = pd.to_datetime(clin_df_kidney["Date of start of 1L systemic therapy"])
clin_df_kidney["GUBB draw minus date start"] = abs(clin_df_kidney["Date of GUBB draw"] - clin_df_kidney["Date of start of 1L systemic therapy"]).dt.days

sex_and_age_df_kidney, cci_df, mets_df, imdc_df, subtype_df, surv_df_kidney = prepare_clinical_data(clin_df_kidney)
surv_df_kidney = surv_df_kidney[surv_df_kidney["Patient_id"].isin(kidney_pts)]
ctDNA_status_kidney = prepare_ctDNA_status(base_kidney_somatic, "Kidney", PATH_sample_information)
chip_status_kidney = prepare_chip_status(base_kidney_chip, "Kidney", PATH_sample_information)
multivars_df_binarized_kidney = prepare_multivars_df(sex_and_age_df_kidney, cci_df, mets_df, ctDNA_status_kidney, subtype_df, chip_status_kidney, imdc_df)

###########################################################################
# Figure 1A. mUC ctDNA OS in 1L patients
df = surv_df_bladder.merge(multivars_df_binarized_bladder, how = "left")
pts_1L_bladder = clin_df_bladder[(clin_df_bladder["Event"] == "1L") & (clin_df_bladder["Patient_id"] != "23-088")]["Patient_id"].reset_index(drop = True)
df = df[df["Patient_id"].isin(pts_1L_bladder)].reset_index(drop = True)
del df["Patient_id"]
del df["CH+"]

os_1L_bladder_ax = make_survival_curve(
    df, 
    "ctDNA+", 
    event_col = "Death",
    duration_col = "OS from cfDNA collection (mo)",
    output_path = None, 
    plot_title = "", 
    print_HRs = False, 
    add_legend = False,
    pos_curve_label_positions = [0.2, 0.4], 
    neg_curve_label_positions =  [0.45, 0.7], 
    print_survival = False,
    xlabel = "OS from cfDNA collection (mo)",
    ax = os_1L_bladder_ax)

# Calculate median survival to print on the plot
muc_ctDNA_dict = return_median_survival(df, stratify_by="ctDNA+", event_col="Death", duration_col="OS from cfDNA collection (mo)")
muc_cox_dict = run_cox_proportional_hazards(df, stratify_by="ctDNA+", event_col="Death", duration_col="OS from cfDNA collection (mo)")

os_1L_bladder_ax.set_title("mUC")
os_1L_bladder_ax.text(0.93, 1, "Median", transform=os_1L_bladder_ax.transAxes, fontsize=8, ha='center', va='center')
os_1L_bladder_ax.text(0.93, 0.93, "mo", transform=os_1L_bladder_ax.transAxes, fontsize=8, ha='center', va='center', fontstyle = "italic")
os_1L_bladder_ax.text(0.90, 0.85, 'ctDNA+', transform=os_1L_bladder_ax.transAxes, fontsize=8, ha='right', va='center')
os_1L_bladder_ax.text(0.90, 0.78, 'ctDNA-', transform=os_1L_bladder_ax.transAxes, fontsize=8, ha='right', va='center')
os_1L_bladder_ax.text(0.95, 0.85, muc_ctDNA_dict["Positive"], transform=os_1L_bladder_ax.transAxes, fontsize=8, ha='center', va='center')
os_1L_bladder_ax.text(0.95, 0.78, muc_ctDNA_dict["Negative"], transform=os_1L_bladder_ax.transAxes, fontsize=8, ha='center', va='center')

# Cox results
hr = round(muc_cox_dict["HR"], 2)
ci_lower = round(muc_cox_dict["CI_lower"], 2)
ci_upper = round(muc_cox_dict["CI_upper"], 2)
p = round(float(muc_cox_dict['p']), 3)
os_1L_bladder_ax.text(0.8, 0.6, f"HR {hr} (95%, {ci_lower}-{ci_upper})", transform=os_1L_bladder_ax.transAxes, fontsize=8, ha='center', va='center')
os_1L_bladder_ax.text(1.02, 0.53, f"p={p}", transform=os_1L_bladder_ax.transAxes, fontsize=8, ha='right', va='center')

###########################################################################
# Figure 1B. mRCC ctDNA OS in 1L patients
df = surv_df_kidney.merge(multivars_df_binarized_kidney, how = "left")
pts_1L_kidney = clin_df_kidney[clin_df_kidney["GUBB draw minus date start"] < 21]["Patient_id"].reset_index(drop = True)
df = df[df["Patient_id"].isin(pts_1L_kidney)].reset_index(drop = True)
del df["Patient_id"]
del df["CH+"]

os_1L_kidney_ax = make_survival_curve(
    df, 
    "ctDNA+",
    event_col = "Death",
    duration_col = "OS from cfDNA collection (mo)",
    output_path = None,
    plot_title = "", 
    print_HRs = True,
    pos_curve_label_positions = [0.27, 0.23], 
    neg_curve_label_positions =  [0.5, 0.72], 
    add_legend = False,
    print_survival = False,
    xlabel = "OS from cfDNA collection (mo)",
    ax = os_1L_kidney_ax)

os_1L_kidney_ax.set_title("mRCC")

mrcc_ctDNA_dict = return_median_survival(df, stratify_by="ctDNA+", event_col="Death", duration_col="OS from cfDNA collection (mo)")
mrcc_cox_dict = run_cox_proportional_hazards(df, duration_col="OS from cfDNA collection (mo)", event_col="Death", stratify_by="ctDNA+", ax = None)

# Log rank test
os_1L_kidney_ax.set_title("mRCC")
os_1L_kidney_ax.text(0.93, 1, "Median", transform=os_1L_kidney_ax.transAxes, fontsize=8, ha='center', va='center')
os_1L_kidney_ax.text(0.93, 0.93, "mo", transform=os_1L_kidney_ax.transAxes, fontsize=8, ha='center', va='center', fontstyle = "italic")
os_1L_kidney_ax.text(0.90, 0.85, 'ctDNA+', transform=os_1L_kidney_ax.transAxes, fontsize=8, ha='right', va='center')
os_1L_kidney_ax.text(0.90, 0.78, 'ctDNA-', transform=os_1L_kidney_ax.transAxes, fontsize=8, ha='right', va='center')
os_1L_kidney_ax.text(0.95, 0.85, mrcc_ctDNA_dict["Positive"], transform=os_1L_kidney_ax.transAxes, fontsize=8, ha='center', va='center')
os_1L_kidney_ax.text(0.95, 0.78, mrcc_ctDNA_dict["Negative"], transform=os_1L_kidney_ax.transAxes, fontsize=8, ha='center', va='center')

# Cox results
hr = round(mrcc_cox_dict["HR"], 2)
ci_lower = round(mrcc_cox_dict["CI_lower"], 2)
ci_upper = round(mrcc_cox_dict["CI_upper"], 2)
p = round(float(mrcc_cox_dict['p']), 3)
os_1L_kidney_ax.text(0.8, 0.6, f"HR {hr} (95%, {ci_lower}-{ci_upper})", transform=os_1L_kidney_ax.transAxes, fontsize=8, ha='center', va='center')
os_1L_kidney_ax.text(1.02, 0.53, f"p={p}", transform=os_1L_kidney_ax.transAxes, fontsize=8, ha='right', va='center')

# ###########################################################################
# # Panel C. mUC ctDNA PFS in 1L
# df = surv_df_bladder.merge(multivars_df_binarized_bladder, how = "left")
# df = df[~df["PFS (mo)"].isna()]
# pts_1L_bladder = clin_df_bladder[(clin_df_bladder["Event"] == "1L") & (clin_df_bladder["Patient_id"] != "23-088")]["Patient_id"].reset_index(drop = True)
# df = df[df["Patient_id"].isin(pts_1L_bladder)].reset_index(drop = True)
# del df["Patient_id"]
# del df["CH+"]

# pfs_bladder_ax = make_survival_curve(
#     df, 
#     "ctDNA+", 
#     event_col = "Progression",
#     duration_col = "PFS (mo)",
#     output_path = None, 
#     plot_title = "", 
#     print_HRs = False, 
#     add_legend = False,
#     pos_curve_label_positions = [0.2, 0.4], 
#     neg_curve_label_positions =  [0.45, 0.7], 
#     print_survival = False,
#     xlabel = "PFS (mo)",
#     ax = pfs_bladder_ax)

# # Calculate median survival to print on the plot
# muc_ctDNA_dict = return_median_survival(df, stratify_by="ctDNA+", event_col="Progression", duration_col="PFS (mo)")
# muc_cox_dict = run_cox_proportional_hazards(df, stratify_by="ctDNA+", event_col="Progression", duration_col="PFS (mo)")

# pfs_bladder_ax.set_title("mUC")
# pfs_bladder_ax.text(0.93, 1, "Median", transform=pfs_bladder_ax.transAxes, fontsize=8, ha='center', va='center')
# pfs_bladder_ax.text(0.93, 0.93, "mo", transform=pfs_bladder_ax.transAxes, fontsize=8, ha='center', va='center', fontstyle = "italic")
# pfs_bladder_ax.text(0.90, 0.85, 'ctDNA+', transform=pfs_bladder_ax.transAxes, fontsize=8, ha='right', va='center')
# pfs_bladder_ax.text(0.90, 0.78, 'ctDNA-', transform=pfs_bladder_ax.transAxes, fontsize=8, ha='right', va='center')
# pfs_bladder_ax.text(0.95, 0.85, muc_ctDNA_dict["Positive"], transform=pfs_bladder_ax.transAxes, fontsize=8, ha='center', va='center')
# pfs_bladder_ax.text(0.95, 0.78, muc_ctDNA_dict["Negative"], transform=pfs_bladder_ax.transAxes, fontsize=8, ha='center', va='center')

# # Cox results
# hr = round(muc_cox_dict["HR"], 2)
# ci_lower = round(muc_cox_dict["CI_lower"], 2)
# ci_upper = round(muc_cox_dict["CI_upper"], 2)
# p = round(float(muc_cox_dict['p']), 3)
# pfs_bladder_ax.text(0.8, 0.6, f"HR {hr} (95%, {ci_lower}-{ci_upper})", transform=pfs_bladder_ax.transAxes, fontsize=8, ha='center', va='center')
# pfs_bladder_ax.text(1.02, 0.53, f"p={p}", transform=pfs_bladder_ax.transAxes, fontsize=8, ha='right', va='center')

# ###########################################################################
# # Panel D. mRCC ctDNA PFS in 1L
# df = surv_df_kidney.merge(multivars_df_binarized_kidney, how = "left")
# pts_1L_kidney = clin_df_kidney[clin_df_kidney["GUBB draw minus date start"] < 21]["Patient_id"].reset_index(drop = True)
# df = df[df["Patient_id"].isin(pts_1L_kidney)].reset_index(drop = True)
# df = df[~df["PFS (mo)"].isna()]

# del df["Patient_id"]
# del df["CH+"]

# pfs_kidney_ax = make_survival_curve(
#     df, 
#     "ctDNA+",
#     event_col = "Progression",
#     duration_col = "PFS (mo)",
#     output_path = None,
#     plot_title = "", 
#     print_HRs = True,
#     pos_curve_label_positions = [0.27, 0.23], 
#     neg_curve_label_positions =  [0.5, 0.72], 
#     add_legend = False,
#     print_survival = False,
#     xlabel = "PFS (mo)",
#     ax = pfs_kidney_ax)

# pfs_kidney_ax.set_title("mRCC")

# mrcc_ctDNA_dict = return_median_survival(df, stratify_by="ctDNA+", event_col="Progression", duration_col="PFS (mo)")
# mrcc_cox_dict = run_cox_proportional_hazards(df, duration_col="PFS (mo)", event_col="Progression", stratify_by="ctDNA+", ax = None)

# # Log rank test
# pfs_kidney_ax.set_title("mRCC")
# pfs_kidney_ax.text(0.93, 1, "Median", transform=pfs_kidney_ax.transAxes, fontsize=8, ha='center', va='center')
# pfs_kidney_ax.text(0.93, 0.93, "mo", transform=pfs_kidney_ax.transAxes, fontsize=8, ha='center', va='center', fontstyle = "italic")
# pfs_kidney_ax.text(0.90, 0.85, 'ctDNA+', transform=pfs_kidney_ax.transAxes, fontsize=8, ha='right', va='center')
# pfs_kidney_ax.text(0.90, 0.78, 'ctDNA-', transform=pfs_kidney_ax.transAxes, fontsize=8, ha='right', va='center')
# pfs_kidney_ax.text(0.95, 0.85, mrcc_ctDNA_dict["Positive"], transform=pfs_kidney_ax.transAxes, fontsize=8, ha='center', va='center')
# pfs_kidney_ax.text(0.95, 0.78, mrcc_ctDNA_dict["Negative"], transform=pfs_kidney_ax.transAxes, fontsize=8, ha='center', va='center')

# # Cox results
# hr = round(mrcc_cox_dict["HR"], 2)
# ci_lower = round(mrcc_cox_dict["CI_lower"], 2)
# ci_upper = round(mrcc_cox_dict["CI_upper"], 2)
# p = round(float(mrcc_cox_dict['p']), 3)
# pfs_kidney_ax.text(0.8, 0.6, f"HR {hr} (95%, {ci_lower}-{ci_upper})", transform=pfs_kidney_ax.transAxes, fontsize=8, ha='center', va='center')
# pfs_kidney_ax.text(1.02, 0.53, f"p={p}", transform=pfs_kidney_ax.transAxes, fontsize=8, ha='right', va='center')

###########################################################################
# Panel E. mUC ctDNA and CH OS
genomics = ctDNA_status_bladder.merge(chip_status_bladder, how = "inner")
ctDNA_neg_bladder = genomics[genomics["ctDNA positive"] == False]["Patient_id"]
ctDNA_pos_CH_neg_bladder = genomics[(genomics["ctDNA positive"] == True) & (genomics["CHIP positive"] == False)]["Patient_id"]
ctDNA_pos_CH_pos_bladder = genomics[(genomics["ctDNA positive"] == True) & (genomics["CHIP positive"] == True)]["Patient_id"]

# Assign 0s and 1s to the columns based on the membership in each group
all_patients_bladder = pd.DataFrame({
    "Patient_id": pd.concat([ctDNA_neg_bladder, ctDNA_pos_CH_neg_bladder, ctDNA_pos_CH_pos_bladder]).unique()
})

all_patients_bladder["ctDNA-"] = all_patients_bladder["Patient_id"].isin(ctDNA_neg_bladder).astype(int)
all_patients_bladder["ctDNA+ CH-"] = all_patients_bladder["Patient_id"].isin(ctDNA_pos_CH_neg_bladder).astype(int)
all_patients_bladder["ctDNA+ CH+"] = all_patients_bladder["Patient_id"].isin(ctDNA_pos_CH_pos_bladder).astype(int)

df = surv_df_bladder.merge(all_patients_bladder)
del df["Patient_id"]

os_with_ch_bladder_ax.set_title("mUC")
os_with_ch_bladder_ax = make_survival_curve_with_forest_3_groups(
    df,
    os_with_ch_bladder_ax, 
    diagnosis = "Bladder",
    event_col = "Death",
    duration_col = "OS from cfDNA collection (mo)",
    pos_curve_label_positions = None, 
    neg_curve_label_positions = None, 
    add_legend = False)

# Run logrank test and print survival on plot
df_modified = df[df["ctDNA-"]==0] # Exclude patients that are ctDNA+ only without CH

muc_dict = return_median_survival(df_modified, stratify_by="ctDNA+ CH+", event_col="Death", duration_col="OS from cfDNA collection (mo)")
muc_cox = run_cox_proportional_hazards(df_modified, stratify_by="ctDNA+ CH+", event_col="Death", duration_col="OS from cfDNA collection (mo)")

hr = round(muc_cox["HR"], 2)
ci_upper = round(muc_cox["CI_upper"], 2)
ci_lower = round(muc_cox["CI_lower"], 2)
p = round(float(muc_cox["p"]), 5)

os_with_ch_bladder_ax.text(1, 0.95, "Median", transform=os_with_ch_bladder_ax.transAxes, fontsize=8, ha='center', va='center')
os_with_ch_bladder_ax.text(1, 0.88, "mo", transform=os_with_ch_bladder_ax.transAxes, fontsize=8, ha='center', va='center', fontstyle = "italic")
os_with_ch_bladder_ax.text(0.97, 0.81, 'ctDNA+ CH+', transform=os_with_ch_bladder_ax.transAxes, fontsize=8, ha='right', va='center')
os_with_ch_bladder_ax.text(0.97, 0.74, 'ctDNA+ CH-', transform=os_with_ch_bladder_ax.transAxes, fontsize=8, ha='right', va='center')
os_with_ch_bladder_ax.text(1, 0.81, muc_dict["Positive"], transform=os_with_ch_bladder_ax.transAxes, fontsize=8, ha='center', va='center')
os_with_ch_bladder_ax.text(1, 0.74, muc_dict["Negative"], transform=os_with_ch_bladder_ax.transAxes, fontsize=8, ha='center', va='center')
os_with_ch_bladder_ax.text(0.89, 0.67, f"HR: {hr} (95% CI: {ci_upper}-{ci_lower})", transform=os_with_ch_bladder_ax.transAxes, fontsize=8, ha='center', va='center')
os_with_ch_bladder_ax.text(0.89, 0.60, f"p={p}", transform=os_with_ch_bladder_ax.transAxes, fontsize=8, ha='center', va='center')

###########################################################################
# Panel F. mRCC ctDNA and CH OS
genomics = ctDNA_status_kidney.merge(chip_status_kidney, how = "inner")
ctDNA_neg_kidney = genomics[genomics["ctDNA positive"] == False]["Patient_id"]
ctDNA_pos_CH_neg_kidney = genomics[(genomics["ctDNA positive"] == True) & (genomics["CHIP positive"] == False)]["Patient_id"]
ctDNA_pos_CH_pos_kidney = genomics[(genomics["ctDNA positive"] == True) & (genomics["CHIP positive"] == True)]["Patient_id"]

# Assign 0s and 1s to the columns based on the membership in each group
all_patients_kidney = pd.DataFrame({
    "Patient_id": pd.concat([ctDNA_neg_kidney, ctDNA_pos_CH_neg_kidney, ctDNA_pos_CH_pos_kidney]).unique()
})

all_patients_kidney["ctDNA-"] = all_patients_kidney["Patient_id"].isin(ctDNA_neg_kidney).astype(int)
all_patients_kidney["ctDNA+ CH-"] = all_patients_kidney["Patient_id"].isin(ctDNA_pos_CH_neg_kidney).astype(int)
all_patients_kidney["ctDNA+ CH+"] = all_patients_kidney["Patient_id"].isin(ctDNA_pos_CH_pos_kidney).astype(int)

df = surv_df_kidney.merge(all_patients_kidney)
del df["Patient_id"]

os_with_ch_kidney_ax.set_title("mRCC")
os_with_ch_kidney_ax = make_survival_curve_with_forest_3_groups(
    df, 
    os_with_ch_kidney_ax, 
    diagnosis = "Kidney",
    event_col = "Death",
    duration_col = "OS from cfDNA collection (mo)",
    pos_curve_label_positions = None, 
    neg_curve_label_positions = None, 
    add_legend = False)

df_modified = df[df["ctDNA-"]==0]

mrcc_dict = return_median_survival(df_modified, stratify_by="ctDNA+ CH+", event_col="Death", duration_col="OS from cfDNA collection (mo)")
mrcc_cox = run_cox_proportional_hazards(df_modified, stratify_by="ctDNA+ CH+", event_col="Death", duration_col="OS from cfDNA collection (mo)")

hr = round(mrcc_cox["HR"], 2)
ci_upper = round(mrcc_cox["CI_upper"], 2)
ci_lower = round(mrcc_cox["CI_lower"], 2)
p = round(float(mrcc_cox["p"]), 5)

os_with_ch_kidney_ax.text(1, 0.95, "Median", transform=os_with_ch_kidney_ax.transAxes, fontsize=8, ha='center', va='center')
os_with_ch_kidney_ax.text(1, 0.88, "mo", transform=os_with_ch_kidney_ax.transAxes, fontsize=8, ha='center', va='center', fontstyle = "italic")
os_with_ch_kidney_ax.text(0.97, 0.81, 'ctDNA+ CH+', transform=os_with_ch_kidney_ax.transAxes, fontsize=8, ha='right', va='center')
os_with_ch_kidney_ax.text(0.97, 0.74, 'ctDNA+ CH-', transform=os_with_ch_kidney_ax.transAxes, fontsize=8, ha='right', va='center')
os_with_ch_kidney_ax.text(1, 0.81, mrcc_dict["Positive"], transform=os_with_ch_kidney_ax.transAxes, fontsize=8, ha='center', va='center')
os_with_ch_kidney_ax.text(1, 0.74, mrcc_dict["Negative"], transform=os_with_ch_kidney_ax.transAxes, fontsize=8, ha='center', va='center')
os_with_ch_kidney_ax.text(0.89, 0.67, f"HR: {hr} (95% CI: {ci_upper}-{ci_lower})", transform=os_with_ch_kidney_ax.transAxes, fontsize=8, ha='center', va='center')
os_with_ch_kidney_ax.text(0.89, 0.60, f"p={p}", transform=os_with_ch_kidney_ax.transAxes, fontsize=8, ha='center', va='center')

gs.tight_layout(fig)
fig.savefig(os.path.join(figure_dir, "SUPP_ctDNA_KM_1L.png"))
fig.savefig(os.path.join(figure_dir, "SUPP_ctDNA_KM_1L.pdf"))