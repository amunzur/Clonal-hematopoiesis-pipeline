
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
source_functions = os.path.join(DIR_working, "workflow/scripts/visualization/UTILITIES_make_chip_plots.py")
sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])

color_dict = {"Bladder": "deepskyblue", "Kidney": "orangered"}

with open(source_functions, 'r') as file:
    script_code = file.read()

exec(script_code)

with open("/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/visualization/survival_analysis/kidney_functions.py", 'r') as file:
    script_code = file.read()

exec(script_code)


PATH_clinical_bladder = os.path.join(DIR_working, "resources/clinical_data/bladder/clinical_data.csv")
PATH_kidney_clinical = '/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/RCC clinical data - mRCC clinical Data.csv'

kidney_pts = sample_info[sample_info["Diagnosis"] == "Kidney"]["Patient_id"].unique().tolist()
bladder_pts = sample_info[sample_info["Diagnosis"] == "Bladder"]["Patient_id"].unique().tolist()

# LOAD CHIP DATASETS
all_vars_chip = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/chip.csv"))
all_vars_chip = all_vars_chip[all_vars_chip["Dependent"] == False]
all_vars_chip = all_vars_chip[all_vars_chip["VAF_n"] > 2]
base_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
prog_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)
prog_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)

# LOAD SOMATIC DATASETS
all_vars_somatic = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/somatic.csv"))
all_vars_somatic = all_vars_somatic[all_vars_somatic["Dependent"] == False]
all_vars_somatic = all_vars_somatic[~all_vars_somatic["Patient_id"].isin(["20-313", "21-184", "21-430"])] # exclude some samples due to oxidative damage
base_kidney_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Kidney")].reset_index(drop = True)
prog_kidney_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "During treatment") & (all_vars_somatic["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Bladder")].reset_index(drop = True)
prog_bladder_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "During treatment") & (all_vars_somatic["Diagnosis"] == "Bladder")].reset_index(drop = True)

fig = plt.figure(figsize=(8, 4))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1], hspace = 0, wspace = 0) # outer most gs with 3 rows

kidney_ax = plt.subplot(gs[0])
bladder_ax = plt.subplot(gs[1])


#############################################################################
# FIGURE 3A . OS. CHIP+ vs CHIP- in kidney.
clin_df = pd.read_csv(PATH_kidney_clinical)

sex_and_age_df, cci_df, mets_df, imdc_df, subtype_df, surv_df = prepare_clinical_data(clin_df)
surv_df = surv_df[surv_df["Patient_id"].isin(kidney_pts)]
ctDNA_status = prepare_ctDNA_status(base_kidney_somatic, "Kidney", PATH_sample_information)
chip_status = prepare_chip_status(base_kidney_chip, "Kidney", PATH_sample_information)
multivars_df_binarized = prepare_multivars_df(sex_and_age_df, cci_df, mets_df, ctDNA_status, subtype_df, chip_status, imdc_df)
merged_df_kidney = surv_df.merge(multivars_df_binarized)
del merged_df_kidney["Patient_id"]

kidney_ax = make_survival_curve(
    merged_df_kidney, 
    "CH+",
    event_col = "Death",
    duration_col = "OS from cfDNA collection (mo)",
    output_path = None,
    plot_title = "", 
    print_HRs = True,
    add_legend = False, 
    print_survival = None, 
    xlabel = "OS from cfDNA collection (mo)",
    ax = kidney_ax)

kidney_ax.set_title("mRCC (CH>2%)")

#############################################################################
# FIGURE 3B . OS. CH+ vs CH- in bladder.
clin_df = pd.read_csv(PATH_clinical_bladder)
clin_df = clin_df[clin_df["First sample?"] == True]

sex_and_age_df_bladder = clin_df[["Patient_id", "Sex", "Age at blood draw"]]
smoking_status = clin_df[["Patient_id", "Previous smoking history"]]
ctDNA_status = prepare_ctDNA_status(base_bladder_somatic, "Bladder", PATH_sample_information)
chip_status = prepare_chip_status(base_bladder_chip, "Bladder", PATH_sample_information)
surv_df = prepare_survival_data(clin_df)
surv_df = surv_df[surv_df["Patient_id"].isin(bladder_pts)]
multivars_df_binarized = prepare_multivars_df_bladder(sex_and_age_df_bladder, smoking_status, ctDNA_status, chip_status)
merged_df_bladder = surv_df.merge(multivars_df_binarized, how = "left")
merged_df_bladder["Death"] = merged_df_bladder["Death"].map({True: True, False: False, 'True': True, 'False': False, 'Lost to follow-up': False}) # Map values to boolean, treating 'Lost to follow-up' as False
del merged_df_bladder["Patient_id"]

bladder_ax = make_survival_curve(
    merged_df_bladder, 
    "CH+", 
    event_col = "Death",
    duration_col = "OS from cfDNA collection (mo)",
    output_path = None, 
    plot_title = "", 
    print_HRs = True, 
    add_legend = False, 
    print_survival = None, 
    xlabel = "OS from cfDNA collection (mo)",
    ax = bladder_ax)

bladder_ax.set_title("mUC (CH>2%)")

# Add curve annotations
kidney_ax.text(0.4, 0.7, "CHIP+", transform=kidney_ax.transAxes, fontsize=10, ha='left', va='center', color = "red")
kidney_ax.text(0.2, 0.43, "CHIP-", transform=kidney_ax.transAxes, fontsize=10, ha='left', va='center', color = "black")

bladder_ax.text(0.15, 0.3, "CHIP+", transform=bladder_ax.transAxes, fontsize=10, ha='left', va='center', color = "red")
bladder_ax.text(0.35, 0.5, "CHIP-", transform=bladder_ax.transAxes, fontsize=10, ha='left', va='center', color = "black")

# Add median surv
# Calculate median surv values
kidney_dict = return_median_survival(merged_df_kidney, "CH+", "Death", "OS from cfDNA collection (mo)")
bladder_dict = return_median_survival(merged_df_bladder, "CH+", "Death", "OS from cfDNA collection (mo)")

# annotate the median surv values
for ax, dic, keyword in zip([kidney_ax, bladder_ax], [kidney_dict, bladder_dict], ["CHIP", "CHIP"]):
    ax.text(0.9, 0.9, "Median", transform=ax.transAxes, fontsize=8, ha='center', va='center')
    ax.text(0.9, 0.83, "months", transform=ax.transAxes, fontsize=8, ha='center', va='center', fontstyle = "italic")
    ax.text(0.7, 0.76, f"{keyword}+", transform=ax.transAxes, fontsize=8, ha='left', va='center', color = "red")
    ax.text(0.7, 0.69, f"{keyword}-", transform=ax.transAxes, fontsize=8, ha='left', va='center', color = "black")
    pos_median = dic["Positive"]
    neg_median = dic["Negative"]
    ax.text(0.9, 0.76, pos_median, transform=ax.transAxes, fontsize=8, ha='center', va='center')
    ax.text(0.9, 0.69, neg_median, transform=ax.transAxes, fontsize=8, ha='center', va='center')
    if "logrank p" in dic.keys(): 
        p = dic["logrank p"]
        ax.text(0.65, 0.62, f"Logrank p={p}", transform=ax.transAxes, fontsize=8, ha='left', va='center')


fig.text(0.07, 0.98, 'A', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
fig.text(0.55, 0.98, 'B', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)

gs.tight_layout(fig)
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/SUPP_OS_CHIP.pdf")
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/SUPP_OS_CHIP.png")


