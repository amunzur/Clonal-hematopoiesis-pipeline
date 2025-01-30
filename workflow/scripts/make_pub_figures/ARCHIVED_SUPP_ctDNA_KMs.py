
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
PATH_sample_information = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"
PATH_kidney_clinical = '/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/RCC clinical data - mRCC clinical Data.csv'
PATH_clinical_bladder = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/bladder/clinical_data.csv"

figure_dir = os.path.join(DIR_working, "results/figures/pub_figures")
source_functions = os.path.join(DIR_working, "workflow/scripts/visualization/UTILITIES_make_chip_plots.py")
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
prog_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)
prog_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)

# LOAD SOMATIC DATASETS
all_vars_somatic = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/somatic.csv"))
all_vars_somatic = all_vars_somatic[all_vars_somatic["Dependent"] == False]
# all_vars_somatic = all_vars_somatic[all_vars_somatic["VAF_t"] > 1]

all_vars_somatic = all_vars_somatic[~all_vars_somatic["Patient_id"].isin(["20-313", "21-184", "21-430"])] # exclude some samples due to oxidative damage
base_kidney_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Kidney")].reset_index(drop = True)
prog_kidney_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "During treatment") & (all_vars_somatic["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Bladder")].reset_index(drop = True)
prog_bladder_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "During treatment") & (all_vars_somatic["Diagnosis"] == "Bladder")].reset_index(drop = True)

all_pts = pd.read_csv(PATH_sample_information, sep="\t", names=["Patient_id", "Date", "Diagnosis", "Timepoint"])
kidney_pts = all_pts[all_pts["Diagnosis"] == "Kidney"]["Patient_id"].unique()
bladder_pts = all_pts[all_pts["Diagnosis"] == "Bladder"]["Patient_id"].unique()

fig = plt.figure(figsize=(8, 11))
outer_gs = gridspec.GridSpec(4, 2, height_ratios=[1, 1, 1, 1], width_ratios = [1, 0.4], hspace=0.01, wspace=0.01)  # Adjusted hspace
inner_gs0_km = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=[1, 0.3], subplot_spec=outer_gs[0, 0], wspace=0.01, hspace = 0.2) # OS - km and text for kidney
inner_gs1_km = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=[1, 0.3], subplot_spec=outer_gs[1, 0], wspace=0.01, hspace = 0.2) # OS - km and text for bladder
inner_gs2_km = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=[1, 0.3], subplot_spec=outer_gs[2, 0], wspace=0.01, hspace = 0.2) # PFS - km and text for kidney
inner_gs3_km = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=[1, 0.3], subplot_spec=outer_gs[3, 0], wspace=0.01, hspace = 0.2) # PFS - km and text for bladder


# inner_gs1 = gridspec.GridSpecFromSubplotSpec(1, 4, width_ratios=[1, 0.5, 1, 0.5], subplot_spec=outer_gs[1], wspace=0.35, hspace = 0.3)
# inner_gs2 = gridspec.GridSpecFromSubplotSpec(1, 4, width_ratios=[1, 0.5, 1, 0.5], subplot_spec=outer_gs[2], wspace=0.35, hspace = 0.3)



with open(source_functions, 'r') as file:
    script_code = file.read()

exec(script_code)


#############################################################################
# FIGURE 1A and B. OS. ctDNA+ vs ctDNA- in kidney.
clin_df = pd.read_csv(PATH_kidney_clinical)

ax0_km = plt.subplot(inner_gs0_km[0])
ax1_survival_printed = plt.subplot(inner_gs0_km[1])
ax2_forest = plt.subplot(outer_gs[0, 1]) # OS kidney forest

sex_and_age_df, cci_df, mets_df, imdc_df, subtype_df, surv_df = prepare_clinical_data(clin_df)
surv_df = surv_df[surv_df["Patient_id"].isin(kidney_pts)]
ctDNA_status = prepare_ctDNA_status(base_kidney_somatic, "Kidney", PATH_sample_information)
chip_status = prepare_chip_status(base_kidney_chip, "Kidney", PATH_sample_information)
multivars_df_binarized = prepare_multivars_df(sex_and_age_df, cci_df, mets_df, ctDNA_status, subtype_df, chip_status, imdc_df)
merged_df_kidney = surv_df.merge(multivars_df_binarized)
del merged_df_kidney["Patient_id"]

ax0_km = make_survival_curve(
    merged_df_kidney, 
    "ctDNA+",
    event_col = "Death",
    duration_col = "OS from cfDNA collection (mo)",
    output_path = None,
    plot_title = "", 
    print_HRs = True,
    pos_curve_label_positions = [0.3, 0.25], 
    neg_curve_label_positions =  [0.75, 0.75], 
    add_legend = False, 
    print_survival = ax1_survival_printed, 
    xlabel = "OS from cfDNA collection (mo)",
    ax = ax0_km)

ax0_km.set_title("mRCC")
ax0_km.set_title("mRCC")

ax1_survival_printed.spines[["right", "top", "bottom", "left"]].set_visible(False)
ax1_survival_printed.set_yticks([])
ax1_survival_printed.set_xticks([])

ax2_forest, cph = plot_cph_forest(merged_df_kidney.drop(["Progression", "PFS (mo)"], axis = 1), ax2_forest, "Death", "OS from cfDNA collection (mo)")
ax2_forest.yaxis.set_tick_params(labelsize=9)


#############################################################################
# FIGURE 1C and D. OS. ctDNA+ vs ctDNA- in bladder.
clin_df = pd.read_csv(PATH_clinical_bladder)
clin_df = clin_df[clin_df["First sample?"] == True]

ax3_km = plt.subplot(inner_gs1_km[0])
ax4_survival_printed = plt.subplot(inner_gs1_km[1])
ax5_forest = plt.subplot(outer_gs[1, 1]) # OS bladder forest

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

ax3_km = make_survival_curve(
    merged_df_bladder, 
    "ctDNA+", 
    event_col = "Death",
    duration_col = "OS from cfDNA collection (mo)",
    output_path = None, 
    plot_title = "", 
    print_HRs = True, 
    pos_curve_label_positions = [0.15, 0.3], 
    neg_curve_label_positions =  [0.5, 0.7], 
    add_legend = False, 
    print_survival = ax4_survival_printed, 
    xlabel = "OS from cfDNA collection (mo)",
    ax = ax3_km)

ax3_km.set_title("mUC")

ax4_survival_printed.spines[["right", "top", "bottom", "left"]].set_visible(False)
ax4_survival_printed.set_yticks([])
ax4_survival_printed.set_xticks([])

ax5_forest, cph = plot_cph_forest(merged_df_bladder.drop(["Progression", "PFS (mo)"], axis = 1), ax5_forest, "Death", "OS from cfDNA collection (mo)")
ax5_forest.yaxis.set_tick_params(labelsize=9)

#############################################################################
# FIGURE 1A and B. PFS. ctDNA+ vs ctDNA- in kidney.
ax6_km = plt.subplot(inner_gs2_km[0])
ax7_survival_printed = plt.subplot(inner_gs2_km[1])
ax8_forest = plt.subplot(outer_gs[2, 1]) # PFS kidney forest

merged_df_kidney_pfs = merged_df_kidney[~merged_df_kidney["PFS (mo)"].isna()]
merged_df_kidney_pfs["PFS (mo)"] = merged_df_kidney_pfs["PFS (mo)"].astype(float)
merged_df_kidney_pfs["Progression"] = merged_df_kidney_pfs["Progression"].astype(bool)

ax6_km = make_survival_curve(
    merged_df_kidney_pfs, 
    "ctDNA+",
    event_col = "Progression",
    duration_col = "PFS (mo)",
    output_path = None,
    plot_title = "", 
    print_HRs = True,
    pos_curve_label_positions = [0.65, 0.1], 
    neg_curve_label_positions =  [0.4, 0.35], 
    add_legend = False, 
    print_survival = ax7_survival_printed, 
    xlabel = "PFS (mo)",
    ax = ax6_km)

ax6_km.set_title("mRCC")

ax7_survival_printed.spines[["right", "top", "bottom", "left"]].set_visible(False)
ax7_survival_printed.set_yticks([])
ax7_survival_printed.set_xticks([])

ax8_forest, cph = plot_cph_forest(merged_df_kidney_pfs.drop(["Death", "OS from cfDNA collection (mo)"], axis = 1), ax8_forest, "Progression", "PFS (mo)")
ax8_forest.yaxis.set_tick_params(labelsize=9)


#############################################################################
# FIGURE 1A and B. PFS. ctDNA+ vs ctDNA- in bladder.
ax9_km = plt.subplot(inner_gs3_km[0])
ax10_survival_printed = plt.subplot(inner_gs3_km[1])
ax11_forest = plt.subplot(outer_gs[3, 1]) # PFS kidney forest

merged_df_bladder_pfs = merged_df_bladder[~merged_df_bladder["PFS (mo)"].isna()]
merged_df_bladder_pfs["PFS (mo)"] = merged_df_bladder_pfs["PFS (mo)"].astype(float)
merged_df_bladder_pfs["Progression"] = merged_df_bladder_pfs["Progression"].astype(bool)

ax9_km = make_survival_curve(
    merged_df_bladder_pfs, 
    "ctDNA+", 
    event_col = "Progression",
    duration_col = "PFS (mo)",
    output_path = None, 
    plot_title = "", 
    print_HRs = True, 
    pos_curve_label_positions = [0.15, 0.1], 
    neg_curve_label_positions =  [0.5, 0.65], 
    add_legend = False, 
    print_survival = ax10_survival_printed, 
    xlabel = "PFS (mo)",
    ax = ax9_km)

ax11_forest, cph = plot_cph_forest(merged_df_bladder_pfs.drop(["Death", "OS from cfDNA collection (mo)"], axis = 1), ax11_forest, "Progression", "PFS (mo)")

ax11_forest.yaxis.set_tick_params(labelsize=9)
ax9_km.set_title("mUC")
ax10_survival_printed.spines[["right", "top", "bottom", "left"]].set_visible(False)
ax10_survival_printed.set_yticks([])
ax10_survival_printed.set_xticks([])

fig.text(0.04, 0.98, 'A', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
fig.text(0.04, 0.74, 'B', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
fig.text(0.04, 0.50, 'C', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
fig.text(0.04, 0.25, 'D', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)

outer_gs.tight_layout(fig)
fig.savefig(os.path.join(figure_dir, "SUPP_ctDNA_KMs.png"))
fig.savefig(os.path.join(figure_dir, "SUPP_ctDNA_KMs.pdf"))
