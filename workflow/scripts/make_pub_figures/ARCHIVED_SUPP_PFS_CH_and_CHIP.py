
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
PATH_kidney_clinical = '/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/RCC clinical data - mRCC clinical Data.csv'
PATH_clinical_bladder = os.path.join(DIR_working, "resources/clinical_data/bladder/clinical_data.csv")

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


PATH_clinical_bladder = os.path.join(DIR_working, "resources/clinical_data/Bladder_enrollment.csv")

# LOAD CHIP DATASETS
all_vars_chip = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/chip.csv"))
all_vars_chip = all_vars_chip[all_vars_chip["Dependent"] == False]
base_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
prog_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)
prog_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)

base_kidney_chip_2_perc = base_kidney_chip[base_kidney_chip["VAF_n"] > 2]
base_bladder_chip_2_perc = base_bladder_chip[base_bladder_chip["VAF_n"] > 2]
# LOAD SOMATIC DATASETS
all_vars_somatic = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/somatic.csv"))
all_vars_somatic = all_vars_somatic[(all_vars_somatic["Dependent"] == False) & (all_vars_chip["Timepoint"] == "Baseline")]
all_vars_somatic = all_vars_somatic[~all_vars_somatic["Patient_id"].isin(["20-313", "21-184", "21-430"])] # exclude some samples due to oxidative damage
base_kidney_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Kidney")].reset_index(drop = True)
prog_kidney_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "During treatment") & (all_vars_somatic["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Bladder")].reset_index(drop = True)
prog_bladder_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "During treatment") & (all_vars_somatic["Diagnosis"] == "Bladder")].reset_index(drop = True)


#############################################################################
# SET UP KIDNEY DFS
clin_df = pd.read_csv(PATH_kidney_clinical)

ax0_km = plt.subplot(inner_gs0[0])
ax1_survival_printed = plt.subplot(inner_gs0[1])

sex_and_age_df, cci_df, mets_df, imdc_df, subtype_df, surv_df = prepare_clinical_data(clin_df)
surv_df = surv_df[surv_df["Patient_id"].isin(kidney_pts)]
ctDNA_status = prepare_ctDNA_status(base_kidney_somatic, "Kidney", PATH_sample_information)

# CH 
chip_status = prepare_chip_status(base_kidney_chip, "Kidney", PATH_sample_information)
multivars_df_binarized = prepare_multivars_df(sex_and_age_df, cci_df, mets_df, ctDNA_status, subtype_df, chip_status, imdc_df)
merged_df_kidney_ch = surv_df.merge(multivars_df_binarized)
del merged_df_kidney_ch["Patient_id"]
merged_df_kidney_ch_pfs = merged_df_kidney_ch[~merged_df_kidney_ch["PFS (mo)"].isna()]
merged_df_kidney_ch_pfs["PFS (mo)"] = merged_df_kidney_ch_pfs["PFS (mo)"].astype(float)
merged_df_kidney_ch_pfs["Progression"] = merged_df_kidney_ch_pfs["Progression"].astype(bool)

# CHIP
chip_status = prepare_chip_status(base_kidney_chip_2_perc, "Kidney", PATH_sample_information)
multivars_df_binarized = prepare_multivars_df(sex_and_age_df, cci_df, mets_df, ctDNA_status, subtype_df, chip_status, imdc_df)
merged_df_kidney_chip = surv_df.merge(multivars_df_binarized)
del merged_df_kidney_chip["Patient_id"]
merged_df_kidney_chip_pfs = merged_df_kidney_chip[~merged_df_kidney_chip["PFS (mo)"].isna()]
merged_df_kidney_chip_pfs["PFS (mo)"] = merged_df_kidney_ch_pfs["PFS (mo)"].astype(float)
merged_df_kidney_chip_pfs["Progression"] = merged_df_kidney_ch_pfs["Progression"].astype(bool)


#############################################################################
# SET UP BLADDER DFS
clin_df = pd.read_csv(PATH_clinical_bladder)
clin_df = clin_df[(clin_df["First sample?"] == True) & (clin_df["Event"] != "Palliative / best supportive care")]

sex_and_age_df_bladder = clin_df[["Patient_id", "Sex", "Age at blood draw"]]
smoking_status = clin_df[["Patient_id", "Previous smoking history"]]
ctDNA_status = prepare_ctDNA_status(base_bladder_somatic, "Bladder", PATH_sample_information)
surv_df = prepare_survival_data(clin_df)
surv_df = surv_df[surv_df["Patient_id"].isin(bladder_pts)]

# CH
chip_status = prepare_chip_status(base_bladder_chip, "Bladder", PATH_sample_information)
multivars_df_binarized = prepare_multivars_df_bladder(sex_and_age_df_bladder, smoking_status, ctDNA_status, chip_status)
merged_df_bladder_ch = surv_df.merge(multivars_df_binarized, how = "left")
merged_df_bladder_ch["Death"] = merged_df_bladder_ch["Death"].map({True: True, False: False, 'True': True, 'False': False, 'Lost to follow-up': False}) # Map values to boolean, treating 'Lost to follow-up' as False
del merged_df_bladder_ch["Patient_id"]
merged_df_bladder_ch_pfs = merged_df_bladder_ch[~merged_df_bladder_ch["PFS (mo)"].isna()]
merged_df_bladder_ch_pfs["PFS (mo)"] = merged_df_bladder_ch_pfs["PFS (mo)"].astype(float)
merged_df_bladder_ch_pfs["Progression"] = merged_df_bladder_ch_pfs["Progression"].astype(bool)

# chip
chip_status = prepare_chip_status(base_bladder_chip_2_perc, "Bladder", PATH_sample_information)
multivars_df_binarized = prepare_multivars_df_bladder(sex_and_age_df_bladder, smoking_status, ctDNA_status, chip_status)
merged_df_bladder_chip = surv_df.merge(multivars_df_binarized, how = "left")
merged_df_bladder_chip["Death"] = merged_df_bladder_chip["Death"].map({True: True, False: False, 'True': True, 'False': False, 'Lost to follow-up': False}) # Map values to boolean, treating 'Lost to follow-up' as False
del merged_df_bladder_chip["Patient_id"]
merged_df_bladder_chip_pfs = merged_df_bladder_chip[~merged_df_bladder_ch["PFS (mo)"].isna()]
merged_df_bladder_chip_pfs["PFS (mo)"] = merged_df_bladder_chip_pfs["PFS (mo)"].astype(float)
merged_df_bladder_chip_pfs["Progression"] = merged_df_bladder_chip_pfs["Progression"].astype(bool)

fig = plt.figure(figsize=(8, 7))
outer_gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[1, 1], hspace = 0, wspace = 0) # outer most gs with 3 rows

kidney_ch = plt.subplot(outer_gs[0])
bladder_ch = plt.subplot(outer_gs[1])
kidney_chip = plt.subplot(outer_gs[2])
bladder_chip = plt.subplot(outer_gs[3])

kidney_ch = make_survival_curve(
    merged_df_kidney_ch_pfs, 
    "CH+",
    event_col = "Progression",
    duration_col = "PFS (mo)",
    output_path = None,
    plot_title = "", 
    print_HRs = True,
    xmax = 50,
    add_legend = False, 
    print_survival = None, 
    xlabel = "PFS (mo)",
    ax = kidney_ch)

bladder_ch = make_survival_curve(
    merged_df_bladder_ch_pfs, 
    "CH+",
    event_col = "Progression",
    duration_col = "PFS (mo)",
    output_path = None,
    plot_title = "", 
    print_HRs = True,
    xmax = 15,
    add_legend = False, 
    print_survival = None, 
    xlabel = "PFS (mo)",
    ax = bladder_ch)

kidney_chip = make_survival_curve(
    merged_df_kidney_chip_pfs, 
    "CH+",
    event_col = "Progression",
    duration_col = "PFS (mo)",
    output_path = None,
    plot_title = "", 
    print_HRs = True,
    xmax = 50,
    add_legend = False, 
    print_survival = None, 
    xlabel = "PFS (mo)",
    ax = kidney_chip)

bladder_chip = make_survival_curve(
    merged_df_bladder_chip_pfs, 
    "CH+",
    event_col = "Progression",
    duration_col = "PFS (mo)",
    output_path = None,
    plot_title = "", 
    print_HRs = True,
    xmax = 15,
    add_legend = False, 
    print_survival = None, 
    xlabel = "PFS (mo)",
    ax = bladder_chip)

# AESTHETICS
kidney_ch.set_title("mRCC (all CH variants)")
bladder_ch.set_title("mUC (all CH variants)")
kidney_chip.set_title("mRCC (CH>2%)")
bladder_chip.set_title("mUC (CH>2%)")

# Add curve annotations
kidney_ch.text(0.43, 0.4, "CH+", transform=kidney_ch.transAxes, fontsize=10, ha='left', va='center', color = "red")
kidney_ch.text(0.22, 0.18, "CH-", transform=kidney_ch.transAxes, fontsize=10, ha='left', va='center', color = "black")

bladder_ch.text(0.26, 0.36, "CH+", transform=bladder_ch.transAxes, fontsize=10, ha='left', va='center', color = "red")
bladder_ch.text(0.5, 0.6, "CH-", transform=bladder_ch.transAxes, fontsize=10, ha='left', va='center', color = "black")

kidney_chip.text(0.45, 0.45, "CHIP+", transform=kidney_chip.transAxes, fontsize=10, ha='left', va='center', color = "red")
kidney_chip.text(0.25, 0.15, "CHIP-", transform=kidney_chip.transAxes, fontsize=10, ha='left', va='center', color = "black")

bladder_chip.text(0.3, 0.7, "CHIP+", transform=bladder_chip.transAxes, fontsize=10, ha='left', va='center', color = "red")
bladder_chip.text(0.25, 0.3, "CHIP-", transform=bladder_chip.transAxes, fontsize=10, ha='left', va='center', color = "black")

# Calculate median surv values
kidney_ch_dict = return_median_survival(merged_df_kidney_ch_pfs, "CH+", "Progression", "PFS (mo)")
bladder_ch_dict = return_median_survival(merged_df_bladder_ch_pfs, "CH+", "Progression", "PFS (mo)")
kidney_chip_dict = return_median_survival(merged_df_kidney_chip_pfs, "CH+", "Progression", "PFS (mo)")
bladder_chip_dict = return_median_survival(merged_df_bladder_chip_pfs, "CH+", "Progression", "PFS (mo)")

# annotate the median surv values

for ax, dic, keyword in zip([kidney_ch, bladder_ch, kidney_chip, bladder_chip], [kidney_ch_dict, bladder_ch_dict, kidney_chip_dict, bladder_chip_dict], ["CH", "CH", "CHIP", "CHIP"]):
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

fig.text(0.05, 0.98, 'A', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
fig.text(0.55, 0.98, 'B', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
fig.text(0.05, 0.47, 'C', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
fig.text(0.55, 0.47, 'D', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)

outer_gs.tight_layout(fig)
fig.savefig(os.path.join(figure_dir, "SUPP_CH_PFS.png"))
fig.savefig(os.path.join(figure_dir, "SUPP_CH_PFS.pdf"))





# sex_and_age_df, cci_df, mets_df, imdc_df, subtype_df = prepare_clinical_data(clin_df)
# ctDNA_status = prepare_ctDNA_status(base_kidney_somatic, diagnosis, PATH_sample_information)
# chip_status = prepare_chip_status(base_kidney_chip, diagnosis, PATH_sample_information)

# genomics_df = chip_status.merge(ctDNA_status)
# genomics_df["ctDNA-"] = (genomics_df["ctDNA negative"] == True)
# genomics_df["ctDNA+"] = (genomics_df["CHIP positive"] == False) & (genomics_df["ctDNA positive"] == True)
# genomics_df["CH and ctDNA"] = (genomics_df["CHIP positive"] == True) & (genomics_df["ctDNA positive"] == True)
# genomics_df = genomics_df[["Patient_id", "CH only", "ctDNA only", "CH and ctDNA"]]

# multivars_df_binarized = prepare_multivars_df_with_3_genomic_vars(sex_and_age_df, cci_df, mets_df, subtype_df, imdc_df, genomics_df)
# surv_df = make_survival_df(PATH_kidney_clean, PATH_sample_information)
# merged_df = surv_df.merge(multivars_df_binarized)
# merged_df = merged_df[['Overall survival', 'Death at last follow up',
#                        'Age at draw', 'CCI>median', 'Male', 'Visceral mets',
#                        "CH+", "ctDNA+", "CH+ and ctDNA+",
#                        'Clear cell', 'IMDC poor/int. risk']]   

# ax0_km = make_survival_curve_with_forest_3_groups(merged_df, ax0_km, ax1_forest, add_legend = False)
# ax1_forest, cph = plot_cph_forest(merged_df, ax1_forest)


# outer_gs.tight_layout(fig)
# fig.savefig(os.path.join(figure_dir, "figure3.png"))
# fig.savefig(os.path.join(figure_dir, "figure3.pdf"))
