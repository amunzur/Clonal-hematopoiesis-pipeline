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
PATH_kidney_clinical = os.path.join(DIR_working, "resources/clinical_data/RCC clinical data - mRCC clinical Data.csv")
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
fig = plt.figure(figsize=(8, 4))
gs = gridspec.GridSpec(1, 2, width_ratios = [1, 1], hspace = 0, wspace = 0)

bladder_ax = plt.subplot(gs[0])
kidney_ax = plt.subplot(gs[1])

############################ BLADDER
clin_df = pd.read_csv(PATH_clinical_bladder)
clin_df = clin_df[(clin_df["First sample?"] == True)]

sex_and_age_df_bladder = clin_df[["Patient_id", "Sex", "Age at blood draw"]]
smoking_status = clin_df[["Patient_id", "Previous smoking history"]]
ctDNA_status = prepare_ctDNA_status(base_bladder_somatic, "Bladder", PATH_sample_information)
chip_status = prepare_chip_status(base_bladder_chip, "Bladder", PATH_sample_information)

genomics = ctDNA_status.merge(chip_status, how = "inner")
ctDNA_neg = genomics[genomics["ctDNA positive"] == False]["Patient_id"]
ctDNA_pos_CH_neg = genomics[(genomics["ctDNA positive"] == True) & (genomics["CHIP positive"] == False)]["Patient_id"]
ctDNA_pos_CH_pos = genomics[(genomics["ctDNA positive"] == True) & (genomics["CHIP positive"] == True)]["Patient_id"]

# Assign 0s and 1s to the columns based on the membership in each group
all_patients = pd.DataFrame({
    "Patient_id": pd.concat([ctDNA_neg, ctDNA_pos_CH_neg, ctDNA_pos_CH_pos]).unique()
})

all_patients["ctDNA-"] = all_patients["Patient_id"].isin(ctDNA_neg).astype(int)
all_patients["ctDNA+ CH-"] = all_patients["Patient_id"].isin(ctDNA_pos_CH_neg).astype(int)
all_patients["ctDNA+ CH+"] = all_patients["Patient_id"].isin(ctDNA_pos_CH_pos).astype(int)

surv_df = prepare_survival_data(clin_df)
surv_df = surv_df[surv_df["Patient_id"].isin(bladder_pts)]
merged_df_bladder = surv_df.merge(all_patients)
merged_df_bladder["Death"] = merged_df_bladder["Death"].map({True: True, False: False, 'True': True, 'False': False, 'Lost to follow-up': False}) # Map values to boolean, treating 'Lost to follow-up' as False
del merged_df_bladder["Patient_id"]

bladder_ax.set_title("mUC")
bladder_ax = make_survival_curve_with_forest_3_groups(
    merged_df_bladder,
    bladder_ax, 
    diagnosis = "Bladder",
    event_col = "Death",
    duration_col = "OS from cfDNA collection (mo)",
    pos_curve_label_positions = None, 
    neg_curve_label_positions = None, 
    add_legend = False)

# Run logrank test and print survival on plot
merged_df_bladder_modified = merged_df_bladder[merged_df_bladder["ctDNA-"]==0]

muc_dict = return_median_survival(merged_df_bladder_modified, stratify_by="ctDNA+ CH+", event_col="Death", duration_col="OS from cfDNA collection (mo)")
muc_cox = run_cox_proportional_hazards(merged_df_bladder_modified, stratify_by="ctDNA+ CH+", event_col="Death", duration_col="OS from cfDNA collection (mo)")

hr = round(muc_cox["HR"], 2)
ci_upper = round(muc_cox["CI_upper"], 2)
ci_lower = round(muc_cox["CI_lower"], 2)
p = round(float(muc_cox["p"]), 5)

bladder_ax.text(1, 0.95, "Median", transform=bladder_ax.transAxes, fontsize=8, ha='center', va='center')
bladder_ax.text(1, 0.88, "mo", transform=bladder_ax.transAxes, fontsize=8, ha='center', va='center', fontstyle = "italic")
bladder_ax.text(0.97, 0.81, 'ctDNA+ CH+', transform=bladder_ax.transAxes, fontsize=8, ha='right', va='center')
bladder_ax.text(0.97, 0.74, 'ctDNA+ CH-', transform=bladder_ax.transAxes, fontsize=8, ha='right', va='center')
bladder_ax.text(1, 0.81, muc_dict["Positive"], transform=bladder_ax.transAxes, fontsize=8, ha='center', va='center')
bladder_ax.text(1, 0.74, muc_dict["Negative"], transform=bladder_ax.transAxes, fontsize=8, ha='center', va='center')
bladder_ax.text(0.89, 0.67, f"HR: {hr} (95% CI: {ci_upper}-{ci_lower})", transform=bladder_ax.transAxes, fontsize=8, ha='center', va='center')
bladder_ax.text(0.89, 0.60, f"p={p}", transform=bladder_ax.transAxes, fontsize=8, ha='center', va='center')

############################ KIDNEY
clin_df = pd.read_csv(PATH_kidney_clinical)
sex_and_age_df, cci_df, mets_df, imdc_df, subtype_df, surv_df = prepare_clinical_data(clin_df)
surv_df = surv_df[surv_df["Patient_id"].isin(kidney_pts)]
ctDNA_status = prepare_ctDNA_status(base_kidney_somatic, "Kidney", PATH_sample_information)
chip_status = prepare_chip_status(base_kidney_chip, "Kidney", PATH_sample_information)

genomics = ctDNA_status.merge(chip_status, how = "inner")
ctDNA_neg = genomics[genomics["ctDNA positive"] == False]["Patient_id"]
ctDNA_pos_CH_neg = genomics[(genomics["ctDNA positive"] == True) & (genomics["CHIP positive"] == False)]["Patient_id"]
ctDNA_pos_CH_pos = genomics[(genomics["ctDNA positive"] == True) & (genomics["CHIP positive"] == True)]["Patient_id"]

# Assign 0s and 1s to the columns based on the membership in each group
all_patients = pd.DataFrame({
    "Patient_id": pd.concat([ctDNA_neg, ctDNA_pos_CH_neg, ctDNA_pos_CH_pos]).unique()
})

all_patients["ctDNA-"] = all_patients["Patient_id"].isin(ctDNA_neg).astype(int)
all_patients["ctDNA+ CH-"] = all_patients["Patient_id"].isin(ctDNA_pos_CH_neg).astype(int)
all_patients["ctDNA+ CH+"] = all_patients["Patient_id"].isin(ctDNA_pos_CH_pos).astype(int)

merged_df_kidney = surv_df.merge(all_patients)
del merged_df_kidney["Patient_id"]

kidney_ax.set_title("mRCC")
kidney_ax = make_survival_curve_with_forest_3_groups(
    merged_df_kidney, 
    kidney_ax, 
    diagnosis = "Kidney",
    event_col = "Death",
    duration_col = "OS from cfDNA collection (mo)",
    pos_curve_label_positions = None, 
    neg_curve_label_positions = None, 
    add_legend = False)

merged_df_kidney_modified = merged_df_kidney[merged_df_kidney["ctDNA-"]==0]

mrcc_dict = return_median_survival(merged_df_kidney_modified, stratify_by="ctDNA+ CH+", event_col="Death", duration_col="OS from cfDNA collection (mo)")
mrcc_cox = run_cox_proportional_hazards(merged_df_kidney_modified, stratify_by="ctDNA+ CH+", event_col="Death", duration_col="OS from cfDNA collection (mo)")

hr = round(mrcc_cox["HR"], 2)
ci_upper = round(mrcc_cox["CI_upper"], 2)
ci_lower = round(mrcc_cox["CI_lower"], 2)
p = round(float(mrcc_cox["p"]), 5)

kidney_ax.text(1, 0.95, "Median", transform=kidney_ax.transAxes, fontsize=8, ha='center', va='center')
kidney_ax.text(1, 0.88, "mo", transform=kidney_ax.transAxes, fontsize=8, ha='center', va='center', fontstyle = "italic")
kidney_ax.text(0.97, 0.81, 'ctDNA+ CH+', transform=kidney_ax.transAxes, fontsize=8, ha='right', va='center')
kidney_ax.text(0.97, 0.74, 'ctDNA+ CH-', transform=kidney_ax.transAxes, fontsize=8, ha='right', va='center')
kidney_ax.text(1, 0.81, mrcc_dict["Positive"], transform=kidney_ax.transAxes, fontsize=8, ha='center', va='center')
kidney_ax.text(1, 0.74, mrcc_dict["Negative"], transform=kidney_ax.transAxes, fontsize=8, ha='center', va='center')
kidney_ax.text(0.89, 0.67, f"HR: {hr} (95% CI: {ci_upper}-{ci_lower})", transform=kidney_ax.transAxes, fontsize=8, ha='center', va='center')
kidney_ax.text(0.89, 0.60, f"p={p}", transform=kidney_ax.transAxes, fontsize=8, ha='center', va='center')

gs.tight_layout(fig)
fig.savefig(os.path.join(figure_dir, "SUPP_ctDNA_and_CH_survival.png"))
fig.savefig(os.path.join(figure_dir, "SUPP_ctDNA_and_CH_survival.pdf"))