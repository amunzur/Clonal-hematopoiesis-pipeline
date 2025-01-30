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
all_vars_chip = all_vars_chip[(all_vars_chip["Dependent"] == False) & (all_vars_chip["Timepoint"] == "Baseline")].reset_index(drop = True)
base_kidney_chip = all_vars_chip[all_vars_chip["Diagnosis"] == "Kidney"].reset_index(drop = True)
base_bladder_chip = all_vars_chip[all_vars_chip["Diagnosis"] == "Bladder"].reset_index(drop = True)

# LOAD SOMATIC DATASETS
all_vars_somatic = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/somatic.csv"))
all_vars_somatic = all_vars_somatic[(all_vars_somatic["Dependent"] == False) & (all_vars_somatic["Timepoint"] == "Baseline")].reset_index(drop = True)
all_vars_somatic = all_vars_somatic[~all_vars_somatic["Patient_id"].isin(["20-313", "21-184", "21-430"])] # exclude some samples due to oxidative damage
base_kidney_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Bladder")].reset_index(drop = True)

all_pts = pd.read_csv(PATH_sample_information, sep="\t", names=["Patient_id", "Date", "Diagnosis", "Timepoint"])
kidney_pts = all_pts[all_pts["Diagnosis"] == "Kidney"]["Patient_id"].unique()
bladder_pts = all_pts[all_pts["Diagnosis"] == "Bladder"]["Patient_id"].unique()

############################ SET UP FIGURE
fig = plt.figure(figsize=(8, 9))
gs = gridspec.GridSpec(3, 1, height_ratios = [1, 1, 1], hspace = 0, wspace = 0) # outer most gs with 3 rows

# ctDNA KMs
ctDNA_km_gs = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=[1, 1], subplot_spec=gs[0], wspace=0.3, hspace = 0.3)
ctDNA_bladder_KM_ax = plt.subplot(ctDNA_km_gs[0])
# ctDNA_bladder_annotations_ax = plt.subplot(ctDNA_km_gs[1])
ctDNA_kidney_KM_ax = plt.subplot(ctDNA_km_gs[1])
# ctDNA_kidney_annotations_ax = plt.subplot(ctDNA_km_gs[3])

# CH KMs
ch_km_gs = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=[1, 1], subplot_spec=gs[1], wspace=0.3, hspace = 0.3)
ch_bladder_KM_ax = plt.subplot(ch_km_gs[0])
# ch_bladder_annotations_ax = plt.subplot(ch_km_gs[1])
ch_kidney_KM_ax = plt.subplot(ch_km_gs[1])
# ch_kidney_annotations_ax = plt.subplot(ch_km_gs[3])

# irAEs
irAE_gs = gridspec.GridSpecFromSubplotSpec(1, 3, width_ratios=[1, 1, 1], subplot_spec=gs[2], wspace=0.35, hspace = 0.3)
############################

############################ FIGURE 4A: ctDNA KM for BLADDER
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

ctDNA_bladder_KM_ax.set_title("mUC")
ctDNA_bladder_KM_ax = make_survival_curve(
    merged_df_bladder, 
    "ctDNA+", 
    event_col = "Death",
    duration_col = "OS from cfDNA collection (mo)",
    output_path = None, 
    plot_title = "mUC", 
    print_HRs = False, 
    pos_curve_label_positions = [0.15, 0.3], 
    neg_curve_label_positions =  [0.5, 0.7], 
    add_legend = False, 
    print_survival = None,
    xlabel = "OS from cfDNA collection (mo)",
    ax = ctDNA_bladder_KM_ax)

# Calculate median survival to print on the plot
muc_ctDNA_dict = return_median_survival(merged_df_bladder, stratify_by="ctDNA+", event_col="Death", duration_col="OS from cfDNA collection (mo)")
logrank_p = round(float(muc_ctDNA_dict["logrank p"]), 4)

muc_ctdna_cox = run_cox_proportional_hazards(merged_df_bladder, stratify_by="ctDNA+", event_col="Death", duration_col="OS from cfDNA collection (mo)")
hr = round(muc_ctdna_cox["HR"], 2)
ci_upper = round(muc_ctdna_cox["CI_upper"], 2)
ci_lower = round(muc_ctdna_cox["CI_lower"], 2)
cox_p = round(float(muc_ctdna_cox["p"]), 5)

ctDNA_bladder_KM_ax.text(1, 0.95,"Median", transform=ctDNA_bladder_KM_ax.transAxes, fontsize=8, ha='center', va='center')
ctDNA_bladder_KM_ax.text(1, 0.88, "mo", transform=ctDNA_bladder_KM_ax.transAxes, fontsize=8, ha='center', va='center', fontstyle = "italic")
ctDNA_bladder_KM_ax.text(0.97, 0.81, 'ctDNA+', transform=ctDNA_bladder_KM_ax.transAxes, fontsize=8, ha='right', va='center')
ctDNA_bladder_KM_ax.text(0.97, 0.74, 'ctDNA-', transform=ctDNA_bladder_KM_ax.transAxes, fontsize=8, ha='right', va='center')
ctDNA_bladder_KM_ax.text(1, 0.81, muc_ctDNA_dict["Positive"], transform=ctDNA_bladder_KM_ax.transAxes, fontsize=8, ha='center', va='center')
ctDNA_bladder_KM_ax.text(1, 0.74, muc_ctDNA_dict["Negative"], transform=ctDNA_bladder_KM_ax.transAxes, fontsize=8, ha='center', va='center')
# ctDNA_bladder_KM_ax.text(1, 0.67, f"Logrank p={logrank_p}", transform=ctDNA_bladder_KM_ax.transAxes, fontsize=8, ha='center', va='center')
ctDNA_bladder_KM_ax.text(0.86, 0.67, f"HR: {hr} (95% CI: {ci_upper}-{ci_lower})", transform=ctDNA_bladder_KM_ax.transAxes, fontsize=8, ha='center', va='center')
ctDNA_bladder_KM_ax.text(0.86, 0.60, f"p={cox_p}", transform=ctDNA_bladder_KM_ax.transAxes, fontsize=8, ha='center', va='center')

############################ FIGURE 4B: ctDNA KM for KIDNEY
clin_df = pd.read_csv(PATH_kidney_clinical)

sex_and_age_df, cci_df, mets_df, imdc_df, subtype_df, surv_df = prepare_clinical_data(clin_df)
surv_df = surv_df[surv_df["Patient_id"].isin(kidney_pts)]
ctDNA_status = prepare_ctDNA_status(base_kidney_somatic, "Kidney", PATH_sample_information)
chip_status = prepare_chip_status(base_kidney_chip, "Kidney", PATH_sample_information)
multivars_df_binarized = prepare_multivars_df(sex_and_age_df, cci_df, mets_df, ctDNA_status, subtype_df, chip_status, imdc_df)
merged_df_kidney = surv_df.merge(multivars_df_binarized)
del merged_df_kidney["Patient_id"]

ctDNA_kidney_KM_ax = make_survival_curve(
    merged_df_kidney, 
    "ctDNA+",
    event_col = "Death",
    duration_col = "OS from cfDNA collection (mo)",
    output_path = None,
    plot_title = "mRCC", 
    print_HRs = True,
    pos_curve_label_positions = [0.3, 0.25], 
    neg_curve_label_positions =  [0.55, 0.80], 
    add_legend = False,
    print_survival = None,
    xlabel = "OS from cfDNA collection (mo)",
    ax = ctDNA_kidney_KM_ax)

mrcc_ctDNA_dict = return_median_survival(merged_df_kidney, stratify_by="ctDNA+", event_col="Death", duration_col="OS from cfDNA collection (mo)")
mrcc_ctdna_cox = run_cox_proportional_hazards(merged_df_kidney, stratify_by="ctDNA+", event_col="Death", duration_col="OS from cfDNA collection (mo)")

hr = round(mrcc_ctdna_cox["HR"], 2)
ci_upper = round(mrcc_ctdna_cox["CI_upper"], 2)
ci_lower = round(mrcc_ctdna_cox["CI_lower"], 2)
p = round(float(mrcc_ctdna_cox["p"]), 5)

ctDNA_kidney_KM_ax.text(1, 0.95, "Median", transform=ctDNA_kidney_KM_ax.transAxes, fontsize=8, ha='center', va='center')
ctDNA_kidney_KM_ax.text(1, 0.88, "mo", transform=ctDNA_kidney_KM_ax.transAxes, fontsize=8, ha='center', va='center', fontstyle = "italic")
ctDNA_kidney_KM_ax.text(0.97, 0.81, 'ctDNA+', transform=ctDNA_kidney_KM_ax.transAxes, fontsize=8, ha='right', va='center')
ctDNA_kidney_KM_ax.text(0.97, 0.74, 'ctDNA-', transform=ctDNA_kidney_KM_ax.transAxes, fontsize=8, ha='right', va='center')
ctDNA_kidney_KM_ax.text(1, 0.81, mrcc_ctDNA_dict["Positive"], transform=ctDNA_kidney_KM_ax.transAxes, fontsize=8, ha='center', va='center')
ctDNA_kidney_KM_ax.text(1, 0.74, mrcc_ctDNA_dict["Negative"], transform=ctDNA_kidney_KM_ax.transAxes, fontsize=8, ha='center', va='center')
ctDNA_kidney_KM_ax.text(0.89, 0.67, f"HR: {hr} (95% CI: {ci_upper}-{ci_lower})", transform=ctDNA_kidney_KM_ax.transAxes, fontsize=8, ha='center', va='center')
ctDNA_kidney_KM_ax.text(0.89, 0.60, f"p={p}", transform=ctDNA_kidney_KM_ax.transAxes, fontsize=8, ha='center', va='center')


############################ FIGURE 2C: CH KM for BLADDER
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

ch_bladder_KM_ax.set_title("mUC")
del merged_df_bladder["Patient_id"]

ch_bladder_KM_ax = make_survival_curve(
    merged_df_bladder, 
    "CH+", 
    event_col = "Death",
    duration_col = "OS from cfDNA collection (mo)",
    output_path = None, 
    plot_title = "mUC", 
    print_HRs = False, 
    pos_curve_label_positions = [0.15, 0.4], 
    neg_curve_label_positions =  [0.3, 0.8], 
    add_legend = False, 
    print_survival = None,
    xlabel = "OS from cfDNA collection (mo)",
    ax = ch_bladder_KM_ax)

muc_ch_dict = return_median_survival(merged_df_bladder, stratify_by="CH+", event_col="Death", duration_col="OS from cfDNA collection (mo)")
muc_ch_cox = run_cox_proportional_hazards(merged_df_bladder, stratify_by="CH+", event_col="Death", duration_col="OS from cfDNA collection (mo)")

hr = round(muc_ch_cox["HR"], 2)
ci_upper = round(muc_ch_cox["CI_upper"], 2)
ci_lower = round(muc_ch_cox["CI_lower"], 2)
cox_p = round(float(muc_ch_cox["p"]), 5)

ch_bladder_KM_ax.text(1, 0.95, "Median", transform=ch_bladder_KM_ax.transAxes, fontsize=8, ha='center', va='center')
ch_bladder_KM_ax.text(1, 0.88, "mo", transform=ch_bladder_KM_ax.transAxes, fontsize=8, ha='center', va='center', fontstyle = "italic")
ch_bladder_KM_ax.text(0.97, 0.81, 'CH+', transform=ch_bladder_KM_ax.transAxes, fontsize=8, ha='right', va='center')
ch_bladder_KM_ax.text(0.97, 0.74, 'CH-', transform=ch_bladder_KM_ax.transAxes, fontsize=8, ha='right', va='center')
ch_bladder_KM_ax.text(1, 0.81, muc_ch_dict["Positive"], transform=ch_bladder_KM_ax.transAxes, fontsize=8, ha='center', va='center')
ch_bladder_KM_ax.text(1, 0.74, muc_ch_dict["Negative"], transform=ch_bladder_KM_ax.transAxes, fontsize=8, ha='center', va='center')
ch_bladder_KM_ax.text(0.89, 0.67, f"HR: {hr} (95% CI: {ci_upper}-{ci_lower})", transform=ch_bladder_KM_ax.transAxes, fontsize=8, ha='center', va='center')
ch_bladder_KM_ax.text(0.89, 0.60, f"p={cox_p}", transform=ch_bladder_KM_ax.transAxes, fontsize=8, ha='center', va='center')

############################ FIGURE 3D: CH KM for KIDNEY
clin_df = pd.read_csv(PATH_kidney_clinical)

sex_and_age_df, cci_df, mets_df, imdc_df, subtype_df, surv_df = prepare_clinical_data(clin_df)
surv_df = surv_df[surv_df["Patient_id"].isin(kidney_pts)]
ctDNA_status = prepare_ctDNA_status(base_kidney_somatic, "Kidney", PATH_sample_information)
chip_status = prepare_chip_status(base_kidney_chip, "Kidney", PATH_sample_information)
multivars_df_binarized = prepare_multivars_df(sex_and_age_df, cci_df, mets_df, ctDNA_status, subtype_df, chip_status, imdc_df)
merged_df_kidney = surv_df.merge(multivars_df_binarized)
# kidney_ch_ax_combined = plot_ctDNA_and_CH_KMs(ctDNA_status, chip_status, merged_df_kidney, PATH_sample_information, kidney_ch_ax_combined, xlabel = None, show_legend = False, xmax = None, plot_title = "")

del merged_df_kidney["Patient_id"]

ch_kidney_KM_ax = make_survival_curve(
    merged_df_kidney, 
    "CH+",
    event_col = "Death",
    duration_col = "OS from cfDNA collection (mo)",
    output_path = None,
    plot_title = "mRCC", 
    print_HRs = False,
    pos_curve_label_positions = [0.47, 0.65], 
    neg_curve_label_positions =  [0.28, 0.42], 
    add_legend = False, 
    print_survival = None,
    xlabel = "OS from cfDNA collection (mo)",
    ax = ch_kidney_KM_ax)

mrcc_ch_dict = return_median_survival(merged_df_kidney, stratify_by="CH+", event_col="Death", duration_col="OS from cfDNA collection (mo)")
mrcc_ch_cox = run_cox_proportional_hazards(merged_df_kidney, stratify_by="CH+", event_col="Death", duration_col="OS from cfDNA collection (mo)")
hr = round(mrcc_ch_cox["HR"], 2)
ci_upper = round(mrcc_ch_cox["CI_upper"], 2)
ci_lower = round(mrcc_ch_cox["CI_lower"], 2)
cox_p = round(float(mrcc_ch_cox["p"]), 5)

ch_kidney_KM_ax.text(1, 0.95, "Median", transform=ch_kidney_KM_ax.transAxes, fontsize=8, ha='center', va='center')
ch_kidney_KM_ax.text(1, 0.88, "mo", transform=ch_kidney_KM_ax.transAxes, fontsize=8, ha='center', va='center', fontstyle = "italic")
ch_kidney_KM_ax.text(0.97, 0.81, 'CH+', transform=ch_kidney_KM_ax.transAxes, fontsize=8, ha='right', va='center')
ch_kidney_KM_ax.text(0.97, 0.74, 'CH-', transform=ch_kidney_KM_ax.transAxes, fontsize=8, ha='right', va='center')
ch_kidney_KM_ax.text(1, 0.81, mrcc_ch_dict["Positive"], transform=ch_kidney_KM_ax.transAxes, fontsize=8, ha='center', va='center')
ch_kidney_KM_ax.text(1, 0.74, mrcc_ch_dict["Negative"], transform=ch_kidney_KM_ax.transAxes, fontsize=8, ha='center', va='center')
ch_kidney_KM_ax.text(0.89, 0.67, f"HR: {hr} (95% CI: {ci_upper}-{ci_lower})", transform=ch_kidney_KM_ax.transAxes, fontsize=8, ha='center', va='center')
ch_kidney_KM_ax.text(0.89, 0.60, f"p={cox_p}", transform=ch_kidney_KM_ax.transAxes, fontsize=8, ha='center', va='center')

#############################################################################

# Now we plot the irAE plots. We look into 3 categories.
############################ SET UP FIGURE
# fig = plt.figure(figsize=(8, 8))
# gs = gridspec.GridSpec(3, 1, height_ratios = [1, 1, 1], hspace = 0, wspace = 0) # outer most gs with 3 rows

# # ctDNA KMs
# ctDNA_km_gs = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=[1, 1], subplot_spec=gs[0], wspace=0.3, hspace = 0.3)
# ctDNA_bladder_KM_ax = plt.subplot(ctDNA_km_gs[0])
# # ctDNA_bladder_annotations_ax = plt.subplot(ctDNA_km_gs[1])
# ctDNA_kidney_KM_ax = plt.subplot(ctDNA_km_gs[1])
# # ctDNA_kidney_annotations_ax = plt.subplot(ctDNA_km_gs[3])

# # CH KMs
# ch_km_gs = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=[1, 1], subplot_spec=gs[1], wspace=0.3, hspace = 0.3)
# ch_bladder_KM_ax = plt.subplot(ch_km_gs[0])
# # ch_bladder_annotations_ax = plt.subplot(ch_km_gs[1])
# ch_kidney_KM_ax = plt.subplot(ch_km_gs[1])
# # ch_kidney_annotations_ax = plt.subplot(ch_km_gs[3])

# irAEs
# irAE_gs = gridspec.GridSpecFromSubplotSpec(1, 4, width_ratios=[1, 1, 1, 0.05], subplot_spec=gs[2], wspace=0.35, hspace = 0.3)
############################




ax1 = plt.subplot(irAE_gs[0])
ax2 = plt.subplot(irAE_gs[1])
ax3 = plt.subplot(irAE_gs[2])
irae_legend_ax = plt.subplot(irAE_gs[3])

# irAEs
ax1, ax2, ax3 = make_irAE_association_plots(all_vars_chip, PATH_kidney_clinical, PATH_sample_information, ax1, ax2, ax3, irae_legend_ax)
# fig.savefig(os.path.join(figure_dir, "figure4.png"))
# fig.savefig(os.path.join(figure_dir, "figure4.pdf"))


fig.text(0.07, 0.97, 'a', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
fig.text(0.52, 0.97, 'b', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
fig.text(0.07, 0.60, 'c', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
fig.text(0.52, 0.60, 'd', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
# fig.text(0.07, 0.25, 'E', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)

gs.tight_layout(fig)
fig.savefig(os.path.join(figure_dir, "figure4.png"))
fig.savefig(os.path.join(figure_dir, "figure4.pdf"))