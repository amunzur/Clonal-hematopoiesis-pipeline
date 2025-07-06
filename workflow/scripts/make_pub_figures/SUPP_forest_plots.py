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
import matplotlib.patches as patches

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

def generate_genomic_subgroups(patients_list1, patients_list2):
	"""
	Checks to see if patients_list1 is in patients_list2. Returns a 2 col df with patient IDs and binary outcomes.
	patients_list1: All patients
	patients_list2: Patients with some sort of mutation
	"""
	df = pd.DataFrame({
		'Patient_id': patients_list1,
		'Status': pd.Series(patients_list1).isin(patients_list2).astype(int)  # 1 if in patients_list2, 0 if not
		})
	return(df)

def calculate_subgroup_HRs(surv_df, genomic_stratification_df, event_col, duration_col, stratify_by = "Status"):
    """
    Runs univariate subgroup analysis, incorporating both PFS and OS. 
    Stratified by: 
    ctDNA, CH, DDR CH, DTA CH, non-DNMT3A CH, DTA CH
    """
    merged = surv_df.merge(genomic_stratification_df, how = "inner")
    cox_dict = run_cox_proportional_hazards(merged, stratify_by=stratify_by, event_col=event_col, duration_col=duration_col)
    return(cox_dict)

def plot_single_HR(ax, ypos, cox_dict, color):
    """
    Given the output of the calculate_subgroup_HRs function, plots it on the ax.
    """
    hr = cox_dict["HR"]
    ci_upper = cox_dict["CI_upper"]
    ci_lower = cox_dict["CI_lower"]
    p = float(cox_dict["p"])
    
    # Calculate the error margins for the confidence intervals
    lower_error = hr - ci_lower
    upper_error = ci_upper - hr
    asymmetric_error = [[lower_error], [upper_error]]
    
    # Plotting
    # ax.scatter(hr, ypos, s=1, edgecolor=None, color=color)
    ax.errorbar(hr, ypos, xerr=asymmetric_error, fmt='o', color=color, capsize=2, elinewidth = 1, markersize = 4)
    
    return ax

def get_count_dict(status_df, clin_df = None):
    # Subset
    if clin_df is not None:
        status_df = status_df[status_df["Patient_id"].isin(clin_df["Patient_id"])].reset_index(drop = True)
    
    status_dict = status_df["Status"].value_counts().reset_index().set_index('index').to_dict()['Status']
    # Map 0 to 'negative' and 1 to 'positive'
    updated_dict = {'positive': status_dict.get(1, 0), 'negative': status_dict.get(0, 0)}
    return(updated_dict)



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

# Generate survival DFs from bladder and kidney
kidney_clin_df = pd.read_csv(PATH_kidney_clinical)
sex_and_age_df, cci_df, mets_df, imdc_df, subtype_df, kidney_surv_df = prepare_clinical_data(kidney_clin_df)
kidney_surv_df = kidney_surv_df[kidney_surv_df["Patient_id"].isin(kidney_pts)]

bladder_clin_df = pd.read_csv(PATH_clinical_bladder)
bladder_clin_df = bladder_clin_df[(bladder_clin_df["First sample?"] == True) & (bladder_clin_df["Patient_id"].isin(bladder_pts))]
bladder_surv_df = prepare_survival_data(bladder_clin_df)

# Generate relevant subgroups genomic dfs
# Bladder
muc_ctdna_status = generate_genomic_subgroups(bladder_pts, patients_list2 = base_bladder_somatic["Patient_id"].unique())
muc_ch_status = generate_genomic_subgroups(bladder_pts, patients_list2 = base_bladder_chip["Patient_id"].unique())
muc_chip_status = generate_genomic_subgroups(bladder_pts, patients_list2 = base_bladder_chip[base_bladder_chip["VAF_n"] >= 2]["Patient_id"].unique())
muc_ddr_ch_status = generate_genomic_subgroups(bladder_pts, patients_list2 = base_bladder_chip[base_bladder_chip["Gene"].isin(["TP53", "BRCA1", "BRCA2", "PPM1D", "ATM", "CHEK2", "ARID1A"])]["Patient_id"].unique())
muc_dta_ch_status = generate_genomic_subgroups(bladder_pts, patients_list2 = base_bladder_chip[base_bladder_chip["Gene"].isin(["DNMT3A", "TET2", "ASXL1"])]["Patient_id"].unique())
muc_splicing_ch_status = generate_genomic_subgroups(bladder_pts, patients_list2= base_bladder_chip[base_bladder_chip["Gene"].isin(["SF3B1", "SRSF2", "U2AF1", "ZRSR2"])]["Patient_id"].unique())
muc_nondta_ch_status = generate_genomic_subgroups(bladder_pts, patients_list2 = base_bladder_chip[~base_bladder_chip["Gene"].isin(["DNMT3A", "TET2", "ASXL1"])]["Patient_id"].unique())

# Kidney
rcc_ctdna_status = generate_genomic_subgroups(kidney_pts, patients_list2 = base_kidney_somatic["Patient_id"].unique())
rcc_ch_status = generate_genomic_subgroups(kidney_pts, patients_list2 = base_kidney_chip["Patient_id"].unique())
rcc_chip_status = generate_genomic_subgroups(kidney_pts, patients_list2 = base_kidney_chip[base_kidney_chip["VAF_n"] >= 2]["Patient_id"].unique())
rcc_ddr_ch_status = generate_genomic_subgroups(kidney_pts, patients_list2 = base_kidney_chip[base_kidney_chip["Gene"].isin(["TP53", "BRCA1", "BRCA2", "PPM1D", "ATM", "CHEK2", "ARID1A"])]["Patient_id"].unique())
rcc_dta_ch_status = generate_genomic_subgroups(kidney_pts, patients_list2 = base_kidney_chip[base_kidney_chip["Gene"].isin(["DNMT3A", "TET2", "ASXL1"])]["Patient_id"].unique())
rcc_splicing_ch_status = generate_genomic_subgroups(kidney_pts, patients_list2= base_kidney_chip[base_kidney_chip["Gene"].isin(["SF3B1", "SRSF2", "U2AF1", "ZRSR2"])]["Patient_id"].unique())
rcc_nondta_ch_status = generate_genomic_subgroups(kidney_pts, patients_list2 = base_kidney_chip[~base_kidney_chip["Gene"].isin(["DNMT3A", "TET2", "ASXL1"])]["Patient_id"].unique())

# Run subgroup analysis
# Bladder
muc_ctdna_cox = calculate_subgroup_HRs(surv_df = bladder_surv_df, genomic_stratification_df = muc_ctdna_status, event_col = "Death", duration_col = "OS from cfDNA collection (mo)")
muc_ch_cox = calculate_subgroup_HRs(surv_df = bladder_surv_df, genomic_stratification_df = muc_ch_status, event_col = "Death", duration_col = "OS from cfDNA collection (mo)")
muc_chip_cox = calculate_subgroup_HRs(surv_df = bladder_surv_df, genomic_stratification_df = muc_chip_status, event_col = "Death", duration_col = "OS from cfDNA collection (mo)")
muc_ddr_ch_cox = calculate_subgroup_HRs(surv_df = bladder_surv_df, genomic_stratification_df = muc_ddr_ch_status, event_col = "Death", duration_col = "OS from cfDNA collection (mo)")
muc_dta_ch_cox = calculate_subgroup_HRs(surv_df = bladder_surv_df, genomic_stratification_df = muc_dta_ch_status, event_col = "Death", duration_col = "OS from cfDNA collection (mo)")
muc_splicing_ch_cox = calculate_subgroup_HRs(surv_df = bladder_surv_df, genomic_stratification_df = muc_splicing_ch_status, event_col = "Death", duration_col = "OS from cfDNA collection (mo)")
muc_nondta_ch_cox = calculate_subgroup_HRs(surv_df = bladder_surv_df, genomic_stratification_df = muc_nondta_ch_status, event_col = "Death", duration_col = "OS from cfDNA collection (mo)")

# Kidney
mrcc_ctdna_cox = calculate_subgroup_HRs(surv_df = kidney_surv_df, genomic_stratification_df = rcc_ctdna_status, event_col = "Death", duration_col = "OS from cfDNA collection (mo)")
mrcc_ch_cox = calculate_subgroup_HRs(surv_df = kidney_surv_df, genomic_stratification_df = rcc_ch_status, event_col = "Death", duration_col = "OS from cfDNA collection (mo)")
mrcc_chip_cox = calculate_subgroup_HRs(surv_df = kidney_surv_df, genomic_stratification_df = rcc_chip_status, event_col = "Death", duration_col = "OS from cfDNA collection (mo)")
mrcc_ddr_ch_cox = calculate_subgroup_HRs(surv_df = kidney_surv_df, genomic_stratification_df = rcc_ddr_ch_status, event_col = "Death", duration_col = "OS from cfDNA collection (mo)")
mrcc_dta_ch_cox = calculate_subgroup_HRs(surv_df = kidney_surv_df, genomic_stratification_df = rcc_dta_ch_status, event_col = "Death", duration_col = "OS from cfDNA collection (mo)")
mrcc_splicing_ch_cox = calculate_subgroup_HRs(surv_df = kidney_surv_df, genomic_stratification_df = rcc_splicing_ch_status, event_col = "Death", duration_col = "OS from cfDNA collection (mo)")
mrcc_nondta_ch_cox = calculate_subgroup_HRs(surv_df = kidney_surv_df, genomic_stratification_df = rcc_nondta_ch_status, event_col = "Death", duration_col = "OS from cfDNA collection (mo)")

# Calculate ns for each subgroup
muc_ctdna_count = get_count_dict(muc_ctdna_status)
muc_ch_count = get_count_dict(muc_ch_status)
muc_chip_count = get_count_dict(muc_chip_status)
muc_ddr_ch_count = get_count_dict(muc_ddr_ch_status)
muc_dta_ch_count = get_count_dict(muc_dta_ch_status)
muc_splicing_ch_count = get_count_dict(muc_splicing_ch_status)
muc_nondta_ch_count = get_count_dict(muc_nondta_ch_status)

rcc_ctdna_count = get_count_dict(rcc_ctdna_status, clin_df = kidney_surv_df)
rcc_ch_count = get_count_dict(rcc_ch_status, clin_df = kidney_surv_df)
rcc_chip_count = get_count_dict(rcc_chip_status, clin_df = kidney_surv_df)
rcc_ddr_ch_count = get_count_dict(rcc_ddr_ch_status, clin_df = kidney_surv_df)
rcc_dta_ch_count = get_count_dict(rcc_dta_ch_status, clin_df = kidney_surv_df)
rcc_splicing_ch_count = get_count_dict(rcc_splicing_ch_status, clin_df = kidney_surv_df)
rcc_nondta_ch_count = get_count_dict(rcc_nondta_ch_status, clin_df = kidney_surv_df)

# Put into lists
cox_list_bladder =  [muc_splicing_ch_cox, muc_ddr_ch_cox, muc_dta_ch_cox, muc_nondta_ch_cox, muc_chip_cox, muc_ch_cox, muc_ctdna_cox]
cox_list_kidney =  [mrcc_splicing_ch_cox, mrcc_ddr_ch_cox, mrcc_dta_ch_cox, mrcc_nondta_ch_cox, mrcc_chip_cox, mrcc_ch_cox, mrcc_ctdna_cox]
    
counts_list_bladder =  [muc_splicing_ch_count, muc_ddr_ch_count, muc_dta_ch_count, muc_nondta_ch_count, muc_chip_count,muc_ch_count,  muc_ctdna_count]
counts_list_kidney =  [rcc_splicing_ch_count, rcc_ddr_ch_count, rcc_dta_ch_count, rcc_nondta_ch_count, rcc_chip_count, rcc_ch_count, rcc_ctdna_count]

############################ SET UP FIGURE
fig = plt.figure(figsize=(8, 4))
gs = gridspec.GridSpec(2, 6, height_ratios = [0.1, 1], width_ratios = [0.25, 0.4, 0.2, 0.1, 0.05, 0.05], hspace = 0, wspace = 0) # outer most gs with 3 rows

forest_ax =  plt.subplot(gs[1, 1])
annotation_ax =  plt.subplot(gs[1, 0], sharey = forest_ax)
hr_ax =  plt.subplot(gs[1, 2], sharey = forest_ax)
p_ax =  plt.subplot(gs[1, 3], sharey = forest_ax)
counts_ax_pos = plt.subplot(gs[1, 4], sharey = forest_ax)
counts_ax_neg = plt.subplot(gs[1, 5], sharey = forest_ax)

# add subplots for titles ax
annot_title_ax = plt.subplot(gs[0, 0])
forest_title_ax = plt.subplot(gs[0, 1])
hr_title_ax = plt.subplot(gs[0, 2])
p_title_ax = plt.subplot(gs[0, 3])
n_pos_title_ax = plt.subplot(gs[0, 4])
n_neg_title_ax = plt.subplot(gs[0, 5])

bladder_color ="deepskyblue"
kidney_color = "orangered"

for i, cox_dict in enumerate(cox_list_kidney):
    forest_ax = plot_single_HR(forest_ax, i-0.1, cox_dict, color = "orangered")

for i, cox_dict in enumerate(cox_list_bladder):
    forest_ax = plot_single_HR(forest_ax, i+0.1, cox_dict, color = "deepskyblue")

# Aes
forest_ax.spines[["top", "right", "left"]].set_visible(False)
forest_ax.axvline(1, linestyle='--', color='black', lw = 0.5)
forest_ax.set_xticks([-1, 0, 1, 2, 3, 4])
forest_ax.set_xticklabels(["-1", "0", "1", "2", "3", "4"])
forest_ax.set_xlabel("HR")

# Background grey rectangles for all axes
for ax in [annotation_ax, forest_ax, hr_ax, p_ax, counts_ax_pos, counts_ax_neg]:
    for i in list(range(len(cox_list_kidney))):
        if i % 2 == 0:  # Check if i is even
            x_limits = ax.get_xlim()  # Get the current x limits of the axis
            rect = patches.Rectangle(
                (x_limits[0], i - 0.5),  # Bottom left corner of the rectangle
                x_limits[1] - x_limits[0],  # Width spans the entire x limits
                1,                          # Height of the rectangle (you may adjust this)
                color='lightgrey',         # Color of the rectangle
                alpha=0.5,                # Transparency of the rectangle
                zorder=0                   # Ensure it's drawn behind other elements
            )
            ax.add_patch(rect)

# Remove all axis elements
for ax in [annotation_ax, hr_ax, p_ax, counts_ax_pos, counts_ax_neg, annot_title_ax, forest_title_ax, hr_title_ax, p_title_ax, n_pos_title_ax, n_neg_title_ax]:
    ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
    ax.tick_params(axis="both", direction="out", which="both", left=False, bottom=False, labelbottom = False, labelleft = False)
    ax.set_xlim((0, 1))

# Add annotations
annots_list = ["Spliceosome CH", "DDR CH", "DTA CH", "non-DTA CH", "CH≥2%", "CH", "ctDNA"]
for i, annot in enumerate(annots_list):
    annotation_ax.text(1, i, annot, ha='right', va='center', fontsize=9)

# Add HRs
for i, cox_dict in enumerate(cox_list_kidney):
    hr = round(cox_dict["HR"], 2)
    ci_upper = round(cox_dict["CI_upper"], 2)
    ci_lower = round(cox_dict["CI_lower"], 2)
    p = round(float(cox_dict["p"]), 5)
    
    HR_text = f"{hr} ({ci_lower}-{ci_upper})"
    p_text = str(p)
    
    hr_ax.text(0, i-0.1, HR_text, ha='left', va='center', fontsize=7, color = "orangered")
    p_ax.text(0, i-0.1, p_text, ha='left', va='center', fontsize=7, color = "orangered")

for i, cox_dict in enumerate(cox_list_bladder):
    hr = round(cox_dict["HR"], 2)
    ci_upper = round(cox_dict["CI_upper"], 2)
    ci_lower = round(cox_dict["CI_lower"], 2)
    p = round(float(cox_dict["p"]), 5)
    
    HR_text = f"{hr} ({ci_lower}-{ci_upper})"
    p_text = str(p)
    
    hr_ax.text(0, i+0.1, HR_text, ha='left', va='center', fontsize=7, color = "deepskyblue")
    p_ax.text(0, i+0.1, p_text, ha='left', va='center', fontsize=7, color = "deepskyblue")

# Annotate the ns
for i, counts_dict in enumerate(counts_list_kidney):
    counts_ax_pos.text(0, i-0.1, counts_dict["positive"], ha='center', va='center', fontsize=7, color = "orangered")
    counts_ax_neg.text(0, i-0.1, counts_dict["negative"], ha='center', va='center', fontsize=7, color = "orangered")

for i, counts_dict in enumerate(counts_list_bladder):
    counts_ax_pos.text(0, i+0.1, counts_dict["positive"], ha='center', va='center', fontsize=7, color = "deepskyblue")
    counts_ax_neg.text(0, i+0.1, counts_dict["negative"], ha='center', va='center', fontsize=7, color = "deepskyblue")

# add titles
titles_list = ["Group", "mUC/mRCC", "HR (95% CI)", "p", "n+ ", "n-"]
for ax, title in zip([annot_title_ax, forest_title_ax, hr_title_ax, p_title_ax, n_pos_title_ax, n_neg_title_ax], titles_list):
    ax.set_ylim((0, 1))
    ax.text(0, 0.5, title, ha='left', va='center', fontsize=9, color = "black")

fig.savefig(os.path.join(figure_dir, "SUPP_forest.png"))
fig.savefig(os.path.join(figure_dir, "SUPP_forest.pdf"))