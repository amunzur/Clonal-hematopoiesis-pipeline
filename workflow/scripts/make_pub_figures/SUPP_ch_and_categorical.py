import pandas as pd
import numpy as np
import os 
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
from scipy.stats import fisher_exact
from scipy.stats import fisher_exact # needs to stay twice
from statsmodels.stats.multitest import multipletests

DIR_working = "/groups/wyattgrp/users/amunzur/pipeline"
PATH_sample_information = os.path.join(DIR_working, "resources/sample_lists/sample_information.tsv")
PATH_kidney_clinical = os.path.join(DIR_working, "resources/clinical_data/Supplementary tables - Clinical data - mRCC.csv")
PATH_clinical_bladder = os.path.join(DIR_working, "resources/clinical_data/bladder/clinical_data.csv")
figure_dir = os.path.join(DIR_working, "results/figures/amazing_figures")
source_functions = os.path.join(DIR_working, "workflow/scripts/visualization/UTILITIES_make_chip_plots.py")
sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])

color_dict = {"Bladder": "deepskyblue", "Kidney": "orangered"}

with open(source_functions, 'r') as file:
    script_code = file.read()

exec(script_code)

# LOAD CHIP DATASETS
all_vars_chip = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/chip.csv"))
all_vars_chip = all_vars_chip[all_vars_chip["Dependent"] == False]
base_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)

kidney_clin = pd.read_csv(PATH_kidney_clinical)
bladder_clin = pd.read_csv(PATH_clinical_bladder)
kidney_pts = sample_info[(sample_info["Diagnosis"] == "Kidney") & (sample_info["Timepoint"] == "Baseline")]["Patient_id"].tolist()
bladder_pts = sample_info[(sample_info["Diagnosis"] == "Bladder") & (sample_info["Timepoint"] == "Baseline")]["Patient_id"].tolist()
kidney_clin = kidney_clin[kidney_clin["Patient_id"].isin(kidney_pts)]
bladder_clin = bladder_clin[bladder_clin["First sample?"]]
bladder_clin = bladder_clin[bladder_clin["Patient_id"].isin(bladder_pts)]

fig = plt.figure(figsize=(6, 5))
outer_gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 1], hspace = 0, wspace = 0)
inner_gs_0 = gridspec.GridSpecFromSubplotSpec(1, 6, width_ratios=[1, 1, 1, 1, 1, 0.2], subplot_spec=outer_gs[0], wspace=0.5, hspace = 0.3)
inner_gs_1 = gridspec.GridSpecFromSubplotSpec(1, 6, width_ratios=[1, 1, 1, 1, 1, 0.2], subplot_spec=outer_gs[1], wspace=0.5, hspace = 0.3)
inner_gs_2 = gridspec.GridSpecFromSubplotSpec(1, 6, width_ratios=[1, 1, 1, 1, 1, 0.2], subplot_spec=outer_gs[2], wspace=0.5, hspace = 0.3)

# Prepare the initial diagnosis df
kidney_initial = kidney_clin[["Patient_id", "Disease stage at initial diagnosis"]]
bladder_initial = bladder_clin[["Patient_id", "Disease stage at initial diagnosis"]]
bladder_initial["Disease stage at initial diagnosis"] = bladder_initial["Disease stage at initial diagnosis"].str.replace(r"\s*MIBC|\s*NMIBC", "", regex=True)

# for bladder prepare df for mibc vs nmibc df.
bladder_staging = bladder_clin[bladder_clin["Disease stage at initial diagnosis"].isin(["Localized NMIBC", "Localized MIBC"])][["Patient_id", "Disease stage at initial diagnosis"]].reset_index(drop = True)
sex_df = pd.concat([kidney_clin[["Patient_id", "Sex"]], bladder_clin[["Patient_id", "Sex"]]]).reset_index(drop = True)
smoking_df = bladder_clin[["Patient_id", "Previous smoking history"]].replace({"Current smoker": "Previous/current smoker", "Previous smoker": "Previous/current smoker"})
initial_diagnosis_df = pd.concat([kidney_initial, bladder_initial]).reset_index(drop = True)

x = pd.concat([base_kidney_chip, base_bladder_chip])
y = x[x["VAF_n"] >= 2]
z = x[x["VAF_n"] >= 10]
inner_gs_0 = ch_presence_absence_bars(x, PATH_sample_information, sex_df, smoking_df, initial_diagnosis_df, bladder_staging, figure_dir, fig_name = "sex_diagnosis_and_CH_presence", gs = inner_gs_0, title = "WBC VAF >= 0.25%", show_x_ticks = False, show_legend = True) # all ch mutations
inner_gs_1 = ch_presence_absence_bars(y, PATH_sample_information, sex_df, smoking_df, initial_diagnosis_df, bladder_staging, figure_dir, fig_name = "sex_diagnosis_and_CH_presence", gs = inner_gs_1, title = "WBC VAF >= 2%", show_ax_titles = False, show_x_ticks = False, show_legend = False) # ch>2
inner_gs_2 = ch_presence_absence_bars(z, PATH_sample_information, sex_df, smoking_df, initial_diagnosis_df, bladder_staging, figure_dir, fig_name = "sex_diagnosis_and_CH_presence", gs = inner_gs_2, title = "WBC VAF >= 10%", show_ax_titles = False, show_legend = False) # ch>3

fig.text(0.02, 0.98, 'A', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
fig.text(0.02, 0.66, 'B', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
fig.text(0.02, 0.33, 'C', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)

outer_gs.tight_layout(fig)
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/SUPP_sex_diagnosis_and_CH_presence.png")
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/SUPP_sex_diagnosis_and_CH_presence.pdf")

# compare ages between ppl originally diagnosed with localized disease vs metastatic disease
age_df = pd.concat([(bladder_clin[["Patient_id", "Age at blood draw"]]), kidney_clin[["Patient_id", "Age at GUBB draw"]].rename(columns = {"Age at GUBB draw": "Age at blood draw"})])
df = pd.concat([kidney_initial, bladder_initial]).merge(age_df)
loc = df[df["Disease stage at initial diagnosis"] == "Localized"]["Age at blood draw"]
mets = df[df["Disease stage at initial diagnosis"] == "Metastatic"]["Age at blood draw"]
mannwhitneyu(loc, mets)
np.median(loc)
np.median(mets)

