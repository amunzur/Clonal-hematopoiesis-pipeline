"""
Divides pts into 3 groups and plots 3 curves.

ctDNA negative
ctDNA positive
ctDNA positive & CH positive

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

def plot_ctDNA_and_CH_KMs(ctDNA_status, chip_status, surv_df_merged, PATH_sample_information, ax, xlabel = None, show_legend = False, xmax = None, plot_title = ""):
    """
    Plots the ctDNA and CH combined KMs. 
    Has three curves: 
    ctDNA negative
    ctDNA positive
    ctDNA and CH positive
    """
    combined = ctDNA_status.merge(chip_status)
    
    # Generate three groups of patients as well
    combined = ctDNA_status.merge(chip_status)
    both_positive_pts = combined[(combined["ctDNA positive"] == True) & (combined["CHIP positive"] == True)]["Patient_id"]
    ctDNA_positive_ch_negative_pts = combined[(combined["ctDNA positive"] == True) & (combined["CHIP positive"] == False)]["Patient_id"]
    ctDNA_negative_pts = combined[combined["ctDNA positive"] == False]["Patient_id"]
    
    both_positive_df = surv_df_merged[surv_df_merged["Patient_id"].isin(both_positive_pts)]
    ctDNA_positive_ch_negative_df = surv_df_merged[surv_df_merged["Patient_id"].isin(ctDNA_positive_ch_negative_pts)]
    ctDNA_negative_df = surv_df_merged[surv_df_merged["Patient_id"].isin(ctDNA_negative_pts)]
    
    del both_positive_df["Patient_id"]
    del ctDNA_positive_ch_negative_df["Patient_id"]
    del ctDNA_negative_df["Patient_id"]
    
    # Fit and plot the KMs:
    kmf_both_positive = KaplanMeierFitter()
    kmf_ctDNA_positive_ch_negative = KaplanMeierFitter()
    kmf_ctDNA_negative = KaplanMeierFitter()
    
    kmf_both_positive.fit(durations=both_positive_df["OS from cfDNA collection (mo)"], event_observed=both_positive_df["Death"], label="ctDNA+ CH+")
    kmf_both_positive.plot_survival_function(ax=ax, ci_show=False, show_censors=True, censor_styles={'marker': '|', 'ms': 7, 'markerfacecolor': "mediumpurple"}, color = "mediumpurple")
    
    kmf_ctDNA_positive_ch_negative.fit(durations=ctDNA_positive_ch_negative_df["OS from cfDNA collection (mo)"], event_observed=ctDNA_positive_ch_negative_df["Death"], label="ctDNA+ CH-")
    kmf_ctDNA_positive_ch_negative.plot_survival_function(ax=ax, ci_show=False, show_censors=True, censor_styles={'marker': '|', 'ms': 7, 'markerfacecolor': "forestgreen"}, color = "forestgreen")
    
    kmf_ctDNA_negative.fit(durations=ctDNA_negative_df["OS from cfDNA collection (mo)"], event_observed=ctDNA_negative_df["Death"], label="ctDNA-")
    kmf_ctDNA_negative.plot_survival_function(ax=ax, ci_show=False, show_censors=True, censor_styles={'marker': '|', 'ms': 7, 'markerfacecolor': "black"}, color = "black")
    
    # AES
    if xlabel is None:
        ax.set_xlabel("Months since cfDNA collection")
    else: 
        ax.set_xlabel(xlabel)
    ax.set_ylabel("Survival fraction")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if show_legend:
        ax.legend(loc='best', frameon=False)
    else:
        ax.legend_.remove()
    ax.set_ylim((0, 1))
    ax.set_yticks((0, 0.2, 0.4, 0.6, 0.8, 1))
    ax.set_yticklabels(("0", "0.2", "0.4", "0.6", "0.8", "1"))
    ax.set_title(plot_title)
    if xmax is None:
        xmax = ax.get_xlim()[1]
    ax.set_xlim((0, xmax))
    ax.set_xticks(np.arange(0, xmax + 10, 10))
    if xmax >= 100:
        ax.set_xticks(np.arange(0, xmax + 10, 20))
    if xmax < 40:
        ax.set_xticks(np.arange(0, xmax + 10, 5))
    
    # add risk counts
    add_at_risk_counts(kmf_both_positive, kmf_ctDNA_positive_ch_negative, kmf_ctDNA_negative, ax=ax, rows_to_show = ["At risk"])
    
    # Annotate the median survival values
    median_both_positive = round(kmf_both_positive.median_survival_time_)
    median_ctDNA_positive_ch_negative = round(kmf_ctDNA_positive_ch_negative.median_survival_time_)
    median_ctDNA_negative = round(kmf_ctDNA_negative.median_survival_time_)
    
    ax.text(1, 0.9, "Median", transform=ax.transAxes, fontsize=8, ha='center', va='center')
    ax.text(1, 0.83, "months", transform=ax.transAxes, fontsize=8, ha='center', va='center', fontstyle = "italic")
    
    ax.text(0.7, 0.76, f"ctDNA+ CH+", transform=ax.transAxes, fontsize=8, ha='left', va='center', color = "mediumpurple")
    ax.text(0.7, 0.69, f"ctDNA+ CH-", transform=ax.transAxes, fontsize=8, ha='left', va='center', color = "forestgreen")
    ax.text(0.7, 0.62, f"ctDNA-", transform=ax.transAxes, fontsize=8, ha='left', va='center', color = "black")
    
    ax.text(1, 0.76, median_both_positive, transform=ax.transAxes, fontsize=8, ha='center', va='center')
    ax.text(1, 0.69, median_ctDNA_positive_ch_negative, transform=ax.transAxes, fontsize=8, ha='center', va='center')
    ax.text(1, 0.62, median_ctDNA_negative, transform=ax.transAxes, fontsize=8, ha='center', va='center')
   
    return(ax)

# LOAD CHIP DATASETS
all_vars_chip = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/chip.csv"))
all_vars_chip = all_vars_chip[all_vars_chip["Dependent"] == False]
base_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
prog_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)
prog_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)

# LOAD SOMATIC DATASETS
all_vars_somatic = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/somatic.csv"))
all_vars_chip = all_vars_chip[all_vars_chip["Dependent"] == False]
all_vars_somatic = all_vars_somatic[~all_vars_somatic["Patient_id"].isin(["20-313", "21-184", "21-430"])] # exclude some samples due to oxidative damage
base_kidney_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Kidney")].reset_index(drop = True)
prog_kidney_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "During treatment") & (all_vars_somatic["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Bladder")].reset_index(drop = True)
prog_bladder_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "During treatment") & (all_vars_somatic["Diagnosis"] == "Bladder")].reset_index(drop = True)

all_pts = pd.read_csv(PATH_sample_information, sep="\t", names=["Patient_id", "Date", "Diagnosis", "Timepoint"])
kidney_pts = all_pts[all_pts["Diagnosis"] == "Kidney"]["Patient_id"].unique()
bladder_pts = all_pts[all_pts["Diagnosis"] == "Bladder"]["Patient_id"].unique()

fig = plt.figure(figsize=(2.5, 11))
gs = gridspec.GridSpec(1, 2, width_ratios = [1, 1], hspace = 0, wspace = 0) # outer most gs with 3 rows
kidney_ch_ax_combined = plt.subplot(gs[0])
bladder_ch_ax_combined = plt.subplot(gs[1])

