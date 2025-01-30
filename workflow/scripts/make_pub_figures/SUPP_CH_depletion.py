"""
Shows the prevalence of CH in the cohort across various minimum detection thresholds.
"""

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
from sklearn.linear_model import LinearRegression


def plot_ch_and_ctDNA_coocurrence(chip_muts_df, ctDNA_muts_df, ax, diagnosis):
    """
    Barchart, Fishers on co-occurence of CH and ctDNA
    """
    all_pts = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])    
    all_pts = all_pts[(all_pts["Timepoint"] == "Baseline") & (all_pts["Diagnosis"] == diagnosis)]["Patient_id"].tolist()
    n_pts_total = len(all_pts)
    
    pts_with_ch = chip_muts_df["Patient_id"].unique()
    pts_with_ctDNA = ctDNA_muts_df["Patient_id"].unique()
    pts_without_ch = [p for p in all_pts if p not in pts_with_ch]
    pts_without_ctDNA = [p for p in all_pts if p not in pts_with_ctDNA]
    
    pts_with_ch_with_ctDNA = list(set(pts_with_ch) & set(pts_with_ctDNA))
    pts_with_ch_without_ctDNA = list(set(pts_with_ch) & set(pts_without_ctDNA))
    
    pts_without_ch_with_ctDNA = list(set(pts_without_ch) & set(pts_with_ctDNA))
    pts_without_ch_without_ctDNA = list(set(pts_without_ch) & set(pts_without_ctDNA))
    
    n_pts_with_ch_with_ctDNA = len(pts_with_ch_with_ctDNA)
    n_pts_with_ch_without_ctDNA = len(pts_with_ch_without_ctDNA)
    
    n_pts_without_ch_with_ctDNA = len(pts_without_ch_with_ctDNA)
    n_pts_without_ch_without_ctDNA = len(pts_without_ch_without_ctDNA)
    
    # Plotting
    ctDNA_absent_color = (34/255, 139/255, 34/255, 120/255)
    
    ax.bar(0, n_pts_with_ch_without_ctDNA, color = ctDNA_absent_color, edgecolor = "none", width = 0.5)
    ax.bar(0, n_pts_with_ch_with_ctDNA, bottom = n_pts_with_ch_without_ctDNA, color = "forestgreen", edgecolor = "none", width = 0.5)
    
    ax.bar(1, n_pts_without_ch_without_ctDNA, color = ctDNA_absent_color, edgecolor = "none", width = 0.5)
    ax.bar(1, n_pts_without_ch_with_ctDNA, bottom = n_pts_without_ch_without_ctDNA, color = "forestgreen", edgecolor = "none", width = 0.5)
    
    # Annotate n patiens on bars
    ax.text(0, n_pts_with_ch_without_ctDNA, str(n_pts_with_ch_without_ctDNA), ha='center', va='top', color='black', fontsize = 8, zorder = 15)
    ax.text(0, n_pts_with_ch_with_ctDNA+n_pts_with_ch_without_ctDNA-2, str(n_pts_with_ch_with_ctDNA), ha='center', va='top', color='black', fontsize = 8, zorder = 15)
    
    ax.text(1, n_pts_without_ch_without_ctDNA, str(n_pts_without_ch_without_ctDNA), ha='center', va='top', color='black', fontsize = 8, zorder = 15)
    ax.text(1, n_pts_without_ch_with_ctDNA+n_pts_without_ch_without_ctDNA-2, str(n_pts_without_ch_with_ctDNA), ha='center', va='top', color='black', fontsize = 8, zorder = 15)
    
    # Run fishers
    oddsratio, p_value = fisher_exact([[n_pts_with_ch_without_ctDNA, n_pts_with_ch_with_ctDNA], [n_pts_without_ch_without_ctDNA, n_pts_without_ch_with_ctDNA]])
    p_val = round(p_value, 3)
    ax.text(0.5, 1, f"P={p_val}", ha='center', va='center', color='black', fontsize = 8, zorder = 15, transform=ax.transAxes)
    
    ax.set_xlim((-0.5, 1.5))
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["CH+", "CH-"])
    ax.spines[["top", "right"]].set_visible(False)
    
    legend_colors = ["forestgreen", ctDNA_absent_color]
    legend_labels = ["ctDNA+", "ctDNA-"]
    legend_handles = [plt.Line2D([0], [0], marker='s', color=color, label=label, markersize=3, linestyle='') for color, label in zip(legend_colors, legend_labels)]
    ax.legend(handles=legend_handles, loc="best", frameon=False, fontsize = 8, handletextpad=0.1)
    
    if diagnosis == "Bladder":
        ax.set_title("mUC")
        ax.set_yticks([0, 25, 50, 75, 100])
        ax.set_ylabel("Number of patients")
    else:
        ax.set_title("mRCC")
        ax.set_yticks([0, 50, 100, 150])
        
    return(ax)

def fit_linear_reg(x_vals, y_vals, ax):
    """
    Fits linear reg. 
    """
    x = np.array(x_vals).reshape(-1, 1)
    y = np.array(y_vals).reshape(-1, 1)
    
    model = LinearRegression()
    model.fit(x, y)
    
    slope = model.coef_[0][0]
    intercept = model.intercept_[0]
    
    ax.plot(x, model.predict(x), color='red', linewidth=0.5)
    ax.text(0.25, 0.1, f'y={slope:.2f}x+{intercept:.2f}', transform=ax.transAxes, fontsize=8, verticalalignment='top', color='red')
    
    # Now we add the y=x line
    x = np.linspace(0, 100, 1000)  # Generate 100 points between 0 and 10
    y = np.random.random(100) * 100  # Random values for y
    ax.plot(x, x, color='blue', linestyle='--', label='y=x', linewidth=0.5)
    ax.text(1, 0.2, 'y=x', transform=ax.transAxes, fontsize=8, horizontalalignment = "right", verticalalignment='top', color='blue')
    return(ax)

def add_pie(xvals, yvals, ax):
    # Calculate the counts of points above, below, and on the y = x line
    above_line = np.sum(yvals > xvals)  # Points where y > x
    below_line = np.sum(yvals < xvals)  # Points where y < x
    
    # Data for the pie chart
    sizes = [above_line, below_line]
    labels = ['Above y=x', 'Below y=x']
    colors = ['lightcoral', 'lightblue']
    
    # Create the pie chart on the same ax as the scatter plot
    ax = ax.pie(sizes, labels=labels, colors=colors, autopct='%d%%', startangle=90, radius=1.5,
    textprops={'fontsize': 8}, 
    wedgeprops={'edgecolor': 'black', 'linewidth': 0.5, 'facecolor': 'whitesmoke'})
    return(ax)

def plot_CH_depletion(baseline, baseline_somatic, ax3, ax4, ax3_pie, ax4_pie):
    """
    Linear reg between cfDNA vaf and WBC vaf
    """
    # Get the highest ctDNA VAF per pt
    pts_max_ctdna = baseline_somatic.groupby("Patient_id")["VAF_t"].max().reset_index()
    
    all_pts = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])    
    all_pts = all_pts[all_pts["Timepoint"] == "Baseline"]["Patient_id"].tolist()
    
    pts_without_ctDNA = [p for p in all_pts if p not in baseline_somatic["Patient_id"].unique()]
    
    ctneg_df = pd.DataFrame({'Patient_id': pts_without_ctDNA, 'VAF_t': 0})
    pts_max_ctdna = pd.concat([pts_max_ctdna, ctneg_df]).reset_index(drop = True)
    
    # Patients with ct max less than 10
    less_10 = pts_max_ctdna[pts_max_ctdna["VAF_t"] < 10]
    less_10_muts = baseline[baseline["Patient_id"].isin(less_10["Patient_id"])]
    
    # Patients with ct max more than or equal to 10
    more_10 = pts_max_ctdna[pts_max_ctdna["VAF_t"] >= 10]
    more_10_muts = baseline[baseline["Patient_id"].isin(more_10["Patient_id"])]
    
    # Plotting
    ax3.scatter(less_10_muts["VAF_n"], less_10_muts["VAF_t"], s = 3, color = "black", alpha = .6)
    ax4.scatter(more_10_muts["VAF_n"], more_10_muts["VAF_t"], s = 3, color = "black", alpha = .6)
    
    # Fit linear reg
    ax3 = fit_linear_reg(x_vals = less_10_muts["VAF_n"], y_vals = less_10_muts["VAF_t"], ax = ax3)
    ax4 = fit_linear_reg(x_vals = more_10_muts["VAF_n"], y_vals = more_10_muts["VAF_t"], ax = ax4)
    
    # Add the pies for the count of dots above, below and on the y=x line.
    ax3_pie = add_pie(xvals = less_10_muts["VAF_n"], yvals = less_10_muts["VAF_t"], ax=ax3_pie)
    ax4_pie = add_pie(xvals = more_10_muts["VAF_n"], yvals = more_10_muts["VAF_t"], ax=ax4_pie)
    
    # Aes
    ax3.spines[["top", "right"]].set_visible(False)
    ax4.spines[["top", "right"]].set_visible(False)
    ax3.set_ylabel("Plasma cfDNA VAF%")
    ax3.set_xlabel("WBC VAF%")
    ax4.set_xlabel("WBC VAF%")
    
    ax3.set_title("Max ctDNA VAF<10%", fontsize = 10)
    ax4.set_title("Max ctDNA VAF≥10%", fontsize = 10)
    
    ax3.set_ylim((0, 50))
    ax3.set_xlim((0, 50))
    ax3.set_xticks([0, 25, 50])
    ax3.set_yticks([0, 25, 50])
    
    ax4.set_ylim((0, 60))
    ax4.set_xlim((0, 60))
    ax4.set_xticks([0, 30, 60])
    ax4.set_yticks([0, 30, 60])
    
    return(ax3, ax4, ax3_pie, ax4_pie)

def explore_ch_only_patients(baseline, baseline_somatic, ax):
    """
    Explores patients with CH mutations only. Makes a plot.
    """
    sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])    
    sample_info = sample_info[["Patient_id", "Diagnosis"]].drop_duplicates()
    
    pts_with_ch = baseline["Patient_id"].unique()
    pts_with_ctdna = baseline_somatic["Patient_id"].unique()
    
    pts_with_ch_only = [p for p in pts_with_ch if p not in pts_with_ctdna]
    
    pts_with_ch_only_df = pd.DataFrame(pts_with_ch_only, columns=["Patient_id"])    
    merged_df = pts_with_ch_only_df.merge(sample_info, on="Patient_id", how="left")
    
    ch_only_df = baseline.merge(merged_df, how = "inner")[["Patient_id", "Diagnosis", "Gene", "VAF_t"]]
    
    other_genes = ['EZH2', 'GNB1', 'JAK2', 'BCOR', 'GNAS', 'SH2B3', 'CUX1', 'SRSF2', 'RAD21', 'CBL', 'BRCC3', 'GATA2', 'SETD2', 'IDH2', 'MYD88', 'SF3B1', 'RUNX1']
    ch_only_df.loc[ch_only_df["Gene"].isin(other_genes), "Gene"] = "Other"
    patient_df = ch_only_df[["Patient_id", "Diagnosis", "Gene"]].drop_duplicates().reset_index(drop = True)
    
    # Convert to log scale
    ch_only_df["VAF_t"] = np.log10(ch_only_df["VAF_t"])
    
    # Plotting. Barchart and the swarmplot superimposed to show
    diagnosis_bar_offset_value_dict = {"Bladder": -0.2, "Kidney": 0.2}
    barcolor_dict = {"Bladder": "deepskyblue", "Kidney": "orangered"}
    
    ax2 = ax.twinx()
    gene_order = ['DNMT3A', 'TET2', 'ASXL1', 'TP53', 'PPM1D', 'CHEK2', 'BRCA1', 'ATM', 'KDM6A', 'TERT', 'STAG2', 'PBRM1', 'ERBB2', 'MTOR', 'KMT2D', "Other"]
    for i, gene in enumerate(gene_order):
        for diagnosis in ["Bladder", "Kidney"]:
            patient_df_gene = patient_df[(patient_df["Diagnosis"] == diagnosis) & (patient_df["Gene"] == gene)]
            offset_val = diagnosis_bar_offset_value_dict[diagnosis]
            color = barcolor_dict[diagnosis]
            ax.bar(i+offset_val, patient_df_gene.shape[0], color=color, width = 0.4, edgecolor = "None")
            
            # Swarm plot
            vaf_df = ch_only_df[(ch_only_df["Gene"] == gene) & (ch_only_df["Diagnosis"] == diagnosis)]["VAF_t"]
            for vaf in vaf_df:
                jitter = np.random.uniform(-0.08, 0.08, 1)
                ax2.scatter(i+offset_val+jitter[0], vaf, color="black", s = 0.15, alpha = 0.7)
    
    # AES
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylabel("Number of patients")
    ax.set_xlim((-1, i+1))
    ax.set_xticks(np.arange(0, i+1))
    ax.set_xticklabels(gene_order, fontsize = 8, rotation = 90, fontstyle = "italic")
    ax.set_ylim([0, 45])
    ax.set_yticks([0, 10, 20, 30, 40])
    ax.set_title("Patient with CH mutations only", loc = "left", fontsize = 10)
    
    ax2.spines[["top", "left"]].set_visible(False)
    ax2.set_ylabel("Plasma cfDNA VAF%")
    ax2.set_ylim((np.log10(0.25), np.log10(50)))
    ax2.set_yticks([np.log10(0.25), np.log10(1), np.log10(2), np.log10(10), np.log10(50)])
    ax2.set_yticklabels([".25", "1", "2", "10", "50"])
    
    legend_colors = ["deepskyblue", "orangered"]
    legend_labels = ["mUC", "mRCC"]
    legend_handles = [plt.Line2D([0], [0], marker='s', color=color, label=label, markersize=5, linestyle='') for color, label in zip(legend_colors, legend_labels)]
    ax.legend(handles=legend_handles, loc="best", frameon=False, fontsize = 8, handlelength=2, handletextpad = 0.1)
    
    return(ax, ax2)

mpl.rcParams['font.size'] = 10

DIR_working = "/groups/wyattgrp/users/amunzur/pipeline"
PATH_sample_information = os.path.join(DIR_working, "resources/sample_lists/sample_information.tsv")
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

all_vars_somatic = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/somatic.csv"))
all_vars_somatic = all_vars_somatic[all_vars_somatic["Dependent"] == False]

base_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
prog_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)
prog_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)

base_kidney_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Kidney")].reset_index(drop = True)
prog_kidney_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "During treatment") & (all_vars_somatic["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Bladder")].reset_index(drop = True)
prog_bladder_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "During treatment") & (all_vars_somatic["Diagnosis"] == "Bladder")].reset_index(drop = True)

baseline = pd.concat([base_kidney_chip, base_bladder_chip]).reset_index(drop = True)
baseline_somatic = pd.concat([base_kidney_somatic, base_bladder_somatic]).reset_index(drop = True)

fig = plt.figure(figsize=(5, 6))
gs = gridspec.GridSpec(2, 1, height_ratios=[0.6, 2], hspace = 0.3)
gs_top_row = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=[1, 1], subplot_spec=gs[0], wspace=0.4)
gs_bottom_two_rows = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios=[1, 1], subplot_spec=gs[1], hspace=0.9)
gs_bottom_row = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=[1, 1], subplot_spec=gs_bottom_two_rows[1], wspace=0.5)
gs_bottom_row_left = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=[1, 0.5], subplot_spec=gs_bottom_row[0], wspace=0)
gs_bottom_row_right = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=[1, 0.5], subplot_spec=gs_bottom_row[1], wspace=0)

ax1 = plt.subplot(gs_top_row[0])
ax2 = plt.subplot(gs_top_row[1])
ax3 = plt.subplot(gs_bottom_row_left[0])
ax4 = plt.subplot(gs_bottom_row_right[0])
ax3_pie = plt.subplot(gs_bottom_row_left[1])
ax4_pie = plt.subplot(gs_bottom_row_right[1])
ax_middle = plt.subplot(gs_bottom_two_rows[0])

ax3.set_aspect('equal')
ax4.set_aspect('equal')

ax1 = plot_ch_and_ctDNA_coocurrence(chip_muts_df = base_bladder_chip, ctDNA_muts_df= base_bladder_somatic, ax = ax1, diagnosis = "Bladder")
ax2 = plot_ch_and_ctDNA_coocurrence(chip_muts_df = base_kidney_chip, ctDNA_muts_df = base_kidney_somatic, ax = ax2, diagnosis = "Kidney")

ax3, ax4, ax3_pie, ax4_pie = plot_CH_depletion(baseline, baseline_somatic, ax3, ax4, ax3_pie, ax4_pie)

ax_middle, ax_swarm = explore_ch_only_patients(baseline, baseline_somatic, ax = ax_middle)

fig.text(0.02, 0.95, 'a', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
fig.text(0.02, 0.67, 'b', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)
fig.text(0.02, 0.36, 'c', size=20, weight='bold', ha='left', va='top', transform=fig.transFigure)

fig.tight_layout()
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/SUPP_depletion.png")
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/SUPP_depletion.pdf")