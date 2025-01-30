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
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.cm import get_cmap
from datetime import datetime
import matplotlib.cm as cm


def return_chip_df(depth_name, min_alt_reads):
    """
    Given a depth name, return the relevant chip output file.
    """
    path_chip = os.path.join("/groups/wyattgrp/users/amunzur/pipeline/results/wbc_downsampling", depth_name, f"variant_calling/min_alt_reads_{min_alt_reads}_chip.csv")
    print(f"Loading {path_chip}")
    df = pd.read_csv(path_chip)
    df = df.rename(columns = {"Sample_name": "Sample_name_n"})
    df.loc[df["Gene"] == "U2AF1\\x3bU2AF1L5", "Gene"] = "U2AF1"
    return(df)

def compare_to_groundtruth(chip, df_depth):
    """
    Determine how many of the mutations could be detected in the df.
    """
    df_depth = df_depth[["Patient_id", "Sample_name_n", 'Date_collected', 'Timepoint', 'Diagnosis', 
        'Chrom', 'Position', 'Ref', 'Alt', 'Gene', 'Protein_annotation']]
    merged = chip.merge(df_depth, how = "left", indicator = True)
    merged_df = merged["_merge"].value_counts(normalize = True).reset_index()
    perc_value = merged_df[merged_df["index"] == "both"].reset_index()["_merge"][0]
    return(perc_value)

def run_comparison(chip_df):
    results_dict_1 = {}
    for depth in [1500, 1400, 1300, 1200, 1100, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100]:
        print(f"Working on depth {depth}")
        results_dict_2 = {}
        for min_reads in [2, 3, 4, 5]:
            df_depth = return_chip_df(f"depth_{depth}", min_reads)  # Assuming this function returns a dataframe
            perc_value = compare_to_groundtruth(chip_df, df_depth)  # Your comparison function
            results_dict_2[min_reads] = perc_value
        results_dict_1[depth] = results_dict_2
    return results_dict_1


def plotting_wbc_downsampling(results_dict_1, ax, title, set_ylabel = True):
    """
    Plots the percentage of mutations caught for each depth value.
    """
    df = pd.DataFrame(results_dict_1)
    df_long = df.reset_index().melt(id_vars='index', var_name='Depth', value_name='Value')
    
    colors_dict = {5: "#264653", 4: "#2a9d8f", 3: "#e9c46a", 2: "#f4a261"}
    df_long["color"] = df_long["index"].map(colors_dict)
    
    ax.scatter(df_long["Depth"], df_long["Value"], color=df_long["color"], zorder = 100, s = 30)
    
    for i, group in df_long.groupby("index"):
        color = group["color"].unique()[0]
        ax.plot(group["Depth"], group["Value"], linestyle='--', color = color, linewidth = 1)
        
    ax.invert_xaxis()
    ax.set_xlim(1550, 50)
    ax.set_xticks(range(1500, 0, -100))
    ticklabels = [f"{tick}X" for tick in ax.get_xticks()]
    ax.set_xticklabels(ticklabels, rotation = 45, ha = "right")
    ax.set_xlabel("Downsampled WBC sequencing depth")
    
    ax.set_ylim((0, 1.05))
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax.set_yticklabels(["0%", "20%", "40%", "60%", "80%", "100%"])
    if set_ylabel:
        ax.set_ylabel("% of CH calls identified")
    
    # add horizontal lines through ticks
    for tick in ax.get_yticks():
        ax.hlines(y=tick, xmin=0.0, xmax=1.0, color='b')
    
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(title, fontsize = 10)
    return(ax)

def plotting_wbc_downsampling_make_legend(ax_wbc_legend):
    """
    Makes legend.
    """
    colors_dict = {"≥5 reads": "#264653", "≥4 reads": "#2a9d8f", "≥3 reads": "#e9c46a", "≥2 reads": "#f4a261"}
    legend_colors = colors_dict.values()
    legend_labels = colors_dict.keys()
    legend_handles = [plt.Line2D([0], [0], marker='o', color=color, label=label, markersize=5, linestyle='') for color, label in zip(legend_colors, legend_labels)]
    ax_wbc_legend.legend(handles=legend_handles, loc="center", frameon=False, fontsize = 8, handletextpad=0.3, ncol = len(legend_colors))
    ax_wbc_legend.spines[["top", "right", "bottom", "left"]].set_visible(False)
    ax_wbc_legend.set_yticks([])
    ax_wbc_legend.set_xticks([])
    return(ax_wbc_legend)




ax_wbc_legend = plotting_wbc_downsampling_make_legend(ax_wbc_legend)

path_chip = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/chip.csv"

chip = pd.read_csv(path_chip)
chip = chip[chip["Timepoint"] == "Baseline"]
chip = chip[chip["Dependent"] == False]

chip = chip[["Patient_id", "Sample_name_n", 'Date_collected', 'Timepoint', 'Diagnosis', 'Chrom', 'Position', 'Ref', 'Alt', 'Gene', 'AAchange', 'VAF_n']]
chip_1_perc = chip[chip["VAF_n"] > 1]

results_dict_chip = run_comparison(chip) # Run comparison for chip
results_dict_chip_1_perc = run_comparison(chip_1_perc) # Run comparison for chip_1_perc

# Plotting
fig = plt.figure(figsize=(8, 3)) 
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1], hspace=0.4)
ax_all_ch = plt.subplot(gs[0])
ax_ch_1 = plt.subplot(gs[1])

ax_all_ch = plotting(results_dict_1, ax_all_ch, title = "All CH variants")
ax_ch_1 = plotting(results_dict_chip_1_perc, ax_ch_1, title = "CH variants with VAF > 1%", set_ylabel = False)

fig.text(0.5, 0.01, 'WBC sequencing depth', ha='center')
gs.tight_layout(fig)
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/wbc_downsampling.png")