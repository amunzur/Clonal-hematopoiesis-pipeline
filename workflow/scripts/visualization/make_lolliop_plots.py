"""
Given a gene name this script plots the CH mutations called in that gene. 
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

mpl.rcParams['font.size'] = 8
mpl.rcParams['text.color'] = 'k'
mpl.rcParams['legend.fontsize'] = 8
mpl.rcParams['legend.handletextpad'] = '0.8'
mpl.rcParams['legend.labelspacing'] = '0.4'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['axes.linewidth'] = 0.6
mpl.rcParams['xtick.major.width'] = 0.4  # Set the linewidth of major ticks
mpl.rcParams['ytick.major.width'] = 0.4  # Set the linewidth of major ticks
mpl.rcParams['xtick.minor.size'] = 0.6
mpl.rcParams['ytick.minor.size'] = 0.6
mpl.rcParams['xtick.labelsize'] = 6
mpl.rcParams['ytick.labelsize'] = 6
mpl.rcParams['axes.labelsize'] = 6


def plot_domains(df, domain_dict, domain_size, domain_color_dict, ax):
    """
    Given a dict of domain names and start and end positions of the domain, plot a rectangle on the ax to represent it. 
    """
    ax.axhline(y=0.5, color='black', linewidth=0.5, zorder = -1)
    for key, value in domain_dict.items():    
        domain_name = key
        domain_start = value[0]    
        domain_end = value[1]
        rect_height = 0.3  # Height of the rectangles
        rect_width = domain_end - domain_start
        rect_y = 0.5 - rect_height / 2  # Position the rectangles such that the horizontal line is in the middle
        ax.add_patch(plt.Rectangle((domain_start, rect_y), rect_width, rect_height, color=domain_color_dict[domain_name]))
        annotation_x = domain_start + rect_width / 2
        ax.annotate(domain_name, (annotation_x, -0.4), ha='center', va='center', fontsize = 8)
    ax.set_xlim((1, domain_size))
    ax.set_xticks(range(1, domain_size + 1, domain_size // 4))  
    return(ax)

def plot_lollipops(df, ax, mut_color_dict):
    """
    Plot the lollipops on the domain rectangles.
    """
    length_df = df["Protein_annotation"].value_counts().reset_index().rename(columns = {"index": "Protein_annotation", "Protein_annotation": "lollipop_length"})
    df = df.merge(length_df, how = "left").drop_duplicates("Protein_annotation")
    for i, group in df.groupby("Protein_annotation"):
        mut_type = group["Consequence"].values[0]
        mut_aa = int(group["Mut aa"].values[0])
        mut_length = int(group["lollipop_length"].values[0]) + 1
        ax.plot([mut_aa, mut_aa], [0.5, mut_length], color="dimgray", linewidth=0.4, zorder = -2)
        ax.scatter(mut_aa, mut_length, s=10, edgecolor = "white", linewidth = 0.2, color = mut_color_dict[mut_type])
        ax.set_xlabel("")
    ax.set_ylim((-1, max(df["lollipop_length"]) + 2))
    ax.set_yticks(list(range(2, max(df["lollipop_length"]) + 2)))
    yticks = len(ax.get_yticks())
    # max_length = max(df["lollipop_length"])
    # ytick_labels = [1] + list(np.repeat("", len(ax.get_yticks()) - 2)) + [max_length]
    # ax.set_yticklabels(ytick_labels, fontsize = 5)
    ax.set_yticklabels(list(range(1, max(df["lollipop_length"])+1)))
    ax.set_ylabel("n mutations", rotation = 90, fontsize = 10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(axis='x', pad=1)
    return(ax)


PATH_muts_list = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/chip_SSCS2_curated_complete.csv"
DIR_working = "/groups/wyattgrp/users/amunzur/pipeline"
mut_color_dict = {
    "Frameshift indel": '#FFC907',
    "Nonframeshift indel": '#a9a9a9',
    "Missense": '#79B443',
    "Stopgain": '#BD4398',
    "Splicing": "darkorange"}

gene_names = ["DNMT3A", "PPM1D", "TP53"]
gene_size_aa = [912, 605, 393]
xtickdict = {"DNMT3A":[1, 300, 600, 912], "PPM1D":[1, 200, 400, 605], "TP53":[1, 100, 200, 300, 393]}
xticklabeldict = {"DNMT3A":["1", "300", "600", "912aa"], "PPM1D":["1", "200", "400", "605aa"], "TP53":["1", "100", "200", "300", "393aa"]}

domains_list = [{"PWWP": [291, 374], "DNA methylase": [634, 767]}, {"PP2C": [67, 368]}, {"P53 TAD": [6, 29], "P53": [95, 288], "P53 tetramer": [318, 358]}]

# LOAD CHIP DATASETS
all_vars_chip = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/chip.csv"))
all_vars_chip["Consequence"] = all_vars_chip["Consequence"].replace({"Frameshift deletion": "Frameshift indel", "Frameshift Deletion": "Frameshift indel", "Frameshift insertion": "Frameshift indel", "Nonframeshift deletion": "Nonframeshift indel", "Nonframeshift insertion": "Nonframeshift indel", "Startloss": "Missense"})
baseline_df = all_vars_chip[all_vars_chip["Timepoint"] == "Baseline"].reset_index(drop = True)
base_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
prog_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)
prog_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)

domain_color_dict = {"PWWP": "limegreen", "DNA methylase":"aquamarine", "PP2C":"violet", "P53 TAD":"lightcoral", "P53":"cornflowerblue", "P53 tetramer":"orange"}

for gene, domain_size, domain_dict in zip(gene_names, gene_size_aa, domains_list): 
    xticklist = xtickdict[gene]
    xticklabel = xticklabeldict[gene]
    df = baseline_df[(baseline_df["Gene"] == gene) & (baseline_df["Protein_annotation"].str.startswith("p."))]
    df["Mut aa"] = df["Protein_annotation"].str.extract(r'(\d+)')
    fig, ax = plt.subplots(figsize = (5.5, 1.2))
    ax = plot_domains(df, domain_dict, domain_size, domain_color_dict, ax) # plot the outline of the gene as distinct genes
    ax = plot_lollipops(df, ax, mut_color_dict)
    ax.set_xticks(xticklist)
    ax.set_xticklabels(xticklabel, fontsize = 8)
    ax.tick_params(axis='x', which='both', bottom=True, length = 2)
    ax.tick_params(axis='y', which='both', length = 2)
    ax.set_title(gene, loc = "left")
    if gene == "DNMT3A":
        ax.annotate("R882", (882, 8.6), ha='center', va='center', fontsize = 6)
        ax.annotate("Y735C", (735, 7.6), ha='center', va='center', fontsize = 6)
        ax.annotate("L637Q", (637, 6.6), ha='center', va='center', fontsize = 6)        
        ax.annotate("R771X", (771, 6.6), ha='center', va='center', fontsize = 6)
    if gene == "TP53":
        ax.annotate("R273H", (273, 4.5), ha='center', va='center', fontsize = 6)
        ax.annotate("M237I", (237, 4.5), ha='center', va='center', fontsize = 6)
        ax.annotate("R205C", (205, 3.5), ha='center', va='center', fontsize = 6)
    if gene == "PPM1D":
        ax.annotate("C478X", (478, 6.5), ha='center', va='center', fontsize = 6)
        ax.annotate("R572X", (572, 5.6), ha='center', va='center', fontsize = 6)
    fig.tight_layout()
    fig.savefig(f"/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/{gene}.pdf")

