"""
Plots lollipops.
"""

import pandas as pd
import numpy as np
import os 
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

def plot_domains(domain_dict, gene_size, domain_color_dict, ax):
    """
    Given a dict of domain names and start and end positions of the domain, plot a rectangle on the ax to represent it. 
    """
    ax.axhline(y=0, color='black', linewidth=0.5, zorder = -1)
    for key, value in domain_dict.items():    
        domain_name = key
        domain_start = value[0]    
        domain_end = value[1]
        rect_height = 0.3  # Height of the rectangles
        rect_width = domain_end - domain_start
        rect_y = 0 - rect_height / 2  # Position the rectangles such that the horizontal line is in the middle
        ax.add_patch(plt.Rectangle((domain_start, rect_y), rect_width, rect_height, color=domain_color_dict[domain_name]))
        annotation_x = domain_start + rect_width / 2
        if key == "HARE-HTH":
            annotation_x += rect_width /2
        ax.annotate(domain_name, (annotation_x, -0.4), ha='center', va='center', fontsize = 8)
    ax.set_xlim((1, gene_size))
    # ax.set_xticks(range(1, domain_size + 1, domain_size // 4))  
    return(ax)

def plot_lollipops(df, ax, mut_color_dict, plot_what):
    """
    Plot the lollipops on the domain rectangles.
    """
    if plot_what == "chip": 
        factor = 1
    else:
        factor = -1
    length_df = df["Protein_annotation"].value_counts().reset_index().rename(columns = {"index": "Protein_annotation", "Protein_annotation": "lollipop_length"})
    df = df.merge(length_df, how = "left").drop_duplicates("Protein_annotation")
    mut_length_df = []
    for i, group in df.groupby("Protein_annotation"):
        mut_type = group["Consequence"].values[0]
        mut_aa = int(group["Mut aa"].values[0])
        mut_length = int(group["lollipop_length"].values[0])
        mut_length_df.append(mut_length)
        ax.plot([mut_aa, mut_aa], [0, mut_length*factor], color="dimgray", linewidth=0.5, zorder = -2)
        ax.scatter(mut_aa, mut_length*factor, s=20, edgecolor = "white", linewidth=0.2, color = mut_color_dict[mut_type])
        ax.set_xlabel("")
    ax.set_ylabel("n mutations", rotation = 90, fontsize = 10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(axis='x', pad=1)
    return(ax, max(mut_length_df))

PATH_muts_list = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/chip.csv"
PATH_ctdna = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/somatic.csv"
PATH_sample_information = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"
sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])[["Patient_id", "Diagnosis"]].drop_duplicates().reset_index(drop = True)

# LOAD CHIP DATASETS
all_vars_chip = pd.read_csv(PATH_muts_list)
all_vars_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Dependent"] == False)].reset_index(drop = True)
all_vars_chip = all_vars_chip.merge(sample_info, how = "inner")
all_vars_chip["Consequence"] = all_vars_chip["Consequence"].replace({"Frameshift deletion": "Frameshift indel", "Frameshift Deletion": "Frameshift indel", "Frameshift insertion": "Frameshift indel", "Nonframeshift deletion": "Nonframeshift indel", "Nonframeshift insertion": "Nonframeshift indel", "Startloss": "Missense"})

ctdna_muts = pd.read_csv(PATH_ctdna)
ctdna_muts["Consequence"] = ctdna_muts["Consequence"].replace({"Frameshift deletion": "Frameshift indel", "Frameshift Deletion": "Frameshift indel", "Frameshift insertion": "Frameshift indel", "Nonframeshift deletion": "Nonframeshift indel", "Nonframeshift insertion": "Nonframeshift indel", "Startloss": "Missense"})
ctdna_muts = ctdna_muts[(ctdna_muts["Timepoint"] == "Baseline") & (ctdna_muts["Dependent"] == False)].reset_index(drop = True)

mut_color_dict = {
    "Frameshift indel": '#FFC907',
    "Nonframeshift indel": '#a9a9a9',
    "Missense": '#79B443',
    "Stopgain": '#BD4398'}

gene_names = ["DNMT3A", "TET2", "ASXL1", "PPM1D", "TP53"]
gene_size_aa_dict = {"DNMT3A": 912, "TET2": 2002, "ASXL1": 1541, "PPM1D": 605, "TP53": 393}
xtickdict = {"DNMT3A":[1, 300, 600, 912], 
             "TET2": [1, 400, 800, 1200, 1600, 2002],
             "ASXL1": [1, 400, 800, 1200, 1541],
             "PPM1D":[1, 200, 400, 605], 
             "TP53":[1, 100, 200, 300, 393]}

xticklabeldict = {"DNMT3A":["1", "300", "600", "912aa"], 
                  "TET2": ["1", "400", "800", "1200", "1600", "2002aa"],
                  "ASXL1": ["1", "400", "800", "1200", "1541aa"],
                  "PPM1D":["1", "200", "400", "605aa"], 
                  "TP53":["1", "100", "200", "300", "393aa"]}

domains_dict = {"DNMT3A": {"PWWP": [291, 374], "DNA methylase": [634, 767]}, 
                "TET2": {"Tet_JBP": [1290, 1905]},
                "ASXL1": {"HARE-HTH": [12, 83], "ASHX": [234, 362], "PHD_3": [1480, 1539]},
                "PPM1D": {"PP2C": [67, 368]}, 
                "TP53": {"P53 TAD": [6, 29], "P53": [95, 288], "P53 tetramer": [318, 358]}}

domain_color_dict = {"PWWP": "limegreen", 
                     "DNA methylase":"aquamarine", 
                     "PP2C":"violet", 
                     "P53 TAD":"lightcoral", 
                     "P53":"cornflowerblue", 
                     "P53 tetramer":"orange", 
                     "Tet_JBP":"limegreen", 
                     "HARE-HTH":"limegreen",
                     "ASHX":"lightcoral",
                     "PHD_3":"cornflowerblue"}

fig = plt.figure(figsize=(8, 8))
gs_outer = gridspec.GridSpec(2, 1, height_ratios=[2, 0.2], hspace = 0.1, wspace = 0)

gs_genes = gridspec.GridSpecFromSubplotSpec(5, 1, height_ratios=[1,1,1,1,1], subplot_spec=gs_outer[0], wspace=0.5, hspace = 0.5)
gs_legend = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios=[0.1, 1], subplot_spec=gs_outer[1], wspace=0.35, hspace = 0.3)

fontsize_protein_annots = 9
labels = ['A', 'B', 'C', 'D', 'E']

for i, gene in enumerate(gene_names):
    ylims_df = {}
    gene_size = gene_size_aa_dict[gene]
    gene_domains_dict = domains_dict[gene]
    xticklist = xtickdict[gene]
    xticklabel = xticklabeldict[gene]
    
    # Set up domains
    ax = plt.subplot(gs_genes[i])
    ax = plot_domains(gene_domains_dict, gene_size, domain_color_dict, ax)
    
    # Plot the lollipops
    max_lollipop_length_MAX = 0
    for plot_what, mut_df in zip(["chip", "ctdna"], [all_vars_chip, ctdna_muts]):
        df = mut_df[(mut_df["Gene"] == gene) & (mut_df["Protein_annotation"].str.startswith("p."))]
        df["Mut aa"] = df["Protein_annotation"].str.extract(r'(\d+)')
        
        ax, max_lollipop_length = plot_lollipops(df, ax, mut_color_dict, plot_what)
        if max_lollipop_length > max_lollipop_length_MAX:
            max_lollipop_length_MAX = max_lollipop_length
        
        ax.set_xticks(xticklist)
        ax.set_xticklabels(xticklabel, fontsize = 9)
        ax.tick_params(axis='x', which='both', bottom=True, length = 2)
        ax.tick_params(axis='y', which='both', length = 2)
        ax.set_title(gene, loc = "left", fontstyle = "italic")
        
        ylims_df[plot_what] = max_lollipop_length
        
        # Add label to upper left corner
        ax.text(-0.02, 1.3, labels[i], transform=ax.transAxes, fontsize=25, fontweight='bold', va='top', ha='right')
        
        if gene == "DNMT3A":
            ax.annotate("R882", (882, 5.5), ha='center', va='center', fontsize = fontsize_protein_annots)
            ax.annotate("R729W", (729, 5.5), ha='center', va='center', fontsize = fontsize_protein_annots)
            ax.annotate("L637Q", (600, 4.5), ha='center', va='center', fontsize = fontsize_protein_annots)
            ax.annotate("Y735C", (770, 4), ha='center', va='center', fontsize = fontsize_protein_annots)
        if gene == "TP53":
            ax.annotate("Y220C", (220, 1.7), ha='center', va='center', fontsize = fontsize_protein_annots)
            ax.annotate("R175H", (175, -4.7), ha='center', va='center', fontsize = fontsize_protein_annots)
            ax.annotate("R248Q", (248, -3.7), ha='center', va='center', fontsize = fontsize_protein_annots)
        if gene == "PPM1D":
            ax.annotate("C478X", (478, 5.7), ha='center', va='center', fontsize = fontsize_protein_annots)
            ax.annotate("R572X", (572, 4.8), ha='center', va='center', fontsize = fontsize_protein_annots)
        if gene == "ASXL1": 
            ax.annotate("Q546X", (533, 3), ha='right', va='center', fontsize = fontsize_protein_annots)
            ax.annotate("E635Rfs*15", (651, 3), ha='left', va='center', fontsize = fontsize_protein_annots)
        if gene == "TET2":
            ax.annotate("I1873T", (1873, 3.2), ha='center', va='center', fontsize = fontsize_protein_annots)
            ax.annotate("N275Ifs*18", (275, 2.2), ha='center', va='center', fontsize = fontsize_protein_annots)
    
    ax.set_ylim((-ylims_df["ctdna"]-1, ylims_df["chip"]+1))
    ax.set_yticks(list(range(-ylims_df["ctdna"]-1, ylims_df["chip"]+1)))
    ax.set_yticklabels([abs(x) for x in range(-ylims_df["ctdna"]-1, ylims_df["chip"]+1)])

# add legend
legend_ax = plt.subplot(gs_legend[0])
legend_colors = mut_color_dict.values()
legend_labels = mut_color_dict.keys()
legend_handles = [plt.Line2D([0], [0], marker='o', color=color, label=label, markersize=5, linestyle='') for color, label in zip(legend_colors, legend_labels)]
legend_ax.legend(handles=legend_handles, loc="upper left", frameon=False, fontsize = 10, handletextpad=0.3, ncol=len(legend_handles))
legend_ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
legend_ax.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labelleft=False)
legend_ax.set_facecolor('none')

gs_outer.tight_layout(fig)
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/SUPP_lollipops.png")
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/SUPP_lollipops.pdf")