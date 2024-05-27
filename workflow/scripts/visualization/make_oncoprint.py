import os
import pandas as pd
import numpy as np
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
from scipy import stats
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.font_manager import FontProperties
from matplotlib.colors import Normalize, to_rgba
from matplotlib.patches import Rectangle
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import seaborn as sns

def remove_frame(ax, list_to_remove):
    for i in list_to_remove:
        ax.spines[i].set_visible(False)

def enumerate_col(df, old_col, new_col):
    df[new_col] = df[old_col].map({gene: number + 1 for number, gene in enumerate(df[old_col].unique())})
    return(df)

def plot_mut_counts_per_gene(mut_counts_df, vaf_df, ax_mut_counts):
    ax_secondary = ax_mut_counts.twiny()
    values = mut_counts_df.values.T
    bottom = np.zeros_like(mut_counts_df.index)
    
    for i, consequence in enumerate(mut_counts_df.columns):
        color = mut_dict[consequence]
        ax_mut_counts.barh(mut_counts_df.index, values[i], left=bottom, color=color, label=consequence, edgecolor="none")
        bottom += values[i]
    
    remove_frame(ax_mut_counts, ["top", "right"])
    ax_mut_counts.set_xlim((0, 100))
    ax_mut_counts.set_xticks([0, 50, 100])
    ax_mut_counts.set_xticklabels(["0", "50", "100"])
    plt.setp(ax_mut_counts.get_yticklabels(), visible=False)
    ax_mut_counts.set_xlabel("Mutation counts\nper gene", fontsize = 8)
    ax_mut_counts.tick_params(axis="x", direction="out", bottom=True, top = True, colors='k')
    ax_mut_counts.tick_params(axis='y', which='both', left=False, right = False)
    
    # Plot the vafs of mutations
    vaf_df = vaf_df.reset_index()
    for i, row in vaf_df.iterrows():
        jitter = np.random.uniform(-0.2, 0.2)
        ax_secondary.scatter(row['VAF_n'], row["Order on oncoprint"] + jitter, s = 1, color = "black")
    # max_vaf = vaf_df["VAF_n"].max()
    ax_secondary.set_xlim((0, 50))
    ax_secondary.set_xticks([0, 25, 50])
    ax_secondary.set_xticklabels(["0", "25", "50"])
    ax_secondary.tick_params(axis="x", direction="out", which="both", bottom=True, top = True, colors='k')
    plt.setp(ax_secondary.get_yticklabels(), visible=False)
    ax_secondary.tick_params(axis='y', which='both', left=False, right = False)
    ax_secondary.set_xlabel("WBC VAF(%)", fontsize = 8)
    remove_frame(ax_secondary, ["right"]) 
    return [ax_mut_counts, ax_secondary]

########################################################

def plot_mut_counts_per_patient(ax_mut_counts_patient, pt_counts, bar_width):
    # 3. PLOTTING THE MUTATION COUNTS PER PATIENT, ON TOP OF THE ONCOPRINT
    ax_mut_counts_patient.bar(x = pt_counts["Samples_enumerated"], height = pt_counts["Count"], width=bar_width, capstyle='butt', color = "black", edgecolor = "none", linewidth = 0.25, zorder = 10)
    plt.setp(ax_mut_counts_patient.get_xticklabels(), visible=False)
    ax_mut_counts_patient.tick_params(axis='x', which='both', left=False)
    remove_frame(ax_mut_counts_patient, ["top", "right"])
    ax_mut_counts_patient.set_ylabel("Mutation\ncounts", rotation = 90)
    ax_mut_counts_patient.set_ylim(0, 15)
    ax_mut_counts_patient.set_yticks([0, 5, 10])
    ax_mut_counts_patient.set_yticklabels(["0", "5", "10"])
    ax_mut_counts.tick_params(axis='y', which='both', left=False, right = False)
    return ax_mut_counts_patient
########################################################

def plot_ctdna_fractions(ax_tc, pts_ctdna_NO_nan, pts_ctdna_nan, bar_width):
    # 4. PLOT THE CTDNA FRACTIONS
    ax_tc.bar(x = pts_ctdna_NO_nan["Samples_enumerated"], height = pts_ctdna_NO_nan["ctDNA fraction (%)"], width=bar_width, capstyle='butt', color = "black", edgecolor = "none", linewidth = 0.25, zorder = 10)
    ax_tc.bar(x = pts_ctdna_nan["Samples_enumerated"], height = 100, width=bar_width, capstyle='butt', color = "whitesmoke", edgecolor = "none", linewidth = 0.25, zorder = 10)
    ax_tc.set_ylabel("ctDNA %")
    
    plt.setp(ax_tc.get_xticklabels(), visible=False)
    ax_tc.tick_params(axis='x', which='both', bottom=False)
    ax_tc.tick_params(axis='y', which='both', left=True)
    
    remove_frame(ax_tc, ["top", "right"])
    ax_tc.spines['bottom'].set_visible(True)
    return(ax_tc)

def plot_chip_fractions(ax_chip_fr, chip_fr, bar_width):
    # 5. PLOT THE CHIP FRACTIONS
    ax_chip_fr.bar(x = chip_fr["Samples_enumerated"], height = chip_fr["CHIP_fr"], width=bar_width, capstyle='butt', color = "black", edgecolor = "none", linewidth = 0.25, zorder = 10)
    ax_chip_fr.set_ylabel("CH %")
    plt.setp(ax_chip_fr.get_xticklabels(), visible=False)
    ax_chip_fr.tick_params(axis='x', which='both', bottom=False)
    ax_chip_fr.set_ylim((0, 100))
    ax_chip_fr.set_yticks([0, 50])
    ax_chip_fr.set_yticklabels(["0", "50"])
    ax_chip_fr.tick_params(axis='y', which='both', left=True, right = False)
    remove_frame(ax_chip_fr, ["top", "right", "bottom"])
    return ax_chip_fr

def add_legend(mut_dict, bbox_location, loc):
    # LEGEND
    legend_handles = [Patch(color=color, label=effect) for effect, color in mut_dict.items()]
    legend_handles.append(Line2D([0], [0], marker='D', color='black', markersize=3, label='4 or more mutations')) # Add "> 2 mutations" as a black diamond
    ax_mut_counts.legend(handles=legend_handles, loc=loc, bbox_to_anchor=bbox_location, frameon = False)
    return ax_mut_counts

def plot_oncoprint(ax_main, muts_df, bar_height, bar_width):
    """
    If plot_mutations_as_squares is set to False, muts will be plotted as colors.
    """
    pts = muts_df["Samples_enumerated"].unique()
    genes = muts_df["Order on oncoprint"].unique()
    
    for i in genes: 
        for j in pts: 
            ax_main.bar(x = j, bottom = i, height = bar_height, width=bar_width, capstyle='butt', color = "whitesmoke", edgecolor = "none", linewidth = 0.25, zorder = 10)
    
    for index, row in muts_df.iterrows():
        ax_main.bar(x = row["Samples_enumerated"], bottom = row["Order on oncoprint"], height = bar_height, width=bar_width, capstyle='butt', color = row["Mutation_color"], edgecolor = "none", linewidth = 0.25, zorder = 1000)
    else:
        # 1.2 Plot the individual mutations as scatter on top of the grid
        size_scatter = 6
        for index, row in muts_df.iterrows():
            if row["Total_mutations_nonsilent_per_sample_per_gene"] == 1: # if the patient has only 1 mutation in a given gene
                ax_main.scatter(x = row["Samples_enumerated"], y = row["Order on oncoprint"] + bar_height/2 , color = row["Mutation_color"], s = size_scatter, marker = "s",  lw = 0, zorder = 9999)
            elif row["Total_mutations_nonsilent_per_sample_per_gene"] == 2:
                if row["Mutation_order"] == 1: # first mutation in the pair
                    ax_main.scatter(x = row["Samples_enumerated"], y = row["Order on oncoprint"] - 0.12 + bar_height/2 , color = row["Mutation_color"], s = size_scatter, marker = "s", lw = 0, zorder = 9999)
                elif row["Mutation_order"] == 2: # second mutation in the pair
                    ax_main.scatter(x = row["Samples_enumerated"], y = row["Order on oncoprint"] + 0.12 + bar_height/2 , color = row["Mutation_color"], s = size_scatter, marker = "s", lw = 0, zorder = 9999)
            if row["Total_mutations_nonsilent_per_sample_per_gene"] == 3:
                if row["Mutation_order"] == 1: # first mutation in the triplet
                    ax_main.scatter(x=row["Samples_enumerated"], y=row["Order on oncoprint"] - 0.24 + bar_height/2, color=row["Mutation_color"], s=size_scatter, marker="s", lw=0, zorder=9999)
                elif row["Mutation_order"] == 2: # second mutation in the triplet
                    ax_main.scatter(x=row["Samples_enumerated"], y=row["Order on oncoprint"] + bar_height/2, color=row["Mutation_color"], s=size_scatter, marker="s", lw=0, zorder=9999)
                elif row["Mutation_order"] == 3: # third mutation in the triplet
                    ax_main.scatter(x=row["Samples_enumerated"], y=row["Order on oncoprint"] + 0.24 + bar_height/2, color=row["Mutation_color"], s=size_scatter, marker="s", lw=0, zorder=9999)
            elif row["Total_mutations_nonsilent_per_sample_per_gene"] > 3:
                ax_main.scatter(x=row["Samples_enumerated"], y=row["Order on oncoprint"] + bar_height/2, marker="d", lw=1, facecolor="black", edgecolors="none", s = size_scatter, zorder = 9999)
    return ax_main



def plot_oncoprint(ax_main, muts_df, bar_height, bar_width):
    """
    If plot_mutations_as_squares is set to False, muts will be plotted as colors.
    """
    pts = muts_df["Samples_enumerated"].unique()
    genes = muts_df["Order on oncoprint"].unique()
    
    for i in genes: 
        for j in pts: 
            ax_main.bar(x = j, bottom = i, height = bar_height, width=bar_width, capstyle='butt', color = "whitesmoke", edgecolor = "none", linewidth = 0.25, zorder = 10)
    
    # 1.2 Plot the individual mutations as scatter on top of the grid
    size_scatter = 6
    for index, row in muts_df.iterrows():
        if row["Total_mutations_nonsilent_per_sample_per_gene"] == 1: # if the patient has only 1 mutation in a given gene
            ax_main.scatter(x = row["Samples_enumerated"], y = row["Order on oncoprint"] + bar_height/2 , color = row["Mutation_color"], s = size_scatter, marker = "s",  lw = 0, zorder = 9999)
        elif row["Total_mutations_nonsilent_per_sample_per_gene"] == 2:
            if row["Mutation_order"] == 1: # first mutation in the pair
                ax_main.scatter(x = row["Samples_enumerated"], y = row["Order on oncoprint"] - 0.12 + bar_height/2 , color = row["Mutation_color"], s = size_scatter, marker = "s", lw = 0, zorder = 9999)
            elif row["Mutation_order"] == 2: # second mutation in the pair
                ax_main.scatter(x = row["Samples_enumerated"], y = row["Order on oncoprint"] + 0.12 + bar_height/2 , color = row["Mutation_color"], s = size_scatter, marker = "s", lw = 0, zorder = 9999)
        if row["Total_mutations_nonsilent_per_sample_per_gene"] == 3:
            if row["Mutation_order"] == 1: # first mutation in the triplet
                ax_main.scatter(x=row["Samples_enumerated"], y=row["Order on oncoprint"] - 0.24 + bar_height/2, color=row["Mutation_color"], s=size_scatter, marker="s", lw=0, zorder=9999)
            elif row["Mutation_order"] == 2: # second mutation in the triplet
                ax_main.scatter(x=row["Samples_enumerated"], y=row["Order on oncoprint"] + bar_height/2, color=row["Mutation_color"], s=size_scatter, marker="s", lw=0, zorder=9999)
            elif row["Mutation_order"] == 3: # third mutation in the triplet
                ax_main.scatter(x=row["Samples_enumerated"], y=row["Order on oncoprint"] + 0.24 + bar_height/2, color=row["Mutation_color"], s=size_scatter, marker="s", lw=0, zorder=9999)
        elif row["Total_mutations_nonsilent_per_sample_per_gene"] > 3:
            ax_main.scatter(x=row["Samples_enumerated"], y=row["Order on oncoprint"] + bar_height/2, marker="d", lw=1, facecolor="black", edgecolors="none", s = size_scatter, zorder = 9999)
    return ax_main

def beautify_OP_ax(ax_main, gene_order_df, DF_Samples_enumerated, xlabel, fontsize = 10):
    ax_main.set_yticks(gene_order_df["Order on oncoprint"] + bar_height/2)
    ax_main.set_yticklabels(gene_order_df["Gene"], fontsize = fontsize)
    ax_main.set_ylim((gene_order_df["Order on oncoprint"].min() - bar_height/2, gene_order_df["Order on oncoprint"].max() + 1))
    
    ax_main.tick_params(axis='y', which='both', left=False)
    ax_main.tick_params(axis='y', which='both', pad=-2) # brings the y ticklabels closer to the plot
    
    ax_main.set_xticks(DF_Samples_enumerated["Samples_enumerated"])
    ax_main.set_xticklabels(DF_Samples_enumerated["Patient_id"], rotation = 90, fontsize = 6)
    ax_main.set_xlim((DF_Samples_enumerated["Samples_enumerated"].min() - 1, DF_Samples_enumerated["Samples_enumerated"].max() + 1))
    
    ax_main.tick_params(axis='x', which='both', bottom=False)
    ax_main.tick_params(axis='x', which='both', pad=-2) # brings the y ticklabels closer to the plot
    ax_main.set_xlabel(xlabel, fontsize = 10)
    
    remove_frame(ax_main, ["top", "bottom", "right", "left"])
    return ax_main

mut_dict = {
    "Frameshift indel": '#FFC907',
    "Nonframeshift indel": '#a9a9a9',
    "Missense": '#79B443',
    "Stopgain": '#BD4398',
    "Splicing": "darkorange"}

# arial_path = '/home/vpc/.fonts/arial_regular.ttf'
# font_properties = FontProperties(fname=arial_path)
# plt.rcParams['font.family'] = font_properties.get_name()


mpl.rcParams['font.size'] = 8
mpl.rcParams['text.color'] = 'k'
mpl.rcParams['legend.fontsize'] = 10
mpl.rcParams['legend.handletextpad'] = '0.8'
mpl.rcParams['legend.labelspacing'] = '0.4'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['axes.linewidth'] = 0.5
mpl.rcParams['xtick.labelsize'] = 7
mpl.rcParams['ytick.labelsize'] = 7
mpl.rcParams['axes.labelsize'] = 10
mpl.rcParams['xtick.major.size'] = 3
mpl.rcParams['ytick.major.size'] = 3

mpl.rcParams['legend.fontsize'] = 6
mpl.rcParams['legend.handletextpad'] = '0.7'
mpl.rcParams['legend.labelspacing'] = '0.2'
plt.rcParams['legend.handlelength'] = 0.4
plt.rcParams['legend.handleheight'] = 0.70




# LOAD DATASETS
DIR_working = "/groups/wyattgrp/users/amunzur/pipeline"
PATH_sample_information = os.path.join(DIR_working, "resources/sample_lists/sample_information.tsv")
PATH_gene_categories = os.path.join(DIR_working, "resources/panel/chip_panel_gene_categories.tsv")
PATH_mutation_ctfractions = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/ctDNA_fractions.csv"

all_vars_chip = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/chip.csv"))
all_vars_chip["Consequence"] = all_vars_chip["Consequence"].replace({"Frameshift deletion":"Frameshift indel", "Frameshift insertion":"Frameshift indel", "Nonframeshift deletion":"Nonframeshift indel", "Nonframeshift insertion":"Nonframeshift indel"})

bladder_chip = all_vars_chip[(all_vars_chip["Diagnosis"] == "Bladder") & (all_vars_chip["Timepoint"] == "Baseline")]
kidney_chip = all_vars_chip[(all_vars_chip["Diagnosis"] == "Kidney") & (all_vars_chip["Timepoint"] == "Baseline")]

# all_vars_somatic = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/somatic_SSCS2_curated_complete.csv"))
# all_vars_somatic["Consequence"]= all_vars_somatic["Consequence"].replace({"Frameshift deletion":"Frameshift indel", "Frameshift insertion":"Frameshift indel", "Nonframeshift deletion":"Nonframeshift indel", "Nonframeshift insertion":"Nonframeshift indel"})
# bladder_somatic = all_vars_somatic[all_vars_somatic["Diagnosis"] == "Bladder"]
# kidney_somatic = all_vars_somatic[all_vars_somatic["Diagnosis"] == "Kidney"]

def prepare_datasets_for_plotting(main_muts_df, mut_dict, PATH_gene_categories, PATH_mutation_ctfractions, MSK_style = False):
    # DATASET 1: MUTATION INFORMATION FOR THE MAIN ONCOPRINT PANEL
    muts_df = main_muts_df.copy()[["Patient_id", "Sample_name_t", "Timepoint", "Gene", "Consequence", "VAF_n", "Chrom"]]
    muts_df["Consequence"] = muts_df["Consequence"].str.lower().replace({"startloss": "missense", "frameshift indel": "indel", "frameshift deletion": "indel","frameshift insertion": "indel", "non-frameshift deletion": "indel", "non-frameshift insertion": "indel"})
    muts_df = muts_df.replace("U2AF1\\x3bU2AF1L5", "U2AF1")
    muts_df = enumerate_col(muts_df, "Gene", "Genes_enumerated")
    muts_df["Total_mutations_nonsilent_per_sample_per_gene"] = muts_df.groupby(['Sample_name_t', 'Gene'])['Gene'].transform('count')
    muts_df["Mutation_order"] = muts_df["Mutation_order"] = muts_df.groupby(['Sample_name_t', 'Gene']).cumcount() + 1
    
    if MSK_style:
        # determine colors
        mut_points = {'stopgain': 4, 'missense': 1, 'indel': 3, 'splicing': 2}        
        mut_colors = {3: '#FFC907', 1: '#79B443', 4: '#BD4398', 2: "darkorange"}
        muts_df["mut_points"] = muts_df["Consequence"].map(mut_points)
        max_pts = muts_df.groupby("Patient_id")["mut_points"].max().reset_index()
        max_pts["Mutation_color"] = max_pts["mut_points"].map(mut_colors)
        muts_df = muts_df.merge(max_pts[["Patient_id", "Mutation_color"]].drop_duplicates(), on = "Patient_id")
    else:
        # generate dict to to map to mut types colors
        set1_palette = sns.color_palette("Set1", n_colors=len(muts_df["Consequence"].unique()))
        consequence_colors = dict(zip(sorted(muts_df["Consequence"].unique()), set1_palette))
        mut_dict.update((consequence, color) for consequence, color in consequence_colors.items() if consequence not in mut_dict)    
        muts_df["Mutation_color"] = muts_df["Consequence"].map(mut_dict)
    
    # Helps labelling the y axis ticks (gene names)    
    gene_order_df = pd.read_csv(PATH_gene_categories, sep = "\t").rename(columns = {"Panel genes": "Gene"})
    gene_order_df = gene_order_df[~pd.isnull(gene_order_df["Order on oncoprint"])]
    genes_present_in_df = muts_df["Gene"].unique()
    gene_order_df = gene_order_df[gene_order_df["Gene"].isin(genes_present_in_df)]
    gene_order_df["Order on oncoprint"] = range(len(gene_order_df), 0, -1)
    muts_df = muts_df.merge(gene_order_df, how = "left")
    muts_df = muts_df[~pd.isnull(muts_df["Order on oncoprint"])]
    
    # DATASET 5 - CHIP FRACTIONS PER PATIENT
    chip_fr = muts_df.copy()[["Patient_id", "VAF_n", "Chrom"]]
    chip_fr = chip_fr.loc[chip_fr.groupby('Patient_id')['VAF_n'].idxmax()].reset_index(drop = True)
    chip_fr["CHIP_fr"] = chip_fr.apply(lambda row: row["VAF_n"] * 2 if row["Chrom"] != "chrX" else row["VAF_n"], axis=1)
    chip_fr = chip_fr.sort_values(by = "CHIP_fr", inplace = False, ascending = False)
    chip_fr = enumerate_col(chip_fr, "Patient_id", "Samples_enumerated")
    
    DF_Samples_enumerated = chip_fr[["Patient_id", "Samples_enumerated"]].drop_duplicates().reset_index(drop = True).sort_values(by = "Samples_enumerated")
    muts_df = muts_df.merge(DF_Samples_enumerated, how = "left")
    
    if MSK_style:
        muts_df.sort_values("mut_points")
        binary_df = muts_df.pivot_table(index='Patient_id', columns='Gene', values='Consequence', aggfunc='count')

# Fill NaN values with 0 (absence of mutation) and non-NaN values with 1 (presence of mutation)
binary_df = binary_df.notnull().astype(int)

# Reset index to make 'Patient_id' a column again
binary_df.reset_index(inplace=True)
binary_df.sort_values(by=["DNMT3A", "TET2", "ASXL1", "JAK2", "ATM", "BRCA2", "PPM1D", "CHEK2", "TP53"])
sorted_df = muts_df.sort_values()

def custom_sort(row):
    gene_order = ["DNMT3A", "TET2", "ASXL1", "JAK2", "ATM", "BRCA2", "PPM1D", "CHEK2", "TP53"]
    return [gene_order.index(gene) for gene in row.index[2:]]

# Sort the DataFrame using the custom sorting function
sorted_df = binary_df.sort_values(by=binary_df.columns.tolist(), key=custom_sort)


        DF_Samples_enumerated = chip_fr[["Patient_id", "Samples_enumerated"]].drop_duplicates().reset_index(drop = True).sort_values(by = "Samples_enumerated")
        muts_df = muts_df.merge(DF_Samples_enumerated, how = "left")

    
    # DATASET 2: FOR THE SIDEWAYS STACKED BARCHART FOR MUTATION COUNTS
    mut_counts_df = muts_df[["Gene", "Consequence"]].pivot_table(index='Gene', columns='Consequence', aggfunc=len, fill_value=0).reset_index()
    mut_counts_df = gene_order_df.merge(mut_counts_df, how = "left").fillna(0)
    mut_counts_df = mut_counts_df.sort_values(by="Order on oncoprint", inplace=False)
    mut_counts_df["Order on oncoprint"] = mut_counts_df["Order on oncoprint"] + bar_height/2
    mut_counts_df.set_index(["Order on oncoprint"], inplace=True)
    del mut_counts_df["Category"]
    del mut_counts_df["Gene"]
    
    # DATASET 2.5: THE VAFS TO SHOW ON TOP OF THE SIDEWAYS STACKED BAR CHART
    vaf_df = muts_df[["Gene", "VAF_n"]]
    vaf_df = gene_order_df.merge(vaf_df, how = "left").fillna(0)
    vaf_df = vaf_df.sort_values(by="Order on oncoprint", inplace=False)
    vaf_df["Order on oncoprint"] = vaf_df["Order on oncoprint"] + bar_height/2
    vaf_df.set_index(["Order on oncoprint"], inplace=True)
    del vaf_df["Category"]
    del vaf_df["Gene"]  
    
    # DATASET 3: MUTATION COUNTS FOR EACH PATIENT, ON TOP OF THE ONCOPRINT
    pt_counts = muts_df["Patient_id"].value_counts().reset_index().rename(columns = {"index": "Patient_id", "Patient_id": "Count"})
    pt_counts = pt_counts.merge(DF_Samples_enumerated, how = "left")
    
    # DATASET 4: CTDNA FRACTIONS
    ctdna_df = pd.read_csv(PATH_mutation_ctfractions)[["Patient_id", "Date_collected", "Mutation_ctDNA_fraction"]]
    ctdna_df["ctDNA fraction (%)"] = ctdna_df["Mutation_ctDNA_fraction"] * 100
    ctdna_df['Date_collected'] = pd.to_datetime(ctdna_df['Date_collected'])
    ctdna_df["Patient_id"] = ctdna_df["Patient_id"].str.replace("GU-", "")
    pts = DF_Samples_enumerated.merge(main_muts_df[["Patient_id", "Date_collected"]].drop_duplicates(), how = "left")
    pts['Date_collected'] = pd.to_datetime(pts['Date_collected'], format='%Y%b%d')
    pts_ctdna = pts.merge(ctdna_df, how = "left", on = ["Patient_id", "Date_collected"])
    
    pts_ctdna_NO_nan = pts_ctdna[~pd.isnull(pts_ctdna["ctDNA fraction (%)"])]
    pts_ctdna_nan = pts_ctdna[pd.isnull(pts_ctdna["ctDNA fraction (%)"])]
    
    return [muts_df, gene_order_df, chip_fr, mut_counts_df, vaf_df, pt_counts, pts_ctdna_NO_nan, pts_ctdna_nan, DF_Samples_enumerated, mut_dict]

    

##############################################################
# SINGLE ONCOPRINT FOR KIDNEY + BLADDER SAMPLES
fig_width = 8.5
fig_height = 5
tc_height = 2
chip_fraction = 2
mutcounts_height = 2
main_height = all_vars_chip["Gene"].unique().shape[0]

fig = plt.figure(figsize = (fig_width, fig_height))
gs = mpl.gridspec.GridSpec(ncols = 2, nrows = 6, height_ratios = [tc_height, chip_fraction, mutcounts_height, main_height, 2, 3], width_ratios = [1, 0.1], hspace = 0, wspace = 0)

# Add and space out subplots
ax_tc = fig.add_subplot(gs[0, 0]) 
ax_mut_counts_patient = fig.add_subplot(gs[2, 0], sharex = ax_tc) 
ax_main = fig.add_subplot(gs[3, 0], sharex = ax_tc)
ax_chip_fr = fig.add_subplot(gs[1, 0], sharex = ax_main) 
ax_mut_counts = fig.add_subplot(gs[3, 1], sharey = ax_main) 

bar_height = 0.9
bar_width = 0.8

[df, gene_order_df, chip_fr, mut_counts_df, vaf_df, pt_counts, pts_ctdna_NO_nan, pts_ctdna_nan, DF_Samples_enumerated, mut_dict_kidney] = prepare_datasets_for_plotting(all_vars_chip, mut_dict, PATH_gene_categories, PATH_mutation_ctfractions, MSK_style = True)
ax_main = plot_oncoprint_MSK_style(ax_main, df, bar_height, bar_width)
ax_main = beautify_OP_ax(ax_main, gene_order_df, DF_Samples_enumerated, xlabel = "Baseline samples")
[ax_mut_counts, ax_secondary] = plot_mut_counts_per_gene(mut_counts_df, vaf_df, ax_mut_counts)
ax_mut_counts_patient = plot_mut_counts_per_patient(ax_mut_counts_patient, pt_counts, bar_width)
ax_tc = plot_ctdna_fractions(ax_tc, pts_ctdna_NO_nan, pts_ctdna_nan, bar_width)
ax_chip_fr = plot_chip_fractions(ax_chip_fr, chip_fr, bar_width)
ax_mut_counts = add_legend(mut_dict_kidney, bbox_location = (2, 0.05), loc = 'lower right')
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/OP.pdf")
##############################################################












##############################################################
# RUNNING ONCOPRINT FOR BLADDER SAMPLES
fig_width = 8.5
fig_height = 5
tc_height = 2
chip_fraction = 2
mutcounts_height = 2
main_height = all_vars_chip["Gene"].unique().shape[0]

fig = plt.figure(figsize = (fig_width, fig_height))
gs = mpl.gridspec.GridSpec(ncols = 2, nrows = 6, height_ratios = [tc_height, chip_fraction, mutcounts_height, main_height, 2, 3], width_ratios = [1, 0.1], hspace = 0, wspace = 0)

# Add and space out subplots
ax_tc = fig.add_subplot(gs[0, 0]) 
ax_mut_counts_patient = fig.add_subplot(gs[2, 0], sharex = ax_tc) 
ax_main = fig.add_subplot(gs[3, 0], sharex = ax_tc)
ax_chip_fr = fig.add_subplot(gs[1, 0], sharex = ax_main) 
ax_mut_counts = fig.add_subplot(gs[3, 1], sharey = ax_main) 

bar_height = 0.9
bar_width = 0.8

[all_vars_chip, gene_order_df, chip_fr, mut_counts_df, vaf_df, pt_counts, pts_ctdna_NO_nan, pts_ctdna_nan, DF_Samples_enumerated, mut_dict_kidney] = prepare_datasets_for_plotting(bladder_chip, mut_dict, PATH_gene_categories, PATH_mutation_ctfractions)
ax_main = plot_oncoprint(ax_main, all_vars_chip, bar_height, bar_width)
ax_main = beautify_OP_ax(ax_main, gene_order_df, DF_Samples_enumerated, xlabel = "Baseline bladder samples")
[ax_mut_counts, ax_secondary] = plot_mut_counts_per_gene(mut_counts_df, vaf_df, ax_mut_counts)
ax_mut_counts_patient = plot_mut_counts_per_patient(ax_mut_counts_patient, pt_counts, bar_width)
ax_tc = plot_ctdna_fractions(ax_tc, pts_ctdna_NO_nan, pts_ctdna_nan, bar_width)
ax_chip_fr = plot_chip_fractions(ax_chip_fr, chip_fr, bar_width)
ax_mut_counts = add_legend(mut_dict_kidney, bbox_location = (2, 0.05), loc = 'lower right')
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/OP_bladder.pdf")
##############################################################

##############################################################
# RUNNING ONCOPRINT FOR KIDNEY SAMPLES
fig_width = 8.5
fig_height = 5
tc_height = 2
chip_fraction = 2
mutcounts_height = 2
main_height = all_vars_chip["Gene"].unique().shape[0]

fig = plt.figure(figsize = (fig_width, fig_height))
gs = mpl.gridspec.GridSpec(ncols = 2, nrows = 6, height_ratios = [tc_height, chip_fraction, mutcounts_height, main_height, 2, 3], width_ratios = [1, 0.1], hspace = 0, wspace = 0)

# Add and space out subplots
ax_tc = fig.add_subplot(gs[0, 0]) 
ax_mut_counts_patient = fig.add_subplot(gs[2, 0], sharex = ax_tc) 
ax_main = fig.add_subplot(gs[3, 0], sharex = ax_tc)
ax_chip_fr = fig.add_subplot(gs[1, 0], sharex = ax_main) 
ax_mut_counts = fig.add_subplot(gs[3, 1], sharey = ax_main) 

bar_height = 0.8
bar_width = 0.7

[all_vars_chip, gene_order_df, chip_fr, mut_counts_df, vaf_df, pt_counts, pts_ctdna_NO_nan, pts_ctdna_nan, DF_Samples_enumerated, mut_dict_kidney] = prepare_datasets_for_plotting(kidney_chip, mut_dict, PATH_gene_categories, PATH_mutation_ctfractions)
ax_main = plot_oncoprint(ax_main, all_vars_chip, bar_height, bar_width)
ax_main = beautify_OP_ax(ax_main, gene_order_df, DF_Samples_enumerated, xlabel = "Baseline kidney samples")
[ax_mut_counts, ax_secondary] = plot_mut_counts_per_gene(mut_counts_df, vaf_df, ax_mut_counts)
ax_mut_counts_patient = plot_mut_counts_per_patient(ax_mut_counts_patient, pt_counts, bar_width)
ax_tc = plot_ctdna_fractions(ax_tc, pts_ctdna_NO_nan, pts_ctdna_nan, bar_width)
ax_chip_fr = plot_chip_fractions(ax_chip_fr, chip_fr, bar_width)
ax_mut_counts = add_legend(mut_dict_kidney, bbox_location = (2, 0.05), loc = 'lower right')
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/OP_kidney.pdf")


