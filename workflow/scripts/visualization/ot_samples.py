import pandas as pd
import numpy as np
import os 
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from scipy.stats import fisher_exact
from scipy.stats import fisher_exact
from scipy.stats import binom
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


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

def get_baseline_and_ot_vafs(muts_df_main, gene, pts_with_ot, return_difference = False): 
    """
    Given a pt id and a gene name, follows all the mutations in that gene and returns their vaf in all samples of that patient.
    If return_difference = True, then returns the VAF difference in the two timepoints.
    """
    
    muts_df = muts_df_main.copy()
    muts_df.loc[(muts_df["Dependent"] == True), "VAF_t"] = muts_df.loc[(muts_df["Dependent"] == True), "VAF_t"]*100
    muts_df.loc[(muts_df["Dependent"] == True), "VAF_n"] = muts_df.loc[(muts_df["Dependent"] == True), "VAF_n"]*100
    
    # Determine if this is a CHIP df or a ctDNA mutation df
    if "Status" in muts_df.columns:
        col = "VAF_t"
    else:
        col = "VAF_n"
    
    muts_df = muts_df[(muts_df["Patient_id"].isin(pts_with_ot)) & (muts_df["Gene"] == gene)]
    results_dict = {}
    CI_dict = {}
    for pt in pts_with_ot:
        pt_df = muts_df[muts_df["Patient_id"] == pt]
        # Go through each mutation in the gene
        for i, group in pt_df.groupby("Protein_annotation"): 
            mut_name = group["Protein_annotation"].unique()[0]
            group["Date_collected"] = pd.to_datetime(group["Date_collected"], format = '%Y%b%d')
            group = group.sort_values(by = "Date_collected")
            vafs_list = group[col].tolist()[0:2] # only showing the extra timepoint after the baseline
            if return_difference:
                vafs_list = calculate_distances(vafs_list)
            
            # Generate a mutation identifier
            mut_identifier = f"{pt}_{gene}_{mut_name}"
            results_dict[mut_identifier] = vafs_list
            # results_list.append(vafs_list)
            
            # Now calculate the confidence interval for each time point
            depth_list = group["Depth_n"].tolist()[0:2]
            vaf_list_nonperc = [i/100 for i in vafs_list]
            mutant_reads_list = [round(v * d) for v, d in zip(vaf_list_nonperc, depth_list)]
            
            def clopper_pearson_ci(mutant_reads, depth, confidence_level):
                alpha = 1 - confidence_level
                # Binomial confidence interval
                lower_bound, upper_bound = binom.interval(confidence_level, depth, mutant_reads / depth)
                return lower_bound / depth*100, upper_bound / depth*100
            
            # Calculate confidence intervals for each VAF
            mutation_list = []
            for mutant_reads, depth in zip(mutant_reads_list, depth_list):
                if mutant_reads == 0:
                    mutant_reads+=1
                lower, upper = clopper_pearson_ci(mutant_reads, depth, confidence_level = 0.95)
                ci_list = [lower, upper]
                mutation_list.append(ci_list)
            CI_dict[mut_identifier] = mutation_list
    return(results_dict, CI_dict)

def calculate_distances(lst):
    """
    computes the distances (differences) between consecutive elements in a given list.
    """
    distances = [lst[i+1] - lst[i] for i in range(len(lst) - 1)]
    return distances

def order_mutations(mutations):
    """
    Orders mutations based on delta.
    
    Parameters:
        mutations (dict): A dictionary where keys are mutation identifiers 
                          and values are lists with two elements: [value1, value2].
                          
    Returns:
        dict: A dictionary of mutations sorted by their delta (value2 - value1) in descending order.
    """
    # Calculate delta for each mutation
    deltas = {mutation: values[1] - values[0] for mutation, values in mutations.items()}
    
    # Sort mutations based on delta value in descending order
    sorted_mutations = sorted(deltas.items(), key=lambda x: x[1], reverse=True)
    
    # Create a new dictionary maintaining the same format as the input
    sorted_mutations_dict = {mutation: mutations[mutation] for mutation, delta in sorted_mutations}
    
    return sorted_mutations_dict

def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def plot_vaf_changes(sorted_mutations, CIs, ax):
    """
    Plots the VAF change from baseline to the next OT sample.
    
    Parameters:
        mutations (dict): A dictionary where keys are mutation identifiers 
                          and values are lists containing two VAF values [baseline_vaf, timepoint2_vaf].
        ax (matplotlib.axes.Axes): The axes on which to plot the data.
    
    Returns:
        ax: The modified axes object with the plot.
    """
    for i, (mutation_id, vaf_values) in enumerate(sorted_mutations.items()):
        ci_lower_baseline = CIs[mutation_id][0][0]
        ci_upper_baseline = CIs[mutation_id][0][1]
        ci_lower_OT = CIs[mutation_id][1][0]
        ci_upper_OT = CIs[mutation_id][1][1]
        
        baseline_vaf = vaf_values[0]
        timepoint2_vaf = vaf_values[1]
        
        # Check for overlap using absolute values
        overlap_baseline_OT = getOverlap([ci_lower_baseline, ci_upper_baseline], [ci_lower_OT, ci_upper_OT])
        
        if timepoint2_vaf > baseline_vaf:
            marker = "^"  # Upward triangle
            if overlap_baseline_OT > 0:
                color = "#ffb3b3" # light red
            else:
                color = "#ff0000" # dark red
        else:
            marker = "v"  # Downward triangle
            if overlap_baseline_OT > 0:
                color = "#b3b7ff"
            else:
                color = "#4955ff" # dark blue
        
        # Plot vertical lines
        ax.vlines(x=i, ymin=min(baseline_vaf, timepoint2_vaf), ymax=max(baseline_vaf, timepoint2_vaf),
                   colors=color, lw=0.8)
        
        # Plot scatter points for timepoint2_vaf
        ax.scatter(i, timepoint2_vaf, color=color, s=15, marker=marker, zorder=50, edgecolors='none')
        
        # Annotate mutations called in OT only
        if baseline_vaf < 0.25:
            higher_value = max(baseline_vaf, timepoint2_vaf)
            ax.text(i, higher_value, "*", ha='center', va='bottom', color="black")
    
    # Aesthetic settings
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylabel("VAF Δ %")
    ax.set_xlabel("")
    ax.set_xticks([])
    ax.set_xticklabels([])
    return ax

def plot_inset_pie(sorted_mutations, CIs, ax):
    """
    Top right adds a new piechart.
    """
    n_overlapping_incr = 0
    n_nonoverlapping_incr = 0
    n_overlapping_decr = 0
    n_nonoverlapping_decr = 0
    for i, (mutation_id, vaf_values) in enumerate(sorted_mutations.items()):
        ci_lower_baseline = CIs[mutation_id][0][0]
        ci_upper_baseline = CIs[mutation_id][0][1]
        ci_lower_OT = CIs[mutation_id][1][0]
        ci_upper_OT = CIs[mutation_id][1][1]
        
        baseline_vaf = vaf_values[0]
        timepoint2_vaf = vaf_values[1]
        
        # Check for overlap using absolute values
        overlap_baseline_OT = getOverlap([ci_lower_baseline, ci_upper_baseline], [ci_lower_OT, ci_upper_OT])
        
        if timepoint2_vaf > baseline_vaf:
            if overlap_baseline_OT > 0:
                n_overlapping_incr+=1
            else:
                n_nonoverlapping_incr+=1
        else:
            if overlap_baseline_OT:
                n_overlapping_decr+=1
            else:
                n_nonoverlapping_decr+=1
    
    # Add the inset plot to top right
    sizes = [n_overlapping_incr, n_nonoverlapping_incr, n_overlapping_decr, n_nonoverlapping_decr]
    labels = ['Overlapping Increment', 'Non-overlapping Increment', 'Overlapping Decrement', 'Non-overlapping Decrement']
    colors = [
        "#ffb3b3",  # Overlapping & increasing
        "#ff0000",  # Non-overlapping & increasing
        "#b3b7ff",  # Overlapping & decreasing
        "#4955ff" # Non-overlapping & decreasing
    ]
    inset_ax = inset_axes(ax, width=0.4, height=0.4, loc='upper right')
    wedges, _ = inset_ax.pie(sizes, colors=colors, autopct=None, startangle=90)
    
    # Manually add percentages just outside the pie chart
    for i, w in enumerate(wedges):
        if sizes[i] > 0:  # Only show percentages for non-zero sizes
            # Get the angle for the wedge
            angle = (w.theta1 + w.theta2) / 2.0
            # Calculate percentage value
            percentage = f'{int(sizes[i] / sum(sizes) * 100)}%'  
            # Calculate label position
            x = np.cos(np.radians(angle)) * 1.2  # Adjust 1.2 to move labels further out
            y = np.sin(np.radians(angle)) * 1.2
            inset_ax.text(x, y, percentage, ha='center', va='center', fontsize=6)
        
    inset_ax.axis('equal')
    
    return(ax, inset_ax)

def plot_CIs(CIs, sorted_mutations, ax):
    """
    Plots the confidence intervals for mutations.
    """
    
    # Order the CIs the same way as sorted mutations
    ordered_CIs = {}
    
    # Iterate over the keys in sorted_mutations to maintain the order
    for mutation in sorted_mutations.keys():
        if mutation in CIs:  # Ensure the mutation is in the CIs dictionary
            ordered_CIs[mutation] = CIs[mutation]
    
    # Plot the CIs
    for i, (mutname, ci_values) in enumerate(ordered_CIs.items()):
        print(i)
        
        baseline_ci_lower = ci_values[0][0]
        baseline_ci_upper = ci_values[0][1]
        baseline_vaf = sorted_mutations[mutname][0]
        
        ot_ci_lower = ci_values[1][0]
        ot_ci_upper = ci_values[1][1]
        ot_vaf = sorted_mutations[mutname][1]
        
        ax.errorbar(x = i-0.2, y = baseline_vaf, yerr=[[baseline_vaf - baseline_ci_lower], [baseline_ci_upper - baseline_vaf]], fmt='none', color='gray', capsize=1, elinewidth = 0.5, capthick = 0.5, zorder = 0)
        ax.errorbar(x = i+0.2, y = ot_vaf, yerr=[[ot_vaf - ot_ci_lower], [ot_ci_upper - ot_vaf]], fmt='none', color='gray', capsize=1, elinewidth = 0.5, capthick = 0.5, zorder = 0)
    
    return(ax)

def plot_the_bar(chip, ax, ax_legend):
    """
    In the form of barchart plots the change in distribution of CH genes in baseline vs OT samples in bladder.
    """
    df_baseline = chip[(chip["Dependent"] == False) & (chip["Timepoint"] == "Baseline")][["Patient_id", "Timepoint", "Gene"]].reset_index(drop = True)
    df_OT = chip[(chip["Dependent"] == False) & (chip["Timepoint"] == "During treatment")][["Patient_id", "Timepoint", "Gene"]].reset_index(drop = True)
    
    # Limit the analysis to pts with OT samples
    df_baseline = df_baseline[df_baseline["Patient_id"].isin(df_OT["Patient_id"])].reset_index(drop = True)
    
    # Distribution of genes
    df_baseline.loc[df_baseline["Gene"].isin(["STAG2", "CBL", "PBRM1", "MYD88", "JAK2", "EZH2", "TERT", "GNAS", "RAD21", "ERBB2", "RUNX1", "KMT2D", "SETD2"]), "Gene"] = "Other"
    df_OT.loc[df_OT["Gene"].isin(["STAG2", "CBL", "PBRM1", "MYD88", "JAK2", "EZH2", "TERT", "GNAS", "RAD21", "ERBB2", "RUNX1", "KMT2D", "SETD2"]), "Gene"] = "Other"
    
    # Assign shades for each gene in each category
    gene_colors = {
        "Other": "Silver",
        
        "SRSF2": "#87CEFA",   # Lighter dodgerblue
        "SF3B1": "#1E90FF",   # Base dodgerblue
        "SH2B3": "#104E8B",    # Darker dodgerblue
        
        "ATM":"#ECAFAF",
        "CHEK2":"#E08B8B",
        "PPM1D":"#D46666",
        "TP53": "#943434",    # Darker darkred
        
        "ASXL1": "#339999",   # Darker teal     
        "TET2": "#008080",    # Lighter teal
        "DNMT3A": "#004c4c",  # Base teal  
    }
    
    # First calculate the fraction for each gene
    n_total_baseline = df_baseline.shape[0]
    n_total_ot = df_OT.shape[0]
    
    baseline_dict = {}
    ot_dict = {}
    
    baseline_dict_n = {}
    ot_dict_n = {}
    
    for i, gene in enumerate(gene_colors):
        n_gene_baseline = df_baseline[df_baseline["Gene"] == gene].shape[0]
        perc_gene_baseline = n_gene_baseline/n_total_baseline*100
        baseline_dict[gene] = perc_gene_baseline
        baseline_dict_n[gene] = n_gene_baseline
        
        n_gene_ot = df_OT[df_OT["Gene"] == gene].shape[0]
        perc_gene_ot = n_gene_ot/n_total_ot*100
        ot_dict[gene] = perc_gene_ot
        ot_dict_n[gene] = n_gene_ot
    
    # Plotting baseline
    bottom = 0
    for i, (gene_name, y_value) in enumerate(baseline_dict.items()):
        ax.bar(0, y_value, bottom = bottom, color=gene_colors[gene_name], edgecolor = "None") # Plot baseline
        yval_text = bottom+y_value*0.5
        if baseline_dict_n[gene_name] > 0:
            ax.text(0, yval_text, baseline_dict_n[gene_name], ha='center', va='center', color="white", fontsize = 6)
        bottom+=y_value
    
    # Plotting OT
    bottom = 0
    for i, (gene_name, y_value) in enumerate(ot_dict.items()):
        ax.bar(1, y_value, bottom = bottom, color=gene_colors[gene_name], edgecolor = "None") # Plot OT
        yval_text = bottom+y_value*0.5
        if ot_dict_n[gene_name] > 0:
            ax.text(1, yval_text, ot_dict_n[gene_name], ha='center', va='center', color="white", fontsize = 6)
        bottom+=y_value
    
    # AES
    ax.set_xticks([0, 1])
    ax.set_ylim((0, 100))
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xticklabels(["Baseline", "OT"])
    ax.set_yticklabels(["0", "25", "50", "75", "100"])
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlabel("")
    ax.set_ylabel("% of mutations")
    ax.set_xlim((-0.5, 1.5))
    
    # Do the legend
    legend_columns = {
        "Col 1": ["DNMT3A", "TET2", "ASXL1", ""],
        "Col 2": ["TP53", "PPM1D", "CHEK2", "ATM"],
        "Col 3": ["SH2B3", "SF3B1", "SRSF2", "Other"]
    }
    
    legend_handles = []
    for col in legend_columns.values():
        for gene in col:
            if len(gene) == 0:
                color = "white"
            else:
                color = gene_colors[gene]
            handle = plt.Line2D([0], [0], marker='s', color=color, label=gene, markersize=5, linestyle='')
            legend_handles.append(handle)
    
    ax_legend.spines[["top", "right", "left", "bottom"]].set_visible(False)
    ax_legend.tick_params(axis='both', bottom=False, labelbottom=False, left=False, labelleft=False)
    
    # Place the legend with specified columns
    ax_legend.legend(handles=legend_handles, loc="upper right", ncol=3, frameon=False, fontsize=6, handlelength=1, handletextpad=0.1, columnspacing=0.4)
    return(ax)    

PATH_sample_information = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"
path_chip = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/chip.csv"
path_somatic = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/somatic.csv"

sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
pts_with_ot = sample_info[sample_info["Timepoint"] == "During treatment"]["Patient_id"].unique().tolist()

chip = pd.read_csv(path_chip)
somatic = pd.read_csv(path_somatic)

genes_list = ["DNMT3A", "TET2", "ASXL1", "TP53", "PPM1D", "CHEK2"]

###############################################################
# CHIP
fig = plt.figure(figsize=(8, 4))
outer_gs = gridspec.GridSpec(1, 2, width_ratios=[0.5, 3], hspace = 0, wspace = 0.2)

# Generate the grid based on the number of arrows per gene
top_grid_width_ratios = []
bottom_grid_with_ratios = []
for i, gene in enumerate(genes_list):
    mutations, CIs = get_baseline_and_ot_vafs(chip, gene, pts_with_ot, return_difference = False) # get the vafs all mutations in a given gene from baseline to timepoint2
    sorted_mutations = order_mutations(mutations)
    n_arrows = len(sorted_mutations.keys())
    
    if i < 3:
        top_grid_width_ratios.append(n_arrows)
    else:
        bottom_grid_with_ratios.append(n_arrows)

# Set up the grid for VAF change gs
vaf_change_gs = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios=[1, 1], hspace = 0.3, wspace = 0.2, subplot_spec=outer_gs[1])
top_row = gridspec.GridSpecFromSubplotSpec(1, 3, width_ratios=top_grid_width_ratios, hspace = 0.3, wspace = 0.2, subplot_spec=vaf_change_gs[0])
bottom_row_gs = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios=[1, 0.05], hspace = 0, wspace = 0, subplot_spec=vaf_change_gs[1])
bottom_row_vafchange = gridspec.GridSpecFromSubplotSpec(1, 3, width_ratios=bottom_grid_with_ratios, hspace = 0.3, wspace = 0.2, subplot_spec=bottom_row_gs[0])
vaf_change_legend_ax = plt.subplot(bottom_row_gs[1])

dta_genes_muts_chip = []
ddr_genes_muts_chip = []

ax0 = plt.subplot(top_row[0])
ax1 = plt.subplot(top_row[1])
ax2 = plt.subplot(top_row[2])
ax3 = plt.subplot(bottom_row_vafchange[0])
ax4 = plt.subplot(bottom_row_vafchange[1]) 
ax5 = plt.subplot(bottom_row_vafchange[2]) 

ax_list = [ax0, ax1, ax2, ax3, ax4, ax5]

for gene_idx, (ax, gene) in enumerate(zip(ax_list, genes_list)):
    print(gene)
    mutations, CIs = get_baseline_and_ot_vafs(chip, gene, pts_with_ot, return_difference = False) # get the vafs all mutations in a given gene from baseline to timepoint2
    sorted_mutations = order_mutations(mutations)
    
    # Calculate median change
    median_change = round(np.median([m[1] - m[0] for m in mutations.values()]), 3)    
    ax.set_title(f"{gene}\nmedian Δ: {median_change}", loc = "left", fontstyle = "italic", fontsize = 8)
    
    # Append to list for doing stats later
    if gene in ["DNMT3A", "TET2", "ASXL1"]:
        dta_genes_muts_chip.extend(sorted_mutations)
    elif gene in ["TP53", "PPM1D", "CHEK2"]:
        ddr_genes_muts_chip.extend(sorted_mutations)
    
    # Plotting and aes
    ax = plot_vaf_changes(sorted_mutations,CIs, ax)
    ax = plot_CIs(CIs, sorted_mutations, ax)
    if gene not in ["DNMT3A", "TP53"]:
        ax.set_ylabel("")
    if gene == "ASXL1":
        ax.set_yticks([0, 4, 8])
        ax.set_yticklabels(["0", "4", "8"])
    elif gene == "PPM1D":
        ax.set_yticks([0, 10, 20, 30, 40])
        ax.set_yticklabels(["0", "10", "20", "30", "40"])
    elif gene == "CHEK2":
        ax.set_yticks([0, 1, 2, 3])
        ax.set_yticklabels(["0", "1", "2", "3"])
    
    # Now plot the inset pie charts
    ax, inset_ax = plot_inset_pie(sorted_mutations, CIs, ax)

# Add legend for VAF change
legend_dict = {
    "Non-overlapping CIs & increasing VAF": "#ff0000",
    "Overlapping CIs & increasing VAF": "#ffb3b3",
    "Non-overlapping CIs & decreasing VAF": "#4955ff",
    "Overlapping CIs & decreasing VAF": "#b3b7ff"
}

legend_colors = legend_dict.values()
legend_labels =  legend_dict.keys()
legend_handles = [plt.Line2D([0], [0], marker='s', color=color, label=label, markersize=5, linestyle='') for color, label in zip(legend_colors, legend_labels)]
vaf_change_legend_ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
vaf_change_legend_ax.tick_params(axis='both', bottom=False, labelbottom=False, left=False, labelleft = False)
vaf_change_legend_ax.legend(handles=legend_handles, ncol = 2, loc="best", frameon=False, fontsize = 6, handlelength=0.5, handletextpad = 0.5)

# Now do the barchart
barchart_gs = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios=[1, 0.1], hspace = 0.15, wspace = 0, subplot_spec=outer_gs[0])

ax_bar = plt.subplot(barchart_gs[0])
ax_legend = plt.subplot(barchart_gs[1])
ax_bar = plot_the_bar(chip, ax_bar, ax_legend)

outer_gs.tight_layout(fig)
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/OT_samples_vaf_change_chip.png")
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/OT_samples_vaf_change_chip.pdf")