import pandas as pd
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ks_2samp
import seaborn as sns
from scipy.stats import mannwhitneyu
import ast
from itertools import chain
from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests
import matplotlib.gridspec as gridspec
import random

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

def make_density_plot(values, title, color, ax, bw_adjust=0.05, annotate_median = True, alpha = 0.75):
    """
    Given a list of values make histogram.
    """
    sns.kdeplot(values, ax=ax, color=color, fill=True, alpha=alpha, bw_adjust=bw_adjust, linewidth = 3)
    if annotate_median:
        ax.text(median_value+10, ax.get_ylim()[1] * 0.9, f'Median: {int(median_value)}', color='black', ha='left')
        median_value = np.median(values)
        ax.axvline(median_value, color='black', linestyle='dashed', linewidth=1)
    ax.set_xlabel('Fragment length')
    ax.set_ylabel('Frequency')
    ax.set_title(title)
    ax.set_xlim((0, 400))
    ax.set_xticks([0, 100, 200, 300, 400])
    ax.set_xticklabels(["0", "100", "200", "300", "400"])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return(ax)

def make_histogram(values, title, color, color_kde, ax, annotate_median = True, alpha = 0.75, bw_adjust= 0.05, show_hist = True):
    """
    Given a list of values make histogram.
    """
    if show_hist:
        ax.hist(values, bins=250, alpha=0.75, color = color)
        kde_ax = ax.twinx()  # Create a twin axes that shares the same x-axis
        ax.spines["top"].set_visible(False)
        kde_ax.spines["top"].set_visible(False)
    else:
        kde_ax = ax
        kde_ax.spines[["top", "right"]].set_visible(False)
    sns.kdeplot(values, ax=kde_ax, color=color_kde, fill=False, bw_adjust=bw_adjust, linewidth=1, zorder=10)
    kde_ax.set_ylabel('Density')  # Label the secondary y-axis
    if annotate_median:
        median_value = np.median(values)
        ax.axvline(median_value, color='black', linestyle='dashed', linewidth=1)
        ax.text(median_value+10, ax.get_ylim()[1] * 0.9, f'Median: {int(median_value)}', color='black', ha='left')
    ax.set_xlabel('Fragment length')
    ax.set_ylabel('Frequency')
    ax.set_title(title)
    ax.set_xlim((0, 400))
    ax.set_xticks([0, 100, 200, 300, 400])
    ax.set_xticklabels(["0", "100", "200", "300", "400"])
    return(ax)

def MWU_mutated_to_WT_fragments(path_to_df, mutated_col_name, WT_col_name): 
    """
    Runs MWU and corrects for multiple testing when comparing the fragments carrying mutated fragments to WT fragments.
    Returns an array of corrected p values. 
    """
    df = pd.read_csv(path_to_df)
    p_values = []   
    # For each ctDNA mutation, run MWU to test if the medians are different
    for i, row in df.iterrows():
        mutant_fr = ast.literal_eval(row[mutated_col_name])
        WT_fr = ast.literal_eval(row[WT_col_name])
        statistic, p_value = mannwhitneyu(mutant_fr, WT_fr)
        p_values.append(p_value)
    
    # Convert p-values to a numpy array
    p_values = np.array(p_values)   
    _, p_values_corrected_bh, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh') # correct for multiple testing
    
    # Add original and corrected p-values to the DataFrame
    df['p_value'] = p_values
    df['p_value_corrected_bh'] = p_values_corrected_bh
    
    # Determine the number of significant p-values
    significant_pvals_count = np.sum(p_values_corrected_bh < 0.05)
    print(f"Number of significant p-values after correction: {significant_pvals_count}") 
    return(df)  

def violin_box(data):
    df_plot = pd.DataFrame(data)
    
    # Map groups to numerical labels
    group_mapping = {group: i for i, group in enumerate(df_plot['Group'].unique())}
    df_plot['Group_num'] = df_plot['Group'].map(group_mapping)
    
    # Create the plot using the ax method
    fig, ax = plt.subplots(figsize=(10, 6)) 
    groups = df_plot['Group'].unique()
    colors = ['blue', 'green', 'red', 'orange']
    violin_width = 0.8  # Adjust the width of violins as needed
    
    violins = ax.violinplot(
        [df_plot[df_plot['Group_num'] == group_mapping[group]]['Fragment Length'] for group in groups], 
        positions=[group_mapping[group] for group in groups],
        widths=violin_width,
        showmedians=True
    )
    
    # Customize the violins' appearance
    for pc in violins['bodies']:
        pc.set_facecolor('grey')
        pc.set_alpha(1.0)  # Normal opacity
    
    # Annotate medians on the violins
    for group, color in zip(groups, colors):
        median_val = np.median(df_plot[df_plot['Group'] == group]['Fragment Length'])
        ax.text(group_mapping[group], median_val, f'{median_val:.2f}', color='black', ha='center', va='bottom')
    
    # Adding labels and title
    ax.set_xlabel('Group')
    ax.set_ylabel('Fragment Length')
    ax.set_title('Fragment Lengths by Group')
    ax.set_xticks(list(group_mapping.values()))
    ax.set_xticklabels(list(group_mapping.keys()))
    ax.legend()
    
    fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/fragmentomics/violin.png")

def create_box_plot_with_jitter_grouped(data, save_path):
    """
    # Function to create box plots with jittered points and median values
    """
    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(6, 3))
    
    # Define fragment types and diagnoses
    fragment_types = data['Fragment Type'].unique()
    diagnoses = data['Diagnosis'].unique()
    
    # Colors for the boxes and jitter points
    colors = {"Bladder": (143/255, 215/255, 239/255, 255/255), "Kidney": (239/255, 169/255, 143/255, 255/255)}
    jitter_colors = {"Bladder": "deepskyblue", "Kidney": "orangered"}
    
    # Create box plots
    box_data = []
    positions = []
    color_list = []
    
    for i, fragment_type in enumerate(fragment_types):
        for j, diagnosis in enumerate(diagnoses):
            subset = data[(data['Fragment Type'] == fragment_type) & (data['Diagnosis'] == diagnosis)]
            box_data.append(subset['Fragment Count'].values)
            positions.append(i + (j * 0.4) - 0.125)  # Adjust positions for dodge effect
            color_list.append(colors[diagnosis])
    
    boxprops = dict(linestyle='-', linewidth=0.5, color='black')
    whiskerprops = dict(color='black', linewidth=0.5)
    capprops = dict(color='black', linewidth=0.5)
    medianprops = dict(color='black', linewidth=0.5)
    
    bp = ax.boxplot(box_data, positions=positions, widths=0.4, patch_artist=True, boxprops=boxprops, whiskerprops=whiskerprops, capprops=capprops, medianprops=medianprops, showfliers=False)
    
    for patch, color in zip(bp['boxes'], color_list):
        patch.set_facecolor(color)
    
    # Add jittered strip plots
    for i, fragment_type in enumerate(fragment_types):
        for j, diagnosis in enumerate(diagnoses):
            subset = data[(data['Fragment Type'] == fragment_type) & (data['Diagnosis'] == diagnosis)]
            jitter_x = np.random.normal(loc=i + (j * 0.4) - 0.125, scale=0.02, size=subset.shape[0])
            ax.scatter(jitter_x, subset['Fragment Count'], alpha=0.7, color=jitter_colors[diagnosis], s=10, zorder = 10)
    
    # Add median values
    for i, fragment_type in enumerate(fragment_types):
        for j, diagnosis in enumerate(diagnoses):
            median_val = data[(data['Fragment Type'] == fragment_type) & (data['Diagnosis'] == diagnosis)]['Fragment Count'].median()
            pos = i + (j * 0.4) - 0.125
            ax.text(pos, median_val, f"{median_val:.1f}", ha='center', va='bottom', color='black', zorder = 100, fontsize = 7)
    
    # Set plot labels and title
    ax.set_xticks(range(len(fragment_types)))
    ax.set_xticklabels(fragment_types)
    ax.set_xlabel("Fragment Type")
    ax.set_ylabel("Fragment Count")
    ax.set_title("")
    
    # Remove legend
    ax.legend().set_visible(False)
    
    # Adjust layout and save plot
    fig.tight_layout()
    fig.savefig(save_path)
    plt.close(fig)

def create_box_plot_with_jitter(data, diagnosis, ax):
    """
    Function to create box plots with jittered points and median values for a specific diagnosis. Kidney and bladder not together.
    """
    
    # Filter data for the specified diagnosis
    filtered_data = data[data['Diagnosis'] == diagnosis]
    
    # Define fragment types for the filtered data
    fragment_types = filtered_data['Fragment Type'].unique()
    
    # Colors for the boxes and jitter points
    colors = {"Bladder": (143/255, 215/255, 239/255, 255/255), "Kidney": (239/255, 169/255, 143/255, 255/255)}
    jitter_color = {"Bladder": "deepskyblue", "Kidney": "orangered"}
    
    # Create box plots
    box_data = []
    positions = []
    color_list = []
    
    for i, fragment_type in enumerate(fragment_types):
        subset = filtered_data[filtered_data['Fragment Type'] == fragment_type]
        box_data.append(subset['Fragment Count'].values)
        positions.append(i)
        color_list.append(colors[diagnosis])
    
    bp = ax.boxplot(box_data, positions=positions, widths=0.6, patch_artist=True, boxprops=dict(linestyle='-', linewidth=0.5, color='black'), whiskerprops=dict(color='black', linewidth=0.5), capprops=dict(color='black', linewidth=0.5), medianprops=dict(color='black', linewidth=0.5), showfliers=False)
    
    for patch, color in zip(bp['boxes'], color_list):
        patch.set_facecolor(color)
    
    # Add jittered strip plots
    for i, fragment_type in enumerate(fragment_types):
        subset = filtered_data[filtered_data['Fragment Type'] == fragment_type]
        jitter_x = np.random.normal(loc=i, scale=0.1, size=subset.shape[0])
        ax.scatter(jitter_x, subset['Fragment Count'], alpha=0.7, color=jitter_color[diagnosis], s=5, zorder=10)
        # # Update x tick labels to include the number of dots (n) in the jitter
        # ax.set_xticks(range(len(fragment_types)))
        # ax.set_xticklabels([f"{fragment_type}\nn={subset.shape[0]}" for fragment_type in fragment_types], fontsize=8, rotation=0, ha='center')
    
    # Add median values
    for i, fragment_type in enumerate(fragment_types):
        median_val = filtered_data[filtered_data['Fragment Type'] == fragment_type]['Fragment Count'].median()
        formatted_median = f"{int(median_val)}" if median_val.is_integer() else f"{median_val:.1f}"
        ax.text(i, median_val, formatted_median, ha='center', va='bottom', color='black', zorder=100, fontsize=7)
    
    # Mann-Whitney U test and p-value display
    # ctDNA
    group1 = filtered_data[filtered_data['Fragment Type'] == "ctDNA WT"]['Fragment Count']
    group2 = filtered_data[filtered_data['Fragment Type'] == "ctDNA mutant"]['Fragment Count']
    statistic, p_value = mannwhitneyu(group1, group2, alternative='two-sided')
    ax.text(0.25, 0.95, f"p-value:\n{p_value:.2e}", ha='center', va='center', transform=ax.transAxes, fontsize=8)
    # CH
    group1 = filtered_data[filtered_data['Fragment Type'] == "CH WT"]['Fragment Count']
    group2 = filtered_data[filtered_data['Fragment Type'] == "CH mutant"]['Fragment Count']
    statistic, p_value = mannwhitneyu(group1, group2, alternative='two-sided')
    ax.text(0.75, 0.95, f"p-value:\n{p_value:.2f}", ha='center', va='center', transform=ax.transAxes, fontsize=8)
    
    # Set plot labels and title
    ax.set_xticks(range(len(fragment_types)))
    ax.set_xticklabels(fragment_types, fontsize = 8)
    ax.set_xlabel("Fragment Type")
    ax.set_ylabel("Fragment Count")
    ax.set_ylim((95, 355))
    ax.set_title(diagnosis)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    
    # Remove legend
    ax.legend().set_visible(False)
    return(ax)
    
def scatter_n_muts_and_fragment_length(melted_df, muts_df, type_of_muts, ax):
    """
    Is there a relationship between the number of mutations a person may have and the median mutated fragment length?
    melted_df: Df that tells the median fragment length per patient
    muts_df: mutation calls. could be CH or ctDNA muts.
    """
    if type_of_muts == "CH": 
        filtered = melted_df[melted_df["Fragment Type"] == "CH mutant"]
    else: 
        filtered = melted_df[melted_df["Fragment Type"] == "ctDNA mutant"]
    # combined
    n_muts = muts_df["Patient_id"].value_counts().reset_index()
    n_muts.columns = ["Patient_id", "n_muts"]
    combined = filtered.merge(n_muts, how = "inner", on = "Patient_id")

    # Plotting
    ax = ax.scatter(combined['n_muts'], combined['Fragment Count'], 
                     alpha=0.7, c=combined['Diagnosis'].apply(lambda x: {'Kidney': 'orangered', 'Bladder': 'deepskyblue'}[x]), label=combined['Diagnosis'])
    return(ax)

DIR_fragment_counts = "/groups/wyattgrp/users/amunzur/pipeline/results/fragmentomics"
DIR_figures = "/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/fragmentomics"
DIR_insertsizes = "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/insert_size/SSCS2"
PATH_sample_information = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"
PATH_vars_ch = ""

sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])

# Load the mut dfs
all_vars_chip = pd.read_csv(os.path.join("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/chip_SSCS2_curated_complete.csv"))
all_vars_chip = sample_info.merge(all_vars_chip, how = "inner")
all_vars_ctDNA = pd.read_csv(os.path.join("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/somatic_SSCS2_curated_complete.csv"))
all_vars_ctDNA = sample_info.merge(all_vars_ctDNA, how = "inner")

CH_fill_color = 
CH_line_color = 
normal_fill_color = "limegreen"
normal_line_color = "darkgreen"
ctDNA_fill_color = "cornflowerblue"
ctDNA_line_color = "mediumblue"

# Compare mutant vs WT fragments for each mutated locus
# CH vs WT fragments
path_to_df = os.path.join(DIR_fragment_counts, "CH.csv")
df_ch = MWU_mutated_to_WT_fragments(path_to_df, mutated_col_name = "CH fragments", WT_col_name = "WT fragments") # 18 sig
df_ch['CH fragments'] = df_ch['CH fragments'].apply(ast.literal_eval)
df_ch['WT fragments'] = df_ch['WT fragments'].apply(ast.literal_eval)
grouped_ch = df_ch.groupby("Patient_id")[["WT fragments", "CH fragments"]].agg(list).reset_index()
grouped_ch["WT fragments"] = grouped_ch["WT fragments"].apply(lambda x: [item for sublist in x for item in sublist])
grouped_ch["CH fragments"] = grouped_ch["CH fragments"].apply(lambda x: [item for sublist in x for item in sublist])
grouped_ch["Median WT CH fragments"] = grouped_ch["WT fragments"].apply(np.median)
grouped_ch["Median CH fragments"] = grouped_ch["CH fragments"].apply(np.median) # this one has median CH and WT values per patient

# ctDNA vs WT fragments
path_to_df = os.path.join(DIR_fragment_counts, "ctDNA.csv")
df_ctDNA = MWU_mutated_to_WT_fragments(path_to_df, mutated_col_name = "ctDNA fragments", WT_col_name = "WT fragments") # 189 sig
df_ctDNA['ctDNA fragments'] = df_ctDNA['ctDNA fragments'].apply(ast.literal_eval)
df_ctDNA['WT fragments'] = df_ctDNA['WT fragments'].apply(ast.literal_eval)
grouped_ctdna = df_ctDNA.groupby("Patient_id")[["WT fragments", "ctDNA fragments"]].agg(list).reset_index()
grouped_ctdna["WT fragments"] = grouped_ctdna["WT fragments"].apply(lambda x: [item for sublist in x for item in sublist])
grouped_ctdna["ctDNA fragments"] = grouped_ctdna["ctDNA fragments"].apply(lambda x: [item for sublist in x for item in sublist])
grouped_ctdna["Median WT ctDNA fragments"] = grouped_ctdna["WT fragments"].apply(np.median)
grouped_ctdna["Median ctDNA fragments"] = grouped_ctdna["ctDNA fragments"].apply(np.median) # this one has median CH and WT values per patient

patient_df = grouped_ctdna[["Patient_id", "Median WT ctDNA fragments", "Median ctDNA fragments"]].merge(grouped_ch[["Patient_id", "Median WT CH fragments", "Median CH fragments"]], how = "outer")
patient_df = patient_df.drop_duplicates()
patient_df = patient_df.merge(sample_info[["Patient_id", "Diagnosis"]].drop_duplicates(), how = "inner") # Add the diagnosis, we will make separate box plots for kidney and bladder
patient_df.columns = ["Patient_id", "ctDNA WT", "ctDNA mutant", "CH WT", "CH mutant", "Diagnosis"]
melted_df = patient_df.melt(id_vars=["Patient_id", "Diagnosis"], var_name="Fragment Type", value_name="Fragment Count")

# Boxplot. showing kidney and bladder together. 
save_path = "/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/fragmentomics/patient_box_plots_with_jitter_no_outliers.pdf"
create_box_plot_with_jitter(melted_df, save_path)

# Boxplot. Bladder and kidney shown in separate subplots.
melted_df = melted_df.dropna(subset=['Fragment Count'])
fig = plt.figure(figsize=(7, 3))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])
ax1 = fig.add_subplot(gs[0]) # Bladder subplot
ax1 = create_box_plot_with_jitter(melted_df, 'Bladder', ax=ax1)
ax2 = fig.add_subplot(gs[1]) # Kidney subplot
ax2 = create_box_plot_with_jitter(melted_df, 'Kidney', ax=ax2)
fig.tight_layout()
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/fragmentomics/combined_box_plots_with_jitter.pdf")
plt.close(fig)

# Scatter. Does median fragment length correlate with the number of mutations a patient may have? For both ctDNA muts and CH muts.
# Create a figure and axis
fig, ax = plt.subplots(figsize=(8, 6))
ax = scatter_n_muts_and_fragment_length(melted_df, muts_df = all_vars_chip, type_of_muts = "CH", ax = ax)
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/fragmentomics/n_muts_and_fragment_length.png")

def generate_histograms(df, diagnosis, title, fill_color, trace_color, save_path):
    if "CH fragments" in df.columns: 
        col_to_plot = "CH fragments"
    elif "ctDNA fragments" in df.columns:
        col_to_plot = "ctDNA fragments"
    
    if diagnosis == "both": 
        ch_del = list(chain.from_iterable(df[df["Type"] == "Deletion"][col_to_plot]))
        ch_ins = list(chain.from_iterable(df[df["Type"] == "Insertion"][col_to_plot]))
        ch_snv = list(chain.from_iterable(df[df["Type"] == "SNV"][col_to_plot]))
    else:
        ch_del = list(chain.from_iterable(df[(df["Type"] == "Deletion") & (df["Diagnosis"] == diagnosis)][col_to_plot]))
        ch_ins = list(chain.from_iterable(df[(df["Type"] == "Insertion") & (df["Diagnosis"] == diagnosis)][col_to_plot]))
        ch_snv = list(chain.from_iterable(df[(df["Type"] == "SNV") & (df["Diagnosis"] == diagnosis)][col_to_plot]))
    all_ch = ch_del + ch_ins + ch_snv
    
    fig = plt.figure(figsize=(8, 6))
    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])
    ax1 = fig.add_subplot(gs[0, 0])
    ax1 = make_histogram(ch_del, "Deletions", color=fill_color, color_kde=trace_color, ax=ax1)
    ax2 = fig.add_subplot(gs[0, 1])
    ax2 = make_histogram(ch_ins, "Insertions", color=fill_color, color_kde=trace_color, ax=ax2)
    ax3 = fig.add_subplot(gs[1, 0])
    ax3 = make_histogram(ch_snv, "SNVs", color=fill_color, color_kde=trace_color, ax=ax3)
    ax4 = fig.add_subplot(gs[1, 1])
    ax4 = make_histogram(all_ch, "All events", color=fill_color, color_kde=trace_color, ax=ax4)
    
    fig.suptitle(title)
    gs.tight_layout(fig)
    fig.savefig(save_path)
    plt.close(fig)


















def plot_patient_jitter_box(patient_df, save_path):
    # Set up the figure and axes
    fig, ax = plt.subplots(figsize=(10, 8))
    # Create a box plot with jitter using Seaborn
    sns.boxplot(data=patient_df, x="Patient_id", y="Median WT ctDNA fragments", hue="Median ctDNA fragments", dodge=True, linewidth=1.5, ax=ax)
    sns.stripplot(data=patient_df, x="Patient_id", y="Median WT ctDNA fragments", hue="Median ctDNA fragments", dodge=True, size=8, alpha=0.5, ax=ax)
    # Set labels and title
    ax.set_xlabel("Patient ID")
    ax.set_ylabel("Median WT ctDNA Fragments")
    ax.set_title("Patient Fragmentomics Analysis")
    # Add legend
    ax.legend(title="Median ctDNA Fragments", loc="upper right")
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha="right")
    # Add vertical lines for median values of each group
    medians = patient_df[["Median WT ctDNA fragments", "Median ctDNA fragments", "Median WT CH fragments", "Median CH fragments"]].median()
    for median_val in medians:
        ax.axhline(y=median_val, color="r", linestyle="--", linewidth=2)
    # Save the plot
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()

# Define the save path for the plot
save_path = "/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/fragmentomics/patient_jitter_box_plot_ax.png"

# Call the function to create and save the jitter box plot with median lines using ax method
plot_patient_jitter_box(patient_df, save_path)

###################################################
# CH focused histograms
###################################################
CH_mutant_fragments = list(chain.from_iterable(df_ch['CH fragments']))
CH_WT_fragments = list(chain.from_iterable(df_ch['WT fragments']))
ctDNA_mutant_fragments = list(chain.from_iterable(df_ctDNA['ctDNA fragments']))
ctDNA_WT_fragments = list(chain.from_iterable(df_ctDNA['WT fragments']))

# Parameters for each plot
params = [
    {
        "diagnosis": "both",
        "title": "Kidney and Bladder",
        "fill_color": "lightcoral",
        "trace_color": "brown",
        "save_path": os.path.join(DIR_figures, "CH_histograms.png"),
        "df": df_ch
    },
    {
        "diagnosis": "Kidney",
        "title": "Kidney",
        "fill_color": (239/255, 169/255, 143/255, 255/255),
        "trace_color": "orangered",
        "save_path": os.path.join(DIR_figures, "CH_histograms_kidney.png"),
        "df": df_ch
    },
    {
        "diagnosis": "Bladder",
        "title": "Bladder",
        "fill_color": (143/255, 215/255, 239/255, 255/255),
        "trace_color": "deepskyblue",
        "save_path": os.path.join(DIR_figures, "CH_histograms_bladder.png"),
        "df": df_ch
    }
]

# Generate histograms for each set of parameters
for param in params:
    generate_histograms(param["df"], param["diagnosis"], param["title"], param["fill_color"], param["trace_color"], param["save_path"])

# Parameters for each plot
params_ctDNA = [
    {
        "diagnosis": "both",
        "title": "Kidney and Bladder",
        "fill_color": "lightcoral",
        "trace_color": "brown",
        "save_path": os.path.join(DIR_figures, "ctDNA_histograms.png"),
        "df": df_ctDNA
    },
    {
        "diagnosis": "Kidney",
        "title": "Kidney",
        "fill_color": (239/255, 169/255, 143/255, 255/255),
        "trace_color": "orangered",
        "save_path": os.path.join(DIR_figures, "ctDNA_histograms_kidney.png"),
        "df": df_ctDNA
    },
    {
        "diagnosis": "Bladder",
        "title": "Bladder",
        "fill_color": (143/255, 215/255, 239/255, 255/255),
        "trace_color": "deepskyblue",
        "save_path": os.path.join(DIR_figures, "ctDNA_histograms_bladder.png"),
        "df": df_ctDNA
    }
]

# Generate histograms for each set of parameters
for param in params_ctDNA:
    generate_histograms(param["df"], param["diagnosis"], param["title"], param["fill_color"], param["trace_color"], param["save_path"])

# Now make a histogram of healthy normals
vip_normals = [os.path.join(DIR_insertsizes, file) for file in os.listdir(DIR_insertsizes) if "VIP" in file]
values = []
for file in vip_normals:
    df = pd.read_csv(file, skiprows=10, sep = "\t")
    repeated_values = df.loc[df.index.repeat(df['All_Reads.fr_count']), 'insert_size'].tolist()
    values.extend(repeated_values)

values_downsampled = random.sample(values, 50000)
fig, ax = plt.subplots(figsize=(4, 4))
ax = make_histogram(values_downsampled, "Control samples", color = , color_kde = normal_line_color, ax = ax, bw_adjust = 0.2) # density plot from just one sample
fig.tight_layout()
fig.savefig(os.path.join(DIR_figures, "normal_fragments_VIP_downsamples.png"))

# Showing all CH, all ctDNA and controls together in a gridspec setting.
fig = plt.figure(figsize=(7, 3))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])
ax1 = fig.add_subplot(gs[0]) # Kidney
ax2 = fig.add_subplot(gs[1]) # Bladder

bladder_ch = list(chain.from_iterable(df_ch[df_ch["Diagnosis"] == "Bladder"]["CH fragments"]))
bladder_ctDNA = list(chain.from_iterable(df_ctDNA[df_ctDNA["Diagnosis"] == "Bladder"]["ctDNA fragments"]))
ax1 = make_histogram(bladder_ch, title = "", color = "lightcoral", color_kde = "brown", ax = ax1, annotate_median = False, alpha = 0.75, bw_adjust= 0.05, show_hist = False)
ax1 = make_histogram(bladder_ctDNA, title = "", color = "cornflowerblue", color_kde = "mediumblue", ax = ax1, annotate_median = False, alpha = 0.75, bw_adjust= 0.05, show_hist = False)
ax1 = make_histogram(values_downsampled, title = "", color = normal_fill_color, color_kde = normal_line_color, ax = ax1, annotate_median = False, alpha = 0.75, bw_adjust= 0.05, show_hist = False)
ax1.set_title("Bladder")

kidney_ch = list(chain.from_iterable(df_ch[df_ch["Diagnosis"] == "Kidney"]["CH fragments"]))
kidney_ctDNA = list(chain.from_iterable(df_ctDNA[df_ctDNA["Diagnosis"] == "Kidney"]["ctDNA fragments"]))
ax2 = make_histogram(kidney_ch, title = "", color = "lightcoral", color_kde = "brown", ax = ax2, annotate_median = False, alpha = 0.75, bw_adjust= 0.05, show_hist = False)
ax2 = make_histogram(kidney_ctDNA, title = "", color = "cornflowerblue", color_kde = "mediumblue", ax = ax2, annotate_median = False, alpha = 0.75, bw_adjust= 0.05, show_hist = False)
ax2 = make_histogram(values_downsampled, title = "", color = normal_fill_color, color_kde = normal_line_color, ax = ax2, annotate_median = False, alpha = 0.75, bw_adjust= 0.05, show_hist = False)
ax2.set_title("Kidney")

legend_colors = ["brown", "mediumblue", normal_fill_color]
legend_labels = ["CH", "ctDNA", "Control"]
legend_handles = [plt.Line2D([0], [0], color=color, label=label, linewidth=2) for color, label in zip(legend_colors, legend_labels)]
ax1.legend(handles=legend_handles, loc="upper right", frameon=False)

gs.tight_layout(fig)
fig.savefig(os.path.join(DIR_figures, "superimposed_kdes.pdf"))