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

DIR_fragment_counts = "/groups/wyattgrp/users/amunzur/pipeline/results/fragmentomics"
DIR_figures = "/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/fragmentomics"
DIR_insertsizes = "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/insert_size/SSCS2"
PATH_sample_information = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"
PATH_vars_ch = ""

sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])

# Load the mut dfs
all_vars_chip = pd.read_csv(os.path.join("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/chip.csv"))
all_vars_chip = sample_info.merge(all_vars_chip, how = "inner")
all_vars_ctDNA = pd.read_csv(os.path.join("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/somatic.csv"))
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

values_downsampled = random.sample(values, 80000)
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