import pandas as pd
import numpy as np
import os

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec



"""
Collect sequencing stats (Median target coverage, On target rate, Duplicate fraction, Number of reads) in one file. 
"""

DIR_HS_metrics = "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/PICARD_HS_metrics/SSCS2"
DIR_umi_metrics = "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/umi_metrics" # SSCS1 here, although doesn't really matter
PATH_averaged_depths_SSCS2 = "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/averaged_depth/SSCS2/averaged_depths.txt"
DIR_raw_read_counts = "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/read_counts/merged" # from raw fastq, no trimmming
PATH_SAMPLE_INFO = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"
PATH_output = "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/seq_quality_metrics/seq_quality_metrics.tsv"

def return_duplicate_frac(file_path):
    """
    Given a file path to the UMI metrics, return the duplicate perc.
    """
    df = pd.read_csv(file_path, sep = "\t")
    duplicate_frac = df[df["family_size"] == 2]["fraction_gt_or_eq_family_size"].iloc[0]
    sample_name = os.path.basename(file_path).split(".")[0]
    result_dict = {sample_name: duplicate_frac}
    return(result_dict)

def return_HS_metrics(file_path): 
    """
    From Picard's HS metrics, return the following in a list: 
    """
    df = pd.read_csv(os.path.join(DIR_HS_metrics, file_path), sep = "\t", skiprows = 6).iloc[0, :]
    # MEDIAN_TARGET_COVERAGE = df["MEDIAN_TARGET_COVERAGE"]
    ON_TARGET_PERCENTAGE = df["PCT_SELECTED_BASES"]
    sample_name = os.path.basename(file_path).split(".")[0]
    result_dict = {sample_name: ON_TARGET_PERCENTAGE}
    return(result_dict)

def return_raw_read_counts(file_path):
    """
    Raw read counts from FastQ.
    """
    df = pd.read_csv(file_path, sep = "\t", names = ["sample", "value"])
    result_dict = {df.iloc[0, ]["sample"]: df.iloc[0, ]["value"]}
    return(result_dict)

DIR_HS_metrics_list = [os.path.join(DIR_HS_metrics, file) for file in os.listdir(DIR_HS_metrics)]
DIR_umi_metrics_list = [os.path.join(DIR_umi_metrics, file) for file in os.listdir(DIR_umi_metrics)]
DIR_raw_read_counts_list = [os.path.join(DIR_raw_read_counts, file) for file in os.listdir(DIR_raw_read_counts)]

# ON TARGET PERCENTAGE AND THE MEDIAN TARGET COVERAGE
hs_dict_list = []
for file_path in DIR_HS_metrics_list: 
    result = return_HS_metrics(file_path)
    hs_dict_list.append(result)

df_hs = pd.melt(pd.DataFrame(hs_dict_list), var_name='Sample_name', value_name='On target rate')
df_hs = df_hs[~pd.isnull(df_hs["On target rate"])].reset_index(drop = True)

# DUPLICATE PERCENTAGE
umi_dict_list = []
for file_path in DIR_umi_metrics_list: 
    result = return_duplicate_frac(file_path)
    umi_dict_list.append(result)

df_dup = pd.melt(pd.DataFrame(umi_dict_list), var_name='Sample_name', value_name='Duplicate fraction')
df_dup = df_dup[~pd.isna(df_dup["Duplicate fraction"])].reset_index(drop = True)

# NUMBER OF READS
reads_dict_list = []
for file_path in DIR_raw_read_counts_list: 
    result = return_raw_read_counts(file_path)
    reads_dict_list.append(result)

df_reads = pd.melt(pd.DataFrame(reads_dict_list), var_name='Sample_name', value_name='Number of reads')
df_reads = df_reads[~pd.isna(df_reads["Number of reads"])]
df_reads["Sample_name"] = df_reads["Sample_name"].str.replace("GUBB-", "GU-")

# MEDIAN DEPTH IN TARGET REGIONS, no need for function.
df_depth_SSCS2 = pd.read_csv(PATH_averaged_depths_SSCS2, sep = " ", names = ["Sample_name", "Mean depth", "Median depth SSCS2"])[["Sample_name", "Median depth SSCS2"]]

# COMBINE ALL AND ADD SAMPLE INFORMATION
df = df_hs.merge(df_depth_SSCS2, how = "inner").merge(df_reads, how = "inner").merge(df_dup, how = "inner")
df['Patient'] = df['Sample_name'].str.split('_', expand=True)[0]
df["Patient"] = df["Patient"].str.replace("GU-", "")
df['Date'] = df['Sample_name'].str.split('[-_]', expand=True)[6]

df_sample_info = pd.read_csv(PATH_SAMPLE_INFO, sep = "\t", names = ["Patient", "Date", "Diagnosis", "Timepoint"])
sequencing_metrics = df.merge(df_sample_info, how = "left")
sequencing_metrics["Number of reads (millions)"] = (sequencing_metrics["Number of reads"]/1e6)*2 # times 2 because R1 and R2
sequencing_metrics = sequencing_metrics[["Patient", "Sample_name", "Diagnosis", "Date", "Timepoint", "Median depth SSCS2", "On target rate", "Duplicate fraction", "Number of reads (millions)"]]
sequencing_metrics = sequencing_metrics.sort_values(by = ["Diagnosis", "Patient"])
sequencing_metrics.to_csv(PATH_output, sep = "\t", index = False)

##########################################################
# VISUALIZING THE SEQUENCING METRICS
##########################################################

PATH_seq_metrics = "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/seq_quality_metrics/seq_quality_metrics.tsv"
DIR_figures = "/groups/wyattgrp/users/amunzur/pipeline/results/figures/seq_metrics"

def enumerate_col(df, old_col, new_col):
    df[new_col] = df[old_col].map({gene: number + 1 for number, gene in enumerate(df[old_col].unique())})
    return(df)

def remove_frame(ax, list_to_remove):
    for i in list_to_remove:
        ax.spines[i].set_visible(False)
    return(ax)

def return_sequencing_metrics(df, consensus): 
    """
    This function prints some sequencing metrics. Useful when making presentations.
    """
    cfDNA_data = df[df['Sample_type'] == 'cfDNA']
    WBC_data = df[df['Sample_type'] == 'WBC']
    # stats
    median_cfDNA_depth = cfDNA_data[f'Median depth {consensus}'].median()
    median_WBC_depth = WBC_data[f'Median depth {consensus}'].median()
    min_cfDNA_depth = cfDNA_data[f'Median depth {consensus}'].min()
    min_WBC_depth = WBC_data[f'Median depth {consensus}'].min()
    max_cfDNA_depth = cfDNA_data[f'Median depth {consensus}'].max()
    max_WBC_depth = WBC_data[f'Median depth {consensus}'].max()
    s = f"{consensus}: Median cfDNA depth: {median_cfDNA_depth}, Median WBC depth: {median_WBC_depth}, min_cfDNA_depth: {min_cfDNA_depth}, min_WBC_depth: {min_WBC_depth}, max_cfDNA_depth: {max_cfDNA_depth}, max_WBC_depth: {max_WBC_depth}"
    return(s)

def make_WBC_cfDNA_depth_plot(df, ax, colname_to_plot, legend = False):
    """
    Makes a bar chart indicating both the cfDNA and WBC depths and saves the plot.
    """
    # generate separate cfDNA and WBC dfs
    cfDNA_data = df[df['Sample_type'] == 'cfDNA'].sort_values(colname_to_plot).reset_index(drop = True)
    cfDNA_data = enumerate_col(cfDNA_data, "Patient_time_point_identifier", "Patient_time_point_enumerated") # helps with plotting
    WBC_data = df[df['Sample_type'] == 'WBC'].merge(cfDNA_data[["Patient_time_point_identifier", "Patient_time_point_enumerated"]], how = "left")
    
    ax.bar(cfDNA_data['Patient_time_point_enumerated'], cfDNA_data[colname_to_plot], color='blue', label='cfDNA')
    ax.bar(WBC_data['Patient_time_point_enumerated'], -WBC_data[colname_to_plot], color='orange', label='WBC')
    ax.set_ylim((-7000, 5000))
    ax.set_yticks(range(-7000, 6000, 1000))
    ax.tick_params(axis='y', labelsize=12)
    ax.set_yticklabels([abs(label) for label in ax.get_yticks()])
    ax = remove_frame(ax, ["top", "right", "bottom"])
    ax.set_xticks([])   
    
    ax.axhline(y=cfDNA_data[colname_to_plot].median(), linestyle='--', color='black', label=f"Median")
    ax.axhline(y=-WBC_data[colname_to_plot].median(), linestyle='--', color='black')   
    
    # Customize the plot
    ax.set_xlabel('cfDNA/WBC pair')
    ax.set_ylabel(colname_to_plot)
    ax.set_title(colname_to_plot)
    if legend:
        ax.legend()
    ax.grid(True, axis='y')
    return ax


def make_WBC_cfDNA_depth_plot_DIFFERENCE(df, DIR_figures, plot_name):
    """
    Makes a bar chart indicating both the cfDNA and WBC depths and saves the plot.
    """
    # Generate separate cfDNA and WBC dfs
    cfDNA_data = df[df['Sample_type'] == 'cfDNA']
    cfDNA_data['Depth_difference'] = cfDNA_data['Median depth SSCS1'] - cfDNA_data['Median depth SSCS2']
    cfDNA_data = cfDNA_data.sort_values("Depth_difference").reset_index(drop = True)
    cfDNA_data = enumerate_col(cfDNA_data, "Patient_time_point_identifier", "Patient_time_point_enumerated") # helps with plotting
    WBC_data = df[df['Sample_type'] == 'WBC'].merge(cfDNA_data[["Patient_time_point_identifier", "Patient_time_point_enumerated"]], how = "left")
    
    diagnosis = cfDNA_data["Diagnosis"].unique()[0]
    # Calculate the difference between Median depth SSCS1 and Median depth SSCS2
    
    WBC_data['Depth_difference'] = WBC_data['Median depth SSCS1'] - WBC_data['Median depth SSCS2']
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot the difference for cfDNA on the positive y-axis
    ax.bar(cfDNA_data['Patient_time_point_enumerated'], cfDNA_data['Depth_difference'], color='blue', label='cfDNA')
    # Plot the difference for WBC on the negative y-axis
    ax.bar(WBC_data['Patient_time_point_enumerated'], -WBC_data['Depth_difference'], color='orange', label='WBC')
    ax.set_ylim((-1500, 1500))
    ax.set_yticks(range(-1500, 1500, 500))
    ax.set_yticklabels([abs(label) for label in ax.get_yticks()])
    ax = remove_frame(ax, ["top", "right", "bottom"])
    ax.set_xticks([])
    # Plot median lines for cfDNA and WBC
    ax.axhline(y=cfDNA_data['Depth_difference'].median(), linestyle='--', color='black', label=f"Median cfDNA")
    ax.axhline(y=-WBC_data['Depth_difference'].median(), linestyle='--', color='black', label=f"Median WBC")   
    # Customize the plot
    ax.set_xlabel('cfDNA/WBC pair')
    ax.set_ylabel('Depth Difference (SSCS1 - SSCS2)')
    ax.set_title(f'Depth Difference (SSCS1 - SSCS2) for cfDNA and WBC {diagnosis} Samples')
    ax.legend()
    ax.grid(True, axis='y')
    # Save the plot as both PNG and PDF
    fig.savefig(os.path.join(DIR_figures, plot_name))
    fig.savefig(os.path.join(DIR_figures, plot_name.replace(".png", ".pdf")))

#############################################
# BLADDER
df_main = pd.read_csv(PATH_seq_metrics, sep = "\t")
df_bladder = (
    df_main
    .loc[df_main['Diagnosis'] == "Bladder"]
    .assign(Sample_type = df_main["Sample_name"].str.extract("(WBC|cfDNA)"), 
            Patient_time_point = df_main["Sample_name"].str.extract("(Baseline|OnTreatment)"), 
            Patient_id = df_main["Sample_name"].str.replace("(_WBC.*|_cfDNA.*)", ""))
)
df_bladder["Patient_time_point_identifier"] = df_bladder["Sample_name"].str.replace("(_WBC*|_cfDNA*)", "")

fig, ax = plt.subplots(figsize=(8, 4))
ax = make_WBC_cfDNA_depth_plot(df_bladder, ax, "Median depth SSCS2", legend = True)
ax.set_title("Median depth SSCS2 in bladder samples")
fig.savefig(os.path.join(DIR_figures, "bladder_consensus_depth.png"))
fig.savefig(os.path.join(DIR_figures, "bladder_consensus_depth.pdf"))

# make_WBC_cfDNA_depth_plot_DIFFERENCE(df_bladder, DIR_figures, "SSCS1 and SSCS2 depth difference bladder.png")
# return_sequencing_metrics(df_bladder, "SSCS1")
# return_sequencing_metrics(df_bladder, "SSCS2")
#############################################

#############################################
# KIDNEY
df_kidney = (
    df_main
    .loc[df_main['Diagnosis'] == "Kidney"]
    .assign(Sample_type = df_main["Sample_name"].str.extract("(WBC|cfDNA)"), 
            Patient_time_point = df_main["Sample_name"].str.extract("(Baseline|OnTreatment)"), 
            Patient_id = df_main["Sample_name"].str.replace("(_WBC.*|_cfDNA.*)", ""))
)
df_kidney["Patient_time_point_identifier"] = df_kidney["Sample_name"].str.replace("(_WBC*|_cfDNA*)", "")
fig, ax = plt.subplots(figsize=(8, 4))

ax = make_WBC_cfDNA_depth_plot(df_kidney, ax, "Median depth SSCS2", legend = True)
ax.set_title("Median depth SSCS2 in kidney samples")
fig.savefig(os.path.join(DIR_figures, "kidney_consensus_depth.png"))
fig.savefig(os.path.join(DIR_figures, "kidney_consensus_depth.pdf"))

# make_WBC_cfDNA_depth_plot_DIFFERENCE(df_kidney, DIR_figures, "SSCS1 and SSCS2 depth difference kidney.png")
# return_sequencing_metrics(df_kidney, "SSCS1")
# return_sequencing_metrics(df_kidney, "SSCS2")
#############################################