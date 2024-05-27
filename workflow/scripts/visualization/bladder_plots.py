import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import os
import seaborn as sns
import numpy as np
import scipy.stats as stats

PATH_bladder_ct_fractions = "/groups/wyattgrp/users/amunzur/pipeline/resources/ct_fractions/bladder_ct_fractions.tsv"
PATH_mutations = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/somatic_curated.csv"
PATH_mutations_SSCS1 = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/somatic_SSCS1.csv"

DIR_figures = "/groups/wyattgrp/users/amunzur/pipeline/results/figures/somatic"

curated = pd.read_csv(PATH_mutations)[['Patient_id', 'Sample_name_t_x', "Gene_x", "Chrom_x", "Position"]]
curated.columns = ['Patient_id', 'Sample_name_t', "Gene", "Chrom", "Position"]
sscs1_raw = pd.read_csv(PATH_mutations_SSCS1)

# Perform anti-join
sscs1 = sscs1_raw.merge(curated, on=['Patient_id', 'Sample_name_t', 'Gene', 'Chrom', 'Position'], how='left', indicator=True)
sscs1 = sscs1[sscs1['_merge'] == 'both'].drop('_merge', axis=1)

kidney = sscs1[sscs1["Diagnosis"] == "Kidney"].reset_index(drop = True)
bladder = sscs1[sscs1["Diagnosis"] == "Bladder"].reset_index(drop = True)

def get_gene_frequency(df):
    """
    Given a df with mutation calls, get the frequency each gene.
    """
    gene_counts = df['Gene'].value_counts()
    total_rows = len(df)
    gene_fraction = gene_counts / total_rows
    gene_frequency_df = pd.DataFrame({'Gene': gene_counts.index, 'Frequency': gene_fraction})
    gene_frequency_df.loc[gene_frequency_df["Gene"].str.contains("U2AF1"), "Gene"] = "U2AF1"
    gene_frequency_df = gene_frequency_df.reset_index(drop = True)
    
    return(gene_frequency_df)

def make_VAF_plot(df, path_figure, p_title, x_column, y_column, x_title, y_title):
    """
    Given a df, make a scatter plot correlating the cfDNA and WBC VAFs.
    """
    fig, ax = plt.subplots()
    ax.scatter(df[x_column]*100, df[y_column]*100, color="black")
    ax.set_xlim(-1, 60)  # Set x-axis limits from 0 to 10
    ax.set_ylim(-1, 60)  # Set y-axis limits from -5 to 5
    ax.set_xlabel(x_title)
    ax.set_ylabel(y_title)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title(p_title)
    
    # Add correlation line
    x = df[x_column]*100
    y = df[y_column]*100
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    ax.plot(x, slope*x + intercept, color='black')
    
    # Add R-squared value
    r_squared = r_value**2
    ax.text(0.5, 0.9, f"R-squared = {r_squared:.2f}", transform=ax.transAxes)
    
    plt.tight_layout()
    plt.savefig(path_figure)

kidney_gene_freq = get_gene_frequency(kidney)
kidney_gene_freq["Diagnosis"] = "Kidney"

bladder_gene_freq = get_gene_frequency(bladder)
bladder_gene_freq["Diagnosis"] = "Bladder"

# A dodge barplot, colored according to diagnosis
merged_df = pd.concat([kidney_gene_freq, bladder_gene_freq])
sns.catplot(x = "Gene", y = 'Frequency', hue = "Diagnosis", kind='bar', data = merged_df, legend = True)
ax = plt.gca()
fig = plt.gcf()
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_title("Frequency of genes with somatic mutations in kidney and bladder")
fig.tight_layout()
plt.savefig(os.path.join(DIR_figures, "Gene_frequency_dodge.png"))

# cfDNA vs WBC VAF 
make_VAF_plot(kidney, os.path.join(DIR_figures, "kidney_VAF_baselines.png"), "Kidney baseline samples", "VAF_t", "VAF_n", "cfDNA VAF (%)", "WBC VAF (%)")
make_VAF_plot(bladder, os.path.join(DIR_figures, "bladder_VAF_baselines.png"), "Bladder baseline samples", "VAF_t", "VAF_n", "cfDNA VAF (%)", "WBC VAF (%)")
make_VAF_plot(bladder, os.path.join(DIR_figures, "bladder_VAF_baselines_ct_fraction.png"), "Bladder baseline samples", "ct_fraction", "VAF_t", "Tumor (%)", "WBC VAF (%)")







def plot_mutation_bubble_chart(df, dir_figures, filename="mutation_bubble_chart.png"):
    """
    Generate a mutation bubble chart based on the provided DataFrame.
    Parameters:
    - df: DataFrame containing mutation data.
    - dir_figures: Directory to save the plot.
    - filename: Name of the file to save the plot. Default is "mutation_bubble_chart.png".
    """
    subset_columns = ['Patient_id', 'Sample_name_t', 'Gene', 'VAF_t']
    subset_df = df[subset_columns]
    heatmap_data = subset_df.pivot_table(values='VAF_t', index='Gene', columns='Sample_name_t', aggfunc='max', fill_value=0)
    plt.figure(figsize=(12, 8))
    for gene in heatmap_data.index:
        for sample in heatmap_data.columns:
            vaf_value = heatmap_data.loc[gene, sample]
            if vaf_value > 0:
                # Adjust the size of the bubbles based on the VAF_t value
                bubble_size = vaf_value * 100
                plt.scatter(sample, gene, s=bubble_size, alpha=0.7, label=f'{gene} ({sample})', edgecolors='black')
    #
    plt.xlabel('Sample Name')
    plt.ylabel('Gene')
    plt.title('Mutation Bubble Chart')
    plt.legend()
    # Save the plot
    file_path = os.path.join(dir_figures, filename)
    plt.savefig(file_path)
    # Show the plot
    # Print the path where the plot is saved
    print(f"Plot saved at: {file_path}")

