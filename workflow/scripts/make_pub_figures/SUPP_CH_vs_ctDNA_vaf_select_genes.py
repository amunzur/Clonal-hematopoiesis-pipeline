"""
For select genes compares the VAF of ctDNA vs CH mutations.
"""
import pandas as pd
import numpy as np
import os 
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import matplotlib.gridspec as gridspec


DIR_working = "/groups/wyattgrp/users/amunzur/pipeline"
figure_dir = os.path.join(DIR_working, "results/figures/pub_figures")
source_functions = os.path.join(DIR_working, "workflow/scripts/visualization/UTILITIES_make_chip_plots.py")
ctDNA_fractions_path = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/ctDNA_fractions.csv"
PATH_sample_information = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"

with open(source_functions, 'r') as file:
    script_code = file.read()

exec(script_code)

# LOAD CHIP DATASETS
select_genes_list = ["CHEK2", "ATM", "TP53", "ASXL1", "BRCA1", "BRCA2", "ERBB2", "ARID1A", "KMT2D", "DNMT3A", "SETD2", "PBRM1"]
all_vars_chip = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/chip.csv"))
all_vars_somatic = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/somatic.csv"))

baseline = all_vars_chip[
    (all_vars_chip["Dependent"] == False) & 
    (all_vars_chip["Timepoint"] == "Baseline")]

baseline_somatic = all_vars_somatic[
    (all_vars_somatic["Dependent"] == False) & 
    (all_vars_somatic["Timepoint"] == "Baseline")]

fig = plt.figure(figsize=(8, 8))
gs = gridspec.GridSpec(4, 3, height_ratios=[1, 1, 1, 1], width_ratios = [1, 1, 1], hspace = 0.6, wspace = 1)

for i, gene in enumerate(select_genes_list):
    print(f"Working on {gene}")
    row = i // 3  # Determine the row index
    col = i % 3   # Determine the column index
    ax = fig.add_subplot(gs[row, col])
    
    # Subset to the gene of interest
    gene_df_ch = baseline[baseline["Gene"] == gene]
    gene_df_somatic = baseline_somatic[baseline_somatic["Gene"] == gene]
    
    # Plotting    
    if i == 0:
        show_legend = True
    else:
        show_legend = False
        
    ax = plot_ctDNA_and_CH_VAF_box(
        ch_df = gene_df_ch, 
        ctdna_df = gene_df_somatic, 
        PATH_sample_information = PATH_sample_information, 
        ax = ax, 
        ax_title = gene,
        hide_xticks = False, 
        all_ctdna_muts_df = baseline_somatic, 
        scatter_zorder = 5,
        scatter_alpha = 1,
        scatter_size = 4, 
        show_legend = show_legend)
    
    if row != 3:
        ax.set_xlabel("")
    
    if col != 0:
        ax.set_ylabel("")


gs.tight_layout(fig)
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/SUPP_CH_vs_ctDNA_vaf_select_genes.png")
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/SUPP_CH_vs_ctDNA_vaf_select_genes.pdf")
