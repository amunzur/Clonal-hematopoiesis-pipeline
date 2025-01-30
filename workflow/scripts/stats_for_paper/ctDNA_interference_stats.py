"""
All code to generate stats relevant for the ctDNA part of the MS.
"""

import pandas as pd
import numpy as np
import os 
import re
import matplotlib as mpl
mpl.use('Agg')
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from scipy.stats import spearmanr
from scipy.stats import chi2_contingency
from scipy.stats import fisher_exact

def make_ch_ctDNA_prevalence_df(baseline_ctDNA, baseline_chip, min_vaf):
    """
    Generates the table to make figure 2A. 
    """
    baseline_ctDNA_gene = baseline_ctDNA[["Gene"]].sort_values("Gene").value_counts().reset_index()
    baseline_ctDNA_gene.columns = ["Gene", "Counts_ctDNA"]
    
    baseline_chip = baseline_chip[baseline_chip["VAF_n"] > min_vaf]
    baseline_chip_gene = baseline_chip[["Gene"]].sort_values("Gene").value_counts().reset_index()
    baseline_chip_gene.columns = ["Gene", "Counts_ch"]
    
    combined = baseline_chip_gene.merge(baseline_ctDNA_gene, on = ["Gene"], how = "outer").fillna(0)
    combined["summed"] = combined["Counts_ch"] + combined["Counts_ctDNA"]
    combined["perc_ch"] = (combined["Counts_ch"]/combined["summed"])*100
    combined["perc_ctDNA"] = (combined["Counts_ctDNA"]/combined["summed"])*100
    combined = combined.sort_values(by = "perc_ch", ascending = False).reset_index(drop = True)
    
    return(combined)

source_functions = os.path.join("/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/visualization/UTILITIES_make_chip_plots.py")
with open(source_functions, 'r') as file:
    script_code = file.read()

exec(script_code)

PATH_sample_information = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"
PATH_ctDNA_fraction = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/ctDNA_fractions.csv"

# LOAD CHIP DATASETS
all_vars_chip = pd.read_csv(os.path.join("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/chip.csv"))
all_vars_chip = all_vars_chip[all_vars_chip["Dependent"] == False]
base_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
prog_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)
prog_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)

# LOAD SOMATIC DATASETS
all_vars_somatic = pd.read_csv(os.path.join("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/somatic.csv"))
all_vars_somatic = all_vars_somatic[~all_vars_somatic["Patient_id"].isin(["20-313", "21-184", "21-430"])] # exclude some samples due to oxidative damage
all_vars_somatic = all_vars_somatic[all_vars_somatic["Dependent"] == False]
base_kidney_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Kidney")].reset_index(drop = True)
prog_kidney_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "During treatment") & (all_vars_somatic["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Bladder")].reset_index(drop = True)
prog_bladder_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "During treatment") & (all_vars_somatic["Diagnosis"] == "Bladder")].reset_index(drop = True)

baseline_chip = pd.concat([base_kidney_chip, base_bladder_chip]).reset_index(drop = True)
baseline_ctDNA = pd.concat([base_kidney_somatic, base_bladder_somatic]).reset_index(drop = True)

###################################################################################
# 1. Genes mutated in ctDNA for bladder and kidney separately. Calculating # of the cohort mutated
sample_info = pd.read_csv(PATH_sample_information, sep="\t", names=["Patient_id", "Date", "Diagnosis", "Timepoint"])
n_bladder = sample_info[sample_info["Diagnosis"] == "Bladder"].drop_duplicates("Patient_id").shape[0]
n_kidney = sample_info[sample_info["Diagnosis"] == "Kidney"].drop_duplicates("Patient_id").shape[0]

# fraction mutated, using ctDNA+ as denominator, not the full cohort
n_bladder_ctDNA_pos = base_bladder_somatic["Patient_id"].unique().shape[0]
n_kidney_ctDNA_pos = base_kidney_somatic["Patient_id"].unique().shape[0]

# subset
bladder_gene_df = base_bladder_somatic[["Patient_id", "Gene"]].drop_duplicates()["Gene"].value_counts().reset_index().rename(columns = {"index": "Gene", "Gene": "Counts_bladder"})
kidney_gene_df = base_kidney_somatic[["Patient_id", "Gene"]].drop_duplicates()["Gene"].value_counts().reset_index().rename(columns = {"index": "Gene", "Gene": "Counts_kidney"})
bladder_gene_df["Cohort fraction"] = round((bladder_gene_df["Counts_bladder"]/n_bladder_ctDNA_pos)*100, 2)
kidney_gene_df["Cohort fraction"] = round((kidney_gene_df["Counts_kidney"]/n_kidney_ctDNA_pos)*100, 2)

###################################################################################
# 2. % of the mRCC and mUC cohorts where we detected ctDNA mutations.
ctDNA_kidney_annotated = annotate_mutation_status(base_kidney_somatic, "Kidney", PATH_sample_information, annotate_what = "ctDNA")
ctDNA_kidney_annotated = ctDNA_kidney_annotated[ctDNA_kidney_annotated["Timepoint"] == "Baseline"]
ctDNA_kidney_annotated["ctDNA status"].value_counts(normalize = True)
ctDNA_kidney_annotated["ctDNA status"].value_counts()

ctDNA_bladder_annotated = annotate_mutation_status(base_bladder_somatic, "Bladder", PATH_sample_information, annotate_what = "ctDNA")
ctDNA_bladder_annotated = ctDNA_bladder_annotated[ctDNA_bladder_annotated["Timepoint"] == "Baseline"]
ctDNA_bladder_annotated["ctDNA status"].value_counts(normalize = True)
ctDNA_bladder_annotated["ctDNA status"].value_counts()

# Showing the entire cohort
pd.concat([ctDNA_kidney_annotated, ctDNA_bladder_annotated])["ctDNA status"].value_counts(normalize = True)
pd.concat([ctDNA_kidney_annotated, ctDNA_bladder_annotated])["ctDNA status"].value_counts()

###################################################################################
# 3. VAF statistics on ctDNA mutations
round(np.median(baseline_ctDNA["VAF_t"]), 2) # median
[min(baseline_ctDNA["VAF_t"]), max(baseline_ctDNA["VAF_t"])] # range
q75, q25 = np.percentile(baseline_ctDNA["VAF_t"], [75 ,25])
iqr = q75 - q25

# 3. VAF statistics on ctDNA mutations - kidney
round(np.median(base_kidney_somatic["VAF_t"]), 2) # median
[round(min(base_kidney_somatic["VAF_t"]), 2), round(max(base_kidney_somatic["VAF_t"]), 2)] # range
q75, q25 = np.percentile(base_kidney_somatic["VAF_t"], [75 ,25])
q25, q75

# 3. VAF statistics on ctDNA mutations - bladder
round(np.median(base_bladder_somatic["VAF_t"]), 2) # median
[min(base_bladder_somatic["VAF_t"]), max(base_bladder_somatic["VAF_t"])] # range
q75, q25 = np.percentile(base_bladder_somatic["VAF_t"], [75 ,25])
q25, q75

# mann whitney u test comparing the ctDNA mutation vaf in bladder vs kidney
mannwhitneyu(base_kidney_somatic["VAF_t"], base_bladder_somatic["VAF_t"])

# mann whitney u test comparing the ctDNA vs CH in kidney vs bladder
mannwhitneyu(base_kidney_somatic["VAF_t"], base_kidney_chip["VAF_t"])
round(np.median(base_kidney_somatic["VAF_t"]), 2) # median
round(np.median(base_kidney_chip["VAF_t"]), 2) # median

mannwhitneyu(base_bladder_somatic["VAF_t"], base_bladder_chip["VAF_t"])
round(np.median(base_bladder_somatic["VAF_t"]), 2) # median
round(np.median(base_bladder_chip["VAF_t"]), 2) # median




###################################################################################
# 4. On ctDNA fraction
df_ct_fraction = pd.read_csv(PATH_ctDNA_fraction)

df_ct_fraction['Date_collected'] = pd.to_datetime(df_ct_fraction['Date_collected'], format='%Y-%m-%d')
sample_info['Date'] = pd.to_datetime(sample_info['Date'], format='%Y%b%d')

# Merge the dataframes on 'Patient_id' and 'Date_collected' / 'Date'
df_ct_fraction = pd.merge(df_ct_fraction, sample_info, how='inner', left_on=['Patient_id', 'Date_collected', 'Diagnosis'], right_on=['Patient_id', 'Date', 'Diagnosis'])
df_ct_fraction = df_ct_fraction.drop(columns=['Date'])

kidney_ctDNA_fr = df_ct_fraction[(df_ct_fraction["Diagnosis"] == "Kidney") & (df_ct_fraction["Timepoint"] == "Baseline")]
bladder_ctDNA_fr = df_ct_fraction[(df_ct_fraction["Diagnosis"] == "Bladder") & (df_ct_fraction["Timepoint"] == "Baseline")]

# Stats - KIDNEY
kidney_ctDNA_pos = kidney_ctDNA_fr[kidney_ctDNA_fr["ctDNA_status"] == "ctDNA positive"]
round(np.median(kidney_ctDNA_pos["Mutation_ctDNA_fraction"])*100, 3) # median
[round(min(kidney_ctDNA_pos["Mutation_ctDNA_fraction"]*100), 2), round(max(kidney_ctDNA_pos["Mutation_ctDNA_fraction"])*100, 2)] # range
q75, q25 = np.percentile(kidney_ctDNA_pos["Mutation_ctDNA_fraction"], [75 ,25])
iqr = round((q75 - q25)*100, 2)

# Stats - BLADDER
bladder_ctDNA_pos = bladder_ctDNA_fr[bladder_ctDNA_fr["ctDNA_status"] == "ctDNA positive"]
round(np.median(bladder_ctDNA_pos["Mutation_ctDNA_fraction"])*100, 3) # median
[round(min(bladder_ctDNA_pos["Mutation_ctDNA_fraction"]*100), 3), round(max(bladder_ctDNA_pos["Mutation_ctDNA_fraction"])*100, 2)] # range
q75, q25 = np.percentile(bladder_ctDNA_pos["Mutation_ctDNA_fraction"], [75 ,25])
iqr = round((q75 - q25)*100, 2)

###################################################################################
# 5. Association between CH detection and ctDNA detection
bladder_ctDNA_status = annotate_mutation_status(base_bladder_somatic, "Bladder", PATH_sample_information, annotate_what = "ctDNA", annotate_gene = False)
bladder_ctDNA_status = bladder_ctDNA_status[bladder_ctDNA_status["Timepoint"] == "Baseline"].replace({"Positive": "ctDNA positive", "Negative": "ctDNA negative"})

kidney_ctDNA_status = annotate_mutation_status(base_kidney_somatic, "Kidney", PATH_sample_information, annotate_what = "ctDNA", annotate_gene = False)
kidney_ctDNA_status = kidney_ctDNA_status[kidney_ctDNA_status["Timepoint"] == "Baseline"].replace({"Positive": "ctDNA positive", "Negative": "ctDNA negative"})

bladder_ch_status = annotate_mutation_status(base_bladder_chip, "Bladder", PATH_sample_information, annotate_what = "CHIP", annotate_gene = False)
bladder_ch_status = bladder_ch_status[bladder_ch_status["Timepoint"] == "Baseline"].replace({"Positive": "CH positive", "Negative": "CH negative"})

kidney_ch_status = annotate_mutation_status(base_kidney_chip, "Kidney", PATH_sample_information, annotate_what = "CHIP", annotate_gene = False)
kidney_ch_status = kidney_ch_status[kidney_ch_status["Timepoint"] == "Baseline"].replace({"Positive": "CH positive", "Negative": "CH negative"})

bladder = bladder_ctDNA_status.merge(bladder_ch_status)
kidney = kidney_ctDNA_status.merge(kidney_ch_status)
merged = pd.concat([bladder, kidney])

# Create the contingency table
contingency_table = pd.crosstab(merged['ctDNA status'], merged['CHIP status'])
odds_ratio, p_value = fisher_exact(contingency_table)

###################################################################################
# 6. CH interference stats
xx = make_ch_ctDNA_prevalence_df(baseline_ctDNA, baseline_chip, min_vaf = 0)
xx_2_perc = make_ch_ctDNA_prevalence_df(baseline_ctDNA, baseline_chip, min_vaf = 2)

# Number of patients with CH mutations in 13 genes commonly mutated in GU malignancies
GU_malignancy_13_genes = ["CHEK2", "ATM", "TERT", "TP53", "BRCA1", "BRCA2", "ARID1A", "ERBB2", "AR", "KMT2D", "PBRM1", "BAP1", "SETD2"]
n_pts_GU = baseline_chip[baseline_chip["Gene"].isin(GU_malignancy_13_genes)]["Patient_id"].unique().size
round(n_pts_GU/(n_bladder + n_kidney)*100, 2)

# when subsetting to 2 perc
baseline_chip_2_perc = baseline_chip[baseline_chip["VAF_n"] > 2]
n_pts_GU_2_perc = baseline_chip_2_perc[baseline_chip_2_perc["Gene"].isin(GU_malignancy_13_genes)]["Patient_id"].unique().size
round(n_pts_GU_2_perc/(n_bladder + n_kidney)*100, 2)

# Origin of mutations in select genes
df_somatic = baseline_ctDNA[["Gene"]].sort_values("Gene").value_counts().reset_index()
df_somatic.columns = ["Gene", "Counts_ctDNA"]
df_ch = baseline_chip[["Gene"]].sort_values("Gene").value_counts().reset_index()
df_ch.columns = ["Gene", "Counts_ch"]   

# Get percentage values
combined = df_ch.merge(df_somatic, on = ["Gene"], how = "outer").fillna(0)
combined["summed"] = combined["Counts_ch"] + combined["Counts_ctDNA"]
combined["perc_ch"] = (combined["Counts_ch"]/combined["summed"])*100
combined["perc_ctDNA"] = (combined["Counts_ctDNA"]/combined["summed"])*100
combined.sort_values(by = "perc_ch", inplace = True, ascending = False)

# Run MWU in ctDNA vs CH mutation VAF for all genes, correct for multiple testing
# Store p-values and genes
p_values = []
genes = combined_filtered["Gene"].tolist()

for gene in genes:
    ctdna_vafs = baseline_ctDNA[baseline_ctDNA["Gene"] == gene]["VAF_t"].tolist()
    ch_vafs = baseline_chip[baseline_chip["Gene"] == gene]["VAF_t"].tolist()
    
    if ctdna_vafs and ch_vafs:  # Ensure both lists are non-empty
        stat, p_value = mannwhitneyu(ctdna_vafs, ch_vafs)
        p_values.append(p_value)
    else:
        p_values.append(1.0)  # Append a non-significant p-value if one of the lists is empty

# Correct for multiple testing using Benjamini-Hochberg
_, corrected_p_values, _, _ = multipletests(p_values, method='fdr_bh')

# Print genes with significant p-values and their median VAFs
for gene, corrected_p_value in zip(genes, corrected_p_values):
    if corrected_p_value < 0.05:
        ctdna_vafs = baseline_ctDNA[baseline_ctDNA["Gene"] == gene]["VAF_t"].tolist()
        ch_vafs = baseline_chip[baseline_chip["Gene"] == gene]["VAF_t"].tolist()
        
        median_ctdna_vaf = np.median(ctdna_vafs) if ctdna_vafs else float('nan')
        median_ch_vaf = np.median(ch_vafs) if ch_vafs else float('nan')
        
        print(f"Gene: {gene}")
        print(f"Median ctDNA VAF: {median_ctdna_vaf}")
        print(f"Median CH VAF: {median_ch_vaf}\n")
