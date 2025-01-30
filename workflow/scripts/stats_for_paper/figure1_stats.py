import os
import pandas as pd
import numpy as np
import re
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from scipy.stats import f_oneway
from scipy.stats import fisher_exact
from scipy.stats import mannwhitneyu

# LOAD DATASETS
DIR_working = "/groups/wyattgrp/users/amunzur/pipeline"
PATH_sample_information = os.path.join(DIR_working, "resources/sample_lists/sample_information.tsv")
PATH_mutation_ctfractions = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/ctDNA_fractions.csv"
PATH_bladder_clinical_data = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/bladder/clinical_data.csv"
PATH_kidney_clinical_data = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/RCC clinical data - mRCC clinical Data.csv"

source_functions = os.path.join(DIR_working, "workflow/scripts/visualization/UTILITIES_make_chip_plots.py")
with open(source_functions, 'r') as file:
    script_code = file.read()

exec(script_code)

# ctDNA mutations
all_vars_ctDNA = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/somatic.csv"))
all_vars_ctDNA = all_vars_ctDNA[(all_vars_ctDNA["Dependent"] == False) & (all_vars_ctDNA["Timepoint"] == "Baseline")]

# CH mutations
all_vars_chip = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/chip.csv"))
all_vars_chip = all_vars_chip[(all_vars_chip["Dependent"] == False) & (all_vars_chip["Timepoint"] == "Baseline")]
bladder_chip = all_vars_chip[all_vars_chip["Diagnosis"] == "Bladder"]
kidney_chip = all_vars_chip[all_vars_chip["Diagnosis"] == "Kidney"]

# 1. In these select genes, are the median vafs different between kidney and bladder?
genes_list = ["DNMT3A", "TET2", "PPM1D", "ASXL1", "TP53", "ATM", "CHEK2", "GNAS", "KMT2D", "STAG2", "CBL", "RAD21", "BRCC3", "JAK2", "SF3B1"]
p_values = []

# Perform Mann-Whitney U test for each gene and collect p-values
for gene in genes_list: 
    bladder_vafs = bladder_chip[bladder_chip["Gene"] == gene]["VAF_n"].tolist()
    kidney_vafs = kidney_chip[kidney_chip["Gene"] == gene]["VAF_n"].tolist()
    stat, p_value = mannwhitneyu(bladder_vafs, kidney_vafs)
    p_values.append(p_value)

# Adjust p-values using Benjamini-Hochberg correction
adjusted_p_values = multipletests(p_values, method="fdr_bh")[1]
for gene, adj_p_value in zip(genes_list, adjusted_p_values):
    print(f"{gene} adjusted p-value: {adj_p_value}")

# For DNMT3A print the medians
bladder_chip[bladder_chip["Gene"] == "DNMT3A"]["VAF_n"].median()
kidney_chip[kidney_chip["Gene"] == "DNMT3A"]["VAF_n"].median()

# 2. The nature of mutations in TET2
all_vars_chip[all_vars_chip["Gene"] == "TET2"]["Consequence"].value_counts(normalize = True)

# 3. Relationship between number of mutations and VAF
mut_counts = all_vars_chip["Patient_id"].value_counts().reset_index()
mut_counts.columns = ["Patient_id", "count"]
mut_counts.loc[mut_counts["count"] > 5, "count"] = ">5"
vaf_df = all_vars_chip[["Patient_id", "VAF_n"]].merge(mut_counts)
groups = vaf_df.groupby('count')['VAF_n'].apply(list)

# Output the results
f_stat, p_value = f_oneway(*groups) # Run ANOVA
print(f"ANOVA F-statistic: {f_stat}")
print(f"P-value: {p_value}")

# 4. About the prevalence of CH when mandating different thresholds
# 2 percent 
n_pts_2_perc = all_vars_chip[all_vars_chip["VAF_n"]>=2].shape[0]
round(n_pts_2_perc/301*100, 2)
# 10 percent
n_pts_10_perc = all_vars_chip[all_vars_chip["VAF_n"]>=10].shape[0]
round(n_pts_10_perc/301*100, 2)

# 5. Fishers exact test - prevalence of CH and disease stage at initial diagnosis
bladder_clinical_path = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/bladder/clinical_data.csv"
kidney_clinical_path = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/CHIP Supplementary tables - Clinical_data_mRCC.csv"

bladder_df_clin_main = pd.read_csv(bladder_clinical_path)
bladder_df_clin_main = bladder_df_clin_main[bladder_df_clin_main["First sample?"]].reset_index(drop = True)
bladder_df = bladder_df_clin_main[["Patient_id", "Disease stage at initial diagnosis"]].reset_index(drop = True)
bladder_df["Disease stage at initial diagnosis"] = bladder_df["Disease stage at initial diagnosis"].str.extract(r'(Localized|Metastatic)')

kidney_df_clin_main = pd.read_csv(kidney_clinical_path)
kidney_df = kidney_df_clin_main[["Patient_id", "Disease stage at initial diagnosis"]]
disease_stage_df = pd.concat([bladder_df, kidney_df]).reset_index(drop = True)

ch_status = annotate_mutation_status(all_vars_chip, "Both", PATH_sample_information, annotate_what = "CHIP")
ch_status = ch_status[ch_status["Timepoint"] == "Baseline"].reset_index(drop = True)[["Patient_id", "CHIP status"]]
disease_stage_df = disease_stage_df.merge(ch_status, how = "inner")

# Creating a contingency table
contingency_table = pd.crosstab(disease_stage_df['Disease stage at initial diagnosis'], disease_stage_df['CHIP status'])
print(contingency_table)
odds_ratio, p_value = fisher_exact(contingency_table)
print(f"Odds Ratio: {odds_ratio}")
print(f"P-value: {p_value}")

# 6. Check age between Disease stage at initial diagnosis localized and metastatic patients
age_df = pd.concat([
    bladder_df_clin_main[["Patient_id", "Age at blood draw"]], 
    kidney_df_clin_main[["Patient_id", "Age at GUBB draw"]].rename(columns = {"Age at GUBB draw": "Age at blood draw"})
    ])
disease_stage_and_age = disease_stage_df.merge(age_df)

# Separate age data for Localized and Metastatic groups
localized_ages = disease_stage_and_age[disease_stage_and_age['Disease stage at initial diagnosis'] == 'Localized']['Age at blood draw']
metastatic_ages = disease_stage_and_age[disease_stage_and_age['Disease stage at initial diagnosis'] == 'Metastatic']['Age at blood draw']
stat, p_value = mannwhitneyu(localized_ages, metastatic_ages, alternative='two-sided')
np.median(localized_ages)
np.median(metastatic_ages)

# Perform Mann-Whitney U test
# Kidney
kidney_age_stage_df = kidney_df_clin_main[["Patient_id", "Disease stage at initial diagnosis"]].merge(age_df)
kidney_localized_ages = kidney_age_stage_df[kidney_age_stage_df['Disease stage at initial diagnosis'] == 'Localized']['Age at blood draw']
kidney_metastatic_ages = kidney_age_stage_df[kidney_age_stage_df['Disease stage at initial diagnosis'] == 'Metastatic']['Age at blood draw']
stat, p_value = mannwhitneyu(localized_ages, metastatic_ages, alternative='two-sided')
np.median(kidney_localized_ages)
np.median(kidney_metastatic_ages)

# Bladder
bladder_age_stage_df = bladder_df_clin_main[["Patient_id", "Disease stage at initial diagnosis"]].merge(age_df)
bladder_localized_ages = bladder_age_stage_df[bladder_age_stage_df['Disease stage at initial diagnosis'].isin(["Localized MIBC", "Localized NMIBC"])]['Age at blood draw']
bladder_metastatic_ages = bladder_age_stage_df[bladder_age_stage_df['Disease stage at initial diagnosis'] == 'Metastatic']['Age at blood draw']
stat, p_value = mannwhitneyu(bladder_localized_ages, bladder_metastatic_ages, alternative='two-sided')
np.median(bladder_localized_ages)
np.median(bladder_metastatic_ages)

# 7. Are ASXL1 mutations more common in smokers? In the blader cohort
bladder_smokers = bladder_df_clin_main[["Patient_id", "Previous smoking history"]].replace({"Previous smoker": "smoker", "Current smoker": "smoker"})
asxl1_status = annotate_mutation_status(bladder_chip, "Bladder", PATH_sample_information, annotate_what = "CHIP", annotate_gene = "ASXL1")
asxl1_status = asxl1_status[asxl1_status["Timepoint"] == "Baseline"].reset_index(drop = True)[["Patient_id", "CHIP status"]]
asxl1_status = asxl1_status.merge(bladder_smokers)

contingency_table_asxl1 = pd.crosstab(asxl1_status['Previous smoking history'], asxl1_status['CHIP status'])
print(contingency_table_asxl1)
odds_ratio, p_value = fisher_exact(contingency_table_asxl1)
print(f"Odds Ratio: {odds_ratio}")
print(f"P-value: {p_value}")

