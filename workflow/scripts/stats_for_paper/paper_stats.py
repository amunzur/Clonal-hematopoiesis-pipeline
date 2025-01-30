import pandas as pd
import numpy as np
import os 
import re
from scipy.stats import ttest_ind
from lifelines import KaplanMeierFitter
from datetime import datetime
from lifelines.statistics import logrank_test
from scipy.stats import fisher_exact
from lifelines import KaplanMeierFitter
from scipy.stats import mannwhitneyu
from lifelines import CoxPHFitter
from statsmodels.stats.multitest import multipletests
from scipy.stats import spearmanr
from scipy.stats import chi2_contingency

def compare_ages(chip_df, gene, PATH_sample_information, age_df, colname = "Gene"):
    """
    Runs MWU to see if patients with CH mutations in a certain gene are older than patients with no mutations in that gene.
    """
    # Ensure genes is a list
    if not isinstance(gene, list):
        gene = [gene]
    mut_pos = chip_df[chip_df[colname].isin(gene)][["Patient_id"]].drop_duplicates().assign(Stat = "+")
    sample_info = pd.read_csv(PATH_sample_information, sep="\t", names=["Patient_id", "Date", "Diagnosis", "Timepoint"])
    sample_info = sample_info[sample_info["Timepoint"] == "Baseline"]
    a = sample_info.merge(mut_pos, how = "left")
    a["Stat"] = a["Stat"].fillna("-")
    a = age_df.merge(a, how = "left", on = "Patient_id")
    # Run MWU
    g1 = a[a['Stat'] == '+']['Age']
    g2 = a[a['Stat'] == '-']['Age']
    stat, p_value = mannwhitneyu(g1, g2)
    pos_median_age = np.median(g1)
    neg_median_age = np.median(g2)
    return gene, p_value, pos_median_age, neg_median_age

def run_multiple_testing(chip_df, genes, PATH_sample_information, age_df):
    results = []
    
    for gene in genes:
        gene, p_value, pos_median_age, neg_median_age = compare_ages(chip_df, gene, PATH_sample_information, age_df)
        results.append((gene, p_value, pos_median_age, neg_median_age))
    
    # Extract p-values
    p_values = [result[1] for result in results]
    
    # Apply Benjamini-Hochberg correction
    reject, pvals_corrected, _, _ = multipletests(p_values, method='fdr_bh')
    
    # Print results with corrected p-values
    for i, (gene, p_value, pos_median_age, neg_median_age) in enumerate(results):
        if reject[i]:
            print(f"{gene}, {pvals_corrected[i]:.2E}, + median age: {pos_median_age}, - median age: {neg_median_age}")


source_functions = os.path.join(DIR_working, "workflow/scripts/visualization/UTILITIES_make_chip_plots.py")
with open(source_functions, 'r') as file:
    script_code = file.read()

exec(script_code)

# number of kidney patients we profiled
DIR_working = "/groups/wyattgrp/users/amunzur/pipeline"
PATH_sample_information = os.path.join(DIR_working, "resources/sample_lists/sample_information.tsv")
sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])

kidney_df = sample_info[(sample_info["Diagnosis"] == "Kidney") ]
bladder_df = sample_info[(sample_info["Diagnosis"] == "Bladder")]

# Number of baseline samples & number of patients
n_baseline_kidney = kidney_df[kidney_df["Timepoint"] == "Baseline"].shape[0]
n_baseline_bladder = bladder_df[bladder_df["Timepoint"] == "Baseline"].shape[0]

# total number of patients we included in the study / genomic analysis, including only mRCC and mUC patients
n_baseline_kidney + n_baseline_bladder

# Number of OT samples
n_ot_kidney = kidney_df[kidney_df["Timepoint"] == "During treatment"].shape[0]
n_ot_bladder = bladder_df[bladder_df["Timepoint"] == "During treatment"].shape[0]

# the total number of samples that come from the metastatic patients
n_baseline_kidney + n_baseline_bladder + n_ot_kidney + n_ot_bladder

# number of pts with n OT sample - kidney
x = kidney_df[kidney_df["Timepoint"] == "During treatment"].groupby("Patient_id").size().value_counts().sort_index().reset_index()
x.columns = ["Number of OT samples", "Number of kidney pts"]
x

# number of pts with n OT sample - bladder
y = bladder_df[bladder_df["Timepoint"] == "During treatment"].groupby("Patient_id").size().value_counts().sort_index().reset_index()
y.columns = ["Number of OT samples", "Number of bladder pts"]
y

# CHIP mutations stats
path_muts = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/chip.csv"
all_vars_chip = pd.read_csv(path_muts)
all_vars_chip = all_vars_chip[all_vars_chip["Dependent"] == False].reset_index(drop = True)

baseline_chip = all_vars_chip[all_vars_chip["Timepoint"] == "Baseline"]
ot_chip = all_vars_chip[all_vars_chip["Timepoint"] == "During treatment"]

base_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
prog_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)
prog_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)

# total number of mutations at baseline
baseline_chip.shape[0]

# Number of unique genes that had a mutation
baseline_chip["Gene"].unique().shape[0]

# Number of baseline samples that had a CH mutation
baseline_chip["Patient_id"].unique().shape[0]

# median VAF of baseline CH mutations
round(np.median(baseline_chip["VAF_t"]), 2)
q1 = np.percentile(baseline_chip["VAF_t"], 25)
q3 = np.percentile(baseline_chip["VAF_t"], 75)
iqr = q3 - q1

# median VAF of baseline CH mutations
round(np.median(baseline_chip["VAF_n"]), 2)
q1 = np.percentile(baseline_chip["VAF_n"], 25)
q3 = np.percentile(baseline_chip["VAF_n"], 75)
iqr = q3 - q1

# the number of kidney patients with a baseline mut and %CH positive in the kidney cohort
a = base_kidney_chip["Patient_id"].unique().shape[0]
round(a/n_baseline_kidney*100, 1) # fraction CH positive in kidneys

# the number of bladder patients with a baseline mut
b = base_bladder_chip["Patient_id"].unique().shape[0]
round(b/n_baseline_bladder*100, 1) # fraction CH positive in bladders

# entire cohort the fraction of patients that are CH+
(a+b)/(n_baseline_kidney+n_baseline_bladder)

# in the kidney cohort the number of pts and the number of muts
q = base_kidney_chip["Patient_id"].value_counts().reset_index()["Patient_id"].value_counts().reset_index()
q.columns = ["n_muts", "count in kidney baseline"]
q.sort_values(by = "n_muts", inplace = True)
q["kidney frac"] = round((q["count in kidney baseline"]/n_baseline_kidney)*100, 2)

w = base_bladder_chip["Patient_id"].value_counts().reset_index()["Patient_id"].value_counts().reset_index()
w.columns = ["n_muts", "count in bladder baseline"]
w.sort_values(by = "n_muts", inplace = True)
w["bladder frac"] = round((w["count in bladder baseline"]/n_baseline_bladder)*100, 2)

# percentage of ppl that had more than one mutation 
q[q["n_muts"] > 1]["kidney frac"].sum() # kidney
w[w["n_muts"] > 1]["bladder frac"].sum() # kidney

# Number of genes that had a CH mutation
baseline_chip["Gene"].unique().shape[0]

# Top genes mutated - number and percentage
gene_counts = baseline_chip["Gene"].value_counts()
gene_percentages = round(gene_counts / gene_counts.sum() * 100, 2) # Calculate the percentages

# # Top genes mutated KIDNEY - number and percentage in the entire cohort
# gene_counts_kidney = base_kidney_chip["Gene"].value_counts()
# gene_percentages_kidney = round(gene_counts_kidney / gene_counts_kidney.sum() * 100, 2) # Calculate the percentages among the CH+ ppl
# gene_percentages_kidney_all = round(gene_counts_kidney / n_baseline_kidney * 100, 2) # Calculate the percentages in the entire kidney group

# # Top genes mutated BLADDER - number and percentage
# gene_counts_bladder = base_bladder_chip["Gene"].value_counts()
# gene_percentages_bladder = round(gene_counts_bladder / gene_counts_bladder.sum() * 100, 2) # Calculate the percentages among the CH+ ppl
# gene_percentages_bladder_all = round(gene_counts_bladder / n_baseline_bladder * 100, 2) # Calculate the percentages in the entire kidney group

# Top genes mutated KIDNEY - number and percentage in the CH+ cohort
n_kidney_denom = base_kidney_chip["Patient_id"].unique().shape[0]
gene_counts_kidney = base_kidney_chip[["Patient_id", "Gene"]].drop_duplicates()["Gene"].value_counts().reset_index()
gene_counts_kidney["Perc in CH+"] = gene_counts_kidney["Gene"]/n_kidney_denom
gene_counts_kidney["Perc in CH+"] = round(gene_counts_kidney["Perc in CH+"]*100, 2)

# Top genes mutated BLADDER - number and percentage in the CH+ cohort
n_bladder_denom = base_bladder_chip["Patient_id"].unique().shape[0]
gene_counts_bladder = base_bladder_chip[["Patient_id", "Gene"]].drop_duplicates()["Gene"].value_counts().reset_index()
gene_counts_bladder["Perc in CH+"] = gene_counts_bladder["Gene"]/n_bladder_denom
gene_counts_bladder["Perc in CH+"] = round(gene_counts_bladder["Perc in CH+"]*100, 2)

# Percent of pts that had 1, 2 and 3 mutations
q = baseline_chip["Patient_id"].value_counts().reset_index()["Patient_id"].value_counts().reset_index().sort_values(by = "index")
q.columns = ["n_muts", "n_patients"]
total_patients = q['n_patients'].sum()
q['percentage'] = (q['n_patients'] / total_patients) * 100

# Percent and number of pts with more than 1 mutation
q[q["n_muts"] > 1]["n_patients"].sum()
round(q[q["n_muts"] > 1]["percentage"].sum(), 2)

# Percent and number of pts with more than 5 mutation
q[q["n_muts"] > 5]["n_patients"].sum()
round(q[q["n_muts"] > 5]["percentage"].sum(), 2)

# median VAF of mutations in patients that have more than 1 CH mutation
pts_list = 

aa = baseline_chip["Patient_id"].value_counts().reset_index()
aa = aa[aa["Patient_id"] > 1]

baseline_chip[baseline_chip["Patient_id"].isin(aa["index"])]["VAF_n"].median()

# Fraction of bladder and kidney cohorts that is CH+ at 2% and 10% VAFs
kidney_2_perc = base_kidney_chip[base_kidney_chip["VAF_n"] >= 2]["Patient_id"].unique().shape[0]/n_baseline_kidney
kidney_10_perc = base_kidney_chip[base_kidney_chip["VAF_n"] >= 10]["Patient_id"].unique().shape[0]/n_baseline_kidney

bladder_2_perc = base_bladder_chip[base_bladder_chip["VAF_n"] >= 2]["Patient_id"].unique().shape[0]/n_baseline_bladder
bladder_10_perc = base_bladder_chip[base_bladder_chip["VAF_n"] >= 10]["Patient_id"].unique().shape[0]/n_baseline_bladder


###################################################
# ASSOCIATIONS WITH AGE
###################################################
PATH_clinical_bladder = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/bladder/clinical_data.csv"
PATH_clinical_kidney = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/Supplementary tables - Clinical data - mRCC.csv"

# Load and rename columns
kidney_clin = pd.read_csv(PATH_clinical_kidney).rename(columns={"GUBB ID": "Patient_id", "Age at GUBB draw": "Age"})
bladder_clin = pd.read_csv(PATH_clinical_bladder).rename(columns={"GUBB ID": "Patient_id", "Age at blood draw": "Age"})
bladder_clin = bladder_clin[bladder_clin["First sample?"]].reset_index(drop = True)

# Subset to samples we are including in the study
kidney_pts = sample_info[(sample_info["Diagnosis"] == "Kidney") & (sample_info["Timepoint"] == "Baseline")]["Patient_id"].tolist()
bladder_pts = sample_info[(sample_info["Diagnosis"] == "Bladder") & (sample_info["Timepoint"] == "Baseline")]["Patient_id"].tolist()

# Filter the clinical data based on the patients included in the study
kidney_clin = kidney_clin[kidney_clin["Patient_id"].isin(kidney_pts)]
bladder_clin = bladder_clin[bladder_clin["Patient_id"].isin(bladder_pts)]

age_df = pd.concat([kidney_clin[["Patient_id", "Age"]], bladder_clin[["Patient_id", "Age"]]]).reset_index(drop = True)

x = kidney_clin[["Patient_id", "Age"]].merge(base_kidney_chip[["Patient_id", "Diagnosis", "VAF_n", 'Gene', 'Chrom', 'Position', 'Ref', 'Alt', "Consequence", "Protein_annotation"]], how = "inner")
y = bladder_clin[["Patient_id", "Age"]].merge(base_bladder_chip[["Patient_id", "Diagnosis", "VAF_n", 'Gene', 'Chrom', 'Position', 'Ref', 'Alt', "Consequence", "Protein_annotation"]], how = "inner")
z = pd.concat([x, y]).reset_index(drop = True)

# Are CH+ patients older than CH- patients?
ch_status = annotate_mutation_status(baseline_chip, "Both", PATH_sample_information, annotate_what = "CHIP")
ch_status = ch_status[ch_status["Timepoint"] == "Baseline"]
ch_pos = ch_status[ch_status["CHIP status"] == "Positive"].merge(age_df)
ch_neg = ch_status[ch_status["CHIP status"] == "Negative"].merge(age_df)
np.median(ch_pos["Age"])
np.median(ch_neg["Age"])
mannwhitneyu(ch_pos["Age"], ch_neg["Age"])


# Among CH+ patien at baseline, are patients with splicing mutations older than patients with DTA mutations?
splicing_df = annotate_mutation_status(baseline_chip, "Both", PATH_sample_information, "CHIP", annotate_gene = ["SF3B1", "SRSF2", "U2AF1", "SH2B3"], drop_dependent = True)
dnmt3a_df = annotate_mutation_status(baseline_chip, "Both", PATH_sample_information, "CHIP", annotate_gene = ["DNMT3A"], drop_dependent = True)

splicing_df = splicing_df[(splicing_df["Timepoint"] == "Baseline") & (splicing_df["CHIP status"] == "Positive")].merge(age_df).rename(columns = {"CHIP status": "Splicing status"})
dnmt3a_df = dnmt3a_df[(dnmt3a_df["Timepoint"] == "Baseline") & (dnmt3a_df["CHIP status"] == "Positive")].merge(age_df).rename(columns = {"CHIP status": "DTA status"})
merged = splicing_df.merge(dnmt3a_df, how = "outer", on = ["Patient_id", "Age", "Diagnosis", "Timepoint"], indicator = True)

only_splicing = merged[merged["_merge"] == "left_only"]
only_dnmt3a = merged[merged["_merge"] == "right_only"]
both = merged[merged["_merge"] == "both"]

stat, p_value = mannwhitneyu(both["Age"], only_dnmt3a["Age"])

np.median(only_splicing["Age"])
np.median(only_dnmt3a["Age"])
np.median(both["Age"])

# Run associations. Are patients with XX genes are significantly older? 
gene_groups = [["DNMT3A", "ASXL1", "TET2"], ["TP53", "PPM1D", "CHEK2", "BRCA1", "BRCA2"], ["SF3B1", "SRSF2", "U2AF1", "SH2B3"]]
run_multiple_testing(baseline_chip, gene_groups, PATH_sample_information, age_df)

# CH+ vs CH- age comparison
ch_status_df = annotate_mutation_status(baseline_chip, "Both", PATH_sample_information, "CHIP", annotate_gene = False, drop_dependent = True)
ch_status_df = ch_status_df.merge(age_df, how =  "inner")
ch_status_df = ch_status_df[ch_status_df["Timepoint"] == "Baseline"]
g1 = ch_status_df[ch_status_df["CHIP status"] == "Positive"]["Age"]
g2 = ch_status_df[ch_status_df["CHIP status"] == "Negative"]["Age"]
stat, p_value = mannwhitneyu(g1, g2)
pos_median_age = np.median(g1)
neg_median_age = np.median(g2)
f"{p_value:.2E}"
ch_status_df["CHIP status"].value_counts()

# comparing the median age at baseline, bladder vs kidney
kidney_cohort = sample_info[(sample_info["Diagnosis"] == "Kidney") & (sample_info["Timepoint"] == "Baseline")].reset_index(drop = True).merge(age_df)
bladder_cohort = sample_info[(sample_info["Diagnosis"] == "Bladder") & (sample_info["Timepoint"] == "Baseline")].reset_index(drop = True).merge(age_df)
stat, p_value = mannwhitneyu(kidney_cohort["Age"], bladder_cohort["Age"])
np.median(kidney_cohort["Age"])
np.median(bladder_cohort["Age"])

# Are older patients more likely to have larger VAF mutations? PAN CANCER
max_vaf_df = z.groupby("Patient_id")["VAF_n"].max().reset_index()
max_vaf_df = max_vaf_df.merge(age_df, how = "inner")
# Calculate Spearman's rank correlation
correlation, p_value = spearmanr(max_vaf_df['Age'], max_vaf_df['VAF_n'])
print(f"Spearman's rank correlation: {correlation:.2f}")
print(f"P-value: {p_value:.2E}")

# Are older patients more likely to have larger VAF mutations? Kidney and bladder separate
max_vaf_df_kidney = base_kidney_chip.groupby("Patient_id")["VAF_n"].max().reset_index()
max_vaf_df_kidney = max_vaf_df_kidney.merge(age_df, how = "inner")
correlation, p_value = spearmanr(max_vaf_df_kidney['Age'], max_vaf_df_kidney['VAF_n'])

max_vaf_df_bladder = base_bladder_chip.groupby("Patient_id")["VAF_n"].max().reset_index()
max_vaf_df_bladder = max_vaf_df_bladder.merge(age_df, how = "inner")
correlation, p_value = spearmanr(max_vaf_df_bladder['Age'], max_vaf_df_bladder['VAF_n'])

# Do older patients have more mutations? PAN CANCER
g = z["Patient_id"].value_counts().reset_index().rename(columns = {"index": "Patient_id", "Patient_id": "n_muts"})
g = g.merge(age_df, how = "inner")
correlation, p_value = spearmanr(g['Age'], g['n_muts'])

# Do older patients have more mutations? Kidney and bladder separate
grouped = z.groupby('Diagnosis')['Patient_id'].value_counts().reset_index(name='n_muts')
grouped = grouped.merge(age_df, how = "inner")
for diagnosis, group in grouped.groupby("Diagnosis"):
    correlation, p_value = spearmanr(group['Age'], group['n_muts'])
    correlation = round(correlation, 4)
    p_value = round(p_value, 4)
    print(f"{diagnosis}, {correlation}, {p_value}")

# Portion of the cohort that is CH+
ch_status_df = annotate_mutation_status(baseline_chip, "Both", PATH_sample_information, "CHIP", annotate_gene = False, drop_dependent = True)
baseline_ch_status = ch_status_df[ch_status_df["Timepoint"] == "Baseline"]
baseline_ch_status= baseline_ch_status["CHIP status"].value_counts().reset_index().rename(columns = {"index": "Status", "CHIP status": "Count"})
n_total = baseline_ch_status["Count"].sum()
baseline_ch_status["Perc"] = round(baseline_ch_status["Count"] /n_total*100, 2)

# Portion of the cohort that is CH+ when we increase the detection threshold to 2% and 10%.
ch_status_2_perc = annotate_mutation_status(baseline_chip[baseline_chip["VAF_n"] >= 2], "Both", PATH_sample_information, "CHIP", annotate_gene = False, drop_dependent = True)
ch_status_2_perc = ch_status_2_perc[ch_status_2_perc["Timepoint"] == "Baseline"]
ch_status_2_perc = ch_status_2_perc["CHIP status"].value_counts().reset_index()
ch_status_2_perc["Perc"] = ch_status_2_perc["CHIP status"]/(ch_status_2_perc["CHIP status"].sum())

ch_status_10_perc = annotate_mutation_status(baseline_chip[baseline_chip["VAF_n"] >= 10], "Both", PATH_sample_information, "CHIP", annotate_gene = False, drop_dependent = True)
ch_status_10_perc = ch_status_10_perc[ch_status_10_perc["Timepoint"] == "Baseline"]
ch_status_10_perc= ch_status_10_perc["CHIP status"].value_counts().reset_index()
ch_status_10_perc["Perc"] = ch_status_10_perc["CHIP status"]/(ch_status_10_perc["CHIP status"].sum())

# Is there a difference in CH prevalence per cancer type?
ch_status_df = annotate_mutation_status(baseline_chip, "Both", PATH_sample_information, "CHIP", annotate_gene = False, drop_dependent = True)
baseline_ch_status = ch_status_df[ch_status_df["Timepoint"] == "Baseline"]
contingency_table = pd.crosstab(baseline_ch_status['Diagnosis'], baseline_ch_status['CHIP status']) # Create a contingency table of Diagnosis vs CHIP status
chi2, p_val, _, _ = chi2_contingency(contingency_table)
baseline_ch_status[["Diagnosis", "CHIP status"]].value_counts()

# Is there a difference in CH prevalence by biological sex?
ch_status_df = annotate_mutation_status(baseline_chip, "Both", PATH_sample_information, "CHIP", annotate_gene = False, drop_dependent = True)
baseline_ch_status = ch_status_df[ch_status_df["Timepoint"] == "Baseline"]
sex_df = pd.concat([kidney_clin[["Patient_id", "Sex"]], bladder_clin[["Patient_id", "Sex"]]]).reset_index(drop = True)
sex_df = baseline_ch_status.merge(sex_df, how = "inner")
contingency_table = pd.crosstab(sex_df['Sex'], sex_df['CHIP status']) # Create a contingency table of Diagnosis vs CHIP status
chi2, p_val, _, _ = chi2_contingency(contingency_table)
sex_df[["Sex", "CHIP status"]].value_counts().reset_index()

# correlation between WBC and cfDNA VAFs at baseline
correlation, p_value = spearmanr(baseline_chip['VAF_n'], baseline_chip['VAF_t'])
print(f"Spearman's rank correlation: {correlation:.2f}")
print(f"P-value: {p_value:.2E}")

# Stats on large CH mutations
large_muts = baseline_chip[baseline_chip["VAF_n"] > 30][["Patient_id", "Gene", "Protein_annotation", "VAF_n"]]
large_muts["Patient_id"].unique().shape[0] # number of pts with large CH mutations

# DNMT3A mutations stats
dnmt3a_muts = baseline_chip[baseline_chip["Gene"] == "DNMT3A"].reset_index(drop = True)
dnmt3a_missense = dnmt3a_muts[dnmt3a_muts["Effects"] == "missense"]["VAF_n"]
dnmt3a_trunc = dnmt3a_muts[dnmt3a_muts["Effects"].isin(["frameshift", "splicing", "stop_gain"])]["VAF_n"]

np.median(dnmt3a_missense)
np.median(dnmt3a_trunc)

mannwhitneyu(dnmt3a_missense, dnmt3a_trunc)

r882_muts = dnmt3a_muts[(dnmt3a_muts["Effects"] == "missense") & (dnmt3a_muts["Protein_annotation"].str.contains("R882"))]["VAF_n"]
non_r882_muts = dnmt3a_muts[(dnmt3a_muts["Effects"] == "missense") & (~dnmt3a_muts["Protein_annotation"].str.contains("R882"))]["VAF_n"]
mannwhitneyu(r882_muts, non_r882_muts)

# TET2 mutations stats
tet2_muts = baseline_chip[baseline_chip["Gene"] == "TET2"].reset_index(drop = True)
tet2_trunc = tet2_muts[tet2_muts["Effects"].isin(["frameshift", "splicing", "stop_gain"])]["VAF_n"]

np.median(dnmt3a_missense)
np.median(dnmt3a_trunc)

mannwhitneyu(dnmt3a_missense, dnmt3a_trunc)

r882_muts = dnmt3a_muts[(dnmt3a_muts["Effects"] == "missense") & (dnmt3a_muts["Protein_annotation"].str.contains("R882"))]["VAF_n"]
non_r882_muts = dnmt3a_muts[(dnmt3a_muts["Effects"] == "missense") & (~dnmt3a_muts["Protein_annotation"].str.contains("R882"))]["VAF_n"]
mannwhitneyu(r882_muts, non_r882_muts)

