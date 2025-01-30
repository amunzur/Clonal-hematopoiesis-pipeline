import numpy as np 
import pandas as pd
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import statsmodels.api as sm
import matplotlib.patches as mpatches
from scipy.stats import kruskal
import scipy.stats as stats
from scipy.stats import linregress
from scipy.stats import fisher_exact
from sksurv.nonparametric import kaplan_meier_estimator

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix



with open(source_functions, 'r') as file:
    script_code = file.read()

exec(script_code)

def run_fishers_exact(df, chip_variable, clinical_variable1): 
    """
    Given two variables in a df, make a contingency table and run fishers exact on it. 
    """
    df_per_pt = data[["Patient_id", chip_variable]]
    contingency_table = pd.crosstab(df_per_pt[chip_variable], df[clinical_variable1])
        
    if (contingency_table <= 2).any().any():
        print("Fisher's exact test is not performed due to cells with counts less than or equal to 2.")
        return None
    
    if contingency_table.shape[0] < 2 or contingency_table.shape[1] < 2:
        print("Fisher's exact test is not performed due to the contingency table being less than 2x2.")
        return None
    
    odds_ratio, p_value = fisher_exact(contingency_table)
    return[contingency_table, odds_ratio, p_value]

def run_mwu(df, categorical_variable, numerical_variable):
    g1_value = df[categorical_variable].unique()[0]
    g2_value = df[categorical_variable].unique()[1]
    g1 = df[df[categorical_variable] == g1_value][numerical_variable]
    g2 = df[df[categorical_variable] == g2_value][numerical_variable]
    stat, p_value = mannwhitneyu(g1, g2)
    return(p_value)

PATH_mutations = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/chip.csv"
PATH_ctdna_mutations = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/somatic.csv"
PATH_clin_data = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/Supplementary tables - Clinical data - mRCC.tsv"
PATH_sample_information = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"

sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Diagnosis", "Date_collected", "Timepoint"])
dict = pd.read_csv(PATH_clin_dict, "\t")
muts_main = pd.read_csv(PATH_mutations)
ctdna_main = pd.read_csv(PATH_ctdna_mutations)
clin_main = pd.read_csv(PATH_clin_data, sep = "\t").rename(columns = {"GUBB ID": "Patient_id"})

clin = clin_main[["Patient_id", "irAE", "Age at GUBB draw"]]
clin = clin[~clin["irAE"].isna()].reset_index(drop = True)

# IMMUNE ADVERSE EVENTS
# Does the prevalence of CHIP differ in patients with and without immune adverse events?
clin = clin_main[["Patient_id", "irAE", "Age at GUBB draw", "Charlson_score", "OS from cfDNA collection (mo)"]]
clin = clin[~pd.isnull(clin["irAE"])].reset_index(drop = True)
clin.loc[clin["irAE"] > 1, "irAE"] = 1
# clin[col] = clin[col].astype(int) > 0

# Clean up muts of interest
ch = muts_main[
    (muts_main["Timepoint"] == "Baseline") & 
    (muts_main['Diagnosis'] == 'Kidney') & 
    (muts_main["Dependent"] == False)].reset_index(drop = True)

ctdna = ctdna_main[
    (ctdna_main["Timepoint"] == "Baseline") & 
    (ctdna_main['Diagnosis'] == 'Kidney') & 
    (ctdna_main["Dependent"] == False)].reset_index(drop = True)

# Now annotate CH status based on various criteria
ctdna_status = annotate_mutation_status(ctdna, "Kidney", PATH_sample_information, annotate_what = "ctDNA", annotate_gene = False, drop_dependent = True)
ch_status = annotate_mutation_status(ch, "Kidney", PATH_sample_information, annotate_what = "CHIP", annotate_gene = False, drop_dependent = True).rename(columns = {"CHIP status": "CH status"})

ctdna_status = ctdna_status[ctdna_status["Timepoint"] == "Baseline"][["Patient_id", "ctDNA status"]]
ch_status = ch_status[ch_status["Timepoint"] == "Baseline"][["Patient_id", "CH status"]]

data = clin.merge(ctdna_status, how = "left").merge(ch_status, how = "left")
data = data[~data["ctDNA status"].isna()]

# Any CH in any gene, vs CH in DTA genes
run_fishers_exact(data, "ctDNA status", "irAE")

# LOGISTIC REGRESSION
for colname in ["CH status", "ctDNA status"]:
    data[colname] = data[colname].map({"Negative": 0, "Positive": 1})

del data["CH status"]

# Univar logistic regression
X = data["ctDNA status"]  # Predictor variables
X_incl_const = sm.add_constant(X)  # Add constant for intercept
y = data['irAE']  # Outcome variable

logit_model = sm.Logit(y, X_incl_const)
logit_result = logit_model.fit()
print("Coefficients:\n", logit_result.params) # Coefficients (log odds)
print("P-values:\n", logit_result.pvalues) # P-values
print("Confidence Intervals:\n", logit_result.conf_int()) # Confidence Intervals
print("Log-Likelihood:\n", logit_result.llf) # Log-likelihood value
print("Pseudo R-squared:\n", logit_result.prsquared) # Pseudo R-squared value (McFadden's R-squared)

# Multivar regression
X = data[["ctDNA status", "Age at GUBB draw", "OS from cfDNA collection (mo)"]]  # Predictor variables
X_incl_const = sm.add_constant(X)  # Add constant for intercept
y = data['irAE']  # Outcome variable

logit_model = sm.Logit(y, X_incl_const)
logit_result = logit_model.fit()
print("Coefficients:\n", logit_result.params) # Coefficients (log odds)
print("P-values:\n", logit_result.pvalues) # P-values
print("Confidence Intervals:\n", logit_result.conf_int()) # Confidence Intervals
print("Log-Likelihood:\n", logit_result.llf) # Log-likelihood value
print("Pseudo R-squared:\n", logit_result.prsquared) # Pseudo R-squared value (McFadden's R-squared)
